#include "mmappet.h"
#include "MSNumpress.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <atomic>
#include <thread>
#include <mutex>
#include <memory>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <zlib.h>

// ─── base64 encode ──────────────────────────────────────────────────────────
static const char b64_table[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

// ─── base64 encode into string ─────────────────────────────────────────────
static void base64_encode_into(std::string& out, const uint8_t* data, size_t len) {
    size_t b64_len = ((len + 2) / 3) * 4;
    out.resize(b64_len);
    char* p = out.data();
    for (size_t i = 0; i < len; i += 3) {
        uint32_t n = uint32_t(data[i]) << 16;
        if (i + 1 < len) n |= uint32_t(data[i + 1]) << 8;
        if (i + 2 < len) n |= uint32_t(data[i + 2]);
        *p++ = b64_table[(n >> 18) & 0x3F];
        *p++ = b64_table[(n >> 12) & 0x3F];
        *p++ = (i + 1 < len) ? b64_table[(n >> 6) & 0x3F] : '=';
        *p++ = (i + 2 < len) ? b64_table[n & 0x3F] : '=';
    }
}

// ─── number formatting (writes directly into buffer, returns chars written) ──
static int fmt_float_buf(char* buf, double val, int max_decimals = 6) {
    int n = snprintf(buf, 64, "%.*f", max_decimals, val);
    char* dot = (char*)memchr(buf, '.', n);
    if (dot) {
        char* end = buf + n - 1;
        while (end > dot + 1 && *end == '0') --end;
        n = (int)(end + 1 - buf);
        buf[n] = '\0';
    }
    return n;
}

// ─── spectrum XML template ─────────────────────────────────────────────────
// {} placeholders (22 total, in order):
//  0: spectrum sequential index, 1: precursor id
//  2:   frag_count
//  3-6: run_id, index, index, index
//  7-8: lowest_mz (3 dec), highest_mz (3 dec)
//  9:   tic (1 dec)
// 10-11: base_peak_mz (3 dec), base_peak_int (1 dec)
// 12-13: rt (6 dec), iim (4 dec)
// 14-15: prec_mz (6 dec), charge
// 16: mz_b64_len
// 17: mz encoding cvParam line (varies: 32-bit float or numpress linear)
// 18: mz_b64_data
// 19: int_b64_len
// 20: int encoding cvParam line (varies: 32-bit float or numpress pic)
// 21: int_b64_data
static constexpr int SPECTRUM_N_SLOTS = 22;
static const char SPECTRUM_TPL[] =
"        <spectrum index=\"{}\" id=\"{}\" defaultArrayLength=\"{}\">\n"
"          <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n"
"          <cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"2\"/>\n"
"          <cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\" value=\"\"/>\n"
"          <cvParam cvRef=\"MS\" accession=\"MS:1000796\" name=\"spectrum title\" value=\"{}.{}.{}.2 File:&quot;&quot;, NativeID:&quot;{}&quot;\"/>\n"
"          <cvParam cvRef=\"MS\" accession=\"MS:1000130\" name=\"positive scan\" value=\"\"/>\n"
"          <cvParam cvRef=\"MS\" accession=\"MS:1000528\" name=\"lowest observed m/z\" value=\"{}\"/>\n"
"          <cvParam cvRef=\"MS\" accession=\"MS:1000527\" name=\"highest observed m/z\" value=\"{}\"/>\n"
"          <cvParam cvRef=\"MS\" accession=\"MS:1000285\" name=\"total ion current\" value=\"{}\"/>\n"
"          <cvParam cvRef=\"MS\" accession=\"MS:1000504\" name=\"base peak m/z\" value=\"{}\"/>\n"
"          <cvParam cvRef=\"MS\" accession=\"MS:1000505\" name=\"base peak intensity\" value=\"{}\"/>\n"
"          <scanList count=\"1\">\n"
"            <cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\" value=\"\"/>\n"
"            <scan>\n"
"              <cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"{}\" unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"/>\n"
"              <cvParam cvRef=\"MS\" accession=\"MS:1002815\" name=\"inverse reduced ion mobility\" value=\"{}\" unitCvRef=\"MS\" unitAccession=\"MS:1002814\" unitName=\"volt-second per square centimeter\"/>\n"
"            </scan>\n"
"          </scanList>\n"
"          <precursorList count=\"1\">\n"
"            <precursor>\n"
"              <selectedIonList count=\"1\">\n"
"                <selectedIon>\n"
"                  <cvParam cvRef=\"MS\" accession=\"MS:1000744\" name=\"selected ion m/z\" value=\"{}\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n"
"                  <cvParam cvRef=\"MS\" accession=\"MS:1000041\" name=\"charge state\" value=\"{}\"/>\n"
"                </selectedIon>\n"
"              </selectedIonList>\n"
"              <activation>\n"
"              </activation>\n"
"            </precursor>\n"
"          </precursorList>\n"
"          <binaryDataArrayList count=\"2\">\n"
"            <binaryDataArray encodedLength=\"{}\">\n"
"{}"
"              <cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\" value=\"\"/>\n"
"              <cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\" value=\"\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n"
"              <binary>{}</binary>\n"
"            </binaryDataArray>\n"
"            <binaryDataArray encodedLength=\"{}\">\n"
"{}"
"              <cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\" value=\"\"/>\n"
"              <cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" value=\"\" unitCvRef=\"MS\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\"/>\n"
"              <binary>{}</binary>\n"
"            </binaryDataArray>\n"
"          </binaryDataArrayList>\n"
"        </spectrum>\n";

// Encoding cvParam lines (constant per run, selected by --numpress flag)
static const char CV_32BIT_FLOAT[] =
    "              <cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\" value=\"\"/>\n";
static const char CV_NUMPRESS_LINEAR[] =
    "              <cvParam cvRef=\"MS\" accession=\"MS:1002746\" name=\"MS-Numpress linear prediction compression\" value=\"\"/>\n";
static const char CV_NUMPRESS_PIC[] =
    "              <cvParam cvRef=\"MS\" accession=\"MS:1002747\" name=\"MS-Numpress positive integer compression\" value=\"\"/>\n";

// Parsed template: fixed text segments between {} placeholders
struct SpectrumTemplate {
    struct Seg { const char* data; size_t len; };
    Seg segs[SPECTRUM_N_SLOTS + 1];

    SpectrumTemplate() {
        const char* p = SPECTRUM_TPL;
        int i = 0;
        const char* start = p;
        while (*p) {
            if (p[0] == '{' && p[1] == '}') {
                segs[i++] = {start, (size_t)(p - start)};
                p += 2;
                start = p;
            } else {
                ++p;
            }
        }
        segs[i] = {start, (size_t)(p - start)};
    }
};
static SpectrumTemplate spectrum_tmpl;

// ─── reusable per-thread buffers ────────────────────────────────────────────
struct RenderBufs {
    std::vector<float>   int_f32;      // intensity float32 conversion
    std::vector<float>   mz_filtered;  // compacted filtered m/z values
    std::vector<uint32_t> int_filtered; // compacted filtered intensity values
    std::vector<uint8_t> zlib_out;     // reusable zlib output
    std::string          mz_b64;       // base64 encoded mz data
    std::string          int_b64;      // base64 encoded intensity data
    std::string          output;       // final XML string
    size_t               slot_offsets[SPECTRUM_N_SLOTS]; // byte offset of each slot value in output
    z_stream             zstrm;        // reusable zlib stream (avoids malloc/free per call)
    bool                 zstrm_init = false;
    // Numpress buffers
    std::vector<float>         mz_rounded; // m/z rounded to N decimals
    std::vector<double>        mz_f64;    // float32 → double for numpress
    std::vector<double>        int_f64;   // uint32 → double for numpress
    std::vector<unsigned char> np_buf;    // numpress output (max: 8 + N*5)

    void init_zlib(int level) {
        memset(&zstrm, 0, sizeof(zstrm));
        int ret = deflateInit2(&zstrm, level, Z_DEFLATED, 15, 8, Z_DEFAULT_STRATEGY);
        if (ret != Z_OK) { fprintf(stderr, "deflateInit2 failed: %d\n", ret); exit(1); }
        zstrm_init = true;
    }

    ~RenderBufs() {
        if (zstrm_init) deflateEnd(&zstrm);
    }
};

static void zlib_compress_into(RenderBufs& bufs, const void* data, size_t len) {
    z_stream& s = bufs.zstrm;
    deflateReset(&s);
    uLongf bound = deflateBound(&s, len);
    bufs.zlib_out.resize(bound);
    s.next_in = (Bytef*)data;
    s.avail_in = len;
    s.next_out = bufs.zlib_out.data();
    s.avail_out = bound;
    int ret = deflate(&s, Z_FINISH);
    if (ret != Z_STREAM_END) { fprintf(stderr, "deflate failed: %d\n", ret); exit(1); }
    bufs.zlib_out.resize(s.total_out);
}

static std::vector<uint8_t> unpack_lsb_bits(const std::filesystem::path& path, size_t n_bits) {
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("Cannot open bitpack file: " + path.string());
    std::vector<uint8_t> packed((std::istreambuf_iterator<char>(f)), {});
    size_t need = (n_bits + 7) / 8;
    if (packed.size() < need) throw std::runtime_error("Bitpack file too short: " + path.string());

    std::vector<uint8_t> bits(n_bits);
    for (size_t i = 0; i < n_bits; i++) bits[i] = (packed[i / 8] >> (i % 8)) & 1u;
    return bits;
}

static std::string metadata_value(const std::filesystem::path& path, const std::string& key) {
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("Cannot open filter metadata: " + path.string());
    std::string line;
    while (std::getline(f, line)) {
        size_t pos = line.find('=');
        if (pos == std::string::npos) continue;
        if (line.substr(0, pos) == key) return line.substr(pos + 1);
    }
    throw std::runtime_error("Missing metadata key '" + key + "' in " + path.string());
}

struct TofFilter {
    bool enabled = false;
    std::vector<uint8_t> precursor_keep;
    std::vector<uint8_t> fragment_keep;

    bool keep_fragment(uint64_t idx) const {
        return !enabled || (idx < fragment_keep.size() && fragment_keep[idx]);
    }
};

static TofFilter load_tof_filter(const std::filesystem::path& filter_path, size_t n_fragments) {
    TofFilter filter;
    if (filter_path.empty()) return filter;

    const auto meta_path = filter_path / "metadata.txt";
    const std::string encoding = metadata_value(meta_path, "encoding");
    if (encoding != "bitpacked_lsb_u8") {
        throw std::runtime_error("Unsupported tof filter encoding: " + encoding);
    }

    size_t n_filter_prec = std::stoull(metadata_value(meta_path, "n_precursors"));
    size_t n_filter_frag = std::stoull(metadata_value(meta_path, "n_fragments"));
    if (n_filter_frag != n_fragments) {
        std::ostringstream msg;
        msg << "tof filter fragment count mismatch: filter=" << n_filter_frag
            << " pmsms=" << n_fragments;
        throw std::runtime_error(msg.str());
    }

    filter.precursor_keep = unpack_lsb_bits(filter_path / "precursor_keep.bin", n_filter_prec);
    filter.fragment_keep = unpack_lsb_bits(filter_path / "fragment_keep.bin", n_filter_frag);
    filter.enabled = true;
    return filter;
}

static uint64_t count_kept_fragments(const TofFilter& filter, uint64_t offset, uint64_t count) {
    if (!filter.enabled) return count;
    uint64_t kept = 0;
    for (uint64_t j = 0; j < count; j++) kept += filter.keep_fragment(offset + j) ? 1 : 0;
    return kept;
}

// ─── render one spectrum into bufs.output ───────────────────────────────────
static void render_spectrum(
    RenderBufs& bufs,
    size_t seq_idx,
    size_t prec_idx,
    const char* run_id,
    size_t run_id_len,
    int zlib_level,
    bool use_numpress,
    int max_decimals,
    const char* mz_enc_cv, size_t mz_enc_cv_len,
    const char* int_enc_cv, size_t int_enc_cv_len,
    const float* frag_mz_data,
    const uint32_t* frag_tof_data,
    const float* tof2mz_data,
    bool use_tof2mz,
    const uint32_t* frag_int_data,
    uint64_t frag_offset,
    uint64_t frag_count,
    const TofFilter* tof_filter,
    double rt,
    double prec_mz,
    uint8_t charge,
    double iim,
    int index_width = 0)
{
    // Compute per-spectrum stats
    float lowest_mz = 0, highest_mz = 0;
    double tic = 0;
    float base_peak_mz = 0;
    float base_peak_int = 0;

    auto mz_at = [&](uint64_t idx) -> float {
        return use_tof2mz ? tof2mz_data[frag_tof_data[idx]] : frag_mz_data[idx];
    };

    const bool use_filter = tof_filter && tof_filter->enabled;
    if (use_filter) {
        bufs.mz_filtered.clear();
        bufs.int_filtered.clear();
        bufs.mz_filtered.reserve(frag_count);
        bufs.int_filtered.reserve(frag_count);
        for (uint64_t j = 0; j < frag_count; j++) {
            uint64_t idx = frag_offset + j;
            if (!tof_filter->keep_fragment(idx)) continue;
            bufs.mz_filtered.push_back(mz_at(idx));
            bufs.int_filtered.push_back(frag_int_data[idx]);
        }
        frag_count = bufs.mz_filtered.size();
    }

    auto kept_mz_at = [&](uint64_t j) -> float {
        return use_filter ? bufs.mz_filtered[j] : mz_at(frag_offset + j);
    };
    auto kept_int_at = [&](uint64_t j) -> uint32_t {
        return use_filter ? bufs.int_filtered[j] : frag_int_data[frag_offset + j];
    };

    if (frag_count > 0) {
        lowest_mz  = kept_mz_at(0);
        highest_mz = kept_mz_at(frag_count - 1);
        uint64_t best_j = 0;
        for (uint64_t j = 0; j < frag_count; j++) {
            float fi = (float)kept_int_at(j);
            tic += fi;
            if (fi > base_peak_int) { base_peak_int = fi; best_j = j; }
        }
        base_peak_mz = kept_mz_at(best_j);
    }

    if (!bufs.zstrm_init) bufs.init_zlib(zlib_level);

    // Truncate m/z values to N decimal places if requested (matches Python MGF pipeline)
    const float* mz_src = use_filter ? bufs.mz_filtered.data() : (use_tof2mz ? nullptr : &frag_mz_data[frag_offset]);
    if (use_filter || use_tof2mz || max_decimals > 0) {
        double scale = pow(10.0, max_decimals);
        bufs.mz_rounded.resize(frag_count);
        char tmp[64];
        for (uint64_t j = 0; j < frag_count; j++) {
            float mz = kept_mz_at(j);
            if (max_decimals > 0) {
                double truncated = floor((double)mz * scale) / scale;
                snprintf(tmp, sizeof(tmp), "%.*f", max_decimals, truncated);
                mz = strtof(tmp, nullptr);
            }
            bufs.mz_rounded[j] = mz;
        }
        mz_src = bufs.mz_rounded.data();
    }

    if (use_numpress) {
        using namespace ms::numpress::MSNumpress;
        // mz: float32 → double → numpress linear → zlib → base64
        bufs.mz_f64.resize(frag_count);
        for (uint64_t j = 0; j < frag_count; j++)
            bufs.mz_f64[j] = (double)mz_src[j];
        bufs.np_buf.resize(8 + frag_count * 5);
        double fp = optimalLinearFixedPoint(bufs.mz_f64.data(), frag_count);
        size_t np_len = encodeLinear(bufs.mz_f64.data(), frag_count, bufs.np_buf.data(), fp);
        zlib_compress_into(bufs, bufs.np_buf.data(), np_len);
        base64_encode_into(bufs.mz_b64, bufs.zlib_out.data(), bufs.zlib_out.size());

        // intensity: u32 → double → numpress pic → zlib → base64
        bufs.int_f64.resize(frag_count);
        for (uint64_t j = 0; j < frag_count; j++)
            bufs.int_f64[j] = (double)kept_int_at(j);
        bufs.np_buf.resize(frag_count * 5);
        np_len = encodePic(bufs.int_f64.data(), frag_count, bufs.np_buf.data());
        zlib_compress_into(bufs, bufs.np_buf.data(), np_len);
        base64_encode_into(bufs.int_b64, bufs.zlib_out.data(), bufs.zlib_out.size());
    } else {
        // mz: float32 → zlib → base64
        zlib_compress_into(bufs, mz_src, frag_count * sizeof(float));
        base64_encode_into(bufs.mz_b64, bufs.zlib_out.data(), bufs.zlib_out.size());

        // intensity: u32 → f32 → zlib → base64
        bufs.int_f32.resize(frag_count);
        for (uint64_t j = 0; j < frag_count; j++)
            bufs.int_f32[j] = (float)kept_int_at(j);
        zlib_compress_into(bufs, bufs.int_f32.data(), frag_count * sizeof(float));
        base64_encode_into(bufs.int_b64, bufs.zlib_out.data(), bufs.zlib_out.size());
    }

    // Format all scalar values
    char seq_s[32], idx_s[32], frag_s[32], charge_s[32], mz_b64l_s[32], int_b64l_s[32];
    char low_s[64], high_s[64], tic_s[64], bp_mz_s[64], bp_int_s[64];
    char rt_s[64], iim_s[64], prec_mz_s[64];

    int seq_n     = (index_width > 0)
        ? snprintf(seq_s, sizeof(seq_s), "%0*zu", index_width, seq_idx)
        : snprintf(seq_s, sizeof(seq_s), "%zu", seq_idx);
    int idx_n     = snprintf(idx_s, sizeof(idx_s), "%zu", prec_idx);
    int frag_n    = snprintf(frag_s, sizeof(frag_s), "%zu", (size_t)frag_count);
    int charge_n  = snprintf(charge_s, sizeof(charge_s), "%u", (unsigned)charge);
    int mz_b64l_n = snprintf(mz_b64l_s, sizeof(mz_b64l_s), "%zu", bufs.mz_b64.size());
    int int_b64l_n= snprintf(int_b64l_s, sizeof(int_b64l_s), "%zu", bufs.int_b64.size());
    int low_n     = fmt_float_buf(low_s, lowest_mz, 3);
    int high_n    = fmt_float_buf(high_s, highest_mz, 3);
    int tic_n     = fmt_float_buf(tic_s, tic, 1);
    int bp_mz_n   = fmt_float_buf(bp_mz_s, base_peak_mz, 3);
    int bp_int_n  = fmt_float_buf(bp_int_s, base_peak_int, 1);
    int rt_n      = fmt_float_buf(rt_s, rt);
    int iim_n     = fmt_float_buf(iim_s, iim, 4);
    int prec_mz_n = fmt_float_buf(prec_mz_s, prec_mz);

    // Map template slots to formatted values
    struct V { const char* d; size_t n; };
    // When index_width > 0 (mmap path), use padded seq for id/title slots too
    const char* id_s = (index_width > 0) ? seq_s : idx_s;
    size_t      id_n = (index_width > 0) ? (size_t)seq_n : (size_t)idx_n;
    V vals[SPECTRUM_N_SLOTS] = {
        {seq_s,     (size_t)seq_n},              //  0: spectrum sequential index
        {id_s,      id_n},                       //  1: id (seq when padded, prec_idx otherwise)
        {frag_s,    (size_t)frag_n},             //  2: frag_count
        {run_id,    run_id_len},                 //  3: run_id
        {id_s,      id_n},                       //  4: scan in title
        {id_s,      id_n},                       //  5: scan in title
        {id_s,      id_n},                       //  6: scan in NativeID
        {low_s,     (size_t)low_n},              //  7: lowest_mz
        {high_s,    (size_t)high_n},             //  8: highest_mz
        {tic_s,     (size_t)tic_n},              //  9: tic
        {bp_mz_s,   (size_t)bp_mz_n},           // 10: base_peak_mz
        {bp_int_s,  (size_t)bp_int_n},           // 11: base_peak_int
        {rt_s,      (size_t)rt_n},               // 12: rt
        {iim_s,     (size_t)iim_n},              // 13: iim
        {prec_mz_s, (size_t)prec_mz_n},         // 14: prec_mz
        {charge_s,  (size_t)charge_n},           // 15: charge
        {mz_b64l_s, (size_t)mz_b64l_n},         // 16: mz_b64_len
        {mz_enc_cv, mz_enc_cv_len},             // 17: mz encoding cvParam
        {bufs.mz_b64.data(),  bufs.mz_b64.size()},  // 18: mz_b64_data
        {int_b64l_s,(size_t)int_b64l_n},         // 19: int_b64_len
        {int_enc_cv, int_enc_cv_len},            // 20: int encoding cvParam
        {bufs.int_b64.data(), bufs.int_b64.size()},  // 21: int_b64_data
    };

    // Render: interleave template segments with values
    std::string& s = bufs.output;
    s.clear();
    size_t total = 0;
    for (int j = 0; j <= SPECTRUM_N_SLOTS; j++) total += spectrum_tmpl.segs[j].len;
    for (int j = 0; j < SPECTRUM_N_SLOTS; j++) total += vals[j].n;
    s.reserve(total);
    for (int j = 0; j < SPECTRUM_N_SLOTS; j++) {
        s.append(spectrum_tmpl.segs[j].data, spectrum_tmpl.segs[j].len);
        bufs.slot_offsets[j] = s.size();
        s.append(vals[j].d, vals[j].n);
    }
    s.append(spectrum_tmpl.segs[SPECTRUM_N_SLOTS].data, spectrum_tmpl.segs[SPECTRUM_N_SLOTS].len);
}

// ─── build header string ────────────────────────────────────────────────────
static std::string build_header(const std::string& run_id, size_t n_spectra, bool indexed) {
    std::string s;
    s.reserve(4096);
    s += "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
    if (indexed) {
        s += "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd\">\n";
        s += "  <mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\" id=\"" + run_id + "\" version=\"1.1.0\">\n";
    } else {
        s += "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\" id=\"" + run_id + "\" version=\"1.1.0\">\n";
    }
    s += "    <cvList count=\"2\">\n";
    s += "      <cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=\"4.1.41\" URI=\"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo\"/>\n";
    s += "      <cv id=\"UO\" fullName=\"Unit Ontology\" version=\"09:04:2014\" URI=\"https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo\"/>\n";
    s += "    </cvList>\n";
    s += "    <fileDescription>\n";
    s += "      <fileContent>\n";
    s += "        <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n";
    s += "        <cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\" value=\"\"/>\n";
    s += "      </fileContent>\n";
    s += "      <sourceFileList count=\"1\">\n";
    s += "        <sourceFile id=\"" + run_id + ".mgf\" name=\"" + run_id + ".mgf\" location=\"file:///\">\n";
    s += "          <cvParam cvRef=\"MS\" accession=\"MS:1000774\" name=\"multiple peak list nativeID format\" value=\"\"/>\n";
    s += "          <cvParam cvRef=\"MS\" accession=\"MS:1001062\" name=\"Mascot MGF format\" value=\"\"/>\n";
    s += "        </sourceFile>\n";
    s += "      </sourceFileList>\n";
    s += "    </fileDescription>\n";
    s += "    <softwareList count=\"1\">\n";
    s += "      <software id=\"mgf_to_mzml_python\" version=\"1.0\">\n";
    s += "        <cvParam cvRef=\"MS\" accession=\"MS:1000615\" name=\"ProteoWizard software\" value=\"\"/>\n";
    s += "      </software>\n";
    s += "    </softwareList>\n";
    s += "    <instrumentConfigurationList count=\"1\">\n";
    s += "      <instrumentConfiguration id=\"IC\">\n";
    s += "        <cvParam cvRef=\"MS\" accession=\"MS:1000031\" name=\"instrument model\" value=\"\"/>\n";
    s += "      </instrumentConfiguration>\n";
    s += "    </instrumentConfigurationList>\n";
    s += "    <dataProcessingList count=\"1\">\n";
    s += "      <dataProcessing id=\"python_conversion\">\n";
    s += "        <processingMethod order=\"0\" softwareRef=\"mgf_to_mzml_python\">\n";
    s += "          <cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\" value=\"\"/>\n";
    s += "        </processingMethod>\n";
    s += "      </dataProcessing>\n";
    s += "    </dataProcessingList>\n";
    s += "    <run id=\"" + run_id + "\" defaultInstrumentConfigurationRef=\"IC\">\n";
    char tmp[128];
    snprintf(tmp, sizeof(tmp), "      <spectrumList count=\"%zu\" defaultDataProcessingRef=\"python_conversion\">\n", n_spectra);
    s += tmp;
    return s;
}

// ─── build footer string ────────────────────────────────────────────────────
// footer_start_offset = byte position in the file where this footer begins
static std::string build_footer(const std::vector<long>& offsets, long footer_start_offset, bool indexed) {
    std::string s;
    s += "      </spectrumList>\n";
    s += "    </run>\n";
    if (!indexed) {
        s += "  </mzML>\n";
        return s;
    }
    s += "  </mzML>\n";

    long index_list_offset = footer_start_offset + (long)s.size();

    s += "  <indexList count=\"2\">\n";
    s += "    <index name=\"spectrum\">\n";
    char tmp[128];
    for (size_t i = 0; i < offsets.size(); i++) {
        snprintf(tmp, sizeof(tmp), "      <offset idRef=\"index=%zu\">%ld</offset>\n", i, offsets[i]);
        s += tmp;
    }
    s += "    </index>\n";
    s += "    <index name=\"chromatogram\">\n";
    s += "    </index>\n";
    s += "  </indexList>\n";
    snprintf(tmp, sizeof(tmp), "  <indexListOffset>%ld</indexListOffset>\n", index_list_offset);
    s += tmp;
    s += "</indexedmzML>\n";
    return s;
}

// ─── multi-charge helpers ────────────────────────────────────────────────────
static std::vector<uint8_t> charges_from_int64(int64_t val) {
    if (val <= 0) return {};
    std::string s = std::to_string(val);
    std::vector<uint8_t> out;
    for (char c : s) out.push_back((uint8_t)(c - '0'));
    return out;
}

static bool detect_multicharge(const std::filesystem::path& prec_dir) {
    std::ifstream f(prec_dir / "schema.txt");
    std::string line;
    while (std::getline(f, line))
        if (line == "int64 charges") return true;
    return false;
}

// ─── main ───────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    // Parse CLI args
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <pmsms_dir> <precursors_dir> <output.mzml> [--tof2mz DIR] [--tof_filter_path DIR] [--run-id NAME] [--zlib-level N] [--threads N] [--decimals N] [--dry-run] [--check-mz-sorted] [--numpress] [--indexed] [--used-spectra-cnt N]\n", argv[0]);
        return 1;
    }

    auto runtime_start = std::chrono::steady_clock::now();


    std::filesystem::path pmsms_dir = argv[1];
    std::filesystem::path precursors_dir = argv[2];
    std::filesystem::path output_path = argv[3];
    std::string run_id = output_path.stem().string();
    int zlib_level = 1;
    int thread_cnt = 1;
    bool dry_run = false;
    bool check_mz_sorted = false;
    bool use_numpress = false;
    bool indexed = false;
    int max_decimals = 0;  // 0 = no rounding of binary m/z values
    size_t used_spectra_cnt = 0;  // 0 = use all
    bool use_tof2mz = false;
    std::filesystem::path tof2mz_path;
    std::filesystem::path tof_filter_path;

    for (int i = 4; i < argc; i++) {
        if (strcmp(argv[i], "--tof2mz") == 0 && i + 1 < argc) {
            tof2mz_path = argv[++i];
            use_tof2mz = true;
        } else if (strcmp(argv[i], "--tof_filter_path") == 0 && i + 1 < argc) {
            tof_filter_path = argv[++i];
        } else if (strcmp(argv[i], "--run-id") == 0 && i + 1 < argc) {
            run_id = argv[++i];
        } else if (strcmp(argv[i], "--zlib-level") == 0 && i + 1 < argc) {
            zlib_level = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            thread_cnt = atoi(argv[++i]);
            if (thread_cnt < 1) thread_cnt = 1;
        } else if (strcmp(argv[i], "--dry-run") == 0) {
            dry_run = true;
        } else if (strcmp(argv[i], "--check-mz-sorted") == 0) {
            check_mz_sorted = true;
        } else if (strcmp(argv[i], "--numpress") == 0) {
            use_numpress = true;
        } else if (strcmp(argv[i], "--indexed") == 0) {
            indexed = true;
        } else if (strcmp(argv[i], "--decimals") == 0 && i + 1 < argc) {
            max_decimals = atoi(argv[++i]);
            if (max_decimals < 1) max_decimals = 1;
        } else if (strcmp(argv[i], "--used-spectra-cnt") == 0 && i + 1 < argc) {
            used_spectra_cnt = (size_t)atol(argv[++i]);
        }
    }

    // ─── Schema detection ────────────────────────────────────────────────
    bool multicharge = detect_multicharge(precursors_dir);

    // ─── Step 1: Open datasets ──────────────────────────────────────────
    fprintf(stderr, "Opening datasets...\n");

    std::unique_ptr<MMappedData<float>> frag_mz_col;
    std::unique_ptr<MMappedData<uint32_t>> frag_tof_col;
    std::unique_ptr<MMappedData<float>> tof2mz_col;
    if (use_tof2mz) {
        frag_tof_col = std::make_unique<MMappedData<uint32_t>>(OpenColumn<uint32_t>(pmsms_dir, "tof"));
        tof2mz_col = std::make_unique<MMappedData<float>>(OpenColumn<float>(tof2mz_path, "x"));
    } else {
        frag_mz_col = std::make_unique<MMappedData<float>>(OpenColumn<float>(pmsms_dir, "mz"));
    }
    auto frag_int_col = OpenColumn<uint32_t>(pmsms_dir, "intensity");
    const float*    frag_mz_data  = use_tof2mz ? nullptr : frag_mz_col->data();
    const uint32_t* frag_tof_data = use_tof2mz ? frag_tof_col->data() : nullptr;
    const float*    tof2mz_data   = use_tof2mz ? tof2mz_col->data() : nullptr;
    const uint32_t* frag_int_data = frag_int_col.data();
    size_t n_fragments = frag_int_col.size();
    TofFilter tof_filter = load_tof_filter(tof_filter_path, n_fragments);

    if (use_tof2mz) {
        madvise((void*)frag_tof_data, frag_tof_col->size() * sizeof(uint32_t), MADV_SEQUENTIAL);
        madvise((void*)tof2mz_data,   tof2mz_col->size()   * sizeof(float),    MADV_SEQUENTIAL);
        if (tof2mz_col->size() == 0) {
            fprintf(stderr, "ERROR: tof2mz table is empty: %s\n", tof2mz_path.c_str());
            return 1;
        }
        for (size_t i = 1; i < tof2mz_col->size(); i++) {
            if (tof2mz_data[i] < tof2mz_data[i - 1]) {
                fprintf(stderr, "ERROR: tof2mz table is not non-decreasing at tof %zu: %.9g > %.9g\n",
                    i, (double)tof2mz_data[i - 1], (double)tof2mz_data[i]);
                return 1;
            }
        }
        for (size_t i = 0; i < frag_tof_col->size(); i++) {
            if ((size_t)frag_tof_data[i] >= tof2mz_col->size()) {
                fprintf(stderr, "ERROR: fragment tof %u at row %zu is out of bounds for tof2mz length %zu\n",
                    frag_tof_data[i], i, tof2mz_col->size());
                return 1;
            }
        }
    } else {
        madvise((void*)frag_mz_data,  frag_mz_col->size()  * sizeof(float),    MADV_SEQUENTIAL);
    }
    madvise((void*)frag_int_data, frag_int_col.size()  * sizeof(uint32_t), MADV_SEQUENTIAL);

    // Flat expanded precursor/charge index vectors (length == n_spectra after expansion)
    std::vector<double>   prec_iim_vec, prec_mz_vec, prec_rt_vec;
    std::vector<uint64_t> frag_start_vec, frag_cnt_vec;
    std::vector<size_t>   effective_prec_idx;
    std::vector<uint8_t>  effective_charge;

    {
        auto prec_iim_col   = OpenColumn<double>  (precursors_dir, "inv_ion_mobility");
        auto prec_mz_col    = OpenColumn<double>  (precursors_dir, "mz");
        auto prec_rt_col    = OpenColumn<double>  (precursors_dir, "rt");
        auto frag_start_col = OpenColumn<uint64_t>(precursors_dir, "fragment_spectrum_start");
        auto frag_cnt_col   = OpenColumn<uint64_t>(precursors_dir, "fragment_event_cnt");

        size_t n_prec = prec_iim_col.size();
        prec_iim_vec.assign  (prec_iim_col.data(),   prec_iim_col.data()   + n_prec);
        prec_mz_vec.assign   (prec_mz_col.data(),    prec_mz_col.data()    + n_prec);
        prec_rt_vec.assign   (prec_rt_col.data(),    prec_rt_col.data()    + n_prec);
        frag_start_vec.assign(frag_start_col.data(), frag_start_col.data() + n_prec);
        frag_cnt_vec.assign  (frag_cnt_col.data(),   frag_cnt_col.data()   + n_prec);

        std::vector<uint8_t> precursor_keep_by_row(n_prec, 1);
        if (tof_filter.enabled) {
            auto didx_idx_col = OpenColumn<uint64_t>(pmsms_dir / "dataindex.mmappet", "idx");
            if (didx_idx_col.size() != tof_filter.precursor_keep.size()) {
                std::ostringstream msg;
                msg << "tof filter precursor count mismatch: filter="
                    << tof_filter.precursor_keep.size()
                    << " dataindex=" << didx_idx_col.size();
                throw std::runtime_error(msg.str());
            }

            size_t didx_j = 0;
            for (size_t pi = 0; pi < n_prec; pi++) {
                const uint64_t frag_start = frag_start_vec[pi];
                while (didx_j < didx_idx_col.size() && didx_idx_col[didx_j] < frag_start) didx_j++;
                if (didx_j == didx_idx_col.size() || didx_idx_col[didx_j] != frag_start) {
                    std::ostringstream msg;
                    msg << "No dataindex row for precursor row " << pi
                        << " fragment_spectrum_start=" << frag_start;
                    throw std::runtime_error(msg.str());
                }
                precursor_keep_by_row[pi] = tof_filter.precursor_keep[didx_j];
            }
        }

        size_t n = (used_spectra_cnt > 0 && used_spectra_cnt < n_prec) ? used_spectra_cnt : n_prec;

        if (!multicharge) {
            auto prec_charge_col = OpenColumn<uint8_t>(precursors_dir, "charge");
            for (size_t i = 0; i < n; i++) {
                if (!precursor_keep_by_row[i]) continue;
                effective_prec_idx.push_back(i);
                effective_charge.push_back(prec_charge_col[i]);
            }
        } else {
            auto prec_charges_col = OpenColumn<int64_t>(precursors_dir, "charges");
            for (size_t pi = 0; pi < n; pi++) {
                if (!precursor_keep_by_row[pi]) continue;
                for (uint8_t c : charges_from_int64(prec_charges_col[pi])) {
                    effective_prec_idx.push_back(pi);
                    effective_charge.push_back(c);
                }
            }
        }
    }

    size_t n_spectra = effective_prec_idx.size();

    fprintf(stderr, "Spectra: %zu, Fragments: %zu, Threads: %d\n", n_spectra, n_fragments, thread_cnt);

    // ─── Build header ───────────────────────────────────────────────────
    std::string header = build_header(run_id, n_spectra, indexed);
    size_t header_size = header.size();

    // ─── Progress tracking ──────────────────────────────────────────────
    std::atomic<size_t> progress{0};
    size_t progress_step = n_spectra < 100 ? 1 : n_spectra / 100;

    const char* run_id_cstr = run_id.c_str();
    size_t run_id_len = run_id.size();

    // Select encoding cvParam strings based on --numpress
    const char* mz_enc_cv  = use_numpress ? CV_NUMPRESS_LINEAR : CV_32BIT_FLOAT;
    size_t mz_enc_cv_len   = use_numpress ? sizeof(CV_NUMPRESS_LINEAR) - 1 : sizeof(CV_32BIT_FLOAT) - 1;
    const char* int_enc_cv = use_numpress ? CV_NUMPRESS_PIC : CV_32BIT_FLOAT;
    size_t int_enc_cv_len  = use_numpress ? sizeof(CV_NUMPRESS_PIC) - 1 : sizeof(CV_32BIT_FLOAT) - 1;

    // Optional: verify all spectra have m/z sorted in ascending order.
    // In tof2mz mode, non-decreasing TOF is sufficient because tof2mz is
    // checked non-decreasing once above.
    if (check_mz_sorted) {
        fprintf(stderr, "Checking m/z sort order...\n");
        for (size_t i = 0; i < n_spectra; i++) {
            size_t pi = effective_prec_idx[i];
            uint64_t off = frag_start_vec[pi];
            uint64_t cnt = frag_cnt_vec[pi];
            for (uint64_t j = 1; j < cnt; j++) {
                if (use_tof2mz) {
                    if (frag_tof_data[off + j] < frag_tof_data[off + j - 1]) {
                        fprintf(stderr, "ERROR: spectrum %zu (precursor %zu) has unsorted TOF at fragment %zu: %u > %u\n",
                            i, pi, (size_t)j, frag_tof_data[off + j - 1], frag_tof_data[off + j]);
                        return 1;
                    }
                } else if (frag_mz_data[off + j] < frag_mz_data[off + j - 1]) {
                    fprintf(stderr, "ERROR: spectrum %zu (precursor %zu) has unsorted m/z at fragment %zu: %.6f > %.6f\n",
                        i, pi, (size_t)j, (double)frag_mz_data[off + j - 1], (double)frag_mz_data[off + j]);
                    return 1;
                }
            }
        }
        fprintf(stderr, "All %zu spectra have sorted m/z values.\n", n_spectra);
    }

    if (dry_run) {
        // ─── Dry-run: compute total size, no file I/O ──────────────────
        std::atomic<size_t> total_bytes{header_size};

        auto worker = [&](size_t begin, size_t end) {
            RenderBufs bufs;
            for (size_t i = begin; i < end; i++) {
                size_t pi = effective_prec_idx[i];
                render_spectrum(bufs, i, pi, run_id_cstr, run_id_len, zlib_level,
                    use_numpress, max_decimals, mz_enc_cv, mz_enc_cv_len, int_enc_cv, int_enc_cv_len,
                    frag_mz_data, frag_tof_data, tof2mz_data, use_tof2mz, frag_int_data,
                    frag_start_vec[pi], frag_cnt_vec[pi], &tof_filter,
                    prec_rt_vec[pi], prec_mz_vec[pi], effective_charge[i], prec_iim_vec[pi]);
                total_bytes.fetch_add(bufs.output.size(), std::memory_order_relaxed);

                size_t p = progress.fetch_add(1) + 1;
                if (p % progress_step == 0)
                    fprintf(stderr, "\rRendering spectra: %zu / %zu (%zu%%)", p, n_spectra, p * 100 / n_spectra);
            }
        };

        std::vector<std::thread> threads;
        size_t chunk = n_spectra / thread_cnt;
        size_t remainder = n_spectra % thread_cnt;
        size_t start = 0;
        for (int t = 0; t < thread_cnt; t++) {
            size_t end = start + chunk + (t < (int)remainder ? 1 : 0);
            threads.emplace_back(worker, start, end);
            start = end;
        }
        for (auto& th : threads) th.join();

        fprintf(stderr, "\rRendering spectra: %zu / %zu (100%%)\n", n_spectra, n_spectra);
        fprintf(stderr, "Dry run: %zu bytes (%.2f GiB) for %zu spectra\n",
            total_bytes.load(), total_bytes.load() / (1024.0 * 1024.0 * 1024.0), n_spectra);
        return 0;

    } else if (!indexed) {
        // ─── Multi-threaded mmap path (no indexing, lock-free) ───────────
        // Estimate file size: header + per-spectrum XML + base64(zlib(fragments)) + footer
        // Base64 of zlib: worst case zlib expands to input+overhead, base64 adds 33%
        // Per fragment: 4 bytes (mz) + 4 bytes (int) = 8 bytes raw → ~11 bytes base64
        size_t total_frags = 0;
        for (size_t i = 0; i < n_spectra; i++)
            total_frags += count_kept_fragments(
                tof_filter,
                frag_start_vec[effective_prec_idx[i]],
                frag_cnt_vec[effective_prec_idx[i]]
            );
        size_t est_size = (header_size + n_spectra * 2000 + total_frags * 12 + 4096) * 2;

        std::filesystem::create_directories(output_path);
        std::filesystem::path mzml_path = output_path / "mzml.mzML";

        int fd = ::open(mzml_path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
        if (fd < 0) { fprintf(stderr, "Cannot open output: %s\n", mzml_path.c_str()); return 1; }
        ftruncate(fd, est_size);

        uint8_t* map = (uint8_t*)mmap(nullptr, est_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (map == MAP_FAILED) { fprintf(stderr, "mmap failed: %s\n", strerror(errno)); ::close(fd); return 1; }

        // Write header
        memcpy(map, header.data(), header_size);
        std::mutex order_lock;
        size_t write_pos = header_size;
        size_t seq_counter = 0;
        int index_width = (n_spectra > 0) ? snprintf(nullptr, 0, "%zu", n_spectra - 1) : 1;
        std::vector<size_t> seq_to_prec(n_spectra);  // mapping: sequential index → precursor index

        auto worker = [&](size_t begin, size_t end) {
            RenderBufs bufs;
            for (size_t i = begin; i < end; i++) {
                size_t pi = effective_prec_idx[i];
                render_spectrum(bufs, 0, pi, run_id_cstr, run_id_len, zlib_level,
                    use_numpress, max_decimals, mz_enc_cv, mz_enc_cv_len, int_enc_cv, int_enc_cv_len,
                    frag_mz_data, frag_tof_data, tof2mz_data, use_tof2mz, frag_int_data,
                    frag_start_vec[pi], frag_cnt_vec[pi], &tof_filter,
                    prec_rt_vec[pi], prec_mz_vec[pi], effective_charge[i], prec_iim_vec[pi],
                    index_width);

                size_t sz = bufs.output.size();
                size_t pos, seq;
                {
                    std::lock_guard<std::mutex> lk(order_lock);
                    pos = write_pos;
                    write_pos += sz;
                    seq = seq_counter++;
                }

                // Patch zero-padded sequential index into all id/index slots
                char seq_buf[32];
                snprintf(seq_buf, sizeof(seq_buf), "%0*zu", index_width, seq);
                for (int slot : {0, 1, 4, 5, 6})
                    memcpy(bufs.output.data() + bufs.slot_offsets[slot], seq_buf, index_width);

                memcpy(map + pos, bufs.output.data(), sz);
                seq_to_prec[seq] = pi;

                size_t p = progress.fetch_add(1) + 1;
                if (p % progress_step == 0)
                    fprintf(stderr, "\rWriting spectra: %zu / %zu (%zu%%)", p, n_spectra, p * 100 / n_spectra);
            }
        };

        std::vector<std::thread> threads;
        size_t chunk = n_spectra / thread_cnt;
        size_t remainder = n_spectra % thread_cnt;
        size_t start = 0;
        for (int t = 0; t < thread_cnt; t++) {
            size_t end = start + chunk + (t < (int)remainder ? 1 : 0);
            threads.emplace_back(worker, start, end);
            start = end;
        }
        for (auto& th : threads) th.join();

        fprintf(stderr, "\rWriting spectra: %zu / %zu (100%%)\n", n_spectra, n_spectra);

        // Write footer
        size_t body_end = write_pos;
        std::vector<long> no_offsets;
        std::string footer = build_footer(no_offsets, 0, false);
        memcpy(map + body_end, footer.data(), footer.size());
        size_t actual_size = body_end + footer.size();

        munmap(map, est_size);
        ftruncate(fd, actual_size);
        ::close(fd);

        // Write id mapping as mmappet (seq_idx → prec_idx)
        std::filesystem::path map_dir = output_path / "idmap.mmappet";
        std::filesystem::create_directories(map_dir);
        {
            FILE* sf = fopen((map_dir / "schema.txt").c_str(), "w");
            fprintf(sf, "uint64 prec_idx\n");
            fclose(sf);

            FILE* bf = fopen((map_dir / "0.bin").c_str(), "wb");
            std::vector<uint64_t> col(n_spectra);
            for (size_t i = 0; i < n_spectra; i++) col[i] = (uint64_t)seq_to_prec[i];
            fwrite(col.data(), sizeof(uint64_t), n_spectra, bf);
            fclose(bf);
        }
        fprintf(stderr, "Wrote id mapping to %s\n", map_dir.c_str());
        fprintf(stderr, "Done. Wrote %zu spectra to %s\n", n_spectra, output_path.c_str());

    } else {
        // ─── Multi-threaded indexed path (two-phase: parallel render, ordered write) ─
        int fd = ::open(output_path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
        if (fd < 0) { fprintf(stderr, "Cannot open output: %s\n", output_path.c_str()); return 1; }

        // Write header
        ::pwrite(fd, header.data(), header.size(), 0);

        std::vector<long> spectrum_offsets(n_spectra);

        // Per-thread result: rendered output + local spectrum offsets
        struct ThreadResult {
            std::string wbuf;
            std::vector<std::pair<size_t, size_t>> offsets; // (spectrum_idx, local_offset)
        };
        std::vector<ThreadResult> results(thread_cnt);

        // Phase 1: parallel render
        auto worker = [&](int tid, size_t begin, size_t end) {
            RenderBufs bufs;
            ThreadResult& res = results[tid];
            res.offsets.reserve(end - begin);

            for (size_t i = begin; i < end; i++) {
                size_t pi = effective_prec_idx[i];
                render_spectrum(bufs, i, pi, run_id_cstr, run_id_len, zlib_level,
                    use_numpress, max_decimals, mz_enc_cv, mz_enc_cv_len, int_enc_cv, int_enc_cv_len,
                    frag_mz_data, frag_tof_data, tof2mz_data, use_tof2mz, frag_int_data,
                    frag_start_vec[pi], frag_cnt_vec[pi], &tof_filter,
                    prec_rt_vec[pi], prec_mz_vec[pi], effective_charge[i], prec_iim_vec[pi]);

                res.offsets.push_back({i, res.wbuf.size()});
                res.wbuf.append(bufs.output);

                size_t p = progress.fetch_add(1) + 1;
                if (p % progress_step == 0)
                    fprintf(stderr, "\rRendering spectra: %zu / %zu (%zu%%)", p, n_spectra, p * 100 / n_spectra);
            }
        };

        std::vector<std::thread> threads;
        size_t chunk = n_spectra / thread_cnt;
        size_t remainder = n_spectra % thread_cnt;
        size_t start = 0;
        for (int t = 0; t < thread_cnt; t++) {
            size_t end = start + chunk + (t < (int)remainder ? 1 : 0);
            threads.emplace_back(worker, t, start, end);
            start = end;
        }
        for (auto& th : threads) th.join();

        fprintf(stderr, "\rRendering spectra: %zu / %zu (100%%)\n", n_spectra, n_spectra);

        // Phase 2: ordered write (thread 0 first, then thread 1, ...)
        size_t body_pos = 0;
        for (int t = 0; t < thread_cnt; t++) {
            ThreadResult& res = results[t];
            size_t thread_start = header_size + body_pos;

            // Resolve spectrum offsets
            for (auto& [idx, local_off] : res.offsets)
                spectrum_offsets[idx] = (long)(thread_start + local_off);

            // Write this thread's buffer
            size_t written = 0;
            while (written < res.wbuf.size()) {
                ssize_t ret = ::pwrite(fd, res.wbuf.data() + written,
                                       res.wbuf.size() - written,
                                       thread_start + written);
                if (ret < 0) {
                    fprintf(stderr, "pwrite failed: %s\n", strerror(errno));
                    exit(1);
                }
                written += ret;
            }
            body_pos += res.wbuf.size();

            // Free memory as we go
            res.wbuf.clear();
            res.wbuf.shrink_to_fit();
        }

        // Write footer
        off_t footer_off = header_size + body_pos;
        std::string footer = build_footer(spectrum_offsets, (long)footer_off, true);

        size_t written = 0;
        while (written < footer.size()) {
            ssize_t ret = ::pwrite(fd, footer.data() + written, footer.size() - written,
                                   footer_off + written);
            if (ret < 0) {
                fprintf(stderr, "pwrite failed: %s\n", strerror(errno));
                exit(1);
            }
            written += ret;
        }

        ftruncate(fd, footer_off + footer.size());

        ::close(fd);
        fprintf(stderr, "Done. Wrote %zu spectra to %s\n", n_spectra, output_path.c_str());
    }

    auto runtime_end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(runtime_end - runtime_start);

    fprintf(stdout, "Processing completed successfully in %.2f seconds.\n", duration.count());

    return 0;
}
