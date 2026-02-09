#include "mmappet.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdarg>
#include <string>
#include <vector>
#include <filesystem>
#include <atomic>
#include <thread>
#include <zlib.h>

// ─── base64 encode ──────────────────────────────────────────────────────────
static const char b64_table[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

// (base64 and zlib functions are inlined in render_spectrum via _into variants)

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

// ─── reusable per-thread buffers ────────────────────────────────────────────
struct RenderBufs {
    std::vector<float>   int_f32;      // intensity float32 conversion
    std::vector<uint8_t> zlib_out;     // reusable zlib output
    std::string          output;       // final XML string
    z_stream             zstrm;        // reusable zlib stream (avoids malloc/free per call)
    bool                 zstrm_init = false;

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

// Appends a string literal to the output
static inline void emit(std::string& s, const char* literal) {
    s.append(literal);
}
static inline void emit_fmt(std::string& s, char* tmp, size_t tmpn, const char* fmt, ...) __attribute__((format(printf, 4, 5)));
static inline void emit_fmt(std::string& s, char* tmp, size_t tmpn, const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    int n = vsnprintf(tmp, tmpn, fmt, args);
    va_end(args);
    s.append(tmp, n);
}

// ─── render one spectrum into bufs.output ───────────────────────────────────
static void render_spectrum(
    RenderBufs& bufs,
    size_t i,
    const char* run_id,
    int zlib_level,
    const float* frag_mz_data,
    const uint32_t* frag_int_data,
    uint64_t frag_offset,
    uint64_t frag_count,
    double rt,
    double prec_mz,
    uint8_t charge,
    double iim)
{
    // Compute per-spectrum stats
    float lowest_mz = 0, highest_mz = 0;
    double tic = 0;
    float base_peak_mz = 0;
    float base_peak_int = 0;

    for (uint64_t j = 0; j < frag_count; j++) {
        float mz = frag_mz_data[frag_offset + j];
        float fi = (float)frag_int_data[frag_offset + j];
        if (j == 0 || mz < lowest_mz) lowest_mz = mz;
        if (j == 0 || mz > highest_mz) highest_mz = mz;
        tic += fi;
        if (j == 0 || fi > base_peak_int) {
            base_peak_int = fi;
            base_peak_mz = mz;
        }
    }

    // Encode mz array: float32 bytes → zlib → base64
    // Reuse bufs.zlib_out for both compressions sequentially
    if (!bufs.zstrm_init) bufs.init_zlib(zlib_level);
    zlib_compress_into(bufs, &frag_mz_data[frag_offset], frag_count * sizeof(float));

    // We need both b64 strings alive at the same time for the XML, so we store
    // mz_b64 in the first part of output, but actually it's simpler to just
    // build the XML in two phases. Instead, let's use a small trick:
    // write XML up to the mz binary, append mz b64 inline, then compress intensity.

    std::string& s = bufs.output;
    s.clear();
    // Reserve generous estimate: 2KB XML + b64 sizes (upper bound)
    size_t b64_bound = ((compressBound(frag_count * 4) + 2) / 3) * 4;
    s.reserve(2048 + b64_bound * 2);

    char tmp[512];
    char fbuf[64];

    emit_fmt(s, tmp, sizeof(tmp),
        "<spectrum index=\"%zu\" id=\"index=%zu\" defaultArrayLength=\"%zu\">\n", i, i, (size_t)frag_count);
    emit(s, "<cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n");
    emit(s, "<cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"2\"/>\n");
    emit(s, "<cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\" value=\"\"/>\n");
    emit_fmt(s, tmp, sizeof(tmp),
        "<cvParam cvRef=\"MS\" accession=\"MS:1000796\" name=\"spectrum title\" value=\"%s.%zu.%zu.2 File:&quot;&quot;, NativeID:&quot;index=%zu&quot;\"/>\n",
        run_id, i, i, i);
    emit(s, "<cvParam cvRef=\"MS\" accession=\"MS:1000130\" name=\"positive scan\" value=\"\"/>\n");

    emit(s, "<cvParam cvRef=\"MS\" accession=\"MS:1000528\" name=\"lowest observed m/z\" value=\"");
    s.append(fbuf, fmt_float_buf(fbuf, lowest_mz, 3));
    emit(s, "\"/>\n");

    emit(s, "<cvParam cvRef=\"MS\" accession=\"MS:1000527\" name=\"highest observed m/z\" value=\"");
    s.append(fbuf, fmt_float_buf(fbuf, highest_mz, 3));
    emit(s, "\"/>\n");

    emit(s, "<cvParam cvRef=\"MS\" accession=\"MS:1000285\" name=\"total ion current\" value=\"");
    s.append(fbuf, fmt_float_buf(fbuf, tic, 1));
    emit(s, "\"/>\n");

    emit(s, "<cvParam cvRef=\"MS\" accession=\"MS:1000504\" name=\"base peak m/z\" value=\"");
    s.append(fbuf, fmt_float_buf(fbuf, base_peak_mz, 3));
    emit(s, "\"/>\n");

    emit(s, "<cvParam cvRef=\"MS\" accession=\"MS:1000505\" name=\"base peak intensity\" value=\"");
    s.append(fbuf, fmt_float_buf(fbuf, base_peak_int, 1));
    emit(s, "\"/>\n");

    emit(s, "<scanList count=\"1\">\n"
            "<cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\" value=\"\"/>\n"
            "<scan>\n"
            "<cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"");
    s.append(fbuf, fmt_float_buf(fbuf, rt));
    emit(s, "\" unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"/>\n"
            "<cvParam cvRef=\"MS\" accession=\"MS:1002815\" name=\"inverse reduced ion mobility\" value=\"");
    s.append(fbuf, fmt_float_buf(fbuf, iim, 4));
    emit(s, "\" unitCvRef=\"MS\" unitAccession=\"MS:1002814\" unitName=\"volt-second per square centimeter\"/>\n"
            "</scan>\n"
            "</scanList>\n"
            "<precursorList count=\"1\">\n"
            "<precursor>\n"
            "<selectedIonList count=\"1\">\n"
            "<selectedIon>\n"
            "<cvParam cvRef=\"MS\" accession=\"MS:1000744\" name=\"selected ion m/z\" value=\"");
    s.append(fbuf, fmt_float_buf(fbuf, prec_mz));
    emit(s, "\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n");
    emit_fmt(s, tmp, sizeof(tmp),
        "<cvParam cvRef=\"MS\" accession=\"MS:1000041\" name=\"charge state\" value=\"%u\"/>\n", (unsigned)charge);
    emit(s, "</selectedIon>\n"
            "</selectedIonList>\n"
            "<activation>\n"
            "</activation>\n"
            "</precursor>\n"
            "</precursorList>\n"
            "<binaryDataArrayList count=\"2\">\n");

    // m/z binary — bufs.zlib_out still has the compressed mz data
    size_t mz_b64_len = ((bufs.zlib_out.size() + 2) / 3) * 4;
    emit_fmt(s, tmp, sizeof(tmp), "<binaryDataArray encodedLength=\"%zu\">\n", mz_b64_len);
    emit(s, "<cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\" value=\"\"/>\n"
            "<cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\" value=\"\"/>\n"
            "<cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\" value=\"\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n"
            "<binary>");
    // Encode base64 directly into s
    size_t b64_start = s.size();
    s.resize(b64_start + mz_b64_len);
    {   // inline base64 encode into s
        const uint8_t* d = bufs.zlib_out.data();
        size_t dlen = bufs.zlib_out.size();
        char* p = s.data() + b64_start;
        for (size_t bi = 0; bi < dlen; bi += 3) {
            uint32_t n = uint32_t(d[bi]) << 16;
            if (bi + 1 < dlen) n |= uint32_t(d[bi + 1]) << 8;
            if (bi + 2 < dlen) n |= uint32_t(d[bi + 2]);
            *p++ = b64_table[(n >> 18) & 0x3F];
            *p++ = b64_table[(n >> 12) & 0x3F];
            *p++ = (bi + 1 < dlen) ? b64_table[(n >> 6) & 0x3F] : '=';
            *p++ = (bi + 2 < dlen) ? b64_table[n & 0x3F] : '=';
        }
    }
    emit(s, "</binary>\n"
            "</binaryDataArray>\n");

    // intensity binary — reuse bufs.zlib_out
    bufs.int_f32.resize(frag_count);
    for (uint64_t j = 0; j < frag_count; j++)
        bufs.int_f32[j] = (float)frag_int_data[frag_offset + j];
    zlib_compress_into(bufs, bufs.int_f32.data(), frag_count * sizeof(float));

    size_t int_b64_len = ((bufs.zlib_out.size() + 2) / 3) * 4;
    emit_fmt(s, tmp, sizeof(tmp), "<binaryDataArray encodedLength=\"%zu\">\n", int_b64_len);
    emit(s, "<cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\" value=\"\"/>\n"
            "<cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\" value=\"\"/>\n"
            "<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" value=\"\" unitCvRef=\"MS\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\"/>\n"
            "<binary>");
    b64_start = s.size();
    s.resize(b64_start + int_b64_len);
    {
        const uint8_t* d = bufs.zlib_out.data();
        size_t dlen = bufs.zlib_out.size();
        char* p = s.data() + b64_start;
        for (size_t bi = 0; bi < dlen; bi += 3) {
            uint32_t n = uint32_t(d[bi]) << 16;
            if (bi + 1 < dlen) n |= uint32_t(d[bi + 1]) << 8;
            if (bi + 2 < dlen) n |= uint32_t(d[bi + 2]);
            *p++ = b64_table[(n >> 18) & 0x3F];
            *p++ = b64_table[(n >> 12) & 0x3F];
            *p++ = (bi + 1 < dlen) ? b64_table[(n >> 6) & 0x3F] : '=';
            *p++ = (bi + 2 < dlen) ? b64_table[n & 0x3F] : '=';
        }
    }
    emit(s, "</binary>\n"
            "</binaryDataArray>\n"
            "</binaryDataArrayList>\n"
            "</spectrum>\n");
}

// ─── build header string ────────────────────────────────────────────────────
static std::string build_header(const std::string& run_id, size_t n_spectra) {
    std::string s;
    s.reserve(4096);
    s += "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
    s += "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd\">\n";
    s += "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\" id=\"" + run_id + "\" version=\"1.1.0\">\n";
    s += "<cvList count=\"2\">\n";
    s += "<cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=\"4.1.41\" URI=\"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo\"/>\n";
    s += "<cv id=\"UO\" fullName=\"Unit Ontology\" version=\"09:04:2014\" URI=\"https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo\"/>\n";
    s += "</cvList>\n";
    s += "<fileDescription>\n";
    s += "<fileContent>\n";
    s += "<cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n";
    s += "<cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\" value=\"\"/>\n";
    s += "</fileContent>\n";
    s += "<sourceFileList count=\"1\">\n";
    s += "<sourceFile id=\"" + run_id + ".mgf\" name=\"" + run_id + ".mgf\" location=\"file:///\">\n";
    s += "<cvParam cvRef=\"MS\" accession=\"MS:1000774\" name=\"multiple peak list nativeID format\" value=\"\"/>\n";
    s += "<cvParam cvRef=\"MS\" accession=\"MS:1001062\" name=\"Mascot MGF format\" value=\"\"/>\n";
    s += "</sourceFile>\n";
    s += "</sourceFileList>\n";
    s += "</fileDescription>\n";
    s += "<softwareList count=\"1\">\n";
    s += "<software id=\"mgf_to_mzml_python\" version=\"1.0\">\n";
    s += "<cvParam cvRef=\"MS\" accession=\"MS:1000615\" name=\"ProteoWizard software\" value=\"\"/>\n";
    s += "</software>\n";
    s += "</softwareList>\n";
    s += "<instrumentConfigurationList count=\"1\">\n";
    s += "<instrumentConfiguration id=\"IC\">\n";
    s += "<cvParam cvRef=\"MS\" accession=\"MS:1000031\" name=\"instrument model\" value=\"\"/>\n";
    s += "</instrumentConfiguration>\n";
    s += "</instrumentConfigurationList>\n";
    s += "<dataProcessingList count=\"1\">\n";
    s += "<dataProcessing id=\"python_conversion\">\n";
    s += "<processingMethod order=\"0\" softwareRef=\"mgf_to_mzml_python\">\n";
    s += "<cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\" value=\"\"/>\n";
    s += "</processingMethod>\n";
    s += "</dataProcessing>\n";
    s += "</dataProcessingList>\n";
    s += "<run id=\"" + run_id + "\" defaultInstrumentConfigurationRef=\"IC\">\n";
    char tmp[128];
    snprintf(tmp, sizeof(tmp), "<spectrumList count=\"%zu\" defaultDataProcessingRef=\"python_conversion\">\n", n_spectra);
    s += tmp;
    return s;
}

// ─── build footer string ────────────────────────────────────────────────────
// footer_start_offset = byte position in the file where this footer begins
static std::string build_footer(const std::vector<long>& offsets, long footer_start_offset) {
    std::string s;
    s += "</spectrumList>\n";
    s += "</run>\n";
    s += "</mzML>\n";

    long index_list_offset = footer_start_offset + (long)s.size();

    s += "<indexList count=\"2\">\n";
    s += "<index name=\"spectrum\">\n";
    char tmp[128];
    for (size_t i = 0; i < offsets.size(); i++) {
        snprintf(tmp, sizeof(tmp), "<offset idRef=\"index=%zu\">%ld</offset>\n", i, offsets[i]);
        s += tmp;
    }
    s += "</index>\n";
    s += "<index name=\"chromatogram\">\n";
    s += "</index>\n";
    s += "</indexList>\n";
    snprintf(tmp, sizeof(tmp), "<indexListOffset>%ld</indexListOffset>\n", index_list_offset);
    s += tmp;
    s += "<fileChecksum>0000000000000000000000000000000000000000</fileChecksum>\n";
    s += "</indexedmzML>\n";
    return s;
}

// ─── main ───────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    // Parse CLI args
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <pmsms_dir> <output.mzml> [--precursors-dir DIR] [--run-id NAME] [--zlib-level N] [--threads N]\n", argv[0]);
        return 1;
    }
    std::filesystem::path pmsms_dir = argv[1];
    std::filesystem::path output_path = argv[2];
    std::filesystem::path precursors_dir = pmsms_dir / "filtered_precursors_with_nontrivial_ms2.mmappet";
    std::string run_id = output_path.stem().string();
    int zlib_level = 1;
    int thread_cnt = 1;

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--precursors-dir") == 0 && i + 1 < argc) {
            precursors_dir = argv[++i];
        } else if (strcmp(argv[i], "--run-id") == 0 && i + 1 < argc) {
            run_id = argv[++i];
        } else if (strcmp(argv[i], "--zlib-level") == 0 && i + 1 < argc) {
            zlib_level = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            thread_cnt = atoi(argv[++i]);
            if (thread_cnt < 1) thread_cnt = 1;
        }
    }

    // ─── Step 1: Open datasets ──────────────────────────────────────────
    fprintf(stderr, "Opening datasets...\n");

    Schema<uint32_t, uint32_t, float, float> frag_schema("tof", "intensity", "score", "mz");
    auto frag_ds = frag_schema.open_dataset(pmsms_dir);

    Schema<uint64_t, uint64_t, uint64_t, uint32_t, float> didx_schema("ms1idx", "size", "idx", "max_group_len", "avg_group_len");
    auto didx_ds = didx_schema.open_dataset(pmsms_dir / "dataindex.mmappet");

    Schema<int64_t, int32_t, int32_t, int32_t, double, uint32_t, double, double, uint8_t>
        prec_schema("precursor_id_before_deisotoping", "frame", "scan", "tof",
                    "inv_ion_mobility", "intensity", "mz", "rt", "charge");
    auto prec_ds = prec_schema.open_dataset(precursors_dir);

    size_t n_spectra = prec_ds.size();
    size_t n_dataindex = didx_ds.size();
    if (n_spectra > n_dataindex) {
        fprintf(stderr, "Error: precursor count (%zu) > dataindex count (%zu)\n", n_spectra, n_dataindex);
        return 1;
    }

    auto& frag_mz_col    = frag_ds.get_column<3>();
    auto& frag_int_col    = frag_ds.get_column<1>();
    auto& didx_size_col   = didx_ds.get_column<1>();
    auto& didx_idx_col    = didx_ds.get_column<2>();
    auto& prec_iim_col    = prec_ds.get_column<4>();
    auto& prec_mz_col     = prec_ds.get_column<6>();
    auto& prec_rt_col     = prec_ds.get_column<7>();
    auto& prec_charge_col = prec_ds.get_column<8>();

    fprintf(stderr, "Spectra: %zu, Fragments: %zu, Threads: %d\n", n_spectra, frag_mz_col.size(), thread_cnt);

    // Raw data pointers for render_spectrum
    const float*    frag_mz_data  = frag_mz_col.data();
    const uint32_t* frag_int_data = frag_int_col.data();

    // ─── Build header ───────────────────────────────────────────────────
    std::string header = build_header(run_id, n_spectra);
    size_t header_size = header.size();

    // ─── Progress tracking ──────────────────────────────────────────────
    std::atomic<size_t> progress{0};
    size_t progress_step = n_spectra < 100 ? 1 : n_spectra / 100;

    const char* run_id_cstr = run_id.c_str();

    if (thread_cnt == 1) {
        // ─── Single-threaded path ───────────────────────────────────────
        FILE* out = fopen(output_path.c_str(), "wb");
        if (!out) { fprintf(stderr, "Cannot open output: %s\n", output_path.c_str()); return 1; }

        fwrite(header.data(), 1, header.size(), out);

        std::vector<long> spectrum_offsets(n_spectra);
        RenderBufs bufs;

        for (size_t i = 0; i < n_spectra; i++) {
            if (i % progress_step == 0)
                fprintf(stderr, "\rWriting spectra: %zu / %zu (%zu%%)", i, n_spectra, i * 100 / n_spectra);

            spectrum_offsets[i] = ftell(out);

            render_spectrum(bufs, i, run_id_cstr, zlib_level,
                frag_mz_data, frag_int_data,
                didx_idx_col[i], didx_size_col[i],
                prec_rt_col[i], prec_mz_col[i], prec_charge_col[i], prec_iim_col[i]);

            fwrite(bufs.output.data(), 1, bufs.output.size(), out);
        }
        fprintf(stderr, "\rWriting spectra: %zu / %zu (100%%)\n", n_spectra, n_spectra);

        long footer_start = ftell(out);
        std::string footer = build_footer(spectrum_offsets, footer_start);
        fwrite(footer.data(), 1, footer.size(), out);
        fclose(out);

    } else {
        // ─── Multi-threaded path (pwrite) ───────────────────────────────
        int fd = ::open(output_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (fd < 0) { fprintf(stderr, "Cannot open output: %s\n", output_path.c_str()); return 1; }

        // Write header
        ::pwrite(fd, header.data(), header.size(), 0);

        std::atomic<size_t> write_pos{0};
        std::vector<long> spectrum_offsets(n_spectra);

        auto worker = [&](size_t begin, size_t end) {
            RenderBufs bufs;
            for (size_t i = begin; i < end; i++) {
                render_spectrum(bufs, i, run_id_cstr, zlib_level,
                    frag_mz_data, frag_int_data,
                    didx_idx_col[i], didx_size_col[i],
                    prec_rt_col[i], prec_mz_col[i], prec_charge_col[i], prec_iim_col[i]);

                const std::string& xml = bufs.output;
                size_t pos = write_pos.fetch_add(xml.size());
                spectrum_offsets[i] = (long)(header_size + pos);

                // pwrite may need multiple calls for large buffers
                size_t written = 0;
                while (written < xml.size()) {
                    ssize_t ret = ::pwrite(fd, xml.data() + written, xml.size() - written,
                                           header_size + pos + written);
                    if (ret < 0) {
                        fprintf(stderr, "pwrite failed: %s\n", strerror(errno));
                        exit(1);
                    }
                    written += ret;
                }

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

        // Write footer at the end
        size_t body_size = write_pos.load();
        off_t footer_off = header_size + body_size;
        std::string footer = build_footer(spectrum_offsets, (long)footer_off);

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

        // Truncate file to exact size
        ftruncate(fd, footer_off + footer.size());
        ::close(fd);
    }

    fprintf(stderr, "Done. Wrote %zu spectra to %s\n", n_spectra, output_path.c_str());
    return 0;
}
