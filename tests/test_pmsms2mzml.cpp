#include "mmappet.h"
#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <cstdio>
#include <string>
#include <vector>
#include <cstring>
#include <zlib.h>

namespace fs = std::filesystem;

// ── input data: charge test (single-charge schema) ───────────────────────────

static constexpr int N_FRAG = 10;
static constexpr int N_PREC = 3;

static uint32_t frag_tof[N_FRAG]       = {0,1,2,3,4,5,6,7,8,9};
static uint32_t frag_intensity[N_FRAG] = {100,200,300,400,500,600,700,800,900,1000};
static float    frag_score[N_FRAG]     = {0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f};
// All mz values are exact binary fractions (powers-of-2 denominators)
static float    frag_mz[N_FRAG]        = {100.125f,200.25f,300.5f,
                                           150.125f,250.25f,350.5f,
                                           120.125f,220.25f,320.5f,420.0f};

static int64_t  prec_pid[N_PREC]    = {1, 2, 3};
static int32_t  prec_frame[N_PREC]  = {10, 20, 30};
static int32_t  prec_scan[N_PREC]   = {1, 2, 3};
static int32_t  prec_tof[N_PREC]    = {100, 200, 300};
static double   prec_iim[N_PREC]    = {0.5, 0.6, 0.7};
static uint32_t prec_inten[N_PREC]  = {1000, 2000, 3000};
static double   prec_mz[N_PREC]     = {500.123, 600.456, 700.789};
static double   prec_rt[N_PREC]     = {60.0, 120.0, 180.0};
static uint8_t  prec_charge[N_PREC] = {2, 3, 2};
static uint64_t prec_fstart[N_PREC] = {0, 3, 6};
static uint64_t prec_fcnt[N_PREC]   = {3, 3, 4};

// ── print input tables ───────────────────────────────────────────────────────

static void print_fragments() {
    printf("\n=== Fragment dataset (pmsms.mmappet, %d rows) ===\n", N_FRAG);
    printf("  row │   tof │ intensity │ score │      mz\n");
    printf("──────┼───────┼───────────┼───────┼─────────\n");
    for (int i = 0; i < N_FRAG; ++i)
        printf("  %3d │   %3u │      %4u │   %.1f │ %7.3f\n",
               i, frag_tof[i], frag_intensity[i], frag_score[i], frag_mz[i]);
}

static void print_precursors() {
    printf("\n=== Precursor dataset (precursors.mmappet, %d rows) ===\n", N_PREC);
    printf("  row │ pid │ frame │ scan │  tof │   iim │ intensity │       mz │     rt │ charge │ fstart │ fcnt\n");
    printf("──────┼─────┼───────┼──────┼──────┼───────┼───────────┼──────────┼────────┼────────┼────────┼─────\n");
    for (int i = 0; i < N_PREC; ++i)
        printf("  %3d │ %3lld │   %3d │  %3d │  %3d │   %.1f │      %4u │ %8.3f │ %6.1f │      %u │      %llu │    %llu\n",
               i,
               (long long)prec_pid[i],
               prec_frame[i], prec_scan[i], prec_tof[i],
               prec_iim[i], prec_inten[i], prec_mz[i], prec_rt[i],
               (unsigned)prec_charge[i],
               (unsigned long long)prec_fstart[i],
               (unsigned long long)prec_fcnt[i]);
}

// ── dataset writers ──────────────────────────────────────────────────────────

static void write_fragments(const fs::path& dir) {
    Schema<uint32_t, uint32_t, float, float> s("tof", "intensity", "score", "mz");
    auto w = s.create_writer(dir);
    w.write_rows(N_FRAG, frag_tof, frag_intensity, frag_score, frag_mz);
}

static void write_precursors(const fs::path& dir) {
    Schema<int64_t, int32_t, int32_t, int32_t, double,
           uint32_t, double, double, uint8_t, uint64_t, uint64_t>
        s("precursor_id_before_deisotoping", "frame", "scan", "tof",
          "inv_ion_mobility", "intensity", "mz", "rt", "charge",
          "fragment_spectrum_start", "fragment_event_cnt");
    auto w = s.create_writer(dir);
    w.write_rows(N_PREC, prec_pid, prec_frame, prec_scan, prec_tof,
                 prec_iim, prec_inten, prec_mz, prec_rt, prec_charge,
                 prec_fstart, prec_fcnt);
}

// ── base64 + zlib helpers ────────────────────────────────────────────────────

static std::vector<uint8_t> base64_decode(const std::string& in) {
    static const char* table =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::vector<uint8_t> out;
    int val = 0, bits = -8;
    for (unsigned char c : in) {
        if (c == '=' || c == '\n' || c == '\r') continue;
        const char* p = strchr(table, c);
        if (!p) continue;
        val = (val << 6) + (int)(p - table);
        bits += 6;
        if (bits >= 0) {
            out.push_back((val >> bits) & 0xFF);
            bits -= 8;
        }
    }
    return out;
}

static std::vector<float> zlib_decompress_floats(const std::vector<uint8_t>& compressed) {
    std::vector<uint8_t> buf(compressed.size() * 8);
    uLongf out_size = buf.size();
    while (true) {
        int ret = uncompress(buf.data(), &out_size, compressed.data(), compressed.size());
        if (ret == Z_OK) break;
        if (ret == Z_BUF_ERROR) { buf.resize(buf.size() * 2); out_size = buf.size(); continue; }
        throw std::runtime_error("zlib decompression failed: " + std::to_string(ret));
    }
    size_t n = out_size / sizeof(float);
    std::vector<float> out(n);
    memcpy(out.data(), buf.data(), n * sizeof(float));
    return out;
}

// ── mzML binary extraction ───────────────────────────────────────────────────

struct SpectrumArrays { std::vector<float> mz, intensity; };

static std::string between(const std::string& s, const std::string& open,
                           const std::string& close, size_t from = 0) {
    size_t a = s.find(open, from);
    if (a == std::string::npos) return {};
    a += open.size();
    size_t b = s.find(close, a);
    if (b == std::string::npos) return {};
    return s.substr(a, b - a);
}

static std::vector<SpectrumArrays> extract_arrays(const std::string& xml) {
    std::vector<SpectrumArrays> result;
    size_t pos = 0;
    while (true) {
        size_t s0 = xml.find("<spectrum ", pos);
        if (s0 == std::string::npos) break;
        size_t s1 = xml.find("</spectrum>", s0);
        if (s1 == std::string::npos) break;
        std::string spec = xml.substr(s0, s1 - s0);

        SpectrumArrays sa;
        size_t apos = 0;
        while (true) {
            size_t a0 = spec.find("<binaryDataArray ", apos);
            if (a0 == std::string::npos) break;
            size_t a1 = spec.find("</binaryDataArray>", a0);
            if (a1 == std::string::npos) break;
            std::string blk = spec.substr(a0, a1 - a0);

            std::string b64 = between(blk, "<binary>", "</binary>");
            auto floats = zlib_decompress_floats(base64_decode(b64));

            if (blk.find("name=\"m/z array\"") != std::string::npos)
                sa.mz = floats;
            else if (blk.find("name=\"intensity array\"") != std::string::npos)
                sa.intensity = floats;

            apos = a1 + 18;
        }
        result.push_back(sa);
        pos = s1 + 11;
    }
    return result;
}

// ── binary verification table ─────────────────────────────────────────────────

static void verify_and_print_arrays(const std::vector<SpectrumArrays>& spectra) {
    printf("\n=== Binary array verification ===\n");
    bool all_ok = true;

    for (int s = 0; s < N_PREC; ++s) {
        int n      = (int)prec_fcnt[s];
        int fstart = (int)prec_fstart[s];
        printf("\nSpectrum %d  (prec_mz=%.3f  charge=%u  rt=%.1f  fragments=%d):\n",
               s, prec_mz[s], (unsigned)prec_charge[s], prec_rt[s], n);
        printf("  idx │  mz_expected │   mz_decoded │ int_expected │ int_decoded\n");
        printf("──────┼──────────────┼──────────────┼──────────────┼────────────\n");

        bool spec_ok = true;
        for (int f = 0; f < n; ++f) {
            int fi        = fstart + f;
            float exp_mz  = frag_mz[fi];
            float exp_int = (float)frag_intensity[fi];
            float dec_mz  = (f < (int)spectra[s].mz.size())        ? spectra[s].mz[f]        : -1.f;
            float dec_int = (f < (int)spectra[s].intensity.size())  ? spectra[s].intensity[f] : -1.f;
            bool row_ok   = (exp_mz == dec_mz) && (exp_int == dec_int);
            if (!row_ok) spec_ok = false;
            printf("  %3d │  %12.3f │  %12.3f │  %12.0f │  %10.0f  %s\n",
                   f, exp_mz, dec_mz, exp_int, dec_int, row_ok ? "" : "<-- MISMATCH");
        }
        printf("  %s\n", spec_ok ? "OK" : "FAILED");
        if (!spec_ok) all_ok = false;
    }
    assert(all_ok && "binary array mismatch");
}

// ── input data: charges test (multi-charge schema, int64 charges column) ─────

static constexpr int MC_N_FRAG = 6;
static constexpr int MC_N_PREC = 2;

static uint32_t mc_frag_tof[MC_N_FRAG]       = {0,1,2,3,4,5};
static uint32_t mc_frag_intensity[MC_N_FRAG] = {100,200,300,400,500,600};
static float    mc_frag_score[MC_N_FRAG]     = {0.5f,0.5f,0.5f,0.5f,0.5f,0.5f};
static float    mc_frag_mz[MC_N_FRAG]        = {100.125f,200.25f,300.5f,
                                                 150.125f,250.25f,350.5f};

static int64_t  mc_prec_pid[MC_N_PREC]     = {1, 2};
static int32_t  mc_prec_frame[MC_N_PREC]   = {10, 20};
static int32_t  mc_prec_scan[MC_N_PREC]    = {1, 2};
static int32_t  mc_prec_tof[MC_N_PREC]     = {100, 200};
static double   mc_prec_iim[MC_N_PREC]     = {0.5, 0.6};
static uint32_t mc_prec_inten[MC_N_PREC]   = {1000, 2000};
static double   mc_prec_mz[MC_N_PREC]      = {500.123, 600.456};
static double   mc_prec_rt[MC_N_PREC]      = {60.0, 120.0};
static int64_t  mc_prec_charges[MC_N_PREC] = {23, 2};  // prec0→[2,3], prec1→[2]
static uint64_t mc_prec_fstart[MC_N_PREC]  = {0, 3};
static uint64_t mc_prec_fcnt[MC_N_PREC]    = {3, 3};

// ── file helper ───────────────────────────────────────────────────────────────

static std::string read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("Cannot open: " + p.string());
    return std::string(std::istreambuf_iterator<char>(f), {});
}

static bool has(const std::string& s, const std::string& sub) {
    return s.find(sub) != std::string::npos;
}

// ── charges test ──────────────────────────────────────────────────────────────

static void run_charges_test() {
    printf("\n\n==============================\n");
    printf("=== charges test (multi-charge expansion) ===\n");
    printf("==============================\n");

    printf("\n=== Fragment dataset (mc_pmsms.mmappet, %d rows) ===\n", MC_N_FRAG);
    printf("  row │      mz │ intensity\n");
    printf("──────┼─────────┼──────────\n");
    for (int i = 0; i < MC_N_FRAG; ++i)
        printf("  %3d │ %7.3f │      %4u\n", i, mc_frag_mz[i], mc_frag_intensity[i]);

    printf("\n=== Precursor dataset (mc_precursors.mmappet, %d rows, int64 charges) ===\n", MC_N_PREC);
    printf("  row │ pid │       mz │     rt │ charges │ fstart │ fcnt\n");
    printf("──────┼─────┼──────────┼────────┼─────────┼────────┼─────\n");
    for (int i = 0; i < MC_N_PREC; ++i)
        printf("  %3d │ %3lld │ %8.3f │ %6.1f │     %3lld │      %llu │    %llu\n",
               i,
               (long long)mc_prec_pid[i],
               mc_prec_mz[i], mc_prec_rt[i],
               (long long)mc_prec_charges[i],
               (unsigned long long)mc_prec_fstart[i],
               (unsigned long long)mc_prec_fcnt[i]);

    fs::path tmp = fs::temp_directory_path() / "pmsms2mzml_charges_test";
    fs::remove_all(tmp);

    // Write fragment dataset
    {
        Schema<uint32_t, uint32_t, float, float> s("tof", "intensity", "score", "mz");
        auto w = s.create_writer(tmp / "mc_pmsms.mmappet");
        w.write_rows(MC_N_FRAG, mc_frag_tof, mc_frag_intensity, mc_frag_score, mc_frag_mz);
    }

    // Write precursor dataset with int64 charges column
    {
        Schema<int64_t,int32_t,int32_t,int32_t,double,
               uint32_t,double,double,int64_t,uint64_t,uint64_t>
            s("precursor_id_before_deisotoping","frame","scan","tof",
              "inv_ion_mobility","intensity","mz","rt","charges",
              "fragment_spectrum_start","fragment_event_cnt");
        auto w = s.create_writer(tmp / "mc_precursors.mmappet");
        w.write_rows(MC_N_PREC, mc_prec_pid, mc_prec_frame, mc_prec_scan, mc_prec_tof,
                     mc_prec_iim, mc_prec_inten, mc_prec_mz, mc_prec_rt, mc_prec_charges,
                     mc_prec_fstart, mc_prec_fcnt);
    }

    fs::path out = tmp / "mc_output.mzML";
    std::string cmd =
        "./pmsms2mzml " + (tmp / "mc_pmsms.mmappet").string() +
        " " + (tmp / "mc_precursors.mmappet").string() +
        " " + out.string() +
        " --run-id mc_test --threads 1 --indexed";

    printf("\n=== Invoking converter ===\n%s\n", cmd.c_str());

    if (system(cmd.c_str()) != 0) {
        fprintf(stderr, "pmsms2mzml exited with non-zero status\n");
        exit(1);
    }

    std::string xml = read_file(out);

    // ── text assertions ───────────────────────────────────────────────────────
    assert(has(xml, "spectrumList count=\"3\""));
    assert(has(xml, "name=\"charge state\" value=\"2\""));
    assert(has(xml, "name=\"charge state\" value=\"3\""));
    assert(has(xml, "name=\"selected ion m/z\" value=\"500.123\""));
    assert(has(xml, "name=\"selected ion m/z\" value=\"600.456\""));
    assert(has(xml, "name=\"scan start time\" value=\"60.0\""));
    assert(has(xml, "name=\"scan start time\" value=\"120.0\""));

    // ── binary array assertions ───────────────────────────────────────────────
    auto spectra = extract_arrays(xml);
    assert(spectra.size() == 3);

    static const float exp_mz0[3]  = {100.125f, 200.25f, 300.5f};
    static const float exp_int0[3] = {100.f, 200.f, 300.f};
    static const float exp_mz2[3]  = {150.125f, 250.25f, 350.5f};
    static const float exp_int2[3] = {400.f, 500.f, 600.f};

    printf("\n=== Binary array verification (charges test) ===\n");
    bool all_ok = true;

    const char* spec_labels[3] = {
        "Spectrum 0 (prec 0, charge=2, mz=500.123, rt=60.0)",
        "Spectrum 1 (prec 0, charge=3, mz=500.123, rt=60.0) [same fragment slice]",
        "Spectrum 2 (prec 1, charge=2, mz=600.456, rt=120.0)"
    };
    const float* exp_mz_arr[3]  = {exp_mz0,  exp_mz0,  exp_mz2};
    const float* exp_int_arr[3] = {exp_int0, exp_int0, exp_int2};

    for (int s = 0; s < 3; ++s) {
        printf("\n%s:\n", spec_labels[s]);
        printf("  idx │  mz_expected │  mz_decoded │  int_expected │  int_decoded\n");
        printf("──────┼──────────────┼─────────────┼───────────────┼─────────────\n");
        bool spec_ok = true;
        for (int f = 0; f < 3; ++f) {
            float dm = (f < (int)spectra[s].mz.size())        ? spectra[s].mz[f]        : -1.f;
            float di = (f < (int)spectra[s].intensity.size()) ? spectra[s].intensity[f] : -1.f;
            bool ok = (dm == exp_mz_arr[s][f]) && (di == exp_int_arr[s][f]);
            if (!ok) spec_ok = false;
            printf("  %3d │  %12.3f │  %11.3f │  %13.0f │  %11.0f  %s\n",
                   f, exp_mz_arr[s][f], dm, exp_int_arr[s][f], di, ok ? "" : "<-- MISMATCH");
        }
        printf("  %s\n", spec_ok ? "OK" : "FAILED");
        if (!spec_ok) all_ok = false;
    }

    assert(all_ok && "charges binary array mismatch");

    printf("\ncharges test passed.\n");
    fs::remove_all(tmp);
}

// ── main ──────────────────────────────────────────────────────────────────────

int main() {
    print_fragments();
    print_precursors();

    fs::path tmp = fs::temp_directory_path() / "pmsms2mzml_test";
    fs::remove_all(tmp);

    write_fragments(tmp / "pmsms.mmappet");
    write_precursors(tmp / "precursors.mmappet");

    fs::path out = tmp / "output.mzML";
    std::string cmd =
        "./pmsms2mzml " + (tmp / "pmsms.mmappet").string() +
        " " + (tmp / "precursors.mmappet").string() +
        " " + out.string() +
        " --run-id test --threads 1 --indexed";

    printf("\n=== Invoking converter ===\n%s\n", cmd.c_str());

    if (system(cmd.c_str()) != 0) {
        fprintf(stderr, "pmsms2mzml exited with non-zero status\n");
        return 1;
    }

    std::string xml = read_file(out);

    // ── text assertions ───────────────────────────────────────────────────────
    assert(has(xml, "spectrumList count=\"3\""));
    assert(has(xml, "defaultArrayLength=\"3\""));
    assert(has(xml, "defaultArrayLength=\"4\""));
    assert(has(xml, "name=\"charge state\" value=\"2\""));
    assert(has(xml, "name=\"charge state\" value=\"3\""));
    assert(has(xml, "name=\"selected ion m/z\" value=\"500.123\""));
    assert(has(xml, "name=\"selected ion m/z\" value=\"600.456\""));
    assert(has(xml, "name=\"selected ion m/z\" value=\"700.789\""));
    assert(has(xml, "name=\"scan start time\" value=\"60.0\""));
    assert(has(xml, "name=\"scan start time\" value=\"120.0\""));
    assert(has(xml, "name=\"scan start time\" value=\"180.0\""));

    // ── binary array assertions ───────────────────────────────────────────────
    verify_and_print_arrays(extract_arrays(xml));

    printf("\ncharge test passed.\n");
    fs::remove_all(tmp);

    run_charges_test();
    return 0;
}
