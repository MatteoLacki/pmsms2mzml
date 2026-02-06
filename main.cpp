#include "mmappet.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <filesystem>
#include <zlib.h>

// ─── base64 encode ──────────────────────────────────────────────────────────
static const char b64_table[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static std::string base64_encode(const uint8_t* data, size_t len) {
    std::string out;
    out.reserve(((len + 2) / 3) * 4);
    for (size_t i = 0; i < len; i += 3) {
        uint32_t n = uint32_t(data[i]) << 16;
        if (i + 1 < len) n |= uint32_t(data[i + 1]) << 8;
        if (i + 2 < len) n |= uint32_t(data[i + 2]);
        out.push_back(b64_table[(n >> 18) & 0x3F]);
        out.push_back(b64_table[(n >> 12) & 0x3F]);
        out.push_back((i + 1 < len) ? b64_table[(n >> 6) & 0x3F] : '=');
        out.push_back((i + 2 < len) ? b64_table[n & 0x3F] : '=');
    }
    return out;
}

// ─── zlib compress ──────────────────────────────────────────────────────────
static std::vector<uint8_t> zlib_compress(const void* data, size_t len, int level) {
    uLongf bound = compressBound(len);
    std::vector<uint8_t> out(bound);
    int ret = compress2(out.data(), &bound, (const Bytef*)data, len, level);
    if (ret != Z_OK) {
        fprintf(stderr, "zlib compress2 failed: %d\n", ret);
        exit(1);
    }
    out.resize(bound);
    return out;
}

// ─── number formatting ──────────────────────────────────────────────────────
// Format a floating-point value with enough decimal places, stripping trailing
// zeros but keeping at least one decimal place (e.g., 202839.0, 576.285).
static std::string fmt_float(double val, int max_decimals = 6) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%.*f", max_decimals, val);
    // Strip trailing zeros, but keep at least one after the dot
    char* dot = strchr(buf, '.');
    if (dot) {
        char* end = buf + strlen(buf) - 1;
        while (end > dot + 1 && *end == '0') --end;
        *(end + 1) = '\0';
    }
    return buf;
}

// ─── main ───────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    // Parse CLI args
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <pmsms_dir> <precursors_dir> <output.mzml> [--run-id NAME] [--zlib-level N]\n", argv[0]);
        return 1;
    }
    std::filesystem::path pmsms_dir = argv[1];
    std::filesystem::path precursors_dir = argv[2];
    std::filesystem::path output_path = argv[3];
    std::string run_id = output_path.stem().string();
    int zlib_level = 1;

    for (int i = 4; i < argc; i++) {
        if (strcmp(argv[i], "--run-id") == 0 && i + 1 < argc) {
            run_id = argv[++i];
        } else if (strcmp(argv[i], "--zlib-level") == 0 && i + 1 < argc) {
            zlib_level = atoi(argv[++i]);
        }
    }

    // ─── Step 1: Open datasets ──────────────────────────────────────────
    fprintf(stderr, "Opening datasets...\n");

    // Fragment data: tof(u32), intensity(u32), score(f32), mz(f32)
    Schema<uint32_t, uint32_t, float, float> frag_schema("tof", "intensity", "score", "mz");
    auto frag_ds = frag_schema.open_dataset(pmsms_dir);

    // Dataindex: ms1idx(u64), size(u64), idx(u64), max_group_len(u32), avg_group_len(f32)
    Schema<uint64_t, uint64_t, uint64_t, uint32_t, float> didx_schema("ms1idx", "size", "idx", "max_group_len", "avg_group_len");
    auto didx_ds = didx_schema.open_dataset(pmsms_dir / "dataindex.mmappet");

    // Precursor metadata: precursor_id_before_deisotoping(i64), frame(i32), scan(i32),
    //   tof(i32), inv_ion_mobility(f64), intensity(u32), mz(f64), rt(f64), charge(u8)
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

    // Get column references
    auto& frag_mz_col   = frag_ds.get_column<3>();  // mz (float)
    auto& frag_int_col   = frag_ds.get_column<1>();  // intensity (uint32_t)
    auto& didx_size_col  = didx_ds.get_column<1>();  // size (uint64_t)
    auto& didx_idx_col   = didx_ds.get_column<2>();  // idx (uint64_t)
    auto& prec_id_col    = prec_ds.get_column<0>();  // precursor_id_before_deisotoping (int64_t)
    auto& prec_iim_col   = prec_ds.get_column<4>();  // inv_ion_mobility (double)
    auto& prec_mz_col    = prec_ds.get_column<6>();  // mz (double)
    auto& prec_rt_col    = prec_ds.get_column<7>();  // rt (double)
    auto& prec_charge_col= prec_ds.get_column<8>();  // charge (uint8_t)

    fprintf(stderr, "Spectra: %zu, Fragments: %zu\n", n_spectra, frag_mz_col.size());

    // ─── Step 2: Open output and write mzML header ──────────────────────
    FILE* out = fopen(output_path.c_str(), "wb");
    if (!out) {
        fprintf(stderr, "Cannot open output: %s\n", output_path.c_str());
        return 1;
    }

    // Preallocate offset vector
    std::vector<long> spectrum_offsets;
    spectrum_offsets.reserve(n_spectra);

    // XML preamble
    fprintf(out, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
    fprintf(out, "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd\">\n");
    fprintf(out, "  <mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\" id=\"%s\" version=\"1.1.0\">\n", run_id.c_str());
    fprintf(out, "    <cvList count=\"2\">\n");
    fprintf(out, "      <cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=\"4.1.41\" URI=\"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo\"/>\n");
    fprintf(out, "      <cv id=\"UO\" fullName=\"Unit Ontology\" version=\"09:04:2014\" URI=\"https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo\"/>\n");
    fprintf(out, "    </cvList>\n");
    fprintf(out, "    <fileDescription>\n");
    fprintf(out, "      <fileContent>\n");
    fprintf(out, "        <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n");
    fprintf(out, "        <cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\" value=\"\"/>\n");
    fprintf(out, "      </fileContent>\n");
    fprintf(out, "      <sourceFileList count=\"1\">\n");
    fprintf(out, "        <sourceFile id=\"%s.mgf\" name=\"%s.mgf\" location=\"file:///\">\n", run_id.c_str(), run_id.c_str());
    fprintf(out, "          <cvParam cvRef=\"MS\" accession=\"MS:1000774\" name=\"multiple peak list nativeID format\" value=\"\"/>\n");
    fprintf(out, "          <cvParam cvRef=\"MS\" accession=\"MS:1001062\" name=\"Mascot MGF format\" value=\"\"/>\n");
    fprintf(out, "        </sourceFile>\n");
    fprintf(out, "      </sourceFileList>\n");
    fprintf(out, "    </fileDescription>\n");
    fprintf(out, "    <softwareList count=\"1\">\n");
    fprintf(out, "      <software id=\"mgf_to_mzml_python\" version=\"1.0\">\n");
    fprintf(out, "        <cvParam cvRef=\"MS\" accession=\"MS:1000615\" name=\"ProteoWizard software\" value=\"\"/>\n");
    fprintf(out, "      </software>\n");
    fprintf(out, "    </softwareList>\n");
    fprintf(out, "    <instrumentConfigurationList count=\"1\">\n");
    fprintf(out, "      <instrumentConfiguration id=\"IC\">\n");
    fprintf(out, "        <cvParam cvRef=\"MS\" accession=\"MS:1000031\" name=\"instrument model\" value=\"\"/>\n");
    fprintf(out, "      </instrumentConfiguration>\n");
    fprintf(out, "    </instrumentConfigurationList>\n");
    fprintf(out, "    <dataProcessingList count=\"1\">\n");
    fprintf(out, "      <dataProcessing id=\"python_conversion\">\n");
    fprintf(out, "        <processingMethod order=\"0\" softwareRef=\"mgf_to_mzml_python\">\n");
    fprintf(out, "          <cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\" value=\"\"/>\n");
    fprintf(out, "        </processingMethod>\n");
    fprintf(out, "      </dataProcessing>\n");
    fprintf(out, "    </dataProcessingList>\n");
    fprintf(out, "    <run id=\"%s\" defaultInstrumentConfigurationRef=\"IC\">\n", run_id.c_str());
    fprintf(out, "      <spectrumList count=\"%zu\" defaultDataProcessingRef=\"python_conversion\">\n", n_spectra);

    // ─── Step 3: Write spectra ──────────────────────────────────────────
    fprintf(stderr, "Writing %zu spectra...\n", n_spectra);

    // Temp buffers for float conversion
    std::vector<float> intensity_f32;

    for (size_t i = 0; i < n_spectra; i++) {
        // Record byte offset
        spectrum_offsets.push_back(ftell(out));

        // Read precursor metadata
        double rt          = prec_rt_col[i];
        double prec_mz     = prec_mz_col[i];
        uint8_t charge     = prec_charge_col[i];
        double iim         = prec_iim_col[i];
        int64_t prec_id    = prec_id_col[i];
        (void)prec_id;  // reserved for future use in title

        // Read fragment slice via dataindex
        uint64_t frag_offset = didx_idx_col[i];
        uint64_t frag_count  = didx_size_col[i];

        // Compute per-spectrum stats
        float lowest_mz = 0, highest_mz = 0;
        double tic = 0;
        float base_peak_mz = 0;
        float base_peak_int = 0;

        for (uint64_t j = 0; j < frag_count; j++) {
            float mz = frag_mz_col[frag_offset + j];
            uint32_t intensity = frag_int_col[frag_offset + j];
            float fi = (float)intensity;
            if (j == 0 || mz < lowest_mz) lowest_mz = mz;
            if (j == 0 || mz > highest_mz) highest_mz = mz;
            tic += fi;
            if (j == 0 || fi > base_peak_int) {
                base_peak_int = fi;
                base_peak_mz = mz;
            }
        }

        // Encode mz array: float32 bytes → zlib → base64
        auto mz_compressed = zlib_compress(
            &frag_mz_col[frag_offset], frag_count * sizeof(float), zlib_level);
        auto mz_b64 = base64_encode(mz_compressed.data(), mz_compressed.size());

        // Encode intensity array: uint32 → float32 → bytes → zlib → base64
        intensity_f32.resize(frag_count);
        for (uint64_t j = 0; j < frag_count; j++)
            intensity_f32[j] = (float)frag_int_col[frag_offset + j];
        auto int_compressed = zlib_compress(
            intensity_f32.data(), frag_count * sizeof(float), zlib_level);
        auto int_b64 = base64_encode(int_compressed.data(), int_compressed.size());

        // Write spectrum XML
        fprintf(out, "        <spectrum index=\"%zu\" id=\"index=%zu\" defaultArrayLength=\"%zu\">\n",
                i, i, (size_t)frag_count);
        fprintf(out, "          <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n");
        fprintf(out, "          <cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"2\"/>\n");
        fprintf(out, "          <cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\" value=\"\"/>\n");

        // Spectrum title: run_id.idx.idx.2 File:"", NativeID:"index=idx"
        fprintf(out, "          <cvParam cvRef=\"MS\" accession=\"MS:1000796\" name=\"spectrum title\" value=\"%s.%zu.%zu.2 File:&quot;&quot;, NativeID:&quot;index=%zu&quot;\"/>\n",
                run_id.c_str(), i, i, i);

        fprintf(out, "          <cvParam cvRef=\"MS\" accession=\"MS:1000130\" name=\"positive scan\" value=\"\"/>\n");
        fprintf(out, "          <cvParam cvRef=\"MS\" accession=\"MS:1000528\" name=\"lowest observed m/z\" value=\"%s\"/>\n", fmt_float(lowest_mz, 3).c_str());
        fprintf(out, "          <cvParam cvRef=\"MS\" accession=\"MS:1000527\" name=\"highest observed m/z\" value=\"%s\"/>\n", fmt_float(highest_mz, 3).c_str());
        fprintf(out, "          <cvParam cvRef=\"MS\" accession=\"MS:1000285\" name=\"total ion current\" value=\"%s\"/>\n", fmt_float(tic, 1).c_str());
        fprintf(out, "          <cvParam cvRef=\"MS\" accession=\"MS:1000504\" name=\"base peak m/z\" value=\"%s\"/>\n", fmt_float(base_peak_mz, 3).c_str());
        fprintf(out, "          <cvParam cvRef=\"MS\" accession=\"MS:1000505\" name=\"base peak intensity\" value=\"%s\"/>\n", fmt_float(base_peak_int, 1).c_str());

        fprintf(out, "          <scanList count=\"1\">\n");
        fprintf(out, "            <cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\" value=\"\"/>\n");
        fprintf(out, "            <scan>\n");
        fprintf(out, "              <cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"%s\" unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"/>\n", fmt_float(rt).c_str());
        fprintf(out, "              <cvParam cvRef=\"MS\" accession=\"MS:1002815\" name=\"inverse reduced ion mobility\" value=\"%s\" unitCvRef=\"MS\" unitAccession=\"MS:1002814\" unitName=\"volt-second per square centimeter\"/>\n", fmt_float(iim, 4).c_str());
        fprintf(out, "            </scan>\n");
        fprintf(out, "          </scanList>\n");

        fprintf(out, "          <precursorList count=\"1\">\n");
        fprintf(out, "            <precursor>\n");
        fprintf(out, "              <selectedIonList count=\"1\">\n");
        fprintf(out, "                <selectedIon>\n");
        fprintf(out, "                  <cvParam cvRef=\"MS\" accession=\"MS:1000744\" name=\"selected ion m/z\" value=\"%s\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n", fmt_float(prec_mz).c_str());
        fprintf(out, "                  <cvParam cvRef=\"MS\" accession=\"MS:1000041\" name=\"charge state\" value=\"%u\"/>\n", (unsigned)charge);
        fprintf(out, "                </selectedIon>\n");
        fprintf(out, "              </selectedIonList>\n");
        fprintf(out, "              <activation>\n");
        fprintf(out, "              </activation>\n");
        fprintf(out, "            </precursor>\n");
        fprintf(out, "          </precursorList>\n");

        fprintf(out, "          <binaryDataArrayList count=\"2\">\n");

        // m/z array
        fprintf(out, "            <binaryDataArray encodedLength=\"%zu\">\n", mz_b64.size());
        fprintf(out, "              <cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\" value=\"\"/>\n");
        fprintf(out, "              <cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\" value=\"\"/>\n");
        fprintf(out, "              <cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\" value=\"\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n");
        fprintf(out, "              <binary>%s</binary>\n", mz_b64.c_str());
        fprintf(out, "            </binaryDataArray>\n");

        // intensity array
        fprintf(out, "            <binaryDataArray encodedLength=\"%zu\">\n", int_b64.size());
        fprintf(out, "              <cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\" value=\"\"/>\n");
        fprintf(out, "              <cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\" value=\"\"/>\n");
        fprintf(out, "              <cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" value=\"\" unitCvRef=\"MS\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\"/>\n");
        fprintf(out, "              <binary>%s</binary>\n", int_b64.c_str());
        fprintf(out, "            </binaryDataArray>\n");

        fprintf(out, "          </binaryDataArrayList>\n");
        fprintf(out, "        </spectrum>\n");
    }

    // ─── Step 4: Write footer ───────────────────────────────────────────
    fprintf(out, "      </spectrumList>\n");
    fprintf(out, "    </run>\n");
    fprintf(out, "  </mzML>\n");

    long index_list_offset = ftell(out);

    fprintf(out, "  <indexList count=\"2\">\n");
    fprintf(out, "    <index name=\"spectrum\">\n");
    for (size_t i = 0; i < n_spectra; i++) {
        fprintf(out, "      <offset idRef=\"index=%zu\">%ld</offset>\n", i, spectrum_offsets[i]);
    }
    fprintf(out, "    </index>\n");
    fprintf(out, "    <index name=\"chromatogram\">\n");
    fprintf(out, "    </index>\n");
    fprintf(out, "  </indexList>\n");
    fprintf(out, "  <indexListOffset>%ld</indexListOffset>\n", index_list_offset);
    fprintf(out, "  <fileChecksum>0000000000000000000000000000000000000000</fileChecksum>\n");
    fprintf(out, "</indexedmzML>\n");

    fclose(out);
    fprintf(stderr, "Done. Wrote %zu spectra to %s\n", n_spectra, output_path.c_str());
    return 0;
}
