CXX      := clang++
# -ffp-contract=off: fragment m/z read through --tof2mz must round exactly as SAGE's
# reader does (no fused multiply-add), so exports hold the m/z SAGE searched.
CXXFLAGS := -std=c++23 -O3 -Wall -Isrc -ffp-contract=off
LDFLAGS  := -lz

pmsms2mzml: main.o MSNumpress.o
	$(CXX) $^ -o $@ $(LDFLAGS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $<

MSNumpress.o: src/MSNumpress.cpp
	$(CXX) $(CXXFLAGS) -c $<

THREADS ?= 1

test: pmsms2mzml
	mkdir -p output
	./pmsms2mzml input/testcutF9477.mmappet input/filtered_precursors_with_nontrivial_ms2_44.mmappet output/tiny_test.mzml --run-id tiny --threads $(THREADS)

run: pmsms2mzml
	mkdir -p output
	./pmsms2mzml input/F9477.mmappet input/F9477.mmappet/filtered_precursors_with_nontrivial_ms2.mmappet output/F9477.mzml --threads $(THREADS)

tests/test_pmsms2mzml: tests/test_pmsms2mzml.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

test-unit: pmsms2mzml tests/test_pmsms2mzml
	./tests/test_pmsms2mzml

clean:
	rm -f *.o pmsms2mzml tests/test_pmsms2mzml