CXX      := clang++
CXXFLAGS := -std=c++20 -O3 -Wall -Isrc -Izlib
LDFLAGS  := -Lzlib -lz -lpthread -lcrypto

myprog: main.o MSNumpress.o
	$(CXX) $^ -o $@ $(LDFLAGS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $<

MSNumpress.o: src/MSNumpress.cpp
	$(CXX) $(CXXFLAGS) -c $<

THREADS ?= 1

test: myprog
	mkdir -p output
	./myprog input/testcutF9477.mmappet output/tiny_test.mzml --precursors-dir input/filtered_precursors_with_nontrivial_ms2_44.mmappet --run-id tiny --threads $(THREADS)

run: myprog
	mkdir -p output
	./myprog input/F9477.mmappet output/F9477.mzml --threads $(THREADS)

clean:
	rm -f *.o myprog