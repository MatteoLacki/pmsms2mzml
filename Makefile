CXX      := clang++
CXXFLAGS := -std=c++20 -O3 -Wall -Isrc -Izlib
LDFLAGS  := -Lzlib -lz

myprog: main.o
	$(CXX) $^ -o $@ $(LDFLAGS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $<

test: myprog
	mkdir -p output
	./myprog input/testcutF9477.mmappet input/filtered_precursors_with_nontrivial_ms2_44.mmappet output/tiny_test.mzml --run-id tiny

clean:
	rm -f *.o myprog