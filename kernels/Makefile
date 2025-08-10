CXX      := g++
CXXFLAGS := -O3 -std=c++17 -g
LDFLAGS  := -lz

# Single
build_word_table: build_word_table.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# OpenMP
build_word_table_mt: build_word_table_mt.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $< -o $@ $(LDFLAGS)

# OpenMP Search DOWN
search_word_table_mt: search_word_table_mt.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $< -o $@ $(LDFLAGS)

# OpenMP Search UP
search_word_table_mt_up: search_word_table_mt_up.cpp
	$(CXX) $(CXXFLAGS) -fopenmp $< -o $@ $(LDFLAGS)

all: build_word_table build_word_table_mt search_word_table_mt  search_word_table_mt_up

clean:
	rm -f build_word_table build_word_table_mt search_word_table_mt  search_word_table_mt_up
