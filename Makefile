all: trimmer

CFLAGS += -O3 -std=c++11

# On RHEL7 -- Works inside devtoolset-7
trimmer: trimmer.cpp bounded_levenshtein_distance.cpp
	$(CXX) -g  -o $@ trimmer.cpp $(CFLAGS) -lboost_program_options$(BOOST_LIB_SUFF) -lboost_iostreams$(BOOST_LIB_SUFF) -lz -DVERSION=\"${VERSION}\"

clean:
	rm -f trimmer
