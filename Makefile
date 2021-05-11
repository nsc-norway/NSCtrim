all: NSCtrim

CFLAGS += -O3 -std=c++11

# If building locally, VERSION is not set by anything,
# so this will set it to the current tag.
# If using "docker-build.sh", use the environment variable
# passed through docker.
VERSION ?= $(shell git describe --tags --dirty )

NSCtrim: NSCtrim.cpp bounded_levenshtein_distance.cpp
	$(CXX) -o $@ NSCtrim.cpp $(CFLAGS) -lboost_program_options$(BOOST_LIB_SUFF) -lboost_iostreams$(BOOST_LIB_SUFF) -lz -DVERSION=\"${VERSION}\"

NSCtrim.static: NSCtrim.cpp bounded_levenshtein_distance.cpp
	$(CXX) -o $@ NSCtrim.cpp -static $(CFLAGS) -lboost_program_options$(BOOST_LIB_SUFF) -lboost_iostreams$(BOOST_LIB_SUFF) -lz -DVERSION=\"${VERSION}\"

docker:
	docker build -t nsctrim:${VERSION} --build-arg VERSION="${VERSION}" .

clean:
	rm -f NSCtrim
