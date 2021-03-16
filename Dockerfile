FROM ubuntu:20.04
ARG VERSION
ENV DEBIAN_FRONTEND "noninteractive"
RUN apt-get update && apt-get -y install build-essential libboost-iostreams-dev libboost-program-options-dev libz-dev
RUN mkdir /trimmer

COPY trimmer.cpp /trimmer
COPY bounded_levenshtein_distance.cpp /trimmer
COPY Makefile /trimmer
ARG VERSION
RUN (cd /trimmer && make && mv /trimmer/trimmer /usr/bin)
