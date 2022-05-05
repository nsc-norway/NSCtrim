FROM ubuntu:20.04 AS BUILDER

ENV DEBIAN_FRONTEND "noninteractive"
# hadolint ignore=DL3009
RUN apt-get update && \
    apt-get -y install --no-install-recommends \
    build-essential \
    libboost-iostreams-dev \
    libboost-program-options-dev \
    zlib1g-dev
ARG VERSION
RUN mkdir /trimmer
WORKDIR /trimmer
COPY NSCtrim.cpp .
COPY bounded_levenshtein_distance.cpp .
COPY Makefile .

RUN make NSCtrim.static

FROM alpine:3.13.5 AS RUNNER
RUN apk add bash 
COPY --from=BUILDER /trimmer/NSCtrim.static /usr/bin/NSCtrim
