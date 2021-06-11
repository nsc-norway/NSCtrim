FROM ubuntu:20.04 AS BUILDER

ENV DEBIAN_FRONTEND "noninteractive"
# hadolint ignore=DL3009
RUN apt-get update && \
    apt-get -y install --no-install-recommends \
    build-essential=12.8ubuntu1 \
    libboost-iostreams-dev=1.71.0.0ubuntu2 \
    libboost-program-options-dev=1.71.0.0ubuntu2 \
    zlib1g-dev=1:1.2.11.dfsg-2ubuntu1.2
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
