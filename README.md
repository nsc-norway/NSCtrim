# NSCtrim - FASTQ level primer trimming tool

Reasonably fast tool to trim pairs of primer sequences from sequencing reads.
It supports paired-end short read data in FASTQ format. The reads should contain
a primer sequence at the start of read 1 and read 2. Reads without a valid
primer combination are filtered, not written to the output. It can be used on
data where you know the exact list of primer (or other prefix) seuqences that
should be present in the reads.

It attempts to match all the provided read 1 primer sequences to the beginning of
the read 1 sequence. For the matches, it checks the read 2 primers. It can use a
simple string equality test, or a dynamic programming approach if the number of 
mismatches allowed (`-m`) is greater than zero.



## Usage

```
 ./NSCtrim [options] \
     PRIMER_FILE \
     INPUT_FILE_R1 INPUT_FILE_R2 \
     OUTPUT_FILE_R1 OUTPUT_FILE_R2 

Allowed options:
  -m [ --mismatches-per-primer ] arg (=0)
                                        Maximum allowed mismatches in primer1 
                                        and primer2 (per primer).
  -s [ --swapped-primer-pairs ]         Also search for reverse primer in read 
                                        1 and forward primer in read 2 
                                        (non-polar amplicons).
  -h [ --help ]                         Show this help message.

```

Example:
```
./NSCtrim primers.txt reads_R1.fastq.gz reads_R2.fastq.gz output_R1.fastq.gz output_R2.fastq.gz
```

### Primer file

The primer file should be an 3-column tab separated file without a header. 

| Column | Value        |
|--------|--------------|
| 1     | Read 1 primer sequence. |
| 2     | Read 2 primer sequence. |
| 3     | Name of primer combination. Used in the output table only. |

Example:

```
AAAGGTTTATACCCTTCCCAGG  AGGCAAACTGAGTTGGACGTG amp001
AACCAACTTTCGATCTCTTGTAG AGGCAAACTGAGTTGGACGTG amp002
```

The second column (Read 2 primer sequence) should contain the actual sequence to search for in
the reads, which is the reverse complement of the DNA sequence to amplify.

**Important warning for short amplicons:** 

If the length of any of the amplicons is less than the length of the sequence reads, extra caution
is required. For example, if read 1 and read 2 are both 151 bases long, then there will be an issue
if the sum of the length of one of the primers and its amplified region is less than 151 bases.
Both read 1 and read 2 will then start with a primer sequence, then contain the full
amplified region, and then also contain the other primer at the end of the read -- either in full or
parts of it. NSCtrim makes no attempt to trim primers that are not at the beginning of the reads.

The recommended way to deal with short amplicons is to clean the data with `fastp` 
( https://github.com/OpenGene/fastp ) AFTER running NSCtrim. fastp removes adapter sequences based
on overlap between read 1 and read 2. If the output of NSCtrim is processed with fastp, it will
treat the primer sequences at the end of the reads as adapters, and remove them.


### Bridging amplicons

Depending on the protocol used for preparing the amplicons, you may have longer amplicons
that bridge several primer sites. These have to be explicitly included in the primer file in order to
be written to the output.

For example if you have two amplicons positioned in the genome like this:
```
-----F1------R1-------F2-----R2-----
```
You may possibly get a product of forward primer F1 and reverse primer R2.


### Swapped primer paris

Some amplicon protocols create "non-polar" amplicons, where they can occur in the reverse
order (this is quite unusual).

Use the option `-s` to search for combinations that have swapped the read 1 and read 2 primers.

In such a protocol if you have the following two amplicons:

| Forward primer | Reverse primer        |
|--------|--------------|
| A123   | B123 |
| C123   | D123 |


You would get this in the data:

| Read 1 | Read 2   |
|--------|------|
| A123   | B123 |
| B123   | A123 |
| C123   | D123 |
| D123   | C123 |

The sequences are not reverse-complemented in any of these cases (apart from the fact that
B123 and D123 should already be the reverse-complement of the amplified template DNA, as
described above).


### Output

STDERR includes detalis about the tool version, parameters, progress and error status.

A table of the number of found primer sequences is written to STDOUT. This can optionally
be redirected to a file using the shell.


## Installation

### Binary download

Statically compiled executables for Linux are available on the release pages in Github.

### Docker

You can build the program from source inside docker, creating a docker image with the NSCtrim tool.
Clone this git repository, and then run `make docker`. The docker image `nsctrim:VERSION` will be created locally.
The version number is determined from the currently checked out git tag, but you can override this version by
setting the environment variable VERSION. The command `NSCtrim` is available on the path inside the docker image.

You can even create a docker image more directly using the github URL. You have to give it the version manually, for the purpose of printing the version in the tool's output (suggestions for how to fix it are welcome). Specify any tag to build -- 1.3.1 in the example:

```
docker build -t nsctrim:1.3.1 https://github.com/nsc-norway/NSCtrim.git#1.3.1 --build-arg VERSION=1.3.1
```

### Compile from source

You can build the program using a C++ compiler.

Requirements:

* C++ compiler
* GNU Make
* Development files for the following libraries, from Boost:
    * `iostreams`
    * `program-options`
* zlib (`libz`)

As an example, you can install the following packages with `apt` on Ubuntu 20.04 to install the requirements:

* `build-essential`
* `libboost-iostreams-dev`
* `libboost-program-options-dev`
* `libz-dev`


Run `make` to make a dynamically linked executable.

Run `make NSCtrim.static` to make a statically linked executable. This can be used on Linux computers that don't have Boost installed. (you can rename the file to `NSCtrim`)

The version string is taken automatically from the currently committed git tag. To define your own version, commit the changes (locally) and create a tag with `git tag`.

There is no installation, but you can copy the file to `/usr/bin` to be able to run it on the command line.
