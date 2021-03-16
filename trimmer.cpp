#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <queue>
#include <memory>
#include <chrono>
#include <sstream>
#include <fstream>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/program_options.hpp>

const size_t BATCH_SIZE = 8192; // reads to input, analyse
const size_t QUEUE_MAX = 16;        // number of batches in input queue
const size_t QUEUE_HIGH_LEVEL = 12; // max number of batch before running input
const size_t OUTPUT_QUEUE_SUM = 65535; // max read in ouput queues
const std::streamsize GZIP_INPUT_BUFFER_SIZE = 16*1024;
const std::streamsize GZIP_OUTPUT_BUFFER_SIZE = 4*1024;

using namespace std;
using namespace boost::iostreams;
using namespace std::chrono;
namespace po = boost::program_options;

#include "bounded_levenshtein_distance.cpp"

// General data classes for primer pair
class PrimerPair {
public:
    std::array<string, 2> primers;
    string name;
    unsigned long n_read_pairs = 0;
    
    PrimerPair(const string name, string primer1, string primer2) :
        name(name), primers({primer1, primer2}) {}

};


// Infrastructure for thread pool and work queues

class Batch {
    public:
    size_t num_in_batch;
    std::array<string, 4> data[BATCH_SIZE];
};


class FastqInput {
    unique_ptr<istream> input_stream;

    public:
    string path;
    size_t input_sizes[4] = {1};

    FastqInput(const string& path, bool compressed) : path(path) {
        if (compressed) {
            filtering_istream* fis = new filtering_istream;
            input_stream.reset(fis);
            file_descriptor_source raw(path);
            // Construct with buffer size. 15 is the default value of the first parameter.
            fis->push(gzip_decompressor(15, GZIP_INPUT_BUFFER_SIZE));
            fis->push(raw);
        }
        else {
            input_stream.reset(new ifstream(path, ios::binary));
        }
    }

    bool readBatch(Batch& bat) {
        int i;
        for (i=0; i<BATCH_SIZE; ++i) {
            for (int j=0; j<4; ++j) bat.data[i][j].reserve(input_sizes[j]);
            getline(*input_stream, bat.data[i][0]);
            input_sizes[0] = max(input_sizes[0], bat.data[i][0].capacity());
            if (input_stream->good()) {
                for (int j=1; j<4; ++j) {
                    getline(*input_stream, bat.data[i][j]);
                    input_sizes[j] = max(input_sizes[j], bat.data[i][j].capacity());
                }
            }
            else if (input_stream->eof()) {
                break;
            }
            else {
                return false;
            }
        }
        bat.num_in_batch = i;
        return true;
    }
    
    bool eof() {
        return input_stream->eof();
    }
};


class OutputJob {
    public:
    int trim;
    std::array<string, 4>& data;
    OutputJob(int trim, std::array<string, 4>& data) :
        trim(trim), data(data) {}
};


class FastqOutput {
    shared_ptr<ostream> out_stream;

    public:
    string path;

    bool openFile(string path, bool compressed) {
        if (compressed) {
            filtering_ostream* fos = new filtering_ostream;
            out_stream.reset(fos);
            file_descriptor_sink fds(path);
            fos->push(gzip_compressor(zlib::default_compression, GZIP_OUTPUT_BUFFER_SIZE));
            fos->push(fds);
            return (fds.is_open());
        }
        else {
            out_stream.reset(new ofstream(path, ios::binary));
            return out_stream->good();
        }
    }

    void writeBuffers(vector<OutputJob>& outputs) {
        for(OutputJob & output : outputs) {
            (*out_stream) << output.data[0] << '\n';
            (*out_stream) << output.data[1].substr(output.trim) << '\n';
            (*out_stream) << output.data[2] << '\n';
            (*out_stream) << output.data[3].substr(output.trim) << '\n';
        }
    }
};


class AnalysisResult {
    // Contains a bool to say whether the primers were found, and the lengths to trim.
    public:
    int found;
    int trim[2];
};


int mismatches(const string& templ, const string& test) {
    unsigned int mismatch = 0;
    size_t len = templ.length();
    for (int i = 0; i < len; ++i) {
        if (templ[i] != test[i]) mismatch++;
    }
    return mismatch;
}

class Analysis {

    vector<PrimerPair>& primer_pairs;

    bool use_levens;
    unsigned int primer_mismatches;

    public:

    unsigned long unknown_read_pairs = 0;
    unsigned long n_total_read_pairs = 0;


    Analysis(vector<PrimerPair>& primer_pairs, bool use_levens, unsigned int primer_mismatches) :
        primer_pairs(primer_pairs), use_levens(use_levens),
        primer_mismatches(primer_mismatches)
    {
    }

    vector<AnalysisResult> analyseReads(std::array<Batch, 2>& bat) {

        vector<AnalysisResult> results(bat[0].num_in_batch);

        int odata = 0; // Output record
        for (int idata=0; idata < bat[0].num_in_batch; ++idata) {
            bool unknown = true;
            int n_trim_r[2];

            for (PrimerPair & pp : primer_pairs) {
                int n_mismatches = 0;
                if (primer_mismatches == 0) {
                    // Check if the reads start with the primer sequences, and trim if they do
                    if (bat[0].data[idata][1].rfind(pp.primers[0], 0) == 0&&
                            bat[1].data[idata][1].rfind(pp.primers[1], 0) == 0) {
                        n_trim_r[0] = pp.primers[0].length();
                        n_trim_r[1] = pp.primers[1].length();
                        n_mismatches = 0;
                    }
                    else {
                        n_mismatches = 1;
                    }
                }
                else if (use_levens) {
                    for (int read=0; read<2 && n_mismatches <= primer_mismatches; ++read) {
                        auto result = mismatch_and_alignment(primer_mismatches + 1,
                                                    pp.primers[read], bat[read].data[idata][1]);
                        n_mismatches += result.first;
                        n_trim_r[read] = result.second;
                    }
                }
                else {
                    for (int read=0; read<2 && n_mismatches <= primer_mismatches; ++read) {
                        if (bat[read].data[idata][1].length() < pp.primers[read].length()) {
                            n_mismatches = primer_mismatches + 1; // Too short, dropped instantly
                            break;
                        }
                        n_mismatches += mismatches(pp.primers[read], bat[read].data[idata][1]);
                        n_trim_r[0] = pp.primers[read].length();
                    }
                }
                if (n_mismatches <= primer_mismatches) {
                    results[idata].found = true;
                    results[idata].trim[0] = n_trim_r[0];
                    results[idata].trim[1] = n_trim_r[1];
                    pp.n_read_pairs += 1;
                    unknown = false;
                    break;
                }
            }
            // No match found, output untrimmed FQ
            if (unknown) {
                unknown_read_pairs++;
                results[idata].found = false;
            }
            if (++n_total_read_pairs % 1000000 == 0) {
                cerr << "Processed " << n_total_read_pairs << " read pairs. Unmatched primers: "
                     << (unknown_read_pairs * 100.0 / n_total_read_pairs) << " %."<< endl;
            }
        }
        return results;
    }
};


class TrimmingManager {
 
    std::array<FastqInput, 2>& inputs;
    Analysis & analysis;
    std::array<FastqOutput, 2>& outputs;

    public:
    // Run times
    double input_time=0, matching_time=0, output_time=0;
    TrimmingManager(    std::array<FastqInput, 2>& inputs,
                        Analysis& analysis,
                        std::array<FastqOutput, 2>& outputs) 
        : inputs(inputs), analysis(analysis), outputs(outputs)
    {
    }


    // Main function reads, trims and writes
    bool run() {
        std::array<Batch, 2> bat;
        auto time1 = high_resolution_clock::now();
        while (!inputs[0].eof() && !inputs[1].eof()) {
            for (int read=0; read<2; ++read) {
                if (!inputs[read].readBatch(bat[read])) {
                    cerr << "READ ERROR" << endl;
                    return false;
                }
            }
            auto time2 = high_resolution_clock::now();
            input_time += duration_cast<microseconds>(time2 - time1).count();
            vector<AnalysisResult> results = analysis.analyseReads(bat);
            auto time3 = high_resolution_clock::now();
            matching_time += duration_cast<microseconds>(time3 - time2).count();
            for (int read=0; read<2; ++read) {
                vector<OutputJob> output_jobs;
                for (int i=0; i<results.size(); ++i) {
                    if (results[i].found) {
                        output_jobs.emplace_back(results[i].trim[read], bat[read].data[i]);
                    }
                }
                outputs[read].writeBuffers(output_jobs);
            }
            time1 = high_resolution_clock::now();
            output_time += duration_cast<microseconds>(time1 - time3).count();
        }
        if (inputs[0].eof() && inputs[1].eof()) {
            return true;
        }
        else {
            cerr << "ERROR: Different length input files." << endl;
            return false;
        }
    }
};

vector<PrimerPair> getPrimerPairs(const string& primer_file) {
    ifstream ppfile(primer_file);
    if (!ppfile) {
        cerr << "Error: Unable to open primer file " << primer_file << endl;
        exit(1);
    }
    vector<PrimerPair> result;
    string line;
    while (ppfile) {
        getline(ppfile, line);
        if (!line.empty()) {
            stringstream ss_line(line);
            string seq1, seq2, dummy, name;
            ss_line >> seq1; ss_line >> seq2; ss_line >> dummy; ss_line >> name;
            result.emplace_back(PrimerPair(name, seq1, seq2));
        }
    }
    return result;
}


int main(int argc, char* argv[]) {

    const string usage =
                "usage:\n " + string(argv[0]) + " [options] \\\n"
                "     PRIMER_FILE \\\n" + 
                "     INPUT_FILE_R1 INPUT_FILE_R2 \\\n" + 
                "     OUTPUT_FILE_R1 OUTPUT_FILE_R2 \n\n";


    // Parse command line into these variables
    string    primer_file,
              input_file_r1, input_file_r2,
              output_file_r[2];
    unsigned int primer_mismatches;
    bool use_hamming;

    cerr << "\ntool " << VERSION << "\n" << endl;

    po::options_description visible("Allowed options");
    visible.add_options()
        ("primer-mismatches,b", po::value<unsigned int>(&primer_mismatches)->default_value(0),
            "Total allowed mismatches in primer1 + primer2.")
        ("use-hamming,H", po::bool_switch(&use_hamming),
            "Use Hamming distance instead of Levenshtein distance.")
        ("help,h", "Show this help message.")
    ;
    po::options_description positionals("Positional options(hidden)");
    positionals.add_options()
        ("PRIMER_FILE", po::value<string>(&primer_file)->required(),
            "Primer pair file")
        ("INPUT_FILE_R1", po::value<string>(&input_file_r1)->required(),
            "Input file read 1")
        ("INPUT_FILE_R2", po::value<string>(&input_file_r2)->required(),
            "Input file read 2")
        ("OUTPUT_FILE_R1", po::value<string>(&output_file_r[0])->required(),
            "Output r1.")
        ("OUTPUT_FILE_R2", po::value<string>(&output_file_r[1])->required(),
            "Output r2.")
    ;
    po::options_description all_options("Allowed options");
    all_options.add(visible);
    all_options.add(positionals);

    po::positional_options_description pos_desc;
    pos_desc.add("PRIMER_FILE", 1);
    pos_desc.add("INPUT_FILE_R1", 1);
    pos_desc.add("INPUT_FILE_R2", 1);
    pos_desc.add("OUTPUT_FILE_R1", 1);
    pos_desc.add("OUTPUT_FILE_R2", 1);

    po::variables_map vm;
    try {
        po::store(
                po::command_line_parser(argc, argv).options(all_options).positional(pos_desc).run(),
                vm
                );
        if (vm.count("help") > 0) {
            cerr << usage << visible << endl;
            return 0;
        }
        else {
            po::notify(vm);
        }
    }
    catch(po::error& e) 
    { 
      cerr << "ERROR: " << e.what() << "\n\n";
      cerr << usage << visible << endl; 
      return 1; 
    }

    bool use_levens = !use_hamming;

    vector<PrimerPair> primers = getPrimerPairs(primer_file);
    if (primers.empty()) {
        cerr << "Error: Barcode file is empty." << endl;
        return 1;
    }

    bool compressed =
            input_file_r1.size() > 3 &&
            input_file_r1.substr(input_file_r1.size() - 3) == ".gz" &&
            input_file_r2.size() > 3 &&
            input_file_r2.substr(input_file_r1.size() - 3) == ".gz";

    // Input files
    std::array<FastqInput, 2> inputs = {
        FastqInput(input_file_r1, compressed),
        FastqInput(input_file_r2, compressed)
    };

    // Output files
    std::array<FastqOutput, 2> outputs;

    // Open output files
    for (int read=0; read<2; ++read) {
        if (!outputs[read].openFile(output_file_r[read], compressed)) {
            cerr << "Failed to open output file " << output_file_r[read] << " for writing." << endl;
            exit(1);
        }
    }
    
    // Print information on startup
    cerr.precision(1);
    cerr << fixed;
    cerr << "\nTrimming " << primers.size() << " primer pairs...\n\n";
    cerr << " Allowed primer mismatches:  " << primer_mismatches << '\n';
    cerr << " String distance:            ";
    if (primer_mismatches == 0) cerr << "NoMismatch" << '\n';
    else if (use_levens)        cerr << "Levenshtein" << '\n';
    else                        cerr << "Hamming" << '\n';
    cerr << " Input/output compression:   " << (compressed ? "gzip" : "off") << '\n';
    cerr << endl;

    Analysis analysis(primers, use_levens, primer_mismatches);
    TrimmingManager manager(inputs, analysis, outputs);
    bool success = manager.run();

    if (!success) { // The Manager will print an error message
        cerr << "\nThe program exited due to an error." << endl;
        return 1;
    }
    else {
        bool delete_files = true;
        for (int i=0; i<primers.size(); ++i) { // Deletes empty files, would be corrupted gzip files.
            delete_files = delete_files && (primers[i].n_read_pairs == 0);
        }
        cerr << "\nCompleted trimming " << analysis.n_total_read_pairs << " read pairs.\n" << endl;
        if (delete_files) {
            cerr << "Removing empty output files.\n" << endl;
            unlink(outputs[0].path.c_str());
            unlink(outputs[1].path.c_str());
        }
        cerr << "\nTime (seconds): Input: " << manager.input_time / 1e6
             << ", Matching: " << manager.matching_time / 1e6
             << ", Output: " << manager.output_time / 1e6 << endl;
        cerr << "\nOutput: " << analysis.n_total_read_pairs - analysis.unknown_read_pairs << " read pairs.\n" << endl;
        cout << "PRIMER_NAME\tNUM_READ_PAIRS\tPCT_READS\n";
        for (PrimerPair& pp : primers) {
            cout.precision(2);
            cout << fixed;
            cout << pp.name << '\t'
                << pp.n_read_pairs << '\t'
                << pp.n_read_pairs * 100.0 / max(analysis.n_total_read_pairs, 1ul) << '\t'
                << '\n';
        }
        cout << "Unmatched\t-\t-\t" << analysis.unknown_read_pairs << '\t'
            << analysis.unknown_read_pairs * 100.0 / analysis.n_total_read_pairs << "\n";
    }
    return 0;
}


