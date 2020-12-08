#include "cmd_distance.h"
#include "utilities.h"
#include "bloom_filter.h"
#include "bit_vector.h"
#include "bit_utilities.h"
#include <string>

#define u32 std::uint32_t
#define u64 std::uint64_t

using namespace std;
void DistanceCommand::short_description(std::ostream &s)
{
    s << commandName << "-- Compute BFs distances" << endl;
}

void DistanceCommand::usage(std::ostream &s,
                            const string &message)
{
    if (!message.empty())
    {
        s << message << endl;
        s << endl;
    }

    short_description(s);
    s << "usage: " << commandName << " [options]" << endl;
    //    123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
    s << "  --list=<filename> file containing a list of bloom filters; only" << endl;
    s << "                    filters with uncompressed bit vectors are allowed" << endl;
    s << "  <filename>        same as --list=<filename>" << endl;
    s << "  --out=<filename>  name for the output matrix" << endl;
    s << "                    (by default this is derived from the list filename)" << endl;
    s << "  <start>..<end>    interval of bits to use from each filter; distance" << endl;
    s << "                    is computed only considers this subset of each filter's bits" << endl;
    s << "                    (by default we use the first " << defaultEndPosition << " bits)" << endl;
    s << "  --bits=<N>        number of bits to use from each filter; same as 0..<N>" << endl;
    s << "Merge options:" << endl;
    s << "  --threshold=<N>   threshold, floating point number ]0,1[" << endl;
    s << "  --matrix=<F>      hamming matrix file" << endl;
    s << "  --merge           merge maximal cliques" << endl;
}
void DistanceCommand::parse(int _argc,
                            char **_argv)
{
    int argc;
    char **argv;

    // defaults

    startPosition = 0;
    endPosition = defaultEndPosition;
    full_vec = true;
    mergeAll = false;
    threshold = default_threshold;
    cl_max_size = default_max_size_cl;
    // skip command name

    argv = _argv + 1;
    argc = _argc - 1;
    if (argc <= 0)
        chastise();

    //////////
    // scan arguments
    //////////

    for (int argIx = 0; argIx < argc; argIx++)
    {
        string arg = argv[argIx];
        string argVal;
        if (arg.empty())
            continue;

        string::size_type argValIx = arg.find('=');
        if (argValIx == string::npos)
            argVal = "";
        else
            argVal = arg.substr(argValIx + 1);

        // --help, etc.

        if ((arg == "--help") || (arg == "-help") || (arg == "--h") || (arg == "-h") || (arg == "?") || (arg == "-?") || (arg == "--?"))
        {
            usage(cerr);
            std::exit(EXIT_SUCCESS);
        }

        if ((arg == "--help=debug") || (arg == "--help:debug") || (arg == "?debug"))
        {
            debug_help(cerr);
            std::exit(EXIT_SUCCESS);
        }

        // --list=<filename>

        if (is_prefix_of(arg, "--list="))
        {
            if (not listFilename.empty())
                chastise("unrecognized option: \"" + arg + "\""
                                                           "\nbloom filters list file was already given as \"" +
                         listFilename + "\"");
            listFilename = argVal;
            continue;
        }

        // --out=<filename>, --tree=<filename>, etc.

        if ((is_prefix_of(arg, "--out=")) || (is_prefix_of(arg, "--output=")))
        {
            matFilename = argVal;
            continue;
        }

        // --bits=<N>

        if ((is_prefix_of(arg, "--bits=")) || (is_prefix_of(arg, "B=")) || (is_prefix_of(arg, "--B=")))
        {
            startPosition = 0;
            full_vec = false;
            endPosition = string_to_unitized_u64(argVal);
            continue;
        }

        // --threshold=<N>

        if (is_prefix_of(arg, "--threshold="))
        {
            threshold = atof(argVal.c_str());
            continue;
        }

        if (is_prefix_of(arg, "--matrix="))
        {
            path_matrix = argVal;
            continue;
        }

        if (is_prefix_of(arg, "--max="))
        {
            cl_max_size = string_to_unitized_u32(argVal);
            continue;
        }

        if (arg == "--merge")
        {
            mergeAll = true;
            continue;
        }
        // unrecognized --option

        if (is_prefix_of(arg, "--"))
            chastise("unrecognized option: \"" + arg + "\"");

        // <start>..<end>

        std::size_t separatorIx = arg.find("..");
        if (separatorIx != string::npos)
        {
            startPosition = string_to_unitized_u64(arg.substr(0, separatorIx));
            endPosition = string_to_unitized_u64(arg.substr(separatorIx + 2));
            full_vec = false;
            if (endPosition <= startPosition)
                chastise("bad interval: " + arg + " (end <= start)");
            continue;
        }

        // <filename>

        if (not listFilename.empty())
            chastise("unrecognized option: \"" + arg + "\""
                                                       "\nbloom filters list file was already given as \"" +
                     listFilename + "\"");
        listFilename = arg;
    }

    // sanity checks

    if (startPosition % 8 != 0)
        chastise("the bit interval's start (" + std::to_string(startPosition) + ") has to be a multiple of 8");

    if (listFilename.empty())
        chastise("you have to provide a file, listing the bloom filters for the tree");

    if (matFilename.empty())
    {
        string::size_type dotIx = listFilename.find_last_of(".");
        if (dotIx == string::npos)
            matFilename = listFilename + ".mat";
        else
            matFilename = listFilename.substr(0, dotIx) + ".mat";
    }
    return;
}

void DistanceCommand::get_vectors()
{
    std::ifstream in(listFilename);
    if (!in)
        fatal("Unable to open " + listFilename);

    BloomFilter *f = nullptr;
    bool list = true;
    string bfname;
    int line = 0;

    while (std::getline(in, bfname))
    {
        line++;
        paths.push_back(bfname);
        BloomFilter *bf = new BloomFilter(strip_blank_ends(bfname));
        bf->preload();
        BitVector *bv = bf->get_bit_vector();
        if (bv->compressor() != bvcomp_uncompressed)
            fatal("Distance computation needs uncompressed vectors");

        if (f == nullptr)
        {
            f = bf;
            list = false;
            if (bf->numBits <= startPosition)
                fatal("error: " + bfname + " has only " + std::to_string(bf->numBits) + " bits" + ", so the bit interval " + std::to_string(startPosition) + ".." + std::to_string(endPosition) + " would be empty");
            if (bf->numBits < endPosition)
            {
                endPosition = bf->numBits;
                cerr << "warning: reducing bit interval to " << startPosition << ".." << endPosition << endl;
            }
        }
        else
        {
            bf->is_consistent_with(f, true);
        }
        
        string bvname = bv->filename;
        size_t offset = bv->offset;
        if (bf != f) delete bf;

        size_t start = offset + sdslbitvectorHeaderBytes + startPosition/8;
        string raw = bvname+":raw"+":"+to_string(start)+":"+to_string(endPosition-startPosition);
        
        BitVector* rawBv = BitVector::bit_vector(raw);
        bf_vectors.emplace_back(rawBv);
    }

    if (f != nullptr) delete f;
    in.close();

    if(list)
        fatal("Error, " + listFilename + " contains no bloom filters");
}

void DistanceCommand::compute_hamming()
{   
    n = bf_vectors.size();
    u64 nbits = endPosition - startPosition;
    u64 nbytes = (nbits + 7) / 8;

    matrix = new double*[n]();
    for (int i=0; i<n; i++)
        matrix[i] = new double[n]();

    uint64_t* bv_vectors[n];

    for (u32 v=0; v<n; v++)
    {
        BitVector* bv = bf_vectors[v];
        bv->load();
        bv_vectors[v] = bv->bits->data();
    }

    u64 dham;
    double rham;
    for (u32 v=0; v<n-1; v++)
    {
        for (u32 w=v+1; w<n; w++)
        {
            dham = hamming_distance(bv_vectors[v], bv_vectors[w], nbits);
            rham = 1.0*dham/nbits;
            matrix[v][w] = rham;
            matrix[w][v] = rham;
        }
    }
}
void DistanceCommand::disk_hamming()
{
    std::ifstream in(listFilename);
    if (!in)
        fatal("Unable to open " + listFilename);

    string bfname;
    int line = 0;

    while (std::getline(in, bfname))
    {
        line++;
        paths.push_back(bfname);
        BloomFilter *bf = new BloomFilter(strip_blank_ends(bfname));
        bf->preload();
        dvec.emplace_back(bf);
    }
    in.close();

    n = dvec.size();
    u64 nbits = dvec[0]->numBits;
    u64 nbytes = (nbits+7) / 8;

    matrix = new double*[n]();
    for (int i=0; i<n; i++)
        matrix[i] = new double[n]();

    u64 dham;
    double rham;
    for (u32 v=0; v<n; v++)
    {
        dvec[v]->load();
        u64* first = dvec[v]->bvs[0]->bits->data();
        for(u32 w=v+1; w<n; w++)
        {
            dvec[w]->load();
            u64* second = dvec[w]->bvs[0]->bits->data();
            dham = hamming_distance(first, second, nbits);
            //cout << 
            rham = 1.0*dham/nbits;
            matrix[v][w] = rham;
            matrix[w][v] = rham;
            dvec[w]->preload();
        }
        delete dvec[v];
    }
}


void DistanceCommand::dump_matrix()
{
    ofstream tsv_mat("hamming_matrix.tsv", ios::out);
    ofstream bin_mat("hamming_matrix.bin", ios::out | ios::binary);

    bin_mat.write((char*)&n, sizeof(size_t));

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            tsv_mat << std::fixed << std::setprecision(4) << matrix[i][j] << " ";
            bin_mat.write((char*)&matrix[i][j], sizeof(double));
        }
        tsv_mat << "\n";
    }
    tsv_mat.close();
    bin_mat.close();
}

void DistanceCommand::load_matrix(string path)
{
    ifstream bin_mat(path, ios::in | ios::binary);
    if (!bin_mat.good()) fatal("Unable to open " + path);
    bin_mat.read((char*)&n, sizeof(size_t));

    matrix = new double*[n];
    for (int i=0; i<n; i++)
    {
        matrix[i] = new double[n];
        bin_mat.read((char*)matrix[i], sizeof(double)*n);
    }
}

void DistanceCommand::merge()
{
    map<int, vector<vector<int>>> cliques = graph->cliques;
    map<int, vector<vector<int>>>::iterator it;

    vector<string> labels;
    std::ifstream in(listFilename);
    if (!in)
        fatal("Unable to open " + listFilename);
    string label;
    while (std::getline(in, label))
    {
        size_t last = label.find_last_of("/");
        if (string::npos != last)
            label.erase(0, last + 1);
        size_t ext = label.find('.');
        if (string::npos != ext)
            label.erase(ext);
        labels.push_back(label);
    }
    for (it = cliques.begin(); it != cliques.end(); it++)
    {
        string filename = "";
        for (auto& vec: it->second)
        {
            filename += to_string(vec[0]) + "_";
            BloomFilter *bf = new BloomFilter(paths[vec[0]]);
            bf->load();
            for (int i=1; i<vec.size(); i++)
            {
                BloomFilter *bf_m = new BloomFilter(paths[vec[i]]);
                filename += to_string(vec[i]) + "_";
                bf_m->load();
                bf->union_with(bf_m->bvs[0]);
                delete bf_m;
                cout << filename << endl;
            }
            filename += ".bf";
            bf->filename = filename;
            bf->save();
            delete bf;
        }
    }
}
DistanceCommand::~DistanceCommand()
{
    for (int i=0; i<n; i++)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
}

int DistanceCommand::execute()
{
    if (!mergeAll)
    {
        get_vectors();
        compute_hamming();
        dump_matrix();
    }
    else if (mergeAll && path_matrix.size())
    {
        load_matrix(path_matrix);
        graph = new DistanceGraph(matrix, threshold, n);
        graph->dump_graph(listFilename);
        graph->dump_adj_matrix();
        graph->compute_cliques();
        merge();
    }
    return 0;
}

DistanceGraph::DistanceGraph(double** matrix, double threshold, size_t size)
    : _matrix(matrix), _threshold(threshold), _nbv(size)
{
    carray.reserve(100);
    _build_graph();
}

DistanceGraph::DistanceGraph() {}

DistanceGraph::~DistanceGraph()
{
    for (int i=0; i<_nbv; i++)
    {
        delete[] _graph[i];
    }
    delete[] _graph;
    delete[] _degre;
}

void DistanceGraph::_build_graph()
{
    _degre = new int[_nbv]();
    _graph = new uint8_t*[_nbv]();
    for (int i=0; i<_nbv; i++)
        _graph[i] = new uint8_t[_nbv]();

    for (int i=0; i<_nbv; i++)
    {
        for (int j=i+1; j<_nbv; j++)
        {
            if (_matrix[i][j] <= _threshold)
            {
                _graph[i][j] = 1;
                _graph[j][i] = 1;
                _degre[i]++;
                _degre[j]++;
            }
        }
    }
}

void DistanceGraph::dump_adj_matrix()
{
    ofstream adj_mat("adj_matrix.txt", ios::out);
    //adj_mat << to_string(_nbv) << "\n";
    for (int i=0; i<_nbv; i++)
    {
        for (int j=0; j<_nbv; j++)
        {
            adj_mat << to_string(_graph[i][j]) << " ";
        }
        adj_mat << endl;
    }
}

void DistanceGraph::compute_cliques()
{
    char buffer[1024];
    readlink("/proc/self/exe", buffer, 1024);
    string dir_bin = dirname(buffer);
    string cmd = "python3 " + dir_bin + "/max_cliques.py adj_matrix.txt max_cliques.txt > /dev/null";
    system(cmd.c_str());
    ifstream in_f("max_cliques_wo_overlap.txt", ios::in);
    string line;
    int n = 0;
    while (getline(in_f, line))
    {
        stringstream ss(line);
        istream_iterator<std::string> begin(ss);
        istream_iterator<std::string> end;
        vector<string> vstrings(begin, end);
        vector<int> tmp;
        transform(vstrings.begin(), vstrings.end(), back_inserter(tmp), 
            [](const string& str) {return stoi(str);}
        );
        cliques[n].push_back(tmp);
        n++;
    }
}

void DistanceGraph::dump_graph(string labels_path)
{
    ofstream dot_graph("distance_graph.dot", ios::out);
    ofstream bin_graph("distance_graph.bin", ios::out | ios::binary);

    vector<string> labels;
    std::ifstream in(labels_path);
    if (!in)
        fatal("Unable to open " + labels_path);
    string label;
    while (std::getline(in, label))
    {
        size_t last = label.find_last_of("/");
        if (string::npos != last)
            label.erase(0, last + 1);
        size_t ext = label.find('.');
        if (string::npos != ext)
            label.erase(ext);
        labels.push_back(label);
    }
    bin_graph.write((char*)&_nbv, sizeof(size_t));

    for (int i=0; i<_nbv; i++)
    {
        bin_graph.write((char*)_graph[i], sizeof(uint8_t)*_nbv);
    }
    
    dot_graph << "graph G {" << "\n";

    for (int i=0; i<_nbv; i++)
    {
        dot_graph << to_string(i) << " [label=\"" << labels[i] << "\"]"<< "\n";
    }

    dot_graph << "\n";

    for (int i=0; i<_nbv; i++)
    {
        for (int j=i+1; j<_nbv; j++)
        {
            if (_graph[i][j])
            {
                dot_graph << to_string(i) << " -> " << to_string(j) << " ";
                dot_graph << "[weight=" << to_string(_matrix[i][j]) << "]" << "\n";
            }
        }
    }

    dot_graph << "}" << "\n";
}

//bool DistanceGraph::is_clique(int v)
//{
//    for (int i=1; i<v; i++)
//    {
//        for (int j=i+1; j<v; j++)
//        {
//            if (_graph[carray[i]][carray[j]] == 0)
//                return false;
//        }
//    }
//    return true;
//}
//
//void DistanceGraph::add(int n, int k)
//{
//    //if (cliques.count(k) == 0)
//    vector<int> temp;
//    for (int i=1; i<n; i++)
//        temp.push_back(carray[i]);
//    cliques[k].push_back(temp);
//}
//
//void DistanceGraph::find_cliques_k(int i, int l, int k)
//{
//    for (int j = i; j<_nbv-(k-l); j++)
//    {
//        if (_degre[j] >= k-1)
//        {
//            carray[l] = j;
//            if (is_clique(l+1))
//            {    
//                if (l<k)
//                {
//                    find_cliques_k(j, l+1, k);
//                }
//                else
//                {
//                    add(l+1, k);
//                }
//            }
//        
//        }
//    }
//
//}
