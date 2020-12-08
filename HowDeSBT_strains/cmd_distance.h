#include "bit_vector.h"
#include "bloom_filter.h"
#include <string>
#include <vector>
#include "utilities.h"
#include "support.h"
#include "commands.h"
using namespace std;

class DistanceGraph
{
public:
    DistanceGraph();
    DistanceGraph(double** matrix, double threshold, size_t size);
    ~DistanceGraph();
    void dump_graph(string labels_path);
    void dump_adj_matrix();
    void compute_cliques();
    //void find_cliques_k(int i, int l, int k);
    //void add(int n, int k);
    //bool is_clique(int v);
private:
    void _build_graph();
    
public:
    map<int, vector<vector<int>>> cliques;
private:
    uint8_t**   _graph;
    double**    _matrix;
    double      _threshold;
    size_t      _nbv;
    int*        _degre;
    vector<int> carray;
};
class DistanceCommand : public Command
{
public:
    static const std::uint64_t defaultEndPosition = 0;
    static const int default_max_size_cl = 10;
    static constexpr double default_threshold = 0.5f;

public:
    DistanceCommand(const std::string &name) : Command(name) {}
    virtual ~DistanceCommand();
    virtual void short_description(std::ostream &s);
    virtual void usage(std::ostream &s, const std::string &message = "");
    virtual void parse(int _argc, char **_argv);
    virtual int execute(void);

    virtual void get_vectors();
    virtual void compute_hamming();
    virtual void dump_matrix();
    virtual void load_matrix(string path);
    //virtual void compute_graph();
    //virtual void dump_graph();
    //virtual void merge_close_bf(double threshold);
    virtual void disk_hamming();
    virtual void merge();

    bool full_vec;
    bool mergeAll;
    std::string listFilename;
    std::string path_matrix;
    std::string matFilename;
    std::uint64_t startPosition; // origin-zero, half-open
    std::uint64_t endPosition;
    size_t n;
    std::vector<string> paths;
    std::vector<BitVector*> bf_vectors;
    std::vector<BloomFilter*> dvec;
    double threshold;
    double** matrix;
    int cl_max_size;
    DistanceGraph* graph;
    vector<string> labels;
};

