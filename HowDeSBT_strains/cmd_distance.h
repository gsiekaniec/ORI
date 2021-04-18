#include "bit_vector.h"
#include "bloom_filter.h"
#include <string>
#include <vector>
#include "utilities.h"
#include "bit_utilities.h"
#include "support.h"
#include "commands.h"

#include <emmintrin.h>
#include <immintrin.h>

using namespace std;

using u64 = uint64_t;
using u32 = uint32_t;
using u8 = uint8_t;

static const u8 poplook8[256] =
{
    0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};

template<typename T>
uint64_t sse_hamming(
  T* x,
  T* y,
  size_t nb_bytes)
{
  u64 ones = 0;
  u32 n = nb_bytes;
  __m128i xi, yi, r, h;
  for (u32 i=0; i<nb_bytes; i+=16)
  {
    if (n > 0 && n < 16) break;
    xi = _mm_loadu_si128(reinterpret_cast<__m128i*>(x+i));
    yi = _mm_loadu_si128(reinterpret_cast<__m128i*>(y+i));
    r  = _mm_xor_si128(xi, yi);
    h  = _mm_unpackhi_epi64(r, r); 
    ones += static_cast<u64>(_mm_popcnt_u64(_mm_cvtsi128_si64(r)));
    ones += static_cast<u64>(_mm_popcnt_u64(_mm_cvtsi128_si64(h)));
    n -= 16;
  }
  for (u32 i=nb_bytes-1; i>=nb_bytes-n; i--)
    ones += poplook8[(x[i]^y[i])];
  return ones;
}

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
    virtual void compute_hamming2();
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

