#ifndef query_H
#define query_H

#include <string>
#include <vector>
#include <iostream>

#include "bloom_filter.h"
#include "query.h"

using namespace std;
//----------
//
// classes in this module--
//
//----------

struct querydata
{
  std::uint32_t batchIx; // index of this query within a batch
  std::string name;
  std::string seq; // nucleotide sequence
};

class Query
{
public:
  Query(const querydata &qd, double threshold);
  virtual ~Query();

  virtual void kmerize(BloomFilter *bf, bool distinct = false, bool populateKmers = false);
  virtual void qgramize(BloomFilter *bf, vector<string> seeds, bool distinct=false, bool populateQgrams=false);
  virtual void sort_kmer_positions();
  virtual void dump_kmer_positions(std::uint64_t numUnresolved = -1);
  virtual std::uint64_t kmer_positions_hash(std::uint64_t numUnresolved = -1);

public:
  std::uint32_t batchIx; // index of this query within a batch
  std::string name;
  std::string seq;                          // nucleotide sequence
  double threshold;                         // search threshold
  std::vector<std::uint64_t> kmerPositions; // the kmers (converted to hash
                                            // .. values) corresponding to this
                                            // .. query; the first numUnresolved
                                            // .. entries are the yet-to-be-resolved
                                            // .. kmers; the resolved kmers are
                                            // .. moved to the tail
  std::vector<std::string> kmers;           // the kmers; this is only populated
                                            // .. in special instances (e.g. for
                                            // .. cmd_query_bf), and in those
                                            // .. cases care should be taken to
                                            // .. assure that the kmerPositions[ix]
                                            // .. corresponds to kmers[ix] for each
                                            // .. ix

  std::uint64_t numPositions;                     // total size of kmerPositions
  std::uint64_t neededToPass;                     // number of kmers required, to judge
                                                  // .. the query as a "pass"
  std::uint64_t neededToFail;                     // number of kmers required, to judge
                                                  // .. the query as a "fail"
  std::uint64_t numUnresolved;                    // number of kmers not yet known to be
                                                  // .. present or absent in all leaves
                                                  // .. of the current subtree
  std::uint64_t numPassed;                        // number of kmers known to be present
                                                  // .. in all leaves of the subtree
  std::uint64_t numFailed;                        // number of kmers known to be absent
                                                  // .. in all leaves of the subtree
  std::uint64_t nodesExamined;                    // number of nodes that were "examined"
                                                  // by this query
  bool adjustKmerCounts;                          // true  => populate matchesAdjusted[]
                                                  // false => don't
  std::vector<std::string> matches;               // names of leaves that match this query
  std::vector<std::uint64_t> matchesNumPassed;    // numPassed corresponding to
                                                  // .. each match; only valid if the
                                                  // .. search reached the leaf without
                                                  // .. having been pruned
  std::vector<std::uint64_t> matchesAdjustedHits; // adjusted value of
                                                  // .. numPassed (corresponding to each
                                                  // .. match), to .. account for
                                                  // .. estimated bloom filter false
                                                  // .. positives; only valid if the
                                                  // .. search reached the leaf without
                                                  // .. having been pruned, and only if
                                                  // .. adjustKmerCounts is true

  // stacks to maintain numUnresolved,
  // .. numPassed, and numFailed as we
  // .. move up and down the tree
  std::vector<std::uint64_t> numUnresolvedStack;
  std::vector<std::uint64_t> numPassedStack;
  std::vector<std::uint64_t> numFailedStack;
  std::vector<std::uint64_t> dbgKmerPositionsHashStack;

public:
  bool dbgKmerize = false;
  bool dbgKmerizeAll = false;

public:
  static void read_query_file(std::istream &in, const std::string &filename,
                              double threshold,
                              std::vector<Query *> &queries);
  static void read_query_file_qgram(const string& filename, double threshold, vector<Query*>& queries);
};

#endif // query_H
