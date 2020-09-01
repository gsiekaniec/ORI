#ifndef cmd_make_bf_qgram_H
#define cmd_make_bf_qgram_H

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <vector>

#include "commands.h"

class MakeBFQCommand : public Command
{
public:
  static const std::uint32_t defaultKmerSize = 20;
  static const std::uint32_t defaultNumHashes = 1;
  static const std::uint64_t defaultNumBits = 500 * 1000;
  static const std::uint64_t defaultUNumBits = defaultNumBits * 2;
  static const std::uint64_t defaultHashSeed = 0;

public:
  MakeBFQCommand(const std::string &name) : Command(name) {}
  virtual ~MakeBFQCommand() {}
  virtual void short_description(std::ostream &s);
  virtual void usage(std::ostream &s, const std::string &message = "");
  virtual void debug_help(std::ostream &s);
  virtual void parse(int _argc, char **_argv);
  virtual int execute(void);
  virtual void make_bloom_filter_qgram(void);
  virtual void make_multi_bloom_filter_qgram(void);
  virtual void make_multi_bloom_filter_qgram_intersect(void);
  virtual void report_stats(const std::string &bfOutFilename, std::uint64_t kmersAdded);
  virtual std::string build_output_filename(void);
  virtual std::string build_stats_filename(const std::string &bfOutFilename);

  std::string seed_file;
  std::string listFilename;
  std::vector<std::string> seqFilenames;
  std::string bfFilename;
  std::string asPerFilename;
  std::uint32_t kmerSize;
  std::uint64_t hashSeed1;
  std::uint64_t hashSeed2;
  std::uint64_t hashModulus;
  std::uint64_t numBits; // nota bene: numBits <= hashModulus
  std::uint64_t numUbits;
  std::uint64_t numIbits;
  bool          intersect;
  std::uint32_t compressor;
  bool outputStats;
  std::string statsFilename;
};

#endif