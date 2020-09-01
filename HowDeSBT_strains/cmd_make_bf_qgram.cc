#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <vector>

#include "utilities.h"
#include "hash.h"
//#include "jelly_kmers.h"
#include "bloom_filter.h"
#include "bloom_filter_file.h"

#include "support.h"
#include "commands.h"
#include "qgram.h"
#include "cmd_make_bf_qgram.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

#define u32 std::uint32_t
#define u64 std::uint64_t

void MakeBFQCommand::short_description(std::ostream &s)
{
  s << commandName << "-- convert a sequence file to a bloom filter using qgrams" << endl;
}

void MakeBFQCommand::usage(std::ostream &s,
                           const string &message)
{
  if (!message.empty())
  {
    s << message << endl;
    s << endl;
  }

  short_description(s);
  s << "usage: " << commandName << " <filename> [<filename>..] [options]" << endl;
  s << "  <filename>         (cumulative) a sequence file, e.g. fasta, fastq (support gzip)" << endl;
  s << "  --qgram=<filename> file containing one seed per line" << endl;
  s << "  --out=<filename>   name for bloom filter file" << endl;
  s << "                     (by default this is derived from first sequence filename)" << endl;
  s << "  --list=<filename>  file containing a list of files; one bloom" << endl;
  s << "                     is created for the union/intersection of sequence files; " << endl;
  s << "                     first line corresponds bf filename without suffix; " << endl;
  s << "                     others lines are path to fasta/q/gz file" << endl;
  s << "  --k=<N>            kmer size (number of nucleotides in a kmer)" << endl;
  s << "                     (default is " << defaultKmerSize << ")" << endl;
  s << "  --hashes=<N>       how many hash functions to use for the filter" << endl;
  s << "                     (default is " << defaultNumHashes << ") WARNING: Currently support only one-hash BF" << endl;
  s << "  --seed=<number>    the hash function's 56-bit seed" << endl;
  s << "  --seed=<number>,<number>  both the hash function seeds; the second seed is" << endl;
  s << "                     only used if more than one hash function is being used" << endl;
  s << "                     (by default the second seed is the first seed plus 1)" << endl;
  s << "  --modulus=<M>      set the hash modulus, if larger than the number of bits" << endl;
  s << "                     (by default this is the same as the number of bits)" << endl;
  s << "  --bits=<N>         number of bits in the bloom filter" << endl;
  s << "                     (default is " << defaultNumBits << ")" << endl;
  s << "  --ubits=<N>        number of bits in the bloom filter that corresponds" << endl;
  s << "                     to the union of all fasta in --list" << endl;
  s << "  --ibits=<N>        number of bits in the bloom filter that corresponds" << endl;
  s << "                     to the union of all fasta in --list" << endl;
  s << "  --uncompressed     make the filter with uncompressed bit vector(s)" << endl;
  s << "                     (this is the default)" << endl;
  s << "  --rrr              make the filter with RRR-compressed bit vector(s)" << endl;
  s << "  --roar             make the filter with roar-compressed bit vector(s)" << endl;
  s << "  --stats[=<filename>] write bloom filter stats to a text file" << endl;
  s << "                     (if no filename is given this is derived from the bloom" << endl;
  s << "                     filter filename)" << endl;
  s << endl;
  s << "When --list is used, each line of the file corresponds to a bloom filter. The" << endl;
  s << "format of each line is" << endl;
  s << "  <filename> [<filename>..] [--kmersin] [--out=<filename>]" << endl;
  s << "with meaning the same as on the command line. No other options (e.g. --k or" << endl;
  s << "--bits) are allowed in the file. These are specified on the command line and" << endl;
  s << "will affect all the bloom filters." << endl;
  s << endl;
}

void MakeBFQCommand::debug_help(std::ostream &s)
{
  s << "--debug= options" << endl;
  s << "  settings" << endl;
  s << "  add" << endl;
  s << "  contains" << endl;
  s << "  kmers" << endl;
  s << "  strings" << endl;
  s << "  v1file" << endl;
}

void MakeBFQCommand::parse(int _argc,
                           char **_argv)
{
  int argc;
  char **argv;

  // defaults

  listFilename = "";
  bfFilename = "";
  kmerSize = defaultKmerSize;
  hashSeed1 = defaultHashSeed;
  hashSeed2 = defaultHashSeed;
  numBits = defaultNumBits;
  numUbits = defaultUNumBits;
  numIbits = defaultUNumBits;
  hashModulus = 0;
  bool hashModulusSet = false;
  compressor = bvcomp_uncompressed;
  outputStats = false;
  statsFilename = "";
  intersect = false;

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

    // --out=<filename>

    if ((is_prefix_of(arg, "--out=")) || (is_prefix_of(arg, "--output=")))
    {
      bfFilename = argVal;
      continue;
    }

    // --qgram=<filename>
    if (is_prefix_of(arg, "--qgram="))
    {
      seed_file = argVal;
      continue;
    }

    // --list=<filename>

    if (is_prefix_of(arg, "--list="))
    {
      listFilename = argVal;
      continue;
    }

    // --k=<N>

    if ((is_prefix_of(arg, "K=")) || (is_prefix_of(arg, "--K=")) || (is_prefix_of(arg, "k=")) || (is_prefix_of(arg, "--k=")) || (is_prefix_of(arg, "--kmer=")) || (is_prefix_of(arg, "--kmersize=")))
    {
      kmerSize = string_to_u32(argVal);
      continue;
    }

    // --modulus=<M>

    if ((is_prefix_of(arg, "--modulus=")) || (is_prefix_of(arg, "M=")) || (is_prefix_of(arg, "--M=")))
    {
      hashModulus = string_to_unitized_u64(argVal);
      hashModulusSet = true;
      continue;
    }

    // --bits=<N>

    if ((is_prefix_of(arg, "--bits=")) || (is_prefix_of(arg, "B=")) || (is_prefix_of(arg, "--B=")))
    {
      numBits = string_to_unitized_u64(argVal);
      continue;
    }

    if (is_prefix_of(arg, "--ubits="))
    {
      numUbits = string_to_unitized_u64(argVal);
      continue;
    }

    if (is_prefix_of(arg, "--ibits="))
    {
      numIbits = string_to_unitized_u64(argVal);
      intersect = true;
      continue;
    }
    // bit vector type

    if (arg == "--uncompressed")
    {
      compressor = bvcomp_uncompressed;
      continue;
    }

    if ((arg == "--rrr") || (arg == "--RRR"))
    {
      compressor = bvcomp_rrr;
      continue;
    }

    if ((arg == "--roar") || (arg == "--roaring"))
    {
      compressor = bvcomp_roar;
      continue;
    }

    // (unadvertised) special bit vector types

    if ((arg == "--zeros") || (arg == "--allzeros") || (arg == "--all_zeros") || (arg == "--all-zeros"))
    {
      compressor = bvcomp_zeros;
      continue;
    }

    if ((arg == "--ones") || (arg == "--allones") || (arg == "--all_ones") || (arg == "--all-ones"))
    {
      compressor = bvcomp_ones;
      continue;
    }

    if (arg == "--uncrrr")
    {
      compressor = bvcomp_unc_rrr;
      continue;
    }

    if (arg == "--uncroar")
    {
      compressor = bvcomp_unc_roar;
      continue;
    }

    // --stats[=<filename>]

    if (arg == "--stats")
    {
      outputStats = true;
      continue;
    }

    if ((is_prefix_of(arg, "--stats=")))
    {
      outputStats = true;
      statsFilename = argVal;
      continue;
    }

    // (unadvertised) debug options

    if (arg == "--debug")
    {
      debug.insert("debug");
      continue;
    }

    if (is_prefix_of(arg, "--debug="))
    {
      for (const auto &field : parse_comma_list(argVal))
        debug.insert(to_lower(field));
      continue;
    }

    // unrecognized --option

    if (is_prefix_of(arg, "--"))
      chastise("unrecognized option: \"" + arg + "\"");

    // <filename>

    seqFilenames.emplace_back(strip_blank_ends(arg));
  }

  // sanity checks

  if ((compressor == bvcomp_zeros) || (compressor == bvcomp_ones))
  {
    if (not listFilename.empty())
      chastise("cannot use --list with --zeros or --ones");
    if (not seqFilenames.empty())
      chastise("cannot use sequence files with --zeros or --ones");
    if (bfFilename.empty())
      chastise("--zeros or --ones requires --out");
  }
  else if (listFilename.empty())
  {
    if (seqFilenames.empty())
      chastise("at least one sequence filename is required");
  }
  else
  {
    if (not seqFilenames.empty())
      chastise("cannot use --list with sequence filenames (e.g. " + seqFilenames[0] + ") in the command");

    if (not bfFilename.empty())
      chastise("cannot use --list with a filter filename (" + bfFilename + ") in the command");
    if (not statsFilename.empty())
      chastise("cannot use --list with a stats filename (" + statsFilename + ") in the command");
  }

  if (kmerSize == 0)
    chastise("kmer size cannot be zero");

  if (numBits < 2)
    chastise("number of bits must be at least 2");

  if (not hashModulusSet)
    hashModulus = numBits;
  if (hashModulus < numBits)
    chastise("hash modulus (" + std::to_string(hashModulus) + ")" + " cannot be less than the number of bits" + " (" + std::to_string(numBits) + ")");

  if (contains(debug, "settings"))
  {
    cerr << "kmerSize    = " << kmerSize << endl;
    cerr << "hashSeed1   = " << hashSeed1 << endl;
    cerr << "hashModulus = " << hashModulus << endl;
    cerr << "numBits     = " << numBits << endl;
    cerr << "compressor  = " << compressor << endl;
  }

  return;
}

int MakeBFQCommand::execute()
{
  if (listFilename.empty())
  {
    make_bloom_filter_qgram();
  }

  else
  {
    std::ifstream in(listFilename);
    if (not in)
      fatal("error: failed to open \"" + listFilename + "\"");

    string line;
    int lineNum = 0;
    while (std::getline(in, line))
    {
      if (lineNum == 0)
      {
        bfFilename = line + ".bf";
        lineNum++;
        continue;
      }
      lineNum++;
      vector<string> tokens = tokenize(line);
      for (size_t argIx = 0; argIx < tokens.size(); argIx++)
      {
        string arg = tokens[argIx];
        string argVal;

        string::size_type argValIx = arg.find('=');
        if (argValIx == string::npos)
          argVal = "";
        else
          argVal = arg.substr(argValIx + 1);

        if ((is_prefix_of(arg, "--out=")) || (is_prefix_of(arg, "--output=")))
        {
          bfFilename = argVal;
          continue;
        }

        if (is_prefix_of(arg, "--"))
          fatal("unrecognized field: \"" + arg + "\"" + " at line " + std::to_string(lineNum) + " in " + listFilename);

        seqFilenames.emplace_back(strip_blank_ends(arg));
      }
    }
    in.close();
    if (intersect)
      make_multi_bloom_filter_qgram_intersect();
    else
      make_multi_bloom_filter_qgram();
  }
  
  return EXIT_SUCCESS;
}

void MakeBFQCommand::make_bloom_filter_qgram()
{
  int numHashes = 1;

  vector<string> seeds = get_seeds(seed_file);
  vector<Hash> hashers;
  int q_size = 0;
  for (size_t i = 0; i < seeds.size(); i++)
  {
    q_size = std::count(seeds[i].begin(), seeds[i].end(), '1');
    hashers.push_back(Hash(q_size, hashSeed1));
  }
  vector<vector<int8_t>> false_positions = get_pos_from_seeds(seeds);
  char *tmp = (char *)malloc(kmerSize + 1);
  tmp[kmerSize] = '\0';
  string kmer;
  
  for (const auto& input_file: seqFilenames)
  {
    u64 qgramAdded = 0;
    BloomFilter *bfi = new BloomFilter(input_file+".bf", kmerSize, numHashes, hashSeed1, hashSeed2, numBits, hashModulus);
    bfi->new_bits(compressor);
    gzFile f_in = gzopen(input_file.c_str(), "r");
    kseq_t *ks = kseq_init(f_in);
    while (kseq_read(ks) >= 0)
    {
      for (size_t i = 0; i < ks->seq.l - kmerSize + 1; i++)
      {
        memcpy(tmp, &ks->seq.s[i], kmerSize);
        kmer = tmp;

        int h = 0;
        for (auto &pos : false_positions)
        {
          bfi->add_qgram(apply_seed(cano(kmer), pos), hashers[h]);
          qgramAdded++;
          h++;
        }
      }
    }
    gzclose(f_in);
    kseq_destroy(ks);
    if (not contains(debug, "v1file"))
    {
      bfi->setSizeKnown = true;
      bfi->setSize = qgramAdded;
    }

    if ((compressor == bvcomp_unc_rrr) || (compressor == bvcomp_unc_roar))
    {
      BitVector *bvi = bfi->bvs[0];
      bvi->unfinished();
    }

    bfi->reportSave = true;
    bfi->save();
    delete bfi;
  }
}

void MakeBFQCommand::make_multi_bloom_filter_qgram()
{
  int numHashes = 1;
  //string bfOutFilename = build_output_filename();
  cout << numBits << endl;
  cout << numUbits << endl;
  u64 bhashModulus = numUbits;
  BloomFilter *bf = new BloomFilter(bfFilename, kmerSize, numHashes, hashSeed1, hashSeed2, numUbits, bhashModulus);
  bf->new_bits(compressor);
  u64 qgramAddedAll = 0;
  //string input_file = seqFilenames[0];

  vector<string> seeds = get_seeds(seed_file);
  vector<Hash> hashers;
  int q_size = 0;
  for (size_t i = 0; i < seeds.size(); i++)
  {
    q_size = std::count(seeds[i].begin(), seeds[i].end(), '1');
    hashers.push_back(Hash(q_size, hashSeed1));
  }
  vector<vector<int8_t>> false_positions = get_pos_from_seeds(seeds);
  char *tmp = (char *)malloc(kmerSize + 1);
  tmp[kmerSize] = '\0';
  string kmer;
  string to_add;
  for (const auto& input_file: seqFilenames)
  {
    u64 qgramAdded = 0;
    BloomFilter *bfi = new BloomFilter(input_file+".bf", kmerSize, numHashes, hashSeed1, hashSeed2, numBits, hashModulus);
    bfi->new_bits(compressor);
    gzFile f_in = gzopen(input_file.c_str(), "r");
    kseq_t *ks = kseq_init(f_in);
    while (kseq_read(ks) >= 0)
    {
      for (size_t i = 0; i < ks->seq.l - kmerSize + 1; i++)
      {
        memcpy(tmp, &ks->seq.s[i], kmerSize);
        kmer = tmp;

        int h = 0;
        for (auto &pos : false_positions)
        {
          to_add = apply_seed(cano(kmer), pos);
          bf->add_qgram(to_add, hashers[h]);
          bfi->add_qgram(to_add, hashers[h]);
          qgramAdded++;
          h++;
        }
      }
    }
    gzclose(f_in);
    kseq_destroy(ks);
    qgramAddedAll += qgramAdded;
    if (not contains(debug, "v1file"))
    {
      bfi->setSizeKnown = true;
      bfi->setSize = qgramAdded;
    }

    if ((compressor == bvcomp_unc_rrr) || (compressor == bvcomp_unc_roar))
    {
      BitVector *bvi = bfi->bvs[0];
      bvi->unfinished();
    }

    bfi->reportSave = true;
    bfi->save();
    delete bfi;
  }

  if (not contains(debug, "v1file"))
  {
    bf->setSizeKnown = true;
    bf->setSize = qgramAddedAll;
  }

  if ((compressor == bvcomp_unc_rrr) || (compressor == bvcomp_unc_roar))
  {
    BitVector *bv = bf->bvs[0];
    bv->unfinished();
  }

  bf->reportSave = true;
  bf->save();
  delete bf;
}

void MakeBFQCommand::make_multi_bloom_filter_qgram_intersect()
{
  int numHashes = 1;
  //string bfOutFilename = build_output_filename();

  u64 qgramAddedAll = 0;
  //string input_file = seqFilenames[0];

  vector<string> seeds = get_seeds(seed_file);
  vector<Hash> hashers;
  int q_size = 0;
  for (size_t i = 0; i < seeds.size(); i++)
  {
    q_size = std::count(seeds[i].begin(), seeds[i].end(), '1');
    hashers.push_back(Hash(q_size, hashSeed1));
  }
  vector<vector<int8_t>> false_positions = get_pos_from_seeds(seeds);
  char *tmp = (char *)malloc(kmerSize + 1);
  tmp[kmerSize] = '\0';
  string kmer;
  string to_add;

  u64 bhashModulus = numIbits;

  BloomFilter *bf_intersect = new BloomFilter(bfFilename, kmerSize, numHashes, hashSeed1, hashSeed2, numIbits, bhashModulus);
  bf_intersect->new_bits(compressor);


  gzFile f_in = gzopen(seqFilenames[0].c_str(), "r");
  kseq_t *ks = kseq_init(f_in);
  while (kseq_read(ks) >= 0)
  {
    for (size_t i = 0; i < ks->seq.l - kmerSize + 1; i++)
    {
      memcpy(tmp, &ks->seq.s[i], kmerSize);
      kmer = tmp;

      int h = 0;
      for (auto &pos : false_positions)
      {
        to_add = apply_seed(cano(kmer), pos);
        bf_intersect->add_qgram(to_add, hashers[h]);
        qgramAddedAll++;
        h++;
      }
    }
  }
  gzclose(f_in);
  kseq_destroy(ks);

  int first = true;
  for (const auto& input_file: seqFilenames)
  {
    if (first)
    {
      first = false;
      continue;
    }
    u64 qgramAdded = 0;
    BloomFilter *bfi = new BloomFilter(input_file+".bf", kmerSize, numHashes, hashSeed1, hashSeed2, numIbits, bhashModulus);
    bfi->new_bits(compressor);
    
    gzFile f_in = gzopen(input_file.c_str(), "r");
    kseq_t *ks = kseq_init(f_in);
    while (kseq_read(ks) >= 0)
    {
      for (size_t i = 0; i < ks->seq.l - kmerSize + 1; i++)
      {
        memcpy(tmp, &ks->seq.s[i], kmerSize);
        kmer = tmp;

        int h = 0;
        for (auto &pos : false_positions)
        {
          to_add = apply_seed(cano(kmer), pos);
          bfi->add_qgram(to_add, hashers[h]);
          qgramAdded++;
          h++;
        }
      }
    }
    bf_intersect->intersect_with(bfi->bvs[0]);
    gzclose(f_in);
    kseq_destroy(ks);
    qgramAddedAll += qgramAdded;
    if (not contains(debug, "v1file"))
    {
      bfi->setSizeKnown = true;
      bfi->setSize = qgramAdded;
    }

    if ((compressor == bvcomp_unc_rrr) || (compressor == bvcomp_unc_roar))
    {
      BitVector *bvi = bfi->bvs[0];
      bvi->unfinished();
    }
    delete bfi;
  }

  if (not contains(debug, "v1file"))
  {
    bf_intersect->setSizeKnown = true;
    bf_intersect->setSize = qgramAddedAll;
  }

  if ((compressor == bvcomp_unc_rrr) || (compressor == bvcomp_unc_roar))
  {
    BitVector *bv = bf_intersect->bvs[0];
    bv->unfinished();
  }

  bf_intersect->reportSave = true;
  bf_intersect->save();
  delete bf_intersect;
}

void MakeBFQCommand::report_stats(const string &bfOutFilename, u64 qgramsAdded)
{
  int numHashes = 1;
  double fpRate = BloomFilter::false_positive_rate(numHashes, numBits, qgramsAdded);

  if (outputStats)
  {
    string statsOutFilename = build_stats_filename(bfOutFilename);
    cerr << "writing bloom filter stats to \"" << statsOutFilename << "\"" << endl;

    std::ofstream statsF(statsOutFilename);

    statsF << "#filename"
           << "\tnumHashes"
           << "\tnumBits"
           << "\tqgramsAdded"
           << "\tbfFpRate"
           << endl;
    statsF << bfOutFilename
           << "\t" << numHashes
           << "\t" << numBits
           << "\t" << qgramsAdded
           << "\t" << fpRate
           << endl;
  }

  if (contains(debug, "fprate"))
  {
    cerr << bfOutFilename << " qgrams inserted: " << qgramsAdded << endl;
    cerr << bfOutFilename << " estimated BF false positive rate: " << fpRate << endl;
  }
}

string MakeBFQCommand::build_output_filename()
{
  string bfOutFilename = bfFilename;

  if (bfOutFilename.empty())
  {
    string ext = "." + BitVector::compressor_to_string(compressor) + ".bf";
    if (ext == ".uncompressed.bf")
      ext = ".bf";

    string seqFilename = seqFilenames[0];
    string::size_type dotIx = seqFilename.find_last_of(".");
    if (dotIx == string::npos)
      bfOutFilename = seqFilename + ext;
    else
      bfOutFilename = seqFilename.substr(0, dotIx) + ext;
  }

  return bfOutFilename;
}

string MakeBFQCommand::build_stats_filename(const string &bfOutFilename)
{
  string statsOutFilename = statsFilename;

  if (statsOutFilename.empty())
  {
    string ext = ".stats";

    string::size_type dotIx = bfOutFilename.find_last_of(".");
    if (dotIx == string::npos)
      statsOutFilename = bfOutFilename + ext;
    else
      statsOutFilename = bfOutFilename.substr(0, dotIx) + ext;
  }

  return statsOutFilename;
}