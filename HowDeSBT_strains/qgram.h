#ifndef qgram_H
#define qgram_H

#include "kseq.h"
#include <stdio.h>
#include <zlib.h>
#include <vector>
#include <iostream>
#include <fstream>

KSEQ_INIT(gzFile, gzread)

//inline std::vector<std::string> get_seeds(std::string &in_seeds)
//{
//  std::vector<std::string> seeds;
//  std::string line;
//  std::ifstream ifs(in_seeds);
//  if (! ifs)
//    fatal("Failed to open " + in_seeds);
//
//  while ( getline(ifs, line) )
//    seeds.push_back(line);
//  return seeds;
//}

inline std::vector<std::string> get_seeds(std::string &in_seeds)
{
  std::vector<std::string> seeds;
  std::string line;
  std::ifstream ifs(in_seeds);
  if (! ifs)
    fatal("Failed to open " + in_seeds);

  while ( getline(ifs, line) )
  {
    if (line.find('*') != std::string::npos)
    {
      fatal("Given seeds must be in the form 0/1, found: " + line);
    }
    else
    {
      seeds.push_back(line);
    }
  }
  return seeds;
}

inline std::vector<std::vector<int8_t>> get_pos_from_seeds(std::vector<std::string> &seeds)
{
  std::vector<std::vector<int8_t>> false_positions;

  for (auto& it : seeds)
  {
    std::vector<int8_t> positions;
    int8_t p = 0;
    for (char& c : it)
    {
      if ( !(c & 0x1) ) positions.push_back(p);
      p++;
    }
    false_positions.push_back(positions);
  }
  return false_positions;
}

inline std::string reverse_comp(std::string kmer)
{
  std::reverse(kmer.begin(), kmer.end());
  for (std::size_t i = 0; i < kmer.length(); ++i){
        switch (kmer[i]){
        case 'A':
            kmer[i] = 'T';
            break;
        case 'C':
            kmer[i] = 'G';
            break;
        case 'G':
            kmer[i] = 'C';
            break;
        case 'T':
            kmer[i] = 'A';
            break;
        }
    }
    return kmer;
}

inline std::string cano(std::string k)
{
  std::string comp = reverse_comp(k);
  if (strcmp(k.c_str(), comp.c_str()) > 0)
  {
    return comp;
  }
  return k;
}

inline std::string apply_seed(std::string kmer, std::vector<int8_t> &pos)
{
  int i = 0;
  for (auto& p : pos)
  {
    kmer.erase(p-i, 1);
    i++;
  }
  return kmer;
}

#endif

