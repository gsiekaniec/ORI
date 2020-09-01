// file_manager.cc-- manage associations between bloom filters and files,
//                   including which filers are resident in memory

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <chrono>

#include "utilities.h"
#include "bloom_filter.h"
#include "bloom_tree.h"
#include "file_manager.h"

using std::string;
using std::vector;
using std::pair;
using std::cerr;
using std::endl;
#define u32 std::uint32_t

//----------
//
// initialize class variables
//
//----------

bool           FileManager::trackMemory     = false;
bool           FileManager::dbgContentLoad  = false;

bool           FileManager::reportOpenClose = false;
string         FileManager::openedFilename  = "";
std::ifstream* FileManager::openedFile      = nullptr;

//----------
//
// FileManager--
//
//----------

FileManager::FileManager
   (BloomTree*	root,
	bool		validateConsistency)
	  :	modelBf(nullptr)
	{
	// scan the tree, setting up hashes from (a) node name to a node's filter
	// filename and from (b) filter filename to names of nodes within; we also
	// verify that the node names are distinct

	std::unordered_map<string,string> nameToFile;

	vector<BloomTree*> order;
	root->post_order(order);

	for (const auto& node : order)
		node->manager = this;

	for (const auto& node : order)
		{
		if (nameToFile.count(node->name) > 0)
			fatal ("error: tree contains more than one node named \""
			     + node->name 
			     + " (in \"" + node->bfFilename + "\""
			     + " and \"" + nameToFile[node->name] + "\")");

		nameToFile[node->name] = node->bfFilename;
		nameToNode[node->name] = node;

		if (filenameToNames.count(node->bfFilename) == 0)
			{
			filenameToNames[node->bfFilename] = new vector<string>;
			if (trackMemory)
				cerr << "@+" << filenameToNames[node->bfFilename]
				     << " allocating names vector for FileManager.filenameToNames[" << node->bfFilename << "]" << endl;
			alreadyPreloaded[node->bfFilename] = false;
			}
		filenameToNames[node->bfFilename]->emplace_back(node->name);
		}

	// preload the content headers for every node, file-by-file; this has
	// two side effects -- (1) the bloom filter properties are checked for
	// consistency, and (2) we are installed as every bloom filter's manager

	if (validateConsistency)
		{
		for (auto iter : filenameToNames) 
			{
			string filename = iter.first;
			preload_content (filename);
			}
		}

	if (trackMemory)
		cerr << "@+" << this << " constructor FileManager" << endl;
	}

FileManager::~FileManager()
	{
	if (trackMemory)
		cerr << "@-" << this << " destructor FileManager" << endl;

	if (modelBf != nullptr) delete modelBf;
	for (auto iter : filenameToNames) 
		{
		vector<string>* names = iter.second;
		if (trackMemory)
			{
			string bfFilename = iter.first;
			cerr << "@-" << names
			     << " discarding names vector for FileManager.filenameToNames[" << bfFilename << "]" << endl;
			}
		delete names;
		}
	}

void FileManager::preload_content
   (const string&	filename)
	{
	wall_time_ty startTime;

	if (filenameToNames.count(filename) == 0)
		fatal ("internal error: attempt to preload content from"
		       " unknown file \"" + filename + "\"");

	if (alreadyPreloaded[filename]) return;

	// $$$ add trackMemory for in
	if (BloomFilter::reportLoadTime || BloomFilter::reportTotalLoadTime) startTime = get_wall_time();
	std::ifstream* in = FileManager::open_file(filename,std::ios::binary|std::ios::in,
	                                           /* positionAtStart*/ true);
	if (not *in)
		fatal ("error: FileManager::preload_content()"
			   " failed to open \"" + filename + "\"");
	if (BloomFilter::reportLoadTime || BloomFilter::reportTotalLoadTime)
		{
		double elapsedTime = elapsed_wall_time(startTime);
		if (BloomFilter::reportLoadTime)
			cerr << "[BloomFilter open] " << std::setprecision(6) << std::fixed << elapsedTime << " secs " << filename << endl;
		if (BloomFilter::reportTotalLoadTime)
			BloomFilter::totalLoadTime += elapsedTime;  // $$$ danger of precision error?
		}

	vector<pair<string,BloomFilter*>> content
		= BloomFilter::identify_content(*in,filename);

	if (dbgContentLoad)
		{
		cerr << "FileManager::preload_content, \"" << filename << "\" contains" << endl;
		for (const auto& templatePair : content)
			{
			string       bfName     = templatePair.first;
			BloomFilter* templateBf = templatePair.second;
			cerr << "  \"" << bfName << "\""
			     << " templateBf=" << templateBf << endl;
			}
		}

	vector<string>* nodeNames = filenameToNames[filename];
	if (content.size() != nodeNames->size())
		fatal ("error: \"" + filename + "\""
		     + " contains " + std::to_string(content.size()) + " bloom filters"
		     + ", in conflict with the tree's topology"
		     + " (expected " + std::to_string(nodeNames->size()) + ")");

	for (const auto& templatePair : content)
		{
		string       bfName     = templatePair.first;
		BloomFilter* templateBf = templatePair.second;
		if (not contains (*nodeNames, bfName))
			fatal ("error: \"" + filename + "\""
			     + " contains the bloom filter \"" + bfName + "\""
			     + ", in conflict with the tree's topology");

		BloomTree* node = nameToNode[bfName];
		if (dbgContentLoad)
			cerr << "FileManager::preload_content (\"" << bfName << "\")"
			     << " node=" << node
			     << " node->name=" << node->name
			     << " node->bf=" << node->bf << endl;

		// if the node has already been loaded, leave it be

		if ((node->bf != nullptr) and (node->bf->ready))
			continue;

		// copy the template into the node's filter

		if (node->bf == nullptr)
			{
			if (dbgContentLoad)
				cerr << "  creating new BF for \"" << bfName << "\"" << endl;
			node->bf = BloomFilter::bloom_filter(templateBf);
			node->bf->manager = this;
			}
		else // node exists but is not ready
			{
			if (dbgContentLoad)
				cerr << "  using existing (but not ready) BF for \"" << bfName << "\"" << endl;
			node->bf->copy_properties(templateBf);
			node->bf->setSizeKnown = templateBf->setSizeKnown;
			node->bf->setSize      = templateBf->setSize;
			}

		node->bf->steal_bits(templateBf);
		delete templateBf;

		// make sure all bloom filters in the tree are consistent

		if (modelBf == nullptr)
			modelBf = BloomFilter::bloom_filter(node->bf);
		else
			node->bf->is_consistent_with(modelBf,/*beFatal*/true);
		}

	alreadyPreloaded[filename] = true;

	// $$$ add trackMemory for in
	FileManager::close_file(in);
	}

void FileManager::load_content
   (const string&	filename,
	const string&	_whichNodeName)
	{
	// ……… when we implement a heap, and dirty bits, we'll need to empty the heap here

	if (filenameToNames.count(filename) == 0)
		fatal ("internal error: attempt to load content from"
		       " unknown file \"" + filename + "\"");

	string whichNodeName = _whichNodeName;
	if (not alreadyPreloaded[filename])
		{
		preload_content(filename);
		whichNodeName = "";  // we will load all nodes in this file
		}

//øøø we only need to load this if it hasn't already been loaded

	vector<string>* nodeNames = filenameToNames[filename];
	for (const auto& nodeName : *nodeNames)
		{
		if ((whichNodeName != "") and (nodeName != whichNodeName))
			continue;
		if (dbgContentLoad)
			cerr << "FileManager::load_content nodeName = \"" << nodeName << "\"" << endl;
		BloomTree* node = nameToNode[nodeName];
		node->bf->reportLoad = reportLoad;
		node->bf->load(/*bypassManager*/ true);
		}

	}

//----------
//
// open_file, close_file--
//	Wrapper for input stream open and close, keeping a file open until some
//	other file is needed.  This (hopefully) saves us the overhead of opening a
//	file, closing it, then opening it again.
//
//----------
//
// CAVEAT:	Any command that reads bit vectors needs to call close_file()
//			before exit, to avoid a memory leak.
//
//----------

std::ifstream* FileManager::open_file
   (const string&			filename,
	std::ios_base::openmode	mode, 
	bool                    positionAtStart)
	{
	if ((openedFile != nullptr) && (filename == openedFilename))
		{
		if (positionAtStart)
			openedFile->seekg(0,openedFile->beg);
		return openedFile;
		}

	if (openedFile != nullptr)
		{
		if (reportOpenClose)
			cerr << "[FileManager close_file] " << openedFilename << endl;
		openedFile->close();
		delete openedFile;
		}

	if (reportOpenClose)
			cerr << "[FileManager open_file] " << filename << endl;
	openedFilename = filename;
	openedFile     = new std::ifstream(filename,mode);
	return openedFile;
	}

void FileManager::close_file
   (std::ifstream*	in,
	bool			really)
	{
	if (openedFile == nullptr) return; // $$$ should this be an error?

	if (in == nullptr)
		{ in = openedFile;  really = true; }
	else if (in != openedFile)
		fatal ("error: FileManager::close_file()"
			   " is asked to close the wrong file");

	if (really)
		{
		if (reportOpenClose)
			cerr << "[FileManager close_file] " << openedFilename << endl;
		openedFile->close();
		delete openedFile;
		openedFilename = "";
		openedFile     = nullptr;
		}
	}

