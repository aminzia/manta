// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Ole Schulz-Trieglaff
///


#include "assembly/SmallAssembler.hh"

#include "boost/foreach.hpp"
#include "boost/unordered_map.hpp"

#include <cassert>
#include <vector>
#include <fstream>
#include <limits>

// compile with this macro to get verbose output:
//#define DEBUG_ASBL


// stream used by DEBUG_ASBL:
#ifdef DEBUG_ASBL
#include "blt_util/log.hh"
#include <iostream>
#endif

// TODO: Some ideas:
// * Store list of reference/read end k-mer. Enforce that the contig finishes with an end k-mer.
// * Compute average coverage per contig.


// maps kmers to positions in read
typedef boost::unordered_map<std::string,unsigned> str_uint_map_t;
// keeps track of k-mers
typedef boost::unordered_map<std::string,bool> str_bool_map_t;

static const std::string alphabet("ACGT");

static const unsigned MIN_KMER_FREQ = 1;


/**
 * Adds base @p base to the end (isEnd is true) or start (otherwise) of the contig.
 *
 *	@return The extended contig.
 */
static
std::string
addBase(const std::string& contig,
        const char base,
        const bool isEnd)
{
    return isEnd ? (contig + base) : (base + contig);
}



/**
 * Returns a suffix (isEnd is true) or prefix (otherwise) of @p contig with length @p length.
 *
 *	@return The suffix or prefix.
 */
static
std::string
getEnd(const std::string& contig,
       const unsigned length,
       const bool isEnd)
{

    const unsigned csize(contig.size());
    /*if(length <= csize) {
        std::cout << "contig=" << contig << std::endl;
    }*/
    assert(length <= csize);

    if (isEnd) return contig.substr((csize-length),length);
    else       return contig.substr(0,length);
}

#ifdef DEBUG_ASBL
static
void
dumpHash(const str_uint_map_t& wordCount,
         const unsigned wordLength,
         const unsigned iteration) {

	std::ofstream outDotFile;
	std::ofstream outTxtFile;
    std::stringstream sstr;
    sstr << "debruijn.graph.k";
    sstr << wordLength;
    sstr << ".iter";
    sstr << iteration;

    std::string outFileStem(sstr.str().c_str());
    std::string dotFileName = outFileStem + ".dot";
    std::string txtFileName = outFileStem + ".txt";

	outDotFile.open(dotFileName.c_str());
	outTxtFile.open(txtFileName.c_str());

	// kmer nodes with coverage higher than this get a different color
	//static const int lowCovGraphVisThreshold(3);
	static const unsigned lowCovGraphVisThreshold(MIN_KMER_FREQ);
	static const std::string lowCovNodeColor("red");
	static const std::string highCovNodeColor("green");

	outDotFile << "graph {\n";
	outDotFile << "node [ style = filled ];\n";
    str_uint_map_t aliasH;
    unsigned n(0);
	for (str_uint_map_t::const_iterator ct = wordCount.begin();ct!=wordCount.end();++ct) {
        aliasH[ct->first] = n++;
        //graphviz out
		outDotFile << aliasH[ct->first] << "[label=\"cov" << ct->second << "\" ";
		if (ct->second>lowCovGraphVisThreshold)
		{
			outDotFile << "color=" << highCovNodeColor << "]";
		} else {
			outDotFile << "color=" << lowCovNodeColor << "]";
		}
		outDotFile << "\n";

		// txt output
		outTxtFile << "VT " << aliasH[ct->first] << " " << ct->first << " " << ct->second << "\n"; 
	}
    // need to add edges here
    static const bool isEnd(true);
	for (str_uint_map_t::const_iterator ct = wordCount.begin();ct!=wordCount.end();++ct) {
        std::string tmp(getEnd(ct->first,wordLength-1,isEnd));
        BOOST_FOREACH(const char symbol, alphabet) {
            const std::string newKey(addBase(tmp,symbol,isEnd));
            if (wordCount.find(newKey) != wordCount.end()) {
                outDotFile << aliasH[ct->first] << " -- " << aliasH[newKey] << ";\n";
                //outTxtFile << "ED " << aliasH[ct->first] << " " << aliasH[newKey] << "\n";
                outTxtFile << "ED " << ct->first << " " << newKey << "\n";
            }
        }
	}
	outDotFile << "}\n";

	outDotFile.close();
	outTxtFile.close();
}
#endif

#if 0
static
void
findMaxPaths(const str_uint_map_t& wordCount,
			 Contig& contigSoFar,
			 const unsigned wordLength,
			 str_bool_map_t& seenVertices,
			 std::vector<Contig>& contigs) {

	// we walk only to the left
	static const bool isEnd(true);

	std::string tmp(getEnd(contigSoFar.seq,wordLength-1,isEnd));
	bool neighbourFound(false);
	BOOST_FOREACH(const char symbol, alphabet) {
		const std::string overlap(addBase(tmp,symbol,isEnd));

        //TODO:: this follows all branches regardless of coverage
        // need to see if this is a good thing
		str_uint_map_t::const_iterator ct = wordCount.find(overlap);
	    if (ct != wordCount.end() && !seenVertices[overlap]) {
	        seenVertices[overlap] = true;
	        // copy contig sequence
	    	Contig newCtg(contigSoFar);
	     	newCtg.seq += symbol;
	     	// update running average of coverage
	     	newCtg.avgCoverage = (ct->second+newCtg.numKmers*newCtg.avgCoverage)/(newCtg.numKmers+1);
	     	++newCtg.numKmers;

	    	//std::cerr << "Extending contig " << symbol << " " << newCtg << "\n";
	    	findMaxPath(wordCount,newCtg,wordLength,seenVertices,contigs);
	    	neighbourFound=true;
	    }
	}
	if (neighbourFound) {
		// at the end of a maximum path, save contig
#ifdef DEBUG_ASBL
		std::cerr << "Max path end: ctg=" << contigSoFar.seq << " cov=" << contigSoFar.avgCoverage << " numKmers=" << contigSoFar.numKmers << std::endl;
#endif
		contigs.push_back(contigSoFar);
		contigSoFar.seq="";
		contigSoFar.avgCoverage=0;
		contigSoFar.numKmers=0;
	}
}
#endif


static
void
doDFS(const str_uint_map_t& wordCount,
	  const Contig& contigSoFar,
	  const unsigned wordLength,
	  str_bool_map_t& seenVertices,
	  std::vector<Contig>& contigs) {

	// we walk only to the left
	static const bool isEnd(true);

	std::string tmp(getEnd(contigSoFar.seq,wordLength-1,isEnd));
	bool neighbourFound(false);
	BOOST_FOREACH(const char symbol, alphabet) {
		const std::string overlap(addBase(tmp,symbol,isEnd));
		//std::cerr << "Testing branch " << overlap << " " << symbol << "\n";
        //typedef std::pair<str_uint_map_t::const_iterator,str_uint_map_t::const_iterator> iterPair;
        //iterPair itP = wordCount.equal_range(overlap);
        //std::cerr << "Hits=" << distance(itP.first,itP.second) << "\n";

        //TODO:: this follows all branches regardless of coverage
        // need to see if this is a good thing
		str_uint_map_t::const_iterator ct = wordCount.find(overlap);
	    if (ct != wordCount.end() && !seenVertices[overlap]) {
	        seenVertices[overlap] = true;
	        // copy contig sequence
	    	Contig newCtg(contigSoFar);
	     	newCtg.seq += symbol;
	     	// update running average of coverage
	     	newCtg.avgCoverage = (ct->second+newCtg.numKmers*newCtg.avgCoverage)/(newCtg.numKmers+1);
	     	++newCtg.numKmers;

	    	//std::cerr << "Extending contig " << symbol << " " << newCtg << "\n";
	    	doDFS(wordCount,newCtg,wordLength,seenVertices,contigs);
	    	neighbourFound=true;
	    }
	}
	if (!neighbourFound) {
		// at the end of a branch, save contig
#ifdef DEBUG_ASBL
		std::cerr << "Branch end: ctg=" << contigSoFar.seq << " cov=" << contigSoFar.avgCoverage << " numKmers=" << contigSoFar.numKmers << " len=" << contigSoFar.seq.size() << std::endl;
#endif
		contigs.push_back(contigSoFar);
	}
}

static
void
pruneHash(str_uint_map_t& wordCount) 
{
    // prune hash, remove singleton k-mers
    unsigned pruneCnt(0);
    for (str_uint_map_t::iterator it=wordCount.begin();
         it!=wordCount.end();) {
        if(it->second <= MIN_KMER_FREQ)
        {
            it = wordCount.erase(it); 
            ++pruneCnt;
        } else {
            ++it;
        }
    }
    //std::cerr << "Pruned " << pruneCnt << " k-mers. Remaining " << wordCount.size() << std::endl;
}

static
void
initDFS(str_uint_map_t& wordCount,
		const std::string& startVertex,
		const unsigned wordLength,
		std::vector<Contig>& contigs
#ifdef DEBUG_ASBL
        ,const unsigned iteration
#endif
        )
         {

#ifdef DEBUG_ASBL
    std::cerr << "INIT DFS START " << startVertex << "\n";
#endif


#ifdef DEBUG_ASBL
    dumpHash(wordCount,wordLength,iteration);
#endif

	str_bool_map_t seenVertices;
	Contig ctg;
	ctg.seq = startVertex;
	ctg.avgCoverage = wordCount[startVertex];
    ctg.numKmers = 1;
	doDFS(wordCount, ctg, wordLength, seenVertices, contigs);
}

static
bool
buildHash(str_uint_map_t& wordHash,
		  std::string& firstWord,
		  std::vector<str_uint_map_t>& readWordOffsets,
		  const AssemblyReadOutput& readInfo,
		  const AssemblyReadInput& reads,
          const unsigned wordLength)
{
	int firstWordPos(std::numeric_limits<int>::max());
	const unsigned readCount(reads.size());

	for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
	{
		const AssemblyReadInfo& rinfo(readInfo[readIndex]);
        // skip reads used in a previous iteration
        if (rinfo.isUsed) continue;

        // stores the index of a kmer in a read sequence
	    const std::string& seq(reads[readIndex].second);
	    const unsigned readLen(seq.size());

	    // this read is unusable for assembly:
	    if (readLen < wordLength) continue;

	    str_uint_map_t& readWordOffset(readWordOffsets[readIndex]);

	    for (unsigned j(0); j<=(readLen-wordLength); ++j)
	    {
	    	const std::string word(seq.substr(j,wordLength));
	        if (readWordOffset.find(word) != readWordOffset.end())
	        {
	        	// try again with different k-mer size
	#ifdef DEBUG_ASBL
	            log_os << " word " << word << " repeated in read " << readIndex << "\n";
                log_os << " res = " << readWordOffset.find(word)->second << "\n";
	#endif
	            return false;
	        }

			// record (0-indexed) start point for word in read
	        #ifdef DEBUG_ASBL
			log_os << "Recording " << word << " at " << j << " in " << seq <<  "\n";
	        #endif
			readWordOffset[word]=j;

			// count occurrences
			++wordHash[word];
			/*if (wordCount[word]>maxWordCount)
			{
				maxWordCount  = wordCount[word];
				maxWord = word;
			}*/

			if(reads[readIndex].first<firstWordPos && j==0) {
				firstWord = word;
				firstWordPos = reads[readIndex].first;
			}
		}
	}

/*    if (maxWordCount < opt.minCoverage)
	{
#ifdef DEBUG_ASBL
		log_os << logtag << "Coverage too low : " << maxWordCount << " " << opt.minCoverage << "\n";
#endif
		return false;
	}
*/
	return true;
}

static
bool
buildContigs(
    const SmallAssemblerOptions& /*opt*/,
    const AssemblyReadInput& reads,
    AssemblyReadOutput& readInfo,
    const unsigned wordLength,
    Assembly& finalContigs,
    unsigned& unusedReads
#ifdef DEBUG_ASBL
    , const unsigned iteration
#endif
    )
{
#ifdef DEBUG_ASBL
    static const std::string logtag("buildContigs: ");
    log_os << logtag << "In SmallAssembler::buildContig. word length=" << wordLength << " number of reads " << reads.size() << "\n";
    for (unsigned readIndex(0); readIndex<reads.size(); ++readIndex)
    {
        log_os << reads[readIndex].first << "  " << reads[readIndex].second << " used=" << readInfo[readIndex].isUsed << "\n";
    }
#endif

    const unsigned readCount(reads.size());

    // a set of read hashes; each read hash stores the starting positions of all kmers in the read
    std::vector<str_uint_map_t> readWordOffsets(readCount);

    // counts the number of occurrences for each kmer in all reads
    str_uint_map_t wordHash;
    std::string firstWord("NA");

    buildHash(wordHash,firstWord,readWordOffsets,readInfo,reads,wordLength);

    // most frequent kmer and its number of occurrences
    //unsigned maxWordCount(0);
    //std::string maxWord;

    // check if this is needed
    //log_os << logtag << " pruneHash size before pruning = " << wordCount.size() << "\n";
    pruneHash(wordHash);
    //log_os << logtag << " pruneHash size after pruning = " << wordCount.size() << "\n";


    int firstWordPos(std::numeric_limits<int>::max());
    for (unsigned readIndex(0);readIndex<readCount; ++readIndex)
    {
        const std::string& seq(reads[readIndex].second);
        const unsigned readLen(seq.size());

        for (unsigned j(0); j<=(readLen-wordLength); ++j)
        {
            const std::string word(seq.substr(j,wordLength));
            bool wordInHash(wordHash.find(word) != wordHash.end());
            if(reads[readIndex].first<firstWordPos && j==0 && wordInHash ) {
                firstWord = word;
                firstWordPos = reads[readIndex].first;
            }
        }
    }

    if (firstWord == "NA")
    {
#ifdef DEBUG_ASBL
        log_os << logtag << " No words left after pruning :" << wordHash.size() << "\n";
#endif
        return false;
    }

#ifdef DEBUG_ASBL
    //log_os << logtag << "Seeding kmer : " << maxWord << "\n";
    log_os << logtag << "Seeding kmer : " << firstWord << " cov=" << wordHash[firstWord] << "\n";
#endif

    // start initial assembly with most frequent kmer as seed
    //AssembledContig contig;
    //walk(opt,maxWord,wordLength,wordCount,contig.seq);
    //walk(opt,firstWord,wordLength,wordCount,contig.seq);


    Assembly contigs;
    initDFS(wordHash,firstWord,wordLength,contigs
#ifdef DEBUG_ASBL
    , iteration
#endif
    );

    /*std::cerr << "RESULTS FOR " << firstWord << "\n";
    for (AssembledSequence::const_iterator ct = contigSeq.begin(); ct!=contigSeq.end(); ++ct) {
        std::cerr << "ASSEMBLED CTG " << *ct << "\n";
    }*/

    assert(contigs.size() != 0);
    //contig.seedReadCount = maxWordCount;

    // done with this now:
    wordHash.clear();

    for (Assembly::iterator ctgIter = contigs.begin(); ctgIter != contigs.end(); ++ctgIter)
    {
        //AssembledContig contig;
        //contig.seq = ctgIter->sequence;

        // finally -- set isUsed and decrement unusedReads
        const unsigned contigSize(ctgIter->seq.size());
        for (unsigned readIndex(0); readIndex<readCount; ++readIndex) 
        {
            const str_uint_map_t& readWordOffset(readWordOffsets[readIndex]);
            AssemblyReadInfo& rinfo(readInfo[readIndex]);

            if (rinfo.isUsed) continue;

            // store all reads sharing k-mers of the current word length with the contig
            for (unsigned j(0); j<=(contigSize-wordLength); ++j) 
            {
                const std::string word(ctgIter->seq.substr(j,wordLength));
                //cerr << "Testing word " << word << " " << readNum << "\n";
                //cerr << "with counts : " << wordCount[word] << "\n";
                if (readWordOffset.count(word))
                {
                    rinfo.isUsed = true;
                    rinfo.contigId = finalContigs.size();

                    assert(unusedReads != 0);
                    --unusedReads;

                    ++ctgIter->seedReadCount;
                    break;
                }
            }
        }

        // throw away short stuff
        if (ctgIter->seq.length() < 100) continue;
        finalContigs.push_back(*ctgIter);
    }

    // don't need this anymore:
    readWordOffsets.clear();

    //contigs.push_back(contig);
    if (finalContigs.size() == 0) return false;

    return true;
}



void
runSmallAssembler(
    const SmallAssemblerOptions& opt,
    const AssemblyReadInput& reads,
    AssemblyReadOutput& assembledReadInfo,
    Assembly& contigs)
{
#ifdef DEBUG_ASBL
    static const std::string logtag("runSmallAssembler: ");
    log_os << logtag << "Starting assembly with " << reads.size() << " reads.\n";
#endif
    assert(alphabet.size()>1);

    assembledReadInfo.clear();
    contigs.clear();

    assembledReadInfo.resize(reads.size());

    unsigned wordLength(opt.minWordLength);
    unsigned unusedReads(reads.size());

    for (unsigned iteration(0); iteration < opt.maxAssemblyIterations; ++iteration)
    {
        if (unusedReads < opt.minSeedReads) return;

        const unsigned lastUnusedReads(unusedReads);
        for (; wordLength<=opt.maxWordLength; wordLength+=opt.wordStepSize)
        {
            //std::cerr << "Starting assembly with k=" << wordLength << " iter= " << iteration << "\n";
            const bool isAssemblySuccess = buildContigs(opt, reads, assembledReadInfo, wordLength, contigs, unusedReads
#ifdef DEBUG_ASBL
            , iteration
#endif
            );
            if (isAssemblySuccess) break;
        }

#ifdef DEBUG_ASBL
        log_os << logtag << "iter: " << iteration << " unused readMap now: " << unusedReads << "\n";
#endif

        // stop if no change in number of unused reads
        if (unusedReads == lastUnusedReads)
        {
#ifdef DEBUG_ASBL
            log_os << logtag << "Number of unused reads (" << unusedReads << ") did not change in this iteration. Stopping.\n";
#endif
            return;
        }
    }
#ifdef DEBUG_ASBL
    log_os << logtag << "Reached max number of assembly iterations: " << opt.maxAssemblyIterations << "\n";
#endif
}

