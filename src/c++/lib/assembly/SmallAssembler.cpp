// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
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



// maps kmers to positions in read
typedef boost::unordered_map<std::string,unsigned> str_uint_map_t;
// keeps track of k-mers
typedef boost::unordered_map<std::string,bool> str_bool_map_t;

static const std::string alphabet("ACGT");


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
    if (isEnd) return contig + base;
    else       return base + contig;
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
    assert(length <= csize);

    if (isEnd) return contig.substr((csize-length),length);
    else       return contig.substr(0,length);
}

static
void
dumpHash(const str_uint_map_t& wordCount,
         const unsigned& wordLength) {


	std::ofstream outFile;
	outFile.open("debruijn.graph.dot");
	outFile << "graph {\n";
    str_uint_map_t aliasH;
    unsigned n(0);
	for (str_uint_map_t::const_iterator ct = wordCount.begin();ct!=wordCount.end();++ct) {
        aliasH[ct->first] = n++;
		outFile << aliasH[ct->first] << "[label=\"cov" << ct->second << "\"]\n"; 
	}
    // need to add edges here
    static const bool isEnd(true);
	for (str_uint_map_t::const_iterator ct = wordCount.begin();ct!=wordCount.end();++ct) {
        std::string tmp(getEnd(ct->first,wordLength-1,isEnd));
        BOOST_FOREACH(const char symbol, alphabet) {
            const std::string newKey(addBase(tmp,symbol,isEnd));
            if (wordCount.find(newKey) != wordCount.end()) {
                outFile << aliasH[ct->first] << " -- " << aliasH[newKey] << ";\n";
            }
        }
	}
	outFile << "}\n";
	outFile.close();
}

/* Prototype of a DFS traversal, not used at the moment
static
void
doDFS(const str_uint_map_t& wordCount,
	  const std::string& v,
	  const unsigned wordLength,
	  str_bool_map_t& seenVertices,
	  const std::string& tmpContig,
	  std::vector<std::string>& contigs) {

	static const bool isEnd(true);

	// when to add a new contig to list?
	// if no new neighbours are left!

	seenVertices[v] = true;

	std::string tmp(getEnd(v,wordLength-1,isEnd));
	bool neighbourFound(false);
	BOOST_FOREACH(const char symbol, alphabet) {
		const std::string newKey(addBase(tmp,symbol,isEnd));
	    if (wordCount.find(newKey) != wordCount.end() && !seenVertices[v]) {
	    	std::string newCtg = tmpContig;
	    	newCtg += symbol;
	    	doDFS(wordCount,newKey,wordLength,seenVertices,newCtg,contigs);
	    	neighbourFound=true;
	    }
	}
	if (!neighbourFound) {
		contigs.push_back(tmpContig);
	}
}

static
void
initDFS(const str_uint_map_t& wordCount,
		const std::string& seed,
		const unsigned wordLength,
		std::vector<std::string>& contigs) {

	str_bool_map_t seenVertices;
	std::string tmpContig("");
	doDFS(wordCount, seed, wordLength, seenVertices, tmpContig, contigs);
}
*/

/**
 * Extends the seed contig (aka most frequent k-mer)
 *
 */
static
void
walk(const SmallAssemblerOptions& opt,
     const std::string& seed,
     const unsigned wordLength,
     const str_uint_map_t& wordCount,
     std::string& contig)
{
    // we start with the seed
    contig = seed;

#ifdef DEBUG_ASBL
    dumpHash(wordCount,wordLength);
#endif

    std::set<std::string> seenBefore;	// records k-mers already encountered during extension
    seenBefore.insert(contig);

    const str_uint_map_t::const_iterator wordCountEnd(wordCount.end());

    // 0 => walk to the right, 1 => walk to the left
    for (unsigned mode(0); mode<2; ++mode)
    {
        const bool isEnd(mode==0);

        while (true)
        {
            const std::string tmp(getEnd(contig, wordLength-1, isEnd));

#ifdef DEBUG_ASBL
            log_os << "# current contig : " << contig << " size : " << contig.size() << "\n"
                   << " getEnd : " << tmp << "\n";
#endif

            if (seenBefore.count(tmp))
            {
#ifdef DEBUG_ASBL
                log_os << "Seen word " << tmp << " before on this walk, terminating" << "\n";
#endif
                break;
            }

            seenBefore.insert(tmp);

            unsigned maxBaseCount(0);
            unsigned totalBaseCount(0);
            char maxBase(alphabet[0]);

            BOOST_FOREACH(const char symbol, alphabet)
            {
                const std::string newKey(addBase(tmp, symbol, isEnd));
                const str_uint_map_t::const_iterator wordCountIter(wordCount.find(newKey));
                if (wordCountIter == wordCountEnd) continue;
#ifdef DEBUG_ASBL
                log_os << "Extending end : base " << symbol << " " << newKey << " " << wordCountIter->second  << "\n";
#endif

                const unsigned val(wordCountIter->second);
                totalBaseCount += val;
                if (val > maxBaseCount)
                {
                    maxBaseCount  = val;
                    maxBase = symbol;
                }
            }
#ifdef DEBUG_ASBL
            log_os << "Winner is : " << maxBase << " with " << maxBaseCount << " occurrences." << "\n";
#endif

            if ((maxBaseCount < opt.minCoverage) ||
                (maxBaseCount < (1.- opt.maxError)* totalBaseCount))
            {
#ifdef DEBUG_ASBL
                log_os << "Coverage or error rate below threshold.\n"
                       << "maxBaseCount : " << maxBaseCount << " minCverage: " << opt.minCoverage << "\n"
                       << "error threshold: " << ((1.0-opt.maxError)* totalBaseCount) << "\n";
#endif
                break;
            }

            /// double check that word exists in reads at least once:
            if (maxBaseCount == 0) break;

#ifdef DEBUG_ASBL
            log_os << "Adding base " << contig << " " << maxBase << " " << mode << "\n";
#endif
            contig = addBase(contig,maxBase,isEnd);
#ifdef DEBUG_ASBL
            log_os << "New contig: " << contig << "\n";
#endif
        }
#ifdef DEBUG_ASBL
        log_os << "mode change. Current mode " << mode << "\n";
#endif
    }
}

static
bool
buildContigs(
    const SmallAssemblerOptions& opt,
    const AssemblyReadInput& reads,
    AssemblyReadOutput& readInfo,
    const unsigned wordLength,
    Assembly& contigs,
    unsigned& unusedReads)
{
    const unsigned readCount(reads.size());

#ifdef DEBUG_ASBL
    static const std::string logtag("buildContigs: ");
    log_os << logtag << "In SmallAssembler::buildContig. word length=" << wordLength << " readCount: " << readCount << "\n";
    for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
    {
        log_os << reads[readIndex] << " used=" << readInfo[readIndex].isUsed << "\n";
    }
#endif

    // a set of read hashes; each read hash stores the starting positions of all kmers in the read
    std::vector<str_uint_map_t> readWordOffsets(readCount);

    // counts the number of occurrences for each kmer in all reads
    str_uint_map_t wordCount;

    // most frequent kmer and its number of occurrences
    unsigned maxWordCount(0);
    std::string maxWord;

    int firstWordPos(std::numeric_limits<int>::max());
    std::string firstWord;

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
                log_os << logtag << "word " << word << " repeated in read " << readIndex << "\n";
#endif
                return false;
            }

            // record (0-indexed) start point for word in read
            //cout << "Recording " << word << " at " << j << "\n";
            readWordOffset[word]=j;

            // count occurrences
            ++wordCount[word];
            if (wordCount[word]>maxWordCount)
            {
                //cout << "Setting max word to " << maxWord << " " << maxOcc << "\n";
                maxWordCount  = wordCount[word];
                maxWord = word;
            }

            if(reads[readIndex].first<firstWordPos && j==0) {
            	firstWord = word;
            	firstWordPos = reads[readIndex].first;
            }
        }
    }

    if (maxWordCount < opt.minCoverage)
    {
#ifdef DEBUG_ASBL
        log_os << logtag << "Coverage too low : " << maxWordCount << " " << opt.minCoverage << "\n";
#endif
        return false;
    }

#ifdef DEBUG_ASBL
    //log_os << logtag << "Seeding kmer : " << maxWord << "\n";
    log_os << logtag << "Seeding kmer : " << firstWord << "\n";
#endif

    // start initial assembly with most frequent kmer as seed
    AssembledContig contig;
    //walk(opt,maxWord,wordLength,wordCount,contig.seq);
    walk(opt,firstWord,wordLength,wordCount,contig.seq);

    // done with this now:
    wordCount.clear();

    const unsigned contigSize(contig.seq.size());

#ifdef DEBUG_ASBL
    log_os << logtag << "First pass assembly resulted in "
           << contig.seq << "\n"
           << " with length " << contigSize << ". Input consisted of " << readCount << " reads.\n";
#endif

    // increment number of reads containing the seeding kmer
    for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
    {
        const str_uint_map_t& readWordOffset(readWordOffsets[readIndex]);
        if (readWordOffset.count(maxWord)) ++contig.seedReadCount;
    }

#ifdef DEBUG_ASBL
    log_os << logtag << "final seeding reading count: " << contig.seedReadCount << "\n";
#endif
    if (contig.seedReadCount < opt.minSeedReads)
    {
#ifdef DEBUG_ASBL
        log_os << "\t...which is below minSeedReadCount of " << opt.minSeedReads << " discarding.\n";
#endif
        return false;
    }

    // finally -- set isUsed and decrement unusedReads
    for (unsigned readIndex(0); readIndex<readCount; ++readIndex)
    {
        const str_uint_map_t& readWordOffset(readWordOffsets[readIndex]);
        AssemblyReadInfo& rinfo(readInfo[readIndex]);

        if (rinfo.isUsed) continue;

        // store all reads sharing k-mers of the current word length with the contig
        // TODO: check if we still needs this
        for (unsigned j(0); j<=(contigSize-wordLength); ++j)
        {
            const std::string word(contig.seq.substr(j,wordLength));
            //cout << "Testing word " << word << " " << readNum << "\n";
            //cout << "with counts : " << wordCount[word] << "\n";
            if (readWordOffset.count(word))
            {
                rinfo.isUsed = true;
                rinfo.contigId = contigs.size();

                assert(unusedReads != 0);
                --unusedReads;
                break;
            }
        }
    }

    // don't need this anymore:
    readWordOffsets.clear();

    contigs.push_back(contig);
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
            const bool isAssemblySuccess = buildContigs(opt, reads, assembledReadInfo, wordLength, contigs, unusedReads);
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

