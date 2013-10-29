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
/// \author Chris Saunders
///

#include "SVFinder.hh"

#include "blt_util/bam_streamer.hh"
#include "common/Exceptions.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "blt_util/bam_header_util.hh"

#include "boost/foreach.hpp"

#include <iostream>



static const bool isExcludeUnpaired(true);

//#define DEBUG_SVDATA

#ifdef DEBUG_SVDATA
#include "blt_util/log.hh"
#endif


SVFinder::
SVFinder(const GSCOptions& opt) :
    _scanOpt(opt.scanOpt),
    _readScanner(_scanOpt,opt.statsFilename,opt.alignFileOpt.alignmentFilename),
    _referenceFilename(opt.referenceFilename)
{
    // load in set:
    _set.load(opt.graphFilename.c_str(),true);

    // setup regionless bam_streams:
    // setup all data for main analysis loop:
    BOOST_FOREACH(const std::string& afile, opt.alignFileOpt.alignmentFilename)
    {
        // avoid creating shared_ptr temporaries:
        streamPtr tmp(new bam_streamer(afile.c_str()));
        _bamStreams.push_back(tmp);
    }
}



// test if read supports an SV on this edge, if so, add to SVData
static
void
addSVNodeRead(
    const std::map<std::string, int32_t>& chromToIndex,
    const SVLocusScanner& scanner,
    const SVLocusNode& localNode,
    const SVLocusNode& remoteNode,
    const bam_record& bamRead,
    const unsigned bamIndex,
    const bool isExpectRepeat,
    const std::string& bkptRef,
    SVCandidateSetReadPairSampleGroup& svDataGroup)
{
    if (scanner.isReadFiltered(bamRead)) return;

    // don't rely on the properPair bit to be set correctly:
    const bool isAnomalous(! scanner.isProperPair(bamRead, bamIndex));

    if (! isAnomalous) return;

    const bool isLargeFragment(scanner.isLargeFragment(bamRead, bamIndex));

    if (! isLargeFragment) return;

#if 0
    /// TODO:  move local-assembly and spanning candidate handling together here:
    const bool isLocalAssemblyEvidence(scanner.isLocalAssemblyEvidence(bamRead));

    if (! (isAnomalous || isLocalAssemblyEvidence)) return;
#endif

    // finally, check to see if the svDataGroup is full... for now, we allow a very large
    // number of reads to be stored in the hope that we never reach this limit, but just in
    // case we don't want to exhaust memory in centromere pileups, etc...
    static const unsigned maxDataSize(2000);
    if (svDataGroup.size() >= maxDataSize)
    {
        if (! svDataGroup.isIncomplete())
        {
            svDataGroup.setIncomplete();
        }
        return;
    }

    typedef std::vector<SVLocus> loci_t;
    loci_t loci;
    scanner.getSVLoci(bamRead,bamIndex,bkptRef,chromToIndex,loci);

    BOOST_FOREACH(const SVLocus& locus, loci)
    {
        if (locus.empty()) continue;
        if (2 != locus.size()) continue;

        unsigned readLocalIndex(0);
        unsigned readRemoteIndex(1);
        if (! locus.getNode(readLocalIndex).isOutCount())
        {
            std::swap(readLocalIndex,readRemoteIndex);
        }

        if (! locus.getNode(readLocalIndex).getInterval().isIntersect(localNode.getInterval())) continue;
        if (! locus.getNode(readRemoteIndex).getInterval().isIntersect(remoteNode.getInterval())) continue;

        svDataGroup.add(bamRead,isExpectRepeat);

        // once any loci has achieved the local/remote overlap criteria, there's no reason to keep scanning loci
        // of the same bam record:
        break;
    }
}



void
SVFinder::
addSVNodeData(
    const std::map<std::string, int32_t>& chromToIndex,
    const SVLocus& locus,
    const NodeIndexType localNodeIndex,
    const NodeIndexType remoteNodeIndex,
    const std::string& referenceFilename,
    SVCandidateSetData& svData)
{
    // get full search interval:
    const SVLocusNode& localNode(locus.getNode(localNodeIndex));
    const SVLocusNode& remoteNode(locus.getNode(remoteNodeIndex));
    GenomeInterval searchInterval(localNode.getInterval());

    searchInterval.range.merge_range(localNode.getEvidenceRange());

    bool isExpectRepeat(svData.setNewSearchInterval(searchInterval));

    /// This is a temporary measure to make the qname collision detection much looser
    /// problems have come up where very large deletions are present in a read, and it is therefore
    /// detected as a repeat in two different regions, even though they are separated by a considerable
    /// distance. Solution is temporarily turn off collision detection whenever two regions are on
    /// the same chrom (ie. almost always)
    ///
    /// TODO: restore more precise collision detection
    if (! isExpectRepeat) isExpectRepeat = (localNode.getInterval().tid == remoteNode.getInterval().tid);

#ifdef DEBUG_SVDATA
    log_os << "addSVNodeData: bp_interval: " << localNode.getInterval()
           << " evidenceInterval: " << localNode.getEvidenceRange()
           << " searchInterval: " << searchInterval
           << " isExpectRepeat: " << isExpectRepeat
           << "\n";
#endif

    // iterate through reads, test reads for association and add to svData:
    unsigned bamIndex(0);
    BOOST_FOREACH(streamPtr& bamPtr, _bamStreams)
    {
        SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
        bam_streamer& read_stream(*bamPtr);
        const bam_header_t& header(*(bamPtr->get_header()));
        const bam_header_info bamHeader(header);

        // set bam stream to new search interval:
        read_stream.set_new_region(searchInterval.tid,searchInterval.range.begin_pos(),searchInterval.range.end_pos());

        reference_contig_segment refSegment;
        getIntervalReferenceSegment(referenceFilename,bamHeader,searchInterval,refSegment);

#ifdef DEBUG_SVDATA
        log_os << "addSVNodeData: scanning bamIndex: " << bamIndex << "\n";
#endif
        while (read_stream.next())
        {
            const bam_record& bamRead(*(read_stream.get_record_ptr()));
            std::string ref;
            {
                ALIGNPATH::path_t apath;
                bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), apath);
                // apparently we need the -1 here
                const int alPos(bamRead.pos()-searchInterval.range.begin_pos()-1);
                const int alLen(apath_ref_length(apath));
                const int refSize(refSegment.seq().size());
                if (alPos < 0 || (alPos+alLen) > refSize) continue;
                ref = refSegment.seq().substr(alPos,alLen);
            }

            // test if read supports an SV on this edge, if so, add to SVData
            addSVNodeRead(chromToIndex,_readScanner,localNode,remoteNode,bamRead, bamIndex,isExpectRepeat,ref,svDataGroup);
        }
        ++bamIndex;
    }
}



// sanity check the final result
void
SVFinder::
checkResult(
    const SVCandidateSetData& svData,
    const std::vector<SVCandidate>& svs) const
{
    using namespace illumina::common;

    const unsigned svCount(svs.size());
    if (0 == svCount) return;

    // check that the counts totaled up from the data match those in the sv candidates
    std::map<unsigned,unsigned> readCounts;
    std::map<unsigned,unsigned> pairCounts;

    for (unsigned i(0); i<svCount; ++i)
    {
        readCounts[i] = 0;
        pairCounts[i] = 0;
    }

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        const SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
        BOOST_FOREACH(const SVCandidateSetReadPair& pair, svDataGroup)
        {
            BOOST_FOREACH(const SVPairAssociation& sva, pair.svLink)
            {

                if (sva.index>=svCount)
                {
                    std::ostringstream oss;
                    oss << "Searching for SVIndex: " << sva.index << " with svSize: " << svCount << "\n";
                    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
                }

                if (SVEvidenceType::isPairType(sva.evtype))
                {
                    if (pair.read1.isSet()) readCounts[sva.index]++;
                    if (pair.read2.isSet()) readCounts[sva.index]++;
                    if (pair.read1.isSet() && pair.read2.isSet()) pairCounts[sva.index] += 2;
                }
            }
        }
    }

    for (unsigned svIndex(0); svIndex<svCount; ++svIndex)
    {
        const unsigned svObsReadCount(svs[svIndex].bp1.getLocalPairCount() + svs[svIndex].bp2.getLocalPairCount());
        const unsigned svObsPairCount(svs[svIndex].bp1.getPairCount() + svs[svIndex].bp2.getPairCount());
        assert(svs[svIndex].bp1.getPairCount() == svs[svIndex].bp2.getPairCount());

        const unsigned dataObsReadCount(readCounts[svIndex]);
        const unsigned dataObsPairCount(pairCounts[svIndex]);

        bool isCountException(false);
        if (isExcludeUnpaired)
        {
            if (svObsReadCount > dataObsReadCount) isCountException=true;
        }
        else
        {
            if (svObsReadCount != dataObsReadCount) isCountException=true;
        }
        if (svObsPairCount != dataObsPairCount) isCountException=true;

        if (isCountException)
        {
            std::ostringstream oss;
            oss << "Unexpected difference in sv and data read counts.\n"
                << "\tSVreadCount: " << svObsReadCount << " DataReadCount: " << dataObsReadCount << "\n"
                << "\tSVpaircount: " << svObsPairCount << " DataPaircount: " << dataObsPairCount << "\n"
                << "\tsvIndex: " << svIndex << " SV: " << svs[svIndex];
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));

        }
    }
}



typedef std::map<unsigned,unsigned> movemap_t;



/// local convenience struct, if only I had closures instead... :<
struct svCandDeleter
{
    svCandDeleter(
        std::vector<SVCandidate>& svs,
        movemap_t& moveSVIndex) :
        _shift(0),
        _isLastIndex(false),
        _lastIndex(0),
        _svs(svs),
        _moveSVIndex(moveSVIndex)
    {}

    void
    deleteIndex(
        const unsigned index)
    {
        assert(index <= _svs.size());

        if (_isLastIndex)
        {
            for (unsigned i(_lastIndex+1); i<index; ++i)
            {
                assert(_shift>0);
                assert(i>=_shift);

                _svs[(i-_shift)] = _svs[i];
                // moveSVIndex has already been set for deleted indices, this sets
                // the move for non-deleted positions:
                _moveSVIndex[i] = (i-_shift);
            }
        }
        _lastIndex=index;
        _isLastIndex=true;
        _shift++;
    }

private:
    unsigned _shift;
    bool _isLastIndex;
    unsigned _lastIndex;
    std::vector<SVCandidate>& _svs;
    movemap_t& _moveSVIndex;
};



// check whether any svs have grown to intersect each other
//
// this is also part of the temp hygen hack, so just make this minimally work:
//
static
void
consolidateOverlap(
    const unsigned bamCount,
    SVCandidateSetData& svData,
    std::vector<SVCandidate>& svs)
{
#ifdef DEBUG_SVDATA
    static const std::string logtag("consolidateOverlap: ");
#endif

    movemap_t moveSVIndex;
    std::set<unsigned> deletedSVIndex;

    std::vector<unsigned> innerIndexShift;

    const unsigned svCount(svs.size());
    for (unsigned outerIndex(1); outerIndex<svCount; ++outerIndex)
    {
        const unsigned prevInnerIndexShift( (outerIndex<=1) ? 0 : innerIndexShift[outerIndex-2]);
        innerIndexShift.push_back(prevInnerIndexShift + deletedSVIndex.count(outerIndex-1));
        for (unsigned innerIndex(0); innerIndex<outerIndex; ++innerIndex)
        {
            if (deletedSVIndex.count(innerIndex)) continue;

            if (svs[innerIndex].isIntersect(svs[outerIndex]))
            {
#ifdef DEBUG_SVDATA
                log_os << logtag << "Merging outer:inner: " << outerIndex << " " << innerIndex << "\n";
#endif
                svs[innerIndex].merge(svs[outerIndex]);
                assert(innerIndexShift.size() > innerIndex);
                assert(innerIndexShift[innerIndex] <= innerIndex);
                moveSVIndex[outerIndex] = (innerIndex - innerIndexShift[innerIndex]);
                deletedSVIndex.insert(outerIndex);
                break;
            }
        }
    }

    if (! deletedSVIndex.empty())
    {
#ifdef DEBUG_SVDATA
        BOOST_FOREACH(const unsigned index, deletedSVIndex)
        {
            log_os << logtag << "deleted index: " << index << "\n";
        }
#endif

        {
            svCandDeleter svDeleter(svs,moveSVIndex);

            BOOST_FOREACH(const unsigned index, deletedSVIndex)
            {
                svDeleter.deleteIndex(index);
            }
            svDeleter.deleteIndex(svCount);
        }

        svs.resize(svs.size()-deletedSVIndex.size());

        // fix indices:
        for (unsigned i(0); i<svs.size(); ++i)
        {
            svs[i].candidateIndex = i;
        }
    }

    if (! moveSVIndex.empty())
    {
#ifdef DEBUG_SVDATA
        BOOST_FOREACH(const movemap_t::value_type& val, moveSVIndex)
        {
            log_os << logtag << "Movemap from: " << val.first << " to: " << val.second << "\n";
        }
#endif

        for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
        {
            SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));
            BOOST_FOREACH(SVCandidateSetReadPair& pair, svDataGroup)
            {
                BOOST_FOREACH(SVPairAssociation& sva, pair.svLink)
                {
                    if (moveSVIndex.count(sva.index))
                    {
                        sva.index = moveSVIndex[sva.index];
                    }
                }
            }
        }
    }
}



void
SVFinder::
getCandidatesFromData(
    const std::map<std::string, int32_t>& chromToIndex,
    SVCandidateSetData& svData,
    std::vector<SVCandidate>& svs)
{
#ifdef DEBUG_SVDATA
    static const std::string logtag("getCandidatesFromData: ");
#endif

    std::vector<SVCandidate> readCandidates;

    const unsigned bamCount(_bamStreams.size());

    //unsigned bamIndex(0);
    for (unsigned bamIndex(0);bamIndex<bamCount;++bamIndex)
    //BOOST_FOREACH(streamPtr& bamPtr, _bamStreams)
    {
        SVCandidateSetReadPairSampleGroup& svDataGroup(svData.getDataGroup(bamIndex));

        // OST: how do I get the reference here??
        /*bam_streamer& read_stream(*bamPtr);
        const bam_header_t& header(*(bamPtr->get_header()));
        const bam_header_info bamHeader(header);

        reference_contig_segment refSegment;
        getReference(referenceFilename,bamHeader,refSegment);*/

        BOOST_FOREACH(SVCandidateSetReadPair& pair, svDataGroup)
        {
            SVCandidateSetRead* localReadPtr(&(pair.read1));
            SVCandidateSetRead* remoteReadPtr(&(pair.read2));
            pair.svLink.clear();

            if (isExcludeUnpaired)
            {
                // in this case both sides of the read pair need to be observed (and not filtered for MAPQ, etc)
                if ((! localReadPtr->isSet()) || (! remoteReadPtr->isSet())) continue;
            }
            else
            {
                if (! localReadPtr->isSet())
                {
                    std::swap(localReadPtr,remoteReadPtr);
                }
            }
            const bam_record* remoteBamRecPtr( remoteReadPtr->isSet() ? &(remoteReadPtr->bamrec) : NULL);
            _readScanner.getBreakendPair(localReadPtr->bamrec, remoteBamRecPtr, bamIndex, chromToIndex, "NA",readCandidates);

#ifdef DEBUG_SVDATA
            log_os << "Checking pair: " << pair << "\n";
            log_os << "Translated to candidates:\n";
            BOOST_FOREACH(const SVCandidate& cand, readCandidates)
            {
                log_os << logtag << "cand: " << cand << "\n";
            }
#endif

            // temporary hack hypoth gen method assumes that only one SV exists for each overlapping breakpoint range with
            // the same orientation:
            //
            // we anticipate so few svs from the POC method, that there's no indexing on them
            // OST 26/09/2013: Be careful when re-arranging or rewriting the code below, under g++ 4.1.2
            // this can lead to an infinite loop.
            BOOST_FOREACH(const SVCandidate& readCand, readCandidates)
            {
#ifdef DEBUG_SVDATA
                log_os << logtag << "Starting assignment for read cand: " << readCand << "\n";
#endif
                bool isSVFound(false);
                unsigned svIndex(0);

                BOOST_FOREACH(SVCandidate& sv, svs)
                {
                    if (sv.isIntersect(readCand))
                    {
#ifdef DEBUG_SVDATA
                        log_os << logtag << "Adding to svIndex: " << svIndex << " match_sv: " << sv << "\n";
#endif
                        isSVFound=true;
                        {
                            using namespace SVEvidenceType;

                            index_t evType(UNKNOWN);
                            for (int i(0); i<SIZE; ++i)
                            {
                                if (readCand.bp1.lowresEvidence.getVal(i) != 0)
                                {
                                    evType = static_cast<index_t>(i);
                                    break;
                                }
                            }

                            pair.svLink.push_back(SVPairAssociation(svIndex,evType));
                        }
                        sv.merge(readCand);
                        break;
                    }
                    svIndex++;
                }

                if (isSVFound) continue;

#ifdef DEBUG_SVDATA
                log_os << logtag << "New svIndex: " << svs.size() << "\n";
#endif
                pair.svLink.push_back(SVPairAssociation(svs.size()));
                svs.push_back(readCand);
                svs.back().candidateIndex = pair.svLink.back().index;
            }
        }
        ++bamIndex;
    }

#ifdef DEBUG_SVDATA
    {
        log_os << logtag << "precount: " << svs.size() << "\n";

        unsigned svIndex(0);
        BOOST_FOREACH(SVCandidate& sv, svs)
        {
            log_os << logtag << "PRECOUNT: index: " << svIndex << " " << sv;
            svIndex++;
        }
    }
#endif

    consolidateOverlap(bamCount,svData,svs);

#ifdef DEBUG_SVDATA
    {
        log_os << logtag << "postcount: " << svs.size() << "\n";

        unsigned svIndex(0);
        BOOST_FOREACH(SVCandidate& sv, svs)
        {
            log_os << logtag << "POSTCOUNT: index: " << svIndex << " " << sv;
            svIndex++;
        }
    }
#endif
}



void
SVFinder::
findCandidateSV(
    const std::map<std::string, int32_t>& chromToIndex,
    const EdgeInfo& edge,
    const std::string& referenceFilename,
    SVCandidateSetData& svData,
    std::vector<SVCandidate>& svs
    )
{
    svData.clear();
    svs.clear();

    const SVLocusSet& set(getSet());
    const unsigned minEdgeCount(set.getMinMergeEdgeCount());

#ifdef DEBUG_SVDATA
    log_os << "SVDATA: Evaluating edge: " << edge << "\n";
#endif

    // first determine if this is an edge we're going to evaluate
    //
    // edge must be bidirectional at the noise threshold of the locus set:
    const SVLocus& locus(set.getLocus(edge.locusIndex));

    if ((locus.getEdge(edge.nodeIndex1,edge.nodeIndex2).getCount() < minEdgeCount) ||
        (locus.getEdge(edge.nodeIndex2,edge.nodeIndex1).getCount() < minEdgeCount))
    {
#ifdef DEBUG_SVDATA
        log_os << "SVDATA: Edge failed min edge count.\n";
#endif
        return;
    }

    // if this is a self-edge, then automatically forward it as is to the assembly module:
    /// TODO: move self-edge handling into the regular hygen routine below
    if (edge.nodeIndex1 == edge.nodeIndex2)
    {
        SVCandidate sv;
        SVBreakend& localBreakend(sv.bp1);
        SVBreakend& remoteBreakend(sv.bp2);

        const SVLocusNode& node(locus.getNode(edge.nodeIndex1));

        static const SVEvidenceType::index_t svUnknown(SVEvidenceType::UNKNOWN);
        localBreakend.lowresEvidence.add(svUnknown, node.getEdge(edge.nodeIndex1).getCount());
        localBreakend.state = SVBreakendState::COMPLEX;
        localBreakend.interval = node.getInterval();

        remoteBreakend.state = SVBreakendState::UNKNOWN;

        sv.candidateIndex=svs.size();
        svs.push_back(sv);

        // minimal setup for svData:
        const unsigned bamCount(_bamStreams.size());
        for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
        {
            svData.getDataGroup(bamIndex);
        }

        svData.setSkipped();

        return;
    }


    // start gathering evidence required for hypothesis generation,
    //
    // first step is to scan through each region, and identify all reads supporting
    // some sort of breakend in the target region, then match up pairs so that they
    // can easily be accessed from each other
    //

    // steps:
    // iterate through regions -- for each region walk from evidence range to breakpoint range picking up all reads associated with the breakend
    // make a common code to determine read breakend association shared with Estimation step
    // pair reads in data structure
    // determine the number of breakends from simple orientation and intersecting joint region logic
    // assign data to each breakend candidates
    // come up with an ultra-simple model-free scoring rule: >10 obs = Q60,k else Q0
    //

    addSVNodeData(chromToIndex, locus,edge.nodeIndex1,edge.nodeIndex2,referenceFilename,svData);
    if (edge.nodeIndex1 != edge.nodeIndex2)
    {
        addSVNodeData(chromToIndex, locus,edge.nodeIndex2,edge.nodeIndex1,referenceFilename,svData);
    }

    getCandidatesFromData(chromToIndex, svData,svs);

    //checkResult(svData,svs);
}

