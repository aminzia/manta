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
/// \author Chris Saunders and Xiaoyu Chen
///

#include "SVScorerBreakendInfo.hh"

#include "blt_util/align_path_bam_util.hh"
#include "manta/SVCandidateUtil.hh"



/// add bam alignment to simple short-range vector depth estimate
///
/// \param[in] beginPos this is the begin position of the range covered by the depth array
///
static
void
addReadToDepthEst(
    const bam_record& bamRead,
    const pos_t beginPos,
    std::vector<unsigned>& depth)
{
    using namespace ALIGNPATH;

    const pos_t endPos(beginPos+depth.size());

    // get cigar:
    path_t apath;
    bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), apath);

    pos_t refPos(bamRead.pos()-1);
    BOOST_FOREACH(const path_segment& ps, apath)
    {
        if (refPos>=endPos) return;

        if (is_segment_align_match(ps.type))
        {
            for (pos_t pos(refPos); pos < (refPos+static_cast<pos_t>(ps.length)); ++pos)
            {
                if (pos>=beginPos)
                {
                    if (pos>=endPos) return;
                    depth[pos-beginPos]++;
                }
            }
        }
        if (is_segment_type_ref_length(ps.type)) refPos += ps.length;
    }
}



void
SVScorer::
getBreakendMaxMappedDepthAndMQ0(
    const bool isMaxDepth,
    const double cutoffDepth,
    const SVBreakend& bp,
    unsigned& maxDepth,
    float& MQ0Frac)
{
    /// define a new interval -/+ 50 bases around the center pos
    /// of the breakpoint
    static const pos_t regionSize(50);

    maxDepth=0;
    MQ0Frac=0;

    unsigned totalReads(0);
    unsigned totalMQ0Reads(0);

    const pos_t centerPos(bp.interval.range.center_pos());
    const known_pos_range2 searchRange(std::max((centerPos-regionSize),0), (centerPos+regionSize));

    if (searchRange.size() == 0) return;

    std::vector<unsigned> depth(searchRange.size(),0);

    bool isCutoff(false);
    bool isNormalFound(false);

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        if (_isAlignmentTumor[bamIndex]) continue;
        isNormalFound=true;

        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        // set bam stream to new search interval:
        bamStream.set_new_region(bp.interval.tid, searchRange.begin_pos(), searchRange.end_pos());

        while (bamStream.next())
        {
            const bam_record& bamRead(*(bamStream.get_record_ptr()));

            // turn filtration down to mapped only to match depth estimate method:
            if (bamRead.is_unmapped()) continue;

            const pos_t refPos(bamRead.pos()-1);
            if (refPos >= searchRange.end_pos()) break;

            addReadToDepthEst(bamRead,searchRange.begin_pos(),depth);

            totalReads++;
            if (0 == bamRead.map_qual()) totalMQ0Reads++;

            if (isMaxDepth)
            {
                const pos_t depthOffset(refPos-searchRange.begin_pos());
                if (depthOffset>=0)
                {
                    if (depth[depthOffset] > cutoffDepth)
                    {
                        isCutoff=true;
                        break;
                    }
                }
            }
        }

        if (isCutoff) break;
    }

    assert(isNormalFound);

    maxDepth = *(std::max_element(depth.begin(),depth.end()));
    if (totalReads>=10)
    {
        MQ0Frac = static_cast<float>(totalMQ0Reads)/static_cast<float>(totalReads);
    }
}



static
void
addSVsToNoiseScore(
    const GenomeInterval& region,
    const std::vector<SVCandidate>& svs,
    float& score)
{
    BOOST_FOREACH(const SVCandidate& sv, svs)
    {
        if (isSVBelowMinSize(sv,1000)) continue;
        if (sv.bp1.interval.isIntersect(region) ||
            sv.bp2.interval.isIntersect(region))
        {
            score += 1;
        }
    }
}



void
getBreakendNoiseScore(
    const SVCandidate& sv,
    const std::vector<SVCandidate>& svs,
    const std::vector<SVCandidate>& offEdgeSvs,
    const bool isBp1,
    float& score)
{
    static const pos_t noiseSpan(250);

    const GenomeInterval& bpRegion(sv.getBp(isBp1).interval);
    GenomeInterval bpTestRegion;
    bpTestRegion.tid = bpRegion.tid;
    bpTestRegion.range.set_begin_pos(bpRegion.range.center_pos()-noiseSpan);
    bpTestRegion.range.set_end_pos(bpRegion.range.center_pos()+noiseSpan);

    score=0;
    addSVsToNoiseScore(bpTestRegion,svs,score);
    addSVsToNoiseScore(bpTestRegion,offEdgeSvs,score);
}

