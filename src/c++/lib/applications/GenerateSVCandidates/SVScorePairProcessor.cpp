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

#include "SVScorePairProcessor.hh"



SVScorePairInitParams::
SVScorePairInitParams(
    const SVLocusScanner& readScanner,
    const SVCandidate& sv,
    const bool isBp1)
{
    /// In case of breakend microhomology approximate the breakend as a point event at the center of the possible range:
    const pos_t centerPos1 = (sv.bp1.interval.range.center_pos());
    const pos_t centerPos2 = (sv.bp2.interval.range.center_pos());

    centerPos = ( isBp1 ? centerPos1 : centerPos2 );

    beginCenterPos = centerPos1;
    endCenterPos = centerPos2;
    if (endCenterPos < beginCenterPos)
    {
        std::swap(beginCenterPos, endCenterPos);
    }

    // total impact of the alt allele on template size, assuming a simple indel:
    altShift = (std::abs(centerPos2-centerPos1)-sv.insertSeq.size());

    minMapQ = (readScanner.getMinMapQ());
}



const GenomeInterval&
SVScorePairProcessor::
nextBamIndex(
    const unsigned bamIndex)
{
    bparams.isSet = true;
    bparams.isTumor = (isAlignmentTumor[bamIndex]);

    /// set the search range around centerPos so that we can get any fragments at the Xth percentile length or smaller which could have
    /// min Fragsupport
    const SVLocusScanner::Range& pRange(readScanner.getEvidencePairRange(bamIndex));
    bparams.minFrag = (static_cast<pos_t>(pRange.min));
    bparams.maxFrag = (static_cast<pos_t>(pRange.max));

    const pos_t maxSupportedFrag(bparams.maxFrag-pairOpt.minFragSupport);

    const pos_t beginPos(iparams.centerPos-maxSupportedFrag);
    const pos_t endPos(iparams.centerPos+maxSupportedFrag+1);
#ifdef DEBUG_MEGAPAIR
    log_os << __FUNCTION__ << ": pair scan begin/end: " << beginPos << " " << endPos << "\n";
#endif

    bparams.fragDistroPtr = &(readScanner.getFragSizeDistro(bamIndex));

    // set bam stream to new search interval:
    const SVBreakend bp( isBp1 ? sv.bp1 : sv.bp2 );
    bparams.interval = GenomeInterval(bp.interval.tid, beginPos, endPos);

    return bparams.interval;
}
