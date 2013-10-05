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

#pragma once

#include "blt_util/bam_header_info.hh"
#include "blt_util/reference_contig_segment.hh"
#include "manta/SVCandidate.hh"

#include <string>


/// test whether the two jump aligner regions will overlap
bool
isRefRegionOverlap(
    const bam_header_info& header,
    const pos_t extraRefEdgeSize,
    const SVCandidate& sv);

bool
isRefRegionOverlapByPos(
    const SVCandidate& sv,
    const pos_t beginPos1,
    const pos_t endPos1,
    const pos_t beginPos2,
    const pos_t endPos2);

void
getIntervalReferenceSegment(
    const std::string& referenceFilename,
    const bam_header_info& header,
    const GenomeInterval& interval,
    reference_contig_segment& intervalRef);

/// extract the reference sequence around each breakend into a reference_contig_segment
/// object
///
/// for each region, we extract the hypothetical breakend region + extraRefEdgeSize bases
/// on each side
///
void
getSVReferenceSegments(
    const std::string& referenceFilename,
    const bam_header_info& header,
    const SVCandidate& sv,
    const pos_t beginPos1,
    const pos_t endPos1,
    const pos_t beginPos2,
    const pos_t endPos2,
    reference_contig_segment& bp1ref,
    reference_contig_segment& bp2ref);
