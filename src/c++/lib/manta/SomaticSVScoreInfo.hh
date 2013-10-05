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
/// \author Chris Saunders and Xiaoyu Chen
///

#pragma once


#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVCandidate.hh"

#include <iosfwd>
#include <string>
#include <set>

struct SVAlignmentInfo
{
    SVAlignmentInfo(
        const SVCandidate& sv,
        const SVCandidateAssemblyData& assemblyData);

    std::string contigSeq;
    std::string bp1RefSeq;
    std::string bp2RefSeq;
    unsigned bp1ContigOffset;
    unsigned bp2ContigOffset;
    bool bp1ContigReversed;
    bool bp2ContigReversed;
    unsigned bp1RefOffset;
    unsigned bp2RefOffset;
};

struct SVSampleInfo
{

    SVSampleInfo() :
        bp1SpanReads(0),
        bp2SpanReads(0),
        spanPairs(0),
        contigSRCount(0),
        contigSREvidence(0),
        contigSRMapQ(0),
        refSRCount(0),
        refSREvidence(0),
        refSRMapQ(0)
    {}

    void
    clear()
    {
        bp1SpanReads = 0;
        bp2SpanReads = 0;
        spanPairs = 0;
        contigSRCount = 0;
        contigSREvidence = 0;
        contigSRMapQ = 0;
        refSRCount = 0;
        refSREvidence = 0;
        refSRMapQ = 0;
    }

    unsigned bp1SpanReads;
    unsigned bp2SpanReads;
    unsigned spanPairs;
    unsigned contigSRCount;
    float contigSREvidence;
    float contigSRMapQ;
    unsigned refSRCount;
    float refSREvidence;
    float refSRMapQ;
};

std::ostream&
operator<<(std::ostream& os, const SVSampleInfo& si);

/// consolidate all somatic scoring results applied to an SV candidate
struct SomaticSVScoreInfo
{
    SomaticSVScoreInfo() :
        bp1MaxDepth(0),
        bp2MaxDepth(0),
        somaticScore(0)
    {}

    void
    clear()
    {
        tumor.clear();
        normal.clear();
        filters.clear();

        bp1MaxDepth=0;
        bp2MaxDepth=0;
        somaticScore=0;
    }

    SVSampleInfo tumor;
    SVSampleInfo normal;

    std::set<std::string> filters;

    unsigned bp1MaxDepth;
    unsigned bp2MaxDepth;
    unsigned somaticScore;
};

std::ostream&
operator<<(std::ostream& os, const SomaticSVScoreInfo& ssi);


