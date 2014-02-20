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
/// \author Chris Saunders
///

#pragma once

#include "alignment/AlignmentScores.hh"
#include "options/SmallAssemblerOptions.hh"


/// Options for the SV refiner step
///
/// Note that we have two categories of options for assembly and alignment,
/// one for small events, and one for large events
///
struct SVRefinerOptions
{
    /// match, mismatch, open score ratios taken from bwa defaults (but not extend!) :
    ///
    SVRefinerOptions() :
        smallSVAlignScores(2, -8, -12, 0, -1),
        largeInsertAlignScores(2, -8, -12, -1, -1),
        spanningAlignScores(2, -8, -12, -1, -15, -1, 0),
        jumpScore(-25)
    {
        spanningAssembleOpt.minContigLength=75; ///< For breakend-spanning assemblies we require a larger contig than for small-variant assemblies
    }

    /// parameters for small SV assembly/alignment:
    AlignmentScores<int> smallSVAlignScores;
    AlignmentScores<int> largeInsertAlignScores;
    SmallAssemblerOptions smallSVAssembleOpt;

    // parameters for large SV assembly/alignment:
    AlignmentScores<int> spanningAlignScores;
    const int jumpScore;
    SmallAssemblerOptions spanningAssembleOpt;
};
