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

#pragma once

#include "SVScorer.hh"


void
getBreakendNoiseScore(
    const SVCandidate& sv,
    const std::vector<SVCandidate>& svs,
    const std::vector<SVCandidate>& offEdgeSvs,
    const bool isBp1,
    float& score);
