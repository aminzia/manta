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

#include "manta/Program.hh"


/// generates candidate calls from graph edges
///
struct GenerateSVCandidates : public manta::Program
{

    const char*
    name() const
    {
        return "GenerateSVCandidates";
    }

    void
    runInternal(int argc, char* argv[]) const;
};
