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

#pragma once

#include <vector>
#include <string>


/// information added to each read in the process of assembly
///
struct AssemblyReadInfo
{
    AssemblyReadInfo() :
        isUsed(false),
        contigId(0)
    {}

    bool isUsed;
    unsigned contigId; ///< index of the contig that this read is used in
};


typedef std::vector<std::pair<int,std::string> > AssemblyReadInput;
typedef std::vector<bool> AssemblyReadReversal;
typedef std::vector<AssemblyReadInfo> AssemblyReadOutput;
