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


/// Input parameters for SmallAssembler, a simple local de-bruijn graph assembler
///
struct SmallAssemblerOptions
{
    /// sets reasonable default values for 30x DNA-seq, 100-150bp reads
    SmallAssemblerOptions() :
        minWordLength(31),
        maxWordLength(61),
        wordStepSize(5),
        minContigLength(15),
        minCoverage(1),
        maxError(0.35),
        minSeedReads(3),
        //maxAssemblyIterations(3)
        maxAssemblyIterations(10)
    {}

    unsigned minWordLength; ///< initial word (kmer) length
    unsigned maxWordLength; ///< max word length
    unsigned wordStepSize;
    unsigned minContigLength; ///< min contig size
    unsigned minCoverage; ///< min. coverage required for contig extension
    double maxError; ///< max error rates allowed during contig extension
    unsigned minSeedReads; ///< min. number of reads required to start assembly
    unsigned maxAssemblyIterations; ///< Max. number of assembly iterations for a cluster before we give up
};

