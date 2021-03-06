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

/// \author Chris Saunders
///

#pragma once

#include <iosfwd>
#include <string>


/// pretty print sequence in such a way that it's easy to locate position number
///
void
printSeq(
    const char* seq,
    std::ostream& os);


inline
void
printSeq(
    const std::string& seq,
    std::ostream& os)
{
    printSeq(seq.c_str(),os);
}
