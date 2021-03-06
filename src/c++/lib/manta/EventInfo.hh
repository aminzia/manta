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

#include <string>


/// data related to a multi-junction event
///
struct EventInfo
{
    EventInfo() :
        junctionCount(1)
    {}

    bool
    isEvent() const
    {
        return (! label.empty());
    }

    unsigned junctionCount;
    std::string label;
};
