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

#include "SummarizeAlignmentStats.hh"
#include "SASOptions.hh"

#include "blt_util/log.hh"
#include "manta/ReadGroupStatsSet.hh"

#include "boost/foreach.hpp"

#include <iostream>



static
void
runSAS(const SASOptions& opt)
{
    static const float quantLevel[] = { 0.25, 0.5, 0.75, 0.9, 0.95, 0.99 };
    static const unsigned quantLevelCount(sizeof(quantLevel)/sizeof(float));

    std::ostream& report_os(std::cout);

    ReadGroupStatsSet rgss;
    rgss.load(opt.statsFilename.c_str());

    const unsigned groupCount(rgss.size());
    for (unsigned i(0); i<groupCount; ++i)
    {
        const ReadGroupStatsSet::KeyType& key(rgss.getKey(i));
#ifdef READ_GROUPS
        report_os << "bamFile:\t" << key.bamLabel << '\n';
        report_os << "readGroup:\t" << key.rgLabel << '\n';
#else
        report_os << "group:\t" << key.bamLabel << '\n';
#endif

        const ReadGroupStats& rgs(rgss.getStats(i));
        report_os << "fragment length observations:\t" << rgs.fragStats.totalObservations() << '\n';
        report_os << "fragment length quantiles:\n";
        for (unsigned j(0); j<quantLevelCount; ++j)
        {
            report_os << quantLevel[j] << '\t' << rgs.fragStats.quantile(quantLevel[j]) << '\n';
        }
        report_os << '\n';
    }
}



void
SummarizeAlignmentStats::
runInternal(int argc, char* argv[]) const
{

    SASOptions opt;

    parseSASOptions(*this,argc,argv,opt);
    runSAS(opt);
}
