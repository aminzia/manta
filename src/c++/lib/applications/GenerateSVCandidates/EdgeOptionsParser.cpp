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

#include "EdgeOptionsParser.hh"

#include "blt_util/parse_util.hh"
#include "blt_util/string_util.hh"


namespace {
const char locusIndexKey[] = "locus-index";
}



boost::program_options::options_description
getOptionsDescription(
    EdgeOptions& opt)
{
    namespace po = boost::program_options;

    po::options_description optdesc("edge-selection");
    optdesc.add_options()
    ("bin-count", po::value(&opt.binCount)->default_value(opt.binCount),
     "Specify how many bins the SV candidate problem should be divided into, where bin-index can be used to specify which bin to solve")
    ("bin-index", po::value(&opt.binIndex)->default_value(opt.binIndex),
     "specify which bin to solve when the SV candidate problem is subdivided into bins. Value must bin in [0,bin-count)")
    (locusIndexKey, po::value<std::string>(),
     "Instead of solving for all SV candidates in a bin, solve for candidates of a particular locus or edge."
     " If this argument is specified then bin-index is ignored."
     " Argument can be one of { locusIndex , locusIndex:nodeIndex , locusIndex:nodeIndex:nodeIndex },"
     " which will run an entire locus, all edges connected to one node in a locus or a single edge, respectively.")
    ;
    return optdesc;
}



bool
parseOptions(
    const boost::program_options::variables_map& vm,
    EdgeOptions& opt,
    std::string& errorMsg)
{
    errorMsg.clear();

    if (vm.count(locusIndexKey))
    {
        using namespace illumina::blt_util;

        const std::string& locusString(boost::any_cast<std::string>(vm[locusIndexKey].value()));

        std::vector<std::string> indices;
        split_string(locusString, ':', indices);
        if (indices.size() > 3)
        {
            errorMsg="locus-index argument can have no more than 3 colon separated segments";
        }

        assert(! indices.empty());

        opt.locusIndex = parse_unsigned_str(indices[0]);
        opt.isLocusIndex = true;
        if (indices.size() > 1)
        {
            opt.nodeIndex1 = parse_unsigned_str(indices[1]);
            opt.isNodeIndex1 = true;
            if (indices.size() > 2)
            {
                opt.nodeIndex2 = parse_unsigned_str(indices[2]);
                opt.isNodeIndex2 = true;
            }
        }
    }

    if (errorMsg.empty())
    {
        if (opt.binCount < 1)
        {
            errorMsg="bin-count must be 1 or greater";
        }
        else if (opt.binIndex >= opt.binCount)
        {
            errorMsg="bin-index must be in range [0,bin-count)";
        }
    }

    return (! errorMsg.empty());
}
