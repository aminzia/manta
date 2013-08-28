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
// <https://github.com/downloads/sequencing/licenses/>.
//

#include "ReadGroupStatsSet.hh"

#include "blt_util/log.hh"
#include "blt_util/io_util.hh"
#include "blt_util/parse_util.hh"
#include "blt_util/string_util.hh"

#include <fstream>
#include <iostream>



// Stats file data format
const unsigned HEAD_FILL_IDX = 0;
const unsigned HEAD_SRC_IDX  = 1;
const unsigned HEAD_NAME_IDX = 2;

// Stats file data format
const unsigned STAT_SOURCE_IDX             = 0;
const unsigned STAT_INS_SIZE_SD_IDX        = 1;
const unsigned STAT_INS_SIZE_MEDIAN_IDX    = 2;
const unsigned STAT_REL_ORIENT_IDX         = 3;


//TODO: implement deserizlize from file
void
ReadGroupStatsSet::
read(const char* filename)
{
    assert(NULL != filename);

    std::ifstream ifs;
    open_ifstream(ifs,filename);
    this->read(ifs);
}

void
ReadGroupStatsSet::
read(std::istream& is)
{

    using namespace illumina::blt_util;

    clear();
    std::map<int,std::string> gmap;

    std::string line;
    while (! is.eof())
    {
        std::getline(is, line);
        if (line.length() == 0) continue;

        std::vector<std::string> data;
        split_string(line,'\t', data);
        if (((data[0] == "#") && (data[1] == "index"))
            || (data[0] == "# Bam_Size_To_Use")
            || (data[0] == "# Bam_Orig_Path"))
        {
            continue;
        }

        if (data[0] == "# Bam_Path")
        {
            gmap[parse_int_str(data[HEAD_SRC_IDX])] = data[HEAD_NAME_IDX];
            continue;
        }
        // Get key string
        const int32_t key = parse_int_str(data[STAT_SOURCE_IDX]);

        // Make sure we have a BAM file source mapping!
        if (gmap.count(key) == 0)
        {
            log_os << "[ERROR]: While loading stats file unmapped data found for index: " << key
                   << ", on line: " << line << '\n';
            exit(EXIT_FAILURE);
        }

        //const ReadGroupStats rps(data);
        //setStats(gmap[key],rps);
    }
}



void
ReadGroupStatsSet::
write(std::ostream& os) const
{
    const unsigned n_groups(_group.size());
    for (unsigned i(0); i<n_groups; ++i)
    {
        os << "# Bam_Path\t" << i << "\t" << _group.get_key(i) << '\n';
    }
    // write column header for better readability
    os << "#\tindex"
       << "\tsample-count\tnumber-of-fragment-sizes"
       << "\treadOrientation"
       << '\n';

    for (unsigned i(0); i<n_groups; ++i)
    {
        os << i << '\t';
        getStats(i).write(os);
        os << '\n';
    }
}

