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

#include "EstimateSVLoci.hh"
#include "ESLOptions.hh"
#include "SVLocusSetFinder.hh"

#include "blt_util/bam_header_util.hh"
#include "blt_util/align_path_bam_util.hh"
#include "blt_util/input_stream_handler.hh"
#include "blt_util/log.hh"

#include "common/OutStream.hh"

#include "manta/SVReferenceUtil.hh"

#include "boost/foreach.hpp"
#include "boost/shared_ptr.hpp"

#include <iostream>
#include <vector>



static
void
runESL(const ESLOptions& opt)
{
    {
        // early test that we have permission to write to output file
        OutStream outs(opt.outputFilename);
    }

    typedef boost::shared_ptr<bam_streamer> stream_ptr;
    std::vector<stream_ptr> bamStreams;

    // setup all data for main alignment loop:
    BOOST_FOREACH(const std::string& afile, opt.alignFileOpt.alignmentFilename)
    {
        stream_ptr tmp(new bam_streamer(afile.c_str(),opt.region.c_str()));
        bamStreams.push_back(tmp);
    }

    const unsigned bamCount(bamStreams.size());

    assert(0 != bamCount);

    // check bam header compatibility:
    if (bamCount > 1)
    {
        /// TODO: provide a better error exception for failed bam header check:
        const bam_header_t* compareHeader(bamStreams[0]->get_header());
        for (unsigned bamIndex(1); bamIndex<bamCount; ++bamIndex)
        {
            const bam_header_t* indexHeader(bamStreams[bamIndex]->get_header());
            if (! check_header_compatibility(compareHeader,indexHeader))
            {
                log_os << "ERROR: incompatible bam headers between files:\n"
                       << "\t" << opt.alignFileOpt.alignmentFilename[0] << "\n"
                       << "\t" << opt.alignFileOpt.alignmentFilename[bamIndex] << "\n";
                exit(EXIT_FAILURE);
            }
        }
    }

    // assume headers compatible after this point....

    const bam_header_t& header(*(bamStreams[0]->get_header()));
    const bam_header_info bamHeader(header);
    int32_t tid(0), beginPos(0), endPos(0);
    parse_bam_region(bamHeader,opt.region,tid,beginPos,endPos);

    const GenomeInterval scanRegion(tid,beginPos,endPos);

    reference_contig_segment refSegment;
    getIntervalReferenceSegment(opt.referenceFilename,bamHeader,scanRegion,refSegment);

    SVLocusSetFinder locusFinder(opt,scanRegion);
    locusFinder.setBamHeader(header);

    input_stream_data sdata;
    for (unsigned bamIndex(0); bamIndex<bamCount; ++bamIndex)
    {
        sdata.register_reads(*bamStreams[bamIndex],bamIndex);
    }

    // loop through alignments:
    input_stream_handler sinput(sdata);
    while (sinput.next())
    {
        const input_record_info current(sinput.get_current());

        if       (current.itype != INPUT_TYPE::READ)
        {
            log_os << "ERROR: invalid input condition.\n";
            exit(EXIT_FAILURE);
        }

        const bam_streamer& readStream(*bamStreams[current.sample_no]);
        const bam_record& read(*(readStream.get_record_ptr()));

        std::string ref;
        {
            ALIGNPATH::path_t apath;
            bam_cigar_to_apath(read.raw_cigar(), read.n_cigar(), apath);
        	const int alPos(read.pos()-scanRegion.range.begin_pos());
        	const int alLen(apath_ref_length(apath));
        	if (alPos < 0) continue;
            ref = refSegment.seq().substr((alPos-1),alLen);
        }

        locusFinder.update(read,current.sample_no,ref);
    }

    // finished updating:
    locusFinder.flush();

    locusFinder.getLocusSet().save(opt.outputFilename.c_str());
}



void
EstimateSVLoci::
runInternal(int argc, char* argv[]) const
{

    ESLOptions opt;

    parseESLOptions(*this,argc,argv,opt);
    runESL(opt);
}
