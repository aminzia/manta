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
///
///

#include "blt_util/bam_record_util.hh"



bool
is_innie_pair(
    const bam_record& bam_read)
{
    if (! bam_read.is_paired()) return false;
    if (bam_read.is_unmapped() || bam_read.is_mate_unmapped()) return false;
    if (bam_read.target_id() != bam_read.mate_target_id()) return false;

    if     (bam_read.pos() < bam_read.mate_pos())
    {
        if (! bam_read.is_fwd_strand()) return false;
        if (  bam_read.is_mate_fwd_strand()) return false;
    }
    else if (bam_read.pos() > bam_read.mate_pos())
    {
        if (  bam_read.is_fwd_strand()) return false;
        if (! bam_read.is_mate_fwd_strand()) return false;
    }
    else
    {
        if (bam_read.is_fwd_strand() == bam_read.is_mate_fwd_strand()) return false;
    }

    return true;
}



unsigned
get_avg_quality(
    const bam_record& bam_read)
{
    const unsigned len(bam_read.read_size());
    if (0 == len) return 0;

    const uint8_t* qual(bam_read.qual());
    unsigned sum(0);
    for (unsigned i(0); i<len; ++i)
    {
        sum+=qual[i];
    }
    // this does not capture the decimal remainder but well...
    return (sum/len);
}


