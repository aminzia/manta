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

#pragma once

#include "blt_util/bam_record.hh"


/// is this read part of mapped pair with 'Innie' orientation?
///
/// Note this does not test MAPQ or fragment size, but could
/// be used as the core of a 'proper-pair' predicate
bool
is_innie_pair(
    const bam_record& bam_read);

/// return average basecall qscore for this read
unsigned
get_avg_quality(
    const bam_record& bam_read);


bool
has_large_gap(
	const bam_record& bam_read);
