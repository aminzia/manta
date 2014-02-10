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

#include <iosfwd>
#include <string>
#include <vector>


/// \brief data pertaining to a de-novo assembly contig
///
/// stores for each contig the sequence and the number of reads
/// containing its seeding k-mer
///
struct Contig {

	Contig() :
		seq(""),
		avgCoverage(0),
		numKmers(0),
        seedReadCount(0)
	{}

	Contig(const Contig & c) :
		seq(c.seq),
		avgCoverage(c.avgCoverage),
		numKmers(c.numKmers),
        seedReadCount(c.seedReadCount)
	{}

	std::string seq;
	float avgCoverage;
	int numKmers;
	unsigned seedReadCount;
};

typedef std::vector<Contig> Assembly;

std::ostream& operator<<(std::ostream& os, const Contig& contig);


//typedef std::vector<AssembledContig> Assembly;

