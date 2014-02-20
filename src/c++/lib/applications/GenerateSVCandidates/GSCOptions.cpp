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

#include "GSCOptions.hh"
#include "EdgeOptionsParser.hh"

#include "blt_util/log.hh"
#include "options/AlignmentFileOptionsParser.hh"
#include "options/ReadScannerOptionsParser.hh"
#include "options/optionsUtil.hh"

#include "boost/foreach.hpp"
#include "boost/program_options.hpp"

#include <iostream>



static
void
usage(
    std::ostream& os,
    const manta::Program& prog,
    const boost::program_options::options_description& visible,
    const char* msg = NULL)
{
    os << "\n" << prog.name() << ": call candidates from an SV Locus graph\n\n";
    os << "version: " << prog.version() << "\n\n";
    os << "usage: " << prog.name() << " [options]\n\n";
    os << visible << "\n\n";

    if (NULL != msg)
    {
        os << msg << "\n\n";
    }
    exit(2);
}



static
void
checkStandardizeUsageFile(
    std::ostream& os,
    const manta::Program& prog,
    const boost::program_options::options_description& visible,
    std::string& filename,
    const char* fileLabel)
{
    std::string errorMsg;
    if ( checkStandardizeInputFile(filename, fileLabel, errorMsg))
    {
        usage(os,prog,visible,errorMsg.c_str());
    }
}



void
parseGSCOptions(
    const manta::Program& prog,
    int argc, char* argv[],
    GSCOptions& opt)
{
    namespace po = boost::program_options;
    po::options_description req("configuration");
    req.add_options()
    ("graph-file", po::value(&opt.graphFilename),
     "sv locus graph file (required)")
    ("align-stats", po::value(&opt.statsFilename),
     "pre-computed alignment statistics for the input alignment files (required)")
    ("chrom-depth", po::value(&opt.chromDepthFilename),
     "average depth estimate for each chromosome")
    ("ref", po::value(&opt.referenceFilename),
     "fasta reference sequence (required)")
    ("truth-vcf", po::value(&opt.truthVcfFilename),
     "optional truth VCF file (for testing)")
    ("edge-runtime-log", po::value(&opt.edgeRuntimeFilename),
     "optionally log edge runtime to this file")
    ("candidate-output-file", po::value(&opt.candidateOutputFilename),
     "Write SV candidates to file (required)")
    ("diploid-output-file", po::value(&opt.diploidOutputFilename),
     "Write germline diploid SVs to file (at least one  non-tumor alignment file must be specified)")
    ("somatic-output-file", po::value(&opt.somaticOutputFilename),
     "Write somatic SVs to file (at least one tumor and non-tumor alignment file must be specified)")
    ("verbose", po::value(&opt.isVerbose)->zero_tokens(),
     "Turn on low-detail INFO logging.")
    ("skip-assembly", po::value(&opt.isSkipAssembly)->zero_tokens(),
     "Turn off all breakend and small-variant assembly. Only large, imprecise variants will be reported based on anomalous read pairs.")
    ("rna", po::value(&opt.isRNA)->zero_tokens(),
     "For RNA input. Skip small deletions and modify diploid scoring.")
    ("min-scored-sv-size", po::value(&opt.minScoredVariantSize)->default_value(opt.minScoredVariantSize),
     "minimum size for variants which are scored and output following initial candidate generation")
    ;

    po::options_description alignDesc(getOptionsDescription(opt.alignFileOpt));
    po::options_description edgeDesc(getOptionsDescription(opt.edgeOpt));
    po::options_description scanDesc(getOptionsDescription(opt.scanOpt));
    po::options_description diploidCallDesc(getOptionsDescription(opt.diploidOpt));
    po::options_description somaticCallDesc(getOptionsDescription(opt.somaticOpt));

    po::options_description help("help");
    help.add_options()
    ("help,h","print this message");

    po::options_description visible("options");
    visible.add(alignDesc).add(scanDesc).add(req).add(edgeDesc).add(diploidCallDesc).add(somaticCallDesc).add(help);

    bool po_parse_fail(false);
    po::variables_map vm;
    try
    {
        po::store(po::parse_command_line(argc, argv, visible,
                                         po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
        po::notify(vm);
    }
    catch (const boost::program_options::error& e)     // todo:: find out what is the more specific exception class thrown by program options
    {
        log_os << "\nERROR: Exception thrown by option parser: " << e.what() << "\n";
        po_parse_fail=true;
    }

    if ((argc<=1) || (vm.count("help")) || po_parse_fail)
    {
        usage(log_os,prog,visible);
    }

    std::string errorMsg;
    if (parseOptions(vm, opt.edgeOpt, errorMsg))
    {}
    else if (parseOptions(vm, opt.alignFileOpt, errorMsg))
    {}
    else if (parseOptions(vm, opt.scanOpt, errorMsg))
    {}
    else if (opt.statsFilename.empty())
    {
        errorMsg="Need the alignment stats file";
    }
    else if (opt.referenceFilename.empty())
    {
        errorMsg="Need the FASTA reference file";
    }

    if (! errorMsg.empty()) usage(log_os,prog,visible,errorMsg.c_str());

    checkStandardizeUsageFile(log_os,prog,visible,opt.graphFilename,"SV locus graph");
    checkStandardizeUsageFile(log_os,prog,visible,opt.referenceFilename,"reference fasta");
    checkStandardizeUsageFile(log_os,prog,visible,opt.statsFilename,"alignment statistics");

    if (! opt.chromDepthFilename.empty())
    {
        checkStandardizeUsageFile(log_os,prog,visible,opt.chromDepthFilename,"chromosome depth");
    }
    if (opt.candidateOutputFilename.empty())
    {
        usage(log_os,prog,visible,"Must specify candidate output file");
    }

    if (opt.diploidOutputFilename.empty())
    {
        usage(log_os,prog,visible,"Must specify diploid output file");
    }

    if (! opt.somaticOutputFilename.empty())
    {
        unsigned normalCount(0);
        unsigned tumorCount(0);
        BOOST_FOREACH(const bool value, opt.alignFileOpt.isAlignmentTumor)
        {
            if (value) tumorCount++;
            else      normalCount++;
        }
        if ((normalCount==0) || (tumorCount==0))
        {
            usage(log_os,prog,visible,"Must specify at least one tumor and non-tumor alignment file for somatic output");
        }
    }
}

