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
/// \author Ole Schulz-Trieglaff
/// \author Bret Barnes
///

#pragma once

#include "blt_util/bam_record.hh"
#include "blt_util/bam_record_util.hh"
#include "blt_util/circularCounter.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "manta/SVCandidate.hh"
#include "svgraph/SVLocus.hh"
#include "options/ReadScannerOptions.hh"

#include <string>
#include <vector>


namespace FragmentSizeType
{
    enum index_t
    {
        COMPRESSED,
        NORMAL,
        VERYCLOSE,
        CLOSE,
        DISTANT
    };
}


/// The counts in the SVLocus Graph represent an abstract weight of evidence supporting each edge/node.
///
/// To support large and small-scale evidence in a single graph, we need to allow for different weightings
/// for different evidence types
///
struct SVObservationWeights
{
    // noise reduction:
    static const unsigned observation = 3; ///< 'average' observation weight, this is used to scale noise filtration, but not for any evidence type

    // input evidence:
    static const unsigned readPair = observation;
    static const unsigned closeReadPair = 1;
    static const unsigned veryCloseReadPair = 1;
    static const unsigned internalReadEvent = observation; ///< indels, soft-clip, etc.
};


/// check bam record for soft-clipping which is interesting enough to be used as SV evidence:
///
/// \param[in] minQ
/// \param[in] minQFrac this fraction of bases must have qual>=minQ within the clipped region
///
void
getSVBreakendCandidateClip(
    const bam_record& bamRead,
    const ALIGNPATH::path_t& apath,
    unsigned& leadingClipLen,
    unsigned& trailingClipLen,
    const uint8_t minQ = 20,
    const float minQFrac = 0.75);


bool
isGoodShadow(
    const bam_record& bamRead,
    const uint8_t lastMapq,
    const std::string& lastQname,
    const double minSingletonMapq);


struct ReadScannerDerivOptions
{
    ReadScannerDerivOptions(const ReadScannerOptions& opt) :
        beforeBreakend(opt.minPairBreakendSize/2),
        afterBreakend(opt.minPairBreakendSize-beforeBreakend)
    {}

    const pos_t beforeBreakend;
    const pos_t afterBreakend;
};



/// consolidate functions which process a read to determine its
/// SV evidence value
///
/// In manta, evidence is scanned (at least) twice: once for SVLocus Graph generation
/// and then once again during hygen/scoring. We need to make sure both of these steps
/// are using the same logic to process read pairs into SV evidence. This class is
/// responsible for the shared logic
///
struct SVLocusScanner
{
    SVLocusScanner(
        const ReadScannerOptions& opt,
        const std::string& statsFilename,
        const std::vector<std::string>& alignmentFilename);

    /// this predicate runs any fast tests on the acceptability of a
    /// read for the SVLocus build
    /// Tests also for low mapq
    bool
    isReadFiltered(const bam_record& bamRead) const;

    unsigned
    getMinMapQ() const
    {
        return _opt.minMapq;
    }

    /// custom version of proper pair bit test:
    bool
    isProperPair(
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex) const;

    /// fragments sizes get thrown is serveral pre-defined categories:
    FragmentSizeType::index_t
    getFragmentSizeType(
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex) const;

    /// test whether a fragment is significantly larger than expected
    ///
    /// this function is useful to eliminate reads which fail the ProperPair test
    /// but are still very small
    ///
    bool
    isLargeFragment(
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex) const;

    /// test whether a fragment is significantly larger than expected and subsample borderline cases
    ///
    /// this behaves like isLargeFragment, except that fragments very close the non-amolous size are
    /// fitered unless they occur at high density
    ///
    bool
    isSampledLargeFragment(
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex,
        bool& isVeryClose,
        bool& isLeftMost) const;

    /// return true if the read is anomalous, for any anomaly type besides being a short innie read:
    bool
    isNonCompressedAnomalous(
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex) const;

    /// return true if the read is anomalous, for any anomaly type besides being a short innie fragment,
    /// with subsampling for very short fragments
    bool
    isSampledNonCompressedAnomalous(
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex,
        bool& isVeryCloseInnie,
        bool& isLeftMost) const;

    /// \brief is the read likely to indicate the presence of a small SV?
    ///
    /// this function flags reads which could contribute to a local small-variant assembly
    /// but would not otherwise be caught by the proper pair function
    ///
    /// "small" here is relative -- it means any event at a size where read pair evidence will not be dominant
    ///
    /// Note that the thresholds in this function are more stringent than the equivalent scan used to
    /// pick up reads prior to assembly -- in this case false positives could clog up the graph and
    /// interfere with larger event discovery if not kept under control
    bool
    isLocalAssemblyEvidence(
        const bam_record& bamRead,
        const reference_contig_segment& refSeq) const;

    /// return zero to many SVLocus objects if the read supports any
    /// structural variant(s) (detectable by manta)
    ///
    /// \param defaultReadGroupIndex the read group index to use in the absence of an RG tag
    /// (for now RGs are ignored for the purpose of gathering insert stats)
    ///
    void
    getSVLoci(
        const bam_record& bamRead,
        const unsigned defaultReadGroupIndex,
        const std::map<std::string, int32_t>& chromToIndex,
        const reference_contig_segment& refSeq,
        std::vector<SVLocus>& loci) const;

    /// get local and remote breakends for each SV Candidate which can be extracted from a read pair
    ///
    /// if remote read is not available, set remoteReadPtr to NULL and a best estimate will be generated for the remote breakend
    ///
    /// for all candidates, if one breakend is estimated from localRead and one is estimated from remoteRead, then
    /// the local breakend will be placed in candidate bp1 and the remote breakend will be placed in candidate.bp2
    ///
    void
    getBreakendPair(
        const bam_record& localRead,
        const bam_record* remoteReadPtr,
        const unsigned defaultReadGroupIndex,
        const std::map<std::string, int32_t>& chromToIndex,
        const reference_contig_segment& localRefSeq,
        const reference_contig_segment* remoteRefSeqPtr,
        std::vector<SVObservation>& candidates) const;

    /// provide direct access to the frag distro for
    /// functions which can't be cached
    ///
    const SizeDistribution&
    getFragSizeDistro(
        const unsigned defaultReadGroupIndex) const
    {
        return _rss.getStats(defaultReadGroupIndex).fragStats;
    }

    struct Range
    {
        Range() :
            min(0),
            max(0)
        {}

        double min;
        double max;
    };

    const Range&
    getEvidencePairRange(const unsigned readGroupIndex) const
    {
        return _stats[readGroupIndex].evidencePair;
    }

    struct CachedReadGroupStats
    {
        CachedReadGroupStats() :
            minCloseFragmentSize(0),
            minDistantFragmentSize(0),
            veryCloseFactor(0)
        {}

        /// fragment size range assumed for the purpose of creating SVLocusGraph regions
        Range breakendRegion;

        /// fragment size range used to determine if a read is anomalous
        Range properPair;

        Range evidencePair;

        int minCloseFragmentSize; ///< beyond the properPair anomalous threshold, there is a threshold to distinguish 'really-close' and 'close' pairs for the purpose of evidence weight
        int minDistantFragmentSize; ///< beyond the properPair anomalous threshold, there is a threshold to distinguish close and far pairs for the purpose of evidence weight

        float veryCloseFactor; ///< precomputed value used to scale down breakend size as fragments get smaller
    };

private:

    /////////////////////////////////////////////////
    // data:
    const ReadScannerOptions _opt;
    const ReadScannerDerivOptions _dopt;
    ReadGroupStatsSet _rss;

    std::vector<CachedReadGroupStats> _stats;

    mutable circularCounter _veryClosePairTracker;

//    std::string lastQname;
//    uint8_t lastMapq;
};

