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

#include "format/VcfWriterSV.hh"

#include "blt_util/samtools_fasta_util.hh"
#include "blt_util/string_util.hh"
#include "blt_util/seq_util.hh"
#include "blt_util/vcf_util.hh"
#include "common/Exceptions.hh"
#include "manta/SVCandidateUtil.hh"

#include <iostream>
#include <sstream>


//#define DEBUG_VCF


#ifdef DEBUG_VCF
#include "blt_util/log.hh"
#endif



VcfWriterSV::
VcfWriterSV(
    const std::string& referenceFilename,
    const SVLocusSet& set,
    std::ostream& os) :
    _referenceFilename(referenceFilename),
    _header(set.header),
    _os(os)
{
}



void
VcfWriterSV::
writeHeaderPrefix(
    const char* progName,
    const char* progVersion)
{
    _os << "##fileformat=VCFv4.1\n";
    _os << "##fileData=" << vcf_fileDate << "\n";
    _os << "##source=" << progName << " " << progVersion << "\n";
    _os << "##reference=file://" << _referenceFilename << "\n";

    BOOST_FOREACH(const bam_header_info::chrom_info& cdata, _header.chrom_data)
    {
        _os << "##contig=<ID=" << cdata.label << ",length=" << cdata.length << ">\n";
    }

    /// vcf 4.1 reserved/suggested INFO tags:
    _os << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n";
    _os << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    _os << "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
    _os << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n";
    _os << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">\n";
    _os << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">\n";
    _os << "##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"CIGAR alignment for each alternate indel allele\">\n";
    _os << "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakend\">\n";
    _os << "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"ID of event associated to breakend\">\n";
    _os << "##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical homology at event breakpoints\">\n";
    _os << "##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical homology at event breakpoints\">\n";

    /// custom INFO tags:
    _os << "##INFO=<ID=SVINSLEN,Number=.,Type=Integer,Description=\"Length of insertion\">\n";
    _os << "##INFO=<ID=SVINSSEQ,Number=.,Type=String,Description=\"Sequence of insertion\">\n";
    _os << "##INFO=<ID=LEFT_SVINSSEQ,Number=.,Type=String,Description=\"Known left side of insertion for an insertion of unknown length\">\n";
    _os << "##INFO=<ID=RIGHT_SVINSSEQ,Number=.,Type=String,Description=\"Known right side of insertion for an insertion of unknown length\">\n";
    _os << "##INFO=<ID=PAIR_COUNT,Number=1,Type=Integer,Description=\"Read pairs supporting this variant where both reads are confidently mapped\">\n";
    _os << "##INFO=<ID=BND_PAIR_COUNT,Number=1,Type=Integer,Description=\"Confidently mapped reads supporting this variant at this breakend (mapping may not be confident at remote breakend)\">\n";
    _os << "##INFO=<ID=UPSTREAM_PAIR_COUNT,Number=1,Type=Integer,Description=\"Confidently mapped reads supporting this variant at the upstream breakend (mapping may not be confident at downstream breakend)\">\n";
    _os << "##INFO=<ID=DOWNSTREAM_PAIR_COUNT,Number=1,Type=Integer,Description=\"Confidently mapped reads supporting this variant at this downstream breakend (mapping may not be confident at upstream breakend)\">\n";
    _os << "##INFO=<ID=INV3,Number=0,Type=Flag,Description=\"Inversion breakends open 3' of reported location\">\n";
    _os << "##INFO=<ID=INV5,Number=0,Type=Flag,Description=\"Inversion breakends open 5' of reported location\">\n";

    addHeaderInfo();

    addHeaderFormat();

    addHeaderFilters();

    _os << "##ALT=<ID=BND,Description=\"Translocation Breakend\">\n";
    _os << "##ALT=<ID=INV,Description=\"Inversion\">\n";
    _os << "##ALT=<ID=DEL,Description=\"Deletion\">\n";
    _os << "##ALT=<ID=INS,Description=\"Insertion\">\n";
    _os << "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">\n";
    _os << "##ALT=<ID=COMPLEX,Description=\"Unknown Candidate Type\">\n";
}


static
void
writeHeaderColKeyPrefix(std::ostream& os)
{
    os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
}



void
VcfWriterSV::
writeHeaderColumnKey()
{
    writeHeaderColKeyPrefix(_os);
    addHeaderFormatSampleKey();
    _os << '\n';
}



static
void
makeInfoField(
    const VcfWriterSV::InfoTag_t& info,
    std::ostream& os)
{
    static const char sep(';');
    bool isFirst(true);
    BOOST_FOREACH(const std::string& is, info)
    {
        if (! isFirst) os << sep;
        else           isFirst = false;
        os << is;
    }
}



static
void
makeFormatSampleField(
    const VcfWriterSV::SampleTag_t& sample,
    std::ostream& os)
{
    static const char sep(':');

    if (sample.empty()) return;

    {
        // first write FORMAT field:
        os << '\t';

        bool isFirst(true);
        BOOST_FOREACH(const VcfWriterSV::SampleTag_t::value_type& fs, sample)
        {
            if (! isFirst) os << sep;
            else           isFirst = false;

            assert(! fs.first.empty());
            os << fs.first;
        }
    }

    unsigned nSamples(0);
    BOOST_FOREACH(const VcfWriterSV::SampleTag_t::value_type& fs, sample)
    {
        const unsigned ns(fs.second.size());
        nSamples = std::max(nSamples, ns);
    }

    for (unsigned sampleIndex(0); sampleIndex < nSamples; ++sampleIndex)
    {
        os << '\t';

        // next write SAMPLE field:
        {
            bool isFirst(true);
            BOOST_FOREACH(const VcfWriterSV::SampleTag_t::value_type& fs, sample)
            {
                if (! isFirst) os << sep;
                else           isFirst = false;

                if (fs.second.size() <= sampleIndex)
                {
                    os << '.';
                }
                else if (fs.second[sampleIndex].empty())
                {
                    os << '.';
                }
                else
                {
                    os << fs.second[sampleIndex];
                }
            }
        }
    }
}



#ifdef DEBUG_VCF

static
void
addDebugInfo(
    const SVBreakend& bp1,
    const SVBreakend& bp2,
    const bool isFirstOfPair,
    const SVCandidateAssemblyData& assemblyData,
    InfoTag_t& infotags)
{
    if (! isFirstOfPair) return;

    // store alignment start + cigar string for each section of the jumping alignment.
    // there can be several contigs per breakend, so we iterate over all of them.
    // only the first breakpoint gets the alignments attached to its VCF entry

    if (assemblyData.isSpanning)
    {
        const unsigned numAlign(assemblyData.spanningAlignments.size());
        std::string cigar1;
        std::string cigar2;
        for (unsigned alignIndex(0); alignIndex<numAlign; ++alignIndex)
        {
            const SVCandidateAssemblyData::JumpAlignmentResultType align(assemblyData.spanningAlignments[alignIndex]);
            infotags.push_back( str(boost::format("CTG_JALIGN_%i_POS_A=%d") %
                                    alignIndex %
                                    (bp1.interval.range.begin_pos()+align.align1.beginPos)) );
            infotags.push_back( str(boost::format("CTG_JALIGN_%i_POS_B=%d") %
                                    alignIndex %
                                    (bp2.interval.range.begin_pos()+align.align2.beginPos)) );

            apath_to_cigar(align.align1.apath,cigar1);
            apath_to_cigar(align.align2.apath,cigar2);

            infotags.push_back( str(boost::format("CTG_JALIGN_%i_CIGAR_A=%s") % alignIndex % cigar1) );
            infotags.push_back( str(boost::format("CTG_JALIGN_%i_CIGAR_B=%s") % alignIndex % cigar2) );
        }
    }
}

#endif



static
void
addSharedInfo(
    const EventInfo& event,
    VcfWriterSV::InfoTag_t& infoTags)
{
    if (event.isEvent())
    {
        infoTags.push_back( str(boost::format("EVENT=%i") % event.label));
    }
}



void
VcfWriterSV::
writeTransloc(
    const SVCandidate& sv,
    const SVId& svId,
    const bool isFirstBreakend,
    const SVCandidateSetData& /*svData*/,
    const SVCandidateAssemblyData& /*adata*/,
    const EventInfo& event)
{
    const bool isImprecise(sv.isImprecise());
    const bool isBreakendRangeSameShift(sv.isBreakendRangeSameShift());

    const SVBreakend& bpA( isFirstBreakend ? sv.bp1 : sv.bp2);
    const SVBreakend& bpB( isFirstBreakend ? sv.bp2 : sv.bp1);

    InfoTag_t infotags;
    SampleTag_t sampletags;

    // get CHROM
    const std::string& chrom(_header.chrom_data[bpA.interval.tid].label);
    const std::string& mateChrom(_header.chrom_data[bpB.interval.tid].label);

    const known_pos_range2& bpArange(bpA.interval.range);
    const known_pos_range2& bpBrange(bpB.interval.range);

    if (! isImprecise)
    {
        assert(bpArange.size() == bpBrange.size());
    }

    // get POS
    pos_t pos(bpArange.center_pos()+1);
    pos_t matePos(bpBrange.center_pos()+1);
    if (! isImprecise)
    {
        pos = bpArange.begin_pos()+1;
        if (isBreakendRangeSameShift)
        {
            matePos = bpBrange.begin_pos()+1;
        }
        else
        {
            matePos = bpBrange.end_pos();
        }
    }

    // get ID
    const std::string& localId(isFirstBreakend ? svId.localId : svId.mateId);
    const std::string& mateId(isFirstBreakend ? svId.mateId : svId.localId);

    // get REF
    std::string ref;
    get_standardized_region_seq(_referenceFilename,chrom,pos-1,pos-1,ref);

    assert(1 == ref.size());

    const bool isReverseInsertSeq(! (isFirstBreakend || (bpA.state != bpB.state)));
    std::string tmpString;
    const std::string* insertSeqPtr(&sv.insertSeq);
    if (isReverseInsertSeq)
    {
        tmpString = reverseCompCopyStr(sv.insertSeq);
        insertSeqPtr = &tmpString;
    }
    const std::string& insertSeq(*insertSeqPtr);

    // build alt:
    boost::format altFormat("%4%%3%%1%:%2%%3%%5%");
    {
        std::string altPrefix;
        std::string altSuffix;
        if     (bpA.state == SVBreakendState::RIGHT_OPEN)
        {
            altPrefix = ref + insertSeq;
        }
        else if (bpA.state == SVBreakendState::LEFT_OPEN)
        {
            altSuffix = insertSeq + ref;
        }
        else
        {
            assert(false && "Unexpected bpA.state");
        }


        char altSep('?');
        if     (bpB.state == SVBreakendState::RIGHT_OPEN)
        {
            altSep=']';
        }
        else if (bpB.state == SVBreakendState::LEFT_OPEN)
        {
            altSep='[';
        }
        else
        {
            assert(false && "Unexpected bpB.state");
        }

        altFormat % mateChrom % matePos % altSep % altPrefix % altSuffix;
    }

    // build INFO field
    infotags.push_back("SVTYPE=BND");
    infotags.push_back("MATEID="+mateId);
    infotags.push_back( str(boost::format("BND_PAIR_COUNT=%i") % bpA.getLocalPairCount()) );
    infotags.push_back( str(boost::format("PAIR_COUNT=%i") % bpA.getPairCount()) );
    if (isImprecise)
    {
        infotags.push_back("IMPRECISE");
    }

    if (bpArange.size() > 1)
    {
        infotags.push_back( str( boost::format("CIPOS=%i,%i") % ((bpArange.begin_pos()+1) - pos) % (bpArange.end_pos() - pos) ));
    }

    if (! isImprecise)
    {
        if (bpArange.size() > 1)
        {
            infotags.push_back( str( boost::format("HOMLEN=%i") % (bpArange.size()-1) ));
            std::string homref;
            get_standardized_region_seq(_referenceFilename,chrom,bpArange.begin_pos()+1,bpArange.end_pos()-1,homref);
            infotags.push_back( str( boost::format("HOMSEQ=%s") % (homref) ));
        }
    }

    if (! insertSeq.empty())
    {
        infotags.push_back( str( boost::format("SVINSLEN=%i") % (insertSeq.size()) ));
        infotags.push_back( str( boost::format("SVINSSEQ=%s") % (insertSeq) ));
    }

    addSharedInfo(event, infotags);

    modifyInfo(event, infotags);
    modifyTranslocInfo(isFirstBreakend, infotags);

    modifySample(sv, sampletags);

#ifdef DEBUG_VCF
    addDebugInfo(bpA, bpB, isFirstBreakend, adata, infotags);
#endif

    // write out record:
    _os << chrom
        << '\t' << pos
        << '\t' << localId // ID
        << '\t' << ref // REF
        << '\t' << str( altFormat ) // ALT
        << '\t';
    writeQual();
    _os << '\t';
    writeFilter();
    _os << '\t';
    makeInfoField(infotags,_os); // INFO
    makeFormatSampleField(sampletags, _os); // FORMAT + SAMPLE
    _os << '\n';
}



void
VcfWriterSV::
writeTranslocPair(
    const SVCandidate& sv,
    const SVId& svId,
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& adata,
    const EventInfo& event)
{
    writeTransloc(sv, svId, true, svData, adata, event);
    writeTransloc(sv, svId, false, svData, adata, event);
}



void
VcfWriterSV::
writeInvdel(
    const SVCandidate& sv,
    const SVId& svId,
    const SVCandidateAssemblyData& /*adata*/,
    const bool isIndel,
    const EventInfo& event)
{
    const bool isImprecise(sv.isImprecise());
    const bool isBreakendRangeSameShift(sv.isBreakendRangeSameShift());

    const bool isBp1First(sv.bp1.interval.range.begin_pos()<=sv.bp2.interval.range.begin_pos());

    const SVBreakend& bpA(isBp1First ? sv.bp1 : sv.bp2);
    const SVBreakend& bpB(isBp1First ? sv.bp2 : sv.bp1);

    InfoTag_t infoTags;
    SampleTag_t sampleTags;

    // get CHROM
    const std::string& chrom(_header.chrom_data[sv.bp1.interval.tid].label);

    const known_pos_range2& bpArange(bpA.interval.range);
    const known_pos_range2& bpBrange(bpB.interval.range);

    if (! isImprecise)
    {
        assert(bpArange.size() == bpBrange.size());
    }

    // above this size all records use symbolic alleles (ie. <DEL>):
    static const unsigned maxNonSymbolicRecordSize(1000);

    // if the variant is a combination of simple insertion and deletions, and below
    // a large-event size threshold, it is classified as a small variant. In this case
    // we report the event using full REF and ALT sequences, plus a CIGAR string for
    // complex in/del combinations
    //
    bool isSmallVariant(false);
    if ((! isImprecise) && isIndel && (! sv.isUnknownSizeInsertion))
    {
        const unsigned deleteSize(bpBrange.begin_pos() - bpArange.begin_pos());
        const unsigned insertSize(sv.insertSeq.size());

        const bool isSmallDelete(deleteSize<=maxNonSymbolicRecordSize);
        const bool isSmallInsert(insertSize<=maxNonSymbolicRecordSize);

        isSmallVariant = (isSmallDelete && isSmallInsert);
    }

    // get POS and endPos
    pos_t pos(bpArange.center_pos()+1);
    pos_t endPos(bpBrange.center_pos());
    if (! isImprecise)
    {
        pos = bpArange.begin_pos()+1;
        if (isBreakendRangeSameShift)
        {
            endPos = bpBrange.begin_pos();
        }
        else
        {
            endPos = bpBrange.end_pos()-1;
        }
    }
    else
    {
        /// check against the rare case arising when CIEND is a subset of CIPOS:
        endPos=std::max(endPos,pos+1);
    }

    if (pos<1) return;

    // get REF
    std::string ref;
    {
        const pos_t beginRefPos(pos-1);
        pos_t endRefPos(beginRefPos);
        if (isSmallVariant) endRefPos=endPos-1;

        get_standardized_region_seq(_referenceFilename, chrom, beginRefPos, endRefPos, ref);

        assert(static_cast<unsigned>(1+endRefPos-beginRefPos) == ref.size());
    }

    // build alt:
    std::string alt;
    if (isSmallVariant)
    {
        alt = ref[0] + sv.insertSeq;
    }
    else
    {
        alt = str( boost::format("<%s>") % svId.getLabel());
    }

    // build INFO field
    std::vector<std::string> words;
    split_string(svId.getLabel(),':',words);
    {
        // note that there's a reasonable argument for displaying these tags only when a
        // symbolic allele is used (by a strict reading of the vcf spec) -- we instead
        // print these fields for all variants for uniformity within the manta vcf:
        //
        infoTags.push_back( str(boost::format("END=%i") % endPos));
        infoTags.push_back( str(boost::format("SVTYPE=%s") % words[0]));
        const pos_t refLen(endPos-pos);
        pos_t svLen(refLen);

        if (! sv.isUnknownSizeInsertion)
        {
            if (isIndel)
            {
                const pos_t insertLen(static_cast<pos_t>(sv.insertSeq.size()));
                if ( insertLen > refLen )
                {
                    svLen = insertLen;
                }
                else
                {
                    svLen = -refLen;
                }
            }
            infoTags.push_back( str(boost::format("SVLEN=%i") % (svLen)));
        }
    }
    infoTags.push_back( str(boost::format("UPSTREAM_PAIR_COUNT=%i") % bpA.getLocalPairCount()) );
    infoTags.push_back( str(boost::format("DOWNSTREAM_PAIR_COUNT=%i") % bpB.getLocalPairCount()) );
    infoTags.push_back( str(boost::format("PAIR_COUNT=%i") % bpA.getPairCount()) );

    if (isSmallVariant)
    {
        if (! sv.insertAlignment.empty())
        {
            std::string cigar;
            apath_to_cigar(sv.insertAlignment,cigar);

            // add the 1M to signify the leading reference base:
            infoTags.push_back( str(boost::format("CIGAR=1M%s") % cigar));
        }
    }

    if (isImprecise)
    {
        infoTags.push_back("IMPRECISE");
    }

    if (bpArange.size() > 1)
    {
        infoTags.push_back( str( boost::format("CIPOS=%i,%i") % ((bpArange.begin_pos()+1) - pos) % (bpArange.end_pos() - pos) ));
    }

    if (! isSmallVariant)
    {
        if (bpBrange.size() > 1)
        {
            infoTags.push_back( str( boost::format("CIEND=%i,%i") % (bpBrange.begin_pos() - endPos) % ((bpBrange.end_pos()-1) - endPos) ));
        }
    }

    if (! isImprecise)
    {
        if (bpArange.size() > 1)
        {
            infoTags.push_back( str( boost::format("HOMLEN=%i") % (bpArange.size()-1) ));
            std::string homref;
            get_standardized_region_seq(_referenceFilename,chrom,bpArange.begin_pos()+1,bpArange.end_pos()-1,homref);
            infoTags.push_back( str( boost::format("HOMSEQ=%s") % (homref) ));
        }
    }

    if (! isSmallVariant)
    {
        if (! (sv.insertSeq.empty() || sv.isUnknownSizeInsertion))
        {
            infoTags.push_back( str( boost::format("SVINSLEN=%i") % (sv.insertSeq.size()) ));
            if (isBp1First || (bpA.state != bpB.state))
            {
                infoTags.push_back( str( boost::format("SVINSSEQ=%s") % (sv.insertSeq) ));
            }
            else
            {
                infoTags.push_back( str( boost::format("SVINSSEQ=%s") % reverseCompCopyStr(sv.insertSeq) ));
            }
        }
    }

    if (sv.isUnknownSizeInsertion)
    {
        if (! sv.unknownSizeInsertionLeftSeq.empty())
        {
            infoTags.push_back( str( boost::format("LEFT_SVINSSEQ=%s") % (sv.unknownSizeInsertionLeftSeq) ));
        }

        if (! sv.unknownSizeInsertionRightSeq.empty())
        {
            infoTags.push_back( str( boost::format("RIGHT_SVINSSEQ=%s") % (sv.unknownSizeInsertionRightSeq) ));
        }
    }

    if (svId.svType == EXTENDED_SV_TYPE::INVERSION)
    {
        if (sv.bp1.state == SVBreakendState::RIGHT_OPEN)
        {
            infoTags.push_back("INV3");
        }
        else if (sv.bp1.state == SVBreakendState::LEFT_OPEN)
        {
            infoTags.push_back("INV5");
        }
        else
        {
            assert(false && "Unexpected inversion configuration");
        }
    }

    addSharedInfo(event, infoTags);

    modifyInfo(event, infoTags);
    modifySample(sv, sampleTags);

    // write out record:
    _os << chrom
        << '\t' << pos
        << '\t' << svId.localId // ID
        << '\t' << ref // REF
        << '\t' << alt // ALT
        << '\t';
    writeQual();
    _os << '\t';
    writeFilter();
    _os << '\t';
    makeInfoField(infoTags,_os); // INFO
    makeFormatSampleField(sampleTags, _os); // FORMAT + SAMPLE
    _os << '\n';
}



static
bool
isAcceptedSVType(
    const EXTENDED_SV_TYPE::index_t svType)
{
    using namespace EXTENDED_SV_TYPE;

    switch (svType)
    {
    case INTERTRANSLOC:
    case INTRATRANSLOC:
    case INVERSION:
    case INSERT:
    case DELETE:
    case TANDUP:
        return true;
    default:
        return false;
    }
}



void
VcfWriterSV::
writeSVCore(
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& adata,
    const SVCandidate& sv,
    const SVId& svId,
    const EventInfo& event)
{
    using namespace EXTENDED_SV_TYPE;
    const index_t svType(getExtendedSVType(sv));

#ifdef DEBUG_VCF
    log_os << "VcfWriterSV::writeSVCore svType: " << SV_TYPE::label(svType) << "\n";
#endif

    if (! isAcceptedSVType(svType))
    {
        using namespace illumina::common;

        std::ostringstream oss;
        oss << "ERROR: sv candidate cannot be classified: " << sv << "\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    if      (isSVTransloc(svType))
    {
        writeTranslocPair(sv, svId, svData, adata, event);
    }
    else
    {
        const bool isIndel(isSVIndel(svType));
        writeInvdel(sv, svId, adata, isIndel, event);
    }
}



void
VcfWriterSV::
writeFilters(
    const std::set<std::string>& filters) const
{
    if (filters.empty())
    {
        _os << "PASS";
    }
    else
    {
        bool isFirst(true);
        BOOST_FOREACH(const std::string& filter, filters)
        {
            if (isFirst)
            {
                isFirst=false;
            }
            else
            {
                _os << ';';
            }
            _os << filter;
        }
    }
}


