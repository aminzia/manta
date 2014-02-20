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

#include "boost/test/unit_test.hpp"

#include "GlobalJumpAligner.hh"

#include "blt_util/align_path.hh"

#include <string>



BOOST_AUTO_TEST_SUITE( test_GlobalJumpAligner )

typedef short int score_t;

static const int badVal(-10000);

static
JumpAlignmentResult<score_t>
testAlignScores(
    const std::string& seq,
    const std::string& ref1,
    const std::string& ref2,
    int match, int mismatch, int open, int extend, int spliceOpen, int offEdge, int spliceOffEdge)
{
    AlignmentScores<score_t> scores(match, mismatch, open, extend, spliceOpen, offEdge, spliceOffEdge);
    score_t jumpScore(-3);
    GlobalJumpAligner<score_t> aligner(scores,jumpScore);
    JumpAlignmentResult<score_t> result;
    aligner.align(
        seq.begin(),seq.end(),
        ref1.begin(),ref1.end(),
        ref2.begin(),ref2.end(),
        result);

    return result;
}

static
JumpAlignmentResult<score_t>
testAlign(
    const std::string& seq,
    const std::string& ref1,
    const std::string& ref2)
{
    return testAlignScores(seq, ref1, ref2, 2,-4,-5,-1,badVal,-1,badVal);
}

static
JumpAlignmentResult<score_t>
testAlignSplice(
    const std::string& seq,
    const std::string& ref1,
    const std::string& ref2)
{
    return testAlignScores(seq, ref1, ref2, 2,-4,-5,-1,-4,-1,-1);
}

static
JumpAlignmentResult<score_t>
testAlign2(
    const std::string& seq,
    const std::string& ref1,
    const std::string& ref2)
{
    static const AlignmentScores<score_t> scores(2,-4,-10,-1,-1);
    static const int jumpScore(-20);
    GlobalJumpAligner<score_t> aligner(scores,jumpScore);
    JumpAlignmentResult<score_t> result;
    aligner.align(
        seq.begin(),seq.end(),
        ref1.begin(),ref1.end(),
        ref2.begin(),ref2.end(),
        result);

    return result;
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAligner0 )
{
    static const std::string seq("ABABACDCDC");
    static const std::string ref1("ABABA");
    static const std::string ref2("CDCDC");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,0);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,0);
}

BOOST_AUTO_TEST_CASE( test_GlobalJumpAligner1 )
{
    static const std::string seq("ABABACDCDC");
    static const std::string ref1("ABABAX");
    static const std::string ref2("CDCDC");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,0);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,0);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAligner2 )
{
    static const std::string seq("ABABACDCDC");
    static const std::string ref1("ABABA");
    static const std::string ref2("XCDCDC");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,0);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,1);
}

BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerLong )
{
    static const std::string seq("ABABACDCDC");
    static const std::string ref1("dslfjfkjaslABABAlsjfkdsflsk");
    static const std::string ref2("sdfldsklkjdCDCDCfsdlkjfslk");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,11);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,11);
}

BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerSplice )
{
    static const std::string seq("AAAAABBBBB");
    static const std::string ref1("xAAAAAGTxxxAGBBBBBx");
    static const std::string ref2("xxxx");

    JumpAlignmentResult<score_t> result = testAlignSplice(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5=7N5=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,1);
}
BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerSpliceRef2 )
{
    static const std::string seq("AAAAABBBBB");
    static const std::string ref1("xxxx");
    static const std::string ref2("xAAAAAGTxxxAGBBBBBx");

    JumpAlignmentResult<score_t> result = testAlignSplice(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5=7N5=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,1);
}
BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerSpliceNoSplice )
{
    static const std::string seq("AAAAABBBBB");
    static const std::string ref1("xAAAAAGGxxxAGBBBBBx");
    static const std::string ref2("xxxx");

    JumpAlignmentResult<score_t> result = testAlignSplice(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5=7D5=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,1);
}

BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerSpliceNoSplice2 )
{
    static const std::string seq("AAAAABBBBB");
    static const std::string ref1("xAAAAAGTxxxGGBBBBBx");
    static const std::string ref2("xxxx");

    JumpAlignmentResult<score_t> result = testAlignSplice(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5=7D5=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,1);
}

BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerSimpleIndels )
{
    static const std::string seq("ABABAABABACDCDCDyCDCDC");
    static const std::string ref1("xABABABABABAx");
    static const std::string ref2("xCDCDCDCDCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5=1D5=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,1);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"6=1I5=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,1);
}



BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerBPInsert )
{
    static const std::string seq("ABABABABABA1234CDCDCDCDCDC");
    static const std::string ref1("xABABABABABAx");
    static const std::string ref2("xCDCDCDCDCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"11=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,1);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"11=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,1);
    BOOST_REQUIRE_EQUAL(result.jumpInsertSize,4u);
}


// define behavior when the breakpoint solution is a 1d range
BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerBPRange )
{
    static const std::string seq("xyzxyzxyzABCABCABCxyzxyzxyz");
    static const std::string ref1("xyzxyzxyzxyzABCABCstustu");
    static const std::string ref2("stustuABCABCxyzxyzxyzxyz");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"12=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,3);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"15=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,6);
    BOOST_REQUIRE_EQUAL(result.jumpInsertSize,0u);
    BOOST_REQUIRE_EQUAL(result.jumpRange,3u);
}


// define behavior when the breakpoint solution is a 1d range
BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerBPRange2 )
{
    static const std::string seq("xyzxyzxyzABCABCABCxyzxyzxyz");
    static const std::string ref1("xyzxyzxyzxyzABCABCABCABCABC");
    static const std::string ref2("ABCABCABCABCABCxyzxyzxyzxyz");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"9=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,3);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"18=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,6);
    BOOST_REQUIRE_EQUAL(result.jumpInsertSize,0u);
    BOOST_REQUIRE_EQUAL(result.jumpRange,9u);
}


// define behavior when the breakpoint solution has an insertion with a repeat
BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerInsert )
{
    static const std::string seq("xyzxyzxyzABCABCABCABCABCABCxyzxyzxyz");
    static const std::string ref1("xyzxyzxyzxyzABCABCstustu");
    static const std::string ref2("stustuABCABCxyzxyzxyzxyz");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"15=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,3);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"15=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,6);
    BOOST_REQUIRE_EQUAL(result.jumpInsertSize,6u);
    BOOST_REQUIRE_EQUAL(result.jumpRange,0u);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerOnly1 )
{
    static const std::string seq("ABABA");
    static const std::string ref1("xABABAx");
    static const std::string ref2("xCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,1);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,0);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerOnly2 )
{
    static const std::string seq("CDCDC");
    static const std::string ref1("xABABAx");
    static const std::string ref2("xCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,0);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,1);
}



BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerOffEdge )
{
    static const std::string seq("123456ABABACDCDC123456");
    static const std::string ref1("xABABAx");
    static const std::string ref2("xCDCDCx");

    JumpAlignmentResult<score_t> result = testAlign(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"5S1X5=");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,0);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"5=1X5S");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,1);
}



BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerRef2Clip )
{
    // extracted from production failure case:
    //
    static const std::string seq("GGCAGAAAAGGAAATA");
    static const std::string ref1("TAAAAAGTAGAT");
    static const std::string ref2("AAAGGAAATA");

    JumpAlignmentResult<score_t> result = testAlign2(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,0);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"6S10=");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,0);
    BOOST_REQUIRE_EQUAL(result.jumpInsertSize,0u);
    BOOST_REQUIRE_EQUAL(result.jumpRange,0u);
}


BOOST_AUTO_TEST_CASE( test_GlobalJumpAlignerRef1Clip )
{
    // extracted from production failure case:
    //
    static const std::string seq("TAAAAAGTAGATTTCGT");
    static const std::string ref1("TAAAAAGTAGAT");
    static const std::string ref2("AAAGGAAATA");

    JumpAlignmentResult<score_t> result = testAlign2(seq,ref1,ref2);

    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align1.apath),"12=5S");
    BOOST_REQUIRE_EQUAL(result.align1.beginPos,0);
    BOOST_REQUIRE_EQUAL(apath_to_cigar(result.align2.apath),"");
    BOOST_REQUIRE_EQUAL(result.align2.beginPos,0);
    BOOST_REQUIRE_EQUAL(result.jumpInsertSize,0u);
    BOOST_REQUIRE_EQUAL(result.jumpRange,0u);
}


BOOST_AUTO_TEST_SUITE_END()

