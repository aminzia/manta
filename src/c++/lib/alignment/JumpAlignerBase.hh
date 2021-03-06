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

#pragma once

#include "AlignerBase.hh"
#include "AlignerUtil.hh"
#include "Alignment.hh"

#include "blt_util/basic_matrix.hh"


/// represents alignment of a query sequence which can switch over from reference1 to reference2
///
/// an empty alignment to one reference indicates that the entire alignment is to the other reference
///
template <typename ScoreType>
struct JumpAlignmentResult
{
    JumpAlignmentResult()
    {
        clear();
    }

    void
    clear()
    {
        score=0;
        jumpInsertSize=0;
        jumpRange=0;
        align1.clear();
        align2.clear();
    }


    ScoreType score;
    unsigned jumpInsertSize;
    unsigned jumpRange; ///< length of sequence over which jump would have the same score (left-most on align1 is reported)
    Alignment align1;
    Alignment align2;
};

template <typename ScoreType>
std::ostream& operator<<(std::ostream& os, JumpAlignmentResult<ScoreType>& alignment);



/// \brief a method to align a contig to two references
///
/// the alignment can make a single jump from reference1 to reference2
///
template <typename ScoreType>
struct JumpAlignerBase : public AlignerBase<ScoreType>
{
    JumpAlignerBase(
        const AlignmentScores<ScoreType>& scores,
        const ScoreType jumpScore) :
        AlignerBase<ScoreType>(scores),
        _jumpScore(jumpScore)
    {}

    /// read-only access to the aligner's scores:
    const ScoreType&
    getJumpScore() const
    {
        return _jumpScore;
    }

protected:

    // backtrace logic shared with the intron jump aligner:
    template <typename SymIter, typename MatrixType>
    void
    backTraceAlignment(
        const SymIter queryBegin, const SymIter queryEnd,
        const SymIter ref1Begin, const SymIter ref1End,
        const SymIter ref2Begin, const SymIter ref2End,
        const size_t querySize, const size_t ref1Size, const size_t ref2Size,
        const MatrixType& ptrMatrix1, const MatrixType& ptrMatrix2,
        const BackTrace<ScoreType>& btraceInput,
        JumpAlignmentResult<ScoreType>& result) const;

    static
    uint8_t
    max4(
        ScoreType& max,
        const ScoreType v0,
        const ScoreType v1,
        const ScoreType v2,
        const ScoreType v3)
    {
        max=v0;
        uint8_t ptr=0;
        if (v1>v0)
        {
            max=v1;
            ptr=1;
        }
        if (v2>max)
        {
            max=v2;
            ptr=2;
        }
        if (v3>max)
        {
            max=v3;
            ptr=3;
        }
        return ptr;
    }

    const ScoreType _jumpScore;
};


#include "alignment/JumpAlignerBaseImpl.hh"


/// Convenience function to run alignment from multiple aligner classes:
template <typename JumpAligner, typename SymIter, typename ScoreType>
void
jumpAlign(
    const JumpAligner& jumpAligner,
    const SymIter queryBegin, const SymIter queryEnd,
    const SymIter ref1Begin, const SymIter ref1End,
    const SymIter ref2Begin, const SymIter ref2End,
    JumpAlignmentResult<ScoreType>& result)
{
    jumpAligner.align(queryBegin,queryEnd,ref1Begin,ref1End,ref2Begin,ref2End,result);
}


