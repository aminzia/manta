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

/// \file

/// \author Chris Saunders
///
#pragma once

#include "blt_util/bam_record.hh"

#include "boost/utility.hpp"

#include <string>



/// convenient bam record iterator for whole genome or chromosome segments
///
//
// Example use:
// while (stream.next()) {
//     const bam_record& read(*(stream.get_record_ptr()));
//     if(read.is_unmapped) foo++;
// }
//
struct bam_streamer : public boost::noncopyable
{

    explicit
    bam_streamer(const char* filename,
                 const char* region = NULL);

    ~bam_streamer();

    /// \brief set new or first region for file:
    void
    set_new_region(const char* region);

    /// \brief set new or first region for file:
    ///
    /// \param beg zero-indexed start pos
    /// \param end zero-indexed end pos
    void
    set_new_region(int reg, int beg, int end);

    bool next();

    const bam_record* get_record_ptr() const
    {
        if (_is_record_set) return &_brec;
        else               return NULL;
    }

    const char* name() const
    {
        return _stream_name.c_str();
    }

    unsigned record_no() const
    {
        return _record_no;
    }

    void report_state(std::ostream& os) const;

    const char*
    target_id_to_name(const int32_t tid) const;

    int32_t
    target_name_to_id(const char* seq_name) const;

    const bam_header_t*
    get_header() const
    {
        return _bfp->header;
    }

private:
    void _load_index();

    bool _is_record_set;
    samfile_t* _bfp;
    bam_index_t* _bidx;
    bam_iter_t _biter;
    bam_record _brec;

    // track for debug only:
    unsigned _record_no;
    std::string _stream_name;
    bool _is_region;
    std::string _region;
};

