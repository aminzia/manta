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
/// \author Xiaoyu Chen
///

#include "SizeDistribution.hh"

#include "blt_util/log.hh"

#include "boost/foreach.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

//#define DEBUG_RPS


static
void
populateCdfQuantiles(
    SizeDistribution::map_type& sizeMap,
    const unsigned totalCount,
    std::vector<int>& quantiles)
{
    const unsigned quantileNum(quantiles.size());
    const float pFactor(1./static_cast<float>(totalCount));

    unsigned fillBase(0);
    unsigned cumulativeCount(0);
    BOOST_REVERSE_FOREACH(SizeDistribution::map_type::value_type& val, sizeMap)
    {
        cumulativeCount += (val.second.count);
        assert(cumulativeCount <= totalCount);

        // update the hash map with cumulative prob value
        val.second.cprob = (cumulativeCount * pFactor);

        const unsigned fillNext = static_cast<unsigned>(rint(val.second.cprob * quantileNum));
        for (; fillBase < fillNext; fillBase++)
        {
            quantiles[fillBase] = val.first;
        }
    }
}




void
SizeDistribution::
calcStats() const
{
#ifdef DEBUG_RPS
    log_os << "Calculating stats...\n"
           << "numOfSized=" << _sizeMap.size() << "\n";
#endif
    _isStatsComputed=true;
    if (_sizeMap.empty()) return;

    populateCdfQuantiles(_sizeMap, _totalCount, _quantiles);
}



int
SizeDistribution::
quantile(const float prob) const
{
    assert((prob >= 0.) && (prob <= 1.));

    static const int maxBin(_quantileNum - 1);
    if (! _isStatsComputed) calcStats();

    int bin(static_cast<int>(ceil(prob * _quantileNum) - 1));
    if (bin < 0) bin=0;
    if (bin > maxBin) bin=maxBin;
    return _quantiles[bin];
}



float
SizeDistribution::
cdf(const int size) const
{
    if (! _isStatsComputed) calcStats();

    const map_type::const_iterator sizeIter(_sizeMap.lower_bound(size));
    if (sizeIter == _sizeMap.end()) return 0;
    return sizeIter->second.cprob;
}



void
SizeDistribution::
filterObservationsOverQuantile(const float prob)
{
    const int maxSize(quantile(prob));
    const map_type::iterator sizeBegin(_sizeMap.begin());
    map_type::iterator sizeEnd(_sizeMap.lower_bound(maxSize));

    for (map_type::iterator sizeIter(sizeBegin); sizeIter != sizeEnd; ++sizeIter)
    {
        if (sizeIter->first <= maxSize)
        {
            sizeEnd = sizeIter;
            break;
        }
        _totalCount -= sizeIter->second.count;
    }
    _sizeMap.erase(sizeBegin,sizeEnd);

    _isStatsComputed=false;
}



std::ostream&
operator<<(std::ostream& os, const SizeDistribution& sd)
{
    os << sd.totalObservations() << '\n';
    return os;
}
