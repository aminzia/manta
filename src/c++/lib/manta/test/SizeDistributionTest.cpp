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

#include "manta/SizeDistribution.cpp"



BOOST_AUTO_TEST_SUITE( test_SizeDistribution )

BOOST_AUTO_TEST_CASE( test_EmptySizeDistribution )
{
    SizeDistribution sd;

    BOOST_REQUIRE_EQUAL(sd.cdf(2),0);
    BOOST_REQUIRE_EQUAL(sd.quantile(0.2),0);
}


BOOST_AUTO_TEST_CASE( test_SizeDistribution1 )
{
    SizeDistribution sd;

    sd.addObservation(1);
    sd.addObservation(2);
    sd.addObservation(3);
    sd.addObservation(4);

    BOOST_REQUIRE_EQUAL(sd.cdf(0),0.);
    BOOST_REQUIRE_EQUAL(sd.cdf(1),0.25);
    BOOST_REQUIRE_EQUAL(sd.cdf(2),0.5);
    BOOST_REQUIRE_EQUAL(sd.cdf(3),0.75);
    BOOST_REQUIRE_EQUAL(sd.cdf(4),1);
    BOOST_REQUIRE_EQUAL(sd.quantile(0.0),1);
    BOOST_REQUIRE_EQUAL(sd.quantile(0.25),1);
    BOOST_REQUIRE_EQUAL(sd.quantile(0.5),2);
    BOOST_REQUIRE_EQUAL(sd.quantile(0.74),3);
    BOOST_REQUIRE_EQUAL(sd.quantile(0.75),3);
    BOOST_REQUIRE_EQUAL(sd.quantile(0.76),4);
    BOOST_REQUIRE_EQUAL(sd.quantile(1.0),4);
}

BOOST_AUTO_TEST_CASE( test_SizeDistributionFilter )
{
    SizeDistribution sd;

    sd.addObservation(1);
    sd.addObservation(2);
    sd.addObservation(3);
    sd.addObservation(4);

    sd.filterObservationsOverQuantile(0.5);

    BOOST_REQUIRE_EQUAL(sd.totalObservations(),2u);

    BOOST_REQUIRE_EQUAL(sd.cdf(0),0.);
    BOOST_REQUIRE_EQUAL(sd.cdf(1),0.5);
    BOOST_REQUIRE_EQUAL(sd.cdf(2),1);
    BOOST_REQUIRE_EQUAL(sd.quantile(0.0),1);
    BOOST_REQUIRE_EQUAL(sd.quantile(0.25),1);
    BOOST_REQUIRE_EQUAL(sd.quantile(0.5),1);
    BOOST_REQUIRE_EQUAL(sd.quantile(0.75),2);
    BOOST_REQUIRE_EQUAL(sd.quantile(1.0),2);
}

BOOST_AUTO_TEST_SUITE_END()

