// Copyright (C) 2010 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Common/HourMinSec.hpp"

using namespace std;
using namespace boost;
using namespace CF;
using namespace CF::Common;

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE( HourMinSec_TestSuite )

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( constructors )
{
  // empty constructor
  HourMinSec hms1;
  BOOST_CHECK_EQUAL(hms1.str(), "0 sec");

  // non-empty constructor
  HourMinSec hms2(4567);
  BOOST_CHECK_EQUAL(hms2.str(), "1 h 16 min 7 sec");
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( set )
{
  HourMinSec hms;
  
  hms.set(0);
  BOOST_CHECK_EQUAL(hms.str(), "0 sec");

  hms.set(3600);
  BOOST_CHECK_EQUAL(hms.str(), "1 h 0 sec");
  
  hms.set(56584);
  BOOST_CHECK_EQUAL(hms.str(), "15 h 43 min 4 sec");

  hms.set(333133);
  BOOST_CHECK_EQUAL(hms.str(), "92 h 32 min 13 sec");
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()

////////////////////////////////////////////////////////////////////////////////

