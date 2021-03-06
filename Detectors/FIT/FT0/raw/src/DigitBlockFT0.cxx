// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "FT0Raw/DigitBlockFT0.h"
using namespace o2::ft0;

int DigitBlockFT0::sEventID = 0;
o2::ft0::LookUpTable DigitBlockFT0::sLookupTable = o2::ft0::LookUpTable::readTable();
int DigitBlockFT0ext::sEventID = 0;
o2::ft0::LookUpTable DigitBlockFT0ext::sLookupTable = o2::ft0::LookUpTable::readTable();
