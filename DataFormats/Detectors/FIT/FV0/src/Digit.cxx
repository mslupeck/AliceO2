//
// Created by Arvind Khuntia on 26/08/19.
//

// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "DataFormatsFV0/Digit.h"
#include <iostream>
#include <iosfwd>
#include "Rtypes.h"

using namespace o2::fv0;

void Digit::printStream(std::ostream& stream) const
{
  stream << "FV0 Digit: event time " << mTime << " BC " << mIntRecord.bc << " orbit " << mIntRecord.orbit << std::endl;
}

std::ostream& operator<<(std::ostream& stream, const Digit& digi)
{
  digi.printStream(stream);
  return stream;
}
