// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_FV0_DIGITIZATION_PARAMETERS
#define ALICEO2_FV0_DIGITIZATION_PARAMETERS

namespace o2::fv0
{
    struct DigitizationParameters {
        Int_t NCellsA = 40;
        Int_t  mMCPs=48;            //number of MCPs
        bool mIsFV0= true;          //amplitude FT0(true) or FV0 (false)
        const Float_t mIntTimeRes = 0.91;//FDD 0.44
        const Float_t mPhotoCathodeEfficiency = 0.23;//FDD 0.18
        const Float_t mLightYield = 0.1;
        const Float_t mPmGain = 4.21e6;        //for PMT R5924-70 at 0.5T
        const Float_t mChargePerADC = 0.6e-12;
        const Float_t mPMTransitTime = 9.5;    // PM response time (corresponds to 1.9 ns rise time)
        const Float_t mPMTransparency = 0.25;  // Transparency of the first dynode of the PM
        const Float_t mPMNbOfSecElec = 9.0;    // Number of secondary electrons emitted from first dynode (per ph.e.)
        const Float_t mShapeAlpha = -0.445;
        const Float_t mShapeN = 2.65;
        const Float_t mShapeSigma = 3.25;
        //const Float_t mPedestal = 0;
        const Float_t mTimeShiftCFD = 1.42;
    };
} // namespace o2::fv0
#endif
