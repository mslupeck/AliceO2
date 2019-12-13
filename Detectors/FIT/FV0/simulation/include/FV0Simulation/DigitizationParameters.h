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
  const Int_t NCellsA = 40;
  const Int_t mPmts = 48;                       // number of PMTs
  const Float_t mIntrinsicTimeRes = 0.91;       // FDD 0.44
  const Float_t mPhotoCathodeEfficiency = 0.23; // quantum efficiency = nOfPhotoE_emitted_by_photocathode / nIncidentPhotons
  const Float_t mLightYield = 0.1;              // TODO: light collection efficiency to be tuned using collision data
  const Float_t mPmtGain = 4.21e6;              // for PMT R5924-70 at 0.5T
  const Float_t mChargePerAdc = 0.6e-12;        // charge
  const Float_t mPmtTransitTime = 9.5;          // PMT response time (corresponds to 1.9 ns rise time)
  const Float_t mPmtTransparency = 0.25;        // Transparency of the first dynode of the PM
  const Float_t mPmtNbOfSecElec = 9.0;          // Number of secondary electrons emitted from first dynode (per ph.e.)
  const Float_t mShapeAlpha = -0.445;
  const Float_t mShapeN = 2.65;
  const Float_t mShapeSigma = 3.25;
  //const Float_t mPedestal = 0;
  const Float_t mTimeShiftCfd = 1.42;     // TODO: copied from FDD, adjust after PM design for FV0 is fixed
  const Int_t photoelMin = 0;             // Integration lower limit
  const Int_t photoelMax = 30;            // Integration upper limit
  const Float_t singleMipThreshold = 3.0; // in [MeV] of deposited energy
};
} // namespace o2::fv0
#endif
