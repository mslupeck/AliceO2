// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "FV0Simulation/Digitizer.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include <CommonDataFormat/InteractionRecord.h>

#include "TMath.h"
#include "TRandom.h"
#include <TH1F.h>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <optional>
#include "MathUtils/CachingTF1.h"
#include <boost/format.hpp>

using namespace o2::fv0;
using namespace o2::math_utils;
using boost::format;

ClassImp(Digitizer);

void Digitizer::process(const std::vector<o2::fv0::Hit>* hits, o2::fv0::Digit* digit, std::vector<std::vector<double>>& channel_times)

{

  auto sorted_hits{*hits};
  std::sort(sorted_hits.begin(), sorted_hits.end(), [](o2::fv0::Hit const& a, o2::fv0::Hit const& b) {
    return a.GetTrackID() < b.GetTrackID();
  });

  digit->setTime(mEventTime);
  digit->setInteractionRecord(mIntRecord);

  //Calculating signal time, amplitude in mean_time +- time_gate --------------
  std::vector<o2::fv0::ChannelData>& channel_data = digit->getChDgData();
  if (channel_data.size() == 0) {
    channel_data.reserve(parameters.mPmts);
    for (int i = 0; i < parameters.mPmts; ++i)
      channel_data.emplace_back(o2::fv0::ChannelData{i, 0, 0});
  }
  Int_t parent = -10;
  Float_t pmtTimeIntegral = mPmtResponse->Integral(-parameters.mPmtTransitTime, 2. * parameters.mPmtTransitTime);
  //LOG(INFO)<<"pmtTimeIntegral \t" <<pmtTimeIntegral;
  Float_t meansPhE = mSinglePhESpectrum->Mean(parameters.photoelMin, parameters.photoelMax);
  //LOG(INFO)<<"Single photo electron spectrum \t" <<meansPhE;
  for (Int_t i = 0; i < parameters.mPmts; i++)
    std::fill(mPmtChargeVsTime[i].begin(), mPmtChargeVsTime[i].end(), 0);

  assert(digit->getChDgData().size() == parameters.mPmts);

  for (auto const& hit : sorted_hits) {
    Int_t const pmt = hit.GetDetectorID();
    // LOG(INFO)<<"pmt ========="<< pmt;
    Double_t const hitValue = hit.GetHitValue() * 1e3; //convert to MeV
    if (hitValue < parameters.singleMipThreshold)
      continue;
    Double_t const nPhoton = hitValue * 10400;
    Int_t const nPhE = SimulateLightYield(pmt, nPhoton);
    Float_t const dt_scintillator = gRandom->Gaus(0, parameters.mIntrinsicTimeRes);
    Float_t const t = dt_scintillator + hit.GetTime() * 1e9;
    //LOG(INFO) << "dt_scintillator = "<<dt_scintillator<<" t "<<t;
    Float_t charge = TMath::Qe() * parameters.mPmtGain * mBinSize / pmtTimeIntegral;
    for (Int_t iPhE = 0; iPhE < nPhE; ++iPhE) {
      Float_t const tPhE = t + mSignalShape->GetRandom(0, mBinSize * Float_t(mNBins));
      //std::cout<<mSignalShape->GetRandom(0, mBinSize * Float_t(mNBins))<<std::endl;
      Float_t const gainVar = mSinglePhESpectrum->GetRandom(parameters.photoelMin, parameters.photoelMax) / meansPhE;
      Int_t const firstBin = TMath::Max((UInt_t)0, (UInt_t)((tPhE - parameters.mPmtTransitTime) / mBinSize));
      Int_t const lastBin = TMath::Min(mNBins - 1, (UInt_t)((tPhE + 2. * parameters.mPmtTransitTime) / mBinSize));
      //LOG(INFO) << "firstBin = "<<firstBin<<" lastbin "<<lastBin<<FairLogger::endl;
      for (Int_t iBin = firstBin; iBin <= lastBin; ++iBin) {
        Float_t const tempT = mBinSize * (0.5 + iBin) - tPhE;
        //mPmtChargeVsTime[pmt][iBin] += gainVar * charge * mPmtResponse->Eval(tempT);
        mPmtChargeVsTime[pmt][iBin] += gainVar * charge * PmtResponse(&tempT, &parameters.mPmtTransitTime);
        //LOG(INFO)<<"mPmtChargeVsTime"<<mBinSize*(0.5 + iBin) ;
      }
    } //photo electron loop

    //charge particles in MCLabel
    Int_t parentID = hit.GetTrackID();
    if (parentID != parent) {
      o2::fv0::MCLabel label(hit.GetTrackID(), mEventID, mSrcID, pmt);
      //LOG(INFO)<<"o2::fv0::MCLabel"<< hit.GetTrackID()<<"   "<<mEventID<<"  "<<mSrcID<<"  "<<pmt;

      int lblCurrent;
      if (mMCLabels) {
        lblCurrent = mMCLabels->getIndexedSize(); // this is the size of mHeaderArray;
        mMCLabels->addElement(lblCurrent, label);
      }
      parent = parentID;
    }
  } //hit loop
  for (Int_t ipmt = 0; ipmt < parameters.mPmts; ++ipmt) {
    channel_data[ipmt].mTime = SimulateTimeCfd(ipmt);
    for (Int_t iTimeBin = 0; iTimeBin < mNBins; ++iTimeBin)
      channel_data[ipmt].mChargeAdc += mPmtChargeVsTime[ipmt][iTimeBin] / parameters.mChargePerAdc;
    //LOG(INFO) << "Adc " << channel_data[ipmt].mChargeAdc << " Time " << channel_data[ipmt].mTime << FairLogger::endl;
  }
}

Float_t Digitizer::SimulateTimeCfd(Int_t channel)
{
  //std::fill(mPmtChargeVsTimeCfd.begin(), mPmtChargeVsTimeCfd.end(), 0);
  Float_t timeCfd = -1024;
  Int_t const binShift = TMath::Nint(parameters.mTimeShiftCfd / mBinSize);
  Float_t sigCurrent = 0;
  Float_t sigPrev = -mPmtChargeVsTime[channel][0];
  for (Int_t iTimeBin = 1; iTimeBin < mNBins; ++iTimeBin) {
    //if (mPmtChargeVsTime[channel][iTimeBin] != 0) std::cout << mPmtChargeVsTime[channel][iTimeBin] / parameters.mChargePerAdc << ", ";
    sigCurrent = (iTimeBin >= binShift
                    ? 5.0 * mPmtChargeVsTime[channel][iTimeBin - binShift] - mPmtChargeVsTime[channel][iTimeBin]
                    : -mPmtChargeVsTime[channel][iTimeBin]);
    if (sigPrev < 0 && sigCurrent >= 0) {
      timeCfd = mBinSize * Float_t(iTimeBin);
      //std::cout<<timeCfd<<std::endl;
      break;
    }
    sigPrev = sigCurrent;
  }
  return timeCfd;
}
//_____________________________________________________________________________
Int_t Digitizer::SimulateLightYield(Int_t pmt, Int_t nPhot)
{
  const Float_t epsilon = 0.0001;
  const Float_t p = parameters.mLightYield * parameters.mPhotoCathodeEfficiency;
  if ((p - 1.0f < epsilon) || nPhot == 0)
    return nPhot;
  const Int_t n = Int_t(nPhot < 100
                          ? gRandom->Binomial(nPhot, p)
                          : gRandom->Gaus(p * nPhot + 0.5, TMath::Sqrt(p * (1 - p) * nPhot)));
  return n;
}
//_____________________________________________________________________________
Float_t Digitizer::PmtResponse(const Float_t* x, const Float_t* t)
{
  // this function describes the PMT time response to a single photoelectron
  Float_t const y = x[0] + t[0]; //parameters.mPmtTransitTime, x -- time (PMT), t -- Pmt response time
  return y * y * TMath::Exp(-y * y / (t[0] * t[0]));
}
//_____________________________________________________________________________
Double_t Digitizer::SinglePhESpectrum(Double_t* x, Double_t*)
{
  // x -- number of photo-electrons emitted from the first dynode
  // this function describes the PMT amplitude response to a single photoelectron
  Double_t y = x[0];
  if (y < 0)
    return 0;
  return (TMath::Poisson(y, parameters.mPmtNbOfSecElec) +
          parameters.mPmtTransparency * TMath::Poisson(y, 1.0));
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void Digitizer::setTriggers(o2::fv0::Digit* digit)
{
  // --------------------------------------------------------------------------
} // trigger

void Digitizer::initParameters()
{
  mEventTime = 0;
}
//_______________________________________________________________________
void Digitizer::init()
{
  std::cout << " @@@ V0Digitizer::init " << std::endl;
  const float y = parameters.mPmtTransitTime;
  format pmresponse("(x+%1%)*(x+%2%)*TMath::Exp(-(x+%3%)*(x+%4%)/(%5% * %6%))");
  pmresponse % y % y % y % y % y % y;
  std::string pmtResponseFormula = pmresponse.str();

  mNBins = 2000;           //Will be computed using detector set-up from CDB
  mBinSize = 25.0 / 256.0; //Will be set-up from CDB

  for (Int_t i = 0; i < parameters.mPmts; i++) {
    mPmtChargeVsTime[i].resize(mNBins);
  }
  mPmtChargeVsTimeCfd.resize(mNBins);

  if (!mPmtResponse) {
    mPmtResponse = std::make_unique<o2::base::CachingTF1>("mPmtResponse", pmtResponseFormula.c_str(), -y, 2 * y);
    mPmtResponse->SetNpx(100);
  }

  if (!mSinglePhESpectrum)
    mSinglePhESpectrum = std::make_unique<o2::base::CachingTF1>("mSinglePhESpectrum", this, &Digitizer::SinglePhESpectrum, parameters.photoelMin, parameters.photoelMax, 0);
  if (!mSignalShape) {
    mSignalShape = std::make_unique<o2::base::CachingTF1>("mSignalShape", "crystalball", 0, 300);
    mSignalShape->SetParameters(1, parameters.mShapeSigma,
                                parameters.mShapeSigma, parameters.mShapeAlpha, parameters.mShapeN);
  } //signal shape
}
//_______________________________________________________________________
void Digitizer::finish()
{
  printParameters();
}

void Digitizer::printParameters()
{
}
