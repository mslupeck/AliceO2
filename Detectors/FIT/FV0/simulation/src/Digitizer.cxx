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

using namespace o2::fv0;
//using o2::fit::Geometry;

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
        channel_data.reserve(parameters.mMCPs);
        for (int i = 0; i < parameters.mMCPs; ++i)
            channel_data.emplace_back(o2::fv0::ChannelData{i,0, 0});
    }

    //for (Int_t i=0;i<parameters.NcellsA)
      //  std::fill(mTime[i].begin(),mTime.end(),0);

   // if (channel_times.size() == 0)
     //   channel_times.resize(parameters.mMCPs);
    Int_t parent = -10;
    Float_t integral = mPMResponse->Integral(-parameters.mPMTransitTime, 2. * parameters.mPMTransitTime);
    Float_t meansPhE = mSinglePhESpectrum->Mean(0, 30);
    //LOG(INFO)<<"Single photo electron spectrum \t" <<meansPhE;
    for (Int_t i = 0; i < parameters.mMCPs; i++)
        std::fill(mTime[i].begin(), mTime[i].end(), 0);

    assert(digit->getChDgData().size() == parameters.mMCPs);

    for (auto const& hit : sorted_hits) {
        Int_t const  pmt = hit.GetDetectorID();
       // LOG(INFO)<<"pmt ========="<< pmt;
        Double_t const hitValue = hit.GetHitValue()*1000;//convert to MeV
        if (hitValue < 6.0) continue;
        Double_t const nPhoton = hitValue*10400;
        Int_t const nPhE = SimulateLightYield(pmt,nPhoton);
        Float_t const dt_scintillator = gRandom->Gaus(0, parameters.mIntTimeRes);
        Float_t const t = dt_scintillator + hit.GetTime()*1e9;
        //LOG(INFO) << "dt_scintillator = "<<dt_scintillator<<" t "<<t;
        //LOG(INFO) << "Nphot = "<<nPhoton<<"  nPhE  "<<nPhE;
        //LOG(INFO) << "NphE = " << nPhE << FairLogger::endl;
        Float_t charge = TMath::Qe() * parameters.mPmGain * mBinSize / integral;
        for (Int_t iPhE = 0; iPhE < nPhE; ++iPhE) {
            Float_t const tPhE = t + mSignalShape->GetRandom(0, mBinSize * Float_t(mNBins));
            //std::cout<<mSignalShape->GetRandom(0, mBinSize * Float_t(mNBins))<<std::endl;
            Float_t const gainVar = mSinglePhESpectrum->GetRandom(0, 30) / meansPhE;
            Int_t const firstBin = TMath::Max((UInt_t)0, (UInt_t)((tPhE - parameters.mPMTransitTime) / mBinSize));
            Int_t const lastBin = TMath::Min(mNBins - 1, (UInt_t)((tPhE + 2. * parameters.mPMTransitTime) / mBinSize));
            //LOG(INFO) << "firstBin = "<<firstBin<<" lastbin "<<lastBin<<FairLogger::endl;
            for (Int_t iBin = firstBin; iBin <= lastBin; ++iBin) {
                Float_t const tempT = mBinSize * (0.5 + iBin) - tPhE;
                mTime[pmt][iBin] += gainVar * charge * mPMResponse->Eval(tempT);
                //LOG(INFO)<<"mTime"<<mBinSize*(0.5 + iBin) ;
            }
        }//photo electron loop

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
    }//hit loop
    for (Int_t ipmt = 0; ipmt < parameters.mMCPs; ++ipmt) {
        channel_data[ipmt].mTime = SimulateTimeCFD(ipmt);
        for (Int_t iBin = 0; iBin < mNBins; ++iBin)
            channel_data[ipmt].mChargeADC += mTime[ipmt][iBin] / parameters.mChargePerADC;
        //LOG(INFO) << "ADC " << channel_data[ipmt].mChargeADC << " Time " << channel_data[ipmt].mTime << FairLogger::endl;
    }
}


Float_t Digitizer::SimulateTimeCFD(Int_t channel)
{
    //std::fill(mTimeCFD.begin(), mTimeCFD.end(), 0);
    Float_t timeCFD = -1024;
    Int_t const  binShift = TMath::Nint(parameters.mTimeShiftCFD / mBinSize);
    Float_t sigCurrent = 0;
    Float_t sigPrev = -mTime[channel][0];
    for (Int_t iBin = 1; iBin < mNBins; ++iBin) {
        //if (mTime[channel][iBin] != 0) std::cout << mTime[channel][iBin] / parameters.mChargePerADC << ", ";
        sigCurrent = (iBin >= binShift
                ? 5.0 * mTime[channel][iBin - binShift] - mTime[channel][iBin]
                : -mTime[channel][iBin]);
        if (sigPrev < 0 && sigCurrent>=0) {
            timeCFD = mBinSize * Float_t(iBin);
            //std::cout<<timeCFD<<std::endl;
            break;
        }
        sigPrev = sigCurrent;
    }
    return timeCFD;
}
//_____________________________________________________________________________
Int_t Digitizer::SimulateLightYield(Int_t pmt, Int_t nPhot)
{
    const Float_t p = parameters.mLightYield * parameters.mPhotoCathodeEfficiency;
    if (p == 1.0f || nPhot == 0)
        return nPhot;
    const Int_t n = Int_t(nPhot < 100
            ? gRandom->Binomial(nPhot, p)
            : gRandom->Gaus(p * nPhot + 0.5, TMath::Sqrt(p * (1 - p) * nPhot)));
    return n;
}
//_____________________________________________________________________________
Double_t Digitizer::PMResponse(Double_t* x, Double_t*)
{
    // this function describes the PM time response to a single photoelectron
    Double_t y = x[0] + parameters.mPMTransitTime;
    return y * y * TMath::Exp(-y * y / (parameters.mPMTransitTime * parameters.mPMTransitTime));
}
//_____________________________________________________________________________
Double_t Digitizer::SinglePhESpectrum(Double_t* x, Double_t*)
{
    // this function describes the PM amplitude response to a single photoelectron
    Double_t y = x[0];
    if (y < 0)
        return 0;
    return (TMath::Poisson(y, parameters.mPMNbOfSecElec) +
                              parameters.mPMTransparency * TMath::Poisson(y, 1.0));
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void Digitizer::setTriggers(o2::fv0::Digit* digit)
{
    // --------------------------------------------------------------------------
}// trigger



void Digitizer::initParameters()
{
    mEventTime = 0;
}
//_______________________________________________________________________
void Digitizer::init()
{
    std::cout << " @@@ V0Digitizer::init " << std::endl;
    mNBins = 2000;           //Will be computed using detector set-up from CDB
    mBinSize = 25.0 / 256.0; //Will be set-up from CDB
    for (Int_t i = 0; i < parameters.mMCPs; i++)
        mTime[i].resize(mNBins);
        mTimeCFD.resize(mNBins);
    if (!mPMResponse)
        mPMResponse = std::make_unique<o2::base::CachingTF1>
                ("mPMResponse", this, &Digitizer::PMResponse,
                         -parameters.mPMTransitTime, 2. * parameters.mPMTransitTime, 0);
    if (!mSinglePhESpectrum)
    mSinglePhESpectrum = std::make_unique<o2::base::CachingTF1>
                         ("mSinglePhESpectrum", this, &Digitizer::SinglePhESpectrum, 0, 30, 0);
    if (!mSignalShape) {
      mSignalShape = std::make_unique<o2::base::CachingTF1>
                     ("mSignalShape", "crystalball", 0, 300);
      mSignalShape->SetParameters(1, parameters.mShapeSigma,
                                     parameters.mShapeSigma, parameters.mShapeAlpha, parameters.mShapeN);
    }//signal shape

}
//_______________________________________________________________________
void Digitizer::finish()
{
    printParameters();
}

void Digitizer::printParameters()
{
}
