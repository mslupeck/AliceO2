// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_FV0_DIGITIZER_H
#define ALICEO2_FV0_DIGITIZER_H

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsFV0/Digit.h"
#include "DataFormatsFV0/MCLabel.h"
#include "FV0Simulation/Detector.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "FV0Simulation/DigitizationParameters.h"
#include <TH1F.h>
#include "MathUtils/CachingTF1.h"
#include "MathUtils/RandomRing.h"//added by a.khuntia
#include <TStopwatch.h>

namespace o2
{
    namespace fv0
    {
        class Digitizer
        {
        public:
            Digitizer(const DigitizationParameters& params, Int_t mode = 0)
            :mMode(mode),
            parameters(params),
            mTime(48),
            mPMResponse()
            {
                init();
            };

            ~Digitizer() = default;

            //void process(const std::vector<HitType>* hits, std::vector<Digit>* digits);
            void process(const std::vector<o2::fv0::Hit>* hits, o2::fv0::Digit* digit, std::vector<std::vector<double>>& channel_times);
            void computeAverage(o2::fv0::Digit& digit);

            void initParameters();
            void printParameters();
            void setEventTime(double value) { mEventTime = value; }
            void setEventID(Int_t id) { mEventID = id; }
            void setSrcID(Int_t id) { mSrcID = id; }
            void setInteractionRecord(uint16_t bc, uint32_t orbit)
            {
                mIntRecord.bc = bc;
                mIntRecord.orbit = orbit;
            }
            const o2::InteractionRecord& getInteractionRecord() const { return mIntRecord; }
            o2::InteractionRecord& getInteractionRecord(o2::InteractionRecord& src) { return mIntRecord; }
            void setInteractionRecord(const o2::InteractionRecord& src) { mIntRecord = src; }
            uint32_t getOrbit() const { return mIntRecord.orbit; }
            uint16_t getBC() const { return mIntRecord.bc; }

            void setTriggers(o2::fv0::Digit* digit);
            Int_t SimulateLightYield(Int_t pmt, Int_t nPhot);//fdd
            Float_t SimulateTimeCFD(Int_t channel);//fdd
            static Float_t PMResponse(const Float_t* x, const Float_t* par);//fdd
            Double_t SinglePhESpectrum(Double_t* x, Double_t* par);//fdd
            //void smearCFDtime(o2::fv0::Digit* digit, std::vector<std::vector<double>> const& channel_times);

            void init();
            void finish();

            void setMCLabels(o2::dataformats::MCTruthContainer<o2::fv0::MCLabel>* mclb) {mMCLabels = mclb; }
            //double get_time(const std::vector<double>& times);

        private:
            // digit info
            // parameters
            Int_t mMode;                      //triggered or continuos
            o2::InteractionRecord mIntRecord; // Interaction record (orbit, bc)
            Int_t mEventID;
            Int_t mSrcID;        // signal, background or QED
            Double_t mEventTime; // timestamp
            //int mNoisePeriod;    //low frequency noise period
            //int mBinshift;       // number of bin to shift positive part of CFD signal

            DigitizationParameters parameters;
            o2::dataformats::MCTruthContainer<o2::fv0::MCLabel>* mMCLabels = nullptr;

            std::vector<std::vector<Float_t>> mTime;                  // Charge time series aka analogue signal pulse from PM
            std::vector<Float_t> mTimeCFD;                            // Time series for CFD measurement
            UInt_t mNBins;                                            // Number of bins in pulse series
            Float_t mBinSize;                                         // Time width of the pulse bin - HPTDC resolution
            std::unique_ptr<o2::base::CachingTF1> mPMResponse;        // function which describes the PM time response
            std::unique_ptr<o2::base::CachingTF1> mSinglePhESpectrum; // function which describes the single ph.e. PM response
            std::unique_ptr<o2::base::CachingTF1> mSignalShape;

            TStopwatch mTimer;

            ClassDefNV(Digitizer, 2);//fdd
        };

    } // namespace fv0
} // namespace o2

#endif
