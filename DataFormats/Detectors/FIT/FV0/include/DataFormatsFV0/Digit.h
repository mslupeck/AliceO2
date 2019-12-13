// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_FV0_DIGIT_H
#define ALICEO2_FV0_DIGIT_H

#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/TimeStamp.h"
#include <iosfwd>
#include "Rtypes.h"

namespace o2
{
namespace fv0
{

struct ChannelData {
  //Int_t ChId;       // channel Id
  //Double_t CFDTime; // time in ns, 0 at lhc clk center
  //Double_t QTCAmpl; // Amplitude in mips
  //int numberOfParticles;
  Int_t mPmtNumber;   // PhotoMultiplier number (0 to 47)
  Float_t mTime;      // Time of Flight
  Short_t mChargeAdc; // ADC sample as present in raw data
  //Bit information from FEE
  //Bool_t mIntegrator;
  //Bool_t mDoubleEvent;
  //Bool_t mEvent1TimeLost;
  //Bool_t mEvent2TimeLost;
  //Bool_t mAdcInGate;
  //Bool_t mTimeTooLate;
  //Bool_t mAmpTooHigh;
  //Bool_t mEventInTrigger;
  //Bool_t mTimeLost;

  ClassDefNV(ChannelData, 1);
};

/// \class Digit
/// \brief FIT digit implementation
using DigitBase = o2::dataformats::TimeStamp<double>;
class Digit : public DigitBase
{
 public:
  Digit() = default;

  Digit(std::vector<ChannelData> ChDgDataArr, Double_t time, uint16_t bc, uint32_t orbit, std::vector<Bool_t> const& triggers)
  {
    setChDgData(std::move(ChDgDataArr));
    setTime(time);
    setInteractionRecord(bc, orbit);
    setTriggers(std::move(triggers));
  }

  ~Digit() = default;

  Double_t getTime() const { return mTime; }
  void setTime(Double_t time) { mTime = time; }

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

  const std::vector<Bool_t>& GetTriggers() const { return mTriggers; }
  std::vector<Bool_t>& GetTriggers() { return mTriggers; }
  void setTriggers(const std::vector<Bool_t>& triggers) { mTriggers = triggers; }
  void setTriggers(std::vector<Bool_t>&& triggers) { mTriggers = std::move(triggers); }

  const std::vector<ChannelData>& getChDgData() const { return mChDgDataArr; }
  std::vector<ChannelData>& getChDgData() { return mChDgDataArr; }
  void setChDgData(const std::vector<ChannelData>& ChDgDataArr) { mChDgDataArr = ChDgDataArr; }
  void setChDgData(std::vector<ChannelData>&& ChDgDataArr) { mChDgDataArr = std::move(ChDgDataArr); }

  void printStream(std::ostream& stream) const;

  void cleardigits()
  {
    mTriggers.clear();
    mChDgDataArr.clear();
  }

 private:
  Double_t mTime;                   // time stamp
  o2::InteractionRecord mIntRecord; // Interaction record (orbit, bc)

  //online triggers processed on TCM

  std::vector<Bool_t> mTriggers; // TODO: unsued until triggers of FV0 are defined
  std::vector<ChannelData> mChDgDataArr;

  ClassDefNV(Digit, 1);
};

std::ostream& operator<<(std::ostream& stream, const Digit& digi);
} // namespace fv0
} // namespace o2
#endif
