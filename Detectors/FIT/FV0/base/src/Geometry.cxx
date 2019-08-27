// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   Geometry.cxx
/// \brief  Implementation of FV0+ geometry.
///
/// \author Maciej Slupecki, University of Jyvaskyla, Finland
/// \author Andreas Molander, University of Helsinki, Finland

// #include <iomanip>
#include "FV0Base/Geometry.h"

#include <FairLogger.h>
#include <TGeoBBox.h>
#include <TGeoCompositeShape.h>
#include <TGeoCone.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoMedium.h>
#include <TGeoTube.h>
#include <TGeoVolume.h>

// #include <cmath>
// #include <sstream>

ClassImp(o2::fv0::Geometry);

using namespace o2::fv0;

Geometry::Geometry(EGeoType initType)
{
  mGeometryType = initType;
  initializeGeometry();
}

Geometry::Geometry(const Geometry& geometry)
{
  this->mGeometryType = geometry.mGeometryType;
  this->mEnabledComponents = geometry.mEnabledComponents;
}

const int Geometry::getCurrentCellId(TVirtualMC* fMC)
{
  int ring = -1;
  int sector = -1;

  fMC->CurrentVolOffID(0, ring);
  fMC->CurrentVolOffID(1, sector);

  // // TODO: Remove (this is for debuging purposes during development only)
  // LOG(INFO) << "FV0 Geometry::getCurrentCellId():";
  // LOG(INFO) << "Ring id  :   " << ring;
  // LOG(INFO) << "Ring name:   " << fMC->CurrentVolOffName(0);
  // LOG(INFO) << "Sector id:   " << sector;
  // LOG(INFO) << "Sector name: " << fMC->CurrentVolOffName(1);
  // LOG(INFO) << "Cell id :    " << sector + 8 * (ring - 1);
  // LOG(INFO) << "";

  return sector + 8 * (ring - 1);
}

bool Geometry::enableComponent(EGeoComponent component, bool enable)
{
  if (mEnabledComponents.find(component) == mEnabledComponents.end()) {
    LOG(WARNING) << "FV0 Geometry::enableComponent(): Component not initialized and cannot be enabled/disabled!";
    return false;
  }

  return mEnabledComponents[component] = enable;
}

void Geometry::buildGeometry()
{
  TGeoVolume* vALIC = gGeoManager->GetVolume("cave");
  if (!vALIC) {
    LOG(FATAL) << "Could not find the top volume";
  }

  // Top volume of FIT V0 detector
  TGeoVolumeAssembly* vFV0 = new TGeoVolumeAssembly("FV0");
  LOG(INFO) << "FV0 Geometry::buildGeometry():: FV0 volume name = " << vFV0->GetName();

  assembleSensVols(vFV0);
  assembleNonSensVols(vFV0);

  TGeoTranslation* globalShift = new TGeoTranslation(sXGlobal, sYGlobal, sZGlobal);

  // TODO: What shoud the copy_no be?
  vALIC->AddNode(vFV0, 1, globalShift);
}

void Geometry::initializeGeometry()
{
  initializeMaps();
  initializeVectors();
  initializeSensVols();
  initializeNonSensVols();
}

void Geometry::initializeMaps()
{
  bool isFull = mGeometryType == eFull;
  bool isRough = isFull || mGeometryType == eRough;
  bool hasScint = isRough || mGeometryType == eOnlySensitive;
  mEnabledComponents.insert(std::pair<EGeoComponent, bool>(eScintillator, hasScint));
  mEnabledComponents.insert(std::pair<EGeoComponent, bool>(ePlastics, isRough));
  mEnabledComponents.insert(std::pair<EGeoComponent, bool>(eFibers, isFull));
  mEnabledComponents.insert(std::pair<EGeoComponent, bool>(eScrews, isFull));
  mEnabledComponents.insert(std::pair<EGeoComponent, bool>(eRods, isFull));
  mEnabledComponents.insert(std::pair<EGeoComponent, bool>(eAluminiumContainer, isRough));
}

void Geometry::initializeVectors()
{
  initializeCellRingRadii();
  initializeSectorTransformations();
  initializeFiberVolumeRadii();
  initializeFiberMedium();
  initializeScrewAndRodRadii();
  initializeScrewTypeMedium();
  initializeRodTypeMedium();
  initializeScrewAndRodPositionsAndDimensions();
}

void Geometry::initializeCellRingRadii()
{
  LOG(INFO) << "FV0 Geometry::initializeCellRadii(): Initializing FV0 cell ring radii";

  // Index of mRAvgRing is NOT linked directly to any ring number
  mRAvgRing.assign(sCellRingRadii, sCellRingRadii + sNCellRings + 1);

  // Set real scintillator radii (reduced by paint thickness and separation gap)
  for (int i = 0; i < mRAvgRing.size() - 1; i++) {
    mRMinScint.push_back(mRAvgRing[i] + sDrSeparationScint);
    mRMaxScint.push_back(mRAvgRing[i + 1] - sDrSeparationScint);
    LOG(INFO) << "FV0 Geometry::initializeCellRadii(): Ring " << i + 1 << " min: " << mRMinScint[i] << ", max: " << mRMaxScint[i];
  }
  // Now indices of mRMinScint and mRMaxScint correspond to the same ring
}

void Geometry::initializeSectorTransformations()
{
  LOG(INFO) << "FV0 Geometry::initializeSectorTransformations(): Initializing FV0 sector transformations";

  for (int iSector = 1; iSector <= sNCellSectors; iSector++) {
    // iSector = 1 corresponds to the first sector counter-clockwise from the y-axis + global azimuthal rotation
    // iSector = 2 corresponds to the next sector in counter-clockwise direction

    // if the sector is rotated, this is the angle (plus an extra 45 for "b"-type cells)
    float phi = iSector * 45;

    std::string rotName = "FV0_rotationSector" + std::to_string(iSector);
    std::string reflName = "FV0_reflectSector" + std::to_string(iSector);
    std::string totalTransName = "FV0_totalTransformSector" + std::to_string(iSector);

    TGeoRotation* rot = new TGeoRotation(rotName.c_str());            // rotation of the sector
    TGeoRotation* refl = new TGeoRotation(reflName.c_str());          // reflection of the sector
    TGeoHMatrix* totTrans = new TGeoHMatrix(totalTransName.c_str());  // the combined transformation

    // The reference to "a" and "b" can be understood with the CAD drawings of the detector.
    switch (iSector)
    {
    case 1: // "a"-mirror
      refl->ReflectX(true);
      break;
    case 2: // "b"-mirror
      refl->ReflectX(true);
      break;
    case 3: // "b"
      rot->SetAngles(phi + 45, 0, 0);
      break;
    case 4: // "a"
      rot->SetAngles(phi, 0, 0);
      break;
    case 5: // "a"-mirror
      refl->ReflectY(true);
      break;
    case 6: // "b"-mirror
      refl->ReflectY(true);
      break;
    case 7: // "b"
      rot->SetAngles(phi + 45, 0, 0);
      break;
    case 8: // "a"
      rot->SetAngles(phi, 0, 0);
      break;
    default:
      break;
    }

    // Rotate the sector with the global rotation angle
    rot->SetAngles(rot->GetPhiRotation() + sGlobalPhiRotation, 0, 0);

    // Combine the rotation and reflection
    *totTrans = *rot * *refl;
    totTrans->RegisterYourself();
    mSectorTrans.push_back(totTrans);
  }
}

void Geometry::initializeFiberVolumeRadii()
{
  mRMinFiber.push_back(0);
  mRMinFiber.push_back(sDrMinAluCone + sEpsilon);
  mRMinFiber.push_back(sDrMinAluFront + sEpsilon);

  mRMaxFiber.push_back(sDrMinAluCone);
  mRMaxFiber.push_back(sDrMinAluFront);
  mRMaxFiber.push_back(mRMaxScint.back());
}

void Geometry::initializeFiberMedium()
{
  // TODO: What happens here if the medium is not founf?
  mMediumFiber.push_back(gGeoManager->GetMedium("FV0_FiberInner$"));
  mMediumFiber.push_back(gGeoManager->GetMedium("FV0_FiberMiddle$"));
  mMediumFiber.push_back(gGeoManager->GetMedium("FV0_FiberOuter$"));
}

void Geometry::initializeScrewAndRodRadii()
{
  mRScrewAndRod.push_back(mRAvgRing[1]);
  mRScrewAndRod.push_back(mRAvgRing[2]);
  mRScrewAndRod.push_back(mRAvgRing[3]);
  mRScrewAndRod.push_back(mRAvgRing[4]);
  mRScrewAndRod.push_back((mRAvgRing[4] + mRAvgRing[5]) / 2);
  mRScrewAndRod.push_back(mRAvgRing[5]);
}

void Geometry::initializeScrewTypeMedium()
{
  TGeoMedium* steel = gGeoManager->GetMedium("FV0_Stainless_Steel$");
  for (int i = 0; i < sNScrewTypes; i++) {
    mMediumScrewTypes.push_back(steel);
  }
}

void Geometry::initializeRodTypeMedium()
{
  TGeoMedium* steel = gGeoManager->GetMedium("FV0_Stainless_Steel$");
  for (int i = 0; i < sNRodTypes; i++) {
    mMediumRodTypes.push_back(steel);
  }
}

void Geometry::addScrewProperties(int screwTypeID, int iRing, float phi) {
  float r = mRScrewAndRod[iRing];
  mScrewTypeIDs.push_back(screwTypeID);
  mScrewPos.push_back(std::vector<float> { cosf(phi * M_PI/180) * r,
                                           sinf(phi * M_PI/180) * r,
                                           sPosScint[2] - sDzScint / 2 + sZShiftScrew + sDzMaxScrewTypes[screwTypeID] / 2 });
  mDrMinScrews.push_back(sDrMinScrewTypes[screwTypeID]);
  mDrMaxScrews.push_back(sDrMaxScrewTypes[screwTypeID]);
  mDzMaxScrews.push_back(sDzMaxScrewTypes[screwTypeID]);
  mDzMinScrews.push_back(sDzMinScrewTypes[screwTypeID]);
}

void Geometry::addRodProperties(int rodTypeID, int iRing) {
  mRodTypeIDs.push_back(rodTypeID);
  mRodPos.push_back(std::vector<float>{ sDxMinRodTypes[rodTypeID] / 2,
                                        mRScrewAndRod[iRing],
                                        sPosScint[2] - sDzScint / 2 + sZShiftRod + sDzMaxRodTypes[rodTypeID] / 2 });
  mDxMinRods.push_back(sDxMinRodTypes[rodTypeID]);
  mDzMaxRods.push_back(sDxMaxRodTypes[rodTypeID]);
  mDyMinRods.push_back(sDyMinRodTypes[rodTypeID]);
  mDyMaxRods.push_back(sDyMaxRodTypes[rodTypeID]);
  mDzMaxRods.push_back(sDzMaxRodTypes[rodTypeID]);
  mDzMinRods.push_back(sDzMinRodTypes[rodTypeID]);

  mRodTypeIDs.push_back(rodTypeID);
  mRodPos.push_back(std::vector<float>{ sDxMinRodTypes[rodTypeID] / 2,
                                        -mRScrewAndRod[iRing],
                                        sPosScint[2] - sDzScint / 2 + sZShiftRod + sDzMaxRodTypes[rodTypeID] / 2 });
  mDxMinRods.push_back(sDxMinRodTypes[rodTypeID]);
  mDzMaxRods.push_back(sDxMaxRodTypes[rodTypeID]);
  mDyMinRods.push_back(sDyMinRodTypes[rodTypeID]);
  mDyMaxRods.push_back(sDyMaxRodTypes[rodTypeID]);
  mDzMaxRods.push_back(sDzMaxRodTypes[rodTypeID]);
  mDzMinRods.push_back(sDzMinRodTypes[rodTypeID]);
}

void Geometry::initializeScrewAndRodPositionsAndDimensions()
{
  LOG(INFO) << "FV0 Geometry::initializeScrewAndRodPositionsAndDimensions(): Initializing screw positions and dimensions";

  for (int iRing = 0; iRing < mRScrewAndRod.size(); iRing++) {
    switch (iRing)
    {
    case 0:
      addRodProperties(0, iRing);
      for (float phi = 45; phi >= -45; phi -= 45) {
        addScrewProperties(0, iRing, phi);
      }
      break;
    case 1:
      addRodProperties(0, iRing);
      for (float phi = 45; phi >= -45; phi -= 45) {
        addScrewProperties(1, iRing, phi);
      }
      break;
    case 2:
      addRodProperties(1, iRing);
      for (float phi = 67.5; phi >= -67.5; phi -= 22.5) {
        addScrewProperties(2, iRing, phi);
      }
      break;
    case 3:
      addRodProperties(2, iRing);
      for (float phi = 67.5; phi >= -67.5; phi -= 22.5) {
        addScrewProperties(3, iRing, phi);
      }
      break;
    case 4:
      addRodProperties(3, iRing);
      for (float phi = 45; phi >= -45; phi -= 45) {
        addScrewProperties(4, iRing, phi);
      }
      break;
    case 5:
      addRodProperties(3, iRing);
      for (float phi = 67.5; phi >= -67.5; phi -= 22.5) {
        addScrewProperties(5, iRing, phi);
      }
      break;
    default:
      break;
    }
  }
}

void Geometry::initializeSensVols()
{
  initializeScintCells();
}

void Geometry::initializeNonSensVols()
{
  initializeScrewHoles();
  initializeRodHoles();
  initializePlasticCells();
  initializeFibers();
  initializeScrews();
  initializeRods();
  initializeMetalContainer();
}

void Geometry::initializeScrewHoles()
{
  LOG(INFO) << "FV0 Geometry::initializeScrewHoles(): Initializing screw holes";

  std::string boolFormula = "";

  for (int i = 0; i < mScrewPos.size(); i++) {
    std::string holeShapeName = "FV0SCREWHOLE" + std::to_string(i);
    std::string holeTransName = "FV0SCREWHOLETRANS" + std::to_string(i);

    createScrewShape(holeShapeName, mScrewTypeIDs[i], sEpsilon, sEpsilon);
    createAndRegisterTrans(holeTransName, mScrewPos[i][0] + sXShiftScrews, mScrewPos[i][1], mScrewPos[i][2]);

    if (i != 0) {
      boolFormula += "+";
    }
    boolFormula += holeShapeName + ":" + holeTransName;
  }

  new TGeoCompositeShape(sScrewHolesCSName.c_str(), boolFormula.c_str());

  LOG(INFO) << "FV0 Geometry::initializeScrewHoles(): Screw holes initialized";
}

void Geometry::initializeRodHoles()
{
  LOG(INFO) << "FV0 Geometry::initializeRodHoles(): Initializing rod holes";

  std::string boolFormula = "";

  for (int i = 0; i < mRodPos.size(); i++) {
    std::string holeShapeName = "FV0" + sRodName + "HOLE" + std::to_string(i);
    std::string holeTransName = "FV0" + sRodName + "HOLETRANS" + std::to_string(i);

    createRodShape(holeShapeName, mRodTypeIDs[i], sEpsilon, sEpsilon);
    createAndRegisterTrans(holeTransName, mRodPos[i][0] + sXShiftScrews, mRodPos[i][1], mRodPos[i][2]);

    if (i != 0) {
      boolFormula += "+";
    }
    boolFormula += holeShapeName + ":" + holeTransName;
  }

  new TGeoCompositeShape(sRodHolesCSName.c_str(), boolFormula.c_str());
}

void Geometry::initializeCells(std::string cellType, float zThickness, TGeoMedium* medium, bool isSensitive) {
  // Creating the two types of cells, "a" and "b", for each ring.
  // All sectors can be assembled with these cells.
  //
  // The reference to "a" and "b" can be understood with the CAD drawings of the detector.

  LOG(INFO) << "FV0 Geometry::initializeCells(): Initializing " << cellType << " cells";

  float dxHoleCut = sDxHoleExtensionScint;           // width of extension of hole 1, 2 and 7 in the "a" cell
  float xHole = sDrSeparationScint + dxHoleCut;      // x-placement of holes 1, 2 and 7 in the "a" cell

  // Sector separation gap shape
  std::string secSepShapeName = "FV0_" + cellType + "SectorSeparation";
  new TGeoBBox(secSepShapeName.c_str(), mRMaxScint.back() + sEpsilon, sDrSeparationScint, zThickness / 2);

  // Sector separation gap rotations
  std::string secSepRot45Name = "FV0_" + cellType + "SecSepRot45";
  std::string secSepRot90Name = "FV0_" + cellType + "SecSepRot90";

  createAndRegisterRot(secSepRot45Name, 45, 0, 0);
  createAndRegisterRot(secSepRot90Name, 90, 0, 0);

  // Hole shapes
  std::string holeSmallName = "FV0_" + cellType + "HoleSmall";
  std::string holeLargeName = "FV0_" + cellType + "HoleLarge";
  std::string holeSmallCutName = "FV0_" + cellType + "HoleSmallCut";
  std::string holeLargeCutName = "FV0_" + cellType + "HoleLargeCut";

  new TGeoTube(holeSmallName.c_str(), 0, sDrHoleSmallScint, zThickness / 2);
  new TGeoTube(holeLargeName.c_str(), 0, sDrHoleLargeScint, zThickness / 2);
  new TGeoBBox(holeSmallCutName.c_str(), dxHoleCut, sDrHoleSmallScint, zThickness / 2);
  new TGeoBBox(holeLargeCutName.c_str(), dxHoleCut, sDrHoleLargeScint, zThickness / 2);

  for (int ir = 0; ir < sNCellRings; ir++) {

    // Radii without separation
    float rMin = mRAvgRing[ir];
    float rMax = mRAvgRing[ir + 1];
    float rMid = rMin + (rMax - rMin) / 2;

    // "a"-type cell
    //
    // Initial placement:
    //
    // y
    // ^
    // |  * * *
    // |  * a * *
    // |  * * *
    // |  * *
    // |
    // 0--------------> x

    std::string aCellName = "FV0" + cellType + sCellName + "a" + std::to_string(ir + 1);

    LOG(INFO) << "FV0 Geometry::initializeCells(): Initializing cell " << aCellName;

    // Base shape
    std::string aCellShapeName = aCellName + "Shape";

    // The cells in the innermost ring have a slightly shifted inner radius origin.
    if (ir == 0) {
      // The innermost "a"-type cell
      std::string a1CellShapeFullName = aCellShapeName + "Full";
      std::string a1CellShapeHoleCutName = aCellShapeName + "HoleCut";
      std::string a1CellShapeHoleCutTransName = a1CellShapeHoleCutName + "Trans";

      new TGeoTubeSeg(a1CellShapeFullName.c_str(), 0, mRMaxScint[ir], zThickness / 2 - sEpsilon, 45, 90);
      new TGeoTube(a1CellShapeHoleCutName.c_str(), 0, mRMinScint[ir], zThickness);

      createAndRegisterTrans(a1CellShapeHoleCutTransName, sXShiftInnerRadiusScint, 0, 0);

      std::string a1BoolFormula = a1CellShapeFullName + "-" + a1CellShapeHoleCutName + ":" + a1CellShapeHoleCutTransName;
      new TGeoCompositeShape(aCellShapeName.c_str(), a1BoolFormula.c_str());
    } else {
      // The rest of the "a"-type cells
      new TGeoTubeSeg(aCellShapeName.c_str(), mRMinScint[ir], mRMaxScint[ir], zThickness / 2, 45, 90);
    }

    // Translations for screw holes
    //
    // 1 = outer left
    // 2 = inner left
    // 3 = outer right
    // 4 = inner right
    // 5 = outer middle
    // 6 = inner middle
    // 7 = half-lenght left
    // 8 = half-length right
    //
    // holes 1, 2 and 7 are sligtly shifted along the rim of the cell

    std::string aHole1TransName = aCellName + "Hole1Trans";
    std::string aHole2TransName = aCellName + "Hole2Trans";
    std::string aHole3TransName = aCellName + "Hole3Trans";
    std::string aHole4TransName = aCellName + "Hole4Trans";
    std::string aHole5TransName = aCellName + "Hole5Trans";
    std::string aHole6TransName = aCellName + "Hole6Trans";
    std::string aHole7TransName = aCellName + "Hole7Trans";
    std::string aHole8TransName = aCellName + "Hole8Trans";
    std::string aHole1CutTransName = aCellName + "Hole1CutTrans";
    std::string aHole2CutTransName = aCellName + "Hole2CutTrans";
    std::string aHole7CutTransName = aCellName + "Hole7CutTrans";

    createAndRegisterTrans(aHole1TransName, xHole, cos(asin(xHole/rMax)) * rMax, 0);
    createAndRegisterTrans(aHole2TransName, xHole, cos(asin(xHole/rMin)) * rMin, 0);
    createAndRegisterTrans(aHole3TransName, sin(45 * M_PI/180) * rMax, cos(45 * M_PI/180) * rMax, 0);
    createAndRegisterTrans(aHole4TransName, sin(45 * M_PI/180) * rMin, cos(45 * M_PI/180) * rMin, 0);
    createAndRegisterTrans(aHole5TransName, sin(22.5 * M_PI/180) * rMax, cos(22.5 * M_PI/180) * rMax, 0);
    createAndRegisterTrans(aHole6TransName, sin(22.5 * M_PI/180) * rMin, cos(22.5 * M_PI/180) * rMin, 0);
    createAndRegisterTrans(aHole7TransName, xHole, cos(asin(xHole/rMid)) * rMid, 0);
    createAndRegisterTrans(aHole8TransName, sin(45 * M_PI/180) * rMid, cos(45 * M_PI/180) * rMid, 0);
    createAndRegisterTrans(aHole1CutTransName, 0, cos(asin(xHole/rMax)) * rMax, 0);
    createAndRegisterTrans(aHole2CutTransName, 0, cos(asin(xHole/rMin)) * rMin, 0);
    createAndRegisterTrans(aHole7CutTransName, 0, cos(asin(xHole/rMid)) * rMid, 0);

    // Composite shape
    std::string aBoolFormula = aCellShapeName;

    // sector separation
    aBoolFormula += "-" + secSepShapeName + ":" + secSepRot45Name;
    aBoolFormula += "-" + secSepShapeName + ":" + secSepRot90Name;

    // outer holes
    aBoolFormula += "-" + ((ir < 2) ? holeSmallName : holeLargeName) + ":" + aHole1TransName;
    aBoolFormula += "-" + ((ir < 2) ? holeSmallCutName : holeLargeCutName) + ":" + aHole1CutTransName;
    aBoolFormula += "-" + ((ir < 2) ? holeSmallName : holeLargeName) + ":" + aHole3TransName;

    // inner holes
    if (ir > 0)  {
      std::string screwHoleName = (ir < 3) ? holeSmallName : holeLargeName;
      std::string screwHoleCutName = (ir < 3) ? holeSmallCutName : holeLargeCutName;

      aBoolFormula += "-" + screwHoleName + ":" + aHole2TransName;
      aBoolFormula += "-" + screwHoleCutName + ":" + aHole2CutTransName;
      aBoolFormula += "-" + screwHoleName + ":" + aHole4TransName;
    }

    // outer middle hole
    if (ir > 1) {
      aBoolFormula += "-" + holeLargeName + ":" + aHole5TransName;
    }

    // inner middle hole
    if (ir > 2) {
      aBoolFormula += "-" + holeLargeName + ":" + aHole6TransName;
    }

    // half-length holes
    if (ir == 4) {
      aBoolFormula += "-" + holeLargeName + ":" + aHole7TransName;
      aBoolFormula += "-" + holeLargeCutName + ":" + aHole7CutTransName;
      aBoolFormula += "-" + holeLargeName + ":" + aHole8TransName;
    }

    std::string aCellCSName = aCellName + "CS";
    TGeoCompositeShape* aCellCs = new TGeoCompositeShape(aCellCSName.c_str(), aBoolFormula.c_str());

    // Cell volume
    TGeoVolume* aCell = new TGeoVolume(aCellName.c_str(), aCellCs, medium);

    // "b"-type cells
    //
    // Initial placement:
    //
    // y
    // ^
    // |
    // |        *
    // |      * * *
    // |    * * b * *
    // |      * * * *
    // |
    // 0--------------> x

    std::string bCellName = "FV0" + cellType + sCellName + "b" + std::to_string(ir + 1);

    LOG(INFO) << "FV0 Geometry::initializeCells(): Initializing cell " << bCellName;

    // Base shape
    std::string bCellShapeName = bCellName + "Shape";

    // The cells in the innermost ring are slightly different than the rest
    if (ir == 0) {
      // The innermost "b"-type cell
      std::string b1CellShapeFullName = bCellShapeName + "Full";
      std::string b1CellShapeHoleCutName = bCellShapeName + "Cut";
      std::string b1CellShapeHoleCutTransName = b1CellShapeHoleCutName + "Trans";

      new TGeoTubeSeg(b1CellShapeFullName.c_str(), 0, mRMaxScint[ir], zThickness / 2 - sEpsilon, 0, 45);
      new TGeoTube(b1CellShapeHoleCutName.c_str(), 0, mRMinScint[ir], zThickness);

      createAndRegisterTrans(b1CellShapeHoleCutTransName, sXShiftInnerRadiusScint, 0, 0);

      std::string b1BoolFormula = b1CellShapeFullName + "-" + b1CellShapeHoleCutName + ":" + b1CellShapeHoleCutTransName;
      new TGeoCompositeShape(bCellShapeName.c_str(), b1BoolFormula.c_str());
    } else {
      // The rest of the "b"-type cells
      new TGeoTubeSeg(bCellShapeName.c_str(), mRMinScint[ir], mRMaxScint[ir], zThickness / 2, 0, 45);
    }

    // Translations for holes
    //
    // 1 = outer left
    // 2 = inner left
    // 3 = outer right
    // 4 = inner right
    // 5 = outer middle
    // 6 = inner middle
    // 7 = half-lenght left
    // 8 = half-length right

    std::string bHole1TransName = bCellName + "Hole1Trans";
    std::string bHole2TransName = bCellName + "Hole2Trans";
    std::string bHole3TransName = bCellName + "Hole3Trans";
    std::string bHole4TransName = bCellName + "Hole4Trans";
    std::string bHole5TransName = bCellName + "Hole5Trans";
    std::string bHole6TransName = bCellName + "Hole6Trans";
    std::string bHole7TransName = bCellName + "Hole7Trans";
    std::string bHole8TransName = bCellName + "Hole8Trans";

    createAndRegisterTrans(bHole1TransName, sin(45 * M_PI/180) * rMax, cos(45 * M_PI/180) * rMax, 0);
    createAndRegisterTrans(bHole2TransName, sin(45 * M_PI/180) * rMin, cos(45 * M_PI/180) * rMin, 0);
    createAndRegisterTrans(bHole3TransName, rMax, 0, 0);
    createAndRegisterTrans(bHole4TransName, rMin, 0, 0);
    createAndRegisterTrans(bHole5TransName, cos(22.5 * M_PI/180) * rMax, sin(22.5 * M_PI/180) * rMax, 0);
    createAndRegisterTrans(bHole6TransName, cos(22.5 * M_PI/180) * rMin, sin(22.5 * M_PI/180) * rMin, 0);
    createAndRegisterTrans(bHole7TransName, sin(45 * M_PI/180) * rMid, cos(45 * M_PI/180) * rMid, 0);
    createAndRegisterTrans(bHole8TransName, rMid, 0, 0);

    // Composite shape
    std::string bBoolFormula = bCellShapeName;

    // sector separation
    bBoolFormula += "-" + secSepShapeName;
    bBoolFormula += "-" + secSepShapeName + ":" + secSepRot45Name;

    // outer holes
    bBoolFormula += "-" + ((ir < 2) ? holeSmallName : holeLargeName) + ":" + bHole1TransName;
    bBoolFormula += "-" + ((ir < 2) ? holeSmallName : holeLargeName) + ":" + bHole3TransName;

    // inner holes
    if (ir > 0) {
      std::string holeName = (ir < 3) ? holeSmallName : holeLargeName;

      bBoolFormula += "-" + holeName + ":" + bHole2TransName;
      bBoolFormula += "-" + holeName + ":" + bHole4TransName;
    }

    // outer middle hole
    if (ir > 1) {
      bBoolFormula += "-" + holeLargeName + ":" + bHole5TransName;
    }

    // inner middle hole
    if (ir > 2) {
      bBoolFormula += "-" + holeLargeName + ":" + bHole6TransName;
    }

    // half-lenght holes
    if (ir == 4) {
      bBoolFormula += "-" + holeLargeName + ":" + bHole7TransName;
      bBoolFormula += "-" + holeLargeName + ":" + bHole8TransName;
    }

    std::string bCellCSName = bCellName + "CS";
    TGeoCompositeShape* bCellCs = new TGeoCompositeShape(bCellCSName.c_str(), bBoolFormula.c_str());

    // Cell volume
    TGeoVolume* bCell = new TGeoVolume(bCellName.c_str(), bCellCs, medium);

    if (isSensitive) {
      mSensitiveVolumeNames.push_back(aCell->GetName());
      mSensitiveVolumeNames.push_back(bCell->GetName());
    }
  }
}

void Geometry::initializeScintCells()
{
  TGeoMedium* medium = gGeoManager->GetMedium("FV0_Scintillator$");
  initializeCells(sScintName, sDzScint, medium, true);
}

void Geometry::initializePlasticCells()
{
  TGeoMedium* medium = gGeoManager->GetMedium("FV0_Plastic$");
  initializeCells(sPlastName, sDzPlast, medium, false);
}

void Geometry::initializeFibers()
{
  LOG(INFO) << "FVO Geometry::initializeFibers(): Initializing fibers";

  float dzFibers = sDzAlu - sDzAluBack - sDzAluFront - sDzScint - sDzPlast - 2 * sEpsilon;  // depth of the fiber volumes

  std::string fiberName = "FV0_Fibers";   // No volume with this name

  std::string fiberSepCutName = fiberName + "SepCut";
  std::string fiberConeCutName = fiberName + "ConeCut";
  std::string fiberHoleCutName = fiberName + "HoleCut";

  std::string fiberTransName = fiberName + "Trans";
  std::string fiberConeCutTransName = fiberConeCutName + "Trans";
  std::string fiberHoleCutTransName = fiberHoleCutName + "Trans";

  new TGeoBBox(fiberSepCutName.c_str(), sDrSeparationScint, mRMaxFiber.back() + sEpsilon, dzFibers / 2 + sEpsilon );
  new TGeoConeSeg(fiberConeCutName.c_str(), sDzAluCone / 2 + sEpsilon, 0, sDrMinAluCone + sXYThicknessAluCone + sEpsilon, 0, sDrMinAluFront + sEpsilon, -90, 90);
  new TGeoTube(fiberHoleCutName.c_str(), 0, mRMinScint.front(), dzFibers / 2 + sEpsilon);

  createAndRegisterTrans(fiberTransName, sPosScint[0] , 0, sZFiber);
  createAndRegisterTrans(fiberConeCutTransName, sPosScint[0], 0, sZCone);
  createAndRegisterTrans(fiberHoleCutTransName, sPosScint[0] + sXShiftInnerRadiusScint, 0, sZFiber);

  for (int i = 0; i < mRMinFiber.size(); i++) {
    std::string fiberShapeName = fiberName + std::to_string(i + 1);

    LOG(INFO) << "FV0 Geometry::initializeFibers(): Initializing fiber volume " << fiberShapeName;

    new TGeoTubeSeg(fiberShapeName.c_str(), mRMinFiber[i], mRMaxFiber[i] - sEpsilon, dzFibers / 2, -90, 90);

    // Composite shape
    std::string boolFormula = "";
    boolFormula += fiberShapeName + ":" + fiberTransName;
    boolFormula += "-" + fiberSepCutName + ":" + fiberTransName;
    boolFormula += "-" + fiberConeCutName + ":" + fiberConeCutTransName;

    if (i == 0) {
      // Cut out the hole in the innermost fiber volume
      boolFormula += "-" + fiberHoleCutName + ":" + fiberHoleCutTransName;
    }

    // Remove holes for screws and rods
    boolFormula += "-" + sScrewHolesCSName;
    boolFormula += "-" + sRodHolesCSName;

    TGeoCompositeShape* fiberCS = new TGeoCompositeShape((fiberShapeName + "CS").c_str(), boolFormula.c_str());

    // Volume
    // TODO: Add check for medium?
    new TGeoVolume(("FV0" + sFiberName + std::to_string(i + 1)).c_str(), fiberCS, mMediumFiber[i]);
  }
}

void Geometry::initializeScrews()
{
  LOG(INFO) << "FV0 Geometry::initializeScrews(): Initializing screws";

  for (int i = 0; i < sNScrewTypes; i++) {
    std::string screwName = "FV0" + sScrewName + std::to_string(i);

    TGeoShape* screwShape = createScrewShape(screwName + "Shape", i, 0, 0, 0);

    // TODO: Add check for medium?

    TGeoMedium* med = mMediumScrewTypes[i];
    if (!med) {
      LOG(INFO) << "ANDREAS med not found";
    }
    new TGeoVolume(screwName.c_str(), screwShape, mMediumScrewTypes[i]);
  }
}

void Geometry::initializeRods()
{
  LOG(INFO) << "FV0 Geometry::initializeRods(): Initializing rods";

  for (int i = 0; i < sNRodTypes; i++) {
    std::string rodName = "FV0" + sRodName + std::to_string(i);

    TGeoShape* rodShape = createRodShape(rodName + "Shape", i, -sEpsilon, -sEpsilon);

    // TODO: Add check for medium?
    new TGeoVolume(rodName.c_str(), rodShape, mMediumScrewTypes[i]);
  }
}

void Geometry::initializeMetalContainer()
{
  // The metal container is constructed starting from the backplate. The backplate is positioned first, relative to
  // the scintillator cells. The rest of the container parts are positioned relative to the backplate.

  // TODO: Make position variables consistent, some are now global coordinates, and some are relative to some other part of the container

  // Backplate
  std::string backPlateName = "FV0_BackPlate";                        // the full backplate
  std::string backPlateStandName = backPlateName + "Stand";           // the stand part of the backplate
  std::string backPlateHoleName = backPlateName + "Hole";             // the hole in the middle of the backplate
  std::string backPlateHoleCutName = backPlateHoleName + "Cut";       // extension of the hole
  std::string backPlateStandTransName = backPlateStandName + "Trans"; // shift of the backplate stand
  std::string backPlateHoleTransName = backPlateHoleName + "Trans";   // shift of the backplate inner radius

  new TGeoTubeSeg(backPlateName.c_str(), 0, sDrMaxAluBack, sDzAluBack / 2, -90, 90);
  new TGeoBBox(backPlateStandName.c_str(), sDxAluStand / 2, (sDrMaxAluBack + sDyAluStand) / 2, sDzAluBack / 2);
  new TGeoTubeSeg(backPlateHoleName.c_str(), 0, sDrAluHole, sDzAluBack / 2, -90, 90);
  new TGeoBBox(backPlateHoleCutName.c_str(), -sXShiftAluHole, sDrAluHole, sDzAluBack);

  createAndRegisterTrans(backPlateStandTransName, sDxAluStand / 2, - (sDrMaxAluBack + sDyAluStand) / 2, 0);
  createAndRegisterTrans(backPlateHoleTransName, sXShiftAluHole, 0, 0);

  // Backplate composite shape
  std::string backPlateBoolFormula = "";
  backPlateBoolFormula += backPlateName;
  backPlateBoolFormula += "+" + backPlateStandName + ":" + backPlateStandTransName;
  backPlateBoolFormula += "-" + backPlateHoleName + ":" + backPlateHoleTransName;
  backPlateBoolFormula += "-" + backPlateHoleCutName;

  std::string backPlateCSName = backPlateName + "CompositeShape";
  std::string backPlateCSTransName = backPlateCSName + "Trans";

  new TGeoCompositeShape(backPlateCSName.c_str(), backPlateBoolFormula.c_str());
  createAndRegisterTrans(backPlateCSTransName, 0, 0, sZAluBack);

  // Frontplate
  float zPosFrontPlate = sZAluFront;                                                    // the z-position o the frontplate
  float dyFrontPlateStand = sDyAluStand + (sDrMaxAluFront - sDrMinAluFront) / 2;        // the height of the total stand overlapping with the rest of the plate
  float yPosFrontPlateStand = -sDrMaxAluFront - sDyAluStand + dyFrontPlateStand / 2;    // the y-position of the total stand

  std::string frontPlateName = "FV0_FrontPlate";
  std::string frontPlateStandName = frontPlateName + "Stand";
  std::string frontPlateTransName = frontPlateName + "Trans";
  std::string frontPlateStandTransName = frontPlateStandName + "Trans";

  new TGeoTubeSeg(frontPlateName.c_str(), sDrMinAluFront, sDrMaxAluFront, sDzAluFront / 2 , -90, 90);
  new TGeoBBox(frontPlateStandName.c_str(), sDxAluStand / 2, dyFrontPlateStand / 2, sDzAluBack / 2);

  createAndRegisterTrans(frontPlateTransName, 0, 0, zPosFrontPlate);
  createAndRegisterTrans(frontPlateStandTransName, sDxAluStand / 2, yPosFrontPlateStand, 0);

  // Frontplate cone composite shape
  std::string frontPlateBoolFormula = "";
  frontPlateBoolFormula += frontPlateName;
  frontPlateBoolFormula += "+" + frontPlateStandName + ":" + frontPlateStandTransName;

  std::string frontPlateCSName = frontPlateName + "CompositeName";

  new TGeoCompositeShape(frontPlateCSName.c_str(), frontPlateBoolFormula.c_str());

  // Frontplate cone
  float thicknessFrontPlateCone = sXYThicknessAluCone;      // radial thickness of frontplate cone in the xy-plane
  float zPosCone = sDzAluFront / 2 - sDzAluCone / 2;       // z-position of the frontplate cone relative to the frontplate

  std::string frontPlateConeName = "FV0_FrontPlateCone";                            // no volume with this name
  std::string frontPlateConeShieldName = frontPlateConeName + "Shield";             // the "sides" of the cone
  std::string frontPlateConeShieldTransName = frontPlateConeShieldName + "Trans";

  new TGeoConeSeg(frontPlateConeShieldName.c_str(),
                  sDzAluCone / 2,
                  sDrMinAluCone,
                  sDrMinAluCone + thicknessFrontPlateCone,
                  sDrMinAluFront - thicknessFrontPlateCone,
                  sDrMinAluFront,
                  -90,
                  90);
  createAndRegisterTrans(frontPlateConeShieldTransName, 0, 0, zPosCone);

  // Frontplate cone "bottom"
  float zPosConePlate = sDzAluFront / 2 - sDzAluCone + thicknessFrontPlateCone / 2;   // z-position of the cone bottom relative to the frontplate
  std::string frontPlateConePlateName = frontPlateConeName + "Plate";                  // the bottom of the cone

  new TGeoTubeSeg(frontPlateConePlateName.c_str(), 0, sDrMinAluCone + thicknessFrontPlateCone,
                  thicknessFrontPlateCone / 2, -90, 90);

  // Frontplate cone bottom composite shape
  std::string frontPlateConePlateCSBoolFormula;
  frontPlateConePlateCSBoolFormula += frontPlateConePlateName;
  frontPlateConePlateCSBoolFormula += "-" + backPlateHoleName + ":" + backPlateHoleTransName;

  std::string frontPlateConePlateCSName = frontPlateConePlateName + "CompositeShape";
  std::string frontPlateConePlateCSTransName = frontPlateConePlateCSName + "Trans";
  new TGeoCompositeShape(frontPlateConePlateCSName.c_str(), frontPlateConePlateCSBoolFormula.c_str());
  createAndRegisterTrans(frontPlateConePlateCSTransName, 0, 0, zPosConePlate);

  // Frontplate cone composite shape
  std::string frontPlateConeCSBoolFormula = "";
  frontPlateConeCSBoolFormula += frontPlateConeShieldName + ":" + frontPlateConeShieldTransName;
  frontPlateConeCSBoolFormula += "+" + frontPlateConePlateCSName + ":" + frontPlateConePlateCSTransName;

  std::string frontPlateConeCSName = frontPlateConeName + "CompositeShape";
  new TGeoCompositeShape(frontPlateConeCSName.c_str(), frontPlateConeCSBoolFormula.c_str());

  // Shields
  float dzShieldGap = 0.7;                      // z-distance between the shields and the front- and backplate outer edges (in z-direction)
  float dzShield = sDzAlu - 2 * dzShieldGap;    // depth of the shields

  // Outer shield
  float zPosOuterShield = (sZAluBack + sZAluFront) / 2;   // z-position of the outer shield

  std::string outerShieldName = "FV0_OuterShield";
  std::string outerShieldTransName = outerShieldName + "Trans";

  new TGeoTubeSeg(outerShieldName.c_str(), sDrMinAluOuterShield, sDrMaxAluOuterShield, dzShield / 2, -90, 90);
  createAndRegisterTrans(outerShieldTransName, 0, 0, zPosOuterShield);

  // Inner shield
  float dzInnerShield = sDzAlu - sDzAluCone - dzShieldGap;                              // depth of the inner shield
  float zPosInnerShield = sZAluBack - sDzAluBack / 2 + dzShieldGap + dzInnerShield / 2; // z-position of the inner shield relative to the backplate

  std::string innerShieldName = "FV0_InnerShield";
  std::string innerShieldCutName = innerShieldName + "Cut";

  new TGeoTubeSeg(innerShieldName.c_str(), sDrMinAluInnerShield, sDrMaxAluInnerShield, dzInnerShield / 2, -90, 90);
  new TGeoBBox(innerShieldCutName.c_str(), fabs(sXShiftAluHole), sDrMaxAluInnerShield, dzInnerShield / 2);

  // Inner shield composite shape
  std::string innerShieldCSBoolFormula;
  innerShieldCSBoolFormula = innerShieldName;
  innerShieldCSBoolFormula += "-" + innerShieldCutName;

  std::string innerShieldCSName = innerShieldName + "CS";
  std::string innerShieldCSTransName = innerShieldCSName + "Trans";
  new TGeoCompositeShape(innerShieldCSName.c_str(), innerShieldCSBoolFormula.c_str());
  createAndRegisterTrans(innerShieldCSTransName, sXShiftAluHole, 0, zPosInnerShield);

  // Cover
  float dzCover = sDzAlu;                             // Depth of the covers
  float zPosCoverConeCut = zPosFrontPlate + zPosCone; // Set the cone cut relative to the frontplate so that the exact position of the aluminium cone part can be used.

  std::string coverName = "FV0_Cover";
  std::string coverConeCutName = coverName + "ConeCut";
  std::string coverHoleCutName = coverName + "HoleCut";

  new TGeoBBox(coverName.c_str(), sDxAluCover / 2, sDrMaxAluOuterShield, dzCover / 2);
  new TGeoCone(coverConeCutName.c_str(), sDzAluCone / 2, 0, sDrMinAluCone + thicknessFrontPlateCone, 0, sDrMinAluFront);
  new TGeoTubeSeg(coverHoleCutName.c_str(), 0, sDrMinAluInnerShield, dzCover / 2, 0, 360);

  std::string coverTransName = coverName + "Trans";
  std::string coverConeCutTransName = coverConeCutName + "Trans";
  std::string coverHoleCutTransName = coverHoleCutName + "Trans";

  createAndRegisterTrans(coverTransName, sDxAluCover / 2, 0, zPosOuterShield);
  createAndRegisterTrans(coverConeCutTransName, 0, 0, zPosCoverConeCut);
  createAndRegisterTrans(coverHoleCutTransName.c_str(), sXShiftAluHole, 0, zPosOuterShield);

  // Cover composite shape
  std::string coverCSBoolFormula = "";
  coverCSBoolFormula += coverName + ":" + coverTransName;
  coverCSBoolFormula += "-" + coverConeCutName + ":" + coverConeCutTransName;
  coverCSBoolFormula += "-" + coverHoleCutName + ":" + coverHoleCutTransName;

  std::string coverCSName = coverName + "CS";
  new TGeoCompositeShape(coverCSName.c_str(), coverCSBoolFormula.c_str());

  // Stand bottom
  float dzStandBottom = sDzAlu - sDzAluBack - sDzAluFront;
  float dyStandBottomGap = 0.5;                             // This bottom part is not vertically aligned with the "front and backplate stands"
  float dxStandBottomHole = 9.4;
  float dzStandBottomHole = 20.4;
  float dxStandBottomHoleSpacing = 3.1;

  std::string standName = "FV0_StandBottom";
  std::string standHoleName = standName + "Hole";

  new TGeoBBox(standName.c_str(), sDxAluStandBottom / 2, sDyAluStandBottom / 2, dzStandBottom / 2);
  new TGeoBBox(standHoleName.c_str(), dxStandBottomHole / 2, sDyAluStandBottom / 2 + sEpsilon, dzStandBottomHole / 2);

  std::string standHoleTrans1Name = standHoleName + "Trans1";
  std::string standHoleTrans2Name = standHoleName + "Trans2";
  std::string standHoleTrans3Name = standHoleName + "Trans3";

  createAndRegisterTrans(standHoleTrans1Name, -dxStandBottomHoleSpacing - dxStandBottomHole, 0, 0);
  createAndRegisterTrans(standHoleTrans2Name, 0, 0, 0);
  createAndRegisterTrans(standHoleTrans3Name, dxStandBottomHoleSpacing + dxStandBottomHole, 0, 0);

  // Stand bottom composite shape
  std::string standCSName = standName + "CS";

  std::string standBoolFormula = "";
  standBoolFormula += standName;
  standBoolFormula += "-" + standHoleName + ":" + standHoleTrans1Name;
  standBoolFormula += "-" + standHoleName + ":" + standHoleTrans2Name;
  standBoolFormula += "-" + standHoleName + ":" + standHoleTrans3Name;

  new TGeoCompositeShape(standCSName.c_str(), standBoolFormula.c_str());

  std::string standCSTransName = standCSName + "Trans";

  createAndRegisterTrans(standCSTransName.c_str(),
                        sDxAluStand - sDxAluStandBottom / 2,
                        -(sDrMaxAluBack + sDyAluStand) + sDyAluStandBottom / 2 + dyStandBottomGap,
                        sZAluMid);

  // Composite shape
  std::string boolFormula = "";
  boolFormula += backPlateCSName + ":" + backPlateCSTransName;
  boolFormula += "+" + frontPlateCSName + ":" + frontPlateTransName;
  boolFormula += "+" + frontPlateConeCSName + ":" + frontPlateTransName;
  boolFormula += "+" + outerShieldName + ":" + outerShieldTransName;
  boolFormula += "+" + innerShieldCSName + ":" + innerShieldCSTransName;
  boolFormula += "+" + coverCSName;
  boolFormula += "+" + standCSName + ":" + standCSTransName;
  boolFormula += "-" + sScrewHolesCSName;                                 // Remove holes for screws
  boolFormula += "-" + sRodHolesCSName;                                   // Remove holes for rods

  std::string aluContCSName = "FV0_AluContCS";
  TGeoCompositeShape* aluContCS = new TGeoCompositeShape(aluContCSName.c_str(), boolFormula.c_str());

  // Volume
  std::string aluContName = "FV0" + sContainerName;
  TGeoMedium* medium = gGeoManager->GetMedium("FV0_Aluminium$");
  new TGeoVolume(aluContName.c_str(), aluContCS, medium);
}

void Geometry::assembleSensVols(TGeoVolume* vFV0)
{
  if (mEnabledComponents[eScintillator]) {
    assembleScintSectors(vFV0);
  }
}

void Geometry::assembleNonSensVols(TGeoVolume* vFV0)
{
  if (mEnabledComponents[ePlastics]) {
    assemblePlasticSectors(vFV0);
  }
  if (mEnabledComponents[eFibers]) {
    assembleFibers(vFV0);
  }
  if (mEnabledComponents[eScrews]) {
    assembleScrews(vFV0);
  }
  if (mEnabledComponents[eRods]) {
    assembleRods(vFV0);
  }
  if (mEnabledComponents[eAluminiumContainer]) {
    assembleMetalContainer(vFV0);
  }
}

void Geometry::assembleScintSectors(TGeoVolume* vFV0)
{
  TGeoVolumeAssembly* sectors = buildSectorAssembly(sScintName);
  vFV0->AddNode(sectors, 1);
}

void Geometry::assemblePlasticSectors(TGeoVolume* vFV0)
{
  TGeoVolumeAssembly* sectors = buildSectorAssembly(sPlastName);
  vFV0->AddNode(sectors, 1, new TGeoTranslation(0, 0, sZPlast));
}

void Geometry::assembleFibers(TGeoVolume* vFV0)
{
  TGeoVolumeAssembly* fibers = new TGeoVolumeAssembly("FV0FIBERS");

  for (int i = 0; i < mRMinFiber.size(); i++) {
    TGeoVolume* fiber = gGeoManager->GetVolume(("FV0" + sFiberName + std::to_string(i + 1)).c_str());
    if (!fiber) {
      LOG(WARNING) << "FV0 Geometry::assembleFibers(): Fiber volume no. " << i + 1 << " not found.";
    }
    fibers->AddNode(fiber, i);
  }

  TGeoRotation* reflection = new TGeoRotation();
  reflection->ReflectX(true);

  vFV0->AddNode(fibers, 1);
  vFV0->AddNode(fibers, 2, reflection);
}

void Geometry::assembleScrews(TGeoVolume* vFV0)
{
  LOG(INFO) << "FV0 Geometry::assembleScrews(): Assembling screws";

  TGeoVolumeAssembly* screws = new TGeoVolumeAssembly("FV0SCREWS");

  // TODO: Add check for screw initialization?

  for (int i = 0; i < mScrewPos.size(); i++) {
    TGeoVolume* screw = gGeoManager->GetVolume(("FV0" + sScrewName + std::to_string(mScrewTypeIDs[i])).c_str());
    if (!screw) {
      LOG(INFO) << "FV0 Geometry::assembleScrews(): Screw not found";
    } else {
      screws->AddNode(screw, i, new TGeoTranslation(mScrewPos[i][0] + sXShiftScrews, mScrewPos[i][1], mScrewPos[i][2]));
    }
  }

  TGeoRotation* reflection = new TGeoRotation();
  reflection->ReflectX(true);

  vFV0->AddNode(screws, 1);
  vFV0->AddNode(screws, 2, reflection);
}

void Geometry::assembleRods(TGeoVolume* vFV0)
{
  LOG(INFO) << "FV0 Geometry::assembleRods(): Assembling rods";

  TGeoVolumeAssembly* rods = new TGeoVolumeAssembly("FV0RODS");

  // TODO: Add check for rod initialization?

  for (int i = 0; i < mRodPos.size(); i++) {
    TGeoVolume* rod = gGeoManager->GetVolume(("FV0" + sRodName + std::to_string(mRodTypeIDs[i])).c_str());

    if (!rod) {
      LOG(INFO) << "FV0 Geometry::assembleRods(): Rod not found";
    } else {
      rods->AddNode(rod, i, new TGeoTranslation(mRodPos[i][0] + sXShiftScrews, mRodPos[i][1], mRodPos[i][2]));
    }
  }

  TGeoRotation* reflection = new TGeoRotation();
  reflection->ReflectX(true);

  vFV0->AddNode(rods, 1);
  vFV0->AddNode(rods, 2, reflection);
}

void Geometry::assembleMetalContainer(TGeoVolume* volV0)
{
  std::string containerName = "FV0" + sContainerName;
  TGeoVolume* container = gGeoManager->GetVolume(containerName.c_str());
  if (!container) {
    LOG(WARNING) << "FV0: Couldn't find volume " << containerName;
  } else {
    TGeoRotation* reflection = new TGeoRotation();
    reflection->ReflectX(true);

    volV0->AddNode(container, 1);
    volV0->AddNode(container, 2, reflection);
  }
}

TGeoVolumeAssembly* Geometry::buildSectorAssembly(std::string cellName)
{
  std::string assemblyName = "FV0" + cellName;
  TGeoVolumeAssembly* assembly = new TGeoVolumeAssembly(assemblyName.c_str());
  TGeoVolumeAssembly* left = new TGeoVolumeAssembly((assemblyName + "LEFT").c_str());
  TGeoVolumeAssembly* right = new TGeoVolumeAssembly((assemblyName + "RIGHT").c_str());

  for (int iSector = 0; iSector < ceil(mSectorTrans.size() / 2); iSector++) {
    TGeoVolumeAssembly* sector = buildSector(cellName, iSector);
    left->AddNode(sector, iSector + 1, mSectorTrans.at(iSector));
  }

  for (int iSector = ceil(mSectorTrans.size() / 2); iSector < mSectorTrans.size(); iSector++) {
    TGeoVolumeAssembly* sector = buildSector(cellName, iSector);
    right->AddNode(sector, iSector + 1, mSectorTrans.at(iSector));
  }

  assembly->AddNode(left, 1, new TGeoTranslation(-sDxAluCover, 0, 0));
  assembly->AddNode(right, 2, new TGeoTranslation(sDxAluCover, 0, 0));

  return assembly;
}

TGeoVolumeAssembly* Geometry::buildSector(std::string cellType, int iSector)
{
  std::stringstream ssSectorName;
  ssSectorName << "FV0" << cellType << sSectorName << iSector + 1;

  LOG(DEBUG) << "FV0 Geometry::buildSector(): building sector " << ssSectorName.str();

  TGeoVolumeAssembly* sector = new TGeoVolumeAssembly(ssSectorName.str().c_str());

  for (int i = 0; i < sNCellRings; i++) {
    std::stringstream ssCellName;
    ssCellName << "FV0" << cellType << sCellName << sCellTypes[iSector] << i + 1;

    TGeoVolume* cell = gGeoManager->GetVolume(ssCellName.str().c_str());

    if (!cell) {
      LOG(WARNING) << "FV0 Geometry::buildSector(): Couldn't find cell volume " << ssCellName.str();
    } else {
      LOG(DEBUG) << "FV0 Geometry::buildSector(): adding cell volume " << ssCellName.str();
      sector->AddNode(cell, i + 1);
    }
  }

  return sector;
}

TGeoShape* Geometry::createScrewShape(std::string shapeName, int screwTypeID, float xEpsilon, float yEpsilon, float zEpsilon)
{
  float xyEpsilon = (fabs(xEpsilon) > fabs(yEpsilon)) ? xEpsilon : yEpsilon;
  float dzMax = sDzMaxScrewTypes[screwTypeID] / 2 + zEpsilon;
  float dzMin = sDzMinScrewTypes[screwTypeID] / 2 + zEpsilon;

  std::string thinPartName = shapeName + "Thin";
  std::string thickPartName = shapeName + "Thick";
  std::string thickPartTransName = thickPartName + "Trans";

  TGeoTube* thinPart = new TGeoTube(thinPartName.c_str(), 0, sDrMinScrewTypes[screwTypeID] + xyEpsilon, dzMax);
  TGeoTube* thickPart = new TGeoTube(thickPartName.c_str(), 0, sDrMaxScrewTypes[screwTypeID] + xyEpsilon, dzMin);
  createAndRegisterTrans(thickPartTransName, 0, 0, -dzMax - sZShiftScrew + sDzScint + sDzPlast + dzMin);

  std::string boolFormula = thinPartName;
  boolFormula += "+" + thickPartName + ":" + thickPartTransName;

  TGeoCompositeShape* screwShape = new TGeoCompositeShape(shapeName.c_str(), boolFormula.c_str());
  return screwShape;
}

TGeoShape* Geometry::createRodShape(std::string shapeName, int rodTypeID, float xEpsilon, float yEpsilon, float zEpsilon)
{
  float dxMin = sDxMinRodTypes[rodTypeID] / 2 + xEpsilon;
  float dxMax = sDxMaxRodTypes[rodTypeID] / 2 + xEpsilon;
  float dyMin = sDyMinRodTypes[rodTypeID] / 2 + yEpsilon;
  float dyMax = sDyMaxRodTypes[rodTypeID] / 2 + yEpsilon;
  float dzMax = sDzMaxRodTypes[rodTypeID] / 2 + zEpsilon;
  float dzMin = sDzMinRodTypes[rodTypeID] / 2 + zEpsilon;

  std::string thinPartName = shapeName + "Thin";
  std::string thickPartName = shapeName + "Thick";
  std::string thickPartTransName = thickPartName + "Trans";

  TGeoBBox* thinPart = new TGeoBBox(thinPartName.c_str(), dxMin, dyMin, dzMax);
  TGeoBBox* thickPart = new TGeoBBox(thickPartName.c_str(), dxMax, dyMax, dzMin);
  createAndRegisterTrans(thickPartTransName, dxMax - dxMin, 0, -dzMax - sZShiftRod + sDzScint + sDzPlast + dzMin);

  std::string boolFormula = thinPartName;
  boolFormula += "+" + thickPartName + ":" + thickPartTransName;
  TGeoCompositeShape* rodShape = new TGeoCompositeShape(shapeName.c_str(), boolFormula.c_str());
  return rodShape;
}

TGeoTranslation* Geometry::createAndRegisterTrans(std::string name, double dx, double dy, double dz)
{
  TGeoTranslation* trans = new TGeoTranslation(name.c_str(), dx, dy, dz);
  trans->RegisterYourself();
  return trans;
}

TGeoRotation* Geometry::createAndRegisterRot(std::string name, double phi, double theta, double psi)
{
  TGeoRotation* rot = new TGeoRotation(name.c_str(), phi, theta, psi);
  rot->RegisterYourself();
  return rot;
}