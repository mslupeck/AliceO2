// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include <iomanip>
//#include <TVector3.h>
#include "V0Base/Geometry.h"

#include <TGeoManager.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoCone.h>
#include <TGeoCompositeShape.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <FairLogger.h>
#include <cmath>
#include <sstream>

ClassImp(o2::v0::Geometry);

using namespace o2::v0;

Geometry::Geometry(EGeoType initType)
{
  mGeometryType = initType;
  initializeGeometry();
  buildGeometry();
}

Geometry::Geometry(const Geometry& geom)
{
  this->mGeometryType = geom.mGeometryType;
}

const int Geometry::getCurrentCellId(TVirtualMC* fMC) {
  int ring = -1;
  int sector = -1;

  fMC->CurrentVolOffID(0, ring);
  fMC->CurrentVolOffID(1, sector);

  return sector + 8 * (ring - 1);
}

void Geometry::initializeGeometry()
{
  initializeVectors();
  initializeSensVols();
  if (mGeometryType == eFull) {
    initializeNonSensVols();
  }
  initializeLuts();
}

void Geometry::initializeVectors()
{
  // RADII

  // Scintillator radii
  LOG(INFO) << "FV0 Geometry::initializeVectors(): Initializing V0 scintillator ring radii.";

  // Index of mvrAvgScint is NOT linked directly to any ring number
  mvrAvgScint.assign(sRingRadiiScint, sRingRadiiScint + sNumberOfRings + 1);

  // Set real scintillator radii (reduced by paint thickness and separation gap)
  for (uint16_t ir = 0; ir < mvrAvgScint.size() - 1; ir++) {
    mvrMinScint.push_back(mvrAvgScint.at(ir) + sDrSeparationScint);
    mvrMaxScint.push_back(mvrAvgScint.at(ir + 1) - sDrSeparationScint);
    LOG(INFO) << "FV0 Geometry::initializeVectors(): Scintillator ring " << ir << " min: " << mvrMinScint.at(ir) << ", max: " << mvrMaxScint.at(ir);
  }
  // Now indices of mvrMinScint and mvrMaxScint correspond to the same ring

  // SECTOR TRANSFORMATIONS
  LOG(INFO) << "FV0 Geometry::initializeVectors(): Initializing V0 sector transformations.";

  for (uint16_t isector = 1; isector <= sBaseNumberOfSectors; isector++) {
    // isector = 1 corresponds to the first sector counter-clockwise from the y-axis + global azimuthal rotation
    // isector = 2 corresponds to the next sector in counter-clockwise direction

    // if the sector is rotated, this is the angle (plus an extra 45 for "b"-type cells)
    float phi = isector * 45;

    std::stringstream ssRotName, ssReflName, ssTotalTransName;
    ssRotName << "FV0_rotationSector" << isector;
    ssReflName << "FV0_reflectSector" << isector;
    ssTotalTransName << "FV0_totalTransformSector" << isector;

    TGeoRotation* rot = new TGeoRotation(ssRotName.str().c_str());            // rotation of the sector
    TGeoRotation* refl = new TGeoRotation(ssReflName.str().c_str());          // reflection of the sector
    TGeoHMatrix* totTrans = new TGeoHMatrix(ssTotalTransName.str().c_str());  // the combined transformation

    // The reference to "a" and "b" can be understood with the CAD drawings of the detector.
    switch (isector)
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
    mvSectorTrans.push_back(totTrans);
  }

  // Fiber volume radii
  mvrMinFiber.push_back(0);
  mvrMinFiber.push_back(sDrMinAluCone + sEpsilon);
  mvrMinFiber.push_back(sDrMinAluFront + sEpsilon);

  mvrMaxFiber.push_back(sDrMinAluCone);
  mvrMaxFiber.push_back(sDrMinAluFront);
  mvrMaxFiber.push_back(mvrMaxScint.back());
}

void Geometry::initializeSensVols()
{
  initializeScintCells();
}

void Geometry::initializeNonSensVols()
{
  initializePlasticCells();
  initializeFibers();
  initializeMetalContainer();
}

void Geometry::initializeCells(std::string cellType, float zThickness, TGeoMedium* medium) {
  // Creating the two types of cells, "a" and "b", for each ring.
  // All sectors can be assembled with these cells.
  //
  // The reference to "a" and "b" can be understood with the CAD drawings of the detector.

  LOG(INFO) << "FV0 Geometry::initializeCells(): Initializing " << cellType << " cells.";
  
  float dz = zThickness / 2;                          // half depth of the cells
  float rHoleSmall = 0.265;                           // radius of the small holes
  float rHoleLarge = 0.415;                           // radius of the large holes
  float xHoleCut = 0.2;                               // width of extension of hole 1, 2 and 7 in the "a" cell
  float dxHole = sDrSeparationScint + xHoleCut;       // x-placement of holes 1, 2 and 7 in the "a" cell

  // Sector separation gap shape
  std::string secSepShapeName = "FV0_" + cellType + "SectorSeparation";
  new TGeoBBox(secSepShapeName.c_str(), mvrMaxScint.back(), sDrSeparationScint, dz);

  // Sector separation gap rotations
  std::string secSepRot45Name = "FV0_" + cellType + "SecSepRot45";
  std::string secSepRot90Name = "FV0_" + cellType + "SecSepRot90";

  TGeoRotation* secSepRot45 = new TGeoRotation(secSepRot45Name.c_str(), 45, 0, 0);  
  TGeoRotation* secSepRot90 = new TGeoRotation(secSepRot90Name.c_str(), 90, 0, 0);

  secSepRot45->RegisterYourself();
  secSepRot90->RegisterYourself();

  // Hole shapes
  std::string holeSmallName = "FV0_" + cellType + "HoleSmall";
  std::string holeLargeName = "FV0_" + cellType + "HoleLarge";
  std::string holeSmallCutName = "FV0_" + cellType + "HoleSmallCut";
  std::string holeLargeCutName = "FV0_" + cellType + "HoleLargeCut";

  new TGeoTube(holeSmallName.c_str(), 0, rHoleSmall, dz);  
  new TGeoTube(holeLargeName.c_str(), 0, rHoleLarge, dz);
  new TGeoBBox(holeSmallCutName.c_str(), xHoleCut, rHoleSmall, dz);
  new TGeoBBox(holeLargeCutName.c_str(), xHoleCut, rHoleLarge, dz);
  
  for (int ir = 0; ir < sNumberOfRings; ir++) {

    // Radii without separation
    float rMin = mvrAvgScint.at(ir);
    float rMax = mvrAvgScint.at(ir + 1);
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

    std::stringstream aCellName;
    aCellName << "FV0" << cellType << "a" << ir + 1;

    LOG(INFO) << "FV0 Geometry::initializeCells(): Initializing cell " << aCellName.str();

    // Base shape
    std::string aCellShapeName = aCellName.str() + "Shape";

    // The cells in the innermost ring have a slightly shifted inner radius origin.
    if (ir == 0) {
      // The innermost "a"-type cell
      std::string a1CellShapeFullName = aCellShapeName + "Full";
      std::string a1CellShapeCutName = aCellShapeName + "Cut";
      std::string a1CellShapeCutTransName = a1CellShapeCutName + "Trans";
      
      new TGeoTubeSeg(a1CellShapeFullName.c_str(), 0, mvrMaxScint.at(ir), dz, 45, 90);
      new TGeoTubeSeg(a1CellShapeCutName.c_str(), 0, mvrMinScint.at(ir), dz, 0, 360);
      TGeoTranslation* a1CellShapeCutTrans = new TGeoTranslation(a1CellShapeCutTransName.c_str(),
                                                                 sRingInnerRadiusDx, 0, 0);
      a1CellShapeCutTrans->RegisterYourself();

      std::string a1BoolFormula = a1CellShapeFullName + "-" + a1CellShapeCutName + ":" + a1CellShapeCutTransName;
      new TGeoCompositeShape(aCellShapeName.c_str(), a1BoolFormula.c_str());
    } else {
      // The rest of the "a"-type cells
      new TGeoTubeSeg(aCellShapeName.c_str(), mvrMinScint.at(ir), mvrMaxScint.at(ir), dz, 45, 90);
    }

    // Translations for screw holes
    std::string aHole1TransName = aCellName.str() + "Hole1Trans";
    std::string aHole2TransName = aCellName.str() + "Hole2Trans";
    std::string aHole3TransName = aCellName.str() + "Hole3Trans";
    std::string aHole4TransName = aCellName.str() + "Hole4Trans";
    std::string aHole5TransName = aCellName.str() + "Hole5Trans";
    std::string aHole6TransName = aCellName.str() + "Hole6Trans";
    std::string aHole7TransName = aCellName.str() + "Hole7Trans";
    std::string aHole8TransName = aCellName.str() + "Hole8Trans";
    std::string aHole1CutTransName = aCellName.str() + "Hole1CutTrans";
    std::string aHole2CutTransName = aCellName.str() + "Hole2CutTrans";
    std::string aHole7CutTransName = aCellName.str() + "Hole7CutTrans";

    // 1 = outer left
    // 2 = inner left
    // 3 = outer right
    // 4 = inner right
    // 5 = outer middle
    // 6 = inner middle
    // 7 = half-lenght left
    // 8 = half-length right

    // holes 1, 2 and 7 are sligtly shifted along the rim of the cell

    TGeoTranslation* aHole1Trans = new TGeoTranslation(aHole1TransName.c_str(), dxHole,
                                                       cos(asin(dxHole/rMax)) * rMax, 0);
    TGeoTranslation* aHole2Trans = new TGeoTranslation(aHole2TransName.c_str(), dxHole,
                                                       cos(asin(dxHole/rMin)) * rMin, 0);
    TGeoTranslation* aHole3Trans = new TGeoTranslation(aHole3TransName.c_str(), sin(45 * M_PI/180) * rMax,
                                                       cos(45 * M_PI/180) * rMax, 0);
    TGeoTranslation* aHole4Trans = new TGeoTranslation(aHole4TransName.c_str(), sin(45 * M_PI/180) * rMin,
                                                       cos(45 * M_PI/180) * rMin, 0);
    TGeoTranslation* aHole5Trans = new TGeoTranslation(aHole5TransName.c_str(), sin(22.5 * M_PI/180) * rMax,
                                                       cos(22.5 * M_PI/180) * rMax, 0);
    TGeoTranslation* aHole6Trans = new TGeoTranslation(aHole6TransName.c_str(), sin(22.5 * M_PI/180) * rMin,
                                                       cos(22.5 * M_PI/180) * rMin, 0);
    TGeoTranslation* aHole7Trans = new TGeoTranslation(aHole7TransName.c_str(), dxHole,
                                                       cos(asin(dxHole/rMid)) * rMid, 0);
    TGeoTranslation* aHole8Trans = new TGeoTranslation(aHole8TransName.c_str(), sin(45 * M_PI/180) * rMid,
                                                       cos(45 * M_PI/180) * rMid, 0);
    TGeoTranslation* aHole1CutTrans = new TGeoTranslation(aHole1CutTransName.c_str(), 0,
                                                          cos(asin(dxHole/rMax)) * rMax, 0);
    TGeoTranslation* aHole2CutTrans = new TGeoTranslation(aHole2CutTransName.c_str(), 0,
                                                          cos(asin(dxHole/rMin)) * rMin, 0);
    TGeoTranslation* aHole7CutTrans = new TGeoTranslation(aHole7CutTransName.c_str(), 0,
                                                          cos(asin(dxHole/rMid)) * rMid, 0);

    aHole1Trans->RegisterYourself();
    aHole2Trans->RegisterYourself();
    aHole3Trans->RegisterYourself();
    aHole4Trans->RegisterYourself();
    aHole5Trans->RegisterYourself();
    aHole6Trans->RegisterYourself();
    aHole7Trans->RegisterYourself();
    aHole8Trans->RegisterYourself();
    aHole1CutTrans->RegisterYourself();
    aHole2CutTrans->RegisterYourself();
    aHole7CutTrans->RegisterYourself();

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

      aBoolFormula.append("-" + screwHoleName + ":" + aHole2TransName);
      aBoolFormula.append("-" + screwHoleCutName + ":" + aHole2CutTransName);
      aBoolFormula.append("-" + screwHoleName + ":" + aHole4TransName);
    }

    // outer middle hole
    if (ir > 1) {
      aBoolFormula.append("-" + holeLargeName + ":" + aHole5TransName);
    }

    // inner middle hole
    if (ir > 2) {
      aBoolFormula.append("-" + holeLargeName + ":" + aHole6TransName);
    }

    // half-length holes
    if (ir == 4) {
      aBoolFormula.append("-" + holeLargeName + ":" + aHole7TransName);
      aBoolFormula.append("-" + holeLargeCutName + ":" + aHole7CutTransName);
      aBoolFormula.append("-" + holeLargeName + ":" + aHole8TransName);
    }

    LOG(INFO) << "FV0 Geometry::initializeCells(): Composite shape formula: " << aBoolFormula;

    std::string aCellCSName = aCellName.str() + "CS";
    TGeoCompositeShape* aCellCs = new TGeoCompositeShape(aCellCSName.c_str(), aBoolFormula.c_str());

    // Cell volume
    TGeoVolume* aCell = new TGeoVolume(aCellName.str().c_str(), aCellCs, medium);
    aCell->RegisterYourself();
    mvSensitiveVolumeNames.push_back(aCell->GetName());

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

    std::stringstream bCellName;
    bCellName << "FV0" << cellType << "b" << ir + 1;

    LOG(INFO) << "FV0 Geometry::initializeCells(): Initializing cell " << bCellName.str();

    // Base shape
    std::string bCellShapeName = bCellName.str() + "Shape";

    // The cells in the innermost ring are slightly different than the rest
    if (ir == 0) {
      // The innermost "b"-type cell
      std::string b1CellShapeFullName = bCellShapeName + "Full";
      std::string b1CellShapeCutName = bCellShapeName + "Cut";
      std::string b1CellShapeCutTransName = b1CellShapeCutName + "Trans";

      TGeoTubeSeg* b1CellShapeFull = new TGeoTubeSeg(b1CellShapeFullName.c_str(), 0, mvrMaxScint.at(ir), dz, 0, 45);
      TGeoTubeSeg* b1CellShapeCut = new TGeoTubeSeg(b1CellShapeCutName.c_str(), 0, mvrMinScint.at(ir), dz, 0, 360);
      TGeoTranslation* b1CellShapeCutTrans = new TGeoTranslation(b1CellShapeCutTransName.c_str(),
                                                                 sRingInnerRadiusDx, 0, 0);
      b1CellShapeCutTrans->RegisterYourself();

      std::string b1BoolFormula = b1CellShapeFullName + "-" + b1CellShapeCutName + ":" + b1CellShapeCutTransName;
      TGeoCompositeShape* b1CellShape = new TGeoCompositeShape(bCellShapeName.c_str(), b1BoolFormula.c_str());
    } else {
      // The rest of the "b"-type cells
      TGeoTubeSeg* bCellShape = new TGeoTubeSeg(bCellShapeName.c_str(), mvrMinScint.at(ir), mvrMaxScint.at(ir), dz, 0,
                                                45);
    }

    // Translations for holes
    std::string bHole1TransName = bCellName.str() + "Hole1Trans";
    std::string bHole2TransName = bCellName.str() + "Hole2Trans";
    std::string bHole3TransName = bCellName.str() + "Hole3Trans";
    std::string bHole4TransName = bCellName.str() + "Hole4Trans";
    std::string bHole5TransName = bCellName.str() + "Hole5Trans";
    std::string bHole6TransName = bCellName.str() + "Hole6Trans";
    std::string bHole7TransName = bCellName.str() + "Hole7Trans";
    std::string bHole8TransName = bCellName.str() + "Hole8Trans";

    // 1 = outer left
    // 2 = inner left
    // 3 = outer right
    // 4 = inner right
    // 5 = outer middle
    // 6 = inner middle
    // 7 = half-lenght left
    // 8 = half-length right

    TGeoTranslation* bHole1Trans = new TGeoTranslation(bHole1TransName.c_str(), sin(45 * M_PI/180) * rMax,
                                                       cos(45 * M_PI/180) * rMax, 0);
    TGeoTranslation* bHole2Trans = new TGeoTranslation(bHole2TransName.c_str(), sin(45 * M_PI/180) * rMin,
                                                       cos(45 * M_PI/180) * rMin, 0);
    TGeoTranslation* bHole3Trans = new TGeoTranslation(bHole3TransName.c_str(), rMax, 0, 0);
    TGeoTranslation* bHole4Trans = new TGeoTranslation(bHole4TransName.c_str(), rMin, 0, 0);
    TGeoTranslation* bHole5Trans = new TGeoTranslation(bHole5TransName.c_str(), cos(22.5 * M_PI/180) * rMax,
                                                       sin(22.5 * M_PI/180) * rMax, 0);
    TGeoTranslation* bHole6Trans = new TGeoTranslation(bHole6TransName.c_str(), cos(22.5 * M_PI/180) * rMin,
                                                       sin(22.5 * M_PI/180) * rMin, 0);
    TGeoTranslation* bHole7Trans = new TGeoTranslation(bHole7TransName.c_str(), sin(45 * M_PI/180) * rMid,
                                                       cos(45 * M_PI/180) * rMid, 0);
    TGeoTranslation* bHole8Trans = new TGeoTranslation(bHole8TransName.c_str(), rMid, 0, 0);
    
    bHole1Trans->RegisterYourself();
    bHole2Trans->RegisterYourself();
    bHole3Trans->RegisterYourself();
    bHole4Trans->RegisterYourself();
    bHole5Trans->RegisterYourself();
    bHole6Trans->RegisterYourself();
    bHole7Trans->RegisterYourself();
    bHole8Trans->RegisterYourself();

    // Composite shape
    std::string bBoolFormula = bCellShapeName;

    // sector separation
    bBoolFormula += "-" + secSepShapeName;
    bBoolFormula += "-" + secSepShapeName + ":" + secSepRot45Name;

    // outer holes
    bBoolFormula += "-" + holeSmallName + ":" + bHole1TransName;
    bBoolFormula += "-" + holeSmallName + ":" + bHole3TransName;

    // inner holes
    if (ir > 0) {
      std::string holeName = (ir < 3) ? holeSmallName : holeLargeName;

      bBoolFormula.append("-" + holeName + ":" + bHole2TransName);
      bBoolFormula.append("-" + holeName + ":" + bHole4TransName);
    }

    // outer middle hole
    if (ir > 1) {
      bBoolFormula.append("-" + holeLargeName + ":" + bHole5TransName);
    }

    // inner middle hole
    if (ir > 2) {
      bBoolFormula.append("-" + holeLargeName + ":" + bHole6TransName);
    }

    // half-lenght holes
    if (ir == 4) {
      bBoolFormula.append("-" + holeLargeName + ":" + bHole7TransName);
      bBoolFormula.append("-" + holeLargeName + ":" + bHole8TransName);
    }

    std::string bCellCSName = bCellName.str() + "CS";
    TGeoCompositeShape* bCellCs = new TGeoCompositeShape(bCellCSName.c_str(), bBoolFormula.c_str());

    // Cell volume
    TGeoVolume* bCell = new TGeoVolume(bCellName.str().c_str(), bCellCs, medium);
    bCell->RegisterYourself();
    mvSensitiveVolumeNames.push_back(bCell->GetName());
  }
}

void Geometry::initializeScintCells()
{
  TGeoMedium* kMed = gGeoManager->GetMedium("V0_Scintillator$");
  initializeCells(sScintCellName, sDzScint, kMed);
}

void Geometry::initializePlasticCells()
{
  TGeoMedium* kMed = gGeoManager->GetMedium("V0_Plastic$");
  initializeCells(sPlastCellName, sDzPlast, kMed);
}

void Geometry::initializeFibers()
{
  LOG(INFO) << "FVO Geometry::initializeFibers(): Initializing fibers";

  float dzFibers = sDzAlu - sDzAluBack - sDzAluFront - sDzScint - sDzPlast - 2 * sEpsilon;  // depth of the fiber volumes
  float zFibers = (sZPlast + sZAluFront) / 2;                                               // z-position of the fiber volumes
  int numberOfFiberVols = mvrMinFiber.size();

  TGeoMedium* medFiberInner = gGeoManager->GetMedium("FiberInner$");
  TGeoMedium* medFiberMiddle = gGeoManager->GetMedium("FiberMiddle$");
  TGeoMedium* medFiberOuter = gGeoManager->GetMedium("FiberOuter$");
  TGeoMedium* medFiber[numberOfFiberVols] = { medFiberInner, medFiberMiddle, medFiberOuter };

  std::string fiberName = "FV0_Fibers";                               // No volume with this name

  std::string fiberSepCutName = fiberName + "SepCut";
  std::string fiberConeCutName = fiberName + "ConeCut";
  std::string fiberHoleCutName = fiberName + "HoleCut";

  std::string fiberTransName = fiberName + "Trans";
  std::string fiberConeCutTransName = fiberConeCutName + "Trans";
  std::string fiberHoleCutTransName = fiberHoleCutName + "Trans";

  new TGeoBBox(fiberSepCutName.c_str(), sDrSeparationScint, mvrMaxFiber.back(), dzFibers / 2);
  new TGeoConeSeg(fiberConeCutName.c_str(), sDzAluCone / 2, 0, sDrMinAluFront, 0, sDrMinAluCone + sXYThicknessAluCone, -90, 90);
  new TGeoTube(fiberHoleCutName.c_str(), 0, mvrMinScint.front(), dzFibers / 2);
  
  TGeoTranslation* fiberTrans = new TGeoTranslation(fiberTransName.c_str(), 0 , 0, zFibers);
  TGeoTranslation* fiberConeCutTrans = new TGeoTranslation(fiberConeCutTransName.c_str(), 0, 0, sZCone);
  TGeoTranslation* fibersHoleCutTrans = new TGeoTranslation(fiberHoleCutTransName.c_str(), sRingInnerRadiusDx, 0, zFibers);
  
  fiberTrans->RegisterYourself();
  fiberConeCutTrans->RegisterYourself();
  fibersHoleCutTrans->RegisterYourself();

  for (int i = 0; i < numberOfFiberVols; i++) {
    std::stringstream fiberShapeName;
    fiberShapeName << fiberName << i + 1;

    new TGeoTubeSeg(fiberShapeName.str().c_str(), mvrMinFiber.at(i), mvrMaxFiber.at(i), dzFibers / 2, -90, 90);

    // Composite shape
    std::string boolFormula = "";
    boolFormula += fiberShapeName.str() + ":" + fiberTransName;
    boolFormula += "-" + fiberSepCutName + ":" + fiberTransName;
    boolFormula += "-" + fiberConeCutName + ":" + fiberConeCutTransName;

    if (i == 0) {
      boolFormula += "-" + fiberHoleCutName + ":" + fiberHoleCutTransName;
    }

    std::string fiberCSName = fiberShapeName.str() + "CS";
    TGeoCompositeShape* fiberCS = new TGeoCompositeShape(fiberCSName.c_str(), boolFormula.c_str());

    // Volume
    std::stringstream fiberName;
    fiberName << "FV0" << sFiberName << i + 1;
    if (!medFiber[i]) {
      LOG(WARN) << "FV0 Geometry::initializeFibers(): Medium for fiber volume no. " << i << " not found!";
    }
    new TGeoVolume(fiberName.str().c_str(), fiberCS, medFiber[i]);
  }
}

void Geometry::initializeMetalContainer()
{
  // The metal container is constructed starting from the backplate. The backplate is positioned first, relative to
  // the scintillator cells. The rest of the container parts are positioned relative to the backplate.

  // TODO: Make position variables consistent, some are now global coordinates, and some are relative to some other part of the container

  // Backplate
  float zPosBackPlate = sZAluBack;                                    // z-position of the backplate

  std::string backPlateName = "FV0_BackPlate";                        // the full backplate
  std::string backPlateHoleName = backPlateName + "Hole";             // the hole in the middle of the backplate
  std::string backPlateHoleCutName = backPlateHoleName + "Cut";       // extension of the hole
  std::string backPlateHoleTransName = backPlateHoleName + "Trans";   // shift of backplate inner radius

  new TGeoTubeSeg(backPlateName.c_str(), 0, sDrMaxAluBack, sDzAluBack / 2, -90, 90);
  new TGeoTubeSeg(backPlateHoleName.c_str(), 0, sDrAluHole, sDzAluBack / 2, -90, 90);
  new TGeoBBox(backPlateHoleCutName.c_str(), -sXAluHole, sDrAluHole, sDzAluBack);
  
  TGeoTranslation* backPlateHoleTrans = new TGeoTranslation(backPlateHoleTransName.c_str(), sXAluHole, 0, 0);
  backPlateHoleTrans->RegisterYourself();

  // Backplate composite shape
  std::string backPlateBoolFormula = "";
  backPlateBoolFormula += backPlateName;
  backPlateBoolFormula += "-" + backPlateHoleName + ":" + backPlateHoleTransName;
  backPlateBoolFormula += "-" + backPlateHoleCutName;

  std::string backPlateCSName = backPlateName + "CompositeShape";
  std::string backPlateCSTransName = backPlateCSName + "Trans";

  new TGeoCompositeShape(backPlateCSName.c_str(), backPlateBoolFormula.c_str());
  TGeoTranslation* backPlateTrans = new TGeoTranslation(backPlateCSTransName.c_str(), 0, 0, zPosBackPlate);
  backPlateTrans->RegisterYourself();

  // Frontplate
  float zPosFrontPlate = sZAluFront;   // the z-position o the frontplate

  std::string frontPlateName = "FV0_FrontPlate";
  std::string frontPlateTransName = frontPlateName + "Trans";

  new TGeoTubeSeg(frontPlateName.c_str(), sDrMinAluFront, sDrMaxAluFront, sDzAluFront / 2 , -90, 90);
  TGeoTranslation* frontPlateTrans = new TGeoTranslation(frontPlateTransName.c_str(), 0, 0, zPosFrontPlate);
  frontPlateTrans->RegisterYourself();

  // Frontplate cone
  // TODO: the cone has a thickness of 0.6, but now 0.6 is used as its thickness in the xy-plane. Calculate the real thickness the in xy-plane.
  float thicknessFrontPlateCone = sThicknessAluCone;        // thickness of frontplate cone
  float zPosCone = -sDzAluFront / 2 + sDzAluCone / 2;       // z-position of the frontplate cone relative to the frontplate

  std::string frontPlateConeName = "FV0_FrontPlateCone";                            // no volume with this name
  std::string frontPlateConeShieldName = frontPlateConeName + "Shield";             // the "sides" of the cone
  std::string frontPlateConeShieldTransName = frontPlateConeShieldName + "Trans";

  new TGeoConeSeg(frontPlateConeShieldName.c_str(),
                  sDzAluCone / 2,
                  sDrMinAluFront - thicknessFrontPlateCone,
                  sDrMinAluFront,
                  sDrMinAluCone,
                  sDrMinAluCone + thicknessFrontPlateCone,
                  -90,
                  90);
  TGeoTranslation* frontPlateConeShieldTrans = new TGeoTranslation(frontPlateConeShieldTransName.c_str(),
                                                                   0,
                                                                   0,
                                                                   zPosCone);
  frontPlateConeShieldTrans->RegisterYourself();

  // Frontplate cone "bottom"
  float zPosConePlate = -sDzAluFront / 2 + sDzAluCone - thicknessFrontPlateCone / 2;   // z-position of the cone bottom relative to the frontplate
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
  TGeoTranslation* frontPlateConePlateCsTrans = new TGeoTranslation(frontPlateConePlateCSTransName.c_str(),
                                                                    0,
                                                                    0,
                                                                    zPosConePlate);
  frontPlateConePlateCsTrans->RegisterYourself();

  // Frontplate cone composite shape
  std::string frontPlateConeCSBoolFormula;
  frontPlateConeCSBoolFormula += frontPlateConeShieldName + ":" + frontPlateConeShieldTransName;
  frontPlateConeCSBoolFormula += "+" + frontPlateConePlateCSName + ":" + frontPlateConePlateCSTransName;

  std::string frontPlateConeCSName = frontPlateConeName + "CompositeShape";
  new TGeoCompositeShape(frontPlateConeCSName.c_str(), frontPlateConeCSBoolFormula.c_str());

  // Shields
  float dzShieldGap = 0.7;                      // z-distance between the shields and the front- and backplate outer edges (in z-direction)
  float dzShield = sDzAlu - 2 * dzShieldGap;    // depth of the shields
  
  // Outer shield
  float zPosOuterShield = zPosBackPlate - (zPosBackPlate - zPosFrontPlate) / 2;   // z-position of the outer shield
  
  std::string outerShieldName = "FV0_OuterShield";
  std::string outerShieldTransName = outerShieldName + "Trans";
  
  new TGeoTubeSeg(outerShieldName.c_str(), sDrMinAluOuterShield, sDrMaxAluOuterShield, dzShield / 2, -90, 90);
  TGeoTranslation* shieldTrans = new TGeoTranslation(outerShieldTransName.c_str(), 0, 0, zPosOuterShield);
  shieldTrans->RegisterYourself();

  // Inner shield
  float dzInnerShield = sDzAlu - sDzAluCone - dzShieldGap;                                  // depth of the inner shield
  float zPosInnerShield = zPosBackPlate + sDzAluBack / 2 - dzShieldGap - dzInnerShield / 2; // z-position of the inner shield relative to the backplate
  
  std::string innerShieldName = "FV0_InnerShield";
  std::string innerShieldCutName = innerShieldName + "Cut";
  
  new TGeoTubeSeg(innerShieldName.c_str(), sDrMinAluInnerShield, sDrMaxAluInnerShield, dzInnerShield / 2, -90, 90);
  new TGeoBBox(innerShieldCutName.c_str(), fabs(sXAluHole), sDrMaxAluInnerShield, dzInnerShield / 2);

  // Inner shield composite shape
  std::string innerShieldCSBoolFormula;
  innerShieldCSBoolFormula = innerShieldName;
  innerShieldCSBoolFormula += "-" + innerShieldCutName;

  std::string innerShieldCSName = innerShieldName + "CS";
  std::string innerShieldCSTransName = innerShieldCSName + "Trans";
  new TGeoCompositeShape(innerShieldCSName.c_str(), innerShieldCSBoolFormula.c_str());
  TGeoTranslation* innerShieldCSTrans = new TGeoTranslation(innerShieldCSTransName.c_str(),
                                                            sXAluHole,
                                                            0,
                                                            zPosInnerShield);
  innerShieldCSTrans->RegisterYourself();

  // Cover
  float dzCover = sDzAlu;                             // Depth of the covers
  float zPosCoverConeCut = zPosFrontPlate + zPosCone; // Set the cone cut relative to the frontplate so that the exact position of the aluminium cone part can be used.
  
  std::string coverName = "FV0_Cover";
  std::string coverConeCutName = coverName + "ConeCut";
  std::string coverHoleCutName = coverName + "HoleCut";
 
  new TGeoBBox(coverName.c_str(), sDxAluCover / 2, sDrMaxAluOuterShield, dzCover / 2);
  new TGeoConeSeg(coverConeCutName.c_str(),
                  sDzAluCone / 2,
                  0,
                  sDrMinAluFront,
                  0,
                  sDrMinAluCone + thicknessFrontPlateCone,
                  -90,
                  90);
  new TGeoTubeSeg(coverHoleCutName.c_str(), 0, sDrMinAluInnerShield, dzCover / 2, 0, 360);  

  std::string coverTransName = coverName + "Trans";
  std::string coverConeCutTransName = coverConeCutName + "Trans";
  std::string coverHoleCutTransName = coverHoleCutName + "Trans";

  TGeoTranslation* coverTrans = new TGeoTranslation(coverTransName.c_str(), sDxAluCover / 2, 0, zPosOuterShield);
  TGeoTranslation* coverConeCutTrans = new TGeoTranslation(coverConeCutTransName.c_str(),
                                                           0,
                                                           0,
                                                           zPosCoverConeCut);
  TGeoTranslation* coverHoleCutTrans = new TGeoTranslation(coverHoleCutTransName.c_str(),
                                                           sXAluHole,
                                                           0,
                                                           zPosOuterShield);

  coverTrans->RegisterYourself();
  coverConeCutTrans->RegisterYourself();
  coverHoleCutTrans->RegisterYourself();

  // Cover composite shape
  std::string coverCSBoolFormula = "";
  coverCSBoolFormula += coverName + ":" + coverTransName;
  coverCSBoolFormula += "-" + coverConeCutName + ":" + coverConeCutTransName;
  coverCSBoolFormula += "-" + coverHoleCutName + ":" + coverHoleCutTransName;

  std::string coverCSName = coverName + "CS";
  new TGeoCompositeShape(coverCSName.c_str(), coverCSBoolFormula.c_str());

  // Composite shape
  std::string boolFormula = "";
  boolFormula += backPlateCSName + ":" + backPlateCSTransName;
  boolFormula += "+" + frontPlateName + ":" + frontPlateTransName;
  boolFormula += "+" + frontPlateConeCSName + ":" + frontPlateTransName;
  boolFormula += "+" + outerShieldName + ":" + outerShieldTransName;
  boolFormula += "+" + innerShieldCSName + ":" + innerShieldCSTransName;
  boolFormula += "+" + coverCSName;

  std::string aluContCSName = "FV0_AluContCS";
  TGeoCompositeShape* aluContCS = new TGeoCompositeShape(aluContCSName.c_str(), boolFormula.c_str());

  // Volume
  std::string aluContName = "FV0" + sContainerName;
  TGeoMedium* kMed = gGeoManager->GetMedium("V0_Aluminium$");
  new TGeoVolume(aluContName.c_str(), aluContCS, kMed);
}

void Geometry::initializeLuts()
{
  // TODO: initialize sth
}

void Geometry::buildGeometry()
{
  TGeoVolume* vALIC = gGeoManager->GetVolume("cave");
  if (!vALIC) {
    LOG(FATAL) << "Could not find the top volume";
  }

  // Top volume of FIT V0 detector
  TGeoVolumeAssembly* vFV0 = new TGeoVolumeAssembly("FITV0");
  LOG(INFO) << "FV0 Geometry::buildGeometry()::Volume name = " << vFV0->GetName();

  assembleSensVols(vFV0);
  if (mGeometryType == eFull) {
    assembleNonSensVols(vFV0);
  }

  TGeoTranslation* trGlobalZshift = new TGeoTranslation(0, 0, sZposition);

  vALIC->AddNode(vFV0, 0, trGlobalZshift);
}

void Geometry::assembleSensVols(TGeoVolumeAssembly* vFV0)
{
  assembleScintSectors(vFV0);
}

void Geometry::assembleNonSensVols(TGeoVolumeAssembly* vFV0)
{
  assemblePlasticSectors(vFV0);
  assembleFibers(vFV0);
  assembleMetalContainer(vFV0);
}

void Geometry::assembleScintSectors(TGeoVolumeAssembly* volV0)
{
  TGeoVolumeAssembly* v0scint = new TGeoVolumeAssembly("FV0SCINT");
  TGeoVolumeAssembly* v0ScintLeft = new TGeoVolumeAssembly("FV0SCINTLEFT");
  TGeoVolumeAssembly* v0ScintRight = new TGeoVolumeAssembly("FV0SCINTRIGHT");

  for (uint16_t isector = 0; isector < ceil(mvSectorTrans.size() / 2); isector++) {
    TGeoVolumeAssembly* sector = buildSector(sScintCellName, isector);
    v0ScintLeft->AddNode(sector, isector + 1, mvSectorTrans.at(isector));
  }
  
  for (uint16_t isector = ceil(mvSectorTrans.size() / 2); isector < mvSectorTrans.size(); isector++) {
    TGeoVolumeAssembly* sector = buildSector(sScintCellName, isector);
    v0ScintRight->AddNode(sector, isector + 1, mvSectorTrans.at(isector));
  }
  
  v0scint->AddNode(v0ScintLeft, 1);
  v0scint->AddNode(v0ScintRight, 1);
  volV0->AddNode(v0scint, 1);
}

void Geometry::assemblePlasticSectors(TGeoVolumeAssembly* volV0) {
  TGeoVolumeAssembly* v0Plast = new TGeoVolumeAssembly("FV0PLAST");
  TGeoVolumeAssembly* v0PlastLeft = new TGeoVolumeAssembly("FV0PLASTLEFT");
  TGeoVolumeAssembly* v0PlastRight = new TGeoVolumeAssembly("FV0PLASTRIGHT");

  for (uint16_t isector = 0; isector < ceil(mvSectorTrans.size() / 2); isector++) {
    TGeoVolumeAssembly* sector = buildSector(sPlastCellName, isector);
    v0PlastLeft->AddNode(sector, isector + 1, mvSectorTrans.at(isector));
  }
  
  for (uint16_t isector = ceil(mvSectorTrans.size() / 2); isector < mvSectorTrans.size(); isector++) {
    TGeoVolumeAssembly* sector = buildSector(sPlastCellName, isector);
    v0PlastRight->AddNode(sector, isector + 1, mvSectorTrans.at(isector));
  }
  
  v0Plast->AddNode(v0PlastLeft, 1);
  v0Plast->AddNode(v0PlastRight, 1);
  volV0->AddNode(v0Plast, 1, new TGeoTranslation(0, 0, sZPlast));
}

void Geometry::assembleFibers(TGeoVolumeAssembly* vFV0)
{
  TGeoVolumeAssembly* fibers = new TGeoVolumeAssembly("FV0FIBERS");

  for (int i = 0; i < mvrMinFiber.size(); i++) {
    std::stringstream ssFiberName;
    ssFiberName << "FV0" << sFiberName << i + 1;
    TGeoVolume* fiber = gGeoManager->GetVolume(ssFiberName.str().c_str());
    fibers->AddNode(fiber, 1);
  }

  TGeoRotation* reflection = new TGeoRotation();
  reflection->ReflectX(true);

  vFV0->AddNode(fibers, 1);
  vFV0->AddNode(fibers, 2, reflection);
}

void Geometry::assembleMetalContainer(TGeoVolumeAssembly* volV0) {
  std::string containerName = "FV0" + sContainerName;
  TGeoVolume* container = gGeoManager->GetVolume(containerName.c_str());
  if (!container) {
    LOG(WARNING) << "FV0: Couldn't find volume " << containerName;
  } else {
    LOG(DEBUG) << "FVO Geometry::assembleMetalContainer(): adding container volume " << containerName;
    TGeoRotation* reflection = new TGeoRotation();
    reflection->ReflectX(true);
    volV0->AddNode(container, 1);
    volV0->AddNode(container, 2, reflection);
  }
}

TGeoVolumeAssembly* Geometry::buildSector(std::string cellType, int iSector)
{
  std::stringstream ssSectorName;
  ssSectorName << "FV0" << cellType << iSector + 1;

  LOG(DEBUG) << "FV0 Geometry::buildSector(): building sector " << ssSectorName.str();

  TGeoVolumeAssembly* sector = new TGeoVolumeAssembly(ssSectorName.str().c_str());
  
  for (int i = 0; i < sNumberOfRings; i++) {
    std::stringstream cellName;
    cellName << "FV0" << cellType << sCellTypes[iSector] << i + 1;

    TGeoVolume* cell = gGeoManager->GetVolume(cellName.str().c_str());

    if (!cell) {
      LOG(WARNING) << "FV0: Couldn't find cell volume " << cellName.str();
    } else {
      LOG(DEBUG) << "FV0 Geometry::buildSector(): adding cell volume " << cellName.str();
      sector->AddNode(cell, i + 1);
    }    
  }

  return sector;
}