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
  initializeVectors();
  initializeScintCells();
  initializeLuts();
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
}

void Geometry::initializeScintCells()
{
  // Creating the two types of cells, "a" and "b", for each ring.
  // All sectors can be assembled with these cells.
  //
  // The reference to "a" and "b" can be understood with the CAD drawings of the detector.

  LOG(INFO) << "FV0 Geometry::initializeScintCells(): Initializing scintillator cells.";
  
  TGeoMedium* kMed = gGeoManager->GetMedium("V0_Scintillator$");
  
  float dz = sDzScint / 2;                            // half depth of the scintillator
  float rHoleSmall = 0.265;                           // radius of the small holes
  float rHoleLarge = 0.415;                           // radius of the large holes
  float xHoleCut = 0.2;                               // width of extension of hole 1, 2 and 7 in the "a" cell
  float dxHole = sDrSeparationScint + xHoleCut;       // x-placement of holes 1, 2 and 7 in the "a" cell

  // Sector separation gap shape
  std::string secSepShapeName = "FV0_scintSectorSeparation";
  new TGeoBBox(secSepShapeName.c_str(), mvrMaxScint.back(), sDrSeparationScint, dz);

  // Sector separation gap rotations
  std::string secSepRot45Name = "FV0_secSepRot45";
  std::string secSepRot90Name = "FV0_secSepRot90";

  TGeoRotation* secSepRot45 = new TGeoRotation(secSepRot45Name.c_str(), 45, 0, 0);  
  TGeoRotation* secSepRot90 = new TGeoRotation(secSepRot90Name.c_str(), 90, 0, 0);

  secSepRot45->RegisterYourself();
  secSepRot90->RegisterYourself();

  // Hole shapes
  std::string holeSmallName = "FV0_scintHoleSmall";
  std::string holeLargeName = "FV0_scintHoleLarge";
  std::string holeSmallCutName = "FV0_scintHoleSmallCut";
  std::string holeLargeCutName = "FV0_scintHoleLargeCut";

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
    aCellName << "FV0" << sScintCellName << "a" << ir + 1;

    LOG(INFO) << "FV0 Geometry::initializeScintCells(): Initializing cell " << aCellName.str();

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

    LOG(INFO) << "FV0 Geometry::initializeScintCells(): Composite shape formula: " << aBoolFormula;

    std::string aCellCSName = aCellName.str() + "CS";
    TGeoCompositeShape* aCellCs = new TGeoCompositeShape(aCellCSName.c_str(), aBoolFormula.c_str());

    // Cell volume
    TGeoVolume* aCell = new TGeoVolume(aCellName.str().c_str(), aCellCs, kMed);
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
    bCellName << "FV0" << sScintCellName << "b" << ir + 1;

    LOG(INFO) << "FV0 Geometry::initializeScintCells(): Initializing cell " << bCellName.str();

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
    TGeoVolume* bCell = new TGeoVolume(bCellName.str().c_str(), bCellCs, kMed);
    bCell->RegisterYourself();
    mvSensitiveVolumeNames.push_back(bCell->GetName());
  }
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
  TGeoVolumeAssembly* volV0 = new TGeoVolumeAssembly("FITV0");
  LOG(INFO) << "FV0 Geometry::buildGeometry()::Volume name = " << volV0->GetName();

  assembleScintSectors(volV0);

  TGeoTranslation* trGlobalZshift = new TGeoTranslation(0, 0, sZposition);

  vALIC->AddNode(volV0, 0, trGlobalZshift);
}

void Geometry::assembleScintSectors(TGeoVolumeAssembly* volV0)
{
  TGeoVolumeAssembly* v0scint = new TGeoVolumeAssembly("FV0SCINT");
  TGeoVolumeAssembly* v0ScintLeft = new TGeoVolumeAssembly("FV0SCINTLEFT");
  TGeoVolumeAssembly* v0ScintRight = new TGeoVolumeAssembly("FV0SCINTRIGHT");

  for (uint16_t isector = 0; isector < ceil(mvSectorTrans.size() / 2); isector++) {
    TGeoVolumeAssembly* sector = buildScintSector(isector);
    v0ScintLeft->AddNode(sector, isector + 1, mvSectorTrans.at(isector));
  }
  
  for (uint16_t isector = ceil(mvSectorTrans.size() / 2); isector < mvSectorTrans.size(); isector++) {
    TGeoVolumeAssembly* sector = buildScintSector(isector);
    v0ScintRight->AddNode(sector, isector + 1, mvSectorTrans.at(isector));
  }
  
  v0scint->AddNode(v0ScintLeft, 1);
  v0scint->AddNode(v0ScintRight, 2);
  volV0->AddNode(v0scint, 1);
}

TGeoVolumeAssembly* Geometry::buildScintSector(uint16_t iSector)
{
  std::stringstream ssSectorName;
  ssSectorName << "FV0" << sScintSectorName << iSector + 1;

  LOG(INFO) << "FV0 Geometry::buildScintSector(): building sector " << ssSectorName.str();

  TGeoVolumeAssembly* sector = new TGeoVolumeAssembly(ssSectorName.str().c_str());
  
  for (int i = 0; i < sNumberOfRings; i++) {
    std::stringstream cellName;
    cellName << "FV0" << sScintCellName << sCellTypes[iSector] << i + 1;

    TGeoVolume* cell = gGeoManager->GetVolume(cellName.str().c_str());

    if (!cell) {
      LOG(WARNING) << "FV0: Couldn't find scintillator cell " << cellName.str();
    } else {
      LOG(INFO) << "FV0 Geometry::buildScintSector(): adding cell " << cellName.str();
      sector->AddNode(cell, i + 1);
    }    
  }

  return sector;
}
