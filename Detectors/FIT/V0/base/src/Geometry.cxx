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

void Geometry::initializeVectors()
{
  // RADII

  // Scintillator radii
  LOG(INFO) << "Geometry::initializeVectors(): Initializing V0 scintillator ring radii.";

  // Index of mvrAvgScint is NOT linked directly to any ring number
  mvrAvgScint.assign(sRingRadiiScint, sRingRadiiScint + sNumberOfRings + 1);

  // Set real scintillator radii (reduced by paint thickness and separation gap)
  for (uint16_t ir = 1; ir < mvrAvgScint.size(); ir++) { // shift of indices to match index with ring number (starting from 0)
    mvrMaxScint.push_back(mvrAvgScint.at(ir) - sDrSeparationScint);
  }
  for (uint16_t ir = 0; ir < mvrAvgScint.size() - 1; ir++) {
    mvrMinScint.push_back(mvrAvgScint.at(ir) + sDrSeparationScint);
  }
  // Now indices of rMinScint and rMaxScint correspond to the same ring

  // SECTOR TRANSFORMATIONS
  LOG(INFO) << "Geometry::initializeVectors(): Initializing V0 sector transformations.";

  for (uint16_t isector = 0; isector < sBaseNumberOfSectors; isector++) {
    // isector = 0 corresponds to the first sector counter-clockwise from the y-axis + global azimuthal rotation

    float phi = sPhiMinScint + (isector + 1) * sDphiScint;

    std::stringstream ssRotName, ssReflName, ssTotalTransName;
    ssRotName << "rotationSector" << isector + 1;
    ssReflName << "reflectSector" << isector + 1;
    ssTotalTransName << "totalTransformSector" << isector + 1;

    TGeoRotation* rot = new TGeoRotation(ssRotName.str().c_str());            // rotation of the sector
    TGeoRotation* refl = new TGeoRotation(ssReflName.str().c_str());          // reflection of the sector
    TGeoHMatrix* totTrans = new TGeoHMatrix(ssTotalTransName.str().c_str());  // the combined transformation

    switch (isector)
    {
    case 0:
      // "a"-mirror
      refl->ReflectX(true);
      break;
    case 1:
      // "b"-mirror
      refl->ReflectX(true);
      break;
    case 2:
      // "b"
      rot->SetAngles(phi + sDphiScint, 0., 0.);
      break;
    case 3:
      // "a"
      rot->SetAngles(phi, 0., 0.);
      break;
    case 4:
      // "a"-mirror
      refl->ReflectY(true);
      break;
    case 5:
      // "b"-mirror
      refl->ReflectY(true);
      break;
    case 6:
      // "b"
      rot->SetAngles(phi + sDphiScint, 0., 0.);
      break;
    case 7:
      // "a"
      rot->SetAngles(phi, 0., 0.);
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

  LOG(INFO) << "Geometry::initializeScintCells(): Initializing scintillator cells.";
  
  TGeoMedium* kMed = gGeoManager->GetMedium("V0_Scintillator$");
  
  float dz = sDzScint / 2;                                      // half depth of the scintillator

  // Sector separation gap
  std::string secSepShapeName = "V0ScintSeparation";
  TGeoBBox* secSepShape = new TGeoBBox(secSepShapeName.c_str(), mvrMaxScint.back(), sDrSeparationScint, dz);

  // Sector separation gap rotations
  std::string secSepRot45Name, secSepRot90Name;
  secSepRot45Name = "V0secSepRot45";
  secSepRot90Name = "V0secSepRot90";
  TGeoRotation* rot45 = new TGeoRotation(secSepRot45Name.c_str(), 45, 0, 0);
  rot45->RegisterYourself();
  TGeoRotation* rot90 = new TGeoRotation(secSepRot90Name.c_str(), 90, 0, 0);
  rot90->RegisterYourself();
  
  for (int ir = 0; ir < sNumberOfRings; ir++) {
    
    // "a"-type cell
    // 
    // Initial placement:
    //
    // y
    // ^
    // |  * * * * *
    // |  * a * *
    // |  * * * 
    // |  * *
    // |
    // 0--------------> x

    std::stringstream aCellName;
    aCellName << "V0a" << ir + 1;

    LOG(INFO) << "Geometry::initializeScintCells(): Initializing cell " << aCellName.str();

    // Base shape
    std::string aCellShapeName = aCellName.str() + "shape";

    // The cells in the innermost ring are slightly different than the rest.
    // The inner radius is slightly shifted.
    if (ir == 0) {
      // The innermost "a"-type cell

      std::string a1CellShapeFullName = aCellShapeName + "Full";
      std::string a1CellShapeCutName = aCellShapeName + "Cut";
      std::string a1CellShapeCutTransName = a1CellShapeCutName + "Trans";
      
      TGeoTubeSeg* a1CellShapeFull = new TGeoTubeSeg(a1CellShapeFullName.c_str(), 0, mvrMaxScint.at(ir), dz, 45, 90);
      TGeoTubeSeg* a1CellShapeCut = new TGeoTubeSeg(a1CellShapeCutName.c_str(), 0, mvrMinScint.at(ir), dz, 0, 360);
      TGeoTranslation* a1CellShapeCutTrans = new TGeoTranslation(a1CellShapeCutTransName.c_str(),
                                                                 sRingInnerRadiusDx, 0, 0);
      a1CellShapeCutTrans->RegisterYourself();

      std::string a1BoolFormula = a1CellShapeFullName + "-" + a1CellShapeCutName + ":" + a1CellShapeCutTransName;
      TGeoCompositeShape* a1CellShape = new TGeoCompositeShape(aCellShapeName.c_str(), a1BoolFormula.c_str());
    } else {
      // The rest of the "a"-type cells
      TGeoTubeSeg* aCellShape = new TGeoTubeSeg(aCellShapeName.c_str(), mvrMinScint.at(ir), mvrMaxScint.at(ir), dz,
                                                45, 90);
    }

    // Composite shape
    std::string aBoolFormula = aCellShapeName
                               + "-" + secSepShapeName + ":" + secSepRot45Name
                               + "-" + secSepShapeName + ":" + secSepRot90Name;

    TGeoCompositeShape* aCellCs = new TGeoCompositeShape("aCellCs", aBoolFormula.c_str());

    // Cell volume
    TGeoVolume* aCell = new TGeoVolume(aCellName.str().c_str(), aCellCs, kMed);
    aCell->RegisterYourself();

    // "b"-type cells
    // 
    // Initial placement:
    //
    // y
    // ^
    // |          *
    // |        * *
    // |      * * *
    // |    * * b *
    // |      * * *
    // |
    // 0--------------> x

    std::stringstream bCellName;
    bCellName << "V0b" << ir + 1;

    LOG(INFO) << "Geometry::initializeScintCells(): Initializing cell " << bCellName.str();

    // Base shape
    std::string bCellShapeName = bCellName.str() + "shape";

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

    // Composite shape
    std::string bBoolFormula = bCellShapeName
                               + "-" + secSepShapeName
                               + "-" + secSepShapeName + ":" + secSepRot45Name;

    TGeoCompositeShape* bCellCs = new TGeoCompositeShape("bCellCs", bBoolFormula.c_str());

    // Cell volume
    TGeoVolume* bCell = new TGeoVolume(bCellName.str().c_str(), bCellCs, kMed);
    bCell->RegisterYourself();
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
  LOG(INFO) << "Geometry::buildGeometry()::Volume name = " << volV0->GetName();

  assembleScintSectors(volV0);

  TGeoTranslation* trGlobalZshift = new TGeoTranslation(0, 0, sZposition);

  vALIC->AddNode(volV0, 0, trGlobalZshift);
}

void Geometry::assembleScintSectors(TGeoVolumeAssembly* volV0)
{
  TGeoVolumeAssembly* v0scint = new TGeoVolumeAssembly("FITV0_scint");
  TGeoVolumeAssembly* v0ScintLeft = new TGeoVolumeAssembly("FITV0_scint_left");
  TGeoVolumeAssembly* v0ScintRight = new TGeoVolumeAssembly("FITV0_scint_right");

  int nRot = mvSectorTrans.size();

  for (uint16_t isector = 0; isector < ceil(nRot / 2); isector++) {
    TGeoVolumeAssembly* sector = buildScintSector(isector);
    v0ScintLeft->AddNode(sector, isector + 1, mvSectorTrans.at(isector));
  }
  
  for (uint16_t isector = ceil(nRot / 2); isector < nRot; isector++) {
    TGeoVolumeAssembly* sector = buildScintSector(isector);
    v0ScintRight->AddNode(sector, isector + 1, mvSectorTrans.at(isector));
  }
  
  v0scint->AddNode(v0ScintLeft, 1);
  v0scint->AddNode(v0ScintRight, 2);
  volV0->AddNode(v0scint, 1);
}

TGeoVolumeAssembly* Geometry::buildScintSector(uint16_t iSector)
{
  std::stringstream sectorName;
  sectorName << "sector" << iSector + 1;

  LOG(INFO) << "Geometry::buildScintSector(): building sector " << sectorName.str();

  TGeoVolumeAssembly* sector = new TGeoVolumeAssembly(sectorName.str().c_str());
  
  for (int i = 0; i < sNumberOfRings; i++) {
    std::stringstream cellName;
    cellName << "V0" << sCellTypes[iSector] << i + 1;

    TGeoVolume* cell = gGeoManager->GetVolume(cellName.str().c_str());

    LOG(INFO) << "Geometry::buildSector(): adding cell " << cellName.str();
    sector->AddNode(cell, i + 1);
  }

  return sector;
}
