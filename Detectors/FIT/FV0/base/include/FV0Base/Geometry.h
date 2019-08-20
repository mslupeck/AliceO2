// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Geometry.h
/// \brief Base definition of FV0+ geometry.
///
/// \author Maciej Slupecki, University of Jyvaskyla, Finland

#ifndef ALICEO2_FV0_GEOMETRY_H_
#define ALICEO2_FV0_GEOMETRY_H_

#include <vector>
#include <TGeoMatrix.h>
#include <TGeoVolume.h>
#include "TVirtualMC.h"

namespace o2
{
namespace fv0
{
/// FV0+ Geometry
class Geometry
{
 public:
  enum EGeoType {
    eUninitilized,
    eDummy,
    eOnlySensitive,
    eFull
  }; // Geometry type options possible to be initialized

  enum EGeoComponent {
    eScint,
    ePlast,
    eFiber,
    eScrew,
    eRod,
    eAluminum
  };

  ///
  /// Default constructor.
  /// It must be kept public for root persistency purposes,
  /// but should never be called by the outside world
  Geometry() { mGeometryType = eUninitilized; };
  /// Standard constructor
  /// \param initType[in]  The type of geometry, that will be initialized
  ///                       -> initType == 0 => only sensitive detector parts
  ///                       -> initType == 1 => sensitive parts and rough structural elements
  ///                       -> initType == 2 => complete, detailed geometry (including screws, etc.)
  /// \return  -
  Geometry(EGeoType initType);
  /// Copy constructor.
  Geometry(const Geometry& geom);

  static constexpr float sEpsilon = 0.01;                   // variable used to make sure one spatial dimension is infinitesimaly larger than the other
  static constexpr float sDzScint = 4;                      // thickness of scintillator
  static constexpr float sDzPlast = 1;                      // thickness of fiber plastic
  static constexpr float sXGlobal = 0;                      // global x-position of the geometrical center of the sensitive parts of the detector
  static constexpr float sYGlobal = 0;                      // global y-position of the geometrical center of the sensitive parts of the detector
  // TODO: Adjust the sZposition once the simulation geometry is implemented, T0 starts at 328
  // at sZposition==320, there is a gap (to be filled with fibers and support) of 8 cm between the plastic of FV0+ and aluminum covers of FT0+
  static constexpr float sZGlobal = 320 - sDzScint / 2;     // global z-position of the geometrical center of the sensitive parts of the detector
  static constexpr float sGlobalPhiRotation = 0;            // global phi rotation (to enable making two detector halves, possible to separate vertically)
  static constexpr float sDxHalvesSeparation = 0;           // separation between the left and right side of the detector
  static constexpr float sDyHalvesSeparation = 0;           // y-position of the right detector part relative to the left
  static constexpr float sDzHalvesSeparation = 0;           // z-position of the right detector part relative to the left
  static constexpr float sDrSeparationScint = 0.03 + 0.04;  // paint thickness + half of separation gap
  static constexpr float sDxCellHoleExt = 0.2;              // extension of the holes closest to the vertical aluminum cover
  static constexpr float sDySeparationScint = sDrSeparationScint;
  static constexpr int sBaseNumberOfSectors = 8; // number of sectors
  static constexpr int sNumberOfRings = 5;                                                               // number of rings
  static constexpr float sRingRadiiScint[sNumberOfRings + 1] = { 4.01, 7.3, 12.9, 21.25, 38.7, 72.115 }; // average ring radii
  static constexpr float sRingInnerRadiusDx = -0.15;                                                     // shift of the inner radius origin
  static constexpr char sCellTypes[sBaseNumberOfSectors] = { 'a', 'b', 'b', 'a', 'a', 'b', 'b', 'a'};

  /// Get the unique ID of the current scintillator cell during simulation.
  /// The ID is a number from 1 to 40 starting from the first cell left of the y-axis
  /// and continues counterclockwise one ring at a time.
  /// \param  fMC  The virtual Monte Carlo interface.
  const int getCurrentCellId(TVirtualMC* fMC);

  /// Get the names of all the sensitive volumes of the geometry.
  const std::vector<std::string> getSensitiveVolumeNames() { return mvSensitiveVolumeNames; };

  /// Enable or disable a geometry component. To be called before the geometry is built. A disabled component will not
  /// be added to the geometry. The enabled components are by default specified by the geometry type.
  /// \param  component The geometry component to be enabled/disabled.
  /// \param  enable  Boolean setting the enabled state. Default is true.
  /// \return The enabled state of the geometry component.
  bool enableComponent(EGeoComponent component, bool enable = true);

  /// Build the geometry.
  void buildGeometry();

 private:
  static constexpr float sDrHoleLarge = 0.415;              // radius of the large scintillator hole
  static constexpr float sDrHoleSmall = 0.265;              // radius of the small scintillator hole
  // Aluminium container constants
  static constexpr float sDzAlu = 30;                       // depth of aluminium container
  static constexpr float sDrAluHole = 4.05;                 // radius of beam hole
  static constexpr float sXAluHole = -0.15;                 // shift of beam hole
  static constexpr float sDrMaxAluBack = 83.1;              // outer radius of aluminium backplate
  static constexpr float sDzAluBack = 1;                    // thickness of aluminium backplate
  static constexpr float sDrMinAluFront = 45.7;             // inner radius of aluminium frontplate
  static constexpr float sDrMaxAluFront = 83.1;             // outer radius of aluminium frontplate
  static constexpr float sDzAluFront = 1;                   // thickness of aluminium frontplate
  static constexpr float sDxAluStand = 40;                  // the width of the aluminium stand
  static constexpr float sDyAluStand = 3;                   // the height of the aluminium stand at x = 0
  static constexpr float sDrMinAluCone = 24.3;              // inner radius at the bottom of aluminium cone
  static constexpr float sDzAluCone = 16.2;                 // depth of alminium frontplate cone
  static constexpr float sThicknessAluCone = 0.6;           // thickness of aluminium frontplate cone
  static constexpr float sXYThicknessAluCone = 0.975;       // radial thickness in the xy-plane of the aluminium cone.
  static constexpr float sDrMinAluOuterShield = 82.5;       // inner radius of outer aluminium shield
  static constexpr float sDrMaxAluOuterShield = 82.65;      // outer radius of outer aluminium shield
  static constexpr float sDrMinAluInnerShield = 4;          // inner radius of inner aluminium shield
  static constexpr float sDrMaxAluInnerShield = 4.05;       // outer radius of inner aluminium shield
  static constexpr float sDxAluCover = 0.15;                // thickness of aluminium cover

  // Local position constants
  static constexpr float sZScint = 0;                                                         // scintillator z-position
  static constexpr float sZPlast = sZScint + sDzScint / 2 + sDzPlast / 2;                     // plastic z-position
  static constexpr float sZAluBack = sZScint - sDzScint / 2 - sDzAluBack / 2;                 // aluminium backplate z-position
  static constexpr float sZAluFront = sZAluBack - sDzAluBack / 2 + sDzAlu - sDzAluFront / 2;  // aluminium frontplate z-position
  static constexpr float sZAluMid = (sZAluBack + sZAluFront) / 2;                             // middle z of aluminum container
  static constexpr float sZFiber = (sZPlast + sZAluFront) / 2;                                // fiber z-position (plastic and frontplate midpoint)
  static constexpr float sZCone = sZAluFront + sDzAluFront / 2 - sDzAluCone / 2;              // aluminium frontplate cone z-position
  static constexpr float sXShiftScrews = sDxAluCover;                                         // x shift of all screw holes

  /// Initialize the geometry.
  void initializeGeometry();

  /// Initialize maps with geometry information.
  void initializeMaps();

  /// Initialize vectors with geometry information.
  void initializeVectors();

  /// Initialize the cell radii.
  void initializeCellRadii();

  /// Initialize sector transformations.
  void initializeSectorTransformations();

  /// Initialize fiber volume radii.
  void initializeFiberRadii();

  /// Initialize the screw type dimensions;
  void initializeScrewTypeDimensions();

  /// Initialize the rod type dimenstions;
  void initializeRodTypeDimensions();

  /// Initialize the position and dimension for every screw and rod.
  void initializeScrewAndRodPositionsAndDimensions();

  /// Initialize the sensitive volumes.
  void initializeSensVols();

  /// Initialize the non-sensitive volumes.
  void initializeNonSensVols();

  /// Initialize a composite shape of all screw holes as well as a translation for this. This shape is removed from all
  /// volumes that the screws are passing through to avoid overlaps.
  void initializeScrewHoles();

  /// Initialize a composite shape of all rod holes as well as a translation for this. This shape is removed from all
  /// volumes that the rods are passing through to avoid overlaps.
  void initializeRodHoles();

  /// Initialize cell volumes with a specified thickness and medium.
  /// \param  cellType  The type of the cells.
  /// \param  zThicknes The thickness of the cells.
  /// \param  medium    The medium of the cells.
  void initializeCells(std::string cellType, const float zThickness, TGeoMedium* medium);

  /// Initialize scintillator cell volumes.
  void initializeScintCells();

  /// Initialize plastic cell volumes for fiber support.
  void initializePlasticCells();

  /// Initialize volumes equivalent to the fibers.
  void initializeFibers();

  /// Initialize the screw volumes.
  void initializeScrews();

  /// Initialize the rod volumes.
  void initializeRods();

  /// Initialize the metal container volume.
  void initializeMetalContainer();

  /// Assemble the sensitive volumes.
  /// \param  vFV0  The FIT V0 volume.
  void assembleSensVols(TGeoVolumeAssembly* vFV0);

  /// Assemble the non sensitive volumes.
  /// \param  vFV0  The FIT V0 volume.
  void assembleNonSensVols(TGeoVolumeAssembly* vFV0);

  /// Assemble the scintillator sectors.
  /// \param  vFV0  The FIT V0 volume.
  void assembleScintSectors(TGeoVolumeAssembly* vFV0);

  /// Assemble the plastice sectors.
  /// \param  vFV0  The FIT V0 volume.
  void assemblePlasticSectors(TGeoVolumeAssembly* vFV0);

  /// Assemble the fibers.
  /// \param  vFV0  The FIT V0 volume.
  void assembleFibers(TGeoVolumeAssembly* vFV0);

  /// Assemble the screwss.
  /// \param  vFV0  The FIT V0 volume.
  void assembleScrews(TGeoVolume* vFV0);

  /// Assemble the rods.
  /// \param  vFV0  The FIT V0 volume.
  void assembleRods(TGeoVolume* vFV0);

  /// Assemble the metal container.
  /// \param  vFV0  The FIT V0 volume.
  void assembleMetalContainer(TGeoVolumeAssembly* vFV0);

  /// Build sector assembly of specified type.
  /// \param  cellName  The type of the cells in the sector assembly.
  /// \return The sector assembly.
  TGeoVolumeAssembly* buildSectorAssembly(std::string cellName);

  /// Build a sector of specified type and number.
  /// \param  cellType  The type of the cells in the sector.
  /// \param  iSector   The numbering of the sector.
  /// \return The sector.
  TGeoVolumeAssembly* buildSector(std::string cellType, int iSector);

  /// Helper function for creating and registering a TGeoTranslation.
  TGeoTranslation* createAndRegisterTrans(std::string name, double dx, double dy, double dz);

  /// Helper function for creating and registering a TGeoRotation.
  TGeoRotation* createAndRegisterRot(std::string name, double phi, double theta, double psi);

  inline static const std::string sScintName = "SCINT";
  inline static const std::string sPlastName = "PLAST";
  inline static const std::string sSectorName = "SECTOR";
  inline static const std::string sCellName = "CELL";
  inline static const std::string sScintSectorName = sScintName + sSectorName;
  inline static const std::string sScintCellName = sScintName + sCellName;
  inline static const std::string sPlastSectorName = sPlastName + sSectorName;
  inline static const std::string sPlastCellName = sPlastName + sCellName;
  inline static const std::string sFiberName = "FIBER";
  inline static const std::string sScrewName = "SCREW";
  inline static const std::string sScrewHolesCSName = "FV0SCREWHOLES";
  inline static const std::string sScrewHolesCSTransName = sScrewHolesCSName + "Trans";
  inline static const std::string sRodName = "ROD";
  inline static const std::string sRodHolesCSName = "FV0RODHOLES";
  inline static const std::string sRodHolesCSTransName = sRodHolesCSName + "Trans";
  inline static const std::string sContainerName = "ALUCONTAINER";

  std::vector<std::string> mvSensitiveVolumeNames;

  std::vector<float> mRAvgRing;                     // average ring radii (index 0 -> ring 1 min, index 1 -> ring 1 max and ring 2 min, ... index 5 -> ring 5 max)
  // The following radii include separation between scintillator rings
  std::vector<float> mRMinScint;                    // inner radii of a ring (.at(0) -> ring 1, .at(4) -> ring 5)
  std::vector<float> mRMaxScint;                    // outer radii of a ring (.at(0) -> ring 1, .at(4) -> ring 5)
  std::vector<float> mRMinFiber;                    // inner radii of fiber volumes (.at(0) -> fiber 1)
  std::vector<float> mRMaxFiber;                    // outer radii of fiber volumes (.at(0) -> fiber 1)

  std::vector<float> mDzScrews;                     // length of screws (.at(n) -> length of screw no. n)
  std::vector<float> mDzScrewTypes;                 // the different length of the screws
  std::vector<float> mDrScrews;                     // radii of the screws (.at(n) -> radius of screw no. n)
  std::vector<float> mDrScrewTypes;                 // the different radii of the screws
  std::vector<float> mRScrews;                      // radii for the screw locations
  std::vector<int> mScrewTypes;                     // the type no. of each screw (.at(n) -> type no. of screw no. n)

  std::vector<float> mDzRods;                       // length of rods (.at(n) -> length of rod no. n)
  std::vector<float> mDzRodTypes;                   // the different length of the rods
  std::vector<float> mDrRods;                       // radii of the rods (.at(n) -> radius of rod no. n)
  std::vector<float> mDrRodTypes;                   // the different radii of the rods
  std::vector<float> mDxRods;                       // width of the rods (.at(n) -> width of rod no. n)
  std::vector<float> mDxRodTypes;                   // the different width of the rods
  std::vector<int> mRodTypes;                       // the type no. of each rod (.at(n) -> type no. of rod no. n)

  std::vector<TGeoMatrix*> mSectorTrans;            // transformations of sectors (.at(0) -> sector 1)
  std::vector<std::vector<float>> mScrewPos;        // xy-coordinates of all the screws
  std::vector<std::vector<float>> mRodPos;          // xy-coordinates of all the rods

  int mGeometryType;                                // same meaning as initType in constructor
  std::map<EGeoComponent, bool> mEnabledComponents; // map of the enabled state of all geometry components.

  ClassDefNV(Geometry, 1);
};
} // namespace fv0
} // namespace o2
#endif
