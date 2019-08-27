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
/// \author Andreas Molander, University of Helsinki, Finland

#ifndef ALICEO2_FV0_GEOMETRY_H_
#define ALICEO2_FV0_GEOMETRY_H_

#include <vector>
#include <TGeoMatrix.h>
#include <TGeoVolume.h>
#include <TVirtualMC.h>

namespace o2
{
namespace fv0
{
/// FV0+ Geometry
class Geometry
{
 public:
  /// Geometry type options possible to be initialized
  enum EGeoType {
    eUninitilized,
    eDummy,
    eOnlySensitive,
    eRough,
    eFull
  };

  /// Geometry components possible to be enabled/disabled
  enum EGeoComponent {
    eScintillator,
    ePlastics,
    eFibers,
    eScrews,
    eRods,
    eAluminiumContainer
  };

  /// Default constructor.
  /// It must be kept public for root persistency purposes,
  /// but should never be called by the outside world
  Geometry() { mGeometryType = eUninitilized; };

  // TODO: What does eDummy geometry type initialize?

  /// Standard constructor
  /// \param initType[in]  The type of geometry, that will be initialized
  ///                       -> initType == eUnitialized   => no parts
  ///                       -> initType == eDummy         => ?
  ///                       -> initType == eOnlySensitive => only sensitive detector parts
  ///                       -> initType == eRough         => sensitive parts and rough structural elements
  ///                       -> initType == eFull          => complete, detailed geometry (including screws, etc.)
  /// \return  -
  Geometry(EGeoType initType);

  /// Copy constructor.
  Geometry(const Geometry& geometry);

  /// Get the unique ID of the current scintillator cell during simulation.
  /// The ID is a number from 1 to 40 starting from the first cell left of the y-axis
  /// and continues counterclockwise one ring at a time.
  /// \param  fMC  The virtual Monte Carlo interface.
  /// \return The ID of the current scintillator cell during simulation.
  const int getCurrentCellId(TVirtualMC* fMC);

  /// Get the names of all the sensitive volumes of the geometry.
  const std::vector<std::string> getSensitiveVolumeNames() { return mSensitiveVolumeNames; };

  /// Enable or disable a geometry component. To be called before the geometry is built. A disabled component will not
  /// be added to the geometry. The enabled components are by default specified by the geometry type.
  /// \param  component The geometry component to be enabled/disabled.
  /// \param  enable  Boolean setting the enabled state. Default is true.
  /// \return The enabled state of the geometry component.
  bool enableComponent(EGeoComponent component, bool enable = true);

  /// Build the geometry.
  void buildGeometry();

 private:
  // General geometry constants
  static constexpr float sEpsilon = 0.01;                   // variable used to make sure one spatial dimension is infinitesimaly larger than the other
  static constexpr float sDzScint = 4;                      // thickness of scintillator
  static constexpr float sDzPlast = 1;                      // thickness of fiber plastic
  static constexpr float sXGlobal = 0;                      // global x-position of the geometrical center of the sensitive parts of the detector
  static constexpr float sYGlobal = 0;                      // global y-position of the geometrical center of the sensitive parts of the detector
  // TODO: Adjust the sZposition once the simulation geometry is implemented, T0 starts at 328
  // at sZGlobal==320, there is a gap (to be filled with fibers and support) of 8 cm between the plastic of FV0+ and aluminum covers of FT0+
  static constexpr float sZGlobal = 320 - sDzScint / 2;     // global z-position of the geometrical center of the sensitive parts of the detector
  static constexpr float sGlobalPhiRotation = 0;            // global phi rotation (to enable making two detector halves, possible to separate vertically)
  static constexpr float sDxHalvesSeparation = 0;           // separation between the left and right side of the detector
  static constexpr float sDyHalvesSeparation = 0;           // y-position of the right detector part relative to sYGlobal
  static constexpr float sDzHalvesSeparation = 0;           // z-position of the right detector part relative to sZGlobal

  // Cell and scintillator constants
  // static constexpr float sDySeparationScint = sDrSeparationScint;
  static constexpr int sNCellSectors = 8;                   // number of sectors
  static constexpr int sNCellRings = 5;                     // number of rings
  static constexpr float sCellRingRadii[sNCellRings + 1] { 4.01, 7.3, 12.9, 21.25, 38.7, 72.115 }; // average ring radii
  static constexpr char sCellTypes[sNCellSectors] { 'a', 'b', 'b', 'a', 'a', 'b', 'b', 'a'};  // ordered cell types per ring in
  static constexpr float sDrSeparationScint = 0.03 + 0.04;  // paint thickness + half of separation gap
  
  static constexpr float sXShiftInnerRadiusScint = -0.15;   // shift of the inner radius origin of the scintillators
  static constexpr float sDxHoleExtensionScint = 0.2;       // extension of the scintillator holes closest to the vertical aluminum cover
  static constexpr float sDrHoleLargeScint = 0.415;         // radius of the large scintillator hole
  static constexpr float sDrHoleSmallScint = 0.265;         // radius of the small scintillator hole

  // Container constants
  static constexpr float sDzAlu = 30;                       // depth of aluminium container
  static constexpr float sDrAluHole = 4.05;                 // radius of beam hole
  static constexpr float sXShiftAluHole = -0.15;            // shift of beam hole
  static constexpr float sDrMaxAluBack = 83.1;              // outer radius of aluminium backplate
  static constexpr float sDzAluBack = 1;                    // thickness of aluminium backplate
  static constexpr float sDrMinAluFront = 45.7;             // inner radius of aluminium frontplate
  static constexpr float sDrMaxAluFront = 83.1;             // outer radius of aluminium frontplate
  static constexpr float sDzAluFront = 1;                   // thickness of aluminium frontplate
  static constexpr float sDxAluStand = 40;                  // the width of the aluminium stand
  static constexpr float sDyAluStand = 3;                   // the height of the aluminium stand at local x == 0
  static constexpr float sDrMinAluCone = 24.3;              // inner radius at the bottom of aluminium cone
  static constexpr float sDzAluCone = 16.2;                 // depth of alminium frontplate cone
  static constexpr float sThicknessAluCone = 0.6;           // thickness of aluminium frontplate cone
  static constexpr float sXYThicknessAluCone = 0.975;       // radial thickness in the xy-plane of the aluminium cone.
  static constexpr float sDrMinAluOuterShield = 82.5;       // inner radius of outer aluminium shield
  static constexpr float sDrMaxAluOuterShield = 82.65;      // outer radius of outer aluminium shield
  static constexpr float sDrMinAluInnerShield = 4;          // inner radius of inner aluminium shield
  static constexpr float sDrMaxAluInnerShield = 4.05;       // outer radius of inner aluminium shield
  static constexpr float sDxAluCover = 0.15;                // thickness of aluminium cover
  static constexpr float sDxAluStandBottom = 38.5;          // width of aluminum stand bottom
  static constexpr float sDyAluStandBottom = 2;             // thickness of aluminum stand bottom

  // Local position constants
  static constexpr float sPosScint[] { sDxAluCover, 0, 0 };                                   // local position of the right half of the scintillator
  static constexpr float sZPlast = sPosScint[2] + sDzScint / 2 + sDzPlast / 2;                // plastic z-position
  static constexpr float sZAluBack = sPosScint[2] - sDzScint / 2 - sDzAluBack / 2;            // aluminium backplate z-position
  static constexpr float sZAluFront = sZAluBack - sDzAluBack / 2 + sDzAlu - sDzAluFront / 2;  // aluminium frontplate z-position
  static constexpr float sZAluMid = (sZAluBack + sZAluFront) / 2;                             // middle z of aluminum container
  static constexpr float sZFiber = (sZPlast + sZAluFront) / 2;                                // fiber z-position (plastic and frontplate midpoint)
  static constexpr float sZCone = sZAluFront + sDzAluFront / 2 - sDzAluCone / 2;              // aluminium frontplate cone z-position
  static constexpr float sXShiftScrews = sPosScint[0];                                         // x shift of all screw holes

  // Screw and rod dimensions
  static constexpr int sNScrewTypes = 6;                                                      // number of different screw types
  static constexpr float sDrMinScrewTypes[sNScrewTypes] { 0.25, 0.25, 0.4, 0.4, 0.4, 0.4 };   // radius of the thinner part of the screw types
  static constexpr float sDrMaxScrewTypes[sNScrewTypes] { 0, 0.5, 0.6, 0.6, 0.6, 0};          // radius of the thicker part of the screw types
  static constexpr float sDzMaxScrewTypes[sNScrewTypes] { 6.02, 13.09, 13.1, 23.1, 28.3, 5 }; // length of the thinner part of the screw types
  static constexpr float sDzMinScrewTypes[sNScrewTypes] { 0, 6.78, 6.58, 15.98, 21.48, 0 };   // length of the thicker part of the screw types
  static constexpr float sZShiftScrew = 0;                                                    // z shift of screws. 0 means they are aligned with the scintillator.

  static constexpr int sNRodTypes = 4;
  static constexpr float sDxMinRodTypes[sNRodTypes] { 0.366, 0.344, 0.344, 0.344 };           // width of the thinner part of the rod types
  static constexpr float sDxMaxRodTypes[sNRodTypes] { 0.536, 0.566, 0.566, 0.716 };           // width of the thicker part of the rod types
  static constexpr float sDyMinRodTypes[sNRodTypes] { 0.5, 0.8, 0.8, 0.8 };                   // height of the thinner part of the rod types
  static constexpr float sDyMaxRodTypes[sNRodTypes] { 0.9, 1.2, 1.2, 1.2 };                   // height of the thicker part of the rod types
  static constexpr float sDzMaxRodTypes[sNRodTypes] { 12.5, 12.5, 22.5, 27.7 };               // length of the thinner part of the rod types
  static constexpr float sDzMinRodTypes[sNRodTypes] { 7.45, 7.45, 17.45, 22.65 };             // length of the thicker part of the rod types
  static constexpr float sZShiftRod = -0.05;                                                  // z shift of rods. 0 means they are aligned with tht scintillators.

  // Strings for volume names etc.
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
  inline static const std::string sRodName = "ROD";
  inline static const std::string sRodHolesCSName = "FV0RODHOLES";
  inline static const std::string sContainerName = "ALUCONTAINER";

  /// Initialize the geometry.
  void initializeGeometry();

  /// Initialize maps with geometry information.
  void initializeMaps();

  /// Initialize vectors with geometry information.
  void initializeVectors();

  /// Initialize the cell ring radii.
  void initializeCellRingRadii();

  /// Initialize sector transformations.
  void initializeSectorTransformations();

  /// Initialize fiber volume radii.
  void initializeFiberVolumeRadii();

  /// Initialize fiber mediums.
  void initializeFiberMedium();

  /// Initialize the radii of the screw and rod positions.
  void initializeScrewAndRodRadii();

  /// Initialize the screw type medium.
  void initializeScrewTypeMedium();

  /// Initialize the rod type medium.
  void initializeRodTypeMedium();

  /// Add a screw propery set to the collection of total screws.
  /// \param  screwTypeID The screw type ID.
  /// \param  iRing       The ring number.
  /// \param  phi         Azimuthal angle of the screw location.
  void addScrewProperties(int screwTypeID, int iRing, float phi);

  /// Add a rod property set to the collectino of total rods.
  /// \param  rodTypeID The rod type ID.
  /// \param  iRing     The ring number.
  void addRodProperties(int rodTypeID, int iRing);

  /// Initialize the position and dimension for every screw and rod.
  void initializeScrewAndRodPositionsAndDimensions();

  /// Initialize the sensitive volumes.
  void initializeSensVols();

  /// Initialize the non-sensitive volumes.
  void initializeNonSensVols();

  /// Initialize a composite shape of all screw holes. This shape is removed from all volumes that the screws are
  /// passing through to avoid overlaps.
  void initializeScrewHoles();

  /// Initialize a composite shape of all rod holes. This shape is removed from all volumes that the rods are passing
  /// through to avoid overlaps.
  void initializeRodHoles();

  /// Initialize cell volumes with a specified thickness and medium.
  /// \param  cellType    The type of the cells.
  /// \param  zThicknes   The thickness of the cells.
  /// \param  medium      The medium of the cells.
  /// \param  isSensitive Specifies if the cells are sensitive volumes.
  void initializeCells(std::string cellType, const float zThickness, TGeoMedium* medium, bool isSensitive);

  /// Initialize scintillator cell volumes.
  void initializeScintCells();

  /// Initialize plastic cell volumes for optical fiber support.
  void initializePlasticCells();

  /// Initialize volumes equivalent to the optical fibers.
  void initializeFibers();

  /// Initialize the screw volumes.
  void initializeScrews();

  /// Initialize the rod volumes.
  void initializeRods();

  /// Initialize the metal container volume.
  void initializeMetalContainer();

  /// Assemble the sensitive volumes.
  /// \param  vFV0  The FIT V0 volume.
  void assembleSensVols(TGeoVolume* vFV0);

  /// Assemble the non sensitive volumes.
  /// \param  vFV0  The FIT V0 volume.
  void assembleNonSensVols(TGeoVolume* vFV0);

  /// Assemble the scintillator sectors.
  /// \param  vFV0  The FIT V0 volume.
  void assembleScintSectors(TGeoVolume* vFV0);

  /// Assemble the plastice sectors.
  /// \param  vFV0  The FIT V0 volume.
  void assemblePlasticSectors(TGeoVolume* vFV0);

  /// Assemble the optical fibers.
  /// \param  vFV0  The FIT V0 volume.
  void assembleFibers(TGeoVolume* vFV0);

  /// Assemble the screwss.
  /// \param  vFV0  The FIT V0 volume.
  void assembleScrews(TGeoVolume* vFV0);

  /// Assemble the rods.
  /// \param  vFV0  The FIT V0 volume.
  void assembleRods(TGeoVolume* vFV0);

  /// Assemble the metal container.
  /// \param  vFV0  The FIT V0 volume.
  void assembleMetalContainer(TGeoVolume* vFV0);

  /// Build sector assembly of specified type.
  /// \param  cellName  The type of the cells in the sector assembly.
  /// \return The sector assembly.
  TGeoVolumeAssembly* buildSectorAssembly(std::string cellName);

  /// Build a sector of specified type and number.
  /// \param  cellType  The type of the cells in the sector.
  /// \param  iSector   The numbering of the sector.
  /// \return The sector.
  TGeoVolumeAssembly* buildSector(std::string cellType, int iSector);
  
  /// Create the shape for a specified screw.
  /// \param  shapeName   The name of the shape.
  /// \param  screwTypeID The number of the screw type.
  /// \param  xEpsilon    Shrinks or expands the x dimensions of the screw shape.
  /// \param  yEpsilon    Shrinks or expands the y dimensions of the screw shape.
  /// \param  zEpsilon    Shrinks or expands the z dimensions of the screw shape.
  /// \return The screw shape.
  TGeoShape* createScrewShape(std::string shapeName, int screwTypeID, float xEpsilon = 0, float yEpsilon = 0, float zEpsilon = 0);

  /// Create the shape for a specified rod.
  /// \param  shapeName The name of the shape.
  /// \param  rodTypeID The number of the rod type.
  /// \param  xEpsilon  Shrinks or expands the x dimensions of the rod shape.
  /// \param  yEpsilon  Shrinks or expands the y dimensions of the rod shape.
  /// \param  zEpsilon  Shrinks or expands the z dimensions of the rod shape.
  /// \return The rod shape.
  TGeoShape* createRodShape(std::string shapeName, int rodTypeID, float xEpsilon = 0, float yEpsilon = 0, float zEpsilon = 0);

  /// Helper function for creating and registering a TGeoTranslation.
  TGeoTranslation* createAndRegisterTrans(std::string name, double dx, double dy, double dz);

  /// Helper function for creating and registering a TGeoRotation.
  TGeoRotation* createAndRegisterRot(std::string name, double phi, double theta, double psi);

  std::vector<std::string> mSensitiveVolumeNames;   // The names of all the sensitive volumes

  std::vector<float> mRAvgRing;                     // average ring radii (index 0 -> ring 1 min, index 1 -> ring 1 max and ring 2 min, ... index 5 -> ring 5 max)
  std::vector<float> mRMinScint;                    // inner radii of scintillator rings (.at(0) -> ring 1, .at(4) -> ring 5)
  std::vector<float> mRMaxScint;                    // outer radii of scintillator rings (.at(0) -> ring 1, .at(4) -> ring 5)
  std::vector<float> mRMinFiber;                    // inner radii of fiber volumes (.at(0) -> fiber 1)
  std::vector<float> mRMaxFiber;                    // outer radii of fiber volumes (.at(0) -> fiber 1)
  std::vector<TGeoMedium*> mMediumFiber;            // Medium of the fiber volumes (.at(n) -> medium of fiber the n:th fiber starting from the middle)

  std::vector<float> mRScrewAndRod;                 // radii of the screw and rod positions.

  std::vector<float> mDrMinScrews;                  // radii of the thinner part of the screws
  std::vector<float> mDrMaxScrews;                  // radii of the thicker part of the screws
  std::vector<float> mDzMaxScrews;                  // length of the thinner part of the screws
  std::vector<float> mDzMinScrews;                  // length of the thicker part of the screws
  
  std::vector<float> mRScrews;                      // radial distance to the screw locations
  std::vector<int> mScrewTypeIDs;                   // the type ID of each screw (.at(n) -> type ID of screw no. n)

  std::vector<float> mDxMinRods;                    // width of the thinner part of the rods
  std::vector<float> mDxMaxRods;                    // width of the thicker part of the rods
  std::vector<float> mDyMinRods;                    // height of the thinner part of the rods
  std::vector<float> mDyMaxRods;                    // height of the thicker part of the rods
  std::vector<float> mDzMaxRods;                    // length of the thinner part of the rods
  std::vector<float> mDzMinRods;                    // length of the thicker part of the rods
  std::vector<int> mRodTypeIDs;                     // the type ID of each rod (.at(n) -> type ID of rod no. n)

  std::vector<TGeoMatrix*> mSectorTrans;            // transformations of sectors (.at(0) -> sector 1)
  std::vector<std::vector<float>> mScrewPos;        // xyz-coordinates of all the screws
  std::vector<std::vector<float>> mRodPos;          // xyz-coordinates of all the rods
  std::vector<TGeoMedium*> mMediumScrewTypes;       // medium of the screw types
  std::vector<TGeoMedium*> mMediumRodTypes;         // medium of the rod types

  int mGeometryType;                                // same meaning as initType in constructor
  std::map<EGeoComponent, bool> mEnabledComponents; // map of the enabled state of all geometry components.

  ClassDefNV(Geometry, 1);
};
} // namespace fv0
} // namespace o2
#endif
