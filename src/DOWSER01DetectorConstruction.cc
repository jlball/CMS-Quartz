//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
/// \file DOWSER01DetectorConstruction.cc
/// \brief Implementation of the DOWSER01DetectorConstruction class

#include "DOWSER01DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4Para.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4SDManager.hh"
#include "G4SDChargedFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSEnergyDeposit3D.hh"
#include "G4PSTrackLength.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4GenericMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01DetectorConstruction::DOWSER01DetectorConstruction()
: G4VUserDetectorConstruction(),
fCheckOverlaps(true)
{ ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01DetectorConstruction::~DOWSER01DetectorConstruction()
{ ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DOWSER01DetectorConstruction::Construct()
{
  //=====================
  // Material Definitions
  //=====================
  //
  //-------- NIST Materials ----------------------------------------------------
  //  Material Information imported from NIST database.
  //
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;
  // nistManager->FindOrBuildMaterial("G4_Pb", fromIsotopes);
  nistManager->FindOrBuildMaterial("G4_Galactic", fromIsotopes);
  G4Material* vacuum  = G4Material::GetMaterial("G4_Galactic");
  
  // nistManager->FindOrBuildMaterial("G4_B", fromIsotopes);
  // G4Material* boron10 = G4Material::GetMaterial("G4_B");
  // density = 1.0 *g/cm3;
  // a = 10.01294*g/mole;
  // G4Material* boron10 = new G4Material("boron10", 5., 10.01294*g/mole, 1.0*g/cm3);
  // nistManager->FindOrBuildMaterial("boron10", fromIsotopes);
  
  G4int protons = 5, neutrons = 5, nucleons=protons+neutrons, isotopes, ncomponents;
  G4double atomicMass = 10.01294*g/mole;
  G4Isotope* isoB10 = new G4Isotope("isoB10", protons, nucleons, atomicMass);
  G4Element* elmB10 = new G4Element("elmB10", "B10", isotopes = 1);
  elmB10->AddIsotope(isoB10,100*perCent);
  G4Material* boron10 = new G4Material("boron10", 0.94770*g/cm3, ncomponents=1, kStateSolid);
  boron10->AddElement(elmB10,100*perCent);
  
  nistManager->FindOrBuildMaterial("G4_POLYETHYLENE", fromIsotopes);
  G4Material* polyethylene = G4Material::GetMaterial("G4_POLYETHYLENE");
  
  // Liquid argon material
  /*
   G4double a;  // mass of a mole;
   G4double z;  // z=mean number of protons;
   G4double density;
   new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
   // The argon by NIST Manager is a gas with a different density
   */
  
  /*
   // Vacuum
   new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
   kStateGas, 2.73*kelvin, 3.e-18*pascal);
   */
  
  // Helium 3
  // new version full Helium 3, 2013-02-02; J.J. Su
  protons=2;
  neutrons=1;
  nucleons=protons+neutrons;
  atomicMass = 3.016*g/mole;
  G4Isotope* he3 = new G4Isotope("He3", protons, nucleons, atomicMass);
  G4Element* He3 = new G4Element("Helium3", "He3", isotopes=1);
  He3->AddIsotope(he3, 100*perCent);
  G4double pressure = 2.5*atmosphere;
  G4double temperature = 273*kelvin;
  G4double molar_constant = Avogadro*k_Boltzmann;  //from clhep
  G4double he3density = (atomicMass*pressure)/(temperature*molar_constant);
  G4Material* Helium3 = new G4Material("Helium3", he3density, ncomponents=1, kStateGas, temperature, pressure);
  Helium3->AddElement(He3, 100*perCent);
  
  // Cd for SETN container, Cadmium
  G4Material* cadmium = nistManager->FindOrBuildMaterial("G4_Cd", fromIsotopes);
  
  // building Xenon
  protons = 54;
  atomicMass = 131.29*g/mole;
  // density = 5.458*mg/cm3;
  pressure = 2*atmosphere;
  G4double Xedensity = (atomicMass*pressure)/(temperature*molar_constant);
  G4Material* xenon = new G4Material("xenon", protons, atomicMass, Xedensity, kStateGas, temperature ,pressure);
  
  //         FindOrBuildMaterial (const G4String &name, G4bool isotopes=true, G4bool warning=false)
  // BuildMaterialWithNewDensity (const G4String &name, G4double pres=CLHEP::STP_Pressure)
  // Ti fuel tank
  // G4Material* fuelTankMaterial = nistManager->FindOrBuildMaterial("G4_Ti", fromIsotopes);

  // Al for substract, Aluminum
  G4Material* aluminum = nistManager->FindOrBuildMaterial("G4_Al", fromIsotopes);
  
  // Xe
  // G4Material* xenon = nistManager->FindOrBuildMaterial("G4_Xe", fromIsotopes);

  
  G4cout << "****************************************************  "  << G4endl;
  G4cout << "************* Temperautre " << temperature << G4endl;
  G4cout << "************* Avogadro constant " << Avogadro << G4endl;
  G4cout << "************* Boltzmann constant " << k_Boltzmann << G4endl;
  G4cout << "************* Pressure " << pressure/atmosphere << G4endl;
  G4cout << "************* Density of Helium " << he3density << G4endl;
  G4cout << "************* Density of Helium " << he3density/(g/cm3) << " g/cm3" << G4endl;
  G4cout << "****************************************************  "  << G4endl;
  
  
  // BC545 Materials
  // Natural Boron-loaded Premium Plastic Scintillator
  // 5% Boron
  // composition, # of atoms-> C: 4.43, H: 5.18, 10B: 0.0559
  //
  G4String symbol;
  G4double a, z, density;
  G4Element* C = new G4Element("Carbon", symbol="C", z=6., a= 12.01*g/mole);
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* B10 = new G4Element("Boron10",symbol="B10", z= 5., a =10.01294*g/mole);
  G4Material* BC545 = new G4Material("BC545", density=1.026*g/cm3, ncomponents=3, kStateSolid);
  G4double fractionmass;
  BC545->AddElement(C, fractionmass = 0.9018);
  BC545->AddElement(H, fractionmass = 0.0887);
  BC545->AddElement(B10, fractionmass = 0.0095);
  
  nistManager->FindOrBuildMaterial("G4_BGO", fromIsotopes);
  G4Material* BGO  = G4Material::GetMaterial("G4_BGO");
  
  // Boron Cabride
  // B4C
  G4int natoms;
  G4Material* B4C = new G4Material("B4C",density = 2.52*g/cm3, ncomponents = 2, kStateSolid);
  B4C->AddElement(C, natoms = 1);
  B4C->AddElement(B10, natoms = 4);
  
  
  
  
  // Print materials
  G4cout << G4endl << "The materials defined are: " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  G4cout << "************  end of material table  ************" << G4endl;
  
  //============================================================================
  //      Definitions of Solids, Logical Volumes, Physical Volumes
  //============================================================================
  
  //-------------
  // World Volume
  //-------------
  G4RotationMatrix* rot = new G4RotationMatrix();
  rot->rotateZ(0.*deg);
  
  G4ThreeVector worldSize = G4ThreeVector(50*cm, 50*cm, 100*cm);
  G4Box * solidWorld
  = new G4Box("soildWorld", worldSize.x()/2., worldSize.y()/2., worldSize.z()/2.);
  G4LogicalVolume * World
  = new G4LogicalVolume(solidWorld, vacuum, "World", 0, 0, 0);
  
  //
  //  Must place the World Physical volume unrotated at (0,0,0).
  G4VPhysicalVolume * worldPV
  = new G4PVPlacement(0,               // no rotation
                      G4ThreeVector(), // at (0,0,0)
                      World,      // its logical volume
                      "WorldPV",         // its name
                      0,               // its mother  volume
                      false,           // no boolean operations
                      0);              // copy number


  // Define dimensions of each components
  G4double rmin = 0*mm, rmax = 12.70*mm; // 0.5*mm; //
  G4double thicknessB10 = 0.001*mm; // 0.250*mm; // 0.5*mm;
  G4double thicknessAl = 1.0*mm; // 0.250*mm; // 0.5*mm;
  G4double thicknessXe = 4.826*mm;
  G4double thicknessPoly = 50.8*mm; // 5.0*mm; //
  G4double thicknessCd = 1.0*mm; // 1*mm; //
  G4double slitWidth = 5.0*mm; //
  G4double diaAlCasing = 40*mm;
  G4double thickAlCasing = 0.127*mm;
  G4double depthAlCasing = 20*mm;
  G4double gapB10film = 4.826*mm;
  
  // Building geometries of each componets
  // G4Tubs * xenonZone = new G4Tubs("xenonZone", rmin, rmax, thicknessXe*0.5, 0.*deg, 360.*deg);
  // G4Tubs* boronFilm  = new G4Tubs("boronFilm", rmin, rmax, thicknessB10*0.5, 0.*deg, 360.*deg);
  // G4Tubs* alSubstrate  = new G4Tubs("alSubstrate", rmin, rmax, thicknessAl*0.5, 0.*deg, 360.*deg);
  // G4Tubs* cdSlit  = new G4Tubs("cdSlit", rmin, rmax, thicknessCd*0.5, 0.*deg, 360.*deg);
  G4Box * boronFilm = new G4Box ("boronFilm", rmax, rmax, 0.5*thicknessB10);
  G4Box * alSubstrate = new G4Box ("alSubstrate", rmax, rmax, 0.5*thicknessAl);
  G4Box * xenonZone = new G4Box ("xenonZone", rmax, rmax, 0.5*thicknessXe);
  G4Tubs * cdSlit  = new G4Tubs("cdSlit", rmin, 0.5*(diaAlCasing + thickAlCasing), thicknessCd*0.5, 0.*deg, 360.*deg);
  G4Tubs* hdpeBlock  = new G4Tubs("hdpeBlock", rmin, 0.5*(diaAlCasing + thickAlCasing), thicknessPoly*0.5, 0.*deg, 360.*deg);
  G4Box * slitWindow = new G4Box ("slitWindow", 0.5*slitWidth, 0.4*rmax, 1.0*cm);
  
  // Building aluminum casing
  G4Tubs* tube1 = new G4Tubs("tube1", 0, 0.5*diaAlCasing, 0.5*depthAlCasing, 0.*deg, 360.*deg);
  G4Tubs* tube2 = new G4Tubs("tube2", 0, 0.5*diaAlCasing + thickAlCasing, 0.5*depthAlCasing + thickAlCasing, 0.*deg, 360.*deg);

  G4VSolid* alCasing
  = new G4SubtractionSolid ("alCasing", tube2, tube1, 0, G4ThreeVector(0,0,0));
  
  // Building attachement
  // for HDPE
  G4Tubs* tube3a = new G4Tubs("tube3a", 0, 1.25*diaAlCasing + thickAlCasing, 0.5*thicknessPoly, 0.*deg, 360.*deg);
  // Tapping a tube with 30 degree slop
  G4double displacement = 1.25*diaAlCasing + thickAlCasing;
  G4Para * tube3b = new G4Para ("tube3b", displacement, displacement, displacement, 0*deg, 65*deg, 0*deg);
  G4VSolid* tube3
  = new G4SubtractionSolid ("tube3", tube3a, tube3b, 0, G4ThreeVector(0,0,0.5*thicknessPoly));
  
  
  // for Cd slit
  G4Tubs* tube4 = new G4Tubs("tube4", 0, 1.25*diaAlCasing + thickAlCasing, 0.5*thicknessCd, 0.*deg, 360.*deg);
  G4Box* slitGap = new G4Box("slitGsp", 0.5*slitWidth, 0.8*rmax, 10*cm);
  G4VSolid* slit1
  = new G4SubtractionSolid ("slit1", tube4, slitGap, 0, G4ThreeVector(0,0,0));

  


  

  rot->rotateY(0.*deg);
  //.....*******oooooOOO000OOOooooo*******.....
  // WEDGE HALF ANGLE
  //.....*******oooooOOO000OOOooooo*******.....
  // G4double rotatingAngle = 30.0;
  // G4RotationMatrix* rot1 = new G4RotationMatrix();
  // rot1->rotateY(rotatingAngle*deg);
  // G4RotationMatrix* rot2 = new G4RotationMatrix();
  // rot2->rotateY(-rotatingAngle*deg);
  // G4RotationMatrix* rot0 = new G4RotationMatrix();
  // rot0->rotateY(rotatingAngle*deg);
  
  
  // .....*******oooooOOO000OOOooooo*******.....
  // Placing Al substrate
  
  displacement = -0.5*(gapB10film+thicknessB10);

  G4LogicalVolume * boronFilmLV = new G4LogicalVolume(boronFilm, boron10, "boronFilmLV", 0, 0, 0);
  G4ThreeVector position01 = G4ThreeVector(0.0, 0.0, displacement);
  G4VPhysicalVolume * boronFilm1
  = new G4PVPlacement(0,               // no rotation
                    position01,    // shift He3Detector toward the beam source
                    boronFilmLV,       // its logical volume
                    "borobFilm1",      // its name
                    World,         // its mother  volume
                    false,           // no boolean operations
                    0,              // copy number
                    fCheckOverlaps);

  displacement = 0.5*(gapB10film+thicknessB10);
  position01 = G4ThreeVector(0.0, 0.0, displacement);
  G4VPhysicalVolume * boronFilm2
  = new G4PVPlacement(0,               // no rotation
                    position01,    // shift He3Detector toward the beam source
                    boronFilmLV,       // its logical volume
                    "borobFilm2",      // its name
                    World,         // its mother  volume
                    false,           // no boolean operations
                    0,              // copy number
                    fCheckOverlaps);

  // Placing B10 film on Al substrate
  displacement = -0.5*(gapB10film+thicknessAl)+thicknessB10;
  G4LogicalVolume * alSubstrateLV = new G4LogicalVolume(alSubstrate, aluminum, "alSubstrateLV", 0, 0, 0);
  position01 = G4ThreeVector(0.0, 0.0, displacement);
  G4VPhysicalVolume * alSubstrate1
  = new G4PVPlacement(0,               // no rotation
                    position01,    // shift He3Detector toward the beam source
                    alSubstrateLV,       // its logical volume
                    "alSubstrate1",      // its name
                    World,         // its mother  volume
                    false,           // no boolean operations
                    0,              // copy number
                    fCheckOverlaps);

  displacement = 0.5*(gapB10film+thicknessAl)+thicknessB10;
  position01 = G4ThreeVector(0.0, 0.0, displacement);
  G4VPhysicalVolume * alSubstrate2
  = new G4PVPlacement(0,               // no rotation
                    position01,    // shift He3Detector toward the beam source
                    alSubstrateLV,       // its logical volume
                    "alSubstrate2",      // its name
                    World,         // its mother  volume
                    false,           // no boolean operations
                    0,              // copy number
                    fCheckOverlaps);


  // Placing xenon gas
  
  G4LogicalVolume * xenonLV = new G4LogicalVolume(xenonZone, aluminum, "xenonLV", 0, 0, 0);
  position01 = G4ThreeVector(0.0, 0.0, 0.0);
  G4VPhysicalVolume * xenonPV
  = new G4PVPlacement(0,               // no rotation
                    position01,    // shift He3Detector toward the beam source
                    xenonLV,       // its logical volume
                    "xenonPV",      // its name
                    World,         // its mother  volume
                    false,           // no boolean operations
                    0,              // copy number
                    fCheckOverlaps);

  
  // Placing aluminum casing
  
  displacement =  -thicknessXe*0.5 ;
  G4LogicalVolume * casingLV = new G4LogicalVolume(alCasing, aluminum, "xenonLV", 0, 0, 0);
  position01 = G4ThreeVector(0.0, 0.0, displacement);
  G4VPhysicalVolume * casingPV
  = new G4PVPlacement(0,               // no rotation
                    position01,    // shift He3Detector toward the beam source
                    casingLV,       // its logical volume
                    "casingPV",      // its name
                    World,         // its mother  volume
                    false,           // no boolean operations
                    0,              // copy number
                    fCheckOverlaps);


  // Placing cd slit on top of the casing
  
  displacement =  0.5*depthAlCasing - 0.5*thicknessXe + thickAlCasing + 0.5*thicknessCd;
  /*
  G4LogicalVolume * cdSlitLV = new G4LogicalVolume(slit1, cadmium, "cdSlitLV", 0, 0, 0);
  position01 = G4ThreeVector(0.0, 0.0, displacement);
  G4VPhysicalVolume * cdSlitPV1
  = new G4PVPlacement(0,               // no rotation
                    position01,    // shift He3Detector toward the beam source
                    cdSlitLV,       // its logical volume
                    "cdSlitPV1",      // its name
                    World,         // its mother  volume
                    false,           // no boolean operations
                    0,              // copy number
                    fCheckOverlaps);
  */

  // Placing cd slit on top of the casing
  
  displacement = 0.5*depthAlCasing - 0.5*thicknessXe + thickAlCasing + 0.5*(thicknessPoly);
  G4LogicalVolume * hdpeLV = new G4LogicalVolume(tube3, polyethylene, "hdpeLV", 0, 0, 0);
  position01 = G4ThreeVector(0.0, 0.0, displacement);
  G4VPhysicalVolume * hdpePV
  = new G4PVPlacement(0,               // no rotation
                    position01,    // shift He3Detector toward the beam source
                    hdpeLV,       // its logical volume
                    "hdpePV",      // its name
                    World,         // its mother  volume
                    false,           // no boolean operations
                    0,              // copy number
                    fCheckOverlaps);
  
  // Placing 2nd cd slit on top of the casing
  /*
  displacement += 0.5*(thicknessCd + thicknessPoly);
  position01 = G4ThreeVector(0.0, 0.0, displacement);
  G4VPhysicalVolume * cdSlitPV2
  = new G4PVPlacement(0,               // no rotation
                    position01,    // shift He3Detector toward the beam source
                    cdSlitLV,       // its logical volume
                    "cdSlitPV2",      // its name
                    World,         // its mother  volume
                    false,           // no boolean operations
                    0,              // copy number
                    fCheckOverlaps);
  */


  
  

  //
  // Visualization attributes
  //
  // WorldLV->SetVisAttributes (G4VisAttributes::Invisible);
  G4VisAttributes* whiteBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  // G4VisAttributes* redBoxVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  G4VisAttributes* greenBoxVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0,0.25));
  G4VisAttributes* redBoxVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.75));
  G4VisAttributes* blueBoxVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.25));
  G4VisAttributes* yellowBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.25));
  G4VisAttributes* lightblueVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,0.5,0.15));
  G4VisAttributes* grayVisAtt= new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.1));
  
  // world
  // WorldLV->SetVisAttributes(whiteBoxVisAtt);
  // collimatorLV->SetVisAttributes(yellowBoxVisAtt);
  // polyFilmLV1->SetVisAttributes(yellowBoxVisAtt);
  // cadmiumLV1->SetVisAttributes(whiteBoxVisAtt);
  // polyFilmLV2->SetVisAttributes(yellowBoxVisAtt);
  // cadmiumLV2->SetVisAttributes(whiteBoxVisAtt);
  
  boronFilmLV->SetVisAttributes (redBoxVisAtt);
  alSubstrateLV->SetVisAttributes (blueBoxVisAtt);
  xenonLV->SetVisAttributes (greenBoxVisAtt);
  casingLV->SetVisAttributes (grayVisAtt);
  // cdSlitLV->SetVisAttributes (yellowBoxVisAtt);
  hdpeLV->SetVisAttributes (lightblueVisAtt);
  
  

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

