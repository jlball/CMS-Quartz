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

  G4Material* Silicon = nistManager->FindOrBuildMaterial("G4_Si", fromIsotopes);
  G4Material* Xenon = nistManager->FindOrBuildMaterial("G4_Xe", fromIsotopes);
  // Al for substract, Aluminum
  G4Material* Aluminum = nistManager->FindOrBuildMaterial("G4_Al", fromIsotopes);

  /*
   // Vacuum
   new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
   kStateGas, 2.73*kelvin, 3.e-18*pascal);
   */
  
  // Cd for SETN container, Cadmium
  G4Material* cadmium = nistManager->FindOrBuildMaterial("G4_Cd", fromIsotopes);

  // Xe
  // G4Material* xenon = nistManager->FindOrBuildMaterial("G4_Xe", fromIsotopes);
  
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


  //Xenon Gas Dimensions:
  G4double Xenon_x = 7 *mm;
  G4double Xenon_y = 25.4 *mm;
  G4double Xenon_z = 25.4 *mm;

  //Aluminum Substrate:
  G4double Al_x = 0.8 *mm;
  G4double Al_y = 25.4 *mm;
  G4double Al_z = 25.4 *mm;

  //Boron-10 Film:
  G4double Boron_x = 0.0001 *mm;
  G4double Boron_y = 25.4 *mm;
  G4double Boron_z = 25.4 *mm;

  //SiPM:
  G4double SiPM_x = 3*mm; 
  G4double SiPM_y = 1*mm; 
  G4double SiPM_z = 3*mm;

  //SOLIDS:
  G4VSolid* XenonSolid = new G4Box("XenonGas", Xenon_x/ 2, Xenon_y/2, Xenon_z/2);
  G4VSolid* AlSubstrateSolid = new G4Box("AlSubstrate", Al_x/ 2, Al_y/2, Al_z/2);
  G4VSolid* BoronFilmSolid = new G4Box("BoronFilm", Boron_x/ 2, Boron_y/2, Boron_z/2);

  G4VSolid* SiPMSolid = new G4Box("SiPM", SiPM_x/2, SiPM_y/2, SiPM_z/2);

  //LOGICAL VOLUMES:
  G4LogicalVolume* XenonLogical = new G4LogicalVolume(XenonSolid, Xenon, "XenonGas");
  G4LogicalVolume* AlSubstrateLogical = new G4LogicalVolume(AlSubstrateSolid, Aluminum, "AlSubstrate");
  G4LogicalVolume* BoronFilmLogical = new G4LogicalVolume(BoronFilmSolid, boron10, "BoronFilm");

  G4LogicalVolume* SiPMLogical = new G4LogicalVolume(SiPMSolid, Silicon, "SiPM");

  //Place volumes in the world:
  new G4PVPlacement(0, G4ThreeVector(), XenonLogical, "XenonGas", World, false, 0);
  new G4PVPlacement(0, G4ThreeVector(Xenon_x/2 + Al_x/2, 0, 0), AlSubstrateLogical, "AlSubstrate", World, false, 0);
  new G4PVPlacement(0, G4ThreeVector(-Xenon_x/2 - Al_x/2, 0, 0), AlSubstrateLogical, "AlSubstrate", World, false, 1);
  new G4PVPlacement(0, G4ThreeVector(Xenon_x/2 + Al_x + Boron_x/2, 0, 0), BoronFilmLogical, "BoronFilm", World, false, 0);

  new G4PVPlacement(0, G4ThreeVector(0, Xenon_y/2 + SiPM_y/2, Xenon_z/8), SiPMLogical, "SiPM", XenonLogical, false, 0);
  new G4PVPlacement(0, G4ThreeVector(0, -Xenon_y/2 - SiPM_y/2, Xenon_z/8), SiPMLogical, "SiPM", XenonLogical, false, 1);

  new G4PVPlacement(0, G4ThreeVector(0, Xenon_y/2 + SiPM_y/2, Xenon_z/8 + Xenon_z/4), SiPMLogical, "SiPM", XenonLogical, false, 2);
  new G4PVPlacement(0, G4ThreeVector(0, -Xenon_y/2 - SiPM_y/2, Xenon_z/8 + Xenon_z/4), SiPMLogical, "SiPM", XenonLogical, false, 3);

  new G4PVPlacement(0, G4ThreeVector(0, Xenon_y/2 + SiPM_y/2, -Xenon_z/8), SiPMLogical, "SiPM", XenonLogical, false, 4);
  new G4PVPlacement(0, G4ThreeVector(0, -Xenon_y/2 - SiPM_y/2, -Xenon_z/8), SiPMLogical, "SiPM", XenonLogical, false, 5);

  new G4PVPlacement(0, G4ThreeVector(0, Xenon_y/2 + SiPM_y/2, -Xenon_z/8 - Xenon_z/4), SiPMLogical, "SiPM", XenonLogical, false, 6);
  new G4PVPlacement(0, G4ThreeVector(0, -Xenon_y/2 - SiPM_y/2, -Xenon_z/8 - Xenon_z/4), SiPMLogical, "SiPM", XenonLogical, false, 7);

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

  BoronFilmLogical->SetVisAttributes (redBoxVisAtt);
  AlSubstrateLogical->SetVisAttributes (blueBoxVisAtt);
  XenonLogical->SetVisAttributes (greenBoxVisAtt);
  SiPMLogical->SetVisAttributes (yellowBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

