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
#include "G4Tubs.hh"
#include "G4Tubss.hh"
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

#include "G4tgbRotationMatrix.hh"

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


  // Al for substract, Aluminum
  G4Material* Aluminum = nistManager->FindOrBuildMaterial("G4_Al", fromIsotopes);


  // Print materials
  G4cout << G4endl << "The materials defined are: " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  G4cout << "************  end of material table  ************" << G4endl;

  //-------- Custom Materials ----------------------------------------------------
  G4double z, a, fractionmass, density;
  G4String name, symbol;
  G4int ncomponents;

  // Call elements from GEANT4 database
  G4NistManager* man = G4NistManager::Instance();
  G4Element* Al = man->FindOrBuildElement("G4_Al");
  G4Element* As = man->FindOrBuildElement("G4_As");
  G4Element* B = man->FindOrBuildElement("G4_B");
  G4Element* Ca = man->FindOrBuildElement("G4_Ca");
  G4Element* Cd = man->FindOrBuildElement("G4_Cd");
  G4Element* Cr = man->FindOrBuildElement("G4_Cr");
  G4Element* Cu = man->FindOrBuildElement("G4_Cu");
  G4Element* Fe = man->FindOrBuildElement("G4_Fe");
  G4Element* K = man->FindOrBuildElement("G4_K");
  G4Element* Li = man->FindOrBuildElement("G4_Li");
  G4Element* Mg = man->FindOrBuildElement("G4_Mg");
  G4Element* Na = man->FindOrBuildElement("G4_Na");
  G4Element* Ni = man->FindOrBuildElement("G4_Ni");
  G4Element* P = man->FindOrBuildElement("G4_P");
  G4Element* Sb = man->FindOrBuildElement("G4_Sb");
  G4Element* Ti = man->FindOrBuildElement("G4_Ti");
  G4Element* Zr = man->FindOrBuildElement("G4_Zr");
  G4Element* O = man->FindOrBuildElement("G4_O");
  G4Element* H = man->FindOrBuildElement("G4_H");
  G4Element* Si = man->FindOrBuildElement("G4_Si");

  G4Element* Cl = man->FindOrBuildElement("G4_Cl");
  G4Element* Mn = man->FindOrBuildElement("G4_Mn");
  G4Element* Ge = man->FindOrBuildElement("G4_Ge");
  G4Element* La = man->FindOrBuildElement("G4_La");
  G4Element* Sm = man->FindOrBuildElement("G4_Sm");
  G4Element* Eu = man->FindOrBuildElement("G4_Eu");
  G4Element* Dy = man->FindOrBuildElement("G4_Dy");
  G4Element* Yb = man->FindOrBuildElement("G4_Yb");
  G4Element* W = man->FindOrBuildElement("G4_W");
  G4Element* U = man->FindOrBuildElement("G4_U");

  // Creating hydroxyl
  G4double z, a, density;
  G4String name, symbol;
  G4int ncomponents, natoms;
  density = 1.000*g/cm3;
  G4Material* OH = new G4Material(name="OH",density,ncomponents=2);
  OH->AddElement(H, natoms=1);
  OH->AddElement(O, natoms=1);
//  G4double z, a, density;
//  G4String name, symbol;
//  G4int ncomponents, natoms;
//  density = 1.000*g/cm3;
//  G4Material* OH = new G4Material(name="OH",density,ncomponents=2);
//  OH->AddElement(H, natoms=1);
//  OH->AddElement(O, natoms=1);

  // Creating SiO2
  G4double z, a, density;
  G4String name, symbol;
  G4int ncomponents, natoms;
  density = 2.65*g/cm3;
  G4Material* SiO2 = new G4Material(name="SiO2",density,ncomponents=2);
  SiO2->AddElement(Si, natoms=1);
  SiO2->AddElement(O, natoms=2);

  // Creating Fused Quartz
  density = 2200 kg/m3;
  G4Material* GE214quartz = new G4Material(name="GE214quartz",density,ncomponents=20);
  ge214quartz->AddMaterial(SiO2, fractionmass=99.97955197*perCent);
  ge214quartz->AddMaterial(SiO2, fractionmass=99.97862127*perCent);
  ge214quartz->AddElement(Al, fractionmass=0.014*perCent);
  ge214quartz->AddElement(As, fractionmass=0.000002*perCent);
  ge214quartz->AddElement(B, fractionmass=0.00002*perCent);
  ge214quartz->AddElement(Ca, fractionmass=0.00004*perCent);
  ge214quartz->AddElement(Cd, fractionmass=0.000001*perCent);
  ge214quartz->AddElement(Cr, fractionmass=0.000005*perCent);
  ge214quartz->AddElement(Cu, fractionmass=0.000005*perCent);
  ge214quartz->AddElement(Fe, fractionmass=0.00002*perCent);
  ge214quartz->AddElement(K, fractionmass=0.00006*perCent);
//  ge214quartz->AddElement(K, fractionmass=0.00006*perCent);
  ge214quartz->AddMaterial(K, fractionmass=0.000586*perCent);
  ge214quartz->AddElement(Li, fractionmass=0.00006*perCent);
  ge214quartz->AddElement(Mg, fractionmass=0.00001*perCent);
  ge214quartz->AddElement(Mn, fractionmass=0.000005*perCent);
  ge214quartz->AddElement(Na, fractionmass=0.00007*perCent);
  ge214quartz->AddElement(Ni, fractionmass=0.00001*perCent);
  ge214quartz->AddElement(P, fractionmass=0.00002*perCent);
  ge214quartz->AddElement(Sb, fractionmass=0.00000003*perCent);
  ge214quartz->AddElement(Ti, fractionmass=0.0011*perCent);
  ge214quartz->AddElement(Zr, fractionmass=0.00008*perCent);
  ge214quartz->AddMaterial(OH, fractionmass=0.005*perCent);
//  ge214quartz->AddMaterial(OH, fractionmass=0.005*perCent);
  ge214quartz->AddMaterial(Cl, fractionmass=0.0002113*perCent);
  ge214quartz->AddMaterial(Mn, fractionmass=0.0000194*perCent);
  ge214quartz->AddMaterial(Ge, fractionmass=0.0000011*perCent);
  ge214quartz->AddMaterial(La, fractionmass=0.0000158*perCent);
  ge214quartz->AddMaterial(Sm, fractionmass=0.0000053*perCent);
  ge214quartz->AddMaterial(Eu, fractionmass=0.0000008*perCent);
  ge214quartz->AddMaterial(Dy, fractionmass=0.0000366*perCent);
  ge214quartz->AddMaterial(Yb, fractionmass=0.0000125*perCent);
  ge214quartz->AddMaterial(W, fractionmass=0.0000375*perCent);
  ge214quartz->AddMaterial(U, fractionmass=0.0000644*perCent);

  //============================================================================
  //      Definitions of Solids, Logical Volumes, Physical Volumes
  //============================================================================

  //-------------
  // World Volume
  //-------------
  G4RotationMatrix* rot = new G4RotationMatrix();
  rot->rotateZ(0.*deg);

  G4ThreeVector worldSize = G4ThreeVector(150*cm, 150*cm, 150*cm);
  G4Tubs * solidWorld
  = new G4Tubs("soildWorld", worldSize.x()/2., worldSize.y()/2., worldSize.z()/2.);
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


  //SOLIDS:


  G4VSolid* QUartzSolid 
  = new G4Tubs("GE214Quartz". 
  G4VSolid* QUartzSolid
  = new G4Tubs("GE214Quartz".
                0*mm, // Inner radius
                100*mm, // Outer radius
                100*mm, // Half length Z (hz)
                0.*deg, // start angle
                360.*deg); // spanningAngle

  G4VSolid* QUartzSolid = new G4Box("Hello", 10*mm, 10*mm, 10*mm);



  //LOGICAL VOLUMES:

  G4LogicalVolume* QuartzLogical = new G4LogicalVolume(QuartzSolid, GE214Quartz, "GE214Quartz");



  //Place volumes in the world:

  new G4PVPlacement(0,                       // no rotation (?)
                    G4ThreeVector(0, 0, 0),  // translation position
                    QuartzLogical,             // logical volume
                    "GE214Quartz",           // volume name
                    World,                   // mother volume
                    false,                   // boolean operations (none)
                    0);                      // its copy number

  new G4VPlacement()
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

  XenonLogical->SetVisAttributes (redBoxVisAtt);
  QuartzLogical->SetVisAttributes (greenBoxisAtt);
  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

