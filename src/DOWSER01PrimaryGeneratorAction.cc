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
/// \file DOWSER01PrimaryGeneratorAction.cc
/// \brief Implementation of the DOWSER01PrimaryGeneratorAction class

#include "DOWSER01PrimaryGeneratorAction.hh"
#include "DOWSER01Analysis.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "DOWSER01EventInformation.hh" // Event Info

// #include "G4CsvAnalysisManager.hh"  // 2013-02-04

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01PrimaryGeneratorAction::DOWSER01PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);
  // nAngle = 0.;
  G4double Emax = 6; // 1 MeV;
  G4double Emin = -0.39794; // 0.4 eV;
  G4double nEnergyT(0);
  fnArray = 1295;
  fAngleArray = 180;
  G4double fEnergyT = 0;
  pi = 3.14159265358979323846264338328;
  // Moon radius: 1737 km, LRO altitude: orbitH km
  orbitH = 50.0; // ground (km)
  moonR = 1737.0;
  RoverRH = (moonR/(moonR + orbitH));
  // vLRO = sqrt(G * moon_Mass /((orbitH+1737)*1000)); // LEO velocity in Y direction
  vLRO = sqrt(4902808346033.4/((orbitH+moonR)*1000)); //
  nEqvEnergy = 0.5*nMass*vLRO*vLRO/1.6022E-19; // eV
  // neutron potential difference
  // del U = G * n_Mass * moon_Mass * (1 / R_moon - 1 / (R_moon + H))
  nPotDiff = 51252.715*(1./1737000 - 1./(1737000 + orbitH*1000)); // potential difference: U_orbit - U_surf
  nMass = 1.6749E-27; // kg
  moonMass = 7.34767309E+22; // kg
  thetaCr = asin(RoverRH);
  // nMC2 = nMass * 299790000.0 * 299790000.0; // joule
  
  std::ifstream theData("LNeutronIntSpectrum.DAT", std::ios::in);
  G4cout << "Open Neutron Spectrum Data " << G4endl;
  
  

  
  G4ParticleDefinition* particleDefinition
  = G4ParticleTable::GetParticleTable()->FindParticle("electron");
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(0.001032720*eV);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01PrimaryGeneratorAction::~DOWSER01PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DOWSER01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event
  G4double Eneutron(0);
  G4int copyNo(0), iflag(0);
  G4double xx(0), yy(0), zz(0), pPhiH(0), pVerH(0), xDiff(0), yDiff(0), zDiff(0);
  // G4double xx(0), yy(0), zz(0), pPhiH(0), pVerH(0);
  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore
  //
  G4double worldZHalfLength = 0;
  G4LogicalVolume* worlLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* worldBox = 0;
  if ( worlLV) worldBox = dynamic_cast< G4Box*>(worlLV->GetSolid());
  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();
  }
  else  {
    G4cerr << "World volume of box not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }
  
  // begining of particle iteration - energy, position, angle assessment
  
  // NEW - SET EVENT INFORMATION
  // Required files:
  // DOWSER01EventInformation.hh
  
  DOWSER01EventInformation* info = new DOWSER01EventInformation();
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // begining of energy assigment
  
  //Set neutron energy
  Eneutron = 2.45;
  
  fParticleGun->SetParticleEnergy(Eneutron*MeV);
  analysisManager->FillH1(1, std::log10(Eneutron));
  info->SetEnergyN(Eneutron);
  
  G4double x0 = 1.0*m;
  G4double y0 = 0*m;
  G4double z0 = 0*m;

  G4double px0 = -1*m;
  G4double py0 = 0*m;
  G4double pz0 = 0*m;
 
  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px0,py0,pz0));
  //analysisManager->FillH1(2, rotationAngle*180.0/pi);
  fParticleGun->GeneratePrimaryVertex(anEvent);
  anEvent->SetUserInformation(info);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


