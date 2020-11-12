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
  
  if(!theData)
  {
    theData.close();
    G4cout << "LNeutronIntSpectrum.DAT does not exit! " << G4endl;
    G4double dElog10 = (Emax - Emin)/fnArray;
    for (G4int copyNo=0; copyNo<fnArray; copyNo++)
    {
      nEnergy[copyNo]=pow(10,Emin + copyNo * dElog10);
      intSpectrum[copyNo] = pow(nEnergy[copyNo],-0.89)
      *(pow(10,Emin + (copyNo+0.5) * dElog10)-pow(10,Emin + (copyNo-0.5) * dElog10));
      nEnergyT += intSpectrum[copyNo];
      intSpectrum[copyNo] = nEnergyT;
    }
    intSpectrum[0] = 0;
    for (G4int copyNo=0; copyNo<fnArray; copyNo++)
    {
      intSpectrum[copyNo] /=  nEnergyT;
    }
    
  } else
  {
    for (G4int copyNo=0; copyNo<fnArray; copyNo++)
    {
      if (!theData.eof())
      {
        theData >> nEnergy[copyNo] >> intSpectrum[copyNo]  ;
      }
    }
  }
  theData.close();
  
  // Lunar neutron angular distribution, LEND's View
  fEnergyT = 0;
  for (G4int copyNo = 0; copyNo < fAngleArray; copyNo++)
  {
    // nFluxHEE[copyNo] = cos(0.5*copyNo*pi/fAngleArray)*std::sqrt(cos(0.5*copyNo*pi/fAngleArray))*sin(0.5*copyNo*pi/fAngleArray);
    nFluxHEE[copyNo] = std::pow(cos(0.5*copyNo*pi/fAngleArray), 3./2.)*sin(0.5*copyNo*pi/fAngleArray);
    nThetaHEE[copyNo] = copyNo*90./fAngleArray;
    fEnergyT += nFluxHEE[copyNo];
    nFluxHEE[copyNo] = fEnergyT;
  }
  
  for (G4int copyNo = 0; copyNo < fAngleArray; copyNo++)
  {
    nFluxHEE[copyNo] /=  fEnergyT;
  }
  
  
  G4ParticleDefinition* particleDefinition
  = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
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
  G4double Eneutron(0), r0(0), x0(0), y0(0), z0(0), vx(0), vy(0), vz(0), theta(0), groundAngle(0),
  P0(0), px0(0), py0(0), pz0(0), rotationAngle(0), theta0(0), deltheta(0), phi0(0);
  G4int copyNo(0), iflag(0);
  G4double xx(0), yy(0), zz(0), pPhiH(0), pVerH(0);
  // G4double xx(0), yy(0), zz(0), pPhiH(0), pVerH(0);
  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore
  //
  G4double worldZHalfLength = 0;
  G4LogicalVolume* worlLV
  = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
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
  // Force to assign energy = 0.5 MeV
  
  // P0 = G4UniformRand();
  // for (int ii = 0; ii < fnArray; ii++)
  // {
  //   if (intSpectrum[ii] >= P0)
  //   {
  //     Eneutron = nEnergy[ii-1] + (P0 - intSpectrum[ii-1])*(nEnergy[ii] - nEnergy[ii-1])/(intSpectrum[ii] - intSpectrum[ii-1]);
  //     break;
  //   }
  // }
  
  Eneutron = 0.5;
  
  // P0 = G4UniformRand();
  // if (P0 <= 0.8286)
  // {
  //   Eneutron = CLHEP::RandGauss::shoot(0.00103272, 0.0000430);
  // } else
  // {
  //   Eneutron = CLHEP::RandGauss::shoot(0.0092738, 0.00042453);
  // }

  // Eneutron = 0.5E6;
  fParticleGun->SetParticleEnergy(Eneutron*eV);
  analysisManager->FillH1(1, std::log10(Eneutron));
  info->SetEnergyN(Eneutron);
  
  rotationAngle = 0; // 1.570796327 rotating angle of beam source
  // rotationAngle = -60*0.017453293;
  // rotationAngle = -0.25 *3.14159265358979323846264338328;
  phi0 = 0.0; // G4UniformRand() *3.14159265358979323846264338328*2.0; // 1.570796327 rotating angle of beam source
  // rotationAngle = 0.5 *3.14159265358979323846264338328*0.5; // 1.570796327 rotating angle of beam source
  // phi0 = 0.*G4UniformRand() *3.14159265358979323846264338328*2.0; // 1.570796327 rotating angle of beam source
  x0 = 10; // +5.875;
  y0 = -1.25 + 2.5 * G4UniformRand();
  z0 = -1.25 + 2.5 * G4UniformRand();
  xx = x0*cos(rotationAngle)*cos(phi0) + z0*cos(phi0)*sin(rotationAngle) -  y0*sin(phi0);
  yy = y0*cos(phi0) + x0*cos(rotationAngle)*sin(phi0) + z0*sin(rotationAngle)*sin(phi0);
  zz = z0*cos(rotationAngle) - x0*sin(rotationAngle) + 0.0;
  px0 = -1; // *cos(phi0);
  py0 = 0; // py0=-sin(rotationAngle)*sin(phi0);
  pz0 = 0;
  fParticleGun->SetParticlePosition(G4ThreeVector(x0*cm, y0*cm, z0*cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px0,py0,pz0));
  info->SetTheta0(rotationAngle*180.0/pi);
  info->SetTheta1(phi0*180.0/pi);
  analysisManager->FillH1(2, rotationAngle*180.0/pi);
  fParticleGun->GeneratePrimaryVertex(anEvent);
  anEvent->SetUserInformation(info);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


