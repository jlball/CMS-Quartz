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
/// \file DOWSER01PrimaryGeneratorAction.hh
/// \brief Definition of the DOWSER01PrimaryGeneratorAction class

#ifndef DOWSER01PrimaryGeneratorAction_h
#define DOWSER01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the calorimeter
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class
/// (see the macros provided with this example).

class DOWSER01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    DOWSER01PrimaryGeneratorAction();
    virtual ~DOWSER01PrimaryGeneratorAction();
    
    virtual void GeneratePrimaries(G4Event* event);
    
    // set methods
    void SetRandomFlag(G4bool value);
    
private:
    G4ParticleGun*  fParticleGun; // G4 particle gun
    G4int fnArray;
    G4int fAngleArray;
    G4double intSpectrum[1295]; // dimension fnArray
    G4double nEnergy[1295];  // dimension fnArray
    G4double nFluxHEE[180];
    G4double nThetaHEE[180];
    // G4double nFluxLEE[180];
    // G4double nThetaLEE[180];
    G4double pi;
    G4bool ifFileInput; // check if input neutron spectrum file was establised
    G4bool ifAngular;  // check if the neutron angular distribution file was established
    G4double nAngle; // neutron incident angle
    // G4double nEnergy; // neutron incident energy
    G4double RoverRH, orbitH, moonR; // R/(R + H), R: Moon radius, H: altitude
    G4double vLRO, nEqvEnergy; // LRO velocity (in Y direction)
    G4double nPotDiff; // neutron potential difference between the moon surface and LRO altitude
    G4double nMass, moonMass, thetaCr;
    
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


