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
/// \file DOWSER01SteppingAction.hh
/// \brief Definition of the DOWSER01SteppingAction class

#ifndef DOWSER01SteppingAction_h
#define DOWSER01SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4LogicalVolume;

/// Stepping action class
///
/// It holds data member fEnergy for accumulating the energy deposit
/// in a selected volume step by step.
/// The selected volume is set from  the detector construction via the
/// SetVolume() function. The accumulated energy deposit is reset for each
/// new event via the Reset() function from the event action.

class DOWSER01SteppingAction : public G4UserSteppingAction
{
public:
  DOWSER01SteppingAction();
  virtual ~DOWSER01SteppingAction();
  
  // static access method
  static DOWSER01SteppingAction* Instance();
  
  // method from the base class
  virtual void UserSteppingAction(const G4Step*);
  
  // reset accumulated energy
  void Reset();
  
  // set methods
  void SetVolume(G4LogicalVolume* volume) { fVolume = volume; }
  
  // get methods
  G4LogicalVolume* GetVolume() const { return fVolume; }
  G4double GetEnergy() const { return fEnergy; }
  // G4double GetPEnergy() const { return fPEnergy; }
  // G4double GetTEnergy() const { return fTEnergy; }
  
  
private:
  static DOWSER01SteppingAction* fgInstance;
  
  G4LogicalVolume* fVolume;
  G4double  fEnergy;
  G4double  kEnergy;
  G4double nEnergy;
  G4double theta;
  G4double phi;
  G4double px;
  G4double py;
  G4double pz;
  G4double x0, y0, z0, x1, y1, z1, nx0, ny0, nx1, ny1, ntheta, theta0, theta1, Li_x1;
  G4String  fParticleName;
  G4String  fParticleNameOld;
  G4int numOfCapture;
  
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
