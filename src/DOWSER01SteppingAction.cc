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
/// \file DOWSER01SteppingAction.cc
/// \brief Implementation of the DOWSER01SteppingAction class

#include "DOWSER01SteppingAction.hh"
#include "DOWSER01DetectorConstruction.hh"
#include "DOWSER01Analysis.hh"  // 2012-0403, J.J. Su
#include "DOWSER01EventInformation.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4CsvAnalysisManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01SteppingAction* DOWSER01SteppingAction::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01SteppingAction* DOWSER01SteppingAction::Instance()
{
    // Static acces function via G4RunManager
    
    return fgInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01SteppingAction::DOWSER01SteppingAction()
: G4UserSteppingAction(),
fVolume(0),
fEnergy(0.),
kEnergy(0),
nEnergy(0),
theta(0), phi(0),px(0),py(0),pz(0),x0(0), y0(0), z0(0), x1(0), y1(0), z1(0), nx0(0), ny0(0), nx1(0), ny1(0),
ntheta(0), theta0(0), theta1(0), fParticleName(), fParticleNameOld(), numOfCapture(0), Li_x1(0)
{
    fgInstance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01SteppingAction::~DOWSER01SteppingAction()
{
    fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DOWSER01SteppingAction::UserSteppingAction(const G4Step* step)
{
    
    // get volume of the current step
    G4int fStepNo;
    G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    fParticleName = step->GetTrack()->GetDefinition()->GetParticleName();
    // G4double edep = step->GetTotalEnergyDeposit();
    fStepNo = step->GetTrack()->GetCurrentStepNumber();
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    // G4CsvAnalysisManager* analysisManager = G4CsvAnalysisManager::Instance();

    //Getting name of the logical volume the step moves through
    G4String nameLogicVolume = volume->GetName();
    

  if (fStepNo == 1)
  {
    G4EventManager* evtMan = G4EventManager::GetEventManager();
    DOWSER01EventInformation* info = (DOWSER01EventInformation*)evtMan->GetConstCurrentEvent()->GetUserInformation();
    if (info==0){G4cout<<"Warning: No Event information recorded!"<<G4endl;}
    else {
      theta0 = info->GetTheta0();
      nEnergy = info->GetEnergyN();
    }

    kEnergy = step->GetPreStepPoint()->GetKineticEnergy()/CLHEP::eV;

    //G4cout << "launched " << fParticleName << " w/ energy: " << kEnergy << G4endl;
    //G4cout << "in logical volume: " << nameLogicVolume << G4endl;

    if (fParticleName == "Li7")
    {
      if (nameLogicVolume == "Detector1")
      {
        Li_x1 = step->GetPreStepPoint()->GetPosition().z()/CLHEP::mm;

        analysisManager->FillH1(5, theta0-90);
        analysisManager->FillH1(7, Li_x1);
        analysisManager->FillH2(3, theta0-90, Li_x1);

        if (0 < Li_x1 && Li_x1 < 50)
        {
          analysisManager->FillH1(8, theta0-90);
        }
        if (50 <= Li_x1 && Li_x1 < 100)
        {
          analysisManager->FillH1(10, theta0-90);
        }
        if (100 <= Li_x1 && Li_x1 < 150)
        {
          analysisManager->FillH1(12, theta0-90);
        }
        if (150 <= Li_x1 && Li_x1 < 200)
        {
          analysisManager->FillH1(14, theta0-90);
        }

      } 
      else if (nameLogicVolume == "Detector2")
      {
        Li_x1 = step->GetPreStepPoint()->GetPosition().z()/CLHEP::mm;

        analysisManager->FillH1(6, theta0-90);
        analysisManager->FillH1(7, Li_x1);
        analysisManager->FillH2(3, theta0-90, Li_x1);

        if (0 < Li_x1 && Li_x1 < 50)
        {
          analysisManager->FillH1(9, theta0-90);
        }
        if (50 <= Li_x1 && Li_x1 < 100)
        {
          analysisManager->FillH1(11, theta0-90);
        }
        if (100 <= Li_x1 && Li_x1 < 150)
        {
          analysisManager->FillH1(13, theta0-90);
        }
        if (150 <= Li_x1 && Li_x1 < 200)
        {
          analysisManager->FillH1(15, theta0-90);
        }
      }

      numOfCapture += 1;
      G4cout << "New B-10 capture event, total is: " << numOfCapture << G4endl;

      theta0 = info->GetTheta0();
      phi = info->GetTheta1();
      analysisManager->FillH1(3, theta0 - 90);
      analysisManager->FillH1(4, phi);
    }

    if (fParticleName == "alpha")
    {
      analysisManager->FillH1(1, kEnergy);
    }
  } 

  if (nameLogicVolume == "Detector1" && fParticleName == "neutron")
  {
    G4double stepEnergy = step->GetPostStepPoint()->GetKineticEnergy()/CLHEP::eV;

    if (stepEnergy != 0)
    {
        analysisManager->FillH1(2, stepEnergy);
        analysisManager->FillH2(1, theta0-90, stepEnergy);
    }
  }

  if (nameLogicVolume == "Detector2" && fParticleName == "neutron")
  {
    G4double stepEnergy = step->GetPostStepPoint()->GetKineticEnergy()/CLHEP::eV;

    if (stepEnergy != 0)
    {
        analysisManager->FillH1(2, stepEnergy);
        analysisManager->FillH2(2, theta0 - 90, stepEnergy);
    }
  }

  if (fParticleName == "alpha" || fParticleName == "Li7")
  {
    z1 = step->GetPreStepPoint()->GetPosition().z()/CLHEP::mm;
    x1 = step->GetPreStepPoint()->GetPosition().x()/CLHEP::mm;
    fParticleNameOld = fParticleName;    
  }

  if ((fParticleName == "alpha") && (step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) && ((nameLogicVolume == "xenonLV") ))
  {
    kEnergy = step->GetPreStepPoint()->GetKineticEnergy()/CLHEP::MeV;
    //analysisManager->FillH1(9, kEnergy);
  }
  if ((fParticleName == "Li7") && (step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) && ((nameLogicVolume == "xenonLV") ))
  {
    kEnergy = step->GetPreStepPoint()->GetKineticEnergy()/CLHEP::MeV;
    //analysisManager->FillH1(10, kEnergy);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DOWSER01SteppingAction::Reset()
{
    // G4cout << "neutron energy deposit " << fEnergy << "MeV" << G4endl;
    // fEnergy = 0.;
    G4cout << "Stepping Action Reset Called" << G4endl;
}

