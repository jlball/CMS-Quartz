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
/// \file DOWSER01EventAction.cc
/// \brief Implementation of the DOWSER01EventAction class

#include "DOWSER01EventAction.hh"
#include "DOWSER01Analysis.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4GenericMessenger.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01EventAction::DOWSER01EventAction()
 : G4UserEventAction(),
   fMessenger(0),
   fPrintModulo(100000)
{
  // Define /DOWSER01/event commands using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/DOWSER01/event/", "Event control");

  // Define /DOWSER01/event/setPrintModulo command
  /*
    G4GenericMessenger::Command& setPrintModulo
    = fMessenger->DeclareProperty("setPrintModulo", 
                                  fPrintModulo, 
    = fMessenger->DeclareProperty("setPrintModulo",
                                  fPrintModulo,
                                 "Print events modulo n");
  setPrintModulo.SetRange("value>0");
    */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DOWSER01EventAction::~DOWSER01EventAction()
{
  delete fMessenger;
}



void DOWSER01EventAction::BeginOfEventAction(const G4Event* event)
{  
{

  G4int eventID = event->GetEventID();
  if ( eventID % fPrintModulo == 0 )  { 
  if ( eventID % fPrintModulo == 0 )  {
    G4cout << "\n---> Begin of event: " << eventID << G4endl;
    //CLHEP::HepRandom::showEngineStatus();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DOWSER01EventAction::EndOfEventAction(const G4Event* event)
{  
 
{

  G4int eventID = event->GetEventID();

  if ( eventID % fPrintModulo == 0) {
    G4cout << "---> End of event: " << eventID << G4endl;     
    G4cout << "---> End of event: " << eventID << G4endl;

  }

}  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
