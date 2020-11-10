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
// ********************************************************************
// * 2013-02-03, J.J. Su                                              *
// *                                                                  *
// * This module is to simulate CSETN count rate                       *
// * 20 bar pressure He3 detector to detect epithermal neutrons       *
// * detector container is made from Cadmium                          *
// * container dimension 60 mm dia x 113 mm L                         *
// * thickness of Cd container: 1.5 mm                                *
// * front end plate: material - UNKNOWN (try ceremic insulator)      *
// * est 1.5mm thick - NOT IMPLEMENTED YET.                           *
// * end end plate same as front end plate, thickness ~ 4 mm          *
// * a conducting copper plat 0.8 mm - NOT IMPLEMENTED YET.           *
// *                                                                  *
// ********************************************************************
//
// $Id$
//
/// \file DOWSER01.cc
/// \brief Main program of the DOWSER01 example

#include "DOWSER01DetectorConstruction.hh"
#include "DOWSER01PrimaryGeneratorAction.hh"
#include "DOWSER01RunAction.hh"
#include "DOWSER01SteppingAction.hh"
#include "DOWSER01EventAction.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
// #include "FTFP_BERT.hh"
#include "QGSP_BERT_HP.hh"
#include "Randomize.hh"
#include <time.h>

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
 
  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // randomize different seed for each run, 2013-02-02; J.J. Su
  G4long seed = time(0);
  CLHEP::HepRandom::setTheSeed(seed);
  CLHEP::HepRandom::showEngineStatus();


  // Set mandatory initialization classes
  //
  DOWSER01DetectorConstruction* detConstruction = new DOWSER01DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  // or the initialization can be done in this manner
  // runManager->SetUserInitialization(new DOWSER01DetectorConstruction());
  //

  // change physics list to QGSP_BERT_HP; 2013-02-02; J.J. Su
  // G4VModularPhysicsList* physicsList = new FTFP_BERT;
  G4VModularPhysicsList* physicsList = new QGSP_BERT_HP;

  // if step limiter is required then add this statement; 2013-02-02; J.J. Su
  // physicsList->RegisterPhysics(new G4StepLimiterBuilder());
  runManager->SetUserInitialization(physicsList);
    
  // Set user action classes
  //
  runManager
    ->SetUserAction(new DOWSER01PrimaryGeneratorAction());
  //
  runManager->SetUserAction(new DOWSER01RunAction());
  // 2013-0226; turn off SteppingAction tracking
  runManager->SetUserAction(new DOWSER01SteppingAction()); // drop SteppingAction tracking
  //
  runManager->SetUserAction(new DOWSER01EventAction());
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
  
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  else  {  
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac"); 
#else
    UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
    if (ui->IsGUI())
      UImanager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();
    delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
