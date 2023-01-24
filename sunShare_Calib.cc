#include "G4RunManager.hh"
#include "G4UImanager.hh" 
#include "G4UIterminal.hh"
#include "G4ios.hh"  

//#include "LHEP_PRECO_HP.hh"
#include "QGSP_BIC_HP.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SuNBetaDecay.hh"
#include "PhysicsList.hh"

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <random>
#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TRint.h"
#include "TPluginManager.h"

#include "G4VisExecutive.hh"


// GLOBAL VARIABLES

  int runTimes;    //total number of events to simulate
  int nCascades; //number of gamma ray cascades from the excited state 
  int nEnergies;   //number of gamma rays in each cascade

	//if you change the following variables you also need to update all of the include/*.hh files
	int nDetectors = 10;       //number of detectors you want the energy for in your ROOT file
  G4String detectorName[10]; //the name of your detectors (see src/DetectorConstruction.cc)
  G4double energy_MeV[10];   //deposited energy in MeV
  G4double energy_keV[10];   //deposited energy in keV
  G4double sigma[10];        //resolution of your detectors
  G4double energy[10];       //energy that will be saved to the ROOT file
  G4double energy_tot;      //total energy deposited in all of your detectors of interest
  int mult;                 //multiplicity of detectors hit
G4double electronKE_mc2;
G4double electronKE_keV;
  double x,y,z;

  std::vector<std::vector<float>> inputArray; //array of energies to simulate [cascade][energies]
  std::vector<float> inputCutoff;    //keeps track of how many events there should be for each cascade
  int cascade=0;             //keeps track of which cascade the simulation is at

  std::vector<float> input_excitation_energies; // keeps track of the excitation energy for a cascade
  std::mt19937 generator;
  std::normal_distribution<double> x_distribution(-0.819725,6.7653875);
  std::normal_distribution<double> y_distribution(1.4020125,7.5242375);
  std::gamma_distribution<double> z_distribution(4.116,24.627);

  TFile *newfile;
  TTree *t;
  TBranch *ebranch;


// THE MAIN PROGRAM

int main(int argc, char** argv)
{

  // creates the file and tree for the .root output
  G4String output_filepath_and_name = argv[3];
  newfile = new TFile(output_filepath_and_name,"RECREATE");
  t = new TTree("t","output from geant");

  // Construct the default run manager
  G4RunManager *runManager = new G4RunManager;

  // set mandatory initialization classes
  DetectorConstruction *Det = new DetectorConstruction;
  runManager->SetUserInitialization(Det);
  PhysicsList *physicsList = new PhysicsList(Det);
  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new QGSP_BIC_HP);

  #ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
     visManager->Initialize();
  #endif

 // set optional user action class
  runManager->SetUserAction(new RunAction);

  // Read input file and create array of energies to simulate 
  G4String input_file = argv[2];
  std::ifstream fin(input_file);

  fin >> runTimes >> nCascades >> nEnergies;
  G4cout << "runTimes: " <<  runTimes <<  "; nCascades: " << nCascades << "; nEnergies: " << nEnergies << G4endl;
  float inputTemp; 
  for (int i=0; i<nCascades; i++)
  {
    std::vector<float> inputSubArray; //array of energies to simulate per cascade [energies]
    fin >> inputTemp;
    input_excitation_energies.push_back(inputTemp);
//     G4cout << "excitationEnergy: " << input_excitation_energies[i] << G4endl;
     fin >> inputTemp;
     inputSubArray.push_back(inputTemp);
//     input_excitation_energies[i] = inputArray[i][0];
//     G4cout << "Probability: " << inputSubArray[0] << G4endl;
     if (i==0)
       inputCutoff.push_back(inputSubArray[0]*(float(runTimes))/100.0);
     else
       inputCutoff.push_back(inputCutoff[i-1] + inputSubArray[0]*(float(runTimes))/100.0);
//     G4cout << "Total Probability: " << inputCutoff[i] << G4endl;
     for (int j=1; j<=nEnergies; j++)
      {
//        G4cout << "Cascade: " << i << "; Energy: " << j << G4endl;
        fin >> inputTemp;
        inputSubArray.push_back(inputTemp);
//        input_excitation_energies[i] += inputArray[i][j];
//        G4cout << inputSubArray[j] << G4endl;
      }
     inputArray.push_back(inputSubArray);
  }
 
  // set mandatory user action class
  runManager->SetUserAction(new PrimaryGeneratorAction);

  EventAction *eventAction = new EventAction;
  runManager->SetUserAction(eventAction);

  runManager->SetUserAction(new SteppingAction(Det)); 


  // Initialize G4 kernel
  runManager->Initialize();


  // for visualization purposes
  if (argc==1)
  {
    G4UIsession *Session = new G4UIterminal;

    //if you want visualization you have to uncomment the next line
    G4UImanager::GetUIpointer()->ApplyCommand("/control/execute vis.mac");

    Session->SessionStart();
    delete Session;
  }
  else
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];  
    G4UImanager::GetUIpointer()->ApplyCommand(command+fileName);
  }

  // run the program
  runManager->BeamOn(runTimes);


  // job termination
  fin.close();
  newfile->Close();

  G4cout << "Rootfile closed!" << G4endl;

  #ifdef G4VIS_USE
     delete visManager;
  #endif

  delete runManager;
  G4cout << "RunManager closed!" << G4endl;

  return 0;
}
