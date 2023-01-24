////////////////////////////////////////////////////////////////////////////////////////////
// Description: This file tells the program what to do before and after each              //
//   gamma-ray cascade. For now the program does the following:                           //
//                                                                                        //
//   BEFORE: - Calculate the gamma-ray cascade to run.                                    //
//           - Set energies and multiplicities back to zero.                              //
//                                                                                        //
//   AFTER:  - Calculate the energy deposited based on your detector resolution function. //
//           - Calculate the multiplicity of detectors hit.                               //
//           - Fill the ROOT Tree.                                                        //
////////////////////////////////////////////////////////////////////////////////////////////

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include <stdio.h>
#include "Randomize.hh"
#include <vector>


EventAction::EventAction()
{}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  //status bar that prints out run number to screen
  if ((evt->GetEventID()+1) % 10000 == 0)                 
   G4cout << ">>> Event " << evt->GetEventID()+1 << G4endl;

 // figures out what gamma-ray cascade to run
  int j = evt->GetEventID()+1; 
  for (int i=0; i<nCascades; i++) 
  {
     if (i==0 && j<inputCutoff[0])
       cascade=i;
     else if (i>0 && j<inputCutoff[i] && j>=inputCutoff[i-1])
       cascade=i;
  }

  //set everything back to zero
  for(int i=0; i<nDetectors; i++) 
   { 
    energy_MeV[i] = 0.0*CLHEP::eV;
    energy_keV[i] = 0.0*CLHEP::keV;
    sigma[i] = 0.0*CLHEP::keV;
    energy[i] = 0.0*CLHEP::keV;
   }
  energy_tot=0.0*CLHEP::keV;
  mult=0;
}

G4double ratio_between_experiment_and_simulation(G4double energy){

  G4double ratio = -1.06515e-09 * pow(energy, 5.0)
    + 4.40364e-09 * pow(energy, 4.0)
    + 1.10134e-05 * pow(energy, 3.0)
    - 0.000300895 * pow(energy, 2.0)
    + 0.00439406 * energy
    + 0.00875016;
  return ratio;
}

void EventAction::EndOfEventAction(const G4Event* evt)
{

  // calculate energy and multiplicity to save to ROOT file
  for(int i=0; i<nDetectors; i++)
   {

     if (i < 8){ // SuN, so that total_energy and mult only correspond to SuN
     
       double threshold = 0.0; //experimental threshold in keV

       if (energy_MeV[i]*1000.0 > threshold)
	 {	
//        std::cout << energy_MeV[i] << std::endl;
	   energy_keV[i]=1000.0*energy_MeV[i];

	   //detector resolution function from ACDombos Thesis Fig3.23
	   sigma[i] = -7.85938e-15 * pow(energy_keV[i],4.0)
	     +2.42564e-10 * pow(energy_keV[i],3.0)
	     -2.92131e-06 * pow(energy_keV[i],2.0)
	     +0.0244892 * energy_keV[i]
	     +6.49392;

	   //detector resolution function from SQuinn Thesis Fig5.11
//	   sigma[i] = 4.57868e-15 * pow(energy_keV[i],4.0)
//	     -1.09955e-10 * pow(energy_keV[i],3.0)
//	     +4.09992e-07 * pow(energy_keV[i],2.0)
//	     +1.27552e-02 * pow(energy_keV[i],1.0)
//	     +17.78350;

	   //detector resolution function from used by WJ in proposal
//       sigma[i]=-5.59375e-15*pow(energy_keV[i],4.0)
//                +1.85975e-10*pow(energy_keV[i],3.0)
//                -2.47836e-6 *pow(energy_keV[i],2.0)
//                +2.33408e-2 *energy_keV[i]
//                +7.00328;


	   energy[i] = G4RandGauss::shoot(energy_keV[i],sigma[i]);

	   // fix threshold discrepancy between experiment and simulation

	   G4double threshold_discrepancy = 70;
        threshold_discrepancy = 0;
	   G4bool increment_multiplicity = true;
	   
	   if (energy[i] < threshold_discrepancy){
	     G4double ratio = ratio_between_experiment_and_simulation(energy[i]);
	     G4double random_number = ratio_between_experiment_and_simulation(threshold_discrepancy) * G4UniformRand();
//         G4cout<<"ratio: "<<ratio<<" ; rand: "<<random_number<<G4endl;
	     if (random_number > ratio){
	       energy[i] = 0;
	       increment_multiplicity = false;
	     }
	   }

	   energy_tot += energy[i];
	   if (increment_multiplicity){
	     ++ mult;
	   }
	 }
       else energy[i]=0.0; 
     }

     else { // DSSD and Veto, in case later on there are variables or a resolution function that only applies to them

       double threshold = 230.0; //experimental threshold in keV
       
       if (energy_MeV[i]*1000.0 > threshold)
	 {	
	   energy_keV[i]=1000.0*energy_MeV[i];

	   //detector resolution function
	   sigma[i] = 0.0077*energy_keV[i];

	   energy[i] = G4RandGauss::shoot(energy_keV[i],sigma[i]);

	 }
       else energy[i]=0.0; 
     }
   }

  t->Fill();

}
