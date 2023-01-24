#include "globals.hh"

class SuNBetaDecay {
public:

  G4double electronRestMassEnergy; // The rest mass energy of the electron in keV.
  
  SuNBetaDecay(G4int input_zDN, G4int input_aDN, G4double input_Qvalue){
    zDN = input_zDN;
    aDN = input_aDN;
    Qvalue = input_Qvalue;
    electronRestMassEnergy = 510.99891;
  }
  
  ~SuNBetaDecay()
  {}
  
  std::complex<G4double> ComplexGammaFunction(const G4double a, const G4double b);

  G4double GetMaximumElectronKineticEnergy();

  G4double GetRandomElectronKineticEnergy();

  G4double GetPhaseSpaceElectronKineticEnergyDistribution(const G4double KE);

  G4double GetPhaseSpaceElectronKineticEnergy();
           
  G4double GetApproximateFermiFunctionA(const G4double KE);

  G4double GetApproximateFermiFunctionB(const G4double KE);

  G4double GetApproximateFermiFunctionC(const G4double KE);

  G4double GetCorrectedElectronKineticEnergy(const G4String flag);

private:

  G4int zDN; // The atomic number of the daughter nucleus in the beta decay (D = daughter, N = nucleus)
  
  G4int aDN; // The mass number of the daughter nucleus in the beta decay (D = daughter, N = nucleus)       
  
  G4double Qvalue; // The ground-state-to-ground-state Q value in keV for the beta decay. This value is from the NNDC.
  
};
