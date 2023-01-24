//   
//

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include <vector>
#include <random>

extern int nEnergies, cascade;
//extern float inputArray[5100][50];
extern std::vector<std::vector<float>> inputArray;
extern double electronKE_mc2;
extern double electronKE_keV;
//extern double input_excitation_energies[5100];
extern std::vector<float> input_excitation_energies;
extern double x,y,z;
extern std::mt19937 generator;
extern std::normal_distribution<double> x_distribution;
extern std::normal_distribution<double> y_distribution;
extern std::gamma_distribution<double> z_distribution;


class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event *anEvent);

  private:
    G4ParticleGun* particleGun;
  G4ParticleGun* particleGun2;

};

#endif

