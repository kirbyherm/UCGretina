// 
// 

#ifndef EventAction_h
#define EventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"
#include <vector>

extern G4double energy_MeV[10], energy_keV[10], energy[10], sigma[10], energy_tot;
extern int mult, cascade, nCascades, nDetectors;
//extern float inputCutoff[5100];
extern std::vector<float> inputCutoff;


class G4Event;

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    ~EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

  private:
};

#endif

