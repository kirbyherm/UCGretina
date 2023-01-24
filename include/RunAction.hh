//
//

#ifndef RunAction_h
#define RunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"
#include "TTree.h"
#include "TFile.h"


extern TFile *newfile;
extern TTree *t;
extern TBranch *ebranch;

extern G4double energy[10], energy_tot;
extern int mult, nDetectors;
extern G4String detectorName[10];
extern G4double electronKE_mc2;
extern G4double electronKE_keV;
extern G4double x,y,z;

class G4Run;

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

  private:
    G4int runIDcounter;
};

#endif

