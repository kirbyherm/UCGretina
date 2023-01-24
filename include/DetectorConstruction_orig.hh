//
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"


extern  G4String detectorName[10];


class G4VPhysicalVolume;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction
{
   public:
      DetectorConstruction();
      ~DetectorConstruction();

   public:
    G4VPhysicalVolume* Construct();

   private:
      void DefineMaterials();
      G4Material *NaI, *Al, *N78O21Ar1, *Cr20Ni8Fe76, *C2F4, *C5O2H8;
  G4Material* Si;
  G4Material* Cu3Zn2;
  G4Material* Cu;
  G4Material* C2H4_FR4;
  G4Material* C2H4_rigid;
  G4Material* C2H4_flexible;
  G4Material* C12H22N2O2;
  G4Material* vacuum;
  
      G4VPhysicalVolume* ConstructDetector();
};
 
#endif

