#ifdef LHTARGET
#include "Target_LH_Messenger.hh"

Target_Messenger::Target_Messenger(Target* Tar)
:aTarget(Tar)
{ 
 
  TargetDir = new G4UIdirectory("/Target/");
  TargetDir->SetGuidance("Target control.");
  
  CCmd = new G4UIcmdWithAString("/Target/Cell",this);
  CCmd->SetGuidance("Select the target cell (thick: 200 mg/cm^2, thin: 50 mg/cm^2, empty: cell body, no windows), notarget: beam pipe only");
  CCmd->SetParameterName("choice",false);
  CCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  BCmd = new G4UIcmdWithADoubleAndUnit("/Target/Bulge",this);
  BCmd->SetGuidance("Set the thickness of the bulges.)");
  BCmd->SetParameterName("choice",false);
  BCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  WCmd = new G4UIcmdWithoutParameter("/Target/Windows",this);
  WCmd->SetGuidance("Construct the Kapton windows.)");

  ACmd = new G4UIcmdWithADoubleAndUnit("/Target/Angle",this);
  ACmd->SetGuidance("Set the angle of tilt of the target.)");
  ACmd->SetParameterName("choice",false);
  ACmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MatCmd = new G4UIcmdWithAString("/Target/Material",this);
  MatCmd->SetGuidance("Select material for the target.");
  MatCmd->SetParameterName("choice",false);
  MatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SDTarCmd = new G4UIcmdWithADouble("/Target/SetDensity",this);
  SDTarCmd->SetGuidance("Set target density in mg/cc.");
  SDTarCmd->SetParameterName("choice",false);
  SDTarCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PosZCmd = new G4UIcmdWithADoubleAndUnit("/Target/SetPosition_Z",this);
  PosZCmd->SetGuidance("Select the position of the target along the beam axis (Z direction.");
  PosZCmd->SetParameterName("choice",false);
  PosZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NSCmd = new G4UIcmdWithAnInteger("/Target/NStep",this);
  NSCmd->SetGuidance("Select the number of steps in the target");
  NSCmd->SetParameterName("choice",false);
  NSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  sFCmd = new G4UIcmdWithAString("/Target/sourceFrame",this);
  sFCmd->SetGuidance("Select a source frame (eu152_Z2707, cs137_E2879, or co56_2012).");
  sFCmd->SetParameterName("choice",false);
  sFCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  sledCmd = new G4UIcmdWithoutParameter("/Target/Sled",this);
  sledCmd->SetGuidance("Construct the target sled.");  

  GCmd = new G4UIcmdWithoutParameter("/Target/Greta",this);
  GCmd->SetGuidance("Construct the Greta LH target");

  CutCmd = new G4UIcmdWithoutParameter("/Target/Cutaway",this);
  CutCmd->SetGuidance("Build a cutaway view of the beamline tee and flanges. For visualization only!");  

  RepCmd = new G4UIcmdWithoutParameter("/Target/Report",this);
  RepCmd->SetGuidance("Report target parameters");  

}



Target_Messenger::~Target_Messenger()
{

  delete TargetDir;
  delete CCmd;
  delete BCmd;
  delete WCmd;
  delete ACmd;
  delete MatCmd;
  delete SDTarCmd;
  delete PosZCmd;
  delete NSCmd;
  delete sFCmd;
  delete sledCmd;
  delete GCmd;
  delete CutCmd;
  delete RepCmd;
}


void Target_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == CCmd )
    { aTarget->setTargetCell(newValue); }
  if( command == BCmd )
    { aTarget->setBulgeThickness(BCmd->GetNewDoubleValue(newValue)); }
  if( command == WCmd )
    { aTarget->setWindows(); }
  if( command == ACmd )
    { aTarget->setTargetAngle(ACmd->GetNewDoubleValue(newValue)); }
  if( command == MatCmd )
   { aTarget->setMaterial(newValue);} 
  if( command == SDTarCmd )
    { aTarget->SetDensity(SDTarCmd->GetNewDoubleValue(newValue)); }
  if( command == PosZCmd )
    { aTarget->SetPositionZ(PosZCmd->GetNewDoubleValue(newValue));}
  if( command == NSCmd )
    { aTarget->setNStep(NSCmd->GetNewIntValue(newValue)); }
  if( command == sFCmd )
    { aTarget->setSourceFrame(newValue); }
  if( command == sledCmd )
    { aTarget->setSled(); }
  if( command == GCmd )
    { aTarget->SetGreta(); }
  if( command == CutCmd )
    { aTarget->setCutaway(); }
  if( command == RepCmd )
    { aTarget->Report(); }
}
#endif
