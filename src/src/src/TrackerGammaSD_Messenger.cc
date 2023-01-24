#include "TrackerGammaSD_Messenger.hh"


TrackerGammaSD_Messenger::TrackerGammaSD_Messenger(TrackerGammaSD* TGSD)
:tracker(TGSD) 
{ 

  PrtGDir = new G4UIdirectory("/GammaPrint/");
  PrtGDir->SetGuidance("Event by event print control for.");

  PrtGSCmd = new G4UIcmdWithoutParameter("/GammaPrint/Track_Set",this);
  PrtGSCmd->SetGuidance("Sets printing of track gamma results.");

  PrtGUCmd = new G4UIcmdWithoutParameter("/GammaPrint/Track_UnSet",this);
  PrtGUCmd->SetGuidance("Un sets printing of track gamma results.");

  PhdACmd = new G4UIcmdWithADouble("/Mode2/PHDA",this);
  PhdACmd->SetGuidance("Set pulse height defect parameter A.");

  PhdBCmd = new G4UIcmdWithADouble("/Mode2/PHDB",this);
  PhdBCmd->SetGuidance("Set pulse height defect parameter B.");

}



TrackerGammaSD_Messenger::~TrackerGammaSD_Messenger()

{  

  delete PrtGDir;
  delete PrtGSCmd;
  delete PrtGUCmd;
  delete PhdACmd;
  delete PhdBCmd;

}



void TrackerGammaSD_Messenger::SetNewValue(G4UIcommand* command,
					   G4String newValue)
{ 

  if( command == PrtGSCmd )
    { tracker->SetPrint(); }

  if( command == PrtGUCmd )
    { tracker->UnSetPrint(); }

  if( command == PhdACmd )
    { tracker->SetPHDA( PhdACmd->GetNewDoubleValue(newValue) ); }
  
  if( command == PhdBCmd )
    { tracker->SetPHDB( PhdBCmd->GetNewDoubleValue(newValue) ); }

}

