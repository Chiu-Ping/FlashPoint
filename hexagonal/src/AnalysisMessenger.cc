#include "AnalysisMessenger.hh"
#include "AnalysisConfig.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

AnalysisMessenger::AnalysisMessenger() {
  fOutputModeCmd = new G4UIcmdWithAString("/analysis/outputMode", this);

  fOutputModeCmd->SetGuidance("Set ROOT output mode: full or light");
  fOutputModeCmd->SetCandidates("full light");
  fOutputModeCmd->SetDefaultValue("full");
}

AnalysisMessenger::~AnalysisMessenger() { delete fOutputModeCmd; }

void AnalysisMessenger::SetNewValue(G4UIcommand *command, G4String value) {
  if (command == fOutputModeCmd) {

    if (value == "light")
      AnalysisConfig::SetMode(AnalysisConfig::Mode::Light);

    else
      AnalysisConfig::SetMode(AnalysisConfig::Mode::Full);
  }
}
