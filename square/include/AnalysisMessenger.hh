#pragma once

#include "G4UImessenger.hh"

class G4UIcmdWithAString;

class AnalysisMessenger : public G4UImessenger {
public:
    AnalysisMessenger();
    ~AnalysisMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

private:
    G4UIcmdWithAString* fOutputModeCmd = nullptr;
};
