#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4RunManagerFactory.hh"
#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4SystemOfUnits.hh"

int main(int argc, char** argv)
{
    G4UIExecutive* ui = nullptr;

    if (argc == 1)
        ui = new G4UIExecutive(argc, argv);

    auto* runManager =
        G4RunManagerFactory::CreateRunManager();

    auto* detector =
        new DetectorConstruction();

    runManager->SetUserInitialization(detector);

    auto* physics = new FTFP_BERT();
    physics->SetDefaultCutValue(0.1 * mm);
    physics->RegisterPhysics(new G4OpticalPhysics());

    runManager->SetUserInitialization(physics);

    auto* optical =
        G4OpticalParameters::Instance();

    optical->SetScintFiniteRiseTime(true);
    optical->SetScintTrackSecondariesFirst(true);
    optical->SetCerenkovTrackSecondariesFirst(true);
    optical->SetScintVerboseLevel(0);
    optical->SetBoundaryVerboseLevel(0);

    runManager->SetUserInitialization(
        new ActionInitialization(detector)
    );

    auto* vis = new G4VisExecutive();
    vis->Initialize();

    auto* uiManager =
        G4UImanager::GetUIpointer();

    if (ui) {
        uiManager->ApplyCommand(
            "/control/execute macros/vis.mac"
        );
        ui->SessionStart();
        delete ui;
    }
    else {
        uiManager->ApplyCommand(
            G4String("/control/execute ")
            + argv[1]
        );
    }

    delete vis;
    delete runManager;

    return 0;
}