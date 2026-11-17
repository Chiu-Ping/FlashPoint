#include "ActionInitialization.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "StackingAction.hh"
#include "SteppingAction.hh"


ActionInitialization::ActionInitialization(
    const DetectorConstruction* detector
)
    : fDetector(detector)
{
}


void
ActionInitialization::BuildForMaster() const
{
    SetUserAction(
        new RunAction(
            fDetector
        )
    );
}


void
ActionInitialization::Build() const
{
    auto* generator =
        new PrimaryGeneratorAction(
            fDetector
        );

    auto* eventAction =
        new EventAction(
            generator,
            fDetector->GetTotalSiPMs(),
            fDetector->GetSiPMGain()
        );


    SetUserAction(
        generator
    );

    SetUserAction(
        new RunAction(
            fDetector
        )
    );

    SetUserAction(
        eventAction
    );

    SetUserAction(
        new StackingAction(
            eventAction
        )
    );

    SetUserAction(
        new SteppingAction(
            generator,
            eventAction
        )
    );
}
