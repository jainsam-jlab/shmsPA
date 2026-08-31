#include "PAActionInitialization.hh"
#include "PAPrimaryGeneratorAction.hh"
#include "PARunAction.hh"
#include "PAEventAction.hh"
#include "PASteppingAction.hh"
#include "PATrackingAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PAActionInitialization::PAActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PAActionInitialization::~PAActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PAActionInitialization::BuildForMaster() const
{
  PAEventAction* eventAction = new PAEventAction;
  SetUserAction(new PARunAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PAActionInitialization::Build() const
{
  SetUserAction(new PAPrimaryGeneratorAction);

  /*
   * Geant4 takes ownership of every action passed to SetUserAction().
   *
   * The raw pointers below are retained only so the action objects can
   * communicate. They must not be deleted manually by another action.
   */
  auto eventAction = new PAEventAction;
  SetUserAction(eventAction);

  /*
   * PATrackingAction owns the temporary record for the one track currently
   * being transported. It also has a non-owning pointer to eventAction so it
   * can report the primary-pion summary.
   */
  auto trackingAction =
      new PATrackingAction(eventAction);
  SetUserAction(trackingAction);

  /*
   * PASteppingAction communicates:
   *
   *   - event detector response to eventAction;
   *   - current-track crossing flags to trackingAction.
   */
  auto steppingAction =
      new PASteppingAction(eventAction, trackingAction);
  SetUserAction(steppingAction);

  SetUserAction(new PARunAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
