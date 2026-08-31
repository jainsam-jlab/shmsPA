// copied from
/// \file optical/OpNovice2/include/SteppingAction.hh
/// \brief Definition of the SteppingAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PASteppingAction_h
#define PASteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Track.hh"
#include "PAEventAction.hh"
#include "globals.hh"

class PATrackingAction;

class PASteppingAction : public G4UserSteppingAction
{
  public:
        /*
     * Both pointers are non-owning. Geant4 owns the registered action
     * objects and keeps them alive for the run.
     */
    PASteppingAction(PAEventAction* eventAction,
                     PATrackingAction* trackingAction);

    virtual ~PASteppingAction();

    // Geant4 calls this once for every step of every transported track.
    void UserSteppingAction(const G4Step* step) override;
    //void SetVerbose(G4int verbose){fVerbose = verbose;};
    /**
     * virtual void SetNPE(G4int npe){NPE=npe;}
    virtual void SetAGCNPE(G4int npe){AGC_NPE=npe;}
    virtual void SetHGCNPE(G4int npe){HGC_NPE=npe;}
    virtual void SetNGCNPE(G4int npe){NGC_NPE=npe;}
    G4int GetNPE(){return NPE;} 
    G4int GetAGCNPE(){return AGC_NPE;} 
    G4int GetHGCNPE(){return NGC_NPE;} 
    G4int GetNGCNPE(){return HGC_NPE;} 
    **/

  private:
    G4int fVerbose;

    // Non-owning link to event-level detector-response accumulators.
    PAEventAction* fEventAction;

    /*
     * Non-owning link to the current track's temporary truth record.
     * PASteppingAction sets flags; PATrackingAction writes the final row.
     */
    PATrackingAction* fTrackingAction;
    //G4int NGC_NPE;
    //G4int HGC_NPE;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
