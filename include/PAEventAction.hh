#ifndef PAEventAction_h
#define PAEventAction_h 1


#include "G4UserEventAction.hh"
#include "globals.hh"

#include <vector>
#include <array>

// named constants
const G4int kDim = 4; //total number SDs

/// Event action

class PAEventAction : public G4UserEventAction
{
public:
    PAEventAction();
    virtual ~PAEventAction();
    
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    virtual void AddS2YNPE(G4int npe){S2YNPE+=npe;}
    virtual void AddAGCNPE(G4int npe){AGCNPE+=npe;}
    virtual void AddHGCNPE(G4int npe){HGCNPE+=npe;}
    virtual void AddNGCNPE(G4int npe){NGCNPE+=npe;}
        virtual void AddCalEnergy(G4double edep){CalEnergy+=edep;}

    /*
     * PATrackingAction calls this when the primary pi+ finishes.
     *
     * The values remain stored until EndOfEventAction() writes the one PA row.
     */
    void SetPrimaryPionSummary(G4int pdg,
                               G4bool reachedCal,
                               G4bool exitedCal)
    {
        PrimaryPDG = pdg;
        PrimaryReachedCal = reachedCal;
        PrimaryExitedCal = exitedCal;
    }

private:
    // hit collections Ids
    std::array<G4int, kDim> fHodHCID;
    G4int S2YNPE;
    G4int AGCNPE;
    G4int HGCNPE;
    G4int NGCNPE;
        G4double CalEnergy;

    /*
     * Event-level primary-pion truth.
     *
     * These are reset at the start of every event and filled when the primary
     * pi+ finishes transport.
     */
    G4int PrimaryPDG;
    G4bool PrimaryReachedCal;
    G4bool PrimaryExitedCal;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
