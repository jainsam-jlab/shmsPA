#ifndef PATrackingAction_h
#define PATrackingAction_h 1

/*
 * PATrackingAction.hh
 * --------------------------------------------------------------------------
 * PATrackingAction records one summary row for each physical Geant4 track.
 *
 * A "track" is the complete history of one transported particle. A track may
 * contain many individual G4Step objects.
 *
 * Geant4 calls the tracking action in this order:
 *
 *     EVENT
 *       |
 *       +-- begin Track 1
 *       |      |
 *       |      +-- PreUserTrackingAction()
 *       |      |      Save the track's initial state.
 *       |      |
 *       |      +-- many calls to PASteppingAction
 *       |      |      Set detector-crossing flags.
 *       |      |
 *       |      +-- PostUserTrackingAction()
 *       |             Save the final state and write one Tracks row.
 *       |
 *       +-- begin Track 2
 *              ...
 *
 * The class owns a map of TrackRecord objects indexed by TrackID. This is
 * necessary because optical processes may suspend a charged parent, transport
 * its photons, and later resume the same parent track.
 *
 * Optical photons are still transported normally, but are excluded from the
 * Tracks ROOT tree to prevent extremely large output files.
 * --------------------------------------------------------------------------
 */

#include "G4UserTrackingAction.hh"
#include "globals.hh"

#include <unordered_map>

class G4Track;
class PAEventAction;

class PATrackingAction : public G4UserTrackingAction
{
  public:
    /*
     * eventAction is a non-owning pointer.
     *
     * PAActionInitialization creates both action objects. Geant4 owns and
     * eventually deletes registered user actions. PATrackingAction must not
     * delete eventAction.
     */
    explicit PATrackingAction(PAEventAction* eventAction);

    virtual ~PATrackingAction();

    /*
     * "override" asks the compiler to verify that these functions really
     * override virtual functions declared by G4UserTrackingAction.
     *
     * Geant4 calls them automatically. User code does not call them.
     *
     * "const G4Track*" means:
     *   - this is a pointer to a G4Track owned by Geant4;
     *   - we may inspect the track but may not modify it through this pointer;
     *   - we must not delete it.
     */
    void PreUserTrackingAction(const G4Track* track) override;
    void PostUserTrackingAction(const G4Track* track) override;

    /*
     * PASteppingAction calls these functions while the current track is being
     * transported. A flag can change from false to true, but never back to
     * false during that track.
     */
    void MarkReachedS1(G4int trackID);
    void MarkReachedHGC(G4int trackID);
    void MarkReachedAGC(G4int trackID);
    void MarkReachedS2(G4int trackID);
    void MarkReachedCal(G4int trackID);
    void MarkExitedCal(G4int trackID);

    /*
     * This check prevents a stepping callback from accidentally changing the
     * record for a different track. It also returns false for optical photons,
     * because they are deliberately excluded from the Tracks tree.
     */
    G4bool IsRecordingTrack(G4int trackID) const;

  private:
    /*
     * A simple structure is used instead of a separate heap-allocated class.
     *
     * It contains the information that must survive between the pre-tracking
     * and post-tracking callbacks, including any suspension/resume cycles.
     * Final quantities are read directly when the track genuinely terminates.
     */
    struct TrackRecord
    {
        G4bool saveThisTrack;

        G4int eventID;
        G4int trackID;
        G4int parentID;
        G4int pdg;

        G4String creatorProcess;
        G4String startVolume;

        G4double startX;
        G4double startY;
        G4double startZ;

        G4double startPx;
        G4double startPy;
        G4double startPz;
        G4double startKE;

        G4bool reachedS1;
        G4bool reachedHGC;
        G4bool reachedAGC;
        G4bool reachedS2;
        G4bool reachedCal;
        G4bool exitedCal;
    };

    // Initialize every field when Geant4 first begins a new track.
    void ResetRecord(TrackRecord& record) const;

    /*
     * fEventAction is non-owning. It is used only to pass the completed
     * primary-pion summary to the event-level PA tree.
     */
    PAEventAction* fEventAction;

    /*
     * Records indexed by TrackID.
     *
     * Usually Geant4 transports one track continuously. Optical processes,
     * however, may suspend a charged parent, transport its photons, and then
     * resume the parent. The map keeps the parent's original start state and
     * accumulated detector flags across those suspension/resume cycles.
     */
    std::unordered_map<G4int, TrackRecord> fTrackRecords;
};

#endif
