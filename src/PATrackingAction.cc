/*
 * PATrackingAction.cc
 * --------------------------------------------------------------------------
 * Event/track flow:
 *
 * EVENT
 *   |
 *   +-- Geant4 begins Track 1
 *   |      |
 *   |      +-- PreUserTrackingAction()
 *   |      |      Create a record, or recover a suspended track's record.
 *   |      |      Filter optical photons from truth output.
 *   |      |      Save starting position, momentum, volume, IDs and process.
 *   |      |
 *   |      +-- Geant4 transports many G4Step objects
 *   |      |      PASteppingAction sets detector-crossing flags.
 *   |      |
 *   |      +-- PostUserTrackingAction()
 *   |             Read final position, momentum, process and volume.
 *   |             Convert Geant4 units to cm and GeV.
 *   |             Write exactly one row into Tracks.
 *   |
 *   +-- Geant4 begins the next track
 *          ...
 *
 * The G4Track object belongs to Geant4. This class only reads it and never
 * deletes it.
 * --------------------------------------------------------------------------
 */

#include "PATrackingAction.hh"

#include "PAAnalysis.hh"
#include "PAEventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

namespace
{
    /*
     * Return a physical-volume name without dereferencing a null pointer.
     *
     * A null volume can occur when a particle is outside the Geant4 world.
     */
    G4String GetVolumeName(const G4VPhysicalVolume* volume,
                           const G4String& nullName)
    {
        if (volume == nullptr) {
            return nullName;
        }

        return volume->GetName();
    }

    /*
     * Determine the process associated with the final step.
     *
     * GetProcessDefinedStep() may return nullptr. Geometry/world-boundary
     * steps do not always have an ordinary physics process associated with
     * them, so every pointer is checked before it is dereferenced.
     */
    G4String GetEndProcessName(const G4Track* track)
    {
        if (track == nullptr) {
            return "Unknown";
        }

        const G4Step* finalStep = track->GetStep();
        if (finalStep == nullptr) {
            return "None";
        }

        const G4StepPoint* postPoint = finalStep->GetPostStepPoint();
        if (postPoint == nullptr) {
            return "None";
        }

        const G4VProcess* process = postPoint->GetProcessDefinedStep();
        if (process != nullptr) {
            return process->GetProcessName();
        }

        /*
         * A world-boundary crossing has no volume on the far side and may
         * have no defining physics process. Give it a meaningful name.
         */
        if (postPoint->GetStepStatus() == fWorldBoundary) {
            return "WorldBoundary";
        }

        return "None";
    }

    /*
     * Determine the final physical volume.
     *
     * For a track killed inside the detector, the post-step physical volume
     * is normally the containing volume. For a track leaving the world, it is
     * null and is represented by the string "OutOfWorld".
     */
    G4String GetEndVolumeName(const G4Track* track)
    {
        if (track == nullptr) {
            return "Unknown";
        }

        const G4Step* finalStep = track->GetStep();
        if (finalStep == nullptr) {
            return GetVolumeName(track->GetVolume(), "OutOfWorld");
        }

        const G4StepPoint* postPoint = finalStep->GetPostStepPoint();
        if (postPoint == nullptr) {
            return GetVolumeName(track->GetVolume(), "OutOfWorld");
        }

        const G4VPhysicalVolume* postVolume =
            postPoint->GetPhysicalVolume();

        if (postVolume != nullptr) {
            return postVolume->GetName();
        }

        if (postPoint->GetStepStatus() == fWorldBoundary) {
            return "OutOfWorld";
        }

        return "None";
    }
}

PATrackingAction::PATrackingAction(PAEventAction* eventAction)
  : G4UserTrackingAction(),
    fEventAction(eventAction)
{
    /*
     * The initializer list above initializes the base class first and then
     * stores the non-owning PAEventAction pointer in fEventAction.
     */
}

PATrackingAction::~PATrackingAction()
{
    /*
     * Nothing is deleted here:
     *
     *   - Geant4 owns the G4Track.
     *   - Geant4 owns registered user-action objects.
     *   - fEventAction is only a non-owning pointer.
     */
}

void PATrackingAction::ResetRecord(TrackRecord& record) const
{
    record.saveThisTrack = false;

    record.eventID = -1;
    record.trackID = -1;
    record.parentID = -1;
    record.pdg = 0;

    record.creatorProcess = "Unknown";


    record.reachedS1 = false;
    record.reachedHGC = false;
    record.reachedAGC = false;
    record.reachedS2 = false;
    record.reachedCal = false;
    record.exitedCal = false;
}

void PATrackingAction::PreUserTrackingAction(const G4Track* track)
{
    if (track == nullptr) {
        return;
    }

    const G4ParticleDefinition* particle =
        track->GetParticleDefinition();

    if (particle == nullptr) {
        return;
    }

    /*
     * Compare particle-definition pointers rather than spelling the particle
     * name. Geant4 maintains one definition object for optical photons.
     *
     * We do not kill the photon. We only decline to write it to Tracks.
     * Its Cherenkov response is still processed by PASteppingAction.
     */
    if (particle == G4OpticalPhoton::OpticalPhotonDefinition()) {
        return;
    }

    /*
     * TrackID is unique only within one event. EventID and TrackID together
     * form the unique key used by later ROOT/Python analysis.
     */
    const G4Event* event =
        G4RunManager::GetRunManager()->GetCurrentEvent();

    const G4int eventID = event != nullptr ? event->GetEventID() : -1;
    const G4int trackID = track->GetTrackID();

    /*
     * A suspended track is later presented to PreUserTrackingAction again.
     * If its record already exists for this event, preserve the original
     * starting state and all detector flags collected before suspension.
     */
    const auto existing = fTrackRecords.find(trackID);
    if (existing != fTrackRecords.end() &&
        existing->second.eventID == eventID) {
        return;
    }

    TrackRecord record;
    ResetRecord(record);
    record.saveThisTrack = true;
    record.eventID = eventID;

    /*
     * The -> operator accesses a member function through a pointer.
     */
    record.trackID = trackID;
    record.parentID = track->GetParentID();
    record.pdg = particle->GetPDGEncoding();

    /*
     * Primary tracks do not have a creator process, because they were created
     * by the event generator rather than by a Geant4 physics process.
     */
    const G4VProcess* creator = track->GetCreatorProcess();

    if (creator != nullptr) {
        record.creatorProcess = creator->GetProcessName();
    }
    else if (track->GetParentID() == 0) {
        record.creatorProcess = "Primary";
    }
    else {
        record.creatorProcess = "Unknown";
    }


    fTrackRecords[trackID] = record;
}

void PATrackingAction::PostUserTrackingAction(const G4Track* track)
{
    if (track == nullptr) {
        return;
    }

    /*
     * Optical photons and invalid tracks were marked as not recordable in
     * PreUserTrackingAction().
     */
    const G4int trackID = track->GetTrackID();
    const auto found = fTrackRecords.find(trackID);
    if (found == fTrackRecords.end()) {
        return;
    }

    /*
     * fSuspend and fPostponeToNextEvent are not true track endings. Geant4
     * will resume these tracks, so keep their records and do not write rows.
     * fAlive is also treated as non-terminal defensively.
     */
    const G4TrackStatus status = track->GetTrackStatus();
    if (status == fSuspend ||
        status == fPostponeToNextEvent ||
        status == fAlive) {
        return;
    }

    TrackRecord& record = found->second;



    const G4String endVolume = GetEndVolumeName(track);
    const G4String endProcess = GetEndProcessName(track);

    G4AnalysisManager* analysisManager =
        G4AnalysisManager::Instance();

    /*
     * Every fill call names both:
     *
     *   1. the ntuple ID, and
     *   2. the column ID.
     *
     * This is essential now that the output contains PA and Tracks.
     */
    analysisManager->FillNtupleIColumn(
        PAAnalysis::kTracksNtuple,
        PAAnalysis::Tracks::kEventID,
        record.eventID);

    analysisManager->FillNtupleIColumn(
        PAAnalysis::kTracksNtuple,
        PAAnalysis::Tracks::kTrackID,
        record.trackID);

    analysisManager->FillNtupleIColumn(
        PAAnalysis::kTracksNtuple,
        PAAnalysis::Tracks::kParentID,
        record.parentID);

    analysisManager->FillNtupleIColumn(
        PAAnalysis::kTracksNtuple,
        PAAnalysis::Tracks::kPDG,
        record.pdg);

    analysisManager->FillNtupleSColumn(
        PAAnalysis::kTracksNtuple,
        PAAnalysis::Tracks::kCreatorProcess,
        record.creatorProcess);


    analysisManager->FillNtupleSColumn(
        PAAnalysis::kTracksNtuple,
        PAAnalysis::Tracks::kEndVolume,
        endVolume);

    analysisManager->FillNtupleSColumn(
        PAAnalysis::kTracksNtuple,
        PAAnalysis::Tracks::kEndProcess,
        endProcess);



    analysisManager->FillNtupleIColumn(
        PAAnalysis::kTracksNtuple,
        PAAnalysis::Tracks::kReachedS1,
        record.reachedS1 ? 1 : 0);

    analysisManager->FillNtupleIColumn(
        PAAnalysis::kTracksNtuple,
        PAAnalysis::Tracks::kReachedHGC,
        record.reachedHGC ? 1 : 0);

    analysisManager->FillNtupleIColumn(
        PAAnalysis::kTracksNtuple,
        PAAnalysis::Tracks::kReachedAGC,
        record.reachedAGC ? 1 : 0);

    analysisManager->FillNtupleIColumn(
        PAAnalysis::kTracksNtuple,
        PAAnalysis::Tracks::kReachedS2,
        record.reachedS2 ? 1 : 0);

    analysisManager->FillNtupleIColumn(
        PAAnalysis::kTracksNtuple,
        PAAnalysis::Tracks::kReachedCal,
        record.reachedCal ? 1 : 0);

    analysisManager->FillNtupleIColumn(
        PAAnalysis::kTracksNtuple,
        PAAnalysis::Tracks::kExitedCal,
        record.exitedCal ? 1 : 0);

    /*
     * All columns above belong to one completed track. AddNtupleRow commits
     * them as one row in the Tracks tree.
     */
    analysisManager->AddNtupleRow(PAAnalysis::kTracksNtuple);

    /*
     * Update the event-level summary only for a primary positive pion.
     *
     * ParentID == 0 means primary.
     * PDG == 211 means pi+.
     *
     * We do not rely only on TrackID == 1.
     */
    if (record.parentID == 0 &&
        record.pdg == 211 &&
        fEventAction != nullptr) {
        fEventAction->SetPrimaryPionSummary(
            record.pdg,
            record.reachedCal,
            record.exitedCal);
    }

    // The terminal row is complete; this TrackID no longer needs state.
    fTrackRecords.erase(found);
}

G4bool PATrackingAction::IsRecordingTrack(G4int trackID) const
{
    return fTrackRecords.find(trackID) != fTrackRecords.end();
}

void PATrackingAction::MarkReachedS1(G4int trackID)
{
    const auto found = fTrackRecords.find(trackID);
    if (found != fTrackRecords.end()) {
        found->second.reachedS1 = true;
    }
}

void PATrackingAction::MarkReachedHGC(G4int trackID)
{
    const auto found = fTrackRecords.find(trackID);
    if (found != fTrackRecords.end()) {
        found->second.reachedHGC = true;
    }
}

void PATrackingAction::MarkReachedAGC(G4int trackID)
{
    const auto found = fTrackRecords.find(trackID);
    if (found != fTrackRecords.end()) {
        found->second.reachedAGC = true;
    }
}

void PATrackingAction::MarkReachedS2(G4int trackID)
{
    const auto found = fTrackRecords.find(trackID);
    if (found != fTrackRecords.end()) {
        found->second.reachedS2 = true;
    }
}

void PATrackingAction::MarkReachedCal(G4int trackID)
{
    const auto found = fTrackRecords.find(trackID);
    if (found != fTrackRecords.end()) {
        found->second.reachedCal = true;
    }
}

void PATrackingAction::MarkExitedCal(G4int trackID)
{
    const auto found = fTrackRecords.find(trackID);
    if (found != fTrackRecords.end()) {
        /*
         * Enforce the physical implication in the stored data:
         *
         *     ExitedCal == true  implies  ReachedCal == true
         */
        found->second.reachedCal = true;
        found->second.exitedCal = true;
    }
}
