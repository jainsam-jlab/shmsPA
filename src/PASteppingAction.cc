//copied from:
/// \file optical/OpNovice2/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



#include "PAAnalysis.hh"
#include "PATrackingAction.hh"

#include "G4Box.hh"
#include "G4GeometryTolerance.hh"
#include "G4LogicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

#include <algorithm>

#include "PASteppingAction.hh"
#include "PAEventAction.hh"
#include "G4Run.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4ProcessManager.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PASteppingAction::PASteppingAction(
    PAEventAction* eventAction,
    PATrackingAction* trackingAction)
  : G4UserSteppingAction(),
    fVerbose(0),
    fEventAction(eventAction),
    fTrackingAction(trackingAction)
{
    /*
     * The initializer list stores two non-owning pointers.
     *
     * PAActionInitialization creates all three action objects. Geant4 owns
     * them, so PASteppingAction must not delete either pointer.
     */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PASteppingAction::~PASteppingAction()
{}
namespace
{
    /*
     * Safely return the logical-volume name at a step point.
     *
     * A step leaving the Geant4 world has a null post-step physical volume.
     */
    G4String LogicalVolumeName(const G4StepPoint* point)
    {
        if (point == nullptr) {
            return "";
        }

        const G4VPhysicalVolume* physical =
            point->GetPhysicalVolume();

        if (physical == nullptr) {
            return "";
        }

        const G4LogicalVolume* logical =
            physical->GetLogicalVolume();

        if (logical == nullptr) {
            return "";
        }

        return logical->GetName();
    }

    /*
     * Count Cherenkov photons created during exactly this step.
     *
     * The returned photons remain fully managed and transported by Geant4.
     * This helper only counts them for the existing detector-response output.
     */
    G4int CountCerenkovPhotons(const G4Step* step)
    {
        if (step == nullptr) {
            return 0;
        }

        G4int count = 0;

        const std::vector<const G4Track*>* secondaries =
            step->GetSecondaryInCurrentStep();

        if (secondaries == nullptr) {
            return count;
        }

        const G4ParticleDefinition* opticalPhoton =
            G4OpticalPhoton::OpticalPhotonDefinition();

        for (const G4Track* secondary : *secondaries) {
            if (secondary == nullptr) {
                continue;
            }

            if (secondary->GetParticleDefinition() != opticalPhoton) {
                continue;
            }

            const G4VProcess* creator =
                secondary->GetCreatorProcess();

            if (creator != nullptr &&
                creator->GetProcessName() == "Cerenkov") {
                ++count;
            }
        }

        return count;
    }

    /*
     * Determine whether this step crosses the downstream calorimeter face.
     *
     * Required conditions:
     *
     *   1. The step begins inside CalLogical.
     *   2. The post-step point is a geometry boundary.
     *   3. The step ends outside CalLogical.
     *   4. The boundary point is at local +z of the calorimeter box.
     *   5. The track is moving in the local +z direction.
     *
     * The touchable transform converts a global position and direction into
     * coordinates local to this particular calorimeter placement. Therefore,
     * no global number such as z = 492 cm is hard-coded.
     */
    G4bool CrossedDownstreamCalorimeterFace(const G4Step* step)
    {
        if (step == nullptr) {
            return false;
        }

        const G4StepPoint* prePoint =
            step->GetPreStepPoint();

        const G4StepPoint* postPoint =
            step->GetPostStepPoint();

        if (prePoint == nullptr || postPoint == nullptr) {
            return false;
        }

        const G4VPhysicalVolume* preVolume =
            prePoint->GetPhysicalVolume();

        if (preVolume == nullptr ||
            preVolume->GetLogicalVolume() == nullptr ||
            preVolume->GetLogicalVolume()->GetName() != "CalLogical") {
            return false;
        }

        if (postPoint->GetStepStatus() != fGeomBoundary) {
            return false;
        }

        /*
         * Reject a step that remains in the same calorimeter logical volume.
         * For this geometry the calorimeter has no daughter segmentation, but
         * keeping this check makes the meaning explicit.
         */
        const G4VPhysicalVolume* postVolume =
            postPoint->GetPhysicalVolume();

        if (postVolume != nullptr &&
            postVolume->GetLogicalVolume() != nullptr &&
            postVolume->GetLogicalVolume()->GetName() == "CalLogical") {
            return false;
        }

        /*
         * The calorimeter solid is a G4Box in PADetectorConstruction.cc.
         * dynamic_cast verifies that assumption at runtime. If the geometry is
         * later changed to another solid type, this safely returns false.
         */
        const G4Box* calBox =
            dynamic_cast<const G4Box*>(
                preVolume->GetLogicalVolume()->GetSolid());

        if (calBox == nullptr) {
            return false;
        }

        const G4TouchableHandle& touchable =
            prePoint->GetTouchableHandle();

        if (!touchable) {
            return false;
        }

        const G4TouchableHistory* history =
            dynamic_cast<const G4TouchableHistory*>(
                touchable.operator->());

        if (history == nullptr) {
            return false;
        }

        /*
         * GetTopTransform() converts global coordinates to coordinates local
         * to the current top-level touchable volume.
         */
        const G4AffineTransform& globalToLocal =
            history->GetHistory()->GetTopTransform();

        const G4ThreeVector localBoundaryPoint =
            globalToLocal.TransformPoint(postPoint->GetPosition());

        const G4ThreeVector localDirection =
            globalToLocal.TransformAxis(
                postPoint->GetMomentumDirection());

        const G4double halfLengthZ =
            calBox->GetZHalfLength();

        /*
         * Navigation points may differ from the mathematical surface by a
         * tiny floating-point tolerance. Use Geant4's geometry tolerance
         * rather than demanding exact equality.
         */
        const G4double surfaceTolerance =
            G4GeometryTolerance::GetInstance()
                ->GetSurfaceTolerance();

        const G4double comparisonTolerance =
            std::max(10.0 * surfaceTolerance, 1.0e-9 * mm);

        const G4bool isAtPositiveZFace =
            localBoundaryPoint.z() >=
            halfLengthZ - comparisonTolerance;

        const G4bool movingDownstream =
            localDirection.z() > 0.0;

        return isAtPositiveZFace && movingDownstream;
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PASteppingAction::UserSteppingAction(const G4Step* step)
{
    /*
     * Geant4 calls this function once for every step taken by every particle.
     *
     * This function performs two kinds of work:
     *
     *   1. Update detector-crossing flags for the current track.
     *   2. Preserve the existing event-level detector response calculations.
     *
     * It does not write a Tracks row. PATrackingAction writes that row after
     * the complete track has finished.
     */
    if (step == nullptr) {
        return;
    }

    /*
     * GetTrack() returns a pointer to the Geant4-owned track taking this step.
     *
     * We inspect the track, but do not delete it.
     */
    G4Track* track = step->GetTrack();

    if (track == nullptr) {
        return;
    }

    /*
     * Every step has:
     *
     *   prePoint  = state immediately before the step
     *   postPoint = state immediately after the step
     */
    const G4StepPoint* prePoint =
        step->GetPreStepPoint();

    const G4StepPoint* postPoint =
        step->GetPostStepPoint();

    if (prePoint == nullptr || postPoint == nullptr) {
        return;
    }

    const G4ParticleDefinition* particle =
        track->GetParticleDefinition();

    if (particle == nullptr) {
        return;
    }

    /*
     * Compare the particle-definition pointer with Geant4's optical-photon
     * definition. Optical photons remain transported, but are not recorded
     * in the Tracks tree.
     */
    const G4bool isOpticalPhoton =
        particle == G4OpticalPhoton::OpticalPhotonDefinition();

    /*
     * These helpers safely return an empty string when a physical or logical
     * volume is null. A null post-step volume is possible when a track leaves
     * the Geant4 world.
     */
    const G4String preLogicalName =
        LogicalVolumeName(prePoint);

    const G4String postLogicalName =
        LogicalVolumeName(postPoint);

    /*
     * Update the truth flags only if PATrackingAction confirms that its
     * temporary record belongs to this exact TrackID.
     *
     * The check also excludes optical photons because PATrackingAction does
     * not create a truth record for them.
     */
    if (fTrackingAction != nullptr &&
        fTrackingAction->IsRecordingTrack(track->GetTrackID())) {

        const G4int trackID = track->GetTrackID();

        /*
         * Checking both the pre-step and post-step volumes handles:
         *
         *   - a track crossing into a detector;
         *   - a secondary particle created inside a detector.
         */

        // S1 consists of separate X and Y logical volumes.
        if (preLogicalName == "S1XLogical" ||
            preLogicalName == "S1YLogical" ||
            postLogicalName == "S1XLogical" ||
            postLogicalName == "S1YLogical") {
            fTrackingAction->MarkReachedS1(trackID);
        }

        // The heavy-gas Cherenkov detector envelope.
        if (preLogicalName == "HGCLogical" ||
            postLogicalName == "HGCLogical") {
            fTrackingAction->MarkReachedHGC(trackID);
        }

        /*
         * ReachedAGC means reaching the aerogel detector.
         *
         * AGCLogical is the detector envelope.
         * AGCTrayLogical is the aerogel radiator inside that envelope.
         */
        if (preLogicalName == "AGCLogical" ||
            preLogicalName == "AGCTrayLogical" ||
            postLogicalName == "AGCLogical" ||
            postLogicalName == "AGCTrayLogical") {
            fTrackingAction->MarkReachedAGC(trackID);
        }

        // S2 consists of separate X and Y logical volumes.
        if (preLogicalName == "S2XLogical" ||
            preLogicalName == "S2YLogical" ||
            postLogicalName == "S2XLogical" ||
            postLogicalName == "S2YLogical") {
            fTrackingAction->MarkReachedS2(trackID);
        }

        /*
         * A track has reached the calorimeter if:
         *
         *   - it begins this step inside CalLogical, or
         *   - it ends this step inside CalLogical.
         *
         * This also marks a secondary born inside the calorimeter.
         */
        if (preLogicalName == "CalLogical" ||
            postLogicalName == "CalLogical") {
            fTrackingAction->MarkReachedCal(trackID);
        }

        /*
         * This helper performs the more restrictive downstream-face test:
         *
         *   - step begins in CalLogical;
         *   - step ends at a geometry boundary;
         *   - step leaves CalLogical;
         *   - boundary is at calorimeter-local +z;
         *   - motion is in the local +z direction.
         */
        if (CrossedDownstreamCalorimeterFace(step)) {
            fTrackingAction->MarkExitedCal(trackID);
        }
    }

    /*
     * Preserve the original CopyNo behavior.
     *
     * If a primary particle is stopped outside CalLogical, save the copy
     * number of the volume containing its final pre-step point.
     */
    if (track->GetTrackStatus() == fStopAndKill &&
        track->GetParentID() == 0 &&
        preLogicalName != "CalLogical") {

        const G4VPhysicalVolume* preVolume =
            prePoint->GetPhysicalVolume();

        if (preVolume != nullptr) {
            G4AnalysisManager::Instance()->FillNtupleIColumn(
                PAAnalysis::kPANtuple,
                PAAnalysis::PA::kCopyNo,
                preVolume->GetCopyNo());
        }
    }

    /*
     * Preserve the original Cherenkov photon counting.
     *
     * CountCerenkovPhotons(step) examines the secondary particles produced
     * during this particular step and counts optical photons whose creator
     * process is "Cerenkov".
     */
    if (!isOpticalPhoton && fEventAction != nullptr) {
        if (preLogicalName == "S2YLogical") {
            fEventAction->AddS2YNPE(
                CountCerenkovPhotons(step));
        }

        if (preLogicalName == "AGCTrayLogical") {
            fEventAction->AddAGCNPE(
                CountCerenkovPhotons(step));
        }

        if (preLogicalName == "HGCLogical") {
            fEventAction->AddHGCNPE(
                CountCerenkovPhotons(step));
        }

        if (preLogicalName == "NGCLogical") {
            fEventAction->AddNGCNPE(
                CountCerenkovPhotons(step));
        }
    }

    /*
     * Preserve calorimeter energy accumulation.
     *
     * GetTotalEnergyDeposit() returns the energy deposited by this step.
     * PAEventAction adds it to the total for the complete event.
     */
    if (preLogicalName == "CalLogical" &&
        fEventAction != nullptr) {
        fEventAction->AddCalEnergy(
            step->GetTotalEnergyDeposit());
    }
}  


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
