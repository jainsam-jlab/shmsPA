#ifndef PAAnalysis_h
#define PAAnalysis_h 1

/*
 * PAAnalysis.hh
 * --------------------------------------------------------------------------
 * This header assigns readable names to every ROOT ntuple and column number.
 *
 * G4AnalysisManager identifies ntuples and columns using integer IDs. For
 * example, the first ntuple normally has ID 0 and its first column has ID 0.
 *
 * Raw integer calls such as:
 *
 *     FillNtupleIColumn(1, 3, pdg);
 *
 * are difficult to understand and easy to misuse. The named constants below
 * allow us to write:
 *
 *     FillNtupleIColumn(
 *         PAAnalysis::kTracksNtuple,
 *         PAAnalysis::Tracks::kPDG,
 *         pdg);
 *
 * The compiler still passes integers to Geant4, but a human can see exactly
 * which tree and column are being filled.
 *
 * IMPORTANT:
 * The order of each enum must match the order of CreateNtuple...Column()
 * calls in PARunAction.cc.
 * --------------------------------------------------------------------------
 */

#include "globals.hh"

namespace PAAnalysis
{
    /*
     * Ntuple IDs
     * ----------
     * CreateNtuple() assigns IDs in creation order.
     *
     * PARunAction creates:
     *   1. PA     -> ID 0
     *   2. Tracks -> ID 1
     */
    enum NtupleId
    {
        kPANtuple     = 0,
        kTracksNtuple = 1
    };

    /*
     * Columns in the existing PA event tree.
     *
     * The original columns retain IDs 0 through 13. The new columns are
     * appended so existing scripts that depend on the old ordering continue
     * to work.
     */
    namespace PA
    {
        enum ColumnId
        {
            kS1XTime = 0,
            kS1YTime,
            kS2XTime,
            kS2YTime,

            kS1XEnergy,
            kS1YEnergy,
            kS2XEnergy,
            kS2YEnergy,

            kS2YNPE,
            kCopyNo,
            kAGCNPE,
            kHGCNPE,
            kNGCNPE,
            kCalEnergy,

            // New event-level truth columns begin here.
            kEventID,
            kPrimaryPDG,
            kPrimaryReachedCal,
            kPrimaryExitedCal
        };
    }

    /*
     * Columns in the new Tracks tree.
     *
     * One completed, non-optical Geant4 track produces one row.
     */
    namespace Tracks
    {
        enum ColumnId
        {
            kEventID = 0,
            kTrackID,
            kParentID,
            kPDG,
            kCreatorProcess,

            kStartX,
            kStartY,
            kStartZ,

            kStartPx,
            kStartPy,
            kStartPz,
            kStartKE,

            kStartVolume,

            kEndX,
            kEndY,
            kEndZ,

            kEndPx,
            kEndPy,
            kEndPz,
            kEndKE,

            kEndVolume,
            kEndProcess,

            kTrackLength,

            kReachedS1,
            kReachedHGC,
            kReachedAGC,
            kReachedS2,
            kReachedCal,
            kExitedCal
        };
    }
}

#endif