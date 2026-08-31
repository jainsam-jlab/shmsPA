#include "PARunAction.hh"
#include "PAAnalysis.hh"
#include "PAEventAction.hh"
#include "G4AnalysisManager.hh"

#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PARunAction::PARunAction(PAEventAction* eventAction)
 : G4UserRunAction(),
   fEventAction(eventAction)
{ 
  // Create analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Default settings
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName("PAOutData.root"); //default, can overwrite in macro

  // Creating ntuple
  //if ( fEventAction ) {
  /*
   * Define the ROOT schemas.
   *
   * CreateNtuple() creates a tree-like ntuple and returns its integer ID.
   * CreateNtuple...Column() appends a typed column.
   * FinishNtuple() completes that ntuple's schema.
   *
   * PA is created first, so its ID is 0.
   * Tracks is created second, so its ID is 1.
   */
  if (fEventAction != nullptr) {
    // ------------------------------------------------------------------
    // Existing event-level PA tree
    // ------------------------------------------------------------------
    analysisManager->CreateNtuple("PA", "Hits and event summary");

    analysisManager->CreateNtupleDColumn(
        PAAnalysis::kPANtuple, "S1XTime");
    analysisManager->CreateNtupleDColumn(
        PAAnalysis::kPANtuple, "S1YTime");
    analysisManager->CreateNtupleDColumn(
        PAAnalysis::kPANtuple, "S2XTime");
    analysisManager->CreateNtupleDColumn(
        PAAnalysis::kPANtuple, "S2YTime");

    analysisManager->CreateNtupleDColumn(
        PAAnalysis::kPANtuple, "S1XEnergy");
    analysisManager->CreateNtupleDColumn(
        PAAnalysis::kPANtuple, "S1YEnergy");
    analysisManager->CreateNtupleDColumn(
        PAAnalysis::kPANtuple, "S2XEnergy");
    analysisManager->CreateNtupleDColumn(
        PAAnalysis::kPANtuple, "S2YEnergy");

    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kPANtuple, "S2YNPE");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kPANtuple, "CopyNo");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kPANtuple, "AGCNPE");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kPANtuple, "HGCNPE");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kPANtuple, "NGCNPE");
    analysisManager->CreateNtupleDColumn(
        PAAnalysis::kPANtuple, "CalEnergy");

    // New columns are appended after the original 14 columns.
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kPANtuple, "EventID");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kPANtuple, "PrimaryPDG");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kPANtuple, "PrimaryReachedCal");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kPANtuple, "PrimaryExitedCal");

    analysisManager->FinishNtuple(
        PAAnalysis::kPANtuple);

    // ------------------------------------------------------------------
    // New one-row-per-track truth tree
    // ------------------------------------------------------------------
    analysisManager->CreateNtuple(
        "Tracks",
        "One row per non-optical Geant4 track");

    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kTracksNtuple, "EventID");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kTracksNtuple, "TrackID");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kTracksNtuple, "ParentID");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kTracksNtuple, "PDG");
    analysisManager->CreateNtupleSColumn(
        PAAnalysis::kTracksNtuple, "CreatorProcess");

    analysisManager->CreateNtupleSColumn(
        PAAnalysis::kTracksNtuple, "EndVolume");
    analysisManager->CreateNtupleSColumn(
        PAAnalysis::kTracksNtuple, "EndProcess");

    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kTracksNtuple, "ReachedS1");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kTracksNtuple, "ReachedHGC");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kTracksNtuple, "ReachedAGC");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kTracksNtuple, "ReachedS2");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kTracksNtuple, "ReachedCal");
    analysisManager->CreateNtupleIColumn(
        PAAnalysis::kTracksNtuple, "ExitedCal");

    analysisManager->FinishNtuple(
        PAAnalysis::kTracksNtuple);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PARunAction::~PARunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PARunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file 
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PARunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // save histograms & ntuple
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
