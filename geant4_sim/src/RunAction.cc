#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>
#include <iomanip>
G4String RunAction::fsRunTag = "";
RunAction::RunAction()
{
  DefineMessenger();
  auto* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetNtupleMerging(true);
  analysisManager->CreateNtuple("PionEvents", "Per-event pion path/energy data");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->CreateNtupleDColumn("labThetaDeg");
  analysisManager->CreateNtupleDColumn("siliconPathLength_mm");
  analysisManager->CreateNtupleDColumn("totalEdep_MeV");
  analysisManager->CreateNtupleDColumn("entryX_mm");
  analysisManager->CreateNtupleDColumn("entryY_mm");
  analysisManager->CreateNtupleDColumn("entryZ_mm");
  analysisManager->CreateNtupleIColumn("nHitsPrimary");
  analysisManager->CreateNtupleIColumn("hitDetector");
  analysisManager->CreateNtupleSColumn("structureTag");
  analysisManager->CreateNtupleDColumn("foldState");
  analysisManager->CreateNtupleDColumn("localIncidenceDeg");
  analysisManager->CreateNtupleDColumn("vertexX_mm");
  analysisManager->CreateNtupleDColumn("vertexY_mm");
  analysisManager->CreateNtupleDColumn("vertexZ_mm");
  analysisManager->CreateNtupleDColumn("scatterAngleDeg");
  analysisManager->CreateNtupleDColumn("kaptonPathLength_mm");
  analysisManager->CreateNtupleDColumn("totalPathLength_mm");
  analysisManager->FinishNtuple();
}
void RunAction::DefineMessenger()
{
  fMessenger = std::make_unique<G4GenericMessenger>(
      this, "/origami/run/", "Run-level output controls");
  fMessenger->DeclareMethod("tag", &RunAction::SetRunTagCmd)
      .SetGuidance("MARK 3b: explicit identifier for THIS macro/sweep, "
                    "folded into the output ROOT filename (e.g. 'exp2', "
                    "'exp3', 'normal'). Set this once near the top of each "
                    "top-level macro, before any /run/beamOn, so different "
                    "macros landing on the same structureTag/foldState/"
                    "momentum/vertexSmear combination in separate process "
                    "launches still produce distinct filenames. Empty "
                    "string (default) omits the tag segment entirely.");
}
void RunAction::BeginOfRunAction(const G4Run* run)
{
  auto* analysisManager = G4AnalysisManager::Instance();
  std::ostringstream foldStr;
  foldStr << std::fixed << std::setprecision(2)
          << DetectorConstruction::GetFoldStateStatic();
  std::string foldTag = foldStr.str();
  for (char& c : foldTag) { if (c == '.') c = 'p'; }
  std::ostringstream momStr;
  momStr << std::fixed << std::setprecision(2)
         << PrimaryGeneratorAction::GetMomentumGeVStatic();
  std::string momTag = momStr.str();
  for (char& c : momTag) { if (c == '.') c = 'p'; }
  std::ostringstream vtxStr;
  vtxStr << std::fixed << std::setprecision(1)
         << PrimaryGeneratorAction::GetVertexSmearMmStatic();
  std::string vtxTag = vtxStr.str();
  for (char& c : vtxTag) { if (c == '.') c = 'p'; }
  std::string tagSegment;
  if (!RunAction::GetRunTagStatic().empty()) {
    tagSegment = "_tag" + RunAction::GetRunTagStatic();
  }
  std::ostringstream fname;
  fname << "origami_" << DetectorConstruction::GetStructureTag()
        << "_fold" << foldTag
        << "_p" << momTag << "GeV"
        << "_vtx" << vtxTag << "mm"
        << tagSegment
        << "_run" << run->GetRunID()
        << ".root";
  analysisManager->OpenFile(fname.str());
}
void RunAction::EndOfRunAction(const G4Run*)
{
  auto* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}
