#ifndef RUN_ACTION_HH
#define RUN_ACTION_HH
#include "G4UserRunAction.hh"
#include "G4GenericMessenger.hh"
#include "globals.hh"
#include <memory>
class G4Run;
class RunAction : public G4UserRunAction
{
public:
  RunAction();
  ~RunAction() override = default;
  void BeginOfRunAction(const G4Run* run) override;
  void EndOfRunAction(const G4Run* run) override;
  static void SetRunTag(G4String tag) { fsRunTag = tag; }
  static G4String GetRunTagStatic() { return fsRunTag; }
private:
  std::unique_ptr<G4GenericMessenger> fMessenger;
  void DefineMessenger();
  void SetRunTagCmd(G4String tag) { SetRunTag(tag); }
  static G4String fsRunTag;
};
#endif
