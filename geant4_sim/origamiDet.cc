#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#ifdef _WIN32
    #include <io.h>
    #include <process.h>
    #define F_OK 0
    #define access _access
    #define sleep(x) Sleep((x) * 1000)
#else
    #include <unistd.h>
#endif
int main(int argc, char** argv)
{
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4long seed = static_cast<G4long>(time(0)) * 1000 + (getpid() % 1000);
  G4Random::setTheSeed(seed);
  G4UIExecutive* ui = nullptr;
  if (argc == 1 || (argc > 1 && G4String(argv[1]).find("vis") != std::string::npos)) {
    ui = new G4UIExecutive(argc, argv);
  }
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetUserInitialization(new DetectorConstruction());
  auto* physicsList = new FTFP_BERT();
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new ActionInitialization());
  G4VisManager* visManager = new G4VisExecutive();
  visManager->Initialize();
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if (ui) {
    if (argc > 1) {
      UImanager->ApplyCommand(G4String("/control/execute ") + argv[1]);
    } else {
      UImanager->ApplyCommand("/control/execute macros/init_vis.mac");
    }
    ui->SessionStart();
    delete ui;
  } else {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }
  delete visManager;
  delete runManager;
  return 0;
}
