////////////////////////////////////////////////////////////////////////////////
//   RunAct.cc
//
//   Definitions of RunAct class's member functions. Details of the user
// actions are here.
//
//                        - Feb 5th 2025, Hoyong Jeong (hoyong5419@korea.ac.kr)
////////////////////////////////////////////////////////////////////////////////



#include "G4Run.hh"
#include "G4SystemOfUnits.hh"

#include "RunAct.hh"
#include "Ana.hh"



RunAct::RunAct() : G4UserRunAction()
{
	//------------------------------------------------
	// Create analysis manager
	// The choice of analysis technology is done via
	// selectin of a namespace in Ana.hh
	//------------------------------------------------
	G4AnalysisManager* AM = G4AnalysisManager::Instance();
	G4cout << "Using " << AM -> GetType() << G4endl;


	//------------------------------------------------
	// Default settings
	//------------------------------------------------
	AM -> SetVerboseLevel(0);


	//------------------------------------------------
	// Creating ntuple
	//------------------------------------------------
	AM -> CreateNtuple("hits", "data");
	AM -> CreateNtupleIColumn("trackID");         // Column ID =  0
	AM -> CreateNtupleSColumn("parName");         // Column ID =  1
	AM -> CreateNtupleIColumn("parID");           // Column ID =  2
	AM -> CreateNtupleIColumn("parCharge");       // Column ID =  3
	AM -> CreateNtupleSColumn("procHist");        // Column ID =  4
	AM -> CreateNtupleIColumn("hadElastic");      // Column ID =  5
	AM -> CreateNtupleIColumn("hIoni");           // Column ID =  6
	AM -> CreateNtupleIColumn("CoulombScat");     // Column ID =  7
	AM -> CreateNtupleIColumn("protonInelastic"); // Column ID =  8
	AM -> CreateNtupleIColumn("etc");             // Column ID =  9
	AM -> CreateNtupleIColumn("exitTar");         // Column ID = 10
	AM -> CreateNtupleDColumn("posXTar");         // Column ID = 11
	AM -> CreateNtupleDColumn("posYTar");         // Column ID = 12
	AM -> CreateNtupleDColumn("posZTar");         // Column ID = 13
	AM -> CreateNtupleDColumn("posThetaTar");     // Column ID = 14
	AM -> CreateNtupleDColumn("posPhiTar");       // Column ID = 15
	AM -> CreateNtupleDColumn("momXTar");         // Column ID = 16
	AM -> CreateNtupleDColumn("momYTar");         // Column ID = 17
	AM -> CreateNtupleDColumn("momZTar");         // Column ID = 18
	AM -> CreateNtupleDColumn("momThetaTar");     // Column ID = 19
	AM -> CreateNtupleDColumn("momPhiTar");       // Column ID = 20
	AM -> CreateNtupleDColumn("kinEgyTar");       // Column ID = 21
	AM -> CreateNtupleDColumn("eDepTar");         // Column ID = 22
	AM -> CreateNtupleDColumn("posXDet");         // Column ID = 23
	AM -> CreateNtupleDColumn("posYDet");         // Column ID = 24
	AM -> FinishNtuple();
}



RunAct::~RunAct()
{
	delete G4AnalysisManager::Instance();
}



void RunAct::BeginOfRunAction(const G4Run* /*run*/)
{
	//------------------------------------------------
	// Inform the runManager to save random number seed
	//------------------------------------------------
//	G4RunManager::GetRunManager()->SetRandomNumberStore(true);


	//------------------------------------------------
	// Get analysis manager
	//------------------------------------------------
	G4AnalysisManager* AM = G4AnalysisManager::Instance();


	//------------------------------------------------
	// Open an output file
	// The default file name is set in RunAct::RunAct(),
	// it can be overwritten in a macro
	//------------------------------------------------
	G4String fileName = "tutorial";
	AM -> OpenFile(fileName);
}



void RunAct::EndOfRunAction(const G4Run* /*run*/)
{
	//------------------------------------------------
	// Save histograms & ntuple
	//------------------------------------------------
	G4AnalysisManager* AM = G4AnalysisManager::Instance();
	AM -> Write();
	AM -> CloseFile();
}
