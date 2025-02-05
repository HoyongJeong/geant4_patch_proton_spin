////////////////////////////////////////////////////////////////////////////////
//   EveAct.cc
//
//   Definitions of EveAct class's member functions. Details of user actions
// are here.
//
//                        - Feb 5th 2025, Hoyong Jeong (hoyong5419@korea.ac.kr)
////////////////////////////////////////////////////////////////////////////////



#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

#include "EveAct.hh"
#include "Ana.hh"



EveAct::EveAct()
{
}



EveAct::~EveAct()
{
}



void EveAct::BeginOfEventAction(const G4Event* /*anEvent*/)
{
	trackID         . clear();
	parName         . clear();
	parID           . clear();
	parCharge       . clear();
	procHist        . clear();
	hadElastic      . clear();
	hIoni           . clear();
	CoulombScat     . clear();
	protonInelastic . clear();
	etc             . clear();
	exitTar         . clear();
	posTar          . clear();
	momTar          . clear();
	kinEgyTar       . clear();
	eDepTar         . clear();
	posDet          . clear();
}



void EveAct::EndOfEventAction(const G4Event* /*anEvent*/)
{
	//------------------------------------------------
	// Get analysis manager
	//------------------------------------------------
	G4AnalysisManager* AM = G4AnalysisManager::Instance();

	G4int iTrack = 1;
	while ( trackID[iTrack] != 0 )
	{
		G4int ID = trackID[iTrack];

		AM -> FillNtupleIColumn( 0, ID                  );
		AM -> FillNtupleSColumn( 1, parName[ID]         );
		AM -> FillNtupleIColumn( 2, parID[ID]           );
		AM -> FillNtupleIColumn( 3, parCharge[ID]       );
		AM -> FillNtupleSColumn( 4, procHist[ID]        );
		AM -> FillNtupleIColumn( 5, hadElastic[ID]      );
		AM -> FillNtupleIColumn( 6, hIoni[ID]           );
		AM -> FillNtupleIColumn( 7, CoulombScat[ID]     );
		AM -> FillNtupleIColumn( 8, protonInelastic[ID] );
		AM -> FillNtupleIColumn( 9, etc[ID]             );
		AM -> FillNtupleIColumn(10, exitTar[ID]         );
		AM -> FillNtupleDColumn(11, posTar[ID] . x()    );
		AM -> FillNtupleDColumn(12, posTar[ID] . y()    );
		AM -> FillNtupleDColumn(13, posTar[ID] . z()    );
		AM -> FillNtupleDColumn(14, posTar[ID] . theta());
		AM -> FillNtupleDColumn(15, posTar[ID] . phi()  );
		AM -> FillNtupleDColumn(16, momTar[ID] . x()    );
		AM -> FillNtupleDColumn(17, momTar[ID] . y()    );
		AM -> FillNtupleDColumn(18, momTar[ID] . z()    );
		AM -> FillNtupleDColumn(19, momTar[ID] . theta());
		AM -> FillNtupleDColumn(20, momTar[ID] . phi()  );
		AM -> FillNtupleDColumn(21, kinEgyTar[ID]       );
		AM -> FillNtupleDColumn(22, eDepTar[ID]         );
		AM -> FillNtupleDColumn(23, posDet[ID] . x()    );
		AM -> FillNtupleDColumn(24, posDet[ID] . y()    );

		AM -> AddNtupleRow();
		iTrack++;
	}
}



void EveAct::SetTrackID(G4int tID)
{ trackID[tID] = tID; }


void EveAct::SetParName(G4int tID, G4String pN)
{ parName[tID] = pN; }


void EveAct::SetParID(G4int tID, G4int pID)
{ parID[tID] = pID; }


void EveAct::SetParCharge(G4int tID, G4int pC)
{ parCharge[tID] = pC; }


void EveAct::AddProcHist(G4int tID, G4String pN)
{ procHist[tID] += pN + " "; }


void EveAct::AddhadElastic(G4int tID)
{ hadElastic[tID]++; }


void EveAct::AddhIoni(G4int tID)
{ hIoni[tID]++; }


void EveAct::AddCoulombScat(G4int tID)
{ CoulombScat[tID]++; }


void EveAct::AddprotonInelastic(G4int tID)
{ protonInelastic[tID]++; }


void EveAct::AddEtc(G4int tID)
{ etc[tID]++; }



//////////////////////////////////////////////////
// For target
//////////////////////////////////////////////////
void EveAct::SetExitTar(G4int tID)
{ exitTar[tID] = 1; }


void EveAct::SetPosTar(G4int tID, G4ThreeVector pP)
{ posTar[tID] = pP; }


void EveAct::SetMomTar(G4int tID, G4ThreeVector pM)
{ momTar[tID] = pM; }


void EveAct::SetKinEgyTar(G4int tID, G4double pKE)
{ kinEgyTar[tID] = pKE; }


void EveAct::AddEDepTar(G4int tID, G4double eD)
{ eDepTar[tID] += eD; }



//////////////////////////////////////////////////
// For detector
//////////////////////////////////////////////////
void EveAct::SetPosDet(G4int tID, G4ThreeVector pP)
{ posDet[tID] = pP; }
