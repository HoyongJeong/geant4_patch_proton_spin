////////////////////////////////////////////////////////////////////////////////
//   EveAct.hh
//
//   This file is a header for EveAct class. User can add user-defined event
// action in this class. So this class works at every event.
//
//                        - Feb 5th 2025, Hoyong Jeong (hoyong5419@korea.ac.kr)
////////////////////////////////////////////////////////////////////////////////



#ifndef EVEACT_h
#define EVEACT_h 1



#include "G4UserEventAction.hh" 
#include <map>



class G4Event;



class EveAct: public G4UserEventAction
{
  public:
	EveAct();
	virtual ~EveAct();

	virtual void BeginOfEventAction(const G4Event*);
	virtual void EndOfEventAction(const G4Event*);

	void SetTrackID(G4int);
	void SetParName(G4int, G4String);
	void SetParID(G4int, G4int);
	void SetParCharge(G4int, G4int);
	void AddProcHist(G4int, G4String);
	void AddhadElastic(G4int);
	void AddhIoni(G4int);
	void AddCoulombScat(G4int);
	void AddprotonInelastic(G4int);
	void AddEtc(G4int);
	void SetExitTar(G4int);
	void SetPosTar(G4int, G4ThreeVector);
	void SetMomTar(G4int, G4ThreeVector);
	void SetKinEgyTar(G4int, G4double);
	void AddEDepTar(G4int, G4double);
	void SetPosDet(G4int, G4ThreeVector);


  private:
	std::map<G4int, G4int> trackID;
	std::map<G4int, G4String> parName;
	std::map<G4int, G4int> parID;
	std::map<G4int, G4int> parCharge;
	std::map<G4int, G4String> procHist;
	std::map<G4int, G4int> hadElastic;
	std::map<G4int, G4int> hIoni;
	std::map<G4int, G4int> CoulombScat;
	std::map<G4int, G4int> protonInelastic;
	std::map<G4int, G4int> etc;
	std::map<G4int, G4int> exitTar;
	std::map<G4int, G4ThreeVector> posTar;
	std::map<G4int, G4ThreeVector> momTar;
	std::map<G4int, G4double> kinEgyTar;
	std::map<G4int, G4double> eDepTar;
	std::map<G4int, G4ThreeVector> posDet;
};



#endif
