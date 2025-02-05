////////////////////////////////////////////////////////////////////////////////
//   SteAct.cc
//
//   Definitions of SteAct class's member functions. Details of user
// actions are here.
//
//                        - Feb 5th 2025, Hoyong Jeong (hoyong5419@korea.ac.kr)
////////////////////////////////////////////////////////////////////////////////



#include "G4VProcess.hh"

#include "SteAct.hh"
#include "EveAct.hh"



SteAct::SteAct(EveAct* EA) : G4UserSteppingAction(), m_EA(EA)
{
}



SteAct::~SteAct()
{
}



void SteAct::UserSteppingAction(const G4Step* step)
{
	//------------------------------------------------
	// Who am I? Where am I? What am I undergoing?
	//------------------------------------------------
	// Track ID
	G4int trackID = step -> GetTrack() -> GetTrackID();

	// Particle name
	G4String parName = step -> GetTrack() -> GetDefinition() -> GetParticleName();

	// Particle ID
	G4int parID = step -> GetTrack() -> GetDefinition() -> GetPDGEncoding();

	// Particle charge
	G4int parCharge = step -> GetTrack() -> GetDefinition() -> GetPDGCharge();

	// Process name
	const G4VProcess* creProc = step -> GetTrack() -> GetCreatorProcess();
	G4String procName = step -> GetPostStepPoint() -> GetProcessDefinedStep() -> GetProcessName();

	// Physical volume
	G4String namePrePV = step -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName();
	G4String namePostPV;
	G4VPhysicalVolume* postPV = step -> GetPostStepPoint() -> GetPhysicalVolume();
	if ( postPV != 0 ) namePostPV = postPV -> GetName();
	else namePostPV = "outside";

	// Status
	G4StepStatus preStat = step -> GetPreStepPoint() -> GetStepStatus();
	G4StepStatus postStat = step -> GetPostStepPoint() -> GetStepStatus();

	// Position
	G4ThreeVector prePos = step -> GetPreStepPoint() -> GetPosition();
	G4ThreeVector postPos = step -> GetPostStepPoint() -> GetPosition();

	// Momentum
	G4ThreeVector preMom = step -> GetPreStepPoint() -> GetMomentum();
	G4ThreeVector postMom = step -> GetPostStepPoint() -> GetMomentum();

	// Kinetic energy
	G4double preKinEgy = step -> GetPreStepPoint() -> GetKineticEnergy();
	G4double postKinEgy = step -> GetPostStepPoint() -> GetKineticEnergy();

	// Energy deposition
	G4double eDep = step -> GetTotalEnergyDeposit();


	//------------------------------------------------
	// Set track ID
	//------------------------------------------------
	m_EA -> SetTrackID(trackID);


	//------------------------------------------------
	// Set particle name
	//------------------------------------------------
	m_EA -> SetParName(trackID, parName);


	//------------------------------------------------
	// Set particle ID
	//------------------------------------------------
	m_EA -> SetParID(trackID, parID);


	//------------------------------------------------
	// Set particle charge
	//------------------------------------------------
	m_EA -> SetParCharge(trackID, parCharge);


	//------------------------------------------------
	// Add process history inside target
	//------------------------------------------------
	if ( (creProc) && (namePrePV == "tarPV") )
	{
		G4String creProcName = creProc -> GetProcessName();
		m_EA -> AddProcHist(trackID, creProcName);
	}

	if ( (procName != "Transportation") && (namePrePV == "tarPV") )
	{
		m_EA -> AddProcHist(trackID, procName);
		// G4cout << procName << G4endl;
		if      ( procName == "hadElastic"      ) m_EA -> AddhadElastic(trackID);
		else if ( procName == "hIoni"           ) m_EA -> AddhIoni(trackID);
		else if ( procName == "CoulombScat"     ) m_EA -> AddCoulombScat(trackID);
		else if ( procName == "protonInelastic" ) m_EA -> AddprotonInelastic(trackID);
		else                                      m_EA -> AddEtc(trackID);
	}


	//------------------------------------------------
	// Set target exit position
	//------------------------------------------------
	if ( (namePrePV == "tarPV") && (postStat == fGeomBoundary) )
	{
		m_EA -> SetExitTar(trackID);
		m_EA -> SetPosTar(trackID, postPos);
		m_EA -> SetMomTar(trackID, postMom);
		m_EA -> SetKinEgyTar(trackID, postKinEgy);
	}


	//------------------------------------------------
	// Add energy deposition in target
	//------------------------------------------------
	if ( namePrePV == "tarPV" )
	{
		m_EA -> AddEDepTar(trackID, eDep);
	}


	//------------------------------------------------
	// Set detector position
	//------------------------------------------------
	if ( (namePrePV == "labPV") && (namePostPV.contains("det")) )
	{
		m_EA -> SetPosDet(trackID, postPos);
	}
}
