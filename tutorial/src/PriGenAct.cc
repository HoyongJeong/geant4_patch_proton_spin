////////////////////////////////////////////////////////////////////////////////
//   PriGenAct.cc
//
//   Definitions of PriGenAct class's member functions.
//
//                        - Feb 5th 2025, Hoyong Jeong (hoyong5419@korea.ac.kr)
////////////////////////////////////////////////////////////////////////////////



#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include "PriGenAct.hh"



PriGenAct::PriGenAct()
{
	PG = new G4ParticleGun();

	//------------------------------------------------
	// Set misalign
	//------------------------------------------------
	displace = 0.0 * mm;
	rotate   = 0; // in radian


	//------------------------------------------------
	// Gun position
	//------------------------------------------------
	gunPos = G4ThreeVector(displace - 2000*std::tan(rotate)*mm, 0, -2000*mm);
	PG -> SetParticlePosition(gunPos);


	//------------------------------------------------
	// Set particle definition
	//------------------------------------------------
	G4ParticleTable* PT = G4ParticleTable::GetParticleTable();
	par = PT -> FindParticle("proton");
	PG -> SetParticleDefinition(par);


	//------------------------------------------------
	// Momentum
	//------------------------------------------------
	momDir = G4ThreeVector(2000*std::tan(rotate), 0, 2000);
	PG -> SetParticleMomentumDirection(momDir);
	PG -> SetParticleEnergy(200.0*MeV);


	//------------------------------------------------
	// Polarization: up!
	//------------------------------------------------
	pol = G4ThreeVector(0.0, 1.0, 0.0);
	PG -> SetParticlePolarization(pol);
}



PriGenAct::~PriGenAct()
{
	delete PG;
}



void PriGenAct::GeneratePrimaries(G4Event* anEvent)
{
	PG -> GeneratePrimaryVertex(anEvent);
}
