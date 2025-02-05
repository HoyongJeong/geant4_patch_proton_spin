////////////////////////////////////////////////////////////////////////////////
//   PriGenAct.hh
//
//   This file is a header for PriGenAct class. You can set primary beam
// options in this class.
//
//                        - Feb 5th 2025, Hoyong Jeong (hoyong5419@korea.ac.kr)
////////////////////////////////////////////////////////////////////////////////



#ifndef PRIGENACT_h
#define PRIGENACT_h 1



#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"



class G4ParticleGun;



class PriGenAct: public G4VUserPrimaryGeneratorAction
{
  public:
	PriGenAct();
	~PriGenAct();

	virtual void GeneratePrimaries(G4Event* anEvent);


  private:
	G4ParticleGun* PG;

	G4ThreeVector gunPos;
	G4ParticleDefinition* par;
	G4ThreeVector momDir;
	G4double mom;
	G4double kinEgy;
	G4ThreeVector pol;

	G4double displace;
	G4double rotate;
};



#endif
