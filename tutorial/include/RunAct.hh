////////////////////////////////////////////////////////////////////////////////
//   RunAct.hh
//
//   This file is a header for RunAct class. In addition to the basic actions,
// user can add additional user-defined actions.
//
//                        - Feb 5th 2025, Hoyong Jeong (hoyong5419@korea.ac.kr)
////////////////////////////////////////////////////////////////////////////////



#ifndef RUNACT_h
#define RUNACT_h 1



#include "G4UserRunAction.hh"



class G4Run;



class RunAct : public G4UserRunAction
{
  public:
	RunAct();
	virtual ~RunAct();

	virtual void BeginOfRunAction(const G4Run*);
	virtual void   EndOfRunAction(const G4Run*);
};



#endif
