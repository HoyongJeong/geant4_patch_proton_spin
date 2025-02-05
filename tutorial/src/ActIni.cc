////////////////////////////////////////////////////////////////////////////////
//   ActIni.cc
//
//   Definitions of ActIni class's member functions.
// All actions must be initialized here in order to use multi thread.
//
//                        - Feb 5th 2025, Hoyong Jeong (hoyong5419@korea.ac.kr)
////////////////////////////////////////////////////////////////////////////////



#include "ActIni.hh"
#include "PriGenAct.hh"
#include "RunAct.hh"
#include "EveAct.hh"
#include "SteAct.hh"



ActIni::ActIni() : G4VUserActionInitialization()
{
}



ActIni::~ActIni()
{
}



void ActIni::BuildForMaster() const
{
	SetUserAction(new RunAct);
}



void ActIni::Build() const
{
	SetUserAction(new PriGenAct);
	SetUserAction(new RunAct);
  
	EveAct* EA = new EveAct;
	SetUserAction(EA);
  
	SetUserAction(new SteAct(EA));
}  
