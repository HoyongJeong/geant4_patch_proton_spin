////////////////////////////////////////////////////////////////////////////////
//   ActIni.hh
//
//   This file is a header for ActIni class. Every actions are initialized
// through this class.
//
//                        - Feb 5th 2025, Hoyong Jeong (hoyong5419@korea.ac.kr)
////////////////////////////////////////////////////////////////////////////////



#ifndef ACTINI_h
#define ACTINI_h 1



#include "G4VUserActionInitialization.hh"



class ActIni : public G4VUserActionInitialization
{
public:
	ActIni();
	virtual ~ActIni();

	virtual void BuildForMaster() const;
	virtual void Build() const;
};



#endif
