////////////////////////////////////////////////////////////////////////////////
//   DetCon.hh
//
//   This file is a header for DetCon class. It's for construction of whole
// geometry of simulation, which includes detector geometry.
//
//                        - Feb 5th 2025, Hoyong Jeong (hoyong5419@korea.ac.kr)
////////////////////////////////////////////////////////////////////////////////



#ifndef DETCON_h
#define DETCON_h 1



#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Element.hh"
#include "G4Material.hh"



class G4VPhysicalVolume;



class DetCon: public G4VUserDetectorConstruction
{
  public:
	DetCon();
	virtual  ~DetCon();
	virtual G4VPhysicalVolume* Construct();


  private:
	void DefineDimensions();
	void ConstructMaterials();
	void DestructMaterials();

	G4Element* elN;
	G4Element* elO;
	G4Element* elAr;
	G4Element* elC;
	G4Element* elH;

	G4Material* Air;
//	G4Material* Graphite;
	G4Material* PureCarbon;
	G4Material* Scint;

	G4double labX, labY, labZ;
	G4double detID, detOD, detT;
	G4double tarD, tarL;
};



#endif
