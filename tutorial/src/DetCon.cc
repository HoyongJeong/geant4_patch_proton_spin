////////////////////////////////////////////////////////////////////////////////
//   DetCon.cc
//
//   Definitions of DetCon class's member functions. And it describes geometry
// of simulation.
//
//                        - Feb 5th 2025, Hoyong Jeong (hoyong5419@korea.ac.kr)
////////////////////////////////////////////////////////////////////////////////



#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

#include "DetCon.hh"



DetCon::DetCon()
{
	ConstructMaterials();
	DefineDimensions();
}



DetCon::~DetCon()
{
	DestructMaterials();
}



void DetCon::DefineDimensions()
{
	//------------------------------------------------
	// World dimensions: Laboratory size
	//------------------------------------------------
	labX    = 2000.0 * mm; // World x dimension
	labY    = 2000.0 * mm; // World y dimension
	labZ    = 2000.0 * mm; // World z dimension


	//------------------------------------------------
	// Target dimensions (spherical)
	//------------------------------------------------
	tarD    =    1.0 * mm; // Target diameter


	//------------------------------------------------
	// (cylindrical, optional)
	//------------------------------------------------
	tarL    =   60.0 * mm; // Target length


	//------------------------------------------------
	// Detector dimensions
	//------------------------------------------------
	detID   =    156 * mm; // Inner diameter
	detOD   =    700 * mm; // Outer diameter
	detT    =      3 * mm; // Thickness
}



G4VPhysicalVolume* DetCon::Construct()
{
	//------------------------------------------------
	// World
	//------------------------------------------------
	G4Box* labSolid = new G4Box("labSolid", labX / 2., labY / 2., labZ / 2.);
	G4LogicalVolume* labLV = new G4LogicalVolume(labSolid, Air, "labLV");
	labLV -> SetVisAttributes(new G4VisAttributes(0));
	G4VPhysicalVolume* labPV = new G4PVPlacement(0, G4ThreeVector(), "labPV", labLV, 0, false, 0);

	//------------------------------------------------
	// Spherical target
	//------------------------------------------------
	G4Sphere* target = new G4Sphere("target", 0., tarD / 2., 0., 360.*degree, 0., 180.*degree);
	G4LogicalVolume* tarLV = new G4LogicalVolume(target, PureCarbon, "tarLV");
	new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), "tarPV", tarLV, labPV, false, 0);


	//------------------------------------------------
	// Detector
	//------------------------------------------------
	G4Tubs* detSolid = new G4Tubs("detSolid", detID / 2., detOD / 2., detT / 2., 0.*degree, 360.*degree);
	G4LogicalVolume* detLV = new G4LogicalVolume(detSolid, Scint, "detLV");
	new G4PVPlacement(0, G4ThreeVector(0, 0, 900.*mm), "detPV", detLV, labPV, false, 0);


	return labPV;
}



void DetCon::ConstructMaterials()
{
	const G4double labTemp = 300.0 * kelvin;

	elN  = new G4Element("Nitrogen",  "N",  7, 14.00674*g/mole);
	elO  = new G4Element(  "Oxygen",  "O",  8, 15.9994 *g/mole);
	elAr = new G4Element(   "Argon", "Ar", 18, 39.948  *g/mole);
	elC  = new G4Element(  "Carbon",  "C",  6, 12.011  *g/mole);
	elH  = new G4Element("Hydrogen",  "H",  1,  1.00794*g/mole);

	Air = new G4Material("Air", 1.2929e-10*g/cm3, 3, kStateGas, labTemp);
	Air -> AddElement(elN, 75.47/99.95);
	Air -> AddElement(elO, 23.20/99.95);
	Air -> AddElement(elAr, 1.28/99.95);

//	Graphite = new G4Material("Graphite", 1.7*g/cm3, 3, kStateSolid, labTemp);
//	Graphite -> AddElement(elC, 99.0/100.0);
//	Graphite -> AddElement(elN,  0.7/100.0);
//	Graphite -> AddElement(elO,  0.3/100.0);

	PureCarbon = new G4Material("PureCarbon", 1.7*g/cm3, 1, kStateSolid, labTemp);
	PureCarbon -> AddElement(elC, 1);

	Scint = new G4Material("Scintillator", 1.032*g/cm3, 2, kStateSolid, labTemp);
	Scint -> AddElement(elC, 10);
	Scint -> AddElement(elH, 11);
}



void DetCon::DestructMaterials()
{
	delete Scint;
	delete PureCarbon;
//	delete Graphite;
	delete Air;

	delete elH;
	delete elC;
	delete elAr;
	delete elO;
	delete elN;
}
