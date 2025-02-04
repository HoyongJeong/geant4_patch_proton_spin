//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// Geant4 Header : G4HadronElastic
//
// Author : V.Ivanchenko 29 June 2009 (redesign old elastic model)
//  
// Modified:
//
// Class Description
// Default model for elastic scattering; GHEISHA algorithm is used 
// Class Description - End

#ifndef G4HadronElastic_h
#define G4HadronElastic_h 1
 
#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"

class G4ParticleDefinition;

class G4HadronElastic : public G4HadronicInteraction
{
public:

  explicit G4HadronElastic(const G4String& name = "hElasticLHEP");

  ~G4HadronElastic() override;
 
  // implementation of the G4HadronicInteraction interface
  G4HadFinalState* ApplyYourself(const G4HadProjectile & aTrack, 
				 G4Nucleus & targetNucleus) override;

  // sample momentum transfer using Lab. momentum
  G4double SampleInvariantT(const G4ParticleDefinition* p, G4double plab,
			    G4int Z, G4int A) override;
  
  G4double GetSlopeCof( const G4int pdg );


  ////////////////////////////////////////////////////////////////////////////////
  // pEDM v1.0
  //
  // Modified: 15-Jan-16 by H. Jeong
  ////////////////////////////////////////////////////////////////////////////////
  // sample phi with asymmetry
  virtual G4double transformFunction(G4double x, G4double A, G4double r, G4double polPer);
  virtual G4double getAP(G4double costheta);
  virtual G4double getAPFromT(G4double t);
  virtual G4double getBetaFromT(G4double t, G4double eKinetic);
  virtual G4double APSpline(G4double x);
  virtual G4double betaSpline150(G4double x);
  virtual G4double betaSpline160(G4double x);
  virtual G4double betaSpline170(G4double x);
  virtual G4double betaSpline180(G4double x);
  virtual G4double betaSpline190(G4double x);
  virtual G4double betaSpline200(G4double x);
  virtual G4double betaSpline210(G4double x);
  virtual G4double betaSpline220(G4double x);
  virtual G4double betaSpline230(G4double x);
  virtual G4double betaSpline240(G4double x);
  virtual G4double CSSpline(G4double x);
  virtual G4double CSSpline2(G4double x);
  virtual G4double SampleInvariantT2(G4double tm);
  virtual G4double BisectionPhiWithT(G4double t, G4double polPerpen);
  ////////////////////////////////////////////////////////////////////////////////

  inline void SetLowestEnergyLimit(G4double value);

  inline G4double LowestEnergyLimit() const;

  inline G4double ComputeMomentumCMS(const G4ParticleDefinition* p, 
				     G4double plab, G4int Z, G4int A);
  
  void ModelDescription(std::ostream&) const override;

protected:

  G4double pLocalTmax;

private:

  G4ParticleDefinition* theProton;
  G4ParticleDefinition* theNeutron;
  G4ParticleDefinition* theDeuteron;
  G4ParticleDefinition* theAlpha;

  G4double lowestEnergyLimit;
  G4int nwarn;
};

inline void G4HadronElastic::SetLowestEnergyLimit(G4double value)
{
  lowestEnergyLimit = value;
}

inline G4double G4HadronElastic::LowestEnergyLimit() const
{
  return lowestEnergyLimit;
}

inline G4double
G4HadronElastic::ComputeMomentumCMS(const G4ParticleDefinition* p, 
				    G4double plab, G4int Z, G4int A)
{
  G4double m1 = p->GetPDGMass();
  G4double m12= m1*m1;
  G4double mass2 = G4NucleiProperties::GetNuclearMass(A, Z);
  return plab*mass2/std::sqrt(m12 + mass2*mass2 + 2.*mass2*std::sqrt(m12 + plab*plab));
}

#endif
