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

#include "G4HadronElastic.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Alpha.hh"
#include "G4Pow.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4HadronicParameters.hh"


G4HadronElastic::G4HadronElastic(const G4String& name) 
  : G4HadronicInteraction(name)
{
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  lowestEnergyLimit= 1.e-6*eV;
  pLocalTmax  = 0.0;
  nwarn = 0;

  theProton   = G4Proton::Proton();
  theNeutron  = G4Neutron::Neutron();
  theDeuteron = G4Deuteron::Deuteron();
  theAlpha    = G4Alpha::Alpha();
}

G4HadronElastic::~G4HadronElastic()
{}


void G4HadronElastic::ModelDescription(std::ostream& outFile) const
{
  outFile << "G4HadronElastic is the base class for all hadron-nucleus\n" 
          << "elastic scattering models except HP.\n" 
          << "By default it uses the Gheisha two-exponential momentum\n"
	  << "transfer parameterization.  The model is fully relativistic\n"
	  << "as opposed to the original Gheisha model which was not.\n"
	  << "This model may be used for all long-lived hadrons at all\n"
	  << "incident energies but fit the data only for relativistic scattering.\n";
}

G4HadFinalState* G4HadronElastic::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();

  const G4HadProjectile* aParticle = &aTrack;
  G4double ekin = aParticle->GetKineticEnergy();

  // no scattering below the limit
  if(ekin <= lowestEnergyLimit) {
    theParticleChange.SetEnergyChange(ekin);
    theParticleChange.SetMomentumChange(0.,0.,1.);
    return &theParticleChange;
  }

  G4int A = targetNucleus.GetA_asInt();
  G4int Z = targetNucleus.GetZ_asInt();

  // Scattered particle referred to axis of incident particle
  const G4ParticleDefinition* theParticle = aParticle->GetDefinition();
  G4double m1 = theParticle->GetPDGMass();
  G4double plab = std::sqrt(ekin*(ekin + 2.0*m1));

  if (verboseLevel>1) {
    G4cout << "G4HadronElastic: " 
	   << aParticle->GetDefinition()->GetParticleName() 
	   << " Plab(GeV/c)= " << plab/GeV  
	   << " Ekin(MeV) = " << ekin/MeV 
	   << " scattered off Z= " << Z 
	   << " A= " << A 
	   << G4endl;
  }

  G4double mass2 = G4NucleiProperties::GetNuclearMass(A, Z);
  G4double e1 = m1 + ekin;
  G4LorentzVector lv(0.0,0.0,plab,e1+mass2);
  G4ThreeVector bst = lv.boostVector();
  G4double momentumCMS = plab*mass2/std::sqrt(m1*m1 + mass2*mass2 + 2.*mass2*e1);

  pLocalTmax = 4.0*momentumCMS*momentumCMS;

  // Sampling in CM system
  G4double t = SampleInvariantT(theParticle, plab, Z, A);

  if(t < 0.0 || t > pLocalTmax) {
    // For the very rare cases where cos(theta) is greater than 1 or smaller than -1,
    // print some debugging information via a "JustWarning" exception, and resample
    // using the default algorithm
#ifdef G4VERBOSE
    if(nwarn < 2) {
      G4ExceptionDescription ed;
      ed << GetModelName() << " wrong sampling t= " << t << " tmax= " << pLocalTmax
	 << " for " << aParticle->GetDefinition()->GetParticleName() 
	 << " ekin=" << ekin << " MeV" 
	 << " off (Z,A)=(" << Z << "," << A << ") - will be resampled" << G4endl;
      G4Exception( "G4HadronElastic::ApplyYourself", "hadEla001", JustWarning, ed);
      ++nwarn;
    }
#endif
    t = G4HadronElastic::SampleInvariantT(theParticle, plab, Z, A);
  }

//  G4double phi  = G4UniformRand()*CLHEP::twopi;
  G4double cost = 1. - 2.0*t/pLocalTmax;

  if (cost > 1.0) { cost = 1.0; }
  else if(cost < -1.0) { cost = -1.0; } 

  G4double sint = std::sqrt((1.0-cost)*(1.0+cost));

  ////////////////////////////////////////////////////////////////////////////////
  // < pEDM v1.0 >
  // Modified for pEDM from here
  // 17-May-16 by H. Jeong
  ////////////////////////////////////////////////////////////////////////////////
  // This part describes phi sampling in given t.
  
  // Choice 1. Just random phi (spin-independent scattering)
  // G4double phi = G4UniformRadn()*CLHEP::twopi;

  // Choice 2. Spin-dependent scattering
  G4double polPhi = 0.0, polBeta = 0.0;
  if ( (aParticle->GetPolarization().x() != 0.0)
    || (aParticle->GetPolarization().y() != 0.0)
    || (aParticle->GetPolarization().z() != 0.0) )
  {
    polPhi  = aParticle->GetPolarization().phi() - 0.25*CLHEP::twopi;
    polBeta = aParticle->GetPolarization().theta();
  }
  G4double sinPolBeta = std::sin(polBeta);
  G4double phi = BisectionPhiWithT(t, sinPolBeta);
  ////////////////////////////////////////////////////////////////////////////////

  if (verboseLevel>1) {
    G4cout << " t= " << t << " tmax(GeV^2)= " << pLocalTmax/(GeV*GeV) 
	   << " Pcms(GeV)= " << momentumCMS/GeV << " cos(t)=" << cost 
	   << " sin(t)=" << sint << G4endl;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // < pEDM v1.0 >
  // Modified for pEDM from here
  // 17-May-16 by H. Jeong
  ////////////////////////////////////////////////////////////////////////////////
  // Polarization dependence
  G4ThreeVector v1(sint*std::cos(phi), sint*std::sin(phi), cost);
  G4ThreeVector it(0., 0., 1.);
  v1 . rotate(polPhi, it);
  v1 *= momentumCMS;
  G4LorentzVector nlv1(v1 . x(), v1 . y(), v1 . z(), std::sqrt(momentumCMS*momentumCMS + m1*m1));
  nlv1.boost(bst); 

  // It was the original part
//  G4LorentzVector nlv1(momentumCMS*sint*std::cos(phi),
//		       momentumCMS*sint*std::sin(phi),
//                       momentumCMS*cost,
//		       std::sqrt(momentumCMS*momentumCMS + m1*m1));

//  nlv1.boost(bst); 
  ////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  // < pEDM v1.0 >
  // Modified for pEDM from here
  // 17-May-16 by H. Jeong
  ////////////////////////////////////////////////////////////////////////////////
  // Describe polarization change here.
  G4ThreeVector polFinal;
  G4ThreeVector polInit = aParticle -> GetPolarization();
  if ( (polInit.x() != 0.0)
    || (polInit.y() != 0.0)
    || (polInit.z() != 0.0) )
  {
    polInit.rotate(-phi, it);
    G4double polS = polInit.x();
    G4double polN = polInit.y();
    G4double polL = polInit.z();
    if ( polN == 1.0 )
    {
      polInit.rotate(phi, it);
      polFinal = polInit;
    }
    else
    {
      G4double newPolN = (getAPFromT(t)+polN) / (1.0+getAPFromT(t)*polN);
      G4double newPolS = std::sqrt(1.0 - newPolN*newPolN);
      polFinal = polInit - G4ThreeVector(0.0, polN, 0.0);
      polFinal = newPolS * polFinal / std::sqrt(polS*polS+polL*polL);
      G4double beta = CLHEP::twopi * getBetaFromT(t, ekin) / 360.0;
      polFinal.rotate(-beta, G4ThreeVector(0.0, 1.0, 0.0));
      polFinal = polFinal + G4ThreeVector(0.0, newPolN, 0.0);
      polFinal.rotate(phi, it);
    }
  }
  else
    polFinal = polInit;
  ////////////////////////////////////////////////////////////////////////////////

  G4double eFinal = nlv1.e() - m1;
  if (verboseLevel > 1) {
    G4cout <<"G4HadronElastic: m= " << m1 << " Efin(MeV)= " << eFinal 
	   << " 4-M Final: " << nlv1 
	   << G4endl;
  }

  if(eFinal <= 0.0) { 
    theParticleChange.SetMomentumChange(0.0,0.0,1.0);
    theParticleChange.SetEnergyChange(0.0);
  } else {
    theParticleChange.SetMomentumChange(nlv1.vect().unit());
    theParticleChange.SetEnergyChange(eFinal);
  }
  lv -= nlv1;
  G4double erec =  std::max(lv.e() - mass2, 0.0);
  if (verboseLevel > 1) {
    G4cout << "Recoil: " <<" m= " << mass2 << " Erec(MeV)= " << erec
	   << " 4-mom: " << lv 
	   << G4endl;
  }
 
  // the recoil is created if kinetic energy above the threshold
  if(erec > GetRecoilEnergyThreshold()) {
    G4ParticleDefinition * theDef = nullptr;
    if(Z == 1 && A == 1)       { theDef = theProton; }
    else if (Z == 1 && A == 2) { theDef = theDeuteron; }
    else if (Z == 1 && A == 3) { theDef = G4Triton::Triton(); }
    else if (Z == 2 && A == 3) { theDef = G4He3::He3(); }
    else if (Z == 2 && A == 4) { theDef = theAlpha; }
    else {
      theDef = 
	G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Z,A,0.0);
    }
    G4DynamicParticle * aSec = new G4DynamicParticle(theDef, lv.vect().unit(), erec);
    theParticleChange.AddSecondary(aSec);
  } else {
    theParticleChange.SetLocalEnergyDeposit(erec);
  }

  return &theParticleChange;
}

// sample momentum transfer in the CMS system 
G4double 
G4HadronElastic::SampleInvariantT(const G4ParticleDefinition* part,
				  G4double mom, G4int, G4int A)
{
  const G4double plabLowLimit = 400.0*CLHEP::MeV;
  const G4double GeV2 = GeV*GeV;
  const G4double z07in13 = std::pow(0.7, 0.3333333333);

  G4int pdg = std::abs(part->GetPDGEncoding());
  G4double tmax = pLocalTmax/GeV2;

  G4double aa, bb, cc, dd;
  G4Pow* g4pow = G4Pow::GetInstance();
  if (A <= 62) {
    if (pdg == 211){ //Pions
      if(mom >= plabLowLimit){     //High energy
	bb = 14.5*g4pow->Z23(A);/*14.5*/
	dd = 10.;
	cc = 0.075*g4pow->Z13(A)/dd;//1.4
	//aa = g4pow->powZ(A, 1.93)/bb;//1.63
	aa = (A*A)/bb;//1.63
      } else {                       //Low energy
	bb = 29.*z07in13*z07in13*g4pow->Z23(A);
	dd = 15.;
	cc = 0.04*g4pow->Z13(A)/dd;//1.4
	aa = g4pow->powZ(A, 1.63)/bb;//1.63
      }
    } else { //Other particles
      bb = 14.5*g4pow->Z23(A);
      dd = 20.;
      aa = (A*A)/bb;//1.63
      cc = 1.4*g4pow->Z13(A)/dd;          
    }
      //===========================
  } else { //(A>62)
    if (pdg == 211) {
      if(mom >= plabLowLimit){ //high
	bb = 60.*z07in13*g4pow->Z13(A);//60
	dd = 30.;
	aa = 0.5*(A*A)/bb;//1.33
	cc = 4.*g4pow->powZ(A,0.4)/dd;//1:0.4     ---    2: 0.4
      } else { //low
	bb = 120.*z07in13*g4pow->Z13(A);//60
	dd = 30.;
	aa = 2.*g4pow->powZ(A,1.33)/bb;
	cc = 4.*g4pow->powZ(A,0.4)/dd;//1:0.4     ---    2: 0.4
      }
    } else {
      bb = 60.*g4pow->Z13(A);
      dd = 25.;
      aa = g4pow->powZ(A,1.33)/bb;//1.33
      cc = 0.2*g4pow->powZ(A,0.4)/dd;//1:0.4     ---    2: 0.4
    }
  }
  G4double q1 = 1.0 - G4Exp(-bb*tmax);
  G4double q2 = 1.0 - G4Exp(-dd*tmax);
  G4double s1 = q1*aa;
  G4double s2 = q2*cc;
  if((s1 + s2)*G4UniformRand() < s2) {
    q1 = q2;
    bb = dd;
  }
  return -GeV2*G4Log(1.0 - G4UniformRand()*q1)/bb;
}

//////////////////////////////////////////////
//
// Cofs for s-,c-,b-particles ds/dt slopes

G4double G4HadronElastic::GetSlopeCof(const G4int pdg )
{
  // The input parameter "pdg" should be the absolute value of the PDG code
  // (i.e. the same value for a particle and its antiparticle).

  G4double coeff = 1.0;

  // heavy barions

  static const G4double  lBarCof1S  = 0.88;
  static const G4double  lBarCof2S  = 0.76;
  static const G4double  lBarCof3S  = 0.64;
  static const G4double  lBarCof1C  = 0.784378;
  static const G4double  lBarCofSC  = 0.664378;
  static const G4double  lBarCof2SC = 0.544378;
  static const G4double  lBarCof1B  = 0.740659;
  static const G4double  lBarCofSB  = 0.620659;
  static const G4double  lBarCof2SB = 0.500659;
  
  if( pdg == 3122 || pdg == 3222 ||  pdg == 3112 || pdg == 3212  )
  {
    coeff = lBarCof1S; // Lambda, Sigma+, Sigma-, Sigma0

  } else if( pdg == 3322 || pdg == 3312   )
  {
    coeff = lBarCof2S; // Xi-, Xi0
  }
  else if( pdg == 3324)
  {
    coeff = lBarCof3S; // Omega
  }
  else if( pdg == 4122 ||  pdg == 4212 ||   pdg == 4222 ||   pdg == 4112   )
  {
    coeff = lBarCof1C; // LambdaC+, SigmaC+, SigmaC++, SigmaC0
  }
  else if( pdg == 4332 )
  {
    coeff = lBarCof2SC; // OmegaC
  }
  else if( pdg == 4232 || pdg == 4132 )
  {
    coeff = lBarCofSC; // XiC+, XiC0
  }
  else if( pdg == 5122 || pdg == 5222 || pdg == 5112 || pdg == 5212    )
  {
    coeff = lBarCof1B; // LambdaB, SigmaB+, SigmaB-, SigmaB0
  }
  else if( pdg == 5332 )
  {
    coeff = lBarCof2SB; // OmegaB-
  }
  else if( pdg == 5132 || pdg == 5232 ) // XiB-, XiB0
  {
    coeff = lBarCofSB;
  }
  // heavy mesons Kaons?
  static const G4double lMesCof1S = 0.82; // Kp/piP kaons?
  static const G4double llMesCof1C = 0.676568;
  static const G4double llMesCof1B = 0.610989;
  static const G4double llMesCof2C = 0.353135;
  static const G4double llMesCof2B = 0.221978;
  static const G4double llMesCofSC = 0.496568;
  static const G4double llMesCofSB = 0.430989;
  static const G4double llMesCofCB = 0.287557;
  static const G4double llMesCofEtaP = 0.88;
  static const G4double llMesCofEta = 0.76;

  if( pdg == 321 || pdg == 311 || pdg == 310 )
  {
    coeff = lMesCof1S; //K+-0
  }
  else if( pdg == 511 ||  pdg == 521  )
  {
    coeff = llMesCof1B; // BMeson0, BMeson+
  }
  else if(pdg == 421 ||  pdg == 411 )
  {
    coeff = llMesCof1C; // DMeson+, DMeson0
  }
  else if( pdg == 531  )
  {
    coeff = llMesCofSB; // BSMeson0
  }
  else if( pdg == 541 )
  {
    coeff = llMesCofCB; // BCMeson+-
  }
  else if(pdg == 431 ) 
  {
    coeff = llMesCofSC; // DSMeson+-
  }
  else if(pdg == 441 || pdg == 443 )
  {
    coeff = llMesCof2C; // Etac, JPsi
  }
  else if(pdg == 553 )
  {
    coeff = llMesCof2B; // Upsilon
  }
  else if(pdg == 221 )
  {
    coeff = llMesCofEta; // Eta
  }
  else if(pdg == 331 )
  {
    coeff = llMesCofEtaP; // Eta'
  } 
  return coeff;
}

////////////////////////////////////////////////////////////////////////////////
// < pEDM v1.0 >
// Modified for pEDM from here
// 17-May-16 by H. Jeong
////////////////////////////////////////////////////////////////////////////////
// Additional functions for pEDM
G4double G4HadronElastic::transformFunction(G4double x,
  G4double A, G4double r, G4double polPerp)
{
  return x + polPerp * A * std::sin(x) - CLHEP::twopi * r;
}

G4double par[3] = {-8.12966e-5, 7.65177e-3, -2.81392e-1};
G4double slope = 4.35264e-2, intersept = 2.22743e-1, xx1 = 17.7;
G4double yy1 = slope*xx1 + intersept;

G4double G4HadronElastic::getAP(G4double cost)
{
  G4double t = std::acos(cost)*360.0/CLHEP::twopi;
  if ( t <= 17.7)
    return slope*std::acos(cost)*360.0/CLHEP::twopi + intersept;
  else if ( (t>17.7) && (t<26.5877) )
  {
    G4double AP = par[0]*(t-xx1)*(t-xx1)*(t*t+2*t*xx1+3*xx1*xx1)
                  + par[1]*(t-xx1)*(t-xx1)*(t+2*xx1)
                  + par[2]*(t-xx1)*(t-xx1) + slope*(t-xx1) + yy1;
    if ( AP <= 1 ) return AP;
      else return 1.0;
  }
   else return 0.0;
}

G4double G4HadronElastic::getAPFromT(G4double t)
{
  double AP = APSpline(t);
  if ( AP > 1.0 ) return 1.0;
  else if ( AP < - 1.0) return -1.0;
  else return AP;
}

G4double G4HadronElastic::getBetaFromT(G4double t, G4double eKinetic)
{
  G4double beta;
  if ( eKinetic < 155.0 ) beta = betaSpline150(t);
  else if ( eKinetic < 165.0 ) beta = betaSpline160(t);
  else if ( eKinetic < 175.0 ) beta = betaSpline170(t);
  else if ( eKinetic < 185.0 ) beta = betaSpline180(t);
  else if ( eKinetic < 195.0 ) beta = betaSpline190(t);
  else if ( eKinetic < 205.0 ) beta = betaSpline200(t);
  else if ( eKinetic < 215.0 ) beta = betaSpline210(t);
  else if ( eKinetic < 225.0 ) beta = betaSpline220(t);
  else if ( eKinetic < 235.0 ) beta = betaSpline230(t);
  else beta = betaSpline240(t);
  return beta;
}

G4double G4HadronElastic::BisectionPhiWithT(G4double t, G4double polPerp)
{
  if ( t < 0.0 ) t = 0;
  // Can choose how to get AP
  G4double AP = getAPFromT(t);
  G4double r = G4UniformRand();
  G4double x1 = 0.0, x2 = CLHEP::twopi;
  G4double x = 0.5 * (x1 + x2);
  G4double ppp = polPerp;
  G4double xc;
  while ( std::abs(transformFunction(x, AP, r, ppp)) >= 1e-15 )
  {
    xc = 0.5 * ( x1 + x2 );
    if ( transformFunction(x1, AP, r, ppp) == 0.0 ) { x = x1; break; }
    if ( transformFunction(x2, AP, r, ppp) == 0.0 ) { x = x2; break; }
    if ( transformFunction(xc, AP, r, ppp) == 0.0 ) { x = xc; break; }
    if ( (transformFunction(x1, AP, r, ppp)*transformFunction(x2, AP, r, ppp) < 0.0)
      && (transformFunction(xc, AP, r, ppp) < 0.0) ) { x1 = xc, x = xc; }
    if ( (transformFunction(x1, AP, r, ppp)*transformFunction(x2, AP, r, ppp) < 0.0)
      && (transformFunction(xc, AP, r, ppp) > 0.0) ) { x2 = xc; x = xc; }
  }
  return x;
}

// Splines
double G4HadronElastic::APSpline(double x) {
   const int fNp = 46, fKstep = 0;
   const double fDelta = -1, fXmin = 0, fXmax = 1.31943e+06;
   const double fX[46] = { 0, 4685.33, 8260.23, 12834.3, 18400.9,
                        24951.6, 32477, 40965.8, 50405.6, 60782.5,
                        71547.8, 83710.8, 96761.7, 110681, 125449,
                        140317, 156679, 173025, 190889, 208622,
                        227886, 246906, 267463, 287662, 308396,
                        329638, 351358, 374594, 397203, 442434,
                        490059, 537568, 585837, 634618, 683663,
                        731526, 776794, 825182, 869399, 926321,
                        978168, 1.03012e+06, 1.07858e+06, 1.12337e+06, 1.23826e+06,
                        1.31943e+06 };
   const double fY[46] = { 0, 0.507, 0.613, 0.707, 0.807,
                        0.893, 0.993, 0.968, 0.782, 0.451,
                        0.044, -0.336, -0.425, -0.13, 0.435,
                        0.868, 0.96, 0.843, 0.653, 0.451,
                        0.294, 0.087, -0.154, -0.348, -0.579,
                        -0.87, -0.996, -0.922, -0.65, -0.073,
                        0.29, 0.427, 0.436, 0.54, 0.621,
                        0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0,
                        0 };
   const double fB[46] = { 0.00018479, 4.97923e-05, 2.00837e-05, 2.01969e-05, 1.49287e-05,
                        1.40453e-05, 7.59933e-06, -1.20526e-05, -2.6476e-05, -3.62632e-05,
                        -3.74608e-05, -2.1331e-05, 7.60634e-06, 3.3072e-05, 3.78413e-05,
                        1.77495e-05, -3.19116e-06, -9.59402e-06, -1.17299e-05, -9.57709e-06,
                        -8.86862e-06, -1.20845e-05, -1.05826e-05, -9.53473e-06, -1.34926e-05,
                        -1.09412e-05, -1.3737e-06, 8.32765e-06, 1.38639e-05, 1.06157e-05,
                        4.97763e-06, 9.85261e-07, 3.18984e-07, 4.68271e-06, -7.71709e-06,
                        -8.16235e-06, 2.24932e-06, -5.82356e-07, 1.73602e-07, -4.44908e-08,
                        1.1905e-08, -3.08542e-09, 8.22593e-10, -3.14182e-10, 1.30072e-10,
                        -2.21968e-10 };
   const double fC[46] = { -2.02212e-08, -8.59176e-09, 2.81446e-10, -2.56703e-10, -6.89695e-10,
                        5.54837e-10, -1.4114e-09, -9.0363e-10, -6.24302e-10, -3.18876e-10,
                        2.07638e-10, 1.1185e-09, 1.09877e-09, 7.30706e-10, -4.07752e-10,
                        -9.43611e-10, -3.36242e-10, -5.54577e-11, -6.41118e-11, 1.85517e-10,
                        -1.48741e-10, -2.03432e-11, 9.34035e-11, -4.15275e-11, -1.49356e-10,
                        2.69469e-10, 1.71018e-10, 2.46492e-10, -1.62062e-12, -7.01932e-11,
                        -4.81907e-11, -3.58439e-11, 2.20404e-11, 6.74138e-11, -3.20239e-10,
                        3.10936e-10, -8.09343e-11, 2.24146e-11, -5.31809e-12, 1.48662e-12,
                        -3.98884e-13, 1.10365e-13, -2.97201e-14, 4.33697e-15, -4.70323e-16,
                        81172 };
   const double fD[46] = { 8.27362e-13, 8.27362e-13, -3.92173e-14, -2.59283e-14, 6.33275e-14,
                        -8.7094e-14, 1.99388e-14, 9.86345e-15, 9.81113e-15, 1.63028e-14,
                        2.49626e-14, -5.03868e-16, -8.81396e-15, -2.56967e-14, -1.20138e-14,
                        1.23738e-14, 5.72572e-15, -1.61484e-16, 4.69239e-15, -5.7836e-15,
                        2.25028e-15, 1.84437e-15, -2.22669e-15, -1.73351e-15, 6.57245e-15,
                        -1.5109e-15, 1.08271e-15, -3.65805e-15, -5.0535e-16, 1.53997e-16,
                        8.66285e-17, 3.99736e-16, 3.10044e-16, -2.63467e-15, 4.39574e-15,
                        -2.88557e-15, 7.11941e-16, -2.09064e-16, 3.98485e-17, -1.21223e-17,
                        3.26717e-18, -9.6358e-19, 2.53488e-19, -1.3947e-20, -1.3947e-20,
                        55532.1 };
   int klow=0;
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((x-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=x-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
}

double G4HadronElastic::betaSpline150(double x) {
   const int fNp = 81, fKstep = 0;
   const double fDelta = -1, fXmin = 0, fXmax = 418272;
   const double fX[81] = { 0, 77.1068, 308.427, 693.962, 1232.32,
                        1925.94, 2773.77, 3773.39, 4926.52, 6233.19,
                        7689.91, 9299.46, 11061.8, 12972.6, 15035.4,
                        17245.2, 19606.5, 22119.2, 24771, 27579.4,
                        30524.7, 33619.5, 36856.1, 40241.4, 43758.9,
                        47423.8, 51227.1, 55167.8, 59244.9, 63467.4,
                        67814.6, 72295.2, 76908.2, 81663.9, 86539,
                        91543.4, 96664, 101924, 107297, 112809,
                        118418, 124162, 130015, 135987, 142064,
                        148259, 154555, 160966, 167474, 174096,
                        180812, 187621, 194539, 201546, 208659,
                        215841, 223125, 230493, 237941, 245469,
                        253095, 260776, 268532, 276380, 284278,
                        292245, 300300, 308400, 316562, 324786,
                        333047, 341389, 349764, 358170, 366652,
                        375160, 383717, 392296, 400920, 409588,
                        418272 };
   const double fY[81] = { 0, 0.06, 0.62, 2.66, 7.07,
                        12.41, 15.83, 17.18, 17.58, 17.74,
                        17.92, 18.24, 18.73, 19.44, 20.42,
                        21.79, 23.71, 26.49, 30.64, 37.04,
                        47.12, 62.61, 83.07, 103.05, 117.39,
                        125.89, 130.29, 131.93, 131.6, 129.72,
                        126.44, 121.77, 115.58, 107.63, 97.57,
                        85.1, 70.2, 53.45, 36.08, 19.35,
                        3.91, -10.24, -23.5, -36.29, -48.97,
                        -61.78, -74.79, -87.89, -100.87, -113.45,
                        -125.4, -136.58, -146.96, -156.58, -165.54,
                        -173.94, -181.89, -189.49, -196.81, -203.91,
                        -210.85, -217.69, -224.47, -231.25, -238.08,
                        -245.03, -252.19, -259.65, -267.53, -275.97,
                        -285.19, -295.44, -307.08, -320.54, -336.19,
                        -354.1, -373.56, -393.19, -411.62, -428.25,
                        -443.18 };
   const double fB[81] = { 0.000341322, 0.00120843, 0.00357454, 0.00690024, 0.00867813,
                        0.00614989, 0.00239119, 0.000647738, 0.000166088, 0.000103979,
                        0.00015527, 0.000237477, 0.000322151, 0.000418767, 0.000539176,
                        0.000704495, 0.000936038, 0.00130102, 0.00186187, 0.00276366,
                        0.00416016, 0.00583692, 0.00643474, 0.0051063, 0.00311756,
                        0.00164886, 0.000736096, 0.000138556, -0.000276492, -0.000605006,
                        -0.000898252, -0.00118806, -0.00149752, -0.00185765, -0.00227266,
                        -0.00271206, -0.00307873, -0.00325203, -0.00316359, -0.00290087,
                        -0.00260048, -0.00234886, -0.00219228, -0.00210484, -0.00207276,
                        -0.00206777, -0.00205927, -0.00202503, -0.00195382, -0.00184295,
                        -0.00171283, -0.00157033, -0.0014344, -0.00131304, -0.00121173,
                        -0.00112828, -0.0010588, -0.00100563, -0.000961835, -0.00092487,
                        -0.000898727, -0.000881989, -0.000867414, -0.000862724, -0.00086755,
                        -0.000878375, -0.000902516, -0.000941298, -0.000991583, -0.00106711,
                        -0.00116671, -0.00129994, -0.00148895, -0.00171727, -0.0019803,
                        -0.00221159, -0.00231031, -0.00223378, -0.00203017, -0.00181312,
                        -0.00163153 };
   const double fC[81] = { 5.74985e-06, 5.49565e-06, 4.73306e-06, 3.89316e-06, -5.90762e-07,
                        -3.05424e-06, -1.37909e-06, -3.65026e-07, -5.2659e-08, 5.12602e-09,
                        3.00841e-08, 2.09902e-08, 2.70552e-08, 2.35097e-08, 3.48604e-08,
                        3.9951e-08, 5.81077e-08, 8.71459e-08, 1.24357e-07, 1.96747e-07,
                        2.77385e-07, 2.64422e-07, -7.97146e-08, -3.127e-07, -2.5268e-07,
                        -1.4807e-07, -9.19234e-08, -5.97087e-08, -4.20907e-08, -3.5711e-08,
                        -3.17447e-08, -3.29366e-08, -3.41474e-08, -4.15792e-08, -4.35488e-08,
                        -4.4254e-08, -2.73515e-08, -5.59691e-09, 2.20553e-08, 2.5614e-08,
                        2.79431e-08, 1.58565e-08, 1.08977e-08, 3.7437e-09, 1.53536e-09,
                        -7.29636e-10, 2.07969e-09, 3.26094e-09, 7.67979e-09, 9.06323e-09,
                        1.03114e-08, 1.06171e-08, 9.03205e-09, 8.28741e-09, 5.95542e-09,
                        5.66463e-09, 3.87346e-09, 3.3431e-09, 2.53638e-09, 2.37402e-09,
                        1.05438e-09, 1.12462e-09, 7.54745e-10, -1.57266e-10, -4.53716e-10,
                        -9.05007e-10, -2.09189e-09, -2.69654e-09, -3.46402e-09, -5.71959e-09,
                        -6.3364e-09, -9.63494e-09, -1.29338e-08, -1.42282e-08, -1.67839e-08,
                        -1.03995e-08, -1.13727e-09, 1.00566e-08, 1.35526e-08, 1.14895e-08,
                        8683.6 };
   const double fD[81] = { -1.09889e-09, -1.09889e-09, -7.26184e-10, -2.77627e-09, -1.18388e-09,
                        6.58605e-10, 3.3815e-10, 9.02947e-11, 1.47411e-11, 5.71101e-12,
                        -1.8833e-12, 1.14712e-12, -6.18531e-13, 1.83413e-12, 7.67876e-13,
                        2.56313e-12, 3.85215e-12, 4.67757e-12, 8.59214e-12, 9.12588e-12,
                        -1.39624e-12, -3.54423e-11, -2.29408e-11, 5.68771e-12, 9.51469e-12,
                        4.92084e-12, 2.72495e-12, 1.44039e-12, 5.0363e-13, 3.04129e-13,
                        -8.86738e-14, -8.74913e-14, -5.20911e-13, -1.34672e-13, -4.6969e-14,
                        1.1003e-12, 1.37863e-12, 1.71537e-12, 2.15227e-13, 1.38423e-13,
                        -7.01309e-13, -2.82435e-13, -3.99267e-13, -1.21131e-13, -1.2188e-13,
                        1.48741e-13, 6.14192e-14, 2.26301e-13, 6.96455e-14, 6.19506e-14,
                        1.49641e-14, -7.63793e-14, -3.54221e-14, -1.09279e-13, -1.34969e-14,
                        -8.19629e-14, -2.39956e-14, -3.61016e-14, -7.18893e-15, -5.76866e-14,
                        3.04807e-15, -1.58975e-14, -3.87347e-14, -1.25112e-14, -1.8882e-14,
                        -4.91136e-14, -2.48855e-14, -3.13422e-14, -9.14214e-14, -2.48877e-14,
                        -1.31806e-13, -1.31296e-13, -5.13312e-14, -1.00442e-13, 2.50125e-13,
                        3.60803e-13, 4.34917e-13, 1.35121e-13, -7.93439e-14, -7.93439e-14,
                        4023.12 };
   int klow=0;
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((x-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=x-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
}

double G4HadronElastic::betaSpline160(double x) {
   const int fNp = 81, fKstep = 0;
   const double fDelta = -1, fXmin = 0, fXmax = 447690;
   const double fX[81] = { 0, 82.3927, 330.288, 742.609, 1319.72,
                        2061.61, 2968.29, 4037.24, 5273.13, 6670.58,
                        8232.1, 9953.76, 11838.8, 13887.1, 16093.8,
                        18457.8, 20989.5, 23671.6, 26515, 29519.6,
                        32671.1, 35982.5, 39453.6, 43068.2, 46841.2,
                        50763.7, 54834.6, 59053, 63417.7, 67927.7,
                        72582, 77379.6, 82319.3, 87400.2, 92621.3,
                        97981.4, 103467, 109089, 114846, 120738,
                        126748, 132891, 139164, 145551, 152052,
                        158678, 165430, 172289, 179255, 186341,
                        193530, 200820, 208227, 215731, 223330,
                        231024, 238808, 246703, 254685, 262735,
                        270889, 279125, 287422, 295819, 304271,
                        312821, 321421, 330093, 338811, 347620,
                        356471, 365386, 374363, 383375, 392445,
                        401545, 410700, 419907, 429138, 438391,
                        447690 };
   const double fY[81] = { 0, 0.07, 0.73, 3.1, 7.95,
                        13.06, 15.79, 16.66, 16.86, 16.97,
                        17.19, 17.59, 18.21, 19.09, 20.34,
                        22.11, 24.65, 28.41, 34.18, 43.28,
                        57.55, 77.53, 98.68, 114.72, 124.5,
                        129.7, 131.86, 131.9, 130.26, 127.15,
                        122.55, 116.3, 108.07, 97.41, 83.9,
                        67.51, 49.12, 30.45, 12.99, -2.79,
                        -17.23, -30.94, -44.49, -58.25, -72.39,
                        -86.78, -101.05, -114.75, -127.51, -139.15,
                        -149.68, -159.22, -167.91, -175.94, -183.45,
                        -190.55, -197.36, -203.96, -210.41, -216.77,
                        -223.12, -229.52, -236.04, -242.78, -249.85,
                        -257.39, -265.59, -274.71, -285.11, -297.28,
                        -311.89, -329.54, -350.25, -372.65, -394.33,
                        -413.64, -430.38, -445.14, -458.54, -471.01,
                        -482.77 };
   const double fB[81] = { 0.000367351, 0.00132434, 0.00393272, 0.00737825, 0.00843241,
                        0.00504411, 0.00156342, 0.000331095, 7.30508e-05, 9.94088e-05,
                        0.000183928, 0.000279973, 0.00037547, 0.000489886, 0.000648046,
                        0.000857317, 0.00117053, 0.00166487, 0.00244574, 0.00369958,
                        0.00538879, 0.00638519, 0.00544826, 0.00345492, 0.00186297,
                        0.000870186, 0.000239259, -0.000198206, -0.000538769, -0.000837716,
                        -0.00113943, -0.00147314, -0.00186922, -0.00233627, -0.00283585,
                        -0.0032499, -0.00339142, -0.00320528, -0.00285233, -0.00252319,
                        -0.00229994, -0.00218197, -0.00214902, -0.00216458, -0.00218093,
                        -0.00215257, -0.00206482, -0.00192135, -0.00173835, -0.00155029,
                        -0.00138358, -0.00123662, -0.00111692, -0.00102646, -0.000952652,
                        -0.000896658, -0.000853988, -0.000820127, -0.000797762, -0.000783263,
                        -0.00077598, -0.000780128, -0.000792034, -0.000817299, -0.000856096,
                        -0.000913289, -0.000996617, -0.00111513, -0.00127611, -0.00150287,
                        -0.0018084, -0.00215387, -0.00243592, -0.00248007, -0.00227183,
                        -0.00196945, -0.00170239, -0.00151697, -0.00139433, -0.00130361,
                        -0.00122836 };
   const double fC[81] = { 5.94382e-06, 5.67118e-06, 4.85092e-06, 3.50551e-06, -1.67889e-06,
                        -2.88821e-06, -9.50745e-07, -2.0209e-07, -6.70248e-09, 2.5564e-08,
                        2.85625e-08, 2.72238e-08, 2.34372e-08, 3.24203e-08, 3.92516e-08,
                        4.9273e-08, 7.44452e-08, 1.09864e-07, 1.64766e-07, 2.52541e-07,
                        2.83455e-07, 1.7449e-08, -2.87366e-07, -2.64107e-07, -1.57826e-07,
                        -9.5276e-08, -5.97068e-08, -4.39983e-08, -3.4028e-08, -3.22569e-08,
                        -3.25679e-08, -3.69916e-08, -4.31899e-08, -4.87329e-08, -4.6953e-08,
                        -3.02934e-08, 4.49412e-09, 2.86162e-08, 3.26908e-08, 2.31746e-08,
                        1.39676e-08, 5.23738e-09, 1.49798e-11, -2.45078e-09, -6.36775e-11,
                        4.34244e-09, 8.65486e-09, 1.22599e-08, 1.40132e-08, 1.25239e-08,
                        1.06662e-08, 9.49348e-09, 6.66752e-09, 5.38688e-09, 4.32559e-09,
                        2.95281e-09, 2.52826e-09, 1.76112e-09, 1.04043e-09, 7.60926e-10,
                        1.32252e-10, -6.3584e-10, -7.99232e-10, -2.20969e-09, -2.38035e-09,
                        -4.30943e-09, -5.37921e-09, -8.28754e-09, -1.01763e-08, -1.55658e-08,
                        -1.8953e-08, -1.98005e-08, -1.16199e-08, 6.72122e-09, 1.6238e-08,
                        1.69878e-08, 1.2184e-08, 7.95579e-09, 5.32912e-09, 4.47541e-09,
                        9298.98 };
   const double fD[81] = { -1.10298e-09, -1.10298e-09, -1.08767e-09, -2.99448e-09, -5.43349e-10,
                        7.12294e-10, 2.33454e-10, 5.26981e-11, 7.69649e-12, 6.40081e-13,
                        -2.59174e-13, -6.69602e-13, 1.46183e-12, 1.0319e-12, 1.41305e-12,
                        3.31435e-12, 4.40176e-12, 6.43633e-12, 9.73791e-12, 3.26974e-12,
                        -2.67773e-11, -2.92712e-11, 2.14491e-12, 9.38971e-12, 5.31547e-12,
                        2.91244e-12, 1.24129e-12, 7.61428e-13, 1.30906e-13, -2.22752e-14,
                        -3.07357e-13, -4.18266e-13, -3.63649e-13, 1.13639e-13, 1.03601e-12,
                        2.11391e-12, 1.43025e-12, 2.35913e-13, -5.38407e-13, -5.10582e-13,
                        -4.73766e-13, -2.77505e-13, -1.28675e-13, 1.22413e-13, 2.21645e-13,
                        2.12906e-13, 1.75185e-13, 8.39005e-14, -7.0052e-14, -8.61335e-14,
                        -5.36244e-14, -1.27182e-13, -5.68869e-14, -4.65505e-14, -5.94802e-14,
                        -1.81784e-14, -3.23921e-14, -3.00935e-14, -1.15747e-14, -2.57003e-14,
                        -3.10837e-14, -6.56454e-15, -5.59925e-14, -6.73023e-15, -7.52137e-14,
                        -4.14615e-14, -1.11791e-13, -7.2212e-14, -2.03946e-13, -1.27562e-13,
                        -3.1688e-14, 3.03776e-13, 6.78389e-13, 3.49758e-13, 2.74638e-14,
                        -1.74914e-13, -1.53086e-13, -9.48486e-14, -3.07538e-14, -3.07538e-14,
                        4293.93 };
   int klow=0;
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((x-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=x-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
}

double G4HadronElastic::betaSpline170(double x) {
   const int fNp = 81, fKstep = 0;
   const double fDelta = -1, fXmin = 0, fXmax = 477290;
   const double fX[81] = { 0, 87.8537, 352.155, 791.794, 1407.14,
                        2198.19, 3164.95, 4304.83, 5622.64, 7112.82,
                        8774.28, 10614, 12619.9, 14804, 17157.2,
                        19683.9, 22372.3, 25239, 28271.8, 31469.7,
                        34831.4, 38363.8, 42058.6, 45914.8, 49940.1,
                        54116.1, 58459.9, 62951.4, 67609.2, 72422.6,
                        77379.6, 82500.6, 87762.3, 93186.6, 98748.8,
                        104459, 110317, 116308, 122444, 128723,
                        135130, 141678, 148365, 155176, 162108,
                        169175, 176359, 183676, 191107, 198668,
                        206322, 214102, 221990, 230000, 238095,
                        246311, 254606, 263018, 271526, 280106,
                        288799, 297582, 306430, 315386, 324404,
                        333503, 342659, 351915, 361223, 370605,
                        380059, 389557, 399123, 408730, 418399,
                        428104, 437868, 447663, 457513, 467390,
                        477290 };
   const double fY[81] = { 0, 0.08, 0.85, 3.58, 8.76,
                        13.46, 15.51, 16, 16.07, 16.17,
                        16.43, 16.92, 17.67, 18.76, 20.31,
                        22.53, 25.82, 30.83, 38.73, 51.37,
                        70.28, 92.58, 111.04, 122.73, 129.14,
                        132.07, 132.63, 131.39, 128.58, 124.2,
                        118.05, 109.73, 98.65, 84.19, 66.27,
                        46.18, 26.27, 8.24, -7.75, -22.44,
                        -36.69, -51.2, -66.38, -82.2, -98.15,
                        -113.46, -127.49, -139.95, -150.86, -160.46,
                        -169, -176.72, -183.83, -190.49, -196.83,
                        -202.94, -208.92, -214.83, -220.74, -226.74,
                        -232.92, -239.39, -246.3, -253.85, -262.31,
                        -272.09, -283.77, -298.25, -316.6, -339.32,
                        -364.73, -389.05, -409.88, -427.31, -442.38,
                        -455.95, -468.54, -480.36, -491.43, -501.7,
                        -511.1 };
   const double fB[81] = { 0.000378301, 0.00143476, 0.00431815, 0.00777858, 0.00795243,
                        0.00393862, 0.00092191, 0.000128479, 3.04255e-05, 0.000106571,
                        0.000209221, 0.000319074, 0.000430858, 0.000571111, 0.000751995,
                        0.00102547, 0.0014426, 0.00210257, 0.00318116, 0.00480993,
                        0.00626383, 0.00592563, 0.00399947, 0.00220602, 0.00107816,
                        0.00037904, -9.19074e-05, -0.000446343, -0.000755021, -0.00106843,
                        -0.00141837, -0.00184812, -0.00237395, -0.00296212, -0.00343248,
                        -0.00352289, -0.00322813, -0.0027944, -0.00244643, -0.00225992,
                        -0.00220584, -0.00223712, -0.00230277, -0.00232945, -0.00225117,
                        -0.00206966, -0.00182975, -0.00158079, -0.00136195, -0.00118657,
                        -0.00104902, -0.000942449, -0.000862945, -0.000804888, -0.000761449,
                        -0.000730395, -0.000710958, -0.000696316, -0.000695461, -0.000703568,
                        -0.000720786, -0.000756029, -0.000807743, -0.000884424, -0.000997754,
                        -0.00116296, -0.00140144, -0.00174997, -0.00220456, -0.00260891,
                        -0.00268735, -0.00238686, -0.00198035, -0.00166829, -0.00146612,
                        -0.00133772, -0.00124624, -0.00116586, -0.00108234, -0.000996072,
                        -0.000901605 };
   const double fC[81] = { 6.15176e-06, 5.87342e-06, 5.03608e-06, 2.83499e-06, -2.55247e-06,
                        -2.52152e-06, -5.98913e-07, -9.71529e-08, 2.27459e-08, 2.83519e-08,
                        3.34313e-08, 2.62793e-08, 2.945e-08, 3.47641e-08, 4.21026e-08,
                        6.61329e-08, 8.90255e-08, 1.41196e-07, 2.14444e-07, 2.9489e-07,
                        1.37591e-07, -2.33336e-07, -2.87977e-07, -1.77107e-07, -1.03086e-07,
                        -6.43245e-08, -4.40959e-08, -3.48152e-08, -3.14561e-08, -3.36551e-08,
                        -3.69416e-08, -4.69763e-08, -5.29594e-08, -5.54733e-08, -2.90905e-08,
                        1.32594e-08, 3.7058e-08, 3.53391e-08, 2.13753e-08, 8.32662e-09,
                        1.15194e-10, -4.89218e-09, -4.92477e-09, 1.00686e-09, 1.02857e-08,
                        1.53995e-08, 1.79924e-08, 1.60326e-08, 1.34166e-08, 9.77867e-09,
                        8.19244e-09, 5.50547e-09, 4.57487e-09, 2.67248e-09, 2.69361e-09,
                        1.0863e-09, 1.25691e-09, 4.83637e-10, -3.8306e-10, -5.61729e-10,
                        -1.419e-09, -2.59389e-09, -3.25046e-09, -5.31142e-09, -7.25624e-09,
                        -1.08996e-08, -1.51479e-08, -2.2506e-08, -2.63308e-08, -1.6768e-08,
                        8.47114e-09, 2.31631e-08, 1.93324e-08, 1.31527e-08, 7.75492e-09,
                        5.47481e-09, 3.89402e-09, 4.31275e-09, 4.16605e-09, 4.56909e-09,
                        9899.89 };
   const double fD[81] = { -1.05605e-09, -1.05605e-09, -1.66886e-09, -2.91839e-09, 1.30421e-11,
                        6.62904e-10, 1.46729e-10, 3.03279e-11, 1.25399e-12, 1.01906e-12,
                        -1.29583e-12, 5.26906e-13, 8.11012e-13, 1.03951e-12, 3.17018e-12,
                        2.83844e-12, 6.06625e-12, 8.05061e-12, 8.38549e-12, -1.55968e-11,
                        -3.5003e-11, -4.92958e-12, 9.58375e-12, 6.12966e-12, 3.09394e-12,
                        1.55233e-12, 6.88744e-13, 2.40392e-13, -1.52282e-13, -2.21002e-13,
                        -6.53171e-13, -3.79035e-13, -1.54482e-13, 1.58108e-12, 2.47199e-12,
                        1.35418e-12, -9.56407e-14, -7.58632e-13, -6.92708e-13, -4.27208e-13,
                        -2.54914e-13, -1.62449e-15, 2.90313e-13, 4.46205e-13, 2.41211e-13,
                        1.203e-13, -8.92836e-14, -1.1734e-13, -1.60383e-13, -6.90814e-14,
                        -1.15119e-13, -3.93298e-14, -7.91593e-14, 8.70127e-16, -6.52153e-14,
                        6.85594e-15, -3.0641e-14, -3.39571e-14, -6.94088e-15, -3.28729e-14,
                        -4.45919e-14, -2.47336e-14, -7.67062e-14, -7.18899e-14, -1.33468e-13,
                        -1.54671e-13, -2.64976e-13, -1.36968e-13, 3.39762e-13, 8.89941e-13,
                        5.15582e-13, -1.33484e-13, -2.14431e-13, -1.86076e-13, -7.83129e-14,
                        -5.39658e-14, 1.42495e-14, -4.96445e-15, 1.36033e-14, 1.36033e-14,
                        4583.66 };
   int klow=0;
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((x-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=x-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
}

double G4HadronElastic::betaSpline180(double x) {
   const int fNp = 81, fKstep = 0;
   const double fDelta = -1, fXmin = 0, fXmax = 507106;
   const double fX[81] = { 0, 93.49, 373.96, 841.41, 1494.31,
                        2335.34, 3361.06, 4573, 5971.15, 7555.53,
                        9322.31, 11274.5, 13407.7, 15725.5, 18228,
                        20909.5, 23768.9, 26811.4, 30037.1, 33431.6,
                        37007.8, 40757.6, 44680.1, 48782.7, 53056.4,
                        57490.6, 62102.8, 66882.6, 71829.1, 76941.1,
                        82217.4, 87645.4, 93246.9, 98997, 104906,
                        110974, 117198, 123565, 130086, 136760,
                        143571, 150517, 157626, 164852, 172224,
                        179740, 187382, 195148, 203037, 211064,
                        219209, 227471, 235848, 244356, 252956,
                        261684, 270499, 279438, 288460, 297603,
                        306824, 316140, 325551, 335055, 344648,
                        354307, 364051, 373880, 383766, 393731,
                        403774, 413868, 424034, 434246, 444500,
                        454821, 465207, 475601, 486054, 496566,
                        507106 };
   const double fY[81] = { 0, 0.09, 0.98, 4.06, 9.47,
                        13.62, 15.06, 15.25, 15.22, 15.32,
                        15.63, 16.2, 17.09, 18.38, 20.26,
                        23.02, 27.21, 33.81, 44.53, 61.6,
                        84.48, 105.97, 120.41, 128.52, 132.49,
                        133.75, 133.04, 130.66, 126.65, 120.78,
                        112.6, 101.33, 86.1, 66.67, 44.8,
                        23.71, 5.28, -10.77, -25.61, -40.44,
                        -56.15, -73.17, -91.1, -108.76, -124.89,
                        -138.82, -150.58, -160.52, -169.07, -176.59,
                        -183.37, -189.62, -195.51, -201.16, -206.67,
                        -212.13, -217.61, -223.2, -228.99, -235.11,
                        -241.74, -249.11, -257.61, -267.83, -280.69,
                        -297.63, -320.26, -348.12, -376.36, -400.23,
                        -419.42, -435.45, -449.61, -462.59, -474.72,
                        -486.04, -496.49, -505.99, -514.51, -522.07,
                        -528.74 };
   const double fB[81] = { 0.000361975, 0.00155135, 0.00468702, 0.00805, 0.00733322,
                        0.0028957, 0.00046761, -1.50161e-05, 3.43613e-06, 0.000118762,
                        0.000231306, 0.000352542, 0.000480528, 0.000641343, 0.000870048,
                        0.00121088, 0.00175682, 0.0026446, 0.00410993, 0.005925,
                        0.00642338, 0.00478841, 0.00271461, 0.00136559, 0.00056041,
                        4.33274e-05, -0.000334838, -0.000653583, -0.000971464, -0.00133194,
                        -0.00179042, -0.00237864, -0.00307047, -0.00362383, -0.0036691,
                        -0.00323393, -0.00271292, -0.00236486, -0.00222076, -0.00224626,
                        -0.00237929, -0.0025059, -0.00251295, -0.00234041, -0.00202473,
                        -0.00168885, -0.00139916, -0.00117261, -0.001003, -0.000878774,
                        -0.000790339, -0.000726942, -0.000681144, -0.000650414, -0.000631619,
                        -0.000622086, -0.000621803, -0.000631654, -0.000652775, -0.000690503,
                        -0.000750191, -0.00083872, -0.000977143, -0.00118773, -0.00151837,
                        -0.00202043, -0.00262589, -0.0029443, -0.00267046, -0.00213169,
                        -0.00172289, -0.00147437, -0.00132325, -0.00122453, -0.00114074,
                        -0.00105193, -0.0009609, -0.000865151, -0.000765874, -0.0006743,
                        -0.00059316 };
   const double fC[81] = { 6.55372e-06, 6.16825e-06, 5.01182e-06, 2.18248e-06, -3.28031e-06,
                        -1.99599e-06, -3.7121e-07, -2.70173e-08, 4.02149e-08, 3.25745e-08,
                        3.11257e-08, 3.09753e-08, 2.90242e-08, 4.03578e-08, 5.10322e-08,
                        7.60707e-08, 1.14861e-07, 1.76927e-07, 2.77336e-07, 2.57377e-07,
                        -1.18015e-07, -3.17992e-07, -2.10709e-07, -1.18112e-07, -7.02908e-08,
                        -4.63216e-08, -3.5672e-08, -3.10129e-08, -3.32513e-08, -3.72651e-08,
                        -4.9627e-08, -5.87424e-08, -6.47657e-08, -3.14687e-08, 2.38085e-08,
                        4.79124e-08, 3.57899e-08, 1.8879e-08, 3.22005e-09, -7.04203e-09,
                        -1.24887e-08, -5.73957e-09, 4.74787e-09, 1.91292e-08, 2.36965e-08,
                        2.09923e-08, 1.69138e-08, 1.22574e-08, 9.2425e-09, 6.23424e-09,
                        4.62272e-09, 3.05064e-09, 2.41675e-09, 1.19503e-09, 9.90502e-10,
                        1.01734e-10, -6.96796e-11, -1.03216e-09, -1.30905e-09, -2.8173e-09,
                        -3.65609e-09, -5.8461e-09, -8.86263e-09, -1.32966e-08, -2.11673e-08,
                        -3.08135e-08, -3.13201e-08, -1.07708e-09, 2.87775e-08, 2.52858e-08,
                        1.54189e-08, 9.2038e-09, 5.66096e-09, 4.00626e-09, 4.16461e-09,
                        4.44046e-09, 4.32422e-09, 4.88795e-09, 4.60884e-09, 4.10286e-09,
                        10539.8 };
   const double fD[81] = { -1.37439e-09, -1.37439e-09, -2.01757e-09, -2.78897e-09, 5.09027e-10,
                        5.28014e-10, 9.46673e-11, 1.60288e-11, -1.60744e-12, -2.73347e-13,
                        -2.56777e-14, -3.04896e-13, 1.62993e-12, 1.42182e-12, 3.1125e-12,
                        4.52205e-12, 6.7997e-12, 1.03759e-11, -1.95989e-12, -3.49904e-11,
                        -1.77763e-11, 9.11702e-12, 7.52338e-12, 3.72988e-12, 1.80184e-12,
                        7.69681e-13, 3.24912e-13, -1.50846e-13, -2.61721e-13, -7.80962e-13,
                        -5.59782e-13, -3.58434e-13, 1.93022e-12, 3.11806e-12, 1.32421e-12,
                        -6.49178e-13, -8.85387e-13, -8.0046e-13, -5.12549e-13, -2.66568e-13,
                        3.23897e-13, 4.91685e-13, 6.63418e-13, 2.06534e-13, -1.19934e-13,
                        -1.7789e-13, -1.9985e-13, -1.27396e-13, -1.24926e-13, -6.5948e-14,
                        -6.34262e-14, -2.52249e-14, -4.78631e-14, -7.92759e-15, -3.39425e-14,
                        -6.48215e-15, -3.58875e-14, -1.02307e-14, -5.49855e-14, -3.03237e-14,
                        -7.83546e-14, -1.06844e-13, -1.55523e-13, -2.73468e-13, -3.32904e-13,
                        -1.73301e-14, 1.0257e-12, 1.00663e-12, -1.16796e-13, -3.27486e-13,
                        -2.05255e-13, -1.1616e-13, -5.40124e-14, 5.14736e-15, 8.90921e-15,
                        -3.7311e-15, 1.80792e-14, -8.90011e-15, -1.60449e-14, -1.60449e-14,
                        4879.06 };
   int klow=0;
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((x-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=x-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
}

double G4HadronElastic::betaSpline190(double x) {
   const int fNp = 81, fKstep = 0;
   const double fDelta = -1, fXmin = 0, fXmax = 537073;
   const double fX[81] = { 0, 98.9086, 395.634, 890.177, 1582.54,
                        2472.71, 3560.71, 4843.77, 6323.87, 8001,
                        9875.17, 11942, 14200.5, 16659.5, 19303.7,
                        22148.6, 25176.4, 28398.1, 31806.6, 35408.3,
                        39195.4, 43166.5, 47329.3, 51665.7, 56192.2,
                        60889.2, 65774.7, 70838.3, 76078.7, 81483.5,
                        87073.8, 92825.5, 98761.2, 104855, 111119,
                        117536, 124121, 130870, 137783, 144845,
                        152052, 159418, 166942, 174607, 182410,
                        190366, 198457, 206681, 215035, 223536,
                        232163, 240915, 249790, 258785, 267918,
                        277148, 286492, 295948, 305513, 315187,
                        324966, 334849, 344810, 354871, 365028,
                        375257, 385577, 395988, 406462, 417022,
                        427640, 438339, 449091, 459919, 470795,
                        481717, 492708, 503739, 514807, 525910,
                        537073 };
   const double fY[81] = { 0, 0.11, 1.13, 4.57, 10.06,
                        13.59, 14.48, 14.45, 14.36, 14.46,
                        14.81, 15.46, 16.48, 18, 20.24,
                        23.63, 28.96, 37.68, 52.23, 74.37,
                        98.98, 117.17, 127.67, 133.03, 135.17,
                        135.12, 133.29, 129.78, 124.38, 116.56,
                        105.39, 89.66, 68.79, 45.03, 22.78,
                        4.14, -11.73, -26.54, -41.88, -59.01,
                        -78.45, -99.14, -118.71, -135.33, -148.73,
                        -159.48, -168.31, -175.79, -182.34, -188.25,
                        -193.73, -198.94, -203.98, -208.96, -213.97,
                        -219.08, -224.41, -230.09, -236.3, -243.33,
                        -251.62, -261.9, -275.46, -294.47, -321.33,
                        -353.86, -383.4, -406.16, -424.11, -439.43,
                        -453.31, -466.21, -478.22, -489.26, -499.22,
                        -508.07, -515.84, -522.63, -528.6, -533.86,
                        -538.56 };
   const double fB[81] = { 0.000489818, 0.00172421, 0.00505863, 0.00820231, 0.00654987,
                        0.00199323, 0.000130042, -9.34932e-05, -1.06619e-05, 0.000123598,
                        0.000248625, 0.000379373, 0.000526334, 0.000717782, 0.000990196,
                        0.00142577, 0.00214554, 0.00336522, 0.00527566, 0.00671648,
                        0.00578148, 0.00345375, 0.00176285, 0.000794573, 0.000200928,
                        -0.000203532, -0.000535471, -0.000852391, -0.00122029, -0.00169202,
                        -0.00233742, -0.00314635, -0.00381622, -0.00382743, -0.00323432,
                        -0.00261461, -0.00225775, -0.00217074, -0.00229844, -0.00256755,
                        -0.00279809, -0.0027571, -0.00240511, -0.00193472, -0.0015171,
                        -0.00120563, -0.000989785, -0.000839387, -0.000734543, -0.000661088,
                        -0.000612758, -0.000579529, -0.000558905, -0.0005496, -0.000549226,
                        -0.000559968, -0.000583036, -0.00062109, -0.000682012, -0.000778226,
                        -0.000927781, -0.00117336, -0.00158183, -0.00224927, -0.00301897,
                        -0.00314673, -0.00252417, -0.00190532, -0.00155572, -0.00136646,
                        -0.00125286, -0.00116136, -0.0010701, -0.000968189, -0.00086318,
                        -0.000757584, -0.000658657, -0.000575363, -0.000504766, -0.000445104,
                        -0.000399319 };
   const double fC[81] = { 6.39544e-06, 6.08474e-06, 5.15261e-06, 1.20414e-06, -3.59082e-06,
                        -1.52799e-06, -1.84508e-07, 1.02884e-08, 4.5675e-08, 3.43781e-08,
                        3.23327e-08, 3.09258e-08, 3.41468e-08, 4.37093e-08, 5.9311e-08,
                        9.3801e-08, 1.43918e-07, 2.34666e-07, 3.25815e-07, 7.42235e-08,
                        -3.21119e-07, -2.65036e-07, -1.41162e-07, -8.21274e-08, -4.90224e-08,
                        -3.70873e-08, -3.08559e-08, -3.17327e-08, -3.84726e-08, -4.8806e-08,
                        -6.66445e-08, -7.39957e-08, -3.88605e-08, 3.70208e-08, 5.76751e-08,
                        3.88841e-08, 1.53159e-08, -2.42385e-09, -1.60473e-08, -2.20643e-08,
                        -9.92453e-09, 1.54891e-08, 3.12918e-08, 3.00797e-08, 2.34401e-08,
                        1.57063e-08, 1.09711e-08, 7.31754e-09, 5.23266e-09, 3.40834e-09,
                        2.19343e-09, 1.60327e-09, 7.20749e-10, 3.13638e-10, -2.72701e-10,
                        -8.91159e-10, -1.57762e-09, -2.44678e-09, -3.92213e-09, -6.02405e-09,
                        -9.26905e-09, -1.55793e-08, -2.54263e-08, -4.09162e-08, -3.48613e-08,
                        2.23702e-08, 3.79529e-08, 2.14896e-08, 1.18885e-08, 6.03355e-09,
                        4.66589e-09, 3.88583e-09, 4.60176e-09, 4.80995e-09, 4.84484e-09,
                        4.82395e-09, 4.17665e-09, 3.37424e-09, 3.00407e-09, 2.36965e-09,
                        11163.2 };
   const double fD[81] = { -1.04712e-09, -1.04712e-09, -2.66136e-09, -2.30851e-09, 7.72443e-10,
                        4.11607e-10, 5.06072e-11, 7.96945e-12, -2.2453e-12, -3.63783e-13,
                        -2.2689e-13, 4.75397e-13, 1.29626e-12, 1.96673e-12, 4.04124e-12,
                        5.51751e-12, 9.38923e-12, 8.91377e-12, -2.32846e-11, -3.4798e-11,
                        4.70751e-12, 9.91916e-12, 4.53789e-12, 2.43789e-12, 8.46999e-13,
                        4.25157e-13, -5.77162e-14, -4.28716e-13, -6.3729e-13, -1.06366e-12,
                        -4.26033e-13, 1.97312e-12, 4.15058e-12, 1.09922e-12, -9.75967e-13,
                        -1.19317e-12, -8.76124e-13, -6.56866e-13, -2.84039e-13, 5.61487e-13,
                        1.15002e-12, 7.00076e-13, -5.27119e-14, -2.83632e-13, -3.24005e-13,
                        -1.9508e-13, -1.48092e-13, -8.31896e-14, -7.15354e-14, -4.69394e-14,
                        -2.24768e-14, -3.31476e-14, -1.50866e-14, -2.13985e-14, -2.23357e-14,
                        -2.44891e-14, -3.06394e-14, -5.14116e-14, -7.24291e-14, -1.10609e-13,
                        -2.12831e-13, -3.29504e-13, -5.13229e-13, 1.98702e-13, 1.86515e-12,
                        5.03293e-13, -5.2712e-13, -3.05551e-13, -1.84819e-13, -4.29361e-14,
                        -2.43033e-14, 2.21957e-14, 6.40848e-15, 1.06931e-15, -6.37667e-16,
                        -1.96308e-14, -2.42472e-14, -1.11479e-14, -1.90473e-14, -1.90473e-14,
                        5152.26 };
   int klow=0;
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((x-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=x-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
}

double G4HadronElastic::betaSpline200(double x) {
   const int fNp = 81, fKstep = 0;
   const double fDelta = -1, fXmin = 0, fXmax = 567246;
   const double fX[81] = { 0, 104.48, 417.919, 940.318, 1671.68,
                        2612, 3761.27, 5116.69, 6680.26, 8451.98,
                        10427.8, 12611, 15001.6, 17594.2, 20387.8,
                        23387.1, 26592.2, 29996.1, 33597.8, 37396,
                        41397.5, 45594, 49984.2, 54566.9, 59341,
                        64315.3, 69469, 74821.2, 80349.7, 86063.6,
                        91961.8, 98043.2, 104306, 110737, 117361,
                        124135, 131099, 138223, 145521, 152976,
                        160602, 168380, 176309, 184404, 192646,
                        201050, 209598, 218286, 227114, 236097,
                        245196, 254447, 263828, 273318, 282954,
                        292715, 302576, 312578, 322675, 332888,
                        343213, 353649, 364170, 374797, 385528,
                        396336, 407242, 418221, 429293, 440431,
                        451660, 462948, 474322, 485752, 497234,
                        508765, 520372, 532024, 543717, 555449,
                        567246 };
   const double fY[81] = { 0, 0.12, 1.28, 5.07, 10.53,
                        13.38, 13.81, 13.59, 13.44, 13.53,
                        13.92, 14.64, 15.78, 17.52, 20.16,
                        24.31, 31.11, 42.76, 62.49, 89.31,
                        112.48, 126.39, 133.61, 136.86, 137.57,
                        136.38, 133.48, 128.69, 121.48, 110.78,
                        94.96, 72.88, 47.09, 23.59, 4.87,
                        -10.57, -25.01, -40.56, -59.18, -81.8,
                        -106.24, -127.93, -144.7, -157.21, -166.81,
                        -174.5, -180.95, -186.59, -191.69, -196.45,
                        -201.01, -205.5, -210, -214.62, -219.45,
                        -224.62, -230.31, -236.79, -244.53, -254.34,
                        -267.79, -287.87, -318.34, -355.12, -385.45,
                        -407.41, -424.8, -440.05, -454.14, -467.29,
                        -479.41, -490.35, -500.04, -508.5, -515.84,
                        -522.19, -527.72, -532.59, -536.96, -540.96,
                        -544.7 };
   const double fB[81] = { 0.0004497, 0.0018322, 0.00543282, 0.00824061, 0.00570675,
                        0.00120532, -9.61845e-05, -0.000165433, -2.60185e-05, 0.000124735,
                        0.000263439, 0.000396798, 0.000562956, 0.000787607, 0.00112541,
                        0.00168183, 0.00263999, 0.00435355, 0.00658455, 0.00686664,
                        0.00454983, 0.00231735, 0.00109364, 0.000388186, -6.19234e-05,
                        -0.000404404, -0.00072076, -0.00107945, -0.00155234, -0.0022324,
                        -0.00316954, -0.00401199, -0.0040188, -0.00323656, -0.00249091,
                        -0.00212293, -0.00207881, -0.00232871, -0.00280346, -0.00320805,
                        -0.00307953, -0.00245947, -0.00179996, -0.00132784, -0.00102272,
                        -0.00082311, -0.000694993, -0.000608951, -0.000550376, -0.000512832,
                        -0.000491597, -0.00048047, -0.00048162, -0.000492514, -0.000512464,
                        -0.0005502, -0.000606628, -0.000697481, -0.000846032, -0.00109868,
                        -0.00154705, -0.00239015, -0.00334842, -0.00328357, -0.00238435,
                        -0.00175716, -0.00146781, -0.00132287, -0.00122615, -0.00113215,
                        -0.00102572, -0.00091084, -0.000794514, -0.000687832, -0.000592776,
                        -0.000511029, -0.000444656, -0.000393707, -0.00035569, -0.000327606,
                        -0.000307827 };
   const double fC[81] = { 6.83422e-06, 6.39802e-06, 5.08943e-06, 2.85362e-07, -3.74995e-06,
                        -1.03718e-06, -9.52725e-08, 4.41821e-08, 4.49824e-08, 4.01061e-08,
                        3.00937e-08, 3.09907e-08, 3.85158e-08, 4.81331e-08, 7.27873e-08,
                        1.12728e-07, 1.86225e-07, 3.17175e-07, 3.02257e-07, -2.27987e-07,
                        -3.50991e-07, -1.80998e-07, -9.77391e-08, -5.61974e-08, -3.80843e-08,
                        -3.0766e-08, -3.06189e-08, -3.63982e-08, -4.91382e-08, -6.98794e-08,
                        -8.90069e-08, -4.95239e-08, 4.84367e-08, 7.32015e-08, 3.93788e-08,
                        1.49443e-08, -8.60973e-09, -2.64653e-08, -3.85857e-08, -1.56859e-08,
                        3.25408e-08, 4.71736e-08, 3.60007e-08, 2.23247e-08, 1.4696e-08,
                        9.05527e-09, 5.93348e-09, 3.96944e-09, 2.66614e-09, 1.51312e-09,
                        8.20759e-10, 3.82064e-10, -5.04693e-10, -6.43315e-10, -1.42692e-09,
                        -2.43926e-09, -3.28314e-09, -5.80025e-09, -8.91124e-09, -1.58278e-08,
                        -2.75959e-08, -5.31903e-08, -3.78918e-08, 4.39947e-08, 3.98032e-08,
                        1.82293e-08, 8.30006e-09, 4.90317e-09, 3.83193e-09, 4.60693e-09,
                        4.87195e-09, 5.30461e-09, 4.92288e-09, 4.41108e-09, 3.8677e-09,
                        3.22122e-09, 2.49716e-09, 1.8756e-09, 1.3757e-09, 1.01806e-09,
                        11797.8 };
   const double fD[81] = { -1.39165e-09, -1.39165e-09, -3.06539e-09, -1.83919e-09, 9.61651e-10,
                        2.73188e-10, 3.42957e-11, 1.70629e-13, -9.17434e-13, -1.68913e-12,
                        1.36953e-13, 1.04929e-12, 1.23647e-12, 2.94179e-12, 4.43888e-12,
                        7.64389e-12, 1.28233e-11, -1.38066e-12, -4.65348e-11, -1.02464e-11,
                        1.35029e-11, 6.32161e-12, 3.02161e-12, 1.26468e-12, 4.90413e-13,
                        9.51446e-15, -3.59929e-13, -7.68147e-13, -1.20998e-12, -1.08097e-12,
                        2.16417e-12, 5.21349e-12, 1.28363e-12, -1.70222e-12, -1.20238e-12,
                        -1.12742e-12, -8.3539e-13, -5.53583e-13, 1.02392e-12, 2.10816e-12,
                        6.27063e-13, -4.69688e-13, -5.63167e-13, -3.08538e-13, -2.23725e-13,
                        -1.21742e-13, -7.53495e-14, -4.9214e-14, -4.27843e-14, -2.53647e-14,
                        -1.58073e-14, -3.15067e-14, -4.86926e-15, -2.71059e-14, -3.45729e-14,
                        -2.85259e-14, -8.38854e-14, -1.02697e-13, -2.25754e-13, -3.79907e-13,
                        -8.1749e-13, 4.84694e-13, 2.56853e-12, -1.30204e-13, -6.65391e-13,
                        -3.0346e-13, -1.03141e-13, -3.22494e-14, 2.31929e-14, 7.86774e-15,
                        1.27754e-14, -1.11873e-14, -1.49264e-14, -1.57749e-14, -1.86872e-14,
                        -2.0794e-14, -1.77821e-14, -1.42507e-14, -1.01613e-14, -1.01613e-14,
                        5444.32 };
   int klow=0;
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((x-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=x-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
}

double G4HadronElastic::betaSpline210(double x) {
   const int fNp = 81, fKstep = 0;
   const double fDelta = -1, fXmin = 0, fXmax = 597573;
   const double fX[81] = { 0, 110.204, 440.815, 990.591, 1761.6,
                        2750.95, 3962.36, 5391.28, 7036.47, 8904.14,
                        10987.3, 13284.6, 15804.8, 18533, 21478.3,
                        24640.7, 28013.6, 31595.8, 35393.5, 39398.8,
                        43610.4, 48036, 52657.2, 57490.6, 62516.5,
                        67753, 73189.3, 78824.4, 84645.4, 90673.8,
                        96885, 103289, 109885, 116672, 123634,
                        130784, 138106, 145612, 153301, 161156,
                        169191, 177388, 185745, 194278, 202966,
                        211808, 220801, 229962, 239271, 248726,
                        258323, 268061, 277938, 287951, 298099,
                        308378, 318764, 329300, 339938, 350698,
                        361579, 372554, 383644, 394846, 406135,
                        417532, 429009, 440589, 452243, 463996,
                        475818, 487707, 499686, 511725, 523823,
                        535975, 548207, 560490, 572788, 585160,
                        597573 };
   const double fY[81] = { 0, 0.14, 1.44, 5.58, 10.87,
                        13.05, 13.08, 12.7, 12.49, 12.57,
                        12.97, 13.74, 15, 16.97, 20.07,
                        25.16, 33.98, 49.97, 76.5, 105.46,
                        124.36, 134.12, 138.72, 140.3, 139.82,
                        137.59, 133.52, 127.11, 117.23, 101.89,
                        79.03, 51.07, 26.07, 7.31, -7.43,
                        -21.02, -36.14, -55.85, -82.58, -112.48,
                        -136.65, -153.12, -164.44, -172.75, -179.29,
                        -184.74, -189.52, -193.88, -197.99, -201.99,
                        -205.98, -210.05, -214.31, -218.88, -223.93,
                        -229.71, -236.66, -245.56, -258.02, -277.43,
                        -309.38, -350.19, -382.31, -404.39, -421.81,
                        -437.41, -452.13, -465.98, -478.68, -489.96,
                        -499.75, -508.13, -515.29, -521.44, -526.8,
                        -531.57, -535.92, -539.98, -543.85, -547.6,
                        -551.24 };
   const double fB[81] = { 0.000552028, 0.00197549, 0.00576968, 0.00819393, 0.00480341,
                        0.000574153, -0.000267478, -0.00022145, -4.06937e-05, 0.000118961,
                        0.000261486, 0.000409906, 0.00059694, 0.000859893, 0.00127533,
                        0.00200629, 0.00334235, 0.00580791, 0.00768072, 0.00610669,
                        0.00313905, 0.00147054, 0.000604873, 9.09002e-05, -0.000266595,
                        -0.000582089, -0.000922764, -0.00137877, -0.00205775, -0.00310162,
                        -0.00419203, -0.00425544, -0.00326599, -0.00236044, -0.00194409,
                        -0.00191949, -0.00227369, -0.00305226, -0.00381218, -0.00354225,
                        -0.00247428, -0.0016231, -0.00113043, -0.000844349, -0.000674152,
                        -0.000567574, -0.000499886, -0.000455629, -0.000430116, -0.000417823,
                        -0.000415113, -0.000422747, -0.000441569, -0.000473961, -0.000524543,
                        -0.000607424, -0.000739559, -0.000974835, -0.00140759, -0.00231688,
                        -0.00353853, -0.00348892, -0.00235733, -0.00168856, -0.00143201,
                        -0.00132003, -0.00124229, -0.00114699, -0.00102722, -0.000893056,
                        -0.000764532, -0.000648201, -0.000550918, -0.000473911, -0.000415253,
                        -0.000371999, -0.00034126, -0.000321486, -0.000308474, -0.000297962,
                        -0.000288765 };
   const double fC[81] = { 6.63838e-06, 6.27829e-06, 5.19799e-06, -7.88474e-07, -3.60901e-06,
                        -6.6578e-07, -2.89718e-08, 6.11837e-08, 4.86858e-08, 3.67978e-08,
                        3.16215e-08, 3.29844e-08, 4.12295e-08, 5.51535e-08, 8.58951e-08,
                        1.45246e-07, 2.50867e-07, 4.37416e-07, 5.57342e-08, -4.48723e-07,
                        -2.55896e-07, -1.21124e-07, -6.62006e-08, -4.01365e-08, -3.09941e-08,
                        -2.92554e-08, -3.34105e-08, -4.75126e-08, -6.91305e-08, -1.04027e-07,
                        -7.15292e-08, 6.16279e-08, 8.83776e-08, 4.50514e-08, 1.47483e-08,
                        -1.1308e-08, -3.70706e-08, -6.6658e-08, -3.21743e-08, 6.65379e-08,
                        6.63746e-08, 3.74622e-08, 2.14883e-08, 1.20417e-08, 7.54773e-09,
                        4.50597e-09, 3.02055e-09, 1.81022e-09, 9.30494e-10, 3.69694e-10,
                        -8.72792e-11, -6.9666e-10, -1.20904e-09, -2.02595e-09, -2.95879e-09,
                        -5.10433e-09, -7.61748e-09, -1.47136e-08, -2.59664e-08, -5.85365e-08,
                        -5.374e-08, 5.82603e-08, 4.37782e-08, 1.5918e-08, 6.80907e-09,
                        3.01611e-09, 3.75765e-09, 4.472e-09, 5.80519e-09, 5.60967e-09,
                        5.26186e-09, 4.52355e-09, 3.59757e-09, 2.79848e-09, 2.05041e-09,
                        1.50898e-09, 1.00389e-09, 6.06045e-10, 4.51986e-10, 3.97684e-10,
                        12412.6 };
   const double fD[81] = { -1.08919e-09, -1.08919e-09, -3.62964e-09, -1.21941e-09, 9.9164e-10,
                        1.75225e-10, 2.10311e-11, -2.53222e-12, -2.12172e-12, -8.28291e-13,
                        1.97746e-13, 1.09054e-12, 1.70124e-12, 3.47915e-12, 6.25585e-12,
                        1.04381e-11, 1.73589e-11, -3.35015e-11, -4.19827e-11, 1.52613e-11,
                        1.01511e-11, 3.96171e-12, 1.79748e-12, 6.06353e-13, 1.10685e-13,
                        -2.54775e-13, -8.34189e-13, -1.23791e-12, -1.92956e-12, 1.74404e-12,
                        6.93068e-12, 1.35179e-12, -2.12797e-12, -1.45081e-12, -1.21477e-12,
                        -1.17292e-12, -1.31398e-12, 1.49494e-12, 4.18888e-12, -6.77227e-15,
                        -1.17568e-12, -6.37129e-13, -3.69059e-13, -1.72416e-13, -1.14672e-13,
                        -5.50561e-14, -4.40372e-14, -3.15009e-14, -1.97722e-14, -1.58714e-14,
                        -2.08587e-14, -1.72923e-14, -2.71949e-14, -3.06436e-14, -6.95764e-14,
                        -8.06548e-14, -2.24507e-13, -3.52594e-13, -1.00894e-12, 1.46941e-13,
                        3.4018e-12, -4.353e-13, -8.28963e-13, -2.6897e-13, -1.10937e-13,
                        2.15375e-14, 2.05627e-14, 3.81304e-14, -5.54528e-15, -9.80669e-15,
                        -2.07016e-14, -2.57667e-14, -2.21239e-14, -2.06126e-14, -1.48517e-14,
                        -1.37632e-14, -1.07973e-14, -4.17543e-15, -1.46306e-15, -1.46306e-15,
                        5742.39 };
   int klow=0;
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((x-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=x-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
}

double G4HadronElastic::betaSpline220(double x) {
   const int fNp = 81, fKstep = 0;
   const double fDelta = -1, fXmin = 0, fXmax = 628095;
   const double fX[81] = { 0, 115.655, 462.622, 1042.17, 1852.18,
                        2891.39, 4163.59, 5667.11, 7398.55, 9356.63,
                        11548.6, 13966.3, 16608.6, 19479.6, 22579.4,
                        25901.6, 29445, 33215.5, 37205.4, 41413.6,
                        45838.7, 50488.4, 55344.1, 60422.7, 65714,
                        71216.9, 76930.1, 82841, 88970.6, 95306.8,
                        101836, 108568, 115502, 122624, 129944,
                        137461, 145160, 153053, 161140, 169386,
                        177837, 186443, 195235, 204194, 213336,
                        222622, 232087, 241710, 251489, 261421,
                        271526, 281759, 292138, 302662, 313328,
                        324134, 335055, 346109, 357297, 368614,
                        380059, 391604, 403248, 415011, 426892,
                        438861, 450944, 463110, 475356, 487707,
                        500132, 512629, 525223, 537854, 550577,
                        563359, 576228, 589122, 602066, 615058,
                        628095 };
   const double fY[81] = { 0, 0.16, 1.61, 6.07, 11.1,
                        12.62, 12.28, 11.74, 11.46, 11.5,
                        11.89, 12.69, 14.04, 16.23, 19.85,
                        26.15, 37.95, 60.78, 94.65, 121.04,
                        134.51, 140.78, 143.31, 143.53, 141.98,
                        138.68, 133.19, 124.44, 110.23, 87.45,
                        57.32, 30.26, 11.23, -2.7, -14.95,
                        -28.61, -47.93, -79.15, -117.61, -144.83,
                        -160.55, -170.44, -177.43, -182.87, -187.4,
                        -191.39, -195.05, -198.54, -201.98, -205.48,
                        -209.14, -213.08, -217.44, -222.44, -228.45,
                        -236.11, -246.79, -263.65, -293.57, -337.92,
                        -374.01, -397.21, -415.1, -431.43, -447.27,
                        -462.43, -476.28, -488.39, -498.69, -507.37,
                        -514.71, -520.98, -526.4, -531.21, -535.57,
                        -539.67, -543.66, -547.63, -551.65, -555.71,
                        -559.75 };
   const double fB[81] = { 0.000617459, 0.00213262, 0.00607469, 0.00802314, 0.00393304,
                        3.73373e-05, -0.000410892, -0.000276144, -6.44759e-05, 9.97067e-05,
                        0.000251911, 0.000411902, 0.000618209, 0.000926404, 0.00145343,
                        0.00242859, 0.00447021, 0.0077575, 0.0080252, 0.00451241,
                        0.00196149, 0.000854362, 0.000246235, -0.000136298, -0.000443324,
                        -0.000762483, -0.00118213, -0.00182432, -0.00290187, -0.0042859,
                        -0.00456375, -0.00337544, -0.00224824, -0.00174269, -0.00167353,
                        -0.00203296, -0.00316718, -0.00467125, -0.0042663, -0.00245619,
                        -0.00141027, -0.000935855, -0.000682515, -0.00054242, -0.000457144,
                        -0.000405103, -0.000371981, -0.000355253, -0.000350439, -0.000355444,
                        -0.000371518, -0.000400085, -0.000443147, -0.000512441, -0.00062225,
                        -0.000814753, -0.00117802, -0.00197958, -0.0034958, -0.00380623,
                        -0.0025032, -0.00167494, -0.00143616, -0.00135527, -0.00130744,
                        -0.00121452, -0.00107352, -0.000917071, -0.000768186, -0.000642384,
                        -0.000543136, -0.000462644, -0.000402819, -0.000359654, -0.000329221,
                        -0.000313855, -0.000307777, -0.000308896, -0.000312055, -0.000312074,
                        -0.000306831 };
   const double fC[81] = { 6.76771e-06, 6.33293e-06, 5.0286e-06, -1.66661e-06, -3.38282e-06,
                        -3.65928e-07, 1.36053e-08, 7.60162e-08, 4.62338e-08, 3.7615e-08,
                        3.18233e-08, 3.43504e-08, 4.37294e-08, 6.36169e-08, 1.06403e-07,
                        1.87124e-07, 3.89055e-07, 4.82797e-07, -4.15704e-07, -4.19052e-07,
                        -1.57412e-07, -8.06975e-08, -4.45427e-08, -3.07801e-08, -2.72438e-08,
                        -3.07544e-08, -4.2698e-08, -6.59473e-08, -1.09846e-07, -1.08587e-07,
                        6.60295e-08, 1.10481e-07, 5.20747e-08, 1.8916e-08, -9.46733e-09,
                        -3.83442e-08, -1.08978e-07, -8.15728e-08, 1.31651e-07, 8.78668e-08,
                        3.58885e-08, 1.92361e-08, 9.5784e-09, 6.0591e-09, 3.26908e-09,
                        2.33519e-09, 1.16421e-09, 5.74208e-10, -8.19969e-11, -4.21834e-10,
                        -1.16902e-09, -1.62264e-09, -2.52605e-09, -4.0584e-09, -6.23677e-09,
                        -1.15781e-08, -2.16864e-08, -5.08215e-08, -8.47106e-08, 5.72802e-08,
                        5.65722e-08, 1.51639e-08, 5.34445e-09, 1.53175e-09, 2.49431e-09,
                        5.26862e-09, 6.40143e-09, 6.45805e-09, 5.69968e-09, 4.48603e-09,
                        3.50142e-09, 2.93945e-09, 1.81083e-09, 1.60649e-09, 7.85562e-10,
                        4.16557e-10, 5.5797e-11, -1.42633e-10, -1.01425e-10, 1.0001e-10,
                        13036.8 };
   const double fD[81] = { -1.25308e-09, -1.25308e-09, -3.85081e-09, -7.06251e-10, 9.67697e-10,
                        9.94422e-11, 1.38366e-11, -5.73367e-12, -1.46722e-12, -8.80754e-13,
                        3.48412e-13, 1.1832e-12, 2.30899e-12, 4.6009e-12, 8.09907e-12,
                        1.89961e-11, 8.28738e-12, -7.50638e-11, -2.65239e-13, 1.97088e-11,
                        5.4996e-12, 2.48196e-12, 9.03307e-13, 2.22776e-13, -2.12656e-13,
                        -6.9684e-13, -1.3111e-12, -2.38723e-12, 6.62395e-14, 8.91498e-12,
                        2.20094e-12, -2.80763e-12, -1.5521e-12, -1.29249e-12, -1.28041e-12,
                        -3.05816e-12, 1.15732e-12, 8.78943e-12, -1.76997e-12, -2.05005e-12,
                        -6.4498e-13, -3.66151e-13, -1.30941e-13, -1.01732e-13, -3.3523e-14,
                        -4.124e-14, -2.04369e-14, -2.2368e-14, -1.14049e-14, -2.46491e-14,
                        -1.47764e-14, -2.90122e-14, -4.85355e-14, -6.80783e-14, -1.64769e-13,
                        -3.08545e-13, -8.78497e-13, -1.00977e-12, 4.18218e-12, -2.0621e-14,
                        -1.19548e-12, -2.81118e-13, -1.08041e-13, 2.70062e-14, 7.72583e-14,
                        3.12524e-14, 1.5513e-15, -2.06422e-14, -3.27549e-14, -2.64138e-14,
                        -1.49895e-14, -2.98718e-14, -5.39249e-15, -2.15083e-14, -9.62269e-15,
                        -9.34427e-15, -5.13003e-15, 1.06117e-15, 5.1682e-15, 5.1682e-15,
                        6029.66 };
   int klow=0;
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((x-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=x-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
}

double G4HadronElastic::betaSpline230(double x) {
   const int fNp = 81, fKstep = 0;
   const double fDelta = -1, fXmin = 0, fXmax = 658832;
   const double fX[81] = { 0, 121.239, 485.824, 1092.45, 1941.56,
                        3033.14, 4367.2, 5943.74, 7759.28, 9816.43,
                        12110.8, 14646, 17421.9, 20432.9, 23683.8,
                        27168, 30884.3, 34838.8, 39023.6, 43437.5,
                        48079.2, 52956.5, 58050.3, 63377.9, 68929.1,
                        74702.5, 80685.6, 86899.2, 93319.2, 99955.4,
                        106807, 113872, 121149, 128624, 136308,
                        144184, 152267, 160538, 169012, 177671,
                        186529, 195567, 204783, 214175, 223759,
                        233515, 243440, 253532, 263788, 274206,
                        284783, 295540, 306430, 317451, 328643,
                        339961, 351447, 363052, 374773, 386632,
                        398625, 410751, 422981, 435313, 447769,
                        460321, 472991, 485752, 498598, 511556,
                        524594, 537710, 550899, 564159, 577517,
                        590910, 604396, 617909, 631509, 645130,
                        658832 };
   const double fY[81] = { 0, 0.18, 1.79, 6.57, 11.24,
                        12.12, 11.44, 10.75, 10.37, 10.34,
                        10.68, 11.46, 12.86, 15.27, 19.49,
                        27.49, 44.28, 78.25, 115.35, 134.68,
                        143.09, 146.6, 147.48, 146.55, 143.97,
                        139.43, 131.97, 119.46, 97.96, 66.24,
                        36.29, 16.4, 3.14, -7.55, -18.77,
                        -35.14, -68.42, -121.54, -152.85, -167.18,
                        -175.38, -181.02, -185.39, -189.04, -192.29,
                        -195.31, -198.23, -201.16, -204.2, -207.47,
                        -211.1, -215.29, -220.31, -226.67, -235.39,
                        -248.8, -272.86, -315.86, -358.89, -385.38,
                        -404.32, -421.53, -438.91, -456.22, -472.24,
                        -486.02, -497.4, -506.75, -514.53, -521.14,
                        -526.86, -531.93, -536.49, -540.73, -544.79,
                        -548.82, -552.94, -557.23, -561.69, -566.21,
                        -570.62 };
   const double fB[81] = { 0.000679838, 0.00227121, 0.00639526, 0.00783035, 0.00305462,
                        -0.00039951, -0.000531049, -0.000330741, -0.000105166, 6.81867e-05,
                        0.000224582, 0.000393424, 0.00062937, 0.000991531, 0.00168321,
                        0.00302483, 0.00660307, 0.00971582, 0.00688511, 0.00264809,
                        0.00114549, 0.000389706, -1.66335e-05, -0.000321146, -0.00061109,
                        -0.000981925, -0.00155298, -0.00257327, -0.00421883, -0.0049052,
                        -0.00362423, -0.00218804, -0.00154897, -0.0013776, -0.0016114,
                        -0.00278673, -0.00580419, -0.00556782, -0.00233266, -0.00116981,
                        -0.00073712, -0.000533677, -0.000423988, -0.00035948, -0.000321752,
                        -0.000299782, -0.000290575, -0.000291598, -0.000303237, -0.000326177,
                        -0.000362998, -0.000419574, -0.000509777, -0.000654834, -0.000938058,
                        -0.00148267, -0.00296504, -0.00404544, -0.00298284, -0.00174988,
                        -0.00145939, -0.00140885, -0.00142603, -0.00136124, -0.00119903,
                        -0.000995335, -0.000808622, -0.000663376, -0.000553366, -0.000470763,
                        -0.0004103, -0.000364021, -0.000330632, -0.000309958, -0.000300659,
                        -0.000301925, -0.000310901, -0.000323421, -0.000331613, -0.000329466,
                        -0.00031159 };
   const double fC[81] = { 6.78936e-06, 6.33659e-06, 4.97502e-06, -2.60934e-06, -3.01509e-06,
                        -1.49241e-07, 5.06402e-08, 7.64156e-08, 4.78312e-08, 3.64373e-08,
                        3.17264e-08, 3.48736e-08, 5.01246e-08, 7.01545e-08, 1.42612e-07,
                        2.4244e-07, 7.20405e-07, 6.67411e-08, -7.43165e-07, -2.16761e-07,
                        -1.0696e-07, -4.79998e-08, -3.17715e-08, -2.53853e-08, -2.68458e-08,
                        -3.73859e-08, -5.80581e-08, -1.06145e-07, -1.50174e-07, 4.67461e-08,
                        1.40222e-07, 6.3059e-08, 2.47545e-08, -1.82742e-09, -2.86011e-08,
                        -1.20613e-07, -2.52716e-07, 2.81293e-07, 1.00488e-07, 3.38133e-08,
                        1.50348e-08, 7.47403e-09, 4.42746e-09, 2.44083e-09, 1.49575e-09,
                        7.56252e-10, 1.71444e-10, -2.72809e-10, -8.62102e-10, -1.33977e-09,
                        -2.14128e-09, -3.1184e-09, -5.16447e-09, -7.99791e-09, -1.73066e-08,
                        -3.08147e-08, -9.82469e-08, 5.15071e-09, 8.55058e-08, 1.84674e-08,
                        5.75401e-09, -1.58573e-09, 1.80981e-10, 5.07262e-09, 7.94975e-09,
                        8.27898e-09, 6.45678e-09, 4.92592e-09, 3.63744e-09, 2.73742e-09,
                        1.89991e-09, 1.62874e-09, 9.0275e-10, 6.56292e-10, 3.98431e-11,
                        -1.34362e-10, -5.31234e-10, -3.95307e-10, -2.0705e-10, 3.64697e-10,
                        13702.4 };
   const double fD[81] = { -1.24485e-09, -1.24485e-09, -4.1675e-09, -1.59285e-10, 8.75135e-10,
                        4.9943e-11, 5.44979e-12, -5.24811e-12, -1.84622e-12, -6.84407e-13,
                        4.13811e-13, 1.83136e-12, 2.21742e-12, 7.42958e-12, 9.55035e-12,
                        4.28708e-11, -5.50991e-11, -6.45112e-11, 3.97535e-11, 7.88512e-12,
                        4.02963e-12, 1.06196e-12, 3.99562e-13, -8.77011e-14, -6.08541e-13,
                        -1.1517e-12, -2.57965e-12, -2.28607e-12, 9.89116e-12, 4.54788e-12,
                        -3.64063e-12, -1.75446e-12, -1.18543e-12, -1.16147e-12, -3.89384e-12,
                        -5.44806e-12, 2.15204e-11, -7.11226e-12, -2.5668e-12, -7.06667e-13,
                        -2.78843e-13, -1.10186e-13, -7.05072e-14, -3.287e-14, -2.52673e-14,
                        -1.96412e-14, -1.4674e-14, -1.91527e-14, -1.52833e-14, -2.52579e-14,
                        -3.02802e-14, -6.26267e-14, -8.57012e-14, -2.77225e-13, -3.97853e-13,
                        -1.95698e-12, 2.96986e-12, 2.28517e-12, -1.88442e-12, -3.53349e-13,
                        -2.01766e-13, 4.81504e-14, 1.32223e-13, 7.69934e-14, 8.74335e-15,
                        -4.79371e-14, -3.99903e-14, -3.34329e-14, -2.31525e-14, -2.14118e-14,
                        -6.89191e-15, -1.83477e-14, -6.19536e-15, -1.53827e-14, -4.33571e-15,
                        -9.81009e-15, 3.35294e-15, 4.61418e-15, 1.39917e-14, 1.39917e-14,
                        6320.52 };
   int klow=0;
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((x-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=x-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
}

double G4HadronElastic::betaSpline240(double x) {
   const int fNp = 81, fKstep = 0;
   const double fDelta = -1, fXmin = 0, fXmax = 689714;
   const double fX[81] = { 0, 126.954, 508.704, 1143.92, 2033.04,
                        3176.06, 4573, 6220.73, 8121.48, 10275.2,
                        12682, 15336.9, 18238.7, 21391.6, 24789.6,
                        28438, 32329.6, 36470.7, 40853.3, 45476.1,
                        50337.8, 55437, 60772.4, 66352.9, 72157.4,
                        78205.1, 84473.3, 90971.2, 97697.5, 104651,
                        111817, 119220, 126833, 134652, 142690,
                        150946, 159402, 168073, 176940, 186001,
                        195270, 204730, 214376, 224226, 234259,
                        244454, 254845, 265412, 276152, 287063,
                        298142, 309387, 320795, 332342, 344046,
                        355906, 367919, 380059, 392346, 404778,
                        417328, 429991, 442792, 455727, 468767,
                        481909, 495176, 508512, 521996, 535541,
                        549172, 562915, 576738, 590607, 604580,
                        618623, 632732, 646906, 661109, 675400,
                        689714 };
   const double fY[81] = { 0, 0.2, 1.99, 7.05, 11.31,
                        11.54, 10.52, 9.64, 9.12, 8.97,
                        9.2, 9.89, 11.27, 13.82, 18.73,
                        29.32, 55.9, 104.8, 134.61, 145.9,
                        150.34, 151.72, 151.26, 149.28, 145.59,
                        139.44, 128.96, 110.08, 78.47, 44.75,
                        22.83, 9.65, 0.29, -8.28, -19.58,
                        -45.99, -123.32, -161.19, -173.21, -179.48,
                        -183.72, -187.03, -189.85, -192.41, -194.85,
                        -197.26, -199.75, -202.41, -205.36, -208.75,
                        -212.81, -217.91, -224.75, -234.89, -252.22,
                        -285.98, -334.18, -367.64, -388.86, -406.97,
                        -425.77, -445.86, -465.4, -482.23, -495.68,
                        -506.25, -514.72, -521.73, -527.71, -532.96,
                        -537.72, -542.16, -546.46, -550.74, -555.13,
                        -559.66, -564.32, -569.02, -573.63, -578.01,
                        -582.07 };
   const double fB[81] = { 0.000690906, 0.00243304, 0.00670235, 0.00754601, 0.00222554,
                        -0.000792059, -0.000652547, -0.000407231, -0.000162226, 1.46148e-05,
                        0.000172599, 0.000353827, 0.00060917, 0.00104993, 0.00193033,
                        0.00422927, 0.0102387, 0.0104612, 0.00399636, 0.00139278,
                        0.000527243, 6.42468e-05, -0.000223752, -0.000486577, -0.000796205,
                        -0.00127928, -0.00214213, -0.00386463, -0.00517387, -0.00408275,
                        -0.00225704, -0.0014258, -0.00107586, -0.00125054, -0.00141942,
                        -0.00688965, -0.00791657, -0.00210001, -0.000875269, -0.000546311,
                        -0.000390123, -0.00031608, -0.000272832, -0.000249596, -0.000238209,
                        -0.000236447, -0.000244074, -0.000261222, -0.000290058, -0.000334392,
                        -0.000403436, -0.000511221, -0.000710025, -0.00108028, -0.00204203,
                        -0.00372352, -0.00362819, -0.00208064, -0.00150505, -0.00145137,
                        -0.0015533, -0.00158824, -0.00143288, -0.00116448, -0.00090806,
                        -0.000711399, -0.000575191, -0.000480247, -0.000411721, -0.000366219,
                        -0.000333924, -0.000315045, -0.000308389, -0.000310433, -0.000318205,
                        -0.000327037, -0.000332173, -0.00032986, -0.000316968, -0.000295521,
                        -0.000271324 };
   const double fC[81] = { 7.17811e-06, 6.54446e-06, 4.63906e-06, -3.31091e-06, -2.67305e-06,
                        3.30402e-08, 6.68305e-08, 8.20507e-08, 4.68485e-08, 3.52592e-08,
                        3.03818e-08, 3.78798e-08, 5.01172e-08, 8.96755e-08, 1.69418e-07,
                        4.60714e-07, 1.08348e-06, -1.02976e-06, -4.45348e-07, -1.17858e-07,
                        -6.01743e-08, -3.06234e-08, -2.3355e-08, -2.37418e-08, -2.96017e-08,
                        -5.02744e-08, -8.7382e-08, -1.77704e-07, -1.69412e-08, 1.73861e-07,
                        8.09126e-08, 3.13673e-08, 1.4601e-08, -3.69409e-08, 1.59307e-08,
                        -6.78528e-07, 5.5709e-07, 1.13742e-07, 2.43803e-08, 1.19248e-08,
                        4.92439e-09, 2.90294e-09, 1.58032e-09, 7.78685e-10, 3.5621e-10,
                        -1.83312e-10, -5.50675e-10, -1.07216e-09, -1.61276e-09, -2.45056e-09,
                        -3.78133e-09, -5.80386e-09, -1.16223e-08, -2.04438e-08, -6.17257e-08,
                        -8.00522e-08, 8.79873e-08, 3.94937e-08, 7.35164e-09, -3.03409e-09,
                        -5.08832e-09, 2.32981e-09, 9.80678e-09, 1.09435e-08, 8.72059e-09,
                        6.2439e-09, 4.0221e-09, 3.09753e-09, 1.98444e-09, 1.37481e-09,
                        9.94492e-10, 3.79228e-10, 1.02273e-10, -2.49624e-10, -3.06624e-10,
                        -3.22319e-10, -4.16618e-11, 2.04829e-10, 7.02894e-10, 7.97768e-10,
                        14313.5 };
   const double fD[81] = { -1.66374e-09, -1.66374e-09, -4.17182e-09, 2.39136e-10, 7.89159e-10,
                        8.06295e-12, 3.07903e-12, -6.1734e-12, -1.79365e-12, -6.75505e-13,
                        9.414e-13, 1.40576e-12, 4.18213e-12, 7.82246e-12, 2.66144e-11,
                        5.33423e-11, -1.70102e-10, 4.4449e-11, 2.36141e-11, 3.95503e-12,
                        1.93173e-12, 4.54095e-13, -2.31032e-14, -3.36519e-13, -1.13941e-12,
                        -1.97333e-12, -4.63343e-12, 7.96691e-12, 9.14669e-12, -4.32355e-12,
                        -2.23077e-12, -7.34157e-13, -2.19727e-12, 2.19256e-12, -2.80395e-11,
                        4.87056e-11, -1.7044e-11, -3.35932e-12, -4.58212e-13, -2.51729e-13,
                        -7.12324e-14, -4.57023e-14, -2.71279e-14, -1.40361e-14, -1.76417e-14,
                        -1.17844e-14, -1.64503e-14, -1.67784e-14, -2.55954e-14, -4.00387e-14,
                        -5.99539e-14, -1.70004e-13, -2.54662e-13, -1.17567e-12, -5.1508e-13,
                        4.66269e-12, -1.33157e-12, -8.71976e-13, -2.78461e-13, -5.45624e-14,
                        1.9526e-13, 1.94706e-13, 2.92943e-14, -5.68237e-14, -6.28196e-14,
                        -5.58189e-14, -2.31104e-14, -2.75161e-14, -1.50025e-14, -9.30037e-15,
                        -1.49232e-14, -6.67866e-15, -8.45752e-15, -1.35981e-15, -3.72555e-16,
                        6.63029e-15, 5.79688e-15, 1.16897e-14, 2.2128e-15, 2.2128e-15,
                        6634.04 };
   int klow=0;
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((x-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=x-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
}

double G4HadronElastic::CSSpline(double x) {
   const int fNp = 45, fKstep = 0;
   const double fDelta = -1, fXmin = 4685.33, fXmax = 1.31943e+06;
   const double fX[45] = { 4685.33, 8260.23, 12834.3, 18400.9, 24951.6,
                        32477, 40965.8, 50405.6, 60782.5, 71547.8,
                        83710.8, 96761.7, 110681, 125449, 140317,
                        156679, 173025, 190889, 208622, 227886,
                        246906, 267463, 287662, 308396, 329638,
                        351358, 374594, 397203, 442434, 490059,
                        537568, 585837, 634618, 683663, 731526,
                        776794, 825182, 869399, 926321, 978168,
                        1.03012e+06, 1.07858e+06, 1.12337e+06, 1.23826e+06, 1.31943e+06 };
   const double fY[45] = { 316, 261, 226, 149, 76.6,
                        46.4, 24.9, 12.4, 6.61, 3.51,
                        1.86, 1.08, 0.648, 0.468, 0.394,
                        0.362, 0.322, 0.278, 0.222, 0.169,
                        0.118, 0.0766, 0.0483, 0.0287, 0.0166,
                        0.00976, 0.00645, 0.00471, 0.00328, 0.00237,
                        0.00135, 0.000659, 0.000219, 3.67e-05, 9e-06,
                        1.59e-05, 2.59e-05, 3.11e-05, 2.92e-05, 1.8e-05,
                        1.26e-05, 7.7e-06, 4.7e-06, 7.5e-06, 1.75e-05 };
   const double fB[45] = { -0.0239173, -0.00910753, -0.00988741, -0.0145102, -0.007184,
                        -0.00258593, -0.00206924, -0.000805728, -0.000383264, -0.000202824,
                        -8.60331e-05, -4.16556e-05, -2.04103e-05, -7.02272e-06, -3.01059e-06,
                        -1.92349e-06, -2.50482e-06, -2.82554e-06, -3.06196e-06, -2.69933e-06,
                        -2.44118e-06, -1.65937e-06, -1.1612e-06, -7.40747e-07, -4.25353e-07,
                        -2.14149e-07, -9.64993e-08, -5.75414e-08, -1.83769e-08, -2.11489e-08,
                        -1.87674e-08, -1.12292e-08, -6.3403e-09, -1.63688e-09, 3.76581e-11,
                        2.16771e-10, 1.71561e-10, 6.45918e-11, -1.62305e-10, -1.77827e-10,
                        -8.65273e-11, -9.38391e-11, -4.02126e-11, 8.50435e-11, 1.59396e-10 };
   const double fC[45] = { 3.01744e-06, 1.12527e-06, -1.29577e-06, 4.65301e-07, 6.53077e-07,
                        -4.20669e-08, 1.02934e-07, 3.09154e-08, 9.79669e-09, 6.9645e-09,
                        2.63769e-09, 7.62641e-10, 7.63641e-10, 1.42891e-10, 1.26962e-10,
                        -6.05201e-11, 2.49565e-11, -4.291e-11, 2.95775e-11, -1.07542e-11,
                        2.43274e-11, 1.37028e-11, 1.09604e-11, 9.31793e-12, 5.53007e-12,
                        4.19378e-12, 8.69433e-13, 8.53693e-13, 1.21831e-14, -7.0387e-14,
                        1.20514e-13, 3.56577e-14, 6.4563e-14, 3.13372e-14, 3.6492e-15,
                        3.07537e-16, -1.24184e-15, -1.17734e-15, -2.80879e-15, 2.50942e-15,
                        -7.52178e-16, 6.01295e-16, 5.96136e-16, 4.94052e-16, 81172 };
   const double fD[45] = { -1.76431e-10, -1.76431e-10, 1.05456e-10, 9.55491e-12, -3.07912e-11,
                        5.69378e-12, -2.54308e-12, -6.78391e-13, -8.76951e-14, -1.18579e-13,
                        -4.78905e-14, 2.39367e-17, -1.40113e-14, -3.57129e-16, -3.81952e-15,
                        1.74303e-15, -1.26639e-15, 1.36259e-15, -6.97853e-16, 6.14836e-16,
                        -1.72276e-16, -4.52563e-17, -2.64046e-17, -5.94412e-17, -2.05077e-17,
                        -4.76891e-17, -2.32063e-19, -6.20157e-18, -5.77916e-19, 1.33941e-18,
                        -5.86e-19, 1.97515e-19, -2.25819e-19, -1.9283e-19, -2.46066e-20,
                        -1.06732e-20, 4.86233e-22, -9.55378e-21, 3.41919e-20, -2.09253e-20,
                        9.30993e-21, -3.83963e-23, -2.9617e-22, -2.9617e-22, 55532.1 };
   int klow=0;
   if(x<=fXmin) klow=0;
   else if(x>=fXmax) klow=fNp-1;
   else {
     if(fKstep) {
       // Equidistant knots, use histogramming
       klow = int((x-fXmin)/fDelta);
       if (klow < fNp-1) klow = fNp-1;
     } else {
       int khig=fNp-1, khalf;
       // Non equidistant knots, binary search
       while(khig-klow>1)
         if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
         else khig=khalf;
     }
   }
   // Evaluate now
   double dx=x-fX[klow];
   return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
}

double G4HadronElastic::CSSpline2(double x)
{
  if ( x < 1.31943e+06) return CSSpline(x);
  else return 0.0;
}

G4double 
G4HadronElastic::SampleInvariantT2(G4double tm)
{
  static const G4double MeV2 = MeV*MeV;
  double csmax = 512.447;
  double xmax = tm/MeV2, xmin = 0.0;
  loop:
  double r1 = G4UniformRand();
  double x = xmin + r1*(xmax - xmin);
  double y = csmax * G4UniformRand();
  if ( y < CSSpline2(x) ) return x*MeV2;
  else goto loop;
}

  ////////////////////////////////////////////////////////////////////////////////
