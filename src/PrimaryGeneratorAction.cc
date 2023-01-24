////////////////////////////////////////////////////////////////////////////////////////////
// Description: This file tells the program what particle to shoot, where to shoot it out //
//   from, what angle to shoot it at, and at what energy the particle has.                //
//                                                                                        //
// Most of the code added below was written to include the Doppler effect for particles   //
// emitted in-flight. The kinetic energy of the moving source is defined by the beamE     //
// variable. Most of the time beamE=0 and the particle is emitted at rest.                //
//                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////

#include "PrimaryGeneratorAction.hh"
#include "SuNBetaDecay.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh" 
#include "G4ParticleDefinition.hh"      
#include "G4UImanager.hh"

#include "globals.hh"
#include "Randomize.hh"
//#include "CLHEP/Random/Random.h"
//#include "CLHEP/Random/RandomEngine.h"
#include "TRandom.h"

#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/math/distributions.hpp>

using namespace CLHEP;

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 1;

  // Setting the particle to a gamma ray
  G4ParticleTable *pTable = G4ParticleTable::GetParticleTable();
  G4String pName;
  G4ParticleDefinition *particle = pTable->FindParticle(pName="gamma");

  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleDefinition(particle);

  G4ParticleDefinition *particle2 = pTable->FindParticle(pName="e-");
  particleGun2 = new G4ParticleGun(n_particle);
  particleGun2->SetParticleDefinition(particle2);

//  G4ParticleDefinition *particle2 = pTable->FindParticle(pName="alpha");
//  particleGun2 = new G4ParticleGun(n_particle);
//  particleGun2->SetParticleDefinition(particle2);

}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete particleGun2;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)    
{

//  std::default_random_engine generator;
//  CLHEP::HepRandomEngine & aRanecuEngine = createEngine(21390);
//  art::Service<art::RandomNumberGenerator> rng;
//  CLHEP::HepRandomEngine aRanecuEngine = rng->getEngine();
//  CLHEP::HepRandomEngine aRanecuEngine;
//  RanecuEngine aRanecuEngine;
  //variables
//  G4double x, y, z;                       // position of source
  G4double theta ;                // angle in CoM frame
  G4double costheta, phi;                 // angle in CoM frame
  G4double cosLabTheta, sinLabTheta;      // angle in Lab frame
  G4double pxLab, pyLab, pzLab;           // momentum vector in Lab frame
  G4double dopplerEffect, gmma, velocity; // doppler variables
  G4double beamE, mamu;                   // energy variables

  mamu= 931.5 *MeV;
  beamE=0.0; //*MeV/u
  G4double totalE = 0.0;

//  std::cout << nEnergies << std::endl;
  for(int i=1; i<=nEnergies+1; i++)
   {
    if ( i!=-2 ) {
   // Generate random angle to shoot particle at
    costheta=2.0*G4UniformRand()-1.0;         // generates a random value for costheta between -1 and 1
    boost::math::beta_distribution<> beta(0.85, 0.85);
    costheta = boost::math::quantile(beta,G4UniformRand())*2 - 1;
    boost::math::beta_distribution<> beta2(2, 2);
    theta = boost::math::quantile(beta2,G4UniformRand())*pi;

//    theta=2*pi*G4UniformRand()-pi;         // generates a random value for costheta between -1 and 1
    costheta=std::cos(theta);         // generates a random value for costheta between -1 and 1
    phi= twopi*G4UniformRand();               // generates a random angle between 0 and 2pi

    // Now we start necessary calcs for Doppler effects (I believe this was initially developed by Shea Mosby)
     gmma = ((beamE*MeV)/mamu) + 1.;          // E_kinetic = gamma*mass*c*c - mass*c*c
     velocity = sqrt(1. - (1./(gmma*gmma)));  // gamma = 1/sqrt(1-(v/c)^2)

   // Move to the lab frame (we assume here that the beam axis is the z axis)
     cosLabTheta = (costheta + velocity)/(1. + velocity*costheta);   // relativistic effects
     sinLabTheta=std::sqrt(1.0-cosLabTheta*cosLabTheta);             // sin^2 + cos^2 = 1
//     sinLabTheta=std::sin(theta);             // sin^2 + cos^2 = 1
     pxLab=sinLabTheta * std::cos(phi);                              // calculates x = sintheta * cosphi
     pyLab=sinLabTheta * std::sin(phi);                              // calculates y = sintheta * sinphi
     pzLab=cosLabTheta;                                              // calculates z = costheta
     dopplerEffect = 1./(gmma*(1. - velocity*cosLabTheta));          // f = 1/gamma * 1/(1-v/c*costheta) * f0

    } else {
    
   // Generate random angle to shoot particle at
//    std::cout << costheta << std::endl;
//    std::cout << phi << std::endl;

    theta=-theta;         // generates a random value for costheta between -1 and 1
    costheta=std::cos(theta);
    if(phi>pi)
    phi= phi-pi;               // generates a random angle between 0 and 2pi
    else
    phi= pi+phi;

//    std::cout << costheta << std::endl;
//    std::cout << phi << std::endl;
    // Now we start necessary calcs for Doppler effects (I believe this was initially developed by Shea Mosby)
     gmma = ((beamE*MeV)/mamu) + 1.;          // E_kinetic = gamma*mass*c*c - mass*c*c
     velocity = sqrt(1. - (1./(gmma*gmma)));  // gamma = 1/sqrt(1-(v/c)^2)

   // Move to the lab frame (we assume here that the beam axis is the z axis)
     cosLabTheta = (costheta + velocity)/(1. + velocity*costheta);   // relativistic effects
//     sinLabTheta=std::sqrt(1.0-cosLabTheta*cosLabTheta);             // sin^2 + cos^2 = 1
     sinLabTheta=std::sin(theta);             // sin^2 + cos^2 = 1
     pxLab=sinLabTheta * std::cos(phi);                              // calculates x = sintheta * cosphi
     pyLab=sinLabTheta * std::sin(phi);                              // calculates y = sintheta * sinphi
     pzLab=cosLabTheta;                                              // calculates z = costheta
//    pxLab=-pxLab;
//    pyLab=-pyLab;
//    pzLab=-pzLab;
     dopplerEffect = 1./(gmma*(1. - velocity*cosLabTheta));          // f = 1/gamma * 1/(1-v/c*costheta) * f0
   } 
   // Sets the original postion of the particle

//     std::normal_distribution<double> x_distribution(-0.819725,6.7653875);
     x=x_distribution(generator);
     while(std::abs(x) > 10.0)
        x=x_distribution(generator);
//     std::normal_distribution<double> y_distribution(1.4020125,7.5242375);
     y=y_distribution(generator);
     while(std::abs(y) > 10.0)
        y=y_distribution(generator);
//     std::gamma_distribution<double> z_distribution(4.116,24.627);
     z=(z_distribution(generator)+594.579)/1000;
     while(std::abs(z) > 0.921 || std::abs(z) < 0.593)
        z=(z_distribution(generator)+594.579)/1000;
     x=x*mm;
     y=y*mm;
     z=z*mm;
     x=0*cm;
     y=0.0*cm;
     z=0.0*cm;
//     z=0.40*cm;
     y=-0.0*cm;
     z=-1.0*cm;
     if (i<=nEnergies){ // Set the properties for the gamma rays
       // Sets the particle's postition, momentum vector, and energy with the Doppler boost
       particleGun->SetParticlePosition(G4ThreeVector(x, y, z));

       particleGun->SetParticleMomentumDirection(G4ThreeVector(pxLab, pyLab, pzLab));

       particleGun->SetParticleEnergy(dopplerEffect*inputArray[cascade][i]*keV);

       // Generate an event
       particleGun->GeneratePrimaryVertex(anEvent);
    totalE += dopplerEffect*inputArray[cascade][i]*keV;
//    G4cout<<i<<", " <<dopplerEffect*inputArray[cascade][i]*keV<<G4endl;

     }

     if (i==nEnergies+1){ // Set the properties for the electron
//       SuNAlphaDecay alpha_decay(95, 146, 5637.82);
//       SuNBetaDecay beta_decay(95, 146, 5637.82);
       SuNBetaDecay beta_decay(28, 60, 2822.81);
//       SuNBetaDecay beta_decay(55, 137, 1175.63);
//       SuNBetaDecay beta_decay(22, 57, 10640.0);
//       electronKE_mc2 = alpha_decay.GetAlphaKineticEnergy();
       electronKE_mc2 = beta_decay.GetCorrectedElectronKineticEnergy("B");
//       electronKE_keV = electronKE_mc2 * (alpha_decay.alphaRestMassEnergy);
       electronKE_keV = electronKE_mc2 * (beta_decay.electronRestMassEnergy);
       totalE += electronKE_keV*keV;
//    G4cout<<i << ", "<<electronKE_keV*keV<<G4endl;
       particleGun2->SetParticlePosition(G4ThreeVector(x, y, z));
       particleGun2->SetParticleMomentumDirection(G4ThreeVector(pxLab, pyLab, pzLab));
       particleGun2->SetParticleEnergy(electronKE_keV*keV);
       particleGun2->GeneratePrimaryVertex(anEvent);
//    G4cout<<totalE<<G4endl;
     }
              
   }

}
