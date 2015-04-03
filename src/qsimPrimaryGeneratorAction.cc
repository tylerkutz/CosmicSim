
#include "qsimPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "qsimIO.hh"
#include "qsimEvent.hh"
#include "qsimtypes.hh"
#include "globals.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

//
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <TRandom.h>
//

//
static std::ifstream myfile;
//

qsimPrimaryGeneratorAction::qsimPrimaryGeneratorAction() {
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);


  fDefaultEvent = new qsimEvent();

  fXmin =  -7.0*cm;
  fXmax =   7.0*cm;

  fYmin =  -1.75*cm;
  fYmax =   1.75*cm;

  fZmin =  -65.0*cm;
  fZmax =  -47.0*cm;

  fEmin = 220*MeV;
  fEmax = 100000*MeV;

  fThetaMin = 40.0*deg;
  fThetaMax = 50.0*deg;
  fPhiMin = -5.0*deg;
  fPhiMax = 5.0*deg;
}

qsimPrimaryGeneratorAction::~qsimPrimaryGeneratorAction() {

  delete fParticleGun;
  delete fDefaultEvent;
}


void qsimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

    /*  Generate event, set IO data */

    // Use default, static single generator
    // Update this just in case things changed	
    // from the command user interface
    fDefaultEvent->Reset();

    // Set data //////////////////////////////////
    // Magic happens here

    TFile *f1 = new TFile("../build/test.root");
    char histname1[256] = "keMuon";
    TH1F *keMuon_distribution = (TH1F*)f1->Get(histname1);
    keMuon_distribution->GetXaxis()->SetRange(log10(fEmin),log10(fEmax));

    double xPos = CLHEP::RandFlat::shoot( fXmin, fXmax );
    double yPos = CLHEP::RandFlat::shoot( fYmin, fYmax );
    double zPos = CLHEP::RandFlat::shoot( fZmin, fZmax );
    double theta = CLHEP::RandFlat::shoot( fThetaMin, fThetaMax );
    //double xPos = -55.5*cm*sin(theta);
    //double zPos = -55.5*cm*cos(theta);
    double phi = CLHEP::RandFlat::shoot( fPhiMin, fPhiMax );
    double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
    //double E = CLHEP::RandFlat::shoot( fEmin, fEmax );
    double E = pow(keMuon_distribution->GetRandom(),10);
    E = E + 105.66;
   
    f1->Close();
    delete f1;
    
    assert( E > 0.0 );
    assert( E > mass );

    double p = sqrt( E*E - mass*mass );

    double pX = sin(theta)*cos(phi)*p;
    double pY = sin(theta)*sin(phi)*p;
    double pZ = cos(theta)*p;

    fDefaultEvent->ProduceNewParticle(
	    G4ThreeVector(xPos, yPos, zPos),
	    G4ThreeVector(pX, pY, pZ ),
	    fParticleGun->GetParticleDefinition()->GetParticleName() );

    /////////////////////////////////////////////////////////////
    // Register and create event
    // 
    double kinE = sqrt(fDefaultEvent->fPartMom[0].mag()*fDefaultEvent->fPartMom[0].mag()
	    + fDefaultEvent->fPartType[0]->GetPDGMass()*fDefaultEvent->fPartType[0]->GetPDGMass() )
	-  fDefaultEvent->fPartType[0]->GetPDGMass();

      fParticleGun->SetParticleDefinition(fDefaultEvent->fPartType[0]);
      fParticleGun->SetParticleMomentumDirection(fDefaultEvent->fPartMom[0].unit());
      fParticleGun->SetParticleEnergy( kinE  );
      fParticleGun->SetParticlePosition( fDefaultEvent->fPartPos[0] );
      fIO->SetEventData(fDefaultEvent);
      fParticleGun->GeneratePrimaryVertex(anEvent);
}

G4ParticleGun* qsimPrimaryGeneratorAction::GetParticleGun() {
  return fParticleGun;
} 

