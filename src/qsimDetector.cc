#include "qsimDetector.hh"
#include "G4SDManager.hh"
//#include "CLHEP/Random/RandFlat.h"
//#include "CLHEP/Random/RandGauss.h"

//
#include <fstream>
#include <iostream>
//

//
static std::ifstream myfile;
//

qsimDetector::qsimDetector( G4String name, G4int detnum ) : G4VSensitiveDetector(name){
    char colname[255];

    fDetNo = detnum;
    assert( fDetNo > 0 );

    fHCID = -1;

//    fTrackSecondaries = false;
    fTrackSecondaries = true;

    sprintf(colname, "genhit_%s_%d", name.data(), detnum);
    collectionName.insert(G4String(colname));

}

qsimDetector::~qsimDetector(){
}

void qsimDetector::Initialize(G4HCofThisEvent *){

    fHitColl = new qsimDetectorHitsCollection( SensitiveDetectorName, collectionName[0] );
}

///////////////////////////////////////////////////////////////////////

G4bool qsimDetector::ProcessHits( G4Step *step, G4TouchableHistory *){
    G4bool badedep = false;
    G4bool badhit  = false;


    // Get touchable volume info
    G4TouchableHistory *hist = 
	(G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
    G4int  copyID = hist->GetReplicaNumber();

    G4StepPoint *prestep = step->GetPreStepPoint();
    G4Track     *track   = step->GetTrack();

//    G4Material* material = track->GetMaterial();

//    printf("Standard detector %d hit by %s!\n", fDetNo, track->GetParticleDefinition()->GetParticleName().data());

//    G4double edep = step->GetTotalEnergyDeposit();

    //  Make pointer to new hit if it's a valid track
    qsimDetectorHit *thishit;
    if( !badhit ){
	thishit = new qsimDetectorHit(fDetNo, copyID);
	fHitColl->insert( thishit );
    }

    if( !badhit ){
	// Hit
	thishit->f3X = prestep->GetPosition();
	thishit->f3V = track->GetVertexPosition();
	thishit->f3D = track->GetVertexMomentumDirection();
	thishit->f3P = track->GetMomentum();

	thishit->fP = track->GetMomentum().mag();
	thishit->fE = track->GetTotalEnergy();
        //
        thishit->fLambda = 1.239824/track->GetTotalEnergy();
        myfile.open ("/home/bulacarl/solid/qsim/build/PMTQE.txt");
        if (myfile.is_open()) {
          while (!myfile.eof()) {
          myfile >> lambda >> QE;
           if (1.239824/track->GetTotalEnergy() >= 1000*lambda){
           //G4cout << 1.239824/track->GetTotalEnergy() << "\t" << 1000*lambda << "\t" << QE << G4endl;
           thisLambda = 1000*lambda;
           thisQE = QE;  
           }
          }
        }
        myfile.close();
        //G4cout << thisLambda << "\t" << thisQE << G4endl;
        thishit->fEff = thisQE/100;
        //
	thishit->fM = track->GetDefinition()->GetPDGMass();

	thishit->fTrID  = track->GetTrackID();
	thishit->fmTrID = track->GetParentID();
	thishit->fPID   = track->GetDefinition()->GetPDGEncoding();

	// FIXME - Enumerate encodings
	thishit->fGen   = (long int) track->GetCreatorProcess();
    }

    return !badedep && !badhit;
}

///////////////////////////////////////////////////////////////////////

void qsimDetector::EndOfEvent(G4HCofThisEvent*HCE) {
    G4SDManager *sdman = G4SDManager::GetSDMpointer();

    if(fHCID<0){ fHCID = sdman->GetCollectionID(collectionName[0]); }

    HCE->AddHitsCollection( fHCID, fHitColl );
    return;
}


