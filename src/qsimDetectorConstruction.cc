
#include "qsimDetectorConstruction.hh"

#include "qsimDetector.hh"
#include "qsimScintDetector.hh"
#include "G4SDManager.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4GenericTrap.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpBoundaryProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

qsimDetectorConstruction::qsimDetectorConstruction()
{
  det_x = det_y = det_z = 275*cm;
  quartz_x = 1.75*cm; 
  quartz_y = 7.*cm; 
  //Change quartz thickness here. 
  quartz_z = 0.5*cm;

  quartz_zPos = -.0*cm;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

qsimDetectorConstruction::~qsimDetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* qsimDetectorConstruction::Construct()
{

//	------------- Materials -------------

  G4double a, z, density;
  G4int nelements;

// Air
// 
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

// Quartz
// 
  G4Element* Si = new G4Element("Silicon", "Si", z=14 , a=28*g/mole);

  G4Material* Quartz = new G4Material("Quartz", density= 2.203*g/cm3, nelements=2);
  Quartz->AddElement(Si, 1);
  Quartz->AddElement(O, 2);

// Mirror
// 
  G4Element* Al = new G4Element("Aluminum", "Al", z=13 , a=27*g/mole);

  G4Element* Pb = new G4Element("Lead", "Pb", z=82 , a=207.2*g/mole);

  G4Material* Alu_Mat = new G4Material("Alu_Mat", 2.7*g/cm3, nelements=1);
  Alu_Mat->AddElement(Al, 1);

  G4Material* Pb_Mat = new G4Material("Pb_Mat", 11.34*g/cm3, nelements=1);
  Pb_Mat->AddElement(Pb, 1);

	//G4Material* Pb_Mat=Air; // To remove lead bricks, uncomment.
	
  G4Material* Mirror = new G4Material("Mirror", density= 2.7*g/cm3, nelements=1);
  Mirror->AddElement(Al, 1);


// Let us make cathode from a special metal (reflectivity 0, efficiency of photoelectrons 25%)
  G4Material* CATH = new G4Material("CATH", density= 2.7*g/cm3, nelements=1);
  CATH->AddElement(Al, 1);


//
// ------------ Generate & Add Material Properties Table ------------
//



const G4int nEntries = 190;

	G4double PhotonEnergy[nEntries] =
		{  2.4,2.42,2.44,2.46,2.48,2.5,2.52,2.54,2.56,2.58,
		2.6,2.62,2.64,2.66,2.68,2.7,2.72,2.74,2.76,2.78,
		2.8,2.82,2.84,2.86,2.88,2.9,2.92,2.94,2.96,2.98,
		3,3.02,3.04,3.06,3.08,3.1,3.12,3.14,3.16,3.18,
		3.2,3.22,3.24,3.26,3.28,3.3,3.32,3.34,3.36,3.38,
		3.4,3.42,3.44,3.46,3.48,3.5,3.52,3.54,3.56,3.58,
		3.6,3.62,3.64,3.66,3.68,3.7,3.72,3.74,3.76,3.78,
		3.8,3.82,3.84,3.86,3.88,3.9,3.92,3.94,3.96,3.98,
		4,4.02,4.04,4.06,4.08,4.1,4.12 ,4.14,4.16,4.18,  //Glass cuts off above 4.135eV, 87 entries
		4.2,4.22,4.24,4.26,4.28,4.3,4.32,4.34,4.36,4.38,
		4.4,4.42,4.44,4.46,4.48,4.5,4.52,4.54,4.56,4.58,
		4.6,4.62,4.64,4.66,4.68,4.7,4.72,4.74,4.76,4.78,
		4.8,4.82,4.84,4.86,4.88,4.9,4.92,4.94,4.96,4.98, //  Cut off -> 4.96eV ~ 250nm
		5,5.02,5.04   ,   5.06,5.08,5.1,5.12,5.14,5.16,5.18,   // 5.04eV = 246 nm is the 30% cutoff, 133 entries
		5.2,5.22,5.24,5.26,5.28,5.3,5.32,5.34,5.36,5.38,
		5.4,5.42,5.44,5.46,5.48,5.5,5.52,5.54,5.56,5.58,	
		5.6,5.62,5.64,5.66,5.68,5.7,5.72,5.74,5.76,5.78,
		5.8,5.82,5.84,5.86,5.88,5.9,5.92,5.94,5.96,5.98,
		6,6.02,6.04,6.06,6.08,6.1,6.12,6.14,6.16,6.18   };  // 200 nm

	G4double RefractiveIndex1[nEntries];
	G4double Absorption1[nEntries];
	G4double RefractiveIndex2[nEntries];
	G4double RefractiveIndex3[nEntries];
	G4double Reflectivity4[nEntries];
	G4double Efficiency4[nEntries];
	G4double Reflectivity3[nEntries];

 	for (int i = 0; i < nEntries; i++) {
		PhotonEnergy[i] = PhotonEnergy[i]*eV;
		RefractiveIndex1[i]= 1.455 -(.005836*PhotonEnergy[i])+(.003374*PhotonEnergy[i]*PhotonEnergy[i]);

//Aluminum
//		Reflectivity3[i] = 0; //.6;
	
//Aluminum Real
   if (PhotonEnergy[i] < 4.135) Reflectivity3[i] = .75;  // regularly .75, .7 below  .56/.53/.46 tunes to 50 PEs
		else if (PhotonEnergy[i] >= 4.135 && PhotonEnergy[i] < 6.203) Reflectivity3[i] = .7;
   else Reflectivity3[i] = .6;		// .6
		
//ALZAK		
//		if (PhotonEnergy[i] < 3.26*eV) {
//			Reflectivity3[i]=.93; }
//		else { Reflectivity3[i] = 0;}

// No Mirror
//		Reflectivity3[i] = 0;
		
//		Absorption1[i] = 50.*cm;  //Uniform
   
        Absorption1[i] = (exp(4.325)*exp(1.191*PhotonEnergy[i])*exp(-.213*PhotonEnergy[i]*PhotonEnergy[i])*exp(-.04086*PhotonEnergy[i]*PhotonEnergy[i]*PhotonEnergy[i]))*m;

       if (Absorption1[i] > 25*m) {Absorption1[i] = 25*m;}

		RefractiveIndex2[i]=1;
		RefractiveIndex3[i]=0;
		Reflectivity4[i]=0;
		Efficiency4[i]=.25;

		
	}
	
//QUARTZ
	
  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();
  myMPT1->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1,nEntries);
  myMPT1->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption1,     nEntries);
  
  Quartz->SetMaterialPropertiesTable(myMPT1);

//
// Air
//

/*  G4double RefractiveIndex2[nEntries] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };  */

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex2, nEntries);
  
  Air->SetMaterialPropertiesTable(myMPT2);

//
// Mirror (refractive index = 0) 
//


/*  G4double RefractiveIndex3[nEntries] =
            { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00 };   */


 

  G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable();
  myMPT3->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex3, nEntries);
  // myMPT3->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption3, nEntries);  
  //  myMPT3->AddProperty("REFLECTIVITY", PhotonEnergy, Reflectivity3, nEntries);
  //  myMPT3->AddProperty("EFFICIENCY",    PhotonEnergy, Efficiency3, nEntries);  

  Mirror->SetMaterialPropertiesTable(myMPT3);

//
// CATH
//

/*
 G4double Reflectivity4[nEntries] =
            { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00 };

  G4double Efficiency4[nEntries] =
           {0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25,  0.25,  0.25,  0.25, 0.25, 0.25,
           0.25, 0.25 };  */



 
  G4MaterialPropertiesTable* myMPT4 = new G4MaterialPropertiesTable();
  myMPT4->AddProperty("REFLECTIVITY",       PhotonEnergy, Reflectivity4,nEntries);
  myMPT4->AddProperty("EFFICIENCY",    PhotonEnergy, Efficiency4, nEntries);  

  CATH->SetMaterialPropertiesTable(myMPT4);

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

//
//	------------- Volumes --------------

// The detector
//
  G4Box* det_box = new G4Box("World",det_x,det_y,det_z);

  G4LogicalVolume* det_log
    = new G4LogicalVolume(det_box,Air,"World",0,0,0);

  G4VPhysicalVolume* det_phys
    = new G4PVPlacement(0,G4ThreeVector(),det_log,"World",0,false,0);

// The quartz
//	

  G4double q_yLB = quartz_y - (quartz_z);

  G4Trap* quartz_box = new G4Trap("Quartz", 2*quartz_x, 2*quartz_z, 2*quartz_y, 2*q_yLB);

  G4LogicalVolume* quartz_log
    = new G4LogicalVolume(quartz_box,Quartz,"Quartz",0,0,0);

  qsimScintDetector* quartzSD = new qsimScintDetector("QuartzSD", 10);

  SDman->AddNewDetector(quartzSD);
  quartz_log->SetSensitiveDetector(quartzSD);

  G4RotationMatrix* rotQ = new G4RotationMatrix;
    rotQ->rotateX(M_PI/2.*rad);
  
  G4VPhysicalVolume* quartz_phys
    = new G4PVPlacement(rotQ,G4ThreeVector(0.0*mm,0,quartz_zPos),quartz_log,"Quartz",
                        det_log,false,0);

//Quartz Suport Frame (need to be re-located)

// Front Plate
//
  G4Box* front_plate_box = new G4Box("front_plate", 8.89*cm,19.863*cm,0.9525*cm);

  //G4LogicalVolume* front_plate_log
    //= new G4LogicalVolume(front_plate_box,Alu_Mat,"front_plate",0,0,0);

  //G4VPhysicalVolume* front_plate_phys
    //= new G4PVPlacement(0,G4ThreeVector(15.625*cm,0*cm,0*cm),front_plate_log,"front_plate",det_log,false,0);

  G4Tubs* plate_hole = new G4Tubs("plate_hole",0*cm,3.1813*cm,0.9525*cm,0*deg,360*deg);

  G4SubtractionSolid* FrontPlateWithHole
    = new G4SubtractionSolid("FrontPlateWithHole", front_plate_box, plate_hole);

  G4LogicalVolume* FrontPlateWithHole_log
    = new G4LogicalVolume(FrontPlateWithHole,Alu_Mat,"front_plate",0,0,0);

  G4RotationMatrix* rotPlate = new G4RotationMatrix;
    rotPlate->rotateY(90*deg);

  //G4VPhysicalVolume* FrontPlateWithHole_phys
    //= new G4PVPlacement(rotPlate,G4ThreeVector(15.625*cm,0*cm,0*cm),FrontPlateWithHole_log,"front_plate",det_log,false,0);

  //G4LogicalVolume* plate_hole_log
    //= new G4LogicalVolume(plate_hole,Air,"plate_hole",0,0,0);

  //G4RotationMatrix* rotPlate = new G4RotationMatrix;
    //rotPlate->rotateY(90*deg);

  //G4VPhysicalVolume* plate_hole_phys
    //= new G4PVPlacement(rotPlate,G4ThreeVector(15.625*cm,0*cm,0*cm),plate_hole_log,"plate_hole",det_log,false,0);


//quartz holder 

  G4Box* quartz_holder_bar = new G4Box("quartz_holder",7.0*cm,0.30*cm,0.30*cm);

  G4Box* quartz_holder_bar_cut = new G4Box("quartz_holder",7.0*cm,0.15*cm,0.15*cm);

  G4SubtractionSolid* QuartzHolderLeft
    = new G4SubtractionSolid("QuartzHolder", quartz_holder_bar, quartz_holder_bar_cut,0,G4ThreeVector(0*cm,-0.15*cm,0.15*cm));

  G4LogicalVolume* quartz_holder_log_left
    = new G4LogicalVolume(QuartzHolderLeft,Alu_Mat,"quartz_holder",0,0,0);

  G4SubtractionSolid* QuartzHolderRight
    = new G4SubtractionSolid("QuartzHolder", quartz_holder_bar, quartz_holder_bar_cut,0,G4ThreeVector(0*cm,0.15*cm,0.15*cm));

  G4LogicalVolume* quartz_holder_log_right
    = new G4LogicalVolume(QuartzHolderRight,Alu_Mat,"quartz_holder",0,0,0);

  //G4LogicalVolume* quartz_holder_log
    //= new G4LogicalVolume(quartz_holder_bar,Alu_Mat,"quartz_holder",0,0,0);

  //G4VPhysicalVolume* quartz_holder_phys_left
    //= new G4PVPlacement(0,G4ThreeVector(7.5*cm,1.78*cm,-0.5*cm),quartz_holder_log_left,"quartz_holder_left",det_log,false,0);

  //G4VPhysicalVolume* quartz_holder_phys_right
    //= new G4PVPlacement(0,G4ThreeVector(7.5*cm,-1.78*cm,-0.5*cm),quartz_holder_log_right,"quartz_holder_right",det_log,false,0);

// Light Guide 

//Front

G4Box* frontPlate_1 = new G4Box("frontPlate_1", 0.025*cm, 44.73*mm/2, 20.19*mm/2);
   
	G4LogicalVolume* frontPlate_1_log = new G4LogicalVolume(frontPlate_1, Mirror, "FrontPlate_1_log",0,0,0);
	
	G4RotationMatrix* rotF_1 = new G4RotationMatrix;
	rotF_1->rotateY(0*deg);
	 G4VPhysicalVolume* frontPlate_1_phys
		= new G4PVPlacement(rotF_1,G4ThreeVector(-76.10*mm,0.17*mm,1.87*mm),frontPlate_1_log,"FrontPlate_1_phys",
							det_log,false,0);

// Top

   G4Box* topPlate_1 = new G4Box("topPlate_1", 136.53*mm/2, 44.73*mm/2, 0.025*cm);
   
	G4LogicalVolume* topPlate_1_log = new G4LogicalVolume(topPlate_1, Mirror, "TopPlate_1_log",0,0,0);
	
	G4RotationMatrix* rotT_1 = new G4RotationMatrix;
	rotT_1->rotateY(0*deg);
	 G4VPhysicalVolume* topPlate_1_phys
		= new G4PVPlacement(rotT_1,G4ThreeVector(-7.83*mm,0.17*mm,-8.73*mm),topPlate_1_log,"TopPlate_1_phys",
							det_log,false,0);

   G4Box* topPlate_2 = new G4Box("topPlate_2", 12.38*mm/2, 44.73*mm/2, 0.025*cm);
   
	G4LogicalVolume* topPlate_2_log = new G4LogicalVolume(topPlate_2, Mirror, "TopPlate_2_log",0,0,0);
	
	G4RotationMatrix* rotT_2 = new G4RotationMatrix;
	rotT_2->rotateY(-255*deg);
	 G4VPhysicalVolume* topPlate_2_phys
		= new G4PVPlacement(rotT_2,G4ThreeVector(62.09*mm,0.17*mm,-14.49*mm),topPlate_2_log,"TopPlate_2_phys",
							det_log,false,0);

   G4Box* topPlate_3 = new G4Box("topPlate_3", 38.58*mm/2, 44.73*mm/2, 0.025*cm);
   
	G4LogicalVolume* topPlate_3_log = new G4LogicalVolume(topPlate_3, Mirror, "TopPlate_3_log",0,0,0);
	
	G4RotationMatrix* rotT_3 = new G4RotationMatrix;
	rotT_3->rotateY(-45*deg);
	 G4VPhysicalVolume* topPlate_3_phys
		= new G4PVPlacement(rotT_3,G4ThreeVector(77.36*mm,-0.17*mm,-34.26*mm),topPlate_3_log,"TopPlate_3_phys",
							det_log,false,0);

// Bottom 
   G4Box* botPlate_1 = new G4Box("botPlate_1", 162.43*mm/2, 44.73*mm/2, 0.025*cm);
   
	G4LogicalVolume* botPlate_1_log = new G4LogicalVolume(botPlate_1, Mirror, "botPlate_1_log",0,0,0);
	
	G4RotationMatrix* rotB_1 = new G4RotationMatrix;
	rotB_1->rotateY(0*deg);
        rotB_1->rotateY(0*deg); 	

	 G4VPhysicalVolume* botPlate_1_phys
		= new G4PVPlacement(rotB_1,G4ThreeVector(5.12*mm, 0.17*mm,11.97*mm),botPlate_1_log,"botPlate_1_phys",
			det_log,false,0);

   G4Box* botPlate_2 = new G4Box("botPlate_2", 45.62*mm/2, 44.73*mm/2, 0.025*cm);
   
	G4LogicalVolume* botPlate_2_log = new G4LogicalVolume(botPlate_2, Mirror, "botPlate_2_log",0,0,0);
	
	G4RotationMatrix* rotB_2 = new G4RotationMatrix;
	rotB_2->rotateY(-45*deg); 	

	 G4VPhysicalVolume* botPlate_2_phys
		= new G4PVPlacement(rotB_2,G4ThreeVector(102.46*mm, 0.17*mm,-4.16*mm),botPlate_2_log,"botPlate_2_phys",
			det_log,false,0);

//  Laterals

//right		
										
    G4Trap* RPlate_1 = new G4Trap("RPlate_1", 0.05*cm, 20.19*mm, 162.43*mm, 137.03*mm);
	  
	G4LogicalVolume* RPlate_1_log = new G4LogicalVolume(RPlate_1, Mirror, "RPlate_1_log",0,0,0);
	
	G4RotationMatrix* rotR_1 = new G4RotationMatrix;
	rotR_1->rotateY(0*deg);
	rotR_1->rotateX(90*deg);

	 G4VPhysicalVolume* RPlate_1_phys
		= new G4PVPlacement(rotR_1,G4ThreeVector(-1.77*mm,22.53*mm,1.01*mm),RPlate_1_log,"RPlate_1_phys",
							det_log,false,0);

    G4Trap* RPlate_2 = new G4Trap("RPlate_2", 0.05*cm, 38.58*mm, 45.62*mm, 38.43*mm);
	  
	G4LogicalVolume* RPlate_2_log = new G4LogicalVolume(RPlate_2, Mirror, "RPlate_2_log",0,0,0);
	
	G4RotationMatrix* rotR_2 = new G4RotationMatrix;
        rotR_2->rotateX(90*deg);
	rotR_2->rotateY(180*deg);
	rotR_2->rotateZ(45*deg);
        
	 G4VPhysicalVolume* RPlate_2_phys
		= new G4PVPlacement(rotR_2,G4ThreeVector(89.09*mm+0.5*mm,22.53*mm,-20.23*mm+0.5*mm),RPlate_2_log,"RPlate_2_phys",
							det_log,false,0);


    G4int nCVtx = 8;
    std::vector<G4TwoVector> cvtx(nCVtx);
    cvtx[0] = G4TwoVector(  0.0*cm,   0.0*mm);
    cvtx[1] = G4TwoVector(-4.93*mm,  11.35*mm);
    cvtx[2] = G4TwoVector( 32.45*mm,   0.0*mm);
    cvtx[3] = G4TwoVector( 32.45*mm,   0.0*mm);
    cvtx[4] = G4TwoVector(  0.0*mm,   0.0*mm);
    cvtx[5] = G4TwoVector(-4.93*mm,  11.35*mm);
    cvtx[6] = G4TwoVector( 32.45*mm,   0.0*mm);
    cvtx[7] = G4TwoVector( 32.45*mm,   0.0*mm);
     
    G4GenericTrap* RPlate_3 = new G4GenericTrap("RPlate_3",0.025*cm,cvtx);
    
    G4LogicalVolume* RPlate_3_log = new G4LogicalVolume(RPlate_3, Mirror, "RPlate_3_log",0,0,0);
    
    G4RotationMatrix* rotR_3 = new G4RotationMatrix;
        rotR_3->rotateX(90*deg);
	rotR_3->rotateY(0*deg);
	rotR_3->rotateZ(38.48*deg);
        
	 G4VPhysicalVolume* RPlate_3_phys
		= new G4PVPlacement(rotR_3,G4ThreeVector(60.94*mm-0.5*mm,22.53*mm,-8.23*mm-0.5*mm),RPlate_3_log,"RPlate_3_phys",
							det_log,false,0);

//Left
	  
	G4LogicalVolume* LPlate_1_log = new G4LogicalVolume(RPlate_1, Mirror, "LPlate_log",0,0,0);
	
	G4RotationMatrix* rotL_1 = new G4RotationMatrix;
	rotL_1->rotateY(0*deg);
	rotL_1->rotateX(90*deg);
	
	 G4VPhysicalVolume* LPlate_1_phys
		= new G4PVPlacement(rotL_1,G4ThreeVector(-1.77*mm,-22.53*mm,1.01*mm),LPlate_1_log,"RPlate_phys",
							det_log,false,0);

	  
	G4LogicalVolume* LPlate_2_log = new G4LogicalVolume(RPlate_2, Mirror, "LPlate_2_log",0,0,0);
	
	G4RotationMatrix* rotL_2 = new G4RotationMatrix;
        rotL_2->rotateX(90*deg);
	rotL_2->rotateY(180*deg);
	rotL_2->rotateZ(45*deg);
        
	 G4VPhysicalVolume* LPlate_2_phys
		= new G4PVPlacement(rotL_2,G4ThreeVector(89.09*mm+0.5*mm,-22.53*mm,-20.23*mm+0.5*mm),LPlate_2_log,"LPlate_2_phys",
							det_log,false,0);

    G4GenericTrap* LPlate_3 = new G4GenericTrap("LPlate_3",0.025*cm,cvtx);
    
    G4LogicalVolume* LPlate_3_log = new G4LogicalVolume(LPlate_3, Mirror, "LPlate_3_log",0,0,0);
    
    G4RotationMatrix* rotL_3 = new G4RotationMatrix;
        rotL_3->rotateX(90*deg);
	rotL_3->rotateY(0*deg);
	rotL_3->rotateZ(39.61*deg);

        
	 G4VPhysicalVolume* LPlate_3_phys
		= new G4PVPlacement(rotL_3,G4ThreeVector(60.94*mm-0.25*mm,-22.53*mm,-8.23*mm-0.25*mm),LPlate_3_log,"LPlate_3_phys",
							det_log,false,0);


//PMT

  // Rotation

    G4double rtphi = -45.0*deg;
    G4RotationMatrix rm;
    rm.rotateY(rtphi);	
    G4double anini = 0*deg;
    G4double anspan = 360*deg;	

// The photomultiplier window dimemsions

    G4double prin = 0;
    G4double prout = 2.3*cm;
    G4double plngth = 1.5*mm;    
    
  //pmt
  
  G4Tubs* pmt = new G4Tubs("PMT",prin,prout,plngth,anini,anspan);

  G4LogicalVolume* pmt_log
    = new G4LogicalVolume(pmt,Air,"PMT",0,0,0);

  // Make PMT Sensitive
	
  
  G4String DetSDname = "tracker1";

  qsimDetector* trackerSD = new qsimDetector(DetSDname, 1);
  
  SDman->AddNewDetector(trackerSD);
  pmt_log->SetSensitiveDetector(trackerSD);

  G4VPhysicalVolume* pmt_phys;

    pmt_phys = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(110.7*mm,0.04*mm,-38.60*mm)),  
                        pmt_log,"PMT",
                        det_log,false,0);	
	
 // Coincidence volumes **** NOTE: Upper scint is above the quartz (First coincidence w/ e-)
 
	G4Box* upperScint = new G4Box("upperScint",0.25*cm,3.5*cm,20*cm);
	G4LogicalVolume* uScint_log = new G4LogicalVolume(upperScint,Air,"upperScint",0,0,0);

	// Make sensitive
			//G4String 
	DetSDname = "tracker3";

	qsimScintDetector* upScintSD = new qsimScintDetector(DetSDname, 1);
  
	SDman->AddNewDetector(upScintSD);
	uScint_log->SetSensitiveDetector(upScintSD);
	
	G4double scintAngle = 90.0*deg;

	G4RotationMatrix* scintRoll = new G4RotationMatrix;
	scintRoll->rotateY(scintAngle);
		
	G4PVPlacement* uScint_phys;

        uScint_phys 
		= new G4PVPlacement(scintRoll,G4ThreeVector(0.0*cm,0.0*cm,55.0*cm),
							uScint_log,"upperScint",det_log,false,0);
								

 /////////////
 
	G4Box* lowScint = new G4Box("lowScint",0.25*cm,3.5*cm,20*cm);
	G4LogicalVolume* lScint_log = new G4LogicalVolume(lowScint,Air,"lowScint",0,0,0);
	
	
	DetSDname = "tracker2";
	 
	 qsimScintDetector* loScintSD = new qsimScintDetector(DetSDname, 2);
	 
	 SDman->AddNewDetector(loScintSD);
	 lScint_log->SetSensitiveDetector(loScintSD);
	 	
	G4PVPlacement* lScint_phys;

        lScint_phys 
		= new G4PVPlacement(scintRoll,G4ThreeVector(0.0*cm,0.0*cm,-55.0*cm),
	     						lScint_log,"lowerScint",det_log,false,0);


 /////////////
 
	G4Box* Pb_blox = new G4Box("Pb_blox",13.0*cm,10.0*cm,20.0*cm);

	G4LogicalVolume* Pb_log = new G4LogicalVolume(Pb_blox,Pb_Mat,"Lead",0,0,0);
		
	G4PVPlacement* Pb_phys;

        Pb_phys 
		= new G4PVPlacement(scintRoll,G4ThreeVector(0.0*cm,0.0*cm,-41*cm),
							Pb_log,"Pb",det_log,false,0);
  


//	------------- Surfaces --------------
//
// Quartz
//
  G4OpticalSurface* OpQuartzSurface = new G4OpticalSurface("QuartzSurface");
  OpQuartzSurface->SetType(dielectric_dielectric);
  OpQuartzSurface->SetFinish(ground);
  OpQuartzSurface->SetModel(unified);

  //  G4LogicalBorderSurface* QuartzSurface = 
  //                                 new G4LogicalBorderSurface("QuartzSurface",
  //                                 quartz_phys,det_phys,OpQuartzSurface);

  //  if(QuartzSurface->GetVolume1() == quartz_phys) G4cout << "Equal" << G4endl;
  //  if(QuartzSurface->GetVolume2() == det_phys  ) G4cout << "Equal" << G4endl;

// Mirrors and cathode

  G4OpticalSurface* MOpSurface = new G4OpticalSurface("MirrorOpSurface");
  G4OpticalSurface* CTHOpSurface = new G4OpticalSurface("CathodeOpSurface");

  MOpSurface -> SetType(dielectric_metal);
  MOpSurface -> SetFinish(ground);
  MOpSurface -> SetModel(glisur);

  //  G4double polish = 0.8;

/*  G4double Reflectivity3[nEntries] =
            { 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
              0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
              0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
              0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
              0.90, 0.90, 0.90, 0.90 };  */


  G4MaterialPropertiesTable* MOpSurfaceProperty = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable* COpSurfaceProperty = new G4MaterialPropertiesTable();

  MOpSurfaceProperty -> AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity3,nEntries);

  MOpSurface -> SetMaterialPropertiesTable(MOpSurfaceProperty);

  COpSurfaceProperty -> AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity4,nEntries);
  COpSurfaceProperty -> AddProperty("EFFICIENCY",PhotonEnergy,Efficiency4,nEntries);

  CTHOpSurface -> SetMaterialPropertiesTable(COpSurfaceProperty);
 

  //G4LogicalSkinSurface* TSurface = new
               //G4LogicalSkinSurface("TMirrorOpS",tmirror_log,MOpSurface);

	G4LogicalSkinSurface* FrontSurface_1 = new
	G4LogicalSkinSurface("FrontMirrorOpS_1",frontPlate_1_log,MOpSurface);

	G4LogicalSkinSurface* TopSurface_1 = new
	G4LogicalSkinSurface("TopMirrorOpS_1",topPlate_1_log,MOpSurface);

	G4LogicalSkinSurface* TopSurface_2 = new
	G4LogicalSkinSurface("TopMirrorOpS_2",topPlate_2_log,MOpSurface);

	G4LogicalSkinSurface* TopSurface_3 = new
	G4LogicalSkinSurface("TopMirrorOpS_3",topPlate_3_log,MOpSurface);
	
	G4LogicalSkinSurface* BotSurface_1 = new
	G4LogicalSkinSurface("BotMirrorOpS_1",botPlate_1_log,MOpSurface);

	G4LogicalSkinSurface* BotSurface_2 = new
	G4LogicalSkinSurface("BotMirrorOpS_2",botPlate_2_log,MOpSurface);

	G4LogicalSkinSurface* LSurface_1 = new
	G4LogicalSkinSurface("LMirrorOpS_1",LPlate_1_log,MOpSurface);

        G4LogicalSkinSurface* LSurface_2 = new
	G4LogicalSkinSurface("LMirrorOpS_2",LPlate_2_log,MOpSurface);

        G4LogicalSkinSurface* LSurface_3 = new
	G4LogicalSkinSurface("LMirrorOpS_3",LPlate_3_log,MOpSurface);	

	G4LogicalSkinSurface* RSurface_1 = new
	G4LogicalSkinSurface("RMirrorOpS_1",RPlate_1_log,MOpSurface);

        G4LogicalSkinSurface* RSurface_2 = new
	G4LogicalSkinSurface("RMirrorOpS_2",RPlate_2_log,MOpSurface);

        G4LogicalSkinSurface* RSurface_3 = new
	G4LogicalSkinSurface("RMirrorOpS_3",RPlate_3_log,MOpSurface);
	
  //G4LogicalSkinSurface* CSurface = new
               //G4LogicalSkinSurface("CMirrorOpS",cmirror_log,MOpSurface);

  //G4LogicalSkinSurface* CathSurface = new
               //G4LogicalSkinSurface("CathOpS",cath_log,CTHOpSurface);


  //G4LogicalSkinSurface* QuartWinSurface = new 
               //G4LogicalSkinSurface("QuartzWinOpS", QuartzWin_log, OpQuartzSurface);





//
// Generate & Add Material Properties Table attached to the optical surfaces
//
  const G4int num = 2;
  G4double Ephoton[num] = {2.038*eV, 4.144*eV};

  //OpticalQuartzSurface 
  G4double RefractiveIndex[num] = {1.46, 1.46};
  G4double SpecularLobe[num]    = {0.3, 0.3};
  G4double SpecularSpike[num]   = {0.2, 0.2};
  G4double Backscatter[num]     = {0.2, 0.2};

  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();
  
  myST1->AddProperty("RINDEX",                Ephoton, RefractiveIndex, num);
  myST1->AddProperty("SPECULARLOBECONSTANT",  Ephoton, SpecularLobe,    num);
  myST1->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike,   num);
  myST1->AddProperty("BACKSCATTERCONSTANT",   Ephoton, Backscatter,     num);

  OpQuartzSurface->SetMaterialPropertiesTable(myST1);




//always return the physical World
  return det_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
