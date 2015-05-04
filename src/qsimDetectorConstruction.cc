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
#include "G4VisAttributes.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

qsimDetectorConstruction::qsimDetectorConstruction()
{
  det_x = det_y = det_z = 275*cm;
  quartz_x = 1.75*cm; 
  quartz_y = 7.5*cm;  //2.5 
//Half cm
    quartz_z = 0.5*cm;
//One cm
//  quartz_z = 0.5*cm;

	quartz_zPos = 0.*cm;//-1.1*cm; //-.9*cm; //-.6*cm;

  cone_rmin1 = 2.1*cm;
  cone_rmax1 = cone_rmin1+.05*cm;
  cone_rmin2 = 2.5*cm;  // normally 2.5*cm;
  cone_rmax2 = cone_rmin2+.05*cm;
  cone_z = quartz_y+.5*cm;    //3
  cone_sphi = 0.;
  cone_fphi = 2*3.1415;

  rin = cone_rmin2;  // normally 2.5*cm;
  rout = rin+.05*cm;
  lngth = 1.9*cm;  // PMT dist. = 2*lngth +1cm  (10.4 == 4.5, 6.8 == 2.9)


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


/*
// Al plate
	G4Box* Al_box = new G4Box("Al_box", 20*cm, 20*cm, .0176*cm);
	
  G4LogicalVolume* Al_log
    = new G4LogicalVolume(Al_box,Alu_Mat,"Al_log",0,0,0);
	
	G4VPhysicalVolume* Al_phys = new G4PVPlacement(0, G4ThreeVector(0,0,-9.0*cm), Al_log, "Aluminum",
									det_log, false, 0);
*/

// The quartz
//	

  G4double q_yLB = quartz_y - (quartz_z);

  G4Trap* quartz_box = new G4Trap("Quartz", 2*quartz_x, 2*quartz_z, 4*quartz_y, 4*q_yLB);

//  G4Box* quartz_box = new G4Box("Quartz",quartz_x,quartz_y,quartz_z);

  G4LogicalVolume* quartz_log
    = new G4LogicalVolume(quartz_box,Quartz,"Quartz",0,0,0);

  qsimScintDetector* quartzSD = new qsimScintDetector("QuartzSD", 10);

  SDman->AddNewDetector(quartzSD);
  quartz_log->SetSensitiveDetector(quartzSD);

  G4RotationMatrix* rotQ = new G4RotationMatrix;
//	rotQ->rotateZ(0.*rad);
	
    rotQ->rotateX(-90.*deg);
    rotQ->rotateY(0*rad);
    rotQ->rotateZ(135.*deg);

  G4VPhysicalVolume* quartz_phys
    = new G4PVPlacement(rotQ,G4ThreeVector(0.0*cm,0,quartz_zPos),quartz_log,"Quartz",
                        det_log,false,0);  // normally zero vector

// Front Plate
//
  G4Box* front_plate_box = new G4Box("front_plate", 8.89*cm,19.863*cm,0.9525*cm);

  G4Tubs* plate_hole = new G4Tubs("plate_hole",0*cm,3.1813*cm,0.9525*cm,0*deg,360*deg);

  G4SubtractionSolid* FrontPlateWithHole
    = new G4SubtractionSolid("FrontPlateWithHole", front_plate_box, plate_hole);

  G4LogicalVolume* FrontPlateWithHole_log
    = new G4LogicalVolume(FrontPlateWithHole,Alu_Mat,"front_plate",0,0,0);

  G4RotationMatrix* rotPlate = new G4RotationMatrix;
    rotPlate->rotateY(90*deg);

  //G4VPhysicalVolume* FrontPlateWithHole_phys
    //= new G4PVPlacement(rotPlate,G4ThreeVector(15.625*cm,0*cm,0*cm),FrontPlateWithHole_log,"front_plate",det_log,false,0);


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

  //G4VPhysicalVolume* quartz_holder_phys_left
    //= new G4PVPlacement(0,G4ThreeVector(7.5*cm,1.78*cm,-0.5*cm),quartz_holder_log_left,"quartz_holder_left",det_log,false,0);

  //G4VPhysicalVolume* quartz_holder_phys_right
    //= new G4PVPlacement(0,G4ThreeVector(7.5*cm,-1.78*cm,-0.5*cm),quartz_holder_log_right,"quartz_holder_right",det_log,false,0);

	//////////////////////Air
	
/*	  G4Box* ruler_box = new G4Box("Ruler",1.5*cm,1*cm,1*cm);
	 
	 G4LogicalVolume* ruler_log
	 = new G4LogicalVolume(ruler_box,Air,"Ruler",0,0,0);
	 
	 G4VPhysicalVolume* ruler_phys
	 = new G4PVPlacement(0,G4ThreeVector(quartz_y+1.4*cm,0.5*cm,0*cm),ruler_log,"Ruler",
	 det_log,false,0);  
*/	 
	 
	
// The small mirror on the quartz
//
/*
  G4Box* mirr_Box = new G4Box("QuMirror", .05*cm, quartz_x, quartz_z*1.4142);
  
  G4LogicalVolume* mirr_log = new G4LogicalVolume(mirr_Box, Mirror, "QuMirror",0,0,0);
  
    G4RotationMatrix* rotM = new G4RotationMatrix;

	//rotM->rotateZ(M_PI*rad);
    rotM->rotateY(-M_PI/4.*rad);
  
  
  G4VPhysicalVolume* mirr_Phys = new G4PVPlacement(rotM,G4ThreeVector(-1*quartz_y+.05*cm,0,quartz_zPos+.06*cm), mirr_log, "QuMirror",
													det_log,false, 0);  // normally z= .06
*/

// Trapezoid tube plates
// Top

   G4Trd* topPlate = new G4Trd("topPlate", 1.*mm, 1.*mm, 38.*mm/2, 64.*mm/2, 228.*mm);
   
	G4LogicalVolume* topPlate_log = new G4LogicalVolume(topPlate, Mirror, "TopPlate_log",0,0,0);


   G4VisAttributes* topPlateVis = new G4VisAttributes(G4Color(1.,0.,0.));
   topPlate_log->SetVisAttributes(topPlateVis);


	G4RotationMatrix* rotT = new G4RotationMatrix;
	rotT->rotateX(0.0*deg);
	rotT->rotateY(-45.0*deg+5.27*deg);//12.27*deg
	 G4VPhysicalVolume* topPlate_phys
		= new G4PVPlacement(rotT,G4ThreeVector(-32.*mm*sin(45*deg)+5.*mm*sin(45*deg),0.0*mm,-32.*mm*sin(45*deg)-5.*mm*sin(45*deg)),topPlate_log,"TopPlate_phys",
//		= new G4PVPlacement(rotT,G4ThreeVector(-13.44*mm,0.0*mm,-41.96*mm),topPlate_log,"TopPlate_phys",
							det_log,false,0);  // normally zero vector


// Bottom 
   G4Trd* botPlate = new G4Trd("botPlate", 1.*mm, 1.*mm, 38.*mm/2, 64.*mm/2, 225.*mm);
   
	G4LogicalVolume* botPlate_log = new G4LogicalVolume(botPlate, Mirror, "botPlate_log",0,0,0);
	
	G4RotationMatrix* rotB = new G4RotationMatrix;
	rotB->rotateX(0.0*deg);
	rotB->rotateY(-45.0*deg-9.42*deg);	

	 G4VPhysicalVolume* botPlate_phys
		= new G4PVPlacement(rotB,G4ThreeVector(-32.*mm*sin(45*deg)-5.*mm*sin(45*deg),0.0*mm,-32.*mm*sin(45*deg)+5.*mm*sin(45*deg)),botPlate_log,"botPlate_phys",
			det_log,false,0);  // normally zero vector


												
    G4int nCVtx = 8;
    std::vector<G4TwoVector> cvtx(nCVtx);
    cvtx[0] = G4TwoVector(0.0*mm, 0.0*mm);
    cvtx[1] = G4TwoVector(228.43*mm, 0.0*mm);
    cvtx[2] = G4TwoVector(228.43*mm, 10.84*mm);
    cvtx[3] = G4TwoVector(0.0*mm, 60.0*mm);
    cvtx[4] = G4TwoVector(0.0*mm, 0.0*mm);
    cvtx[5] = G4TwoVector(228.43*mm, 0.0*mm);
    cvtx[6] = G4TwoVector(228.43*mm, 10.84*mm);
    cvtx[7] = G4TwoVector(0.0*mm, 60.0*mm);
     
    G4GenericTrap* RPlate = new G4GenericTrap("RPlate",0.025*cm,cvtx);
	  
	G4LogicalVolume* RPlate_log = new G4LogicalVolume(RPlate, Mirror, "RPlate_log",0,0,0);
	
	G4RotationMatrix* rotR = new G4RotationMatrix;
	rotR->rotateX(-90.0*deg+1.5*deg);
	rotR->rotateZ(135.0*deg);
        rotR->rotateY(0.0*deg);

//	 G4VPhysicalVolume* RPlate_phys
//	 = new G4PVPlacement(rotR,G4ThreeVector(42.43*mm,25.44*mm,63.64*mm),RPlate_log,"RPlate_phys",
//							det_log,false,0);  // normally zero vector
	  
	G4LogicalVolume* LPlate_log = new G4LogicalVolume(RPlate, Mirror, "LPlate_log",0,0,0);
	
	G4RotationMatrix* rotL = new G4RotationMatrix;
	rotL->rotateX(-90.0*deg-1.5*deg);
	rotL->rotateZ(135.0*deg);
	
//	 G4VPhysicalVolume* LPlate_phys
//		= new G4PVPlacement(rotL,G4ThreeVector(42.43*mm,-25.44*mm,63.64*mm),LPlate_log,"LPlate_phys",
//							det_log,false,0);  // normally zero vector

// The cone mirror
//
	
  G4Cons* cmirror_cone = new G4Cons("CMirror",cone_rmin1,cone_rmax1,
  cone_rmin2,cone_rmax2,cone_z,cone_sphi,cone_fphi);

  G4LogicalVolume* cmirror_log
    = new G4LogicalVolume(cmirror_cone,Mirror,"CMirror",0,0,0);

  // Rotation

  G4VPhysicalVolume* cmirror_phys;

    G4double rtphi = 45.*deg;
    G4RotationMatrix rm;
    rm.rotateY(rtphi);
    G4double tmpx = 0.5*cm;

//    cmirror_phys
//    = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(tmpx,0., 0.)),  
//                        cmirror_log,"CMirror",
//                        det_log,false,0);


// The tube mirror
//	//	//	//	//	//	//	//	//	//	

    G4double rin = 2.5*cm;
    G4double rout = 2.55*cm;
    G4double lngth = 2.3*cm;  // reg. tube
//    G4double lngth = 3.5*cm;  // long tube	
    G4double anini = 0*deg;
    G4double anspan = 360*deg;
/*


G4Cons* mirror_tube = new G4Cons("TMirror",cone_rmin2,cone_rmax2,
  2.5*cm,2.55*cm,lngth,cone_sphi,cone_fphi);

  G4LogicalVolume* tmirror_log
    = new G4LogicalVolume(mirror_tube,Mirror,"TMirror",0,0,0);


  // Rotation

  G4VPhysicalVolume* tmirror_phys;

 
    G4double tmp = lngth+cone_z+tmpx;

    tmirror_phys      // Only for tube at full distance (7.1cm)
   // = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(11.8*cm,0.,0.)),  // Long tube
    = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(12.4*cm,0.,0.)),  // Long tube
						tmirror_log,"TMirror",
                        det_log,false,0);
*/
 
 /*   
  // Rotation

  G4VPhysicalVolume* cmirror_phys;

 

    rm.rotateY(rtphi);



*/
 //	//	//	//	//	//	//	//	//	//	

// The photomultiplier
//	quartz window

    G4double prin = 0;
    G4double prout = 2.6*cm;
    G4double plngth = 1.5*mm;
    
/*
  G4Tubs* quartz_window = new G4Tubs("QuartzWin",prin,prout,3*mm,anini,anspan);

  G4LogicalVolume* QuartzWin_log
  = new G4LogicalVolume(quartz_window,Quartz,"QuartzWin",0,0,0);

  G4VPhysicalVolume* QuartzWin_phys;

  QuartzWin_phys = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector( 15.150*cm, 0., 0.)), QuartzWin_log,"QuartzWin", det_log, false, 0);
*/

    
  //pmt
  
  G4Tubs* pmt = new G4Tubs("PMT",prin,2.3*cm,plngth,anini,anspan);

  G4LogicalVolume* pmt_log
    = new G4LogicalVolume(pmt,Air,"PMT",0,0,0);

  // Make PMT Sensitive
	
  
  G4String DetSDname = "tracker1";

  qsimDetector* trackerSD = new qsimDetector(DetSDname, 1);
  
  SDman->AddNewDetector(trackerSD);
  pmt_log->SetSensitiveDetector(trackerSD);

  // Rotation


//  G4double sep = 0.5*cm;
//    G4double ptmp = tmp+lngth+sep;


  G4double pmt_z = 14.65*cm*sin(45*deg)		//Get to end of quartz
  	+ 7.*cm*sin(45*deg)		//7 cm PMT-quartz separation	
	+ 8.333*mm*sin(45*deg);		//Shift PMT off quartz center

  G4double pmt_x = 14.65*cm*sin(45*deg)		//Get to end of quartz
  	+ 7.*cm*sin(45*deg)		//7 cm PMT-quartz separation	
	- 8.333*mm*sin(45*deg);		//Shift PMT off quartz center

//    pmt_phys
      //  = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(10.1*cm,0.,0.)),  // Original sim. position
	  //    = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(10.7*cm,0.,0.)),  // Cosmic test sim. position
	
  G4VPhysicalVolume* pmt_phys
	= new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(pmt_x,0.,pmt_z)),  // Final detector-pmt length position
//		  = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(63.64*mm+0.5*cm*sin(45*deg),0.,42.43*mm+0.5*cm*sin(45*deg))),  // Final detector-pmt length position
                        pmt_log,"PMT",
                        det_log,false,0);

 // metal cathode
/*
    G4double cin = 0;
    G4double cout = 2.6*cm;
    G4double clngth = 0.1*mm;


  G4Tubs* cath = new G4Tubs("CATH",cin,cout,clngth,anin,anspan);

  G4LogicalVolume* cath_log
    = new G4LogicalVolume(cath,CATH,"CATH",0,0,0);

  qsimDetector* cathSD = new qsimDetector("cath", 2);
  
  SDman->AddNewDetector(cathSD);
  cath_log->SetSensitiveDetector(cathSD);

  // Rotation

  G4VPhysicalVolume* cath_phys;

    G4double ctmp = 10.1*cm+plngth;

    cath_phys
    //    = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(ctmp,0.,0.)),cath_log,"CATH",det_log,false,0);  // Original sim. position
	//	  = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(10.85*cm,0.,0.)),cath_log,"CATH",det_log,false,0);  // Cosmic test position (3.7cm quartz-pmt)
		  = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(14.25*cm,0.,0.)),cath_log,"CATH",det_log,false,0);  // Cosmic test position (7.1cm quartz-pmt)
*/	
	
 // Coincidence volumes **** NOTE: Upper scint is below the quartz (First coincidence w/ e-)
 
	G4Box* upperScint = new G4Box("upperScint",0.25*cm,3.5*cm,20*cm);
	G4LogicalVolume* uScint_log = new G4LogicalVolume(upperScint,Air,"upperScint",0,0,0);

	// Make sensitive
			//G4String 
	DetSDname = "tracker3";

	qsimScintDetector* upScintSD = new qsimScintDetector(DetSDname, 1);
  
	SDman->AddNewDetector(upScintSD);
	uScint_log->SetSensitiveDetector(upScintSD);
	
	G4double scintAngle = 180.0*deg;

	G4RotationMatrix* scintRoll = new G4RotationMatrix;
	scintRoll->rotateY(scintAngle);
		
	G4double upScint_pos=(-1*quartz_z)+(-2.5*cm-(quartz_y*sin(scintAngle)))*sin(scintAngle);
		
	G4PVPlacement* uScint_phys;
	//uScint_phys 
	//	= new G4PVPlacement(scintRoll,G4ThreeVector(upScint_pos-1.*cm,0.0*cm,upScint_pos-1.*cm),
	//						uScint_log,"upperScint",det_log,false,0);

        //uScint_phys 
	//	= new G4PVPlacement(scintRoll,G4ThreeVector(-55.0*cm,0.0*cm,0.0*cm),
	//						uScint_log,"upperScint",det_log,false,0);
								

 /////////////
 
	G4Box* lowScint = new G4Box("lowScint",0.25*cm,3.5*cm,20*cm);
	G4LogicalVolume* lScint_log = new G4LogicalVolume(lowScint,Air,"lowScint",0,0,0);

    // Make sensitive
/*			//G4String 
	DetSDname = "/tracker3";

	qsimTrackerSD* loScintSD = new qsimTrackerSD(DetSDname);
  
	SDman->AddNewDetector(loScintSD);
	lScint_log->SetSensitiveDetector(loScintSD);
*/	
	
		DetSDname = "tracker2";
	 
	 qsimScintDetector* loScintSD = new qsimScintDetector(DetSDname, 2);
	 
	 SDman->AddNewDetector(loScintSD);
	 lScint_log->SetSensitiveDetector(loScintSD);
	 
	
	G4double loScint_pos=upScint_pos+70.71*cm;

//(-1*quartz_z)+(41.25*cm-(quartz_y*sin(scintAngle)))*sin(scintAngle);
		
	G4PVPlacement* lScint_phys;
	//lScint_phys 
	//	= new G4PVPlacement(scintRoll,G4ThreeVector(loScint_pos,0.0*cm,loScint_pos),
	//						lScint_log,"lowerScint",det_log,false,0);

        //lScint_phys 
	//	= new G4PVPlacement(scintRoll,G4ThreeVector(55.0*cm,0.0*cm,0.0*cm),
	  //   						lScint_log,"lowerScint",det_log,false,0);


 /////////////
 
	G4Box* Pb_blox = new G4Box("Pb_blox",13.0*cm,10.0*cm,20.0*cm);
		// Really 10x5x10 cm half-lengths, expanded to ensure nothing
		//   can hit the scint. w/o the lead.
	G4LogicalVolume* Pb_log = new G4LogicalVolume(Pb_blox,Pb_Mat,"Lead",0,0,0);

	G4double Pb_pos=loScint_pos-11.25*cm; //(-1*quartz_z)+(30.0*cm-(quartz_y*sin(scintAngle)))*sin(scintAngle);
		
	G4PVPlacement* Pb_phys;
	//Pb_phys 
	//	= new G4PVPlacement(scintRoll,G4ThreeVector(Pb_pos,0.0*cm,Pb_pos),
	//						Pb_log,"Pb",det_log,false,0);

        //Pb_phys 
	//	= new G4PVPlacement(scintRoll,G4ThreeVector(41.0*cm,0.0*cm,0.0*cm),
	//						Pb_log,"Pb",det_log,false,0);
  


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

	G4LogicalSkinSurface* TopSurface = new
	G4LogicalSkinSurface("TopMirrorOpS",topPlate_log,MOpSurface);
	
	G4LogicalSkinSurface* BotSurface = new
	G4LogicalSkinSurface("BotMirrorOpS",botPlate_log,MOpSurface);

	G4LogicalSkinSurface* LSurface = new
	G4LogicalSkinSurface("LMirrorOpS",LPlate_log,MOpSurface);
	
	G4LogicalSkinSurface* RSurface = new
	G4LogicalSkinSurface("RMirrorOpS",RPlate_log,MOpSurface);
	
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
