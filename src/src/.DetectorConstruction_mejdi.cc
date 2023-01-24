////////////////////////////////////////////////////////////////////////////////////////////
// Author: Steve Quinn                                                                    //
//                                                                                        //
// Description: This is the file where the SuN detector is built. You can follow the      //
//   syntax below to build and place your own detectors. There is also many useful        //
//   examples online.                                                                     //
//                                                                                        //
// Steps: 1. Define the elements that you need                                            //
//        2. Use these elements to define the materials for the experimental setup        //
//        3. Create an experimental room of air, vacuum, etc.                             //
//        4. Create your experimental setup and place it in the experimental room         //
//        5. Apply the color scheme you want for the optional visualization               //
//                                                                                        //
// Important: The way the code is currently set up, you are required to fill the array    //
//   called detectorName[i] with the name of the detectors you want to save in your       //
//   ROOT file.  In this example the names are "Det1" and "Det2" (see syntax below).      //
////////////////////////////////////////////////////////////////////////////////////////////

#include "DetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4NistManager.hh"

#include "G4UnitsTable.hh"
#include "globals.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include <iostream>
#include <sstream>
#include "G4String.hh"
#include "G4ios.hh"
#include <stdio.h>
#include "G4SystemOfUnits.hh"

DetectorConstruction::DetectorConstruction()
  :  NaI(0), Al(0), N78O21Ar1(0), Cr20Ni8Fe76(0), C2F4(0), C5O2H8(0),
     Si(0), Cu3Zn2(0), Cu(0), C2H4_FR4(0), C2H4_rigid(0), C2H4_flexible(0), C12H22N2O2(0), vacuum(0)
{
}

DetectorConstruction::~DetectorConstruction()
{
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();

  return ConstructDetector();
}

void DetectorConstruction::DefineMaterials()
{

  // define Parameters
     G4String name, symbol;  
     G4double a, z, n, density;           
     G4int ncomponents, natoms;

  // define Elements

     a = 22.99*g/mole;
     //G4Element* Na = new G4Element(name="Sodium" ,symbol="Na" , z= 11., a);
     G4Isotope* isoNa23 = new G4Isotope(name="23Na" ,z= 11. ,n=23, a);
     G4Isotope* isoNa24 = new G4Isotope(name="24Na" ,z= 11. ,n=24, 23.99*g/mole);
     G4Element* Na = new G4Element(name="Sodium" ,symbol="Na" ,100*perCent);
     Na->AddIsotope(isoNa23, 1);

     a = 126.90*g/mole;
     //G4Element* I = new G4Element(name="Iodine" ,symbol="I" , z= 53., a);
     G4Isotope* isoI127 = new G4Isotope(name="127I" ,z= 53. ,n=127, a);
     G4Isotope* isoI128 = new G4Isotope(name="128I" ,z= 53. ,n=128, 127.91*g/mole);
     G4Element* I = new G4Element(name="Iodine" ,symbol="I" ,100*perCent);
     I->AddIsotope(isoI127, 1);

     a = 204.38*g/mole;
     G4Element* Tl = new G4Element(name="Thalium" ,symbol="Tl" , z= 81., a);

     a = 26.982*g/mole;
     G4Element* elAl  = new G4Element(name="element_Aluminum",symbol="elAl" , z= 13., a);

     a = 14.00*g/mole;
     G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

     a = 16.00*g/mole;
     G4Element* O  = new G4Element(name="Oxygen",symbol="O" , z= 8., a);

     a = 39.95*g/mole;
     G4Element* Ar  = new G4Element(name="Argon",symbol="Ar" , z= 18., a);

     a = 51.996*g/mole;
     G4Element* Cr  = new G4Element(name="Chromium"  ,symbol="Cr" , z= 24., a);

     a = 58.69*g/mole;
     G4Element* Ni  = new G4Element(name="Nickel" ,symbol="Ni" , z= 28., a);

     a = 55.847*g/mole;
     G4Element* Fe  = new G4Element(name="Iron"  ,symbol="Fe" , z= 26., a);

     a = 12.011*g/mole;
     G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

     a = 18.998*g/mole;
     G4Element* F  = new G4Element(name="Fluorine"  ,symbol="F" , z= 9., a);

     a = 1.008*g/mole;
     G4Element* H  = new G4Element(name="Hydrogen"  ,symbol="H" , z= 1., a);

     // define more elements
     
     a = 28.085*g/mole;
     G4Element* elSi = new G4Element(name = "element_Silicon", symbol = "elSi", z = 14., a);

     a = 63.546*g/mole;
     G4Element* elCu = new G4Element(name = "element_Copper", symbol = "elCu", z = 29., a);

     a = 65.38*g/mole;
     G4Element* Zn = new G4Element(name = "Zinc", symbol = "Zn", z = 30., a);

  // define Materials

   //..........Stainless Steel..........
     density = 8.0*g/cm3;
     Cr20Ni8Fe76 = new G4Material(name="Stainless_Steel", density, ncomponents=3);
     Cr20Ni8Fe76->AddElement(Cr, natoms=20);
     Cr20Ni8Fe76->AddElement(Fe, natoms=76);
     Cr20Ni8Fe76->AddElement(Ni, natoms=8);

   //..........Polytetrafluorine (PTFE)............
     density = 2.20*g/cm3;
     C2F4 = new G4Material(name="PTFE", density, ncomponents=2);
     C2F4->AddElement(C, natoms=2);
     C2F4->AddElement(F, natoms=4); 

   //..........NaI......................
     density = 3.67*g/cm3;
     NaI = new G4Material(name="Sodium Iodine", density, ncomponents=3);
     NaI->AddElement(Na, natoms=1000);
     NaI->AddElement(I, natoms=1000);
     NaI->AddElement(Tl, natoms=1);

   //..........Al.......................	
     density = 2.698*g/cm3;
     Al = new G4Material (name="Aluminum", density, ncomponents=1);
     Al->AddElement(elAl, natoms=1);
	
	
   //..........Air.....................
     density = 1.2927*mg/cm3;
     N78O21Ar1 = new G4Material (name="Air", density, ncomponents=3);
     N78O21Ar1->AddElement(N, natoms=78);
     N78O21Ar1->AddElement(O, natoms=21);
     N78O21Ar1->AddElement(Ar, natoms=1);

   //..........Acrylic.....................
     density = 1.18*g/cm3;
     C5O2H8 = new G4Material (name="Acrylic", density, ncomponents=3);
     C5O2H8->AddElement(C, natoms=5);
     C5O2H8->AddElement(O, natoms=2);
     C5O2H8->AddElement(H, natoms=8);

     // define more materials

     //..........Silicon..................
     density = 2.3290*g/cm3;
     Si = new G4Material(name = "Silicon", density, ncomponents = 1);
     Si->AddElement(elSi, natoms=1);

     //..........Brass....................
     density = 8.5*g/cm3;
     Cu3Zn2 = new G4Material(name = "Brass", density, ncomponents = 2);
     Cu3Zn2->AddElement(elCu, natoms=3);
     Cu3Zn2->AddElement(Zn, natoms=2);
     
     //..........Copper...................
     density = 8.96*g/cm3;
     Cu = new G4Material(name = "Copper", density, ncomponents = 1);
     Cu->AddElement(elCu, natoms=1);

     //..........Polyethylene...................
     // For the PCB FR4 frame around the silicon of the DSSD (since I cannot find the molecular formula for FR4 anywhere)
     density = 1.850*g/cm3; // density of FR4
     C2H4_FR4 = new G4Material(name = "PE_FR4", density, ncomponents = 2);
     C2H4_FR4->AddElement(C, natoms=2);
     C2H4_FR4->AddElement(H, natoms=4);

     //..........Polyethylene...................
     // For the plastic washers
     density = 0.97*g/cm3; // rigid density of polyethylene
     C2H4_rigid = new G4Material(name = "PE_rigid", density, ncomponents = 2);
     C2H4_rigid->AddElement(C, natoms=2);
     C2H4_rigid->AddElement(H, natoms=4);

     //..........Polyethylene...................
     // For the flexible ribbon cables of the DSSD
     density = 0.93*g/cm3; // flexible density of polyethylene
     C2H4_flexible = new G4Material(name = "PE_flexible", density, ncomponents = 2);
     C2H4_flexible->AddElement(C, natoms=2);
     C2H4_flexible->AddElement(H, natoms=4);
     
     //..........Nylon...................
     // For the silicon detector holder
     density = 1.15*g/cm3;
     C12H22N2O2 = new G4Material(name = "Nylon", density, ncomponents = 4);
     C12H22N2O2->AddElement(C, natoms=12);
     C12H22N2O2->AddElement(H, natoms=22);
     C12H22N2O2->AddElement(N, natoms=2);
     C12H22N2O2->AddElement(O, natoms=2);

     //..........Vacuum...................
     // For inside the experimental room / inside of the beam pipe
     vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
     
  // Print out Elements and Materials
    G4cout << "\n\n ####-------------------------------------------------------#### \n";
    G4cout << "\n\n\n\n\t\t #### List of elements used #### \n";
    G4cout << *(G4Element::GetElementTable());
    G4cout << "\n\n\n\n\t\t #### List of materials used #### \n";
    G4cout << *(G4Material::GetMaterialTable());
    G4cout << "\n\n ####-------------------------------------------------------#### \n";
}


G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{

// EXPERIMENTAL ROOM
  G4Tubs* room_tube = new G4Tubs("room", 0.0*cm, 100.0*cm, 300.0*cm, 0.0*deg, 360.0*deg);
  G4LogicalVolume* room_log = new G4LogicalVolume(room_tube,vacuum,"room",0,0,0); // room is filled with vacuum, not air
  G4VPhysicalVolume* room_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),"room",room_log,NULL,false,0);

// BEAM PIPE
  G4double outerR_beam = 2.06375*cm;                   //edit this to change the radius of the beam pipe
  G4double innerR_beam = outerR_beam - 0.889*mm;      //edit this to change the thickness of the beam pipe
  G4double halflength_beam = (68.58/2)*cm;              //edit this to change the length of the beam pipe
  G4double startAngle_beam = 0.*deg;
  G4double spanAngle_beam = 360.*deg;

  G4Tubs* beam_tube = new G4Tubs("beam_tube",innerR_beam, outerR_beam, halflength_beam, startAngle_beam, spanAngle_beam);

  G4LogicalVolume* beam_log = new G4LogicalVolume(beam_tube,Al,"beam_log",0,0,0); // beam pipe is made out of aluminum, not stainless steel


  // The first group of lengths are based on the Pos_z from SIDES OF SUN, below. For reference:
  // Pos_z = 2.0*length_scint + 2.0*width_Al_vert + 4.0*width_Refl + halflength_Al_side;
  // Only this time it should be length_Al_side, not halflength_Al_side
  // The second group of lengths is based on my diagram
  G4VPhysicalVolume* beam_phys = new G4PVPlacement(0,G4ThreeVector(0.0*mm,0.0*mm,
								   (-1.0)*((halflength_beam) - ((2*(101.5*mm) + 2*(0.50*mm) + 4*(0.25*mm) + 13.0*mm) + ((0.9525*cm) + (0.635*cm) + (0.3175*cm) + (0.635*cm) + (0.3175*cm))))),
						   beam_log,"beam_phys",room_log,false,0);


// DIMENSIONS OF SuN                           **(changing these will scale the whole simulation)**
  G4double innerR_scint = 22.5*mm;               // 45mm borehole
  G4double outerR_scint = 203.0*mm;              // total of 406mm in diameter
  G4double length_scint = 101.5*mm;              // total of 406mm in length
  G4double width_Refl = 0.25*mm;                 // width of reflector 
  G4double width_Al_vert = 0.50*mm;              // width of alumiunum between each segment
  G4double width_Al_horiz = 0.75*mm;             // width of aluminum between top and bottom half
  G4double innerR_Al = 21.5*mm;                  // 43mm in center           
  G4double outerR_Al = 222.5*mm;                 // thick outer casing     


// NaI SCINTILLATOR
  G4double halflength_scint = 0.5*length_scint;
  G4double startAngle_scint = 0.0*deg;
  G4double spanAngle_scint = 180.0*deg;

  G4Tubs* scint_tube = new G4Tubs("scint_tube",innerR_scint,outerR_scint,halflength_scint,startAngle_scint,spanAngle_scint);
 
  G4LogicalVolume* scint_log = new G4LogicalVolume(scint_tube,NaI,"scint_log",0,0,0);

// REFLECTOR
  G4double innerR_Refl = innerR_scint - width_Refl;
  G4double outerR_Refl = outerR_scint + 2.0*width_Refl;
  G4double length_Refl = length_scint + 2.0*width_Refl;
  G4double halflength_Refl = 0.5*length_Refl;

  G4Tubs* refl_tube = new G4Tubs("refl_tube",innerR_Refl,outerR_Refl,halflength_Refl,startAngle_scint,spanAngle_scint);

  G4SubtractionSolid* refl_sub = new G4SubtractionSolid("refl_sub",refl_tube,scint_tube,0,G4ThreeVector(0.0*mm,width_Refl,0.0*mm));

  G4LogicalVolume*  refl_log = new G4LogicalVolume(refl_sub,C2F4,"refl_log",0,0,0);

// ALUMINUM          
  G4double length_Al = length_Refl + width_Al_vert;   
  G4double halflength_Al = 0.5*length_Al;

  G4Tubs* al_tube = new G4Tubs("al_tube",innerR_Al,outerR_Al,halflength_Al,startAngle_scint,spanAngle_scint);

  G4SubtractionSolid* al_sub = new G4SubtractionSolid("al_sub",al_tube,refl_tube,0,G4ThreeVector(0.0*mm,width_Al_horiz,0.0*mm));
 
  G4LogicalVolume*  al_log = new G4LogicalVolume(al_sub,Al,"al_log",0,0,0);


  G4double length_Al_side = 13.0*mm;   
  G4double halflength_Al_side = 0.5*length_Al_side;

  G4Tubs* al_tube_side = new G4Tubs("al_tube_side",innerR_Al,outerR_Al,halflength_Al_side,0.0*deg,360.0*deg);

  G4LogicalVolume*  al_log_side = new G4LogicalVolume(al_tube_side,Al,"al_log_side",0,0,0);


// ====================== Placing the volumes ==============================//

 // Rotation Matrix
    G4RotationMatrix* rot_180 = new G4RotationMatrix();
      rot_180->rotateZ(180*deg);

  G4double Pos_x = 0.0*mm;
  G4double Pos_y_Al = 0.0*mm;
  G4double Pos_y_Refl = width_Al_horiz;
  G4double Pos_y_Scint = Pos_y_Refl + width_Refl;
  G4double Pos_z = -1.5*width_Al_vert - 3.0*width_Refl - 3.0*halflength_scint;
  
  G4String topName;
  G4String botName;

for (int i=1; i<=4; i++)
 {

     if(i==1)
      {
        topName = "T1";
        botName = "B1";
      }
     if(i==2)
      {
        topName = "T2";
        botName = "B2";
      }
     if(i==3)
      {
        topName = "T3";
        botName = "B3";
      }
     if(i==4)
      {
        topName = "T4";
        botName = "B4";
      }

  // TOP OF SUN

     Pos_y_Al = 0.0*mm;
     Pos_y_Refl = width_Al_horiz;
     Pos_y_Scint = Pos_y_Refl + width_Refl;

     //aluminum
     G4VPhysicalVolume* al_top = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y_Al,Pos_z),al_log,"al_top",room_log,false,0);

     //reflector  
     G4VPhysicalVolume* refl_top = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y_Refl,Pos_z),refl_log,"refl_top",room_log,false,0);

     //scintillator
     G4VPhysicalVolume* scint_top = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y_Scint,Pos_z),scint_log,topName,room_log,false,0);


  // BOTTOM OF SUN

     Pos_y_Al = 0.0*mm;
     Pos_y_Refl = -width_Al_horiz;
     Pos_y_Scint = Pos_y_Refl - width_Refl;

     //aluminum
     G4VPhysicalVolume* al_bottom = new G4PVPlacement(rot_180,G4ThreeVector(Pos_x,Pos_y_Al,Pos_z),al_log,"al_bottom",room_log,false,0);

     //reflector  
     G4VPhysicalVolume* refl_bottom = new G4PVPlacement(rot_180,G4ThreeVector(Pos_x,Pos_y_Refl,Pos_z),refl_log,"refl_bottom",room_log,false,0);

     //scintillator
     G4VPhysicalVolume* scint_bottom = new G4PVPlacement(rot_180,G4ThreeVector(Pos_x,Pos_y_Scint,Pos_z),scint_log,botName,room_log,false,0);


     Pos_z = Pos_z + width_Al_vert + 2.0*width_Refl + length_scint;
   }


// SIDES OF SUN

  Pos_x = 0.0*mm;
  Pos_y_Al = 0.0*mm;
  Pos_z = 2.0*length_scint + 2.0*width_Al_vert + 4.0*width_Refl + halflength_Al_side; 
  G4VPhysicalVolume* al_sideA = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y_Al,Pos_z),al_log_side,"al_sideA",room_log,false,0);
  G4VPhysicalVolume* al_sideB = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y_Al,-Pos_z),al_log_side,"al_sideB",room_log,false,0);



  
  // DSSD construction variables
  // Below, cv stands for construction variable
  G4double dssd_cv_frame_diameter = 35.50*mm;
  G4double dssd_cv_frame_f = 23.20*mm;
  G4double dssd_cv_frame_b = 24.20*mm;
  G4double dssd_cv_frame_gv = 1.588*mm;
  G4double dssd_cv_frame_gh1 = 1.0*mm;
  G4double dssd_cv_frame_gh2 = 1.0*mm;
  G4double dssd_cv_frame_dch = 15.0*mm;
  G4double dssd_cv_frame_rh = 0.5*(2.40*mm);
  G4double dssd_cv_frame_dbh = 5.0*mm;
  G4double dssd_cv_frame_thickness = 3.37*mm;

  ///////////////////////////////////////////
  // Black PCB FR4 frame around the silicon of the DSSD
  // This will have portions removed from it to create the final frame
  ///////////////////////////////////////////
  
  G4double dssd_frame_whole_inner_R = 0.0*mm;
  G4double dssd_frame_whole_outer_R = 0.5*(dssd_cv_frame_diameter);
  G4double dssd_frame_whole_length = dssd_cv_frame_thickness;
  G4double dssd_frame_whole_half_length = 0.5*dssd_frame_whole_length;
  G4double dssd_frame_whole_start_angle = 0.0*deg;
  G4double dssd_frame_whole_span_angle = 360.0*deg;

  G4Tubs* dssd_frame_whole_tube = new G4Tubs("dssd_frame_whole_tube",
					     dssd_frame_whole_inner_R,
					     dssd_frame_whole_outer_R,
					     dssd_frame_whole_half_length,
					     dssd_frame_whole_start_angle,
					     dssd_frame_whole_span_angle);
  
  G4LogicalVolume* dssd_frame_whole_log = new G4LogicalVolume(dssd_frame_whole_tube,
							      C2H4_FR4,
							      "dssd_frame_whole_log", 0, 0, 0);

  /*
  // Only for visualization							      
  G4VPhysicalVolume* dssd_frame_whole = new G4PVPlacement(0, G4ThreeVector(0,0,0),
							  dssd_frame_whole_log,
							  "dssd_frame_whole", room_log, false, 0);
  */

  ///////////////////////////////////////////
  // front side cutout (left side, as seen by incoming beam)
  ///////////////////////////////////////////
  
  G4double dssd_frame_cutout_left_side_half_X = 0.5*(dssd_frame_whole_outer_R - 0.5*dssd_cv_frame_b - dssd_cv_frame_gh1 - dssd_cv_frame_gh2);
  G4double dssd_frame_cutout_left_side_half_Y = 0.5*dssd_cv_frame_f;
  G4double dssd_frame_cutout_left_side_half_Z = dssd_frame_whole_half_length + 1.0*mm; // +1.0*mm to correctly visualize the subtraction

  G4Box* dssd_frame_cutout_left_side_box = new G4Box("dssd_frame_cutout_left_side_box",
						     dssd_frame_cutout_left_side_half_X,
						     dssd_frame_cutout_left_side_half_Y,
						     dssd_frame_cutout_left_side_half_Z);

  G4LogicalVolume* dssd_frame_cutout_left_side_log = new G4LogicalVolume(dssd_frame_cutout_left_side_box,
									 C2H4_FR4,
									 "dssd_frame_cutout_left_side_log", 0, 0, 0);

  /*
  // Only for visualization
  G4VPhysicalVolume* dssd_frame_cutout_left_side = new G4PVPlacement(0, G4ThreeVector(0.5*dssd_cv_frame_b + dssd_cv_frame_gh1 + dssd_cv_frame_gh2 + 0.5*(dssd_frame_whole_outer_R - 0.5*dssd_cv_frame_b - dssd_cv_frame_gh1 - dssd_cv_frame_gh2),0,0),
								     dssd_frame_cutout_left_side_log,
								     "dssd_frame_cutout_left_side", room_log, false, 0);
  */

  ///////////////////////////////////////////
  // back side cutout (bottom side, as seen by incoming beam)
  ///////////////////////////////////////////
  
  G4double dssd_frame_cutout_bottom_side_half_X = 0.5*dssd_cv_frame_b;
  G4double dssd_frame_cutout_bottom_side_half_Y = 0.5*(dssd_frame_whole_outer_R - 0.5*dssd_cv_frame_f - dssd_cv_frame_gv);
  G4double dssd_frame_cutout_bottom_side_half_Z = dssd_frame_whole_half_length + 1.0*mm; // +1.0*mm to correctly visualize the subtraction

  G4Box* dssd_frame_cutout_bottom_side_box = new G4Box("dssd_frame_cutout_bottom_side_box",
						       dssd_frame_cutout_bottom_side_half_X,
						       dssd_frame_cutout_bottom_side_half_Y,
						       dssd_frame_cutout_bottom_side_half_Z);

  G4LogicalVolume* dssd_frame_cutout_bottom_side_log = new G4LogicalVolume(dssd_frame_cutout_bottom_side_box,
									   C2H4_FR4,
									   "dssd_frame_cutout_bottom_side_log", 0, 0, 0);

  /*
  // Only for visualization
  G4VPhysicalVolume* dssd_frame_cutout_bottom_side = new G4PVPlacement(0, G4ThreeVector(0,-0.5*dssd_cv_frame_f - dssd_cv_frame_gv - 0.5*(dssd_frame_whole_outer_R - 0.5*dssd_cv_frame_f - dssd_cv_frame_gv),0),
								       dssd_frame_cutout_bottom_side_log,
								       "dssd_frame_cutout_bottom_side", room_log, false, 0);
  */

  ///////////////////////////////////////////
  // silicon cutout
  ///////////////////////////////////////////

  G4double dssd_frame_cutout_silicon_half_X = 0.5*(21.8*mm); // the whole silicon chip, not just the active area
  G4double dssd_frame_cutout_silicon_half_Y = 0.5*(21.8*mm); // the whole silicon chip, not just the active area
  G4double dssd_frame_cutout_silicon_half_Z = dssd_frame_whole_half_length + 1.0*mm; // +1.0*mm to correctly visualize the subtraction
  
  G4Box* dssd_frame_cutout_silicon_box = new G4Box("dssd_frame_cutout_silicon_box",
						   dssd_frame_cutout_silicon_half_X,
						   dssd_frame_cutout_silicon_half_Y,
						   dssd_frame_cutout_silicon_half_Z);

  ///////////////////////////////////////////
  // additional front side cutout (left side, as seen by incoming beam)
  // this is on the front side of the detector where the ribbon cables attach to the silicon chip
  ///////////////////////////////////////////

  G4double dssd_frame_additional_cutout_left_side_half_X = 0.5*dssd_cv_frame_gh1;
  G4double dssd_frame_additional_cutout_left_side_half_Y = 0.5*dssd_cv_frame_f;
  G4double dssd_frame_additional_cutout_left_side_half_Z = 0.5*dssd_frame_whole_half_length + 0.1*mm; // thickness is estimated as a quarter of the thickness of the frame; +0.1*mm to correctly visualize the subtraction

  G4Box* dssd_frame_additional_cutout_left_side_box = new G4Box("dssd_frame_additional_cutout_left_side_box",
								dssd_frame_additional_cutout_left_side_half_X,
								dssd_frame_additional_cutout_left_side_half_Y,
								dssd_frame_additional_cutout_left_side_half_Z);
  
  G4LogicalVolume* dssd_frame_additional_cutout_left_side_log = new G4LogicalVolume(dssd_frame_additional_cutout_left_side_box,
										    C2H4_FR4,
										    "dssd_frame_additional_cutout_left_side_log", 0, 0, 0);
  
  // Only for visualization
  /*
  // Have to use 21.8*mm instead of dssd_cv_frame_b as I initially thought
  //G4VPhysicalVolume* dssd_frame_additional_cutout_left_side = new G4PVPlacement(0, G4ThreeVector(0.5*dssd_cv_frame_b + 0.5*dssd_cv_frame_gh1, 0, (-1.0)*dssd_frame_additional_cutout_left_side_half_Z),
  G4VPhysicalVolume* dssd_frame_additional_cutout_left_side = new G4PVPlacement(0, G4ThreeVector(0.5*(21.8*mm) + 0.5*dssd_cv_frame_gh1, 0, (-1.0)*dssd_frame_additional_cutout_left_side_half_Z),
										dssd_frame_additional_cutout_left_side_log,
										"dssd_frame_additional_cutout_left_side", room_log, false, 0);
  */

  ///////////////////////////////////////////
  // additional back side cutout (bottom side AND back side, as seen by incoming beam)
  // this is on the back side of the detector where the ribbon cables attach to the silicon chip
  ///////////////////////////////////////////
  
  G4double dssd_frame_additional_cutout_bottom_side_half_X = 0.5*dssd_cv_frame_b;
  G4double dssd_frame_additional_cutout_bottom_side_half_Y = 0.5*dssd_cv_frame_gv;
  G4double dssd_frame_additional_cutout_bottom_side_half_Z = 0.5*dssd_frame_whole_half_length + 0.1*mm; // thickness is estimated as a quarter of the thickness of the frame; +0.1*mm to correctly visualize the subtraction

  G4Box* dssd_frame_additional_cutout_bottom_side_box = new G4Box("dssd_frame_additional_cutout_bottom_side_box",
								  dssd_frame_additional_cutout_bottom_side_half_X,
								  dssd_frame_additional_cutout_bottom_side_half_Y,
								  dssd_frame_additional_cutout_bottom_side_half_Z);
  
  G4LogicalVolume* dssd_frame_additional_cutout_bottom_side_log = new G4LogicalVolume(dssd_frame_additional_cutout_bottom_side_box,
										      C2H4_FR4,
										      "dssd_frame_additional_cutout_bottom_side_log", 0, 0, 0);
  
  // Only for visualization
  /*
  // Have to use 21.8*mm instead of dssd_cv_frame_f as I initially thought
  //G4VPhysicalVolume* dssd_frame_additional_cutout_bottom_side = new G4PVPlacement(0, G4ThreeVector(0, -1.0*(0.5*(dssd_cv_frame_f) + 0.5*dssd_cv_frame_gv), dssd_frame_additional_cutout_bottom_side_half_Z),
  G4VPhysicalVolume* dssd_frame_additional_cutout_bottom_side = new G4PVPlacement(0, G4ThreeVector(0, -1.0*(0.5*(21.8*mm) + 0.5*dssd_cv_frame_gv), dssd_frame_additional_cutout_bottom_side_half_Z),
										dssd_frame_additional_cutout_bottom_side_log,
										"dssd_frame_additional_cutout_bottom_side", room_log, false, 0);
  */

  ///////////////////////////////////////////
  // vacuum hole
  ///////////////////////////////////////////
  
  G4double dssd_frame_vacuum_hole_inner_R = 0.0*mm;
  G4double dssd_frame_vacuum_hole_outer_R = dssd_cv_frame_rh;
  G4double dssd_frame_vacuum_hole_length = dssd_cv_frame_thickness;
  G4double dssd_frame_vacuum_hole_half_length = 0.5*dssd_frame_whole_length + 1.0*mm; // +1.0*mm to correctly visualize the subtraction
  G4double dssd_frame_vacuum_hole_start_angle = 0.0*deg;
  G4double dssd_frame_vacuum_hole_span_angle = 360.0*deg;

  G4Tubs* dssd_frame_vacuum_hole_tube = new G4Tubs("dssd_frame_vacuum_hole_tube",
						   dssd_frame_vacuum_hole_inner_R,
						   dssd_frame_vacuum_hole_outer_R,
						   dssd_frame_vacuum_hole_half_length,
						   dssd_frame_vacuum_hole_start_angle,
						   dssd_frame_vacuum_hole_span_angle);
  
  // Subtract volumes from the "whole" frame to create the final frame

  // Subtract the front side cutout (left side, as seen by incoming beam)
  // The placement of the subtraction should match the commented out physical volume for this portion (which is for visualization only)
  G4SubtractionSolid* dssd_frame_sub_1 = new G4SubtractionSolid("dssd_frame_sub_1",
								dssd_frame_whole_tube,
								dssd_frame_cutout_left_side_box,
								0,
								G4ThreeVector(0.5*dssd_cv_frame_b + dssd_cv_frame_gh1 + dssd_cv_frame_gh2 + 0.5*(dssd_frame_whole_outer_R - 0.5*dssd_cv_frame_b - dssd_cv_frame_gh1 - dssd_cv_frame_gh2),0,0));

  // Subtract the back side cutout (bottom side, as seen by incoming beam)
  // The placement of the subtraction should match the commented out physical volume for this portion (which is for visualization only)
  G4SubtractionSolid* dssd_frame_sub_2 = new G4SubtractionSolid("dssd_frame_sub_2",
								dssd_frame_sub_1,
								dssd_frame_cutout_bottom_side_box,
								0,
								G4ThreeVector(0,-0.5*dssd_cv_frame_f - dssd_cv_frame_gv - 0.5*(dssd_frame_whole_outer_R - 0.5*dssd_cv_frame_f - dssd_cv_frame_gv),0));
  
  // Subtract the three vacuum holes on the top side, as seen by incoming beam
  // top, middle
  G4SubtractionSolid* dssd_frame_sub_3 = new G4SubtractionSolid("dssd_frame_sub_3",
								dssd_frame_sub_2,
								dssd_frame_vacuum_hole_tube,
								0,
								G4ThreeVector(0,dssd_cv_frame_dch,0));
  // top, left
  G4SubtractionSolid* dssd_frame_sub_4 = new G4SubtractionSolid("dssd_frame_sub_4",
								dssd_frame_sub_3,
								dssd_frame_vacuum_hole_tube,
								0,
								G4ThreeVector(dssd_cv_frame_dbh, dssd_cv_frame_dch,0));
  // top, right
  G4SubtractionSolid* dssd_frame_sub_5 = new G4SubtractionSolid("dssd_frame_sub_5",
								dssd_frame_sub_4,
								dssd_frame_vacuum_hole_tube,
								0,
								G4ThreeVector(-1.0*dssd_cv_frame_dbh, dssd_cv_frame_dch,0));
  
  // Subtract three vacuum holes on the right side, as seen by incoming beam
  // right, middle
  G4SubtractionSolid* dssd_frame_sub_6 = new G4SubtractionSolid("dssd_frame_sub_6",
								dssd_frame_sub_5,
								dssd_frame_vacuum_hole_tube,
								0,
								G4ThreeVector(-1.0*dssd_cv_frame_dch,0,0));
  // right, top
  G4SubtractionSolid* dssd_frame_sub_7 = new G4SubtractionSolid("dssd_frame_sub_7",
								dssd_frame_sub_6,
								dssd_frame_vacuum_hole_tube,
								0,
								G4ThreeVector(-1.0*dssd_cv_frame_dch, dssd_cv_frame_dbh,0));
  // right, bottom
  G4SubtractionSolid* dssd_frame_sub_8 = new G4SubtractionSolid("dssd_frame_sub_8",
								dssd_frame_sub_7,
								dssd_frame_vacuum_hole_tube,
								0,
								G4ThreeVector(-1.0*dssd_cv_frame_dch, -1.0*dssd_cv_frame_dbh,0));

  // remove the space for the silicon chip
  G4SubtractionSolid* dssd_frame_sub_9 = new G4SubtractionSolid("dssd_frame_sub_9",
								dssd_frame_sub_8,
								dssd_frame_cutout_silicon_box,
								0,
								G4ThreeVector(0,0,0));

  // remove the additional material from the left, front side (as seen by the incoming beam)
  G4SubtractionSolid* dssd_frame_sub_10 = new G4SubtractionSolid("dssd_frame_sub_10",
								 dssd_frame_sub_9,
								 dssd_frame_additional_cutout_left_side_box,
								 0,
								 // Have to use 21.8*mm instead of dssd_cv_frame_b as I initially thought
								 G4ThreeVector(0.5*(21.8*mm) + 0.5*dssd_cv_frame_gh1, 0, (-1.0)*dssd_frame_additional_cutout_left_side_half_Z));

  // remove the additional material from the bottom, back side (as seen by the incoming beam)
  G4SubtractionSolid* dssd_frame_sub_11 = new G4SubtractionSolid("dssd_frame_sub_11",
								 dssd_frame_sub_10,
								 dssd_frame_additional_cutout_bottom_side_box,
								 0,
								 // Have to use 21.8*mm instead of dssd_cv_frame_b as I initially thought
								 G4ThreeVector(0, -1.0*(0.5*(21.8*mm) + 0.5*dssd_cv_frame_gv), dssd_frame_additional_cutout_bottom_side_half_Z));
  
  G4LogicalVolume* dssd_frame_log = new G4LogicalVolume(dssd_frame_sub_11,
							C2H4_FR4,
							"dssd_frame_log", 0, 0, 0);
  
  G4VPhysicalVolume* dssd_frame = new G4PVPlacement(0, G4ThreeVector(0,0,0),
						    dssd_frame_log,
						    "dssd_frame", room_log, false, 0);










  ///////////////////////////////////////////
  // The silicon chip of the DSSD
  // The non-active portion
  ///////////////////////////////////////////

  G4double dssd_silicon_chip_non_active_area_length = 1.03*mm;
  G4double dssd_silicon_chip_non_active_area_half_length = 0.5*dssd_silicon_chip_non_active_area_length;
  G4double dssd_silicon_chip_non_active_area_start_angle = 0.0*deg;
  G4double dssd_silicon_chip_non_active_area_span_angle = 360.0*deg;
  G4int dssd_silicon_chip_non_active_area_num_sides = 4;
  G4int dssd_silicon_chip_non_active_area_num_Z_planes = 2;
  G4double dssd_silicon_chip_non_active_area_z_plane[] = {(-1.0)*(dssd_silicon_chip_non_active_area_half_length), dssd_silicon_chip_non_active_area_half_length};
  G4double dssd_silicon_chip_non_active_area_inner_r[] = {0.5*(20.0*mm), 0.5*(20.0*mm)};
  G4double dssd_silicon_chip_non_active_area_outer_r[] = {0.5*(21.8*mm), 0.5*(21.8*mm)};
  
  G4Polyhedra* dssd_silicon_chip_non_active_area_poly = new G4Polyhedra("dssd_silicon_chip_non_active_area_poly",
									dssd_silicon_chip_non_active_area_start_angle,
									dssd_silicon_chip_non_active_area_span_angle,
									dssd_silicon_chip_non_active_area_num_sides,
									dssd_silicon_chip_non_active_area_num_Z_planes,
									dssd_silicon_chip_non_active_area_z_plane,
									dssd_silicon_chip_non_active_area_inner_r,
									dssd_silicon_chip_non_active_area_outer_r);
  
  G4LogicalVolume* dssd_silicon_chip_non_active_area_log = new G4LogicalVolume(dssd_silicon_chip_non_active_area_poly,
									       Si,
									       "dssd_silicon_chip_non_active_area_log", 0, 0, 0);

  G4RotationMatrix* rot_45 = new G4RotationMatrix();
  rot_45->rotateZ(45*deg);
  
  G4VPhysicalVolume* dssd_silicon_chip_non_active_area = new G4PVPlacement(rot_45, G4ThreeVector(0,0,0),
									   dssd_silicon_chip_non_active_area_log,
									   "dssd_silicon_chip_non_active_area", room_log, false, 0);
  
  ///////////////////////////////////////////
  // The silicon chip of the DSSD
  // The active-area portion
  ///////////////////////////////////////////

  G4double dssd_silicon_chip_active_area_half_X = 0.5*(20.0*mm);
  G4double dssd_silicon_chip_active_area_half_Y = 0.5*(20.0*mm);
  G4double dssd_silicon_chip_active_area_half_Z = 0.5*(1.03*mm);
  
  G4Box* dssd_silicon_chip_active_area_box = new G4Box("dssd_silicon_chip_active_area_box",
						       dssd_silicon_chip_active_area_half_X,
						       dssd_silicon_chip_active_area_half_Y,
						       dssd_silicon_chip_active_area_half_Z);
  
  G4LogicalVolume* dssd_silicon_chip_active_area_log = new G4LogicalVolume(dssd_silicon_chip_active_area_box,
									   Si,
									   "dssd_silicon_chip_active_area_log", 0, 0, 0);
  
  G4VPhysicalVolume* dssd_silicon_chip_active_area = new G4PVPlacement(0, G4ThreeVector(0,0,0),
								       dssd_silicon_chip_active_area_log,
								       "DSSD", room_log, false, 0);
  










  ///////////////////////////////////////////
  // The ribbon cables of the DSSD
  // (left side, as seen by the incoming beam)
  // (this is for the front side of the DSSD)
  ///////////////////////////////////////////

  G4double bending_offset = 1.50*mm;
  G4RotationMatrix* rot_left_side_45 = new G4RotationMatrix();
  rot_left_side_45->rotateY(-45*deg);

  G4double dssd_ribbon_cable_left_side_half_X = (1./50.)*mm;
  G4double dssd_ribbon_cable_left_side_half_Y = 0.5*dssd_cv_frame_f;
  G4double dssd_ribbon_cable_left_side_half_Z = 0.5*((300.0-23.0)*mm);
  
  G4Box* dssd_ribbon_cable_left_side_box = new G4Box("dssd_ribbon_cable_left_side_box",
						     dssd_ribbon_cable_left_side_half_X,
						     dssd_ribbon_cable_left_side_half_Y,
						     dssd_ribbon_cable_left_side_half_Z);

  G4LogicalVolume* dssd_ribbon_cable_left_side_log = new G4LogicalVolume(dssd_ribbon_cable_left_side_box,
									 C2H4_flexible,
									 "dssd_ribbon_cable_left_side_log", 0, 0, 0);
  
  G4VPhysicalVolume* dssd_ribbon_cable_left_side = new G4PVPlacement(0, G4ThreeVector(0.5*dssd_cv_frame_b + dssd_cv_frame_gh1 + dssd_cv_frame_gh2 + bending_offset, 0, dssd_ribbon_cable_left_side_half_Z - dssd_frame_whole_half_length + bending_offset),
								     dssd_ribbon_cable_left_side_log,
								     "dssd_ribbon_cable_left_side", room_log, false, 0);
  
  ///////////////////////////////////////////
  // The ribbon cables of the DSSD
  // (BENT part of the left side, as seen by the incoming beam)
  // (this is for the front side of the DSSD)
  ///////////////////////////////////////////
  
  G4double dssd_ribbon_cable_bent_left_side_half_X = (1./50.)*mm;
  G4double dssd_ribbon_cable_bent_left_side_half_Y = 0.5*dssd_cv_frame_f;
  G4double dssd_ribbon_cable_bent_left_side_half_Z = 0.5*sqrt(pow(bending_offset,2) + pow(bending_offset,2));
  
  G4Box* dssd_ribbon_cable_bent_left_side_box = new G4Box("dssd_ribbon_cable_bent_left_side_box",
							  dssd_ribbon_cable_bent_left_side_half_X,
							  dssd_ribbon_cable_bent_left_side_half_Y,
							  dssd_ribbon_cable_bent_left_side_half_Z);
  
  G4LogicalVolume* dssd_ribbon_cable_bent_left_side_log = new G4LogicalVolume(dssd_ribbon_cable_bent_left_side_box,
									      C2H4_flexible,
									      "dssd_ribbon_cable_bent_left_side_log", 0, 0, 0);
  
  G4VPhysicalVolume* dssd_ribbon_cable_bent_left_side = new G4PVPlacement(rot_left_side_45, G4ThreeVector(0.5*dssd_cv_frame_b + dssd_cv_frame_gh1 + dssd_cv_frame_gh2 + 0.5*bending_offset, 0, (-1.0)*dssd_frame_whole_half_length + 0.5*bending_offset),
									  dssd_ribbon_cable_bent_left_side_log,
									  "dssd_ribbon_cable_bent_left_side", room_log, false, 0);

  ///////////////////////////////////////////
  // The ribbon cables of the DSSD
  // (bottom side, as seen by the incoming beam)
  // (this is for the back side of the DSSD)
  ///////////////////////////////////////////

  G4RotationMatrix* rot_bottom_side_45 = new G4RotationMatrix();
  rot_bottom_side_45->rotateX(-45*deg);

  G4double dssd_ribbon_cable_bottom_side_half_X = 0.5*dssd_cv_frame_b;
  G4double dssd_ribbon_cable_bottom_side_half_Y = (1./50.)*mm;
  G4double dssd_ribbon_cable_bottom_side_half_Z = 0.5*((300.0-23.0)*mm);
  
  G4Box* dssd_ribbon_cable_bottom_side_box = new G4Box("dssd_ribbon_cable_bottom_side_box",
						       dssd_ribbon_cable_bottom_side_half_X,
						       dssd_ribbon_cable_bottom_side_half_Y,
						       dssd_ribbon_cable_bottom_side_half_Z);

  G4LogicalVolume* dssd_ribbon_cable_bottom_side_log = new G4LogicalVolume(dssd_ribbon_cable_bottom_side_box,
									   C2H4_flexible,
									   "dssd_ribbon_cable_bottom_side_log", 0, 0, 0);
  
  G4VPhysicalVolume* dssd_ribbon_cable_bottom_side = new G4PVPlacement(0, G4ThreeVector(0, -0.5*dssd_cv_frame_f - dssd_cv_frame_gv - bending_offset, dssd_ribbon_cable_bottom_side_half_Z + dssd_frame_whole_half_length + bending_offset),
								       dssd_ribbon_cable_bottom_side_log,
								       "dssd_ribbon_cable_bottom_side", room_log, false, 0);
  
  ///////////////////////////////////////////
  // The ribbon cables of the DSSD
  // (BENT part of the bottom side, as seen by the incoming beam)
  // (this is for the back side of the DSSD)
  ///////////////////////////////////////////
  
  G4double dssd_ribbon_cable_bent_bottom_side_half_X = 0.5*dssd_cv_frame_b;
  G4double dssd_ribbon_cable_bent_bottom_side_half_Y = (1./50.)*mm;
  G4double dssd_ribbon_cable_bent_bottom_side_half_Z = 0.5*sqrt(pow(bending_offset,2) + pow(bending_offset,2));
  
  G4Box* dssd_ribbon_cable_bent_bottom_side_box = new G4Box("dssd_ribbon_cable_bent_bottom_side_box",
							    dssd_ribbon_cable_bent_bottom_side_half_X,
							    dssd_ribbon_cable_bent_bottom_side_half_Y,
							    dssd_ribbon_cable_bent_bottom_side_half_Z);
  
  G4LogicalVolume* dssd_ribbon_cable_bent_bottom_side_log = new G4LogicalVolume(dssd_ribbon_cable_bent_bottom_side_box,
										C2H4_flexible,
										"dssd_ribbon_cable_bent_bottom_side_log", 0, 0, 0);
  
  G4VPhysicalVolume* dssd_ribbon_cable_bent_bottom_side = new G4PVPlacement(rot_bottom_side_45, G4ThreeVector(0, -0.5*dssd_cv_frame_f - dssd_cv_frame_gv - 0.5*bending_offset, dssd_frame_whole_half_length + 0.5*bending_offset),
									    dssd_ribbon_cable_bent_bottom_side_log,
									    "dssd_ribbon_cable_bent_bottom_side", room_log, false, 0);












 





  ///////////////////////////////////////////
  // whole nylon veto detector holder
  // This will have portions removed from it to create the final holder
  ///////////////////////////////////////////

  G4double veto_holder_whole_inner_R = 0.0*mm;
  G4double veto_holder_whole_outer_R = ((3.9116)/2)*cm;
  G4double veto_holder_whole_length = 1.5875*cm;
  G4double veto_holder_whole_half_length = 0.5*veto_holder_whole_length;
  G4double veto_holder_whole_start_angle = 0.0*deg;
  G4double veto_holder_whole_span_angle = 360.0*deg;
  
  G4Tubs* veto_holder_whole_tube = new G4Tubs("veto_holder_whole_tube",
					      veto_holder_whole_inner_R,
					      veto_holder_whole_outer_R,
					      veto_holder_whole_half_length,
					      veto_holder_whole_start_angle,
					      veto_holder_whole_span_angle);
  
  G4LogicalVolume* veto_holder_whole_log = new G4LogicalVolume(veto_holder_whole_tube,
							       C12H22N2O2,
							       "veto_holder_whole_log", 0, 0, 0);

  // Distance between the center of the DSSD frame and the veto holder
  G4double db_dssd_and_veto_holder = 2.0*cm;
  
  // Only for visualization
  /*
  G4VPhysicalVolume* veto_holder_whole = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder),
							   veto_holder_whole_log,
							   "veto_holder_whole", room_log, false, 0);
  */

  ///////////////////////////////////////////
  // cutout material for the actual silicon veto detector
  ///////////////////////////////////////////

  G4double veto_holder_silicon_cutout_inner_R = 0.0*mm;
  G4double veto_holder_silicon_cutout_outer_R = ((2.8702)/2)*cm;
  G4double veto_holder_silicon_cutout_length = 1.22936*cm;
  G4double veto_holder_silicon_cutout_half_length = 0.5*veto_holder_silicon_cutout_length + 0.1*mm; // +0.1 mm to correctly visualize the subtraction
  G4double veto_holder_silicon_cutout_start_angle = 0.0*deg;
  G4double veto_holder_silicon_cutout_span_angle = 360.0*deg;
  
  G4Tubs* veto_holder_silicon_cutout_tube = new G4Tubs("veto_holder_silicon_cutout_tube",
						       veto_holder_silicon_cutout_inner_R,
						       veto_holder_silicon_cutout_outer_R,
						       veto_holder_silicon_cutout_half_length,
						       veto_holder_silicon_cutout_start_angle,
						       veto_holder_silicon_cutout_span_angle);
  
  G4LogicalVolume* veto_holder_silicon_cutout_log = new G4LogicalVolume(veto_holder_silicon_cutout_tube,
									C12H22N2O2,
									"veto_holder_silicon_cutout_log", 0, 0, 0);

  // Only for visualization
  // -0.1*mm in order to undo the +0.1*mm above (that was done to correctly visualize the subtraction)
  /*
  G4VPhysicalVolume* veto_holder_silicon_cutout = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length + (veto_holder_silicon_cutout_half_length - 0.1*mm)),
								    veto_holder_silicon_cutout_log,
								    "veto_holder_silicon_cutout", room_log, false, 0);
  */

  ///////////////////////////////////////////
  // cutout material in the back for the microdot to lemo cable
  ///////////////////////////////////////////

  G4double veto_holder_back_cutout_inner_R = 0.0*mm;
  G4double veto_holder_back_cutout_outer_R = ((1.905)/2)*cm;
  G4double veto_holder_back_cutout_length = veto_holder_whole_length - veto_holder_silicon_cutout_length;
  G4double veto_holder_back_cutout_half_length = 0.5*veto_holder_back_cutout_length + 0.1*mm; // +0.1 mm to correctly visualize the subtraction
  G4double veto_holder_back_cutout_start_angle = 0.0*deg;
  G4double veto_holder_back_cutout_span_angle = 360.0*deg;
  
  G4Tubs* veto_holder_back_cutout_tube = new G4Tubs("veto_holder_back_cutout_tube",
						    veto_holder_back_cutout_inner_R,
						    veto_holder_back_cutout_outer_R,
						    veto_holder_back_cutout_half_length,
						    veto_holder_back_cutout_start_angle,
						    veto_holder_back_cutout_span_angle);

  G4LogicalVolume* veto_holder_back_cutout_log = new G4LogicalVolume(veto_holder_back_cutout_tube,
								     C12H22N2O2,
								     "veto_holder_back_cutout_log", 0, 0, 0);

  // Only for visualization
  // -0.1*mm in order to undo the +0.1*mm above (that was done to correctly visualize the subtraction)
  /*
  G4VPhysicalVolume* veto_holder_back_cutout = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder + veto_holder_whole_half_length - (veto_holder_back_cutout_half_length - 0.1*mm)),
								 veto_holder_back_cutout_log,
								 "veto_holder_back_cutout", room_log, false, 0);
  */

  ///////////////////////////////////////////
  // left side cutout, as seen by incoming beam
  ///////////////////////////////////////////
  
  //G4double veto_holder_cutout_left_side_half_X = 0.5*((2.8702/2)*cm - 1.27*cm);
  G4double veto_holder_cutout_left_side_half_X = 1.0*cm; // arbitrary, provided this is longer than a certain length and the placement is correct
  //G4double veto_holder_cutout_left_side_half_Y = 1.27*cm;
  G4double veto_holder_cutout_left_side_half_Y = 2*cm; // arbitrary, provided this is longer than a certain length and the placement is correct
  G4double veto_holder_cutout_left_side_half_Z = veto_holder_whole_half_length + 1.0*mm; // +1.0*mm to correctly visualize the subtraction
  
  G4Box* veto_holder_cutout_left_side_box = new G4Box("veto_holder_cutout_left_side_box",
						      veto_holder_cutout_left_side_half_X,
						      veto_holder_cutout_left_side_half_Y,
						      veto_holder_cutout_left_side_half_Z);
  
  G4LogicalVolume* veto_holder_cutout_left_side_log = new G4LogicalVolume(veto_holder_cutout_left_side_box,
									  C12H22N2O2,
									  "veto_holder_cutout_left_side_log", 0, 0, 0);
  
  // Only for visualization
  /*
  G4VPhysicalVolume* veto_holder_cutout_left_side = new G4PVPlacement(0, G4ThreeVector(1.27*cm + veto_holder_cutout_left_side_half_X, 0, db_dssd_and_veto_holder),
								      veto_holder_cutout_left_side_log,
								      "veto_holder_cutout_left_side", room_log, false, 0);
  */

  ///////////////////////////////////////////
  // bottom side cutout, as seen by incoming beam
  ///////////////////////////////////////////
  
  //G4double veto_holder_cutout_bottom_side_half_X = 1.27*cm;
  G4double veto_holder_cutout_bottom_side_half_X = 2*cm; // arbitrary, provided this is longer than a certain length and the placement is correct
  //G4double veto_holder_cutout_bottom_side_half_Y = 0.5*((2.8702/2)*cm - 1.27*cm);
  G4double veto_holder_cutout_bottom_side_half_Y = 1.0*cm; // arbitrary, provided this is longer than a certain length and the placement is correct
  G4double veto_holder_cutout_bottom_side_half_Z = veto_holder_whole_half_length + 1.0*mm; // +1.0*mm to correctly visualize the subtraction
  
  G4Box* veto_holder_cutout_bottom_side_box = new G4Box("veto_holder_cutout_bottom_side_box",
							veto_holder_cutout_bottom_side_half_X,
							veto_holder_cutout_bottom_side_half_Y,
							veto_holder_cutout_bottom_side_half_Z);
  
  G4LogicalVolume* veto_holder_cutout_bottom_side_log = new G4LogicalVolume(veto_holder_cutout_bottom_side_box,
									    C12H22N2O2,
									    "veto_holder_cutout_bottom_side_log", 0, 0, 0);
  
  // Only for visualization
  /*
  G4VPhysicalVolume* veto_holder_cutout_bottom_side = new G4PVPlacement(0, G4ThreeVector(0, (-1.0)*(1.27*cm + veto_holder_cutout_bottom_side_half_Y), db_dssd_and_veto_holder),
									veto_holder_cutout_bottom_side_log,
									"veto_holder_cutout_bottom_side", room_log, false, 0);
  */

  ///////////////////////////////////////////
  // holes for the threaded rods
  ///////////////////////////////////////////

  G4double veto_holder_rod_hole_inner_R = 0.0*mm;
  G4double veto_holder_rod_hole_outer_R = dssd_cv_frame_rh;
  G4double veto_holder_rod_hole_length = veto_holder_whole_length;
  G4double veto_holder_rod_hole_half_length = 0.5*veto_holder_rod_hole_length + 1.0*mm; // +1.0*mm to correctly visualize the subtraction
  G4double veto_holder_rod_hole_start_angle = 0.0*deg;
  G4double veto_holder_rod_hole_span_angle = 360.0*deg;

  G4Tubs* veto_holder_rod_hole_tube = new G4Tubs("veto_holder_rod_hole_tube",
						 veto_holder_rod_hole_inner_R,
						 veto_holder_rod_hole_outer_R,
						 veto_holder_rod_hole_half_length,
						 veto_holder_rod_hole_start_angle,
						 veto_holder_rod_hole_span_angle);
  
  // Subtract volumes from the "whole" holder to create the final holder

  ///////////////////////////////////////////
  // cutout material for the actual silicon veto detector
  // The placement of the subtraction should match the commented out physical volume for this portion (which is for visualization only), except
  // there is no db_dssd_and_veto_holder because the subtraction is with respect to the first volume
  ///////////////////////////////////////////
  G4SubtractionSolid* veto_holder_sub_1 = new G4SubtractionSolid("veto_holder_sub_1",
								 veto_holder_whole_tube,
								 veto_holder_silicon_cutout_tube,
								 0,
								 G4ThreeVector(0, 0, (-1.0)*veto_holder_whole_half_length + (veto_holder_silicon_cutout_half_length - 0.1*mm)));

  ///////////////////////////////////////////
  // cutout material in the back for the microdot to lemo cable
  // The placement of the subtraction should match the commented out physical volume for this portion (which is for visualization only), except
  // there is no db_dssd_and_veto_holder because the subtraction is with respect to the first volume
  ///////////////////////////////////////////
  G4SubtractionSolid* veto_holder_sub_2 = new G4SubtractionSolid("veto_holder_sub_2",
								 veto_holder_sub_1,
								 veto_holder_back_cutout_tube,
								 0,
								 G4ThreeVector(0, 0, veto_holder_whole_half_length - (veto_holder_back_cutout_half_length - 0.1*mm)));

  ///////////////////////////////////////////
  // left side cutout, as seen by incoming beam
  // The placement of the subtraction should match the commented out physical volume for this portion (which is for visualization only), except
  // there is no db_dssd_and_veto_holder because the subtraction is with respect to the first volume
  ///////////////////////////////////////////
  G4SubtractionSolid* veto_holder_sub_3 = new G4SubtractionSolid("veto_holder_sub_3",
								 veto_holder_sub_2,
								 veto_holder_cutout_left_side_box,
								 0,
								 G4ThreeVector(1.27*cm + veto_holder_cutout_left_side_half_X, 0, 0));

  ///////////////////////////////////////////
  // bottom side cutout, as seen by incoming beam
  // The placement of the subtraction should match the commented out physical volume for this portion (which is for visualization only), except
  // there is no db_dssd_and_veto_holder because the subtraction is with respect to the first volume
  ///////////////////////////////////////////
  G4SubtractionSolid* veto_holder_sub_4 = new G4SubtractionSolid("veto_holder_sub_4",
								 veto_holder_sub_3,
								 veto_holder_cutout_bottom_side_box,
								 0,
								 G4ThreeVector(0, (-1.0)*(1.27*cm + veto_holder_cutout_bottom_side_half_Y), 0));

  ///////////////////////////////////////////
  // top hole cutout for threaded rod, as seen by incoming beam
  ///////////////////////////////////////////
  G4SubtractionSolid* veto_holder_sub_5 = new G4SubtractionSolid("veto_holder_sub_5",
								 veto_holder_sub_4,
								 veto_holder_rod_hole_tube,
								 0,
								 G4ThreeVector(0,dssd_cv_frame_dch,0)); // using the same location as hole on the DSSD frame, even though this is slightly different than 7mdx85_008

  ///////////////////////////////////////////
  // right hole cutout for threaded rod, as seen by incoming beam
  ///////////////////////////////////////////
  G4SubtractionSolid* veto_holder_sub_6 = new G4SubtractionSolid("veto_holder_sub_6",
								 veto_holder_sub_5,
								 veto_holder_rod_hole_tube,
								 0,
								 G4ThreeVector(-1.0*dssd_cv_frame_dch,0,0)); // using the same location as hole on the DSSD frame, even though this is slightly different than 7mdx85_008
  
  G4LogicalVolume* veto_holder_log = new G4LogicalVolume(veto_holder_sub_6,
							 C12H22N2O2,
							 "veto_holder_log", 0, 0, 0);

  G4VPhysicalVolume* veto_holder = new G4PVPlacement(0, G4ThreeVector(0,0,db_dssd_and_veto_holder),
						     veto_holder_log,
						     "veto_holder", room_log, false, 0);

  
  
  










  







  ///////////////////////////////////////////
  // veto housing
  // This will have portions removed from it to create the final detector
  ///////////////////////////////////////////

  G4double veto_whole_inner_R = 0.0*mm;
  G4double veto_whole_outer_R = ((28.6)/2)*mm;
  G4double veto_whole_length = 12.3*mm;
  G4double veto_whole_half_length = 0.5*veto_whole_length;
  G4double veto_whole_start_angle = 0.0*deg;
  G4double veto_whole_span_angle = 360.0*deg;
  
  G4Tubs* veto_whole_tube = new G4Tubs("veto_whole_tube",
				       veto_whole_inner_R,
				       veto_whole_outer_R,
				       veto_whole_half_length,
				       veto_whole_start_angle,
				       veto_whole_span_angle);
  
  G4LogicalVolume* veto_whole_log = new G4LogicalVolume(veto_whole_tube,
							Cu3Zn2,
							"veto_whole_log", 0, 0, 0);

  // Only for visualization
  /*
  G4VPhysicalVolume* veto_whole = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length + veto_whole_half_length),
						    veto_whole_log,
						    "veto_whole", room_log, false, 0);
  */

  ///////////////////////////////////////////
  // interior cutout of the veto detector brass housing
  ///////////////////////////////////////////

  G4double veto_interior_cutout_inner_R = 0.0*mm;
  G4double veto_interior_cutout_outer_R = ((28.6)/2)*mm - 1.5*mm; // this is a radius, so only subtract once width of the wall
  G4double veto_interior_cutout_length = 12.3*mm - 2*(1.5*mm); // this is a length, so subtract twice the width of the wall
  G4double veto_interior_cutout_half_length = 0.5*veto_interior_cutout_length; // no small add on here because you will never visualize this subtraction
  G4double veto_interior_cutout_start_angle = 0.0*deg;
  G4double veto_interior_cutout_span_angle = 360.0*deg;
  
  G4Tubs* veto_interior_cutout_tube = new G4Tubs("veto_interior_cutout_tube",
						 veto_interior_cutout_inner_R,
						 veto_interior_cutout_outer_R,
						 veto_interior_cutout_half_length,
						 veto_interior_cutout_start_angle,
						 veto_interior_cutout_span_angle);
  
  G4LogicalVolume* veto_interior_cutout_log = new G4LogicalVolume(veto_interior_cutout_tube,
								  Cu3Zn2,
								  "veto_interior_cutout_log", 0, 0, 0);

  // Only for visualization
  /*
  G4VPhysicalVolume* veto_interior_cutout = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length + veto_whole_half_length), // + veto_whole_half_length, not + veto_interior_cutout_half_length
							      veto_interior_cutout_log,
							      "veto_interior_cutout", room_log, false, 0);
  */
  
  ///////////////////////////////////////////
  // cutout for the window of the silicon wafer
  ///////////////////////////////////////////

  G4double veto_window_cutout_inner_R = 0.0*mm;
  G4double veto_window_cutout_outer_R = ((19.5)/2)*mm;
  G4double veto_window_cutout_length = 1.5*mm;
  G4double veto_window_cutout_half_length = 0.5*veto_window_cutout_length + 0.1*mm; // +0.1*mm to correctly visualize the subtraction
  G4double veto_window_cutout_start_angle = 0.0*deg;
  G4double veto_window_cutout_span_angle = 360.0*deg;
  
  G4Tubs* veto_window_cutout_tube = new G4Tubs("veto_window_cutout_tube",
					       veto_window_cutout_inner_R,
					       veto_window_cutout_outer_R,
					       veto_window_cutout_half_length,
					       veto_window_cutout_start_angle,
					       veto_window_cutout_span_angle);
  
  G4LogicalVolume* veto_window_cutout_log = new G4LogicalVolume(veto_window_cutout_tube,
								Cu3Zn2,
								"veto_window_cutout_log", 0, 0, 0);
  
  // Only for visualization
  // -0.1*mm in order to undo the +0.1*mm above (that was done to correctly visualize the subtraction)
  /*
  G4VPhysicalVolume* veto_window_cutout = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length + (veto_window_cutout_half_length - 0.1*mm)),
							    veto_window_cutout_log,
							    "veto_window_cutout", room_log, false, 0);
  */


  ///////////////////////////////////////////
  // Subtract the volumes to make the final volume
  ///////////////////////////////////////////
  
  // Remember, subtraction is with respect to the first volume, so the placement here does not match the commented out physical placement above (which is only for visualization)
  G4SubtractionSolid* veto_housing_sub_1 = new G4SubtractionSolid("veto_housing_sub_1",
								  veto_whole_tube,
								  veto_interior_cutout_tube,
								  0,
								  G4ThreeVector(0, 0, 0));

  // Remember, subtraction is with respect to the first volume, so the placement here does not match the commented out physical placement above (which is only for visualization)
  G4SubtractionSolid* veto_housing_sub_2 = new G4SubtractionSolid("veto_housing_sub_2",
								  veto_housing_sub_1,
								  veto_window_cutout_tube,
								  0,
								  G4ThreeVector(0, 0, (-1.0)*veto_whole_half_length + (veto_window_cutout_half_length - 0.1*mm)));

  G4LogicalVolume* veto_housing_log = new G4LogicalVolume(veto_housing_sub_2,
							  Cu3Zn2,
							  "veto_housing_log", 0, 0, 0);

  // The placement should match the commented out placement of veto_whole_tube (which is for visualization only)
  G4VPhysicalVolume* veto_housing = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length + veto_whole_half_length),
						      veto_housing_log,
						      "veto_housing", room_log, false, 0);

  



  ///////////////////////////////////////////
  // silicon wafer in the veto detector
  ///////////////////////////////////////////

  G4double veto_silicon_wafer_inner_R = 0.0*mm;
  G4double veto_silicon_wafer_outer_R = ((19.5)/2)*mm;
  G4double veto_silicon_wafer_length = 0.5*mm;
  G4double veto_silicon_wafer_half_length = 0.5*veto_silicon_wafer_length;
  G4double veto_silicon_wafer_start_angle = 0.0*deg;
  G4double veto_silicon_wafer_span_angle = 360.0*deg;
  
  G4Tubs* veto_silicon_wafer_tube = new G4Tubs("veto_silicon_wafer_tube",
					       veto_silicon_wafer_inner_R,
					       veto_silicon_wafer_outer_R,
					       veto_silicon_wafer_half_length,
					       veto_silicon_wafer_start_angle,
					       veto_silicon_wafer_span_angle);

  G4LogicalVolume* veto_silicon_wafer_log = new G4LogicalVolume(veto_silicon_wafer_tube,
								Si,
								"veto_silicon_wafer_log", 0, 0, 0);

  G4VPhysicalVolume* veto_silicon_wafer = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length + veto_window_cutout_length + veto_silicon_wafer_half_length),
							    veto_silicon_wafer_log,
							    "Veto", room_log, false, 0);
  


  
  ///////////////////////////////////////////
  // Rear microdot mount on the veto detector
  // The hex nut on the back of the veto detector
  ///////////////////////////////////////////

  G4double veto_housing_hex_nut_length = ((7.1)/2)*mm;
  G4double veto_housing_hex_nut_half_length = 0.5*veto_housing_hex_nut_length;
  G4double veto_housing_hex_nut_start_angle = 0.0*deg;
  G4double veto_housing_hex_nut_span_angle = 360.0*deg;
  G4int veto_housing_hex_nut_num_sides = 6;
  G4int veto_housing_hex_nut_num_Z_planes = 2;
  G4double veto_housing_hex_nut_z_plane[] = {(-1.0)*(veto_housing_hex_nut_half_length), veto_housing_hex_nut_half_length};
  G4double veto_housing_hex_nut_inner_r[] = {2*mm, 2*mm}; // same as veto_housing_cable_connection_inner_R
  G4double veto_housing_hex_nut_outer_r[] = {4*mm, 4*mm};
  
  G4Polyhedra* veto_housing_hex_nut_poly = new G4Polyhedra("veto_housing_hex_nut_poly",
							   veto_housing_hex_nut_start_angle,
							   veto_housing_hex_nut_span_angle,
							   veto_housing_hex_nut_num_sides,
							   veto_housing_hex_nut_num_Z_planes,
							   veto_housing_hex_nut_z_plane,
							   veto_housing_hex_nut_inner_r,
							   veto_housing_hex_nut_outer_r);
  
  G4LogicalVolume* veto_housing_hex_nut_log = new G4LogicalVolume(veto_housing_hex_nut_poly,
								  Cr20Ni8Fe76,
								  "veto_housing_hex_nut_log", 0, 0, 0);

  G4VPhysicalVolume* veto_housing_hex_nut = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length + veto_whole_length + veto_housing_hex_nut_half_length),
							      veto_housing_hex_nut_log,
							      "veto_housing_hex_nut", room_log, false, 0);
  
  ///////////////////////////////////////////
  // Rear microdot mount on the veto detector
  // The part that connects to the microdot end of the lemo to microdot cable
  ///////////////////////////////////////////

  G4double veto_housing_cable_connection_inner_R = 2.0*mm;
  G4double veto_housing_cable_connection_outer_R = 3.0*mm;
  G4double veto_housing_cable_connection_length = ((7.1)/2)*mm;
  G4double veto_housing_cable_connection_half_length = 0.5*veto_housing_cable_connection_length;
  G4double veto_housing_cable_connection_start_angle = 0.0*deg;
  G4double veto_housing_cable_connection_span_angle = 360.0*deg;
  
  G4Tubs* veto_housing_cable_connection_tube = new G4Tubs("veto_housing_cable_connection_tube",
							  veto_housing_cable_connection_inner_R,
							  veto_housing_cable_connection_outer_R,
							  veto_housing_cable_connection_half_length,
							  veto_housing_cable_connection_start_angle,
							  veto_housing_cable_connection_span_angle);
  
  G4LogicalVolume* veto_housing_cable_connection_log = new G4LogicalVolume(veto_housing_cable_connection_tube,
									   Cr20Ni8Fe76,
									   "veto_housing_cable_connection_log", 0, 0, 0);
  
  G4VPhysicalVolume* veto_housing_cable_connection = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length + veto_whole_length + veto_housing_hex_nut_length + veto_housing_cable_connection_half_length),
								       veto_housing_cable_connection_log,
								       "veto_housing_cable_connection", room_log, false, 0);
  



























  
  ///////////////////////////////////////////
  // microdot to lemo cable
  // part A
  // (The end that attaches directory to the back of the veto detector)
  ///////////////////////////////////////////
  
  G4double veto_cable_connection_gap = 1.0*mm; // the gap where the microdot to lemo cable attaches to the back of the veto detector
  
  G4double microdot_to_lemo_cable_A_inner_R = 3.0*mm;
  G4double microdot_to_lemo_cable_A_outer_R = 4.0*mm;
  G4double microdot_to_lemo_cable_A_length = 6.0*mm;
  G4double microdot_to_lemo_cable_A_half_length = 0.5*microdot_to_lemo_cable_A_length;
  G4double microdot_to_lemo_cable_A_start_angle = 0.0*deg;
  G4double microdot_to_lemo_cable_A_span_angle = 360.0*deg;
  
  G4Tubs* microdot_to_lemo_cable_A_tube = new G4Tubs("microdot_to_lemo_cable_A_tube",
						     microdot_to_lemo_cable_A_inner_R,
						     microdot_to_lemo_cable_A_outer_R,
						     microdot_to_lemo_cable_A_half_length,
						     microdot_to_lemo_cable_A_start_angle,
						     microdot_to_lemo_cable_A_span_angle);
  
  G4LogicalVolume* microdot_to_lemo_cable_A_log = new G4LogicalVolume(microdot_to_lemo_cable_A_tube,
								      Cu3Zn2,
								      "microdot_to_lemo_cable_A_log", 0, 0, 0);

  // Do not need veto_housing_cable_connection_length in the placement
  G4VPhysicalVolume* microdot_to_lemo_cable_A = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length + veto_whole_length + veto_housing_hex_nut_length + veto_cable_connection_gap + microdot_to_lemo_cable_A_half_length),
								  microdot_to_lemo_cable_A_log,
								  "microdot_to_lemo_cable_A", room_log, false, 0);
  
  ///////////////////////////////////////////
  // microdot to lemo cable
  // part B
  // (The part attached to part A)
  ///////////////////////////////////////////
  
  G4double microdot_to_lemo_cable_B_inner_R = 2.5*mm;
  G4double microdot_to_lemo_cable_B_outer_R = 3.0*mm;
  G4double microdot_to_lemo_cable_B_length = 4.0*mm;
  G4double microdot_to_lemo_cable_B_half_length = 0.5*microdot_to_lemo_cable_B_length;
  G4double microdot_to_lemo_cable_B_start_angle = 0.0*deg;
  G4double microdot_to_lemo_cable_B_span_angle = 360.0*deg;
  
  G4Tubs* microdot_to_lemo_cable_B_tube = new G4Tubs("microdot_to_lemo_cable_B_tube",
						     microdot_to_lemo_cable_B_inner_R,
						     microdot_to_lemo_cable_B_outer_R,
						     microdot_to_lemo_cable_B_half_length,
						     microdot_to_lemo_cable_B_start_angle,
						     microdot_to_lemo_cable_B_span_angle);
  
  G4LogicalVolume* microdot_to_lemo_cable_B_log = new G4LogicalVolume(microdot_to_lemo_cable_B_tube,
								      Cu3Zn2,
								      "microdot_to_lemo_cable_B_log", 0, 0, 0);

  // Do not need veto_housing_cable_connection_length in the placement
  G4VPhysicalVolume* microdot_to_lemo_cable_B = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length + veto_whole_length + veto_housing_hex_nut_length + veto_cable_connection_gap + microdot_to_lemo_cable_A_length + microdot_to_lemo_cable_B_half_length),
								  microdot_to_lemo_cable_B_log,
								  "microdot_to_lemo_cable_B", room_log, false, 0);

  ///////////////////////////////////////////
  // microdot to lemo cable
  // part C
  // (The part attached to part B)
  ///////////////////////////////////////////
  
  G4double microdot_to_lemo_cable_C_inner_R = 2.0*mm;
  G4double microdot_to_lemo_cable_C_outer_R = 2.5*mm;
  G4double microdot_to_lemo_cable_C_length = 5.0*mm;
  G4double microdot_to_lemo_cable_C_half_length = 0.5*microdot_to_lemo_cable_C_length;
  G4double microdot_to_lemo_cable_C_start_angle = 0.0*deg;
  G4double microdot_to_lemo_cable_C_span_angle = 360.0*deg;
  
  G4Tubs* microdot_to_lemo_cable_C_tube = new G4Tubs("microdot_to_lemo_cable_C_tube",
						     microdot_to_lemo_cable_C_inner_R,
						     microdot_to_lemo_cable_C_outer_R,
						     microdot_to_lemo_cable_C_half_length,
						     microdot_to_lemo_cable_C_start_angle,
						     microdot_to_lemo_cable_C_span_angle);
  
  G4LogicalVolume* microdot_to_lemo_cable_C_log = new G4LogicalVolume(microdot_to_lemo_cable_C_tube,
								      Cu3Zn2,
								      "microdot_to_lemo_cable_C_log", 0, 0, 0);

  // Do not need veto_housing_cable_connection_length in the placement
  G4VPhysicalVolume* microdot_to_lemo_cable_C = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length + veto_whole_length + veto_housing_hex_nut_length + veto_cable_connection_gap + microdot_to_lemo_cable_A_length + microdot_to_lemo_cable_B_length + microdot_to_lemo_cable_C_half_length),
								  microdot_to_lemo_cable_C_log,
								  "microdot_to_lemo_cable_C", room_log, false, 0);
  
  ///////////////////////////////////////////
  // microdot to lemo cable
  // part D
  // (The part attached to part C)
  ///////////////////////////////////////////
  
  G4double microdot_to_lemo_cable_D_length = 3*mm;
  G4double microdot_to_lemo_cable_D_half_length = 0.5*microdot_to_lemo_cable_D_length;
  G4double microdot_to_lemo_cable_D_start_angle = 0.0*deg;
  G4double microdot_to_lemo_cable_D_span_angle = 360.0*deg;
  G4int microdot_to_lemo_cable_D_num_sides = 6;
  G4int microdot_to_lemo_cable_D_num_Z_planes = 2;
  G4double microdot_to_lemo_cable_D_z_plane[] = {(-1.0)*(microdot_to_lemo_cable_D_half_length), microdot_to_lemo_cable_D_half_length};
  G4double microdot_to_lemo_cable_D_inner_r[] = {1.75*mm, 1.75*mm};
  G4double microdot_to_lemo_cable_D_outer_r[] = {2.00*mm, 2.00*mm};
  
  G4Polyhedra* microdot_to_lemo_cable_D_poly = new G4Polyhedra("microdot_to_lemo_cable_D_poly",
							       microdot_to_lemo_cable_D_start_angle,
							       microdot_to_lemo_cable_D_span_angle,
							       microdot_to_lemo_cable_D_num_sides,
							       microdot_to_lemo_cable_D_num_Z_planes,
							       microdot_to_lemo_cable_D_z_plane,
							       microdot_to_lemo_cable_D_inner_r,
							       microdot_to_lemo_cable_D_outer_r);

  G4LogicalVolume* microdot_to_lemo_cable_D_log = new G4LogicalVolume(microdot_to_lemo_cable_D_poly,
								      Cu3Zn2,
								      "microdot_to_lemo_cable_D_log", 0, 0, 0);
  
  // Do not need veto_housing_cable_connection_length in the placement
  G4VPhysicalVolume* microdot_to_lemo_cable_D = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length + veto_whole_length + veto_housing_hex_nut_length + veto_cable_connection_gap + microdot_to_lemo_cable_A_length + microdot_to_lemo_cable_B_length + microdot_to_lemo_cable_C_length + microdot_to_lemo_cable_D_half_length),
								  microdot_to_lemo_cable_D_log,
								  "microdot_to_lemo_cable_D", room_log, false, 0);

  ///////////////////////////////////////////
  // microdot to lemo cable
  // part E
  // (The wire, which is attached to part D)
  ///////////////////////////////////////////
  
  G4double microdot_to_lemo_cable_E_inner_R = 0.0*mm;
  G4double microdot_to_lemo_cable_E_outer_R = 1.0*mm;
  G4double microdot_to_lemo_cable_E_length = 18.8*cm;
  G4double microdot_to_lemo_cable_E_half_length = 0.5*microdot_to_lemo_cable_E_length;
  G4double microdot_to_lemo_cable_E_start_angle = 0.0*deg;
  G4double microdot_to_lemo_cable_E_span_angle = 360.0*deg;
  
  G4Tubs* microdot_to_lemo_cable_E_tube = new G4Tubs("microdot_to_lemo_cable_E_tube",
						     microdot_to_lemo_cable_E_inner_R,
						     microdot_to_lemo_cable_E_outer_R,
						     microdot_to_lemo_cable_E_half_length,
						     microdot_to_lemo_cable_E_start_angle,
						     microdot_to_lemo_cable_E_span_angle);
  
  G4LogicalVolume* microdot_to_lemo_cable_E_log = new G4LogicalVolume(microdot_to_lemo_cable_E_tube,
								      Cu,
								      "microdot_to_lemo_cable_E_log", 0, 0, 0);
  
  // Do not need veto_housing_cable_connection_length in the placement
  G4VPhysicalVolume* microdot_to_lemo_cable_E = new G4PVPlacement(0, G4ThreeVector(0, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length + veto_whole_length + veto_housing_hex_nut_length + veto_cable_connection_gap + microdot_to_lemo_cable_A_length + microdot_to_lemo_cable_B_length + microdot_to_lemo_cable_C_length + microdot_to_lemo_cable_D_length + microdot_to_lemo_cable_E_half_length),
								  microdot_to_lemo_cable_E_log,
								  "microdot_to_lemo_cable_E", room_log, false, 0);

















  



  ///////////////////////////////////////////
  // stainless steel threaded rods that support the DSSD and veto detector
  ///////////////////////////////////////////
  
  G4double stainless_steel_threaded_rod_inner_R = 0.0*mm;
  G4double stainless_steel_threaded_rod_outer_R = dssd_cv_frame_rh;
  G4double stainless_steel_threaded_rod_length = 30*cm;
  G4double stainless_steel_threaded_rod_half_length = 0.5*stainless_steel_threaded_rod_length;
  G4double stainless_steel_threaded_rod_start_angle = 0.0*deg;
  G4double stainless_steel_threaded_rod_span_angle = 360.0*deg;
  
  G4Tubs* stainless_steel_threaded_rod_tube = new G4Tubs("stainless_steel_threaded_rod_tube",
							 stainless_steel_threaded_rod_inner_R,
							 stainless_steel_threaded_rod_outer_R,
							 stainless_steel_threaded_rod_half_length,
							 stainless_steel_threaded_rod_start_angle,
							 stainless_steel_threaded_rod_span_angle);
  
  G4LogicalVolume* stainless_steel_threaded_rod_log = new G4LogicalVolume(stainless_steel_threaded_rod_tube,
									  Cr20Ni8Fe76,
									  "stainless_steel_threaded_rod_log", 0, 0, 0);
  
  // top rod (as seen by the incoming beam)
  // 10.79 inches = 27.4066*cm, which is the length that centers the DSSD in SuN
  G4VPhysicalVolume* stainless_steel_threaded_rod_top = new G4PVPlacement(0, G4ThreeVector(0, dssd_cv_frame_dch, stainless_steel_threaded_rod_half_length - (stainless_steel_threaded_rod_length - 27.4066*cm)),
									  stainless_steel_threaded_rod_log,
									  "stainless_steel_threaded_rod_top", room_log, false, 0);
  // right rod (as seen by the incoming beam)
  // 10.79 inches = 27.4066*cm, which is the length that centers the DSSD in SuN
  G4VPhysicalVolume* stainless_steel_threaded_rod_right = new G4PVPlacement(0, G4ThreeVector(-1.0*dssd_cv_frame_dch, 0, stainless_steel_threaded_rod_half_length - (stainless_steel_threaded_rod_length - 27.4066*cm)),
									    stainless_steel_threaded_rod_log,
									    "stainless_steel_threaded_rod_right", room_log, false, 0);
  
























								 




								 



  
  ///////////////////////////////////////////
  // plastic washers on the rods that hold the DSSD and veto detector in place
  ///////////////////////////////////////////
  G4double rod_plastic_washer_inner_R = (2.40/2)*mm;
  G4double rod_plastic_washer_outer_R = 3*mm;
  G4double rod_plastic_washer_length = 2*mm;
  G4double rod_plastic_washer_half_length = 0.5*rod_plastic_washer_length;
  G4double rod_plastic_washer_start_angle = 0.0*deg;
  G4double rod_plastic_washer_span_angle = 360.0*deg;
  
  G4Tubs* rod_plastic_washer_tube = new G4Tubs("rod_plastic_washer_tube",
					       rod_plastic_washer_inner_R,
					       rod_plastic_washer_outer_R,
					       rod_plastic_washer_half_length,
					       rod_plastic_washer_start_angle,
					       rod_plastic_washer_span_angle);
  
  G4LogicalVolume* rod_plastic_washer_log = new G4LogicalVolume(rod_plastic_washer_tube,
								C2H4_rigid,
								"rod_plastic_washer_log", 0, 0, 0);
  
  // top, front side of the DSSD frame (as seen by the incoming beam)
  G4VPhysicalVolume* rod_plastic_washer_top_front_dssd = new G4PVPlacement(0, G4ThreeVector(0, dssd_cv_frame_dch, (-1.0)*(dssd_frame_whole_half_length + rod_plastic_washer_half_length)),
									   rod_plastic_washer_log,
									   "rod_plastic_washer_top_front_dssd", room_log, false, 0);

  // right, front side of the DSSD frame (as seen by the incoming beam)
  G4VPhysicalVolume* rod_plastic_washer_right_front_dssd = new G4PVPlacement(0, G4ThreeVector((-1.0)*dssd_cv_frame_dch, 0, (-1.0)*(dssd_frame_whole_half_length + rod_plastic_washer_half_length)),
									     rod_plastic_washer_log,
									     "rod_plastic_washer_right_front_dssd", room_log, false, 0);

  // top, back side of the DSSD frame (as seen by the incoming beam)
  G4VPhysicalVolume* rod_plastic_washer_top_back_dssd = new G4PVPlacement(0, G4ThreeVector(0, dssd_cv_frame_dch, dssd_frame_whole_half_length + rod_plastic_washer_half_length),
									  rod_plastic_washer_log,
									  "rod_plastic_washer_top_back_dssd", room_log, false, 0);
  
  // right, back side of the DSSD frame (as seen by the incoming beam)
  G4VPhysicalVolume* rod_plastic_washer_right_back_dssd = new G4PVPlacement(0, G4ThreeVector((-1.0)*dssd_cv_frame_dch, 0, dssd_frame_whole_half_length + rod_plastic_washer_half_length),
									    rod_plastic_washer_log,
									    "rod_plastic_washer_right_back_dssd", room_log, false, 0);
  
  // top, front side of the veto holder (as seen by the incoming beam)
  G4VPhysicalVolume* rod_plastic_washer_top_front_veto = new G4PVPlacement(0, G4ThreeVector(0, dssd_cv_frame_dch, db_dssd_and_veto_holder - veto_holder_whole_half_length - rod_plastic_washer_half_length),
									   rod_plastic_washer_log,
									   "rod_plastic_washer_top_front_veto", room_log, false, 0);

  // right, front side of the veto holder (as seen by the incoming beam)
  G4VPhysicalVolume* rod_plastic_washer_right_front_veto = new G4PVPlacement(0, G4ThreeVector((-1.0)*dssd_cv_frame_dch, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length - rod_plastic_washer_half_length),
									     rod_plastic_washer_log,
									     "rod_plastic_washer_right_front_veto", room_log, false, 0);

  // top, back side of the veto holder (as seen by the incoming beam)
  G4VPhysicalVolume* rod_plastic_washer_top_back_veto = new G4PVPlacement(0, G4ThreeVector(0, dssd_cv_frame_dch, db_dssd_and_veto_holder + veto_holder_whole_half_length + rod_plastic_washer_half_length),
									  rod_plastic_washer_log,
									  "rod_plastic_washer_top_back_veto", room_log, false, 0);

  // right, back side of the veto holder (as seen by the incoming beam)
  G4VPhysicalVolume* rod_plastic_washer_right_back_veto = new G4PVPlacement(0, G4ThreeVector((-1.0)*dssd_cv_frame_dch, 0, db_dssd_and_veto_holder + veto_holder_whole_half_length + rod_plastic_washer_half_length),
									    rod_plastic_washer_log,
									    "rod_plastic_washer_right_back_veto", room_log, false, 0);
  







  





  
  ///////////////////////////////////////////
  // brass hex nuts on the rods that hold the DSSD and veto detector in place
  ///////////////////////////////////////////

  G4double rod_hex_nut_length = 1.5*mm;
  G4double rod_hex_nut_half_length = 0.5*rod_hex_nut_length;
  G4double rod_hex_nut_start_angle = 0.0*deg;
  G4double rod_hex_nut_span_angle = 360.0*deg;
  G4int rod_hex_nut_num_sides = 6;
  G4int rod_hex_nut_num_Z_planes = 2;
  G4double rod_hex_nut_z_plane[] = {(-1.0)*(rod_hex_nut_half_length), rod_hex_nut_half_length};
  G4double rod_hex_nut_inner_r[] = {(2.40/2)*mm, (2.40/2)*mm};
  G4double rod_hex_nut_outer_r[] = {2*mm, 2*mm};
  
  G4Polyhedra* rod_hex_nut_poly = new G4Polyhedra("rod_hex_nut_poly",
						  rod_hex_nut_start_angle,
						  rod_hex_nut_span_angle,
						  rod_hex_nut_num_sides,
						  rod_hex_nut_num_Z_planes,
						  rod_hex_nut_z_plane,
						  rod_hex_nut_inner_r,
						  rod_hex_nut_outer_r);
  
  G4LogicalVolume* rod_hex_nut_log = new G4LogicalVolume(rod_hex_nut_poly,
							 Cu3Zn2,
							 "rod_hex_nut_log", 0, 0, 0);
  
  // top, front side of the DSSD frame (as seen by the incoming beam)
  G4VPhysicalVolume* rod_hex_nut_top_front_dssd = new G4PVPlacement(0, G4ThreeVector(0, dssd_cv_frame_dch, (-1.0)*(dssd_frame_whole_half_length + rod_plastic_washer_length + rod_hex_nut_half_length)),
								    rod_hex_nut_log,
								    "rod_hex_nut_top_front_dssd", room_log, false, 0);

  // right, front side of the DSSD frame (as seen by the incoming beam)
  G4VPhysicalVolume* rod_hex_nut_right_front_dssd = new G4PVPlacement(0, G4ThreeVector((-1.0)*dssd_cv_frame_dch, 0, (-1.0)*(dssd_frame_whole_half_length + rod_plastic_washer_length + rod_hex_nut_half_length)),
								      rod_hex_nut_log,
								      "rod_hex_nut_right_front_dssd", room_log, false, 0);  

  // top, back side of the DSSD frame (as seen by the incoming beam)
  G4VPhysicalVolume* rod_hex_nut_top_back_dssd = new G4PVPlacement(0, G4ThreeVector(0, dssd_cv_frame_dch, dssd_frame_whole_half_length + rod_plastic_washer_length + rod_hex_nut_half_length),
								   rod_hex_nut_log,
								   "rod_hex_nut_top_back_dssd", room_log, false, 0);
  
  // right, back side of the DSSD frame (as seen by the incoming beam)
  G4VPhysicalVolume* rod_hex_nut_right_back_dssd = new G4PVPlacement(0, G4ThreeVector((-1.0)*dssd_cv_frame_dch, 0, dssd_frame_whole_half_length + rod_plastic_washer_length + rod_hex_nut_half_length),
								     rod_hex_nut_log,
								     "rod_hex_nut_right_back_dssd", room_log, false, 0);  

  // top, front side of the veto holder (as seen by the incoming beam)
  G4VPhysicalVolume* rod_hex_nut_top_front_veto = new G4PVPlacement(0, G4ThreeVector(0, dssd_cv_frame_dch, db_dssd_and_veto_holder - veto_holder_whole_half_length - rod_plastic_washer_length - rod_hex_nut_half_length),
								    rod_hex_nut_log,
								    "rod_hex_nut_top_front_veto", room_log, false, 0);

  // right, front side of the veto holder (as seen by the incoming beam)
  G4VPhysicalVolume* rod_hex_nut_right_front_veto = new G4PVPlacement(0, G4ThreeVector((-1.0)*dssd_cv_frame_dch, 0, db_dssd_and_veto_holder - veto_holder_whole_half_length - rod_plastic_washer_length - rod_hex_nut_half_length),
								      rod_hex_nut_log,
								      "rod_hex_nut_right_front_veto", room_log, false, 0);

  // top, back side of the veto holder (as seen by the incoming beam)
  G4VPhysicalVolume* rod_hex_nut_top_back_veto = new G4PVPlacement(0, G4ThreeVector(0, dssd_cv_frame_dch, db_dssd_and_veto_holder + veto_holder_whole_half_length + rod_plastic_washer_length + rod_hex_nut_half_length),
								   rod_hex_nut_log,
								   "rod_hex_nut_top_back_veto", room_log, false, 0);

  // right, back side of the veto holder (as seen by the incoming beam)
  G4VPhysicalVolume* rod_hex_nut_right_back_veto = new G4PVPlacement(0, G4ThreeVector((-1.0)*dssd_cv_frame_dch, 0, db_dssd_and_veto_holder + veto_holder_whole_half_length + rod_plastic_washer_length + rod_hex_nut_half_length),
								     rod_hex_nut_log,
								     "rod_hex_nut_right_back_veto", room_log, false, 0);















  
  
  

  
// **************************************************************************************************
//  FOR THE CODE TO PROPERLY SAVE THINGS TO ROOT, YOU NEED TO TELL IT WHAT YOU NAMED YOUR DETECTORS
// **************************************************************************************************
//              The name of the detectors is in the G4VPhysicalVolume command                      
	  detectorName[0] = "T1";
	  detectorName[1] = "T2";
	  detectorName[2] = "T3";
	  detectorName[3] = "T4";
	  detectorName[4] = "B1";
	  detectorName[5] = "B2";
	  detectorName[6] = "B3";
	  detectorName[7] = "B4";
	  detectorName[8] = "DSSD"; // "dssd_silicon_chip_active_area"
	  detectorName[9] = "Veto"; // "veto_silicon_wafer"
	  
// **************************************************************************************************


/*
//____________THE TARGET HOLDER______________________

  G4double innerR_pipe = 3.64*0.5*mm;	 //tubes for air cooling   
  G4double outerR_pipe = 4.76*0.5*mm;               
  G4double halflength_pipe =879.39*0.5*mm;  

  G4double innerR_holder = 0.0*mm;	//holder  
  G4double outerR_holder = 17.15*mm;               
  G4double halflength_holder =4.76*0.5*mm;  

  G4double halfX_top = 11.52*mm;	//remove top part of holder
  G4double halfY_top = 5.0*mm;
  G4double halfZ_top = halflength_holder+1.0*mm;

  G4double halfX_indent = 17.17*0.5*mm;	//indent for target frame
  G4double halfY_indent = 20.73*0.5*mm;
  G4double halfZ_indent = 4.16*0.5*mm;

  G4double halfX_hole = 13.11*0.5*mm;	//hole for beam
  G4double halfY_hole = 16.77*0.5*mm;
  G4double halfZ_hole = halflength_holder+1.0*mm;

  G4double cutout_halfX1=13.465*mm;	//remove material 
  G4double cutout_halfX2=13.185*mm;
  G4double cutout_halfY1=2.38*0.5*mm;
  G4double cutout_halfY2=2.38*0.5*mm;
  G4double cutout_halfZ =2.54*0.5*mm;

  G4double innerR_cutout = 0.0*mm; 	//remove more material
  G4double outerR_cutout = 6.06*mm;               
  G4double halflength_cutout =2.54*0.5*mm;   

  G4double innerR_end = 0.0*mm;		//end where you hold it
  G4double outerR_end = 54.86*0.5*mm;               
  G4double halflength_end =3.05*0.5*mm;     

  G4double innerR1_endcap = 0.0*mm;	//endcap where you hold it
  G4double outerR1_endcap = 54.86*0.5*mm;   
  G4double innerR2_endcap = 0.0*mm;
  G4double outerR2_endcap = 39.70*0.5*mm;              
  G4double halflength_endcap = 2.03*0.5*mm;    

  G4double innerR_endHole1 = 0.0*mm;	//hole for air cooling
  G4double outerR_endHole1 = 3.64*0.5*mm;               
  G4double halflength_endHole1 =halflength_end+1.0*mm;   

  G4double innerR1_endHole2 = 0.0*mm;	//hole for air cooling
  G4double outerR1_endHole2 = 4.85*0.5*mm;   
  G4double innerR2_endHole2 = 0.0*mm;
  G4double outerR2_endHole2 = 6.38*0.5*mm;              
  G4double halflength_endHole2 =halflength_endcap+1.0*mm;      
          
  G4Tubs* pipe_tube = new G4Tubs("pipe_tube",innerR_pipe, outerR_pipe, halflength_pipe, 0.0*deg, 360.0*deg);
  G4Tubs* holder_tube = new G4Tubs("holder_tube",innerR_holder, outerR_holder, halflength_holder, 0.0*deg, 360.0*deg);
  G4Box* top_box = new G4Box("top_box",halfX_top,halfY_top,halfZ_top);
  G4Box* indent_box = new G4Box("indent_box",halfX_indent,halfY_indent,halfZ_indent);
  G4Box* hole_box = new G4Box("hole_box",halfX_hole,halfY_hole,halfZ_hole);
  G4Trd* cutout = new G4Trd("cutout",cutout_halfX1,cutout_halfX2,cutout_halfY1,cutout_halfY2,cutout_halfZ);
  G4Tubs* cutout_tube1 = new G4Tubs("cutout_tube",innerR_cutout, outerR_cutout, halflength_cutout, 90.0*deg, 90.0*deg);
  G4Tubs* cutout_tube2 = new G4Tubs("cutout_tube",innerR_cutout, outerR_cutout, halflength_cutout, 0.0*deg, 90.0*deg);
  G4Tubs* end_tube = new G4Tubs("end_tube",innerR_end, outerR_end, halflength_end, 0.0*deg, 360.0*deg);     
  G4Cons* endcap_cone = new G4Cons("endcap_cone",innerR1_endcap, outerR1_endcap, innerR2_endcap, outerR2_endcap, halflength_end, 0.0*deg, 360.0*deg);        
  G4Tubs* endHole_tube = new G4Tubs("endHole_tube",innerR_endHole1, outerR_endHole1, halflength_endHole1, 0.0*deg, 360.0*deg);    
  G4Cons* endHole_cone = new G4Cons("endHole_cone",innerR1_endHole2, outerR1_endHole2, innerR2_endHole2, outerR2_endHole2, halflength_endHole2, 0.0*deg, 360.0*deg);
 
  G4SubtractionSolid* holder_sub1 = new G4SubtractionSolid("holder_sub1",holder_tube,top_box,0,G4ThreeVector(0.0*mm,12.705*mm+halfY_top,0.0*mm));
  G4SubtractionSolid* holder_sub2 = new G4SubtractionSolid("holder_sub2",holder_sub1,indent_box,0,G4ThreeVector(0.0*mm,0.0*mm,halfZ_indent));
  G4SubtractionSolid* holder_sub3 = new G4SubtractionSolid("holder_sub3",holder_sub2,hole_box,0,G4ThreeVector(0.0*mm,0.0*mm,0.0*mm));           
  G4SubtractionSolid* holder_sub4 = new G4SubtractionSolid("holder_sub4",holder_sub3,cutout,0,G4ThreeVector(0.0*mm,4.705*mm,halflength_holder));    
  G4SubtractionSolid* holder_sub5 = new G4SubtractionSolid("holder_sub5",holder_sub4,cutout_tube1,rot_180,G4ThreeVector(7.035*mm,-3.425*mm,halflength_holder));
  G4SubtractionSolid* holder_sub6 = new G4SubtractionSolid("holder_sub6",holder_sub5,cutout_tube2,rot_180,G4ThreeVector(-7.035*mm,-3.425*mm,halflength_holder));

  G4SubtractionSolid* end_sub1 = new G4SubtractionSolid("end_sub1",end_tube,endHole_tube,0,G4ThreeVector(12.19*mm,0.0*mm,0.0*mm));
  G4SubtractionSolid* end_sub2 = new G4SubtractionSolid("end_sub2",end_sub1,endHole_tube,0,G4ThreeVector(-12.19*mm,0.0*mm,0.0*mm));

  G4SubtractionSolid* endcap_sub1 = new G4SubtractionSolid("endcap_sub1",endcap_cone,endHole_cone,0,G4ThreeVector(12.19*mm,0.0*mm,0.0*mm));
  G4SubtractionSolid* endcap_sub2 = new G4SubtractionSolid("endcap_sub2",endcap_sub1,endHole_cone,0,G4ThreeVector(-12.19*mm,0.0*mm,0.0*mm));

  G4LogicalVolume* pipe_log = new G4LogicalVolume(pipe_tube,Al,"pipe_log",0,0,0);
  G4LogicalVolume* holder_log = new G4LogicalVolume(holder_sub6,Al,"holder_log",0,0,0);
  G4LogicalVolume* end_log = new G4LogicalVolume(end_sub2,Al,"end_log",0,0,0);
  G4LogicalVolume* endcap_log = new G4LogicalVolume(endcap_sub2,Al,"endcap_log",0,0,0);

  G4VPhysicalVolume* holder_phys = new G4PVPlacement(0,G4ThreeVector(0.0*mm,0.0*mm,0.0*mm),holder_log,"holder_phys",room_log,false,0);

  G4VPhysicalVolume* pipe1_phys = new G4PVPlacement(0,G4ThreeVector(12.19*mm,0.0*mm,halflength_holder+halflength_pipe),pipe_log,"pipe1_phys",room_log,false,0);

  G4VPhysicalVolume* pipe2_phys = new G4PVPlacement(0,G4ThreeVector(-12.19*mm,0.0*mm,halflength_holder+halflength_pipe),pipe_log,"pipe2_phys",room_log,false,0);
 
  G4VPhysicalVolume* end_phys = new G4PVPlacement(0,G4ThreeVector(0.0*mm,0.0*mm,halflength_holder+2.0*halflength_pipe+halflength_end),end_log,"end_phys",room_log,false,0);

  G4VPhysicalVolume* endcap_phys = new G4PVPlacement(0,G4ThreeVector(0.0*mm,0.0*mm,halflength_holder+2.0*halflength_pipe+2.0*halflength_end+halflength_endcap),endcap_log,"endcap_phys",room_log,false,0);
*/



// VISUALIZATION STUFF

  room_log->SetVisAttributes (G4VisAttributes::Invisible);

//visualization for scintillators = GREEN
  G4VisAttributes *GreenAttr = new G4VisAttributes(G4Colour(0.,1.,0.));     
  GreenAttr->SetVisibility(true);
  GreenAttr->SetForceSolid(true);

//visualization for reflector = PURPLE
  G4VisAttributes *PurpleAttr = new G4VisAttributes(G4Colour(1.,0.,1.));  
  PurpleAttr->SetVisibility(true);
  PurpleAttr->SetForceSolid(true);

//visualization for aluminum = GREY
  G4VisAttributes *GreyAttr = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  GreyAttr->SetVisibility(true);
  GreyAttr->SetForceSolid(true);

//visualization for BLUE
  G4VisAttributes *BlueAttr = new G4VisAttributes(G4Colour(0.,0.,1.));
  BlueAttr->SetVisibility(true);
  BlueAttr->SetForceSolid(true);

//visualization for RED
  G4VisAttributes *RedAttr = new G4VisAttributes(G4Colour(1.,0.,0.));
  RedAttr->SetVisibility(true);
  RedAttr->SetForceSolid(true);

// applying the color scheme
  //scint_log->SetVisAttributes(GreenAttr);
  //refl_log->SetVisAttributes(PurpleAttr);
  //al_log->SetVisAttributes(GreyAttr);
  //al_log_side->SetVisAttributes(GreyAttr);
  //  beam_log->SetVisAttributes(BlueAttr);
//  holder_log->SetVisAttributes(GreyAttr);
//  pipe_log->SetVisAttributes(GreyAttr);
//  end_log->SetVisAttributes(GreyAttr);
//  endcap_log->SetVisAttributes(GreyAttr);

  

  
  dssd_frame_log->SetVisAttributes(BlueAttr);
  dssd_silicon_chip_non_active_area_log->SetVisAttributes(RedAttr);
  dssd_silicon_chip_active_area_log->SetVisAttributes(GreyAttr);
  //dssd_ribbon_cable_left_side_log->SetVisAttributes(PurpleAttr);
  //dssd_ribbon_cable_bent_left_side_log->SetVisAttributes(PurpleAttr);
  //dssd_ribbon_cable_bottom_side_log->SetVisAttributes(PurpleAttr);
  //dssd_ribbon_cable_bent_bottom_side_log->SetVisAttributes(PurpleAttr);
  
  //veto_holder_log->SetVisAttributes(BlueAttr);
  stainless_steel_threaded_rod_log->SetVisAttributes(GreyAttr);
  rod_plastic_washer_log->SetVisAttributes(RedAttr);
  rod_hex_nut_log->SetVisAttributes(GreenAttr);
  veto_housing_log->SetVisAttributes(GreenAttr);
  veto_silicon_wafer_log->SetVisAttributes(RedAttr);
  veto_housing_hex_nut_log->SetVisAttributes(BlueAttr);
  veto_housing_cable_connection_log->SetVisAttributes(PurpleAttr);
  microdot_to_lemo_cable_A_log->SetVisAttributes(RedAttr);
  microdot_to_lemo_cable_B_log->SetVisAttributes(BlueAttr);
  microdot_to_lemo_cable_C_log->SetVisAttributes(GreyAttr);
  microdot_to_lemo_cable_D_log->SetVisAttributes(BlueAttr);
  microdot_to_lemo_cable_E_log->SetVisAttributes(RedAttr);
  
  // Only for visualization
  dssd_frame_whole_log->SetVisAttributes(BlueAttr);
  dssd_frame_cutout_left_side_log->SetVisAttributes(GreenAttr);
  dssd_frame_cutout_bottom_side_log->SetVisAttributes(PurpleAttr);
  dssd_frame_additional_cutout_left_side_log->SetVisAttributes(RedAttr);
  dssd_frame_additional_cutout_bottom_side_log->SetVisAttributes(PurpleAttr);
  
  veto_holder_whole_log->SetVisAttributes(BlueAttr);
  veto_holder_silicon_cutout_log->SetVisAttributes(RedAttr);
  veto_holder_back_cutout_log->SetVisAttributes(GreenAttr);
  veto_holder_cutout_left_side_log->SetVisAttributes(PurpleAttr);
  veto_holder_cutout_bottom_side_log->SetVisAttributes(GreyAttr);

  veto_whole_log->SetVisAttributes(GreenAttr);
  veto_interior_cutout_log->SetVisAttributes(RedAttr);
  veto_window_cutout_log->SetVisAttributes(BlueAttr);
  
  return room_phys;

}
