// 2018.10.10
//
// E498 S.Nakamura
//

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fDefaultMaterial(0),
 fSolidWorld(0),fLogicWorld(0),fPhysiWorld(0),
 fLaBr3CrystalMaterial(0),fAlCaseMaterial(0),
 fSolidLaBr3Det(0), fLogicLaBr3Det(0), fPhysiLaBr3Det(0),
 fSolidLaBr3Crystal(0), fLogicLaBr3Crystal(0), fPhysiLaBr3Crystal(0),
 fSolidAlCase(0), fLogicAlCase(0), fPhysiAlCase(0),
 fSolidFaceAlCase(0), fLogicFaceAlCase(0), fPhysiFaceAlCase(0),
 fDetectorMessenger(0)
{
  // default parameter values
  ComputeLaBr3Parameters();

  fWorldSizeX = 1000*cm;
  fWorldSizeY = 1000*cm;
  fWorldSizeZ = 1000*cm;

  fLaBr3CrystalLength   = 20.3*cm;
  fLaBr3CrystalDiameter = 8.9*cm;
  fAlCaseDiameter       = 12.0*cm;
  fAlCaseInnerDiameter  = 10.03*cm;
  fAlCaseLength         = 41.2*cm;
  fFaceAlCaseThickness  = 0.5*cm;

  // Definition of materials
  DefineMaterials();
  fDefaultMaterial      = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  fLaBr3CrystalMaterial = G4NistManager::Instance()->FindOrBuildMaterial("LaBr3");
  fAlCaseMaterial       = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

  // create commands for interactive definition
  fDetectorMessenger = new DetectorMessenger(this);
}


DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructLaBr3();
}


void DetectorConstruction::DefineMaterials()
{ 
  // use G4-NIST materials data base
  G4String symbol;
  G4double a, z, density;
  G4double fractionmass;
  G4int ncomp;

  //Definition of LaBr3
  G4Element*  La    = new G4Element ("Lanthanum", symbol="La", z=57, a=138.91*g/mole);
  G4Element*  Br    = new G4Element ("Bromine",   symbol="Br", z=35, a=79.90*g/mole);
  G4Element*  Ce    = new G4Element ("Cerium",    symbol="Ce", z=58, a=140.11*g/mole);
  G4Material* LaBr3 = new G4Material("LaBr3", density=5.1*g/cm3, ncomp=3);
  LaBr3->AddElement(La, fractionmass = 0.34855);
  LaBr3->AddElement(Br, fractionmass = 0.60145);
  LaBr3->AddElement(Ce, fractionmass = 0.05);
  LaBr3->GetIonisation()->SetMeanExcitationEnergy(454.5*eV);

  // print table
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


void DetectorConstruction::ComputeLaBr3Parameters()
{
  fLaBr3DetLength   =  fFaceAlCaseThickness + fAlCaseLength;  //Length of LaBr3 detector system
  fLaBr3DetDiameter =  fAlCaseDiameter;                       //Diameter of LaBr3 detector system
  fZposLaBr3Det     = -0.5*fLaBr3DetLength - 0.5*cm;          //Z position of LaBr3 detector system

  //Dont't touch below sentences! This is relative position of each components of LaBr3 det.   
  fZposFaceAlCase   =  0.5*fLaBr3DetLength - 0.5*fFaceAlCaseThickness;
  fZposAlCase       =  0.5*fLaBr3DetLength - fFaceAlCaseThickness - 0.5*fAlCaseLength;
  fZposLaBr3Crystal =  0.5*fLaBr3DetLength - fFaceAlCaseThickness - 0.1*cm - fLaBr3CrystalLength*0.5;
}

G4VPhysicalVolume* DetectorConstruction::ConstructLaBr3()
{
  // Clean old geometry, if any
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // complete the parameters definition
  ComputeLaBr3Parameters();
     
  // World
  fSolidWorld = new G4Box("World",                                //its name
                   fWorldSizeX/2,fWorldSizeY/2,fWorldSizeZ/2);  //its size
                         
  fLogicWorld = new G4LogicalVolume(fSolidWorld,            //its solid
                                   fDefaultMaterial,        //its material
                                   "World");                //its name
                                   
  fPhysiWorld = new G4PVPlacement(0,                        //no rotation
                                 G4ThreeVector(),           //at (0,0,0)
                                 fLogicWorld,             //its logical volume
                                 "World",                   //its name
                                 0,                         //its mother  volume
                                 false,                  //no boolean operation
                                 0);                        //copy number

  // Milano LaBr3 detector
  G4RotationMatrix rotm   = G4RotationMatrix(0, 0, 0);
  G4ThreeVector position  = G4ThreeVector(0., 0., fZposLaBr3Det);
  G4Transform3D transform = G4Transform3D(rotm, position);          

  fSolidLaBr3Det = new G4Tubs("LaBr3Det",
			      0.,
			      0.5*fLaBr3DetDiameter,
			      0.5*fLaBr3DetLength,
			      0.*deg,
			      360.*deg);
  
  fLogicLaBr3Det = new G4LogicalVolume(fSolidLaBr3Det,        //its solid
				       fDefaultMaterial,   //its material
				       "LaBr3Det");     //its name
  
  fPhysiLaBr3Det = new G4PVPlacement(transform,                    //no rotation
				     fLogicLaBr3Det,           //its logical volume
				     "LaBr3Det",         //its name
				     fLogicWorld,           //its mother  volume
				     false,              //no boolean operation
				     0);                    //copy number

  // Milano LaBr3 Crystal
  fSolidLaBr3Crystal = new G4Tubs("LaBr3Crystal",                //its name
				  0.,
				  0.5*fLaBr3CrystalDiameter,
				  0.5*fLaBr3CrystalLength,
				  0.*deg,
				  360.*deg);
  
  fLogicLaBr3Crystal = new G4LogicalVolume(fSolidLaBr3Crystal,        //its solid
					   fLaBr3CrystalMaterial,   //its material
					   "LaBr3Crystal");     //its name
  
  fPhysiLaBr3Crystal = new G4PVPlacement(0,                    //no rotation
					 G4ThreeVector(0.,0.,fZposLaBr3Crystal),
					 fLogicLaBr3Crystal,           //its logical volume
					 "LaBr3Crystal",         //its name
					 fLogicLaBr3Det,           //its mother  volume
					 false,              //no boolean operation
					 0);                    //copy number

  // Milano LaBr3 Al Case
  fSolidAlCase = new G4Tubs("AlCase",                //its name
			    0.5*fAlCaseInnerDiameter,
			    0.5*fAlCaseDiameter,
			    0.5*fAlCaseLength,
			    0.*deg,
			    360.*deg);
  
  fLogicAlCase = new G4LogicalVolume(fSolidAlCase,        //its solid
				     fAlCaseMaterial,   //its material
				     "AlCase");     //its name
  
  fPhysiAlCase = new G4PVPlacement(0,                    //no rotation
				   G4ThreeVector(0.,0.,fZposAlCase),
				   fLogicAlCase,           //its logical volume
				   "AlCase",         //its name
				   fLogicLaBr3Det,           //its mother  volume
				   false,              //no boolean operation
				   0);                    //copy number

  // Milano LaBr3 Face Al Case
  fSolidFaceAlCase = new G4Tubs("FaceAlCase",                //its name
				0.,
				0.5*fAlCaseDiameter,
				0.5*fFaceAlCaseThickness,
				0.*deg,
				360.*deg);
  
  fLogicFaceAlCase = new G4LogicalVolume(fSolidFaceAlCase,        //its solid
					 fAlCaseMaterial,   //its material
					 "FaceAlCase");     //its name
  
  fPhysiFaceAlCase = new G4PVPlacement(0,                    //no rotation
				       G4ThreeVector(0.,0.,fZposFaceAlCase),
				       fLogicFaceAlCase,           //its logical volume
				       "FaceAlCase",         //its name
				       fLogicLaBr3Det,           //its mother  volume
				       false,              //no boolean operation
				       0);                    //copy number
  
  PrintLaBr3Parameters();     

  // Visualization attributes
  fLogicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());

  G4VisAttributes* simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);

  G4VisAttributes* WireFrameVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.0));
  WireFrameVisAtt->SetForceWireframe(true);
  fLogicLaBr3Det->SetVisAttributes(WireFrameVisAtt);

  G4VisAttributes* RedVisAtt       = new G4VisAttributes(G4Colour(0.75,0.,0.));
  fLogicAlCase->SetVisAttributes(RedVisAtt);
  fLogicFaceAlCase->SetVisAttributes(RedVisAtt);

  G4VisAttributes* BlueVisAtt      = new G4VisAttributes(G4Colour(0.,0.,0.75));
  fLogicLaBr3Crystal->SetVisAttributes(BlueVisAtt);

  //always return the physical World
  return fPhysiWorld;
}


void DetectorConstruction::PrintLaBr3Parameters()
{
}


void DetectorConstruction::SetFaceAlCaseThickness(G4double val)
{
  fFaceAlCaseThickness = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetLaBr3CrystalDiameter(G4double val)
{
  fLaBr3CrystalDiameter = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetLaBr3CrystalLength(G4double val)
{
  fLaBr3CrystalLength = val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
