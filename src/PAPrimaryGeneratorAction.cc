#include "PAPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <iostream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PAPrimaryGeneratorAction::PAPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),     
  fParticleGun(nullptr), fMessenger(nullptr), 
  fMomentum(5.*GeV)
{


  G4int nofParticles = 1;
  fParticleGun  = new G4ParticleGun(nofParticles);
  
  auto particleTable = G4ParticleTable::GetParticleTable();
  fProton = particleTable->FindParticle("proton");
  //start with proton
  auto mass = fProton->GetPDGMass();
  G4double Ekin = sqrt(fMomentum*fMomentum + mass*mass) - mass;
  
  //gun particle and position: does not need to change
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,ztar*m));
  fParticleGun->SetParticleDefinition(fProton);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
  fParticleGun->SetParticleEnergy(Ekin);
  
  // define commands for this class
  DefineCommands();

  //read in event parameters from a file
  //get information at initialization
  ReadGeneratedEvents();
  NEVENTS = delta.size();
  nev=0;

  fUseGenerated = false;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PAPrimaryGeneratorAction::~PAPrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PAPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  	auto mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
	//get information 
	if (!fUseGenerated){
		//leave default, consistent with /gun/ commands
  		fParticleGun->GeneratePrimaryVertex(event);
	} else {
		//and generate all the events read in from the file
		G4double P0 = fMomentum; //read central momentum
		G4double Pmom = P0*(1.0+delta[nev]/100.);
		//calculate gun position, momentum
		G4double posX = xfp[nev]/100. + xpfp[nev]*ztar; //x at ztar
		G4double posY = yfp[nev]/100. +ypfp[nev]*ztar; //y at ztar
		G4double momX = xpfp[nev]; //dx/dz
		G4double momY = ypfp[nev]; //dy/dz
		G4double momZ = 1; //dz/dz
		G4double Ekin = sqrt(Pmom*Pmom + mass*mass)-mass;
		auto Pvec = G4ThreeVector(momX, momY, momZ);
		auto Xvec = G4ThreeVector(posX*m,posY*m,ztar*m); //default z position
		//assign these quantities to the event
		fParticleGun->SetParticleMomentumDirection(Pvec);
		fParticleGun->SetParticleEnergy(Ekin);
		fParticleGun->SetParticlePosition(Xvec);

		fParticleGun->GeneratePrimaryVertex(event);
	   //output progress statement
	   if ((nev+1)%1000==0)G4cout<<"Processed "<<nev+1<<" events"<<G4endl;
	   //if using the entire file (>163333 events)
	   if (nev == NEVENTS-1) {
		   G4cout<< "End of generated events"<<G4endl;
		   return;
	   }
	}
	nev++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PAPrimaryGeneratorAction::DefineCommands()
{
  // Define /PA/generator command directory using generic messenger class
  fMessenger 
    = new G4GenericMessenger(this, 
                             "/PA/generator/", 
                             "Primary generator control");
            
  // momentum command
  auto& momentumCmd
    = fMessenger->DeclarePropertyWithUnit("momentum", "GeV", fMomentum, 
        "SHMS Central Momentum");
  momentumCmd.SetParameterName("p", true);
  momentumCmd.SetRange("p>=0.");                                
  momentumCmd.SetDefaultValue("1.");
  // ok
  //momentumCmd.SetParameterName("p", true);
  //momentumCmd.SetRange("p>=0.");     

  //choose to use generated events or not
  auto& genCmd = fMessenger->DeclareProperty("useGenerated",fUseGenerated);
  genCmd.SetGuidance("Use generated events from file");
  genCmd.SetParameterName("bool", true);
  genCmd.SetDefaultValue("false");

  auto& fileCmd = fMessenger->DeclareProperty("setInFile",fInputFileName);
  fileCmd.SetGuidance("Set file name for input kinematics");
  fileCmd.SetParameterName("string",true);
  fileCmd.SetDefaultValue("run5055kine.txt");
  
}


//..oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PAPrimaryGeneratorAction::ReadGeneratedEvents(){
  
  //loosely based off Bill's code
  //read genkinfile.txt, which has structure xfp, yfp, xpfp, ypfp, delta
  std::string line;
  std::ifstream fin( fInputFileName , std::ios::in );
  if( !fin.is_open() ) { G4cerr << "Can't open file" << G4endl; exit(1);}


  //create variables to hold data
  G4double delta_tmp, xfp_tmp, yfp_tmp, xpfp_tmp, ypfp_tmp;

  NEVENTS=0;

  //for every line in file
  while (std::getline(fin, line)) {
     NEVENTS++; //increment event counter

     //read the line
    sscanf( line.c_str() , "%lf %lf %lf %lf %lf" ,
            &xfp_tmp, &yfp_tmp, &xpfp_tmp,  &ypfp_tmp, &delta_tmp);
    delta.push_back(delta_tmp);
    xfp.push_back(xfp_tmp);
    yfp.push_back(yfp_tmp);
    xpfp.push_back(xpfp_tmp);
    ypfp.push_back(ypfp_tmp);

    }
}



//..oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
