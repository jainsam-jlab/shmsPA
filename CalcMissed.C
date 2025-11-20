//Calculation of missed triggers correction
//ACP 2025
//Takes name of file to analyze, central momentum, and particle type as arguments
//Should also change the analysis cuts: jump to section marked by -----------ooooooOOOOOooooo---------

void CalcMissed(TString filename = "PAOutData.root", Double_t P = 5, TString ptype = "pi+"){
	TString partList = "pi+K+p+e+";
	if (!partList.Contains(ptype)){
		cout<<"Not a valid particle type. Please enter a particle type of pi+, p+, K+, or e+"<<endl;
		cin>>ptype;
	}
	//import file, tree
	TFile *f1 = new TFile(filename);
	TTree *t1 = (TTree*)f1->Get("PA");
	//define variables and assign to tree
	Double_t s1x_time, s1x_energy;
	Double_t s1y_time, s1y_energy;
	Double_t s2x_time, s2x_energy;
	Double_t s2y_time, s2y_energy;
	Int_t NPE, copyNo;
	Int_t AGC_NPE, NGC_NPE, HGC_NPE;
	Double_t cal_energy;
	t1->SetBranchAddress("S1XEnergy",&s1x_energy);
	t1->SetBranchAddress("S1YEnergy",&s1y_energy);
	t1->SetBranchAddress("S2XEnergy",&s2x_energy);
	t1->SetBranchAddress("S2YEnergy",&s2y_energy);
	t1->SetBranchAddress("S1XTime",&s1x_time);
	t1->SetBranchAddress("S1YTime",&s1y_time);
	t1->SetBranchAddress("S2XTime",&s2x_time);
	t1->SetBranchAddress("S2YTime",&s2y_time);
	t1->SetBranchAddress("S2YNPE",&NPE);
	t1->SetBranchAddress("CopyNo",&copyNo);
	t1->SetBranchAddress("AGCNPE",&AGC_NPE);
	t1->SetBranchAddress("HGCNPE",&HGC_NPE);
	t1->SetBranchAddress("NGCNPE",&NGC_NPE);
	t1->SetBranchAddress("CalEnergy",&cal_energy);
	TH1I *h1 = new TH1I("h1","copyNo",60,0,60);
	//other things needed:
	Double_t dt1, dt2, dt3; //time differences
	Double_t energy_threshold = 0.5; //MeV
	Double_t time_window = 20; //ns
	Int_t npe_threshold = 100; //unitless
	Int_t nEvents = t1->GetEntries();
	//initialize counters to zero
	Int_t nMissed = 0; //missed triggers
	Int_t nTrig; //fired triggers per event
	Int_t nProton=0;
        Int_t nPion=0;
        Int_t nPositron=0;
        Int_t nKaon =0;
        Int_t nContam =0;
	Int_t nStopped = 0;
	Int_t nGoodEv = 0;
	Bool_t found; //for good particle
	//loop through tree to simulate 3/4 trigger
	for (int i=0; i<nEvents; i++){
		t1->GetEntry(i);
		//count raw missed triggers
		nTrig=0;
		if (s1x_energy > energy_threshold) nTrig++;
		dt1 = s1y_time - s1x_time;
                dt2 = s2x_time - s1x_time;
                dt3 = s2y_time - s1x_time;
		if (s1y_energy > energy_threshold && abs(dt1)<time_window)nTrig++;
		if (s2x_energy > energy_threshold && abs(dt2)<time_window)nTrig++;
		if (NPE > npe_threshold && abs(dt3)<time_window)nTrig++;
		if (nTrig<3)nMissed++;
		//count stopped particles
		if (copyNo!=0 && copyNo!=4)nStopped++;
		//determine signals of different particles
		if (nTrig>2){
			found = true;
		//USER EDITS: change cuts in if statements
		//--------ooooooOOOOOOooooo---------ooooOOOOOooooo---------ooooOOOOOooooo--------
		//proton: no signal in HGC or AGC
			if (AGC_NPE<20 && HGC_NPE <50){
				nProton++;
				found = true;
			}
		//kaon: signal in AGC, not in HGC
			if (AGC_NPE>20 && HGC_NPE<50){
				nKaon++;
				found = true;
			}
		//electron: signals in all Cherenkovs, lepton in Cal
			if (AGC_NPE>20 && HGC_NPE>50 && NGC_NPE>5 && cal_energy/P > 0.7){
				nPositron++;
				found = true;
			}
		//pion: signals in HGC and AGC
			 if (AGC_NPE>20 && HGC_NPE>50 && NGC_NPE<5){
				 nPion++;
				 found = true;
			 }
		//--------ooooooOOOOOOooooo---------ooooOOOOOooooo---------ooooOOOOOooooo--------
			 if (found==false)nContam++;

		}
	}
	//calculate correction factor
	//for correct particle type only
	if (ptype=="pi+") nGoodEv=nPion;
	if (ptype=="p+") nGoodEv=nProton;
	if (ptype=="K+") nGoodEv=nKaon;
	if (ptype=="e+") nGoodEv=nPositron;

	//calculations
	Double_t abs = (nEvents-nGoodEv)*100/(double)nEvents;
	Double_t abs_err = sqrt(nGoodEv)*100/(double)nEvents;
	Double_t corr = nGoodEv/(double)nEvents;
	Double_t corr_err = sqrt(nGoodEv)/(double)nEvents;
	
	//print results
	//input block
	cout<<endl;
	cout<<"Central momentum: "<<P<<" GeV"<<endl;
	cout<<"Incident particle type: "<<ptype<<endl;
	
	//interpretation block
	cout<<endl;
	cout<<"Total events: "<<nEvents<<endl;
	cout<<"Missed 3/4 triggers: "<<nMissed<<" or ";
	cout<<(double)nMissed*100/nEvents<<"%"<<endl;
	cout<<"Total stopped tracks: "<<nStopped<<" or "<<(double)nStopped*100/nEvents<<"%"<<endl;
	cout<<Form("3/4 was triggered by: %d p+, %d pi+, %d K+, %d positrons, and %d contaminated tracks",nProton,nPion,nKaon,nPositron,nContam)<<endl;
	
	//results block
	cout<<endl;
	cout<<"The fraction of events with a missing trigger or a secondary event of a different particle type is:"<<endl;
	cout<<Form("%5.3f +/- %5.3f %%",abs,abs_err)<<endl;
	cout<<"The correction factor to apply to your experimental yields is:"<<endl;
	cout<<Form("%5.5f +/- %5.5f",corr,corr_err)<<endl;
}


