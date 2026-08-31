#include <TFile.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct PA {
  int id{}, agc{}, hgc{}, ngc{}, s2ynpe{};
  double s1xe{}, s1ye{}, s2xe{}, s2ye{}, s1xt{}, s1yt{}, s2xt{}, s2yt{}, cal{};
};
struct Tr {
  int id{}, parent{}, pdg{};
  std::string creator, endVol, endProc;
  bool s1{}, hgc{}, agc{}, s2{}, cal{}, exit{};
};
enum PID { PION, PROTON, KAON, POSITRON, CONTAM, TRIGFAIL };
static int trig(const PA& p) {
  return (p.s1xe > .5) + (p.s1ye > .5 && std::abs(p.s1yt-p.s1xt) < 20) +
         (p.s2xe > .5 && std::abs(p.s2xt-p.s1xt) < 20) +
         (p.s2ynpe > 100 && std::abs(p.s2yt-p.s1xt) < 20);
}
static int pid(const PA& p, double mom) {
  if (trig(p) < 3) return TRIGFAIL;
  if (p.agc < 20 && p.hgc < 50) return PROTON;
  if (p.agc > 20 && p.hgc < 50) return KAON;
  if (p.agc > 20 && p.hgc > 50 && p.ngc > 5 && p.cal/(mom*1000.) > .7) return POSITRON;
  if (p.agc > 20 && p.hgc > 50) return PION;
  return CONTAM;
}
static int stage(const Tr& t) {
  if (t.exit) return 6; if (t.cal) return 5; if (t.s2) return 4;
  if (t.agc) return 3; if (t.hgc) return 2; if (t.s1) return 1; return 0;
}
static const char* stageName(int s) {
  static const char* n[]={"BeforeS1","S1","HGC","AGC","S2","Cal","ExitedCal"};
  return n[s];
}

int main(int argc, char** argv) {
  if (argc < 5) {
    std::cerr << "usage: " << argv[0]
              << " INPUT.root OUTPUT.root CENTRAL_MOMENTUM_GEV TARGET_PDG [MAX_TRACK_ENTRIES]\n";
    return 2;
  }
  const std::string in=argv[1], out=argv[2];
  const double mom=std::stod(argv[3]);
  int targetPDG=std::stoi(argv[4]);
  const long long max=argc>5 ? std::stoll(argv[5]) : -1;
  TStopwatch total, phase; total.Start(); phase.Start();

  TFile fi(in.c_str(),"READ");
  auto* pt=dynamic_cast<TTree*>(fi.Get("PA"));
  auto* tt=dynamic_cast<TTree*>(fi.Get("Tracks"));
  if (fi.IsZombie() || !pt || !tt) { std::cerr<<"Cannot open PA and Tracks in "<<in<<"\n"; return 1; }

  PA p;
#define PAB(n,v) pt->SetBranchAddress(n,&v)
  PAB("EventID",p.id); PAB("AGCNPE",p.agc); PAB("HGCNPE",p.hgc); PAB("NGCNPE",p.ngc);
  PAB("S2YNPE",p.s2ynpe); PAB("S1XEnergy",p.s1xe); PAB("S1YEnergy",p.s1ye);
  PAB("S2XEnergy",p.s2xe); PAB("S2YEnergy",p.s2ye); PAB("S1XTime",p.s1xt);
  PAB("S1YTime",p.s1yt); PAB("S2XTime",p.s2xt); PAB("S2YTime",p.s2yt);
  PAB("CalEnergy",p.cal);
  std::unordered_map<int,PA> pas; pas.reserve(pt->GetEntries()*2);
  for(long long i=0;i<pt->GetEntries();++i){pt->GetEntry(i);pas[p.id]=p;}
  std::cout<<"TIMING read_PA_s "<<phase.RealTime()<<" events "<<pas.size()<<"\n"; phase.Start(true);

  TFile fo(out.c_str(),"RECREATE");
  TTree es("EventSummary","Selected-particle truth correlated with event detector response");
  TTree dd("DirectDaughters","Direct daughters with the selected PDG");
  TTree at("TargetTracks","Tracks with the selected PDG");

  int eid,primaryID,pS1,pHGC,pAGC,pS2,pCal,pExit,nDirect,hasDirect,dHGC,dAGC,dS2,dCal,dExit;
  int bestID,bestStage,nTarget,nDesc,anyHGC,anyAGC,anyS2,anyCal,descHGC,descAGC,descS2,descCal;
  int primaryPID,directPID,descPID,anyPID,ts1x,ts1y,ts2x,ts2y,ntr,tpass,oldClass,oldPion;
  std::string pEndVol,pEndProc,bestStageName;
#define EB(n,v) es.Branch(n,&v)
  EB("EventID",eid); EB("TargetPDG",targetPDG); EB("PrimaryTrackID",primaryID);
  EB("PrimaryReachedS1",pS1); EB("PrimaryReachedHGC",pHGC); EB("PrimaryReachedAGC",pAGC);
  EB("PrimaryReachedS2",pS2); EB("PrimaryReachedCal",pCal); EB("PrimaryExitedCal",pExit);
  EB("PrimaryEndVolume",pEndVol); EB("PrimaryEndProcess",pEndProc);
  EB("NDirectDaughters",nDirect); EB("HasDirectDaughter",hasDirect);
  EB("DirectDaughterReachedHGC",dHGC); EB("DirectDaughterReachedAGC",dAGC);
  EB("DirectDaughterReachedS2",dS2); EB("DirectDaughterReachedCal",dCal);
  EB("DirectDaughterExitedCal",dExit); EB("BestDirectDaughterTrackID",bestID);
  EB("BestDirectDaughterFurthestStage",bestStage); EB("BestDirectDaughterFurthestStageName",bestStageName);
  EB("NTargetTracks",nTarget); EB("NDescendantTracks",nDesc);
  EB("AnyTargetReachedHGC",anyHGC); EB("AnyTargetReachedAGC",anyAGC);
  EB("AnyTargetReachedS2",anyS2); EB("AnyTargetReachedCal",anyCal);
  EB("AnyDescendantReachedHGC",descHGC); EB("AnyDescendantReachedAGC",descAGC);
  EB("AnyDescendantReachedS2",descS2); EB("AnyDescendantReachedCal",descCal);
  EB("PrimaryInPIDRegion",primaryPID); EB("DirectDaughterInPIDRegion",directPID);
  EB("DescendantInPIDRegion",descPID); EB("AnyTargetInPIDRegion",anyPID);
  EB("TriggerS1X",ts1x); EB("TriggerS1Y",ts1y); EB("TriggerS2X",ts2x); EB("TriggerS2Y",ts2y);
  EB("NTrig",ntr); EB("TriggerPass",tpass); EB("OldPIDClass",oldClass); EB("OldPionPIDPass",oldPion);
  EB("AGCNPE",p.agc); EB("HGCNPE",p.hgc); EB("NGCNPE",p.ngc); EB("S2YNPE",p.s2ynpe);
  EB("S1XEnergy",p.s1xe); EB("S1YEnergy",p.s1ye); EB("S2XEnergy",p.s2xe); EB("S2YEnergy",p.s2ye);
  EB("S1XTime",p.s1xt); EB("S1YTime",p.s1yt); EB("S2XTime",p.s2xt); EB("S2YTime",p.s2yt);
  EB("CalEnergy",p.cal);

  int deid,dprim,did,dparent,dstage,dS1,dH,dA,dS,dC,dE,dTrig,dNTrig,dANPE,dHNPE,dNNPE,dOldPion;
  std::string dCreator,dEndVol,dEndProc,dStageName;
#define DB(n,v) dd.Branch(n,&v)
  DB("EventID",deid); DB("TargetPDG",targetPDG); DB("PrimaryTrackID",dprim);
  DB("DaughterTrackID",did); DB("DaughterParentID",dparent); DB("CreatorProcess",dCreator);
  DB("ReachedS1",dS1); DB("ReachedHGC",dH); DB("ReachedAGC",dA); DB("ReachedS2",dS);
  DB("ReachedCal",dC); DB("ExitedCal",dE); DB("EndVolume",dEndVol); DB("EndProcess",dEndProc);
  DB("FurthestStage",dstage); DB("FurthestStageName",dStageName); DB("TriggerPass",dTrig);
  DB("NTrig",dNTrig); DB("AGCNPE",dANPE); DB("HGCNPE",dHNPE); DB("NGCNPE",dNNPE);
  DB("OldPionPIDPass",dOldPion);

  int aeid,aid,apar,agen,aPrimary,aDirect,aDesc,aS1,aH,aA,aS,aC,aE;
  std::string aCreator,aEndVol,aEndProc;
#define AB(n,v) at.Branch(n,&v)
  AB("EventID",aeid); AB("TargetPDG",targetPDG); AB("TrackID",aid); AB("ParentID",apar);
  AB("GenerationDepth",agen); AB("IsPrimary",aPrimary); AB("IsDirectDaughter",aDirect);
  AB("IsPrimaryDescendant",aDesc); AB("ReachedS1",aS1); AB("ReachedHGC",aH);
  AB("ReachedAGC",aA); AB("ReachedS2",aS); AB("ReachedCal",aC); AB("ExitedCal",aE);
  AB("CreatorProcess",aCreator); AB("EndVolume",aEndVol); AB("EndProcess",aEndProc);

  int te,tid,tpar,tpdg,rs1,rh,ra,rs2,rc,re;
  char creator[256]{},ev[256]{},eproc[256]{};
  tt->SetCacheSize(32LL*1024*1024); tt->SetBranchStatus("*",0);
#define TB(n,v) tt->SetBranchStatus(n,1); tt->SetBranchAddress(n,v)
  TB("EventID",&te); TB("TrackID",&tid); TB("ParentID",&tpar); TB("PDG",&tpdg);
  TB("CreatorProcess",creator); TB("EndVolume",ev); TB("EndProcess",eproc);
  TB("ReachedS1",&rs1); TB("ReachedHGC",&rh); TB("ReachedAGC",&ra);
  TB("ReachedS2",&rs2); TB("ReachedCal",&rc); TB("ExitedCal",&re);

  long long missingPA=0,badPrimary=0,duplicates=0,events=0;
  std::vector<Tr> tracks; int current=-1;
  auto process=[&](){
    if(tracks.empty()) return;
    eid=current; auto pi=pas.find(eid); if(pi==pas.end()){++missingPA;return;} p=pi->second;
    std::unordered_map<int,size_t> byid; byid.reserve(tracks.size()*2);
    for(size_t i=0;i<tracks.size();++i) if(!byid.emplace(tracks[i].id,i).second) ++duplicates;
    std::vector<size_t> prim;
    for(size_t i=0;i<tracks.size();++i) if(tracks[i].parent==0&&tracks[i].pdg==targetPDG) prim.push_back(i);
    if(prim.size()!=1){++badPrimary;return;}
    const Tr& pr=tracks[prim[0]]; primaryID=pr.id; pS1=pr.s1;pHGC=pr.hgc;pAGC=pr.agc;
    pS2=pr.s2;pCal=pr.cal;pExit=pr.exit;pEndVol=pr.endVol;pEndProc=pr.endProc;
    nDirect=dHGC=dAGC=dS2=dCal=dExit=0;bestID=-1;bestStage=-1;
    nTarget=nDesc=anyHGC=anyAGC=anyS2=anyCal=descHGC=descAGC=descS2=descCal=0;
    auto ancestry=[&](const Tr&t,int&depth){depth=0;int par=t.parent;std::unordered_set<int> guard;
      while(par&&guard.insert(par).second){++depth;if(par==primaryID)return true;
        auto q=byid.find(par);if(q==byid.end())return false;par=tracks[q->second].parent;}return false;};
    for(const Tr&t:tracks) if(t.pdg==targetPDG){
      ++nTarget;int depth=0;bool isp=t.parent==0&&t.id==primaryID,direct=t.parent==primaryID;
      bool desc=direct||ancestry(t,depth);if(direct)depth=1;if(desc&&!isp)++nDesc;
      anyHGC|=t.hgc;anyAGC|=t.agc;anyS2|=t.s2;anyCal|=t.cal;
      if(desc&&!isp){descHGC|=t.hgc;descAGC|=t.agc;descS2|=t.s2;descCal|=t.cal;}
      aeid=eid;aid=t.id;apar=t.parent;agen=isp?0:depth;aPrimary=isp;aDirect=direct;aDesc=desc&&!isp;
      aS1=t.s1;aH=t.hgc;aA=t.agc;aS=t.s2;aC=t.cal;aE=t.exit;aCreator=t.creator;
      aEndVol=t.endVol;aEndProc=t.endProc;at.Fill();
      if(direct){++nDirect;dHGC|=t.hgc;dAGC|=t.agc;dS2|=t.s2;dCal|=t.cal;dExit|=t.exit;
        int st=stage(t);if(st>bestStage){bestStage=st;bestID=t.id;}
        deid=eid;dprim=primaryID;did=t.id;dparent=t.parent;dCreator=t.creator;dS1=t.s1;
        dH=t.hgc;dA=t.agc;dS=t.s2;dC=t.cal;dE=t.exit;dEndVol=t.endVol;dEndProc=t.endProc;
        dstage=st;dStageName=stageName(st);dNTrig=trig(p);dTrig=dNTrig>=3;dANPE=p.agc;
        dHNPE=p.hgc;dNNPE=p.ngc;dOldPion=pid(p,mom)==PION;dd.Fill();}
    }
    hasDirect=nDirect>0;directPID=dHGC||dAGC;descPID=descHGC||descAGC;
    anyPID=anyHGC||anyAGC;primaryPID=pHGC||pAGC;bestStageName=bestStage>=0?stageName(bestStage):"None";
    ts1x=p.s1xe>.5;ts1y=p.s1ye>.5&&std::abs(p.s1yt-p.s1xt)<20;
    ts2x=p.s2xe>.5&&std::abs(p.s2xt-p.s1xt)<20;ts2y=p.s2ynpe>100&&std::abs(p.s2yt-p.s1xt)<20;
    ntr=ts1x+ts1y+ts2x+ts2y;tpass=ntr>=3;oldClass=pid(p,mom);oldPion=oldClass==PION;
    es.Fill();++events;
  };
  const long long entries=max>=0?std::min<long long>(max,tt->GetEntries()):tt->GetEntries();
  for(long long i=0;i<entries;++i){tt->GetEntry(i);if(current!=-1&&te!=current){process();tracks.clear();}
    current=te;tracks.push_back({tid,tpar,tpdg,creator,ev,eproc,(bool)rs1,(bool)rh,(bool)ra,
      (bool)rs2,(bool)rc,(bool)re});if(i&&i%50000000==0)std::cout<<"PROGRESS "<<i<<"/"<<entries<<"\n";}
  process();
  std::cout<<"TIMING read_tracks_classify_s "<<phase.RealTime()<<" tracks "<<entries<<" events "<<events<<"\n";
  phase.Start(true);fo.cd();es.Write();dd.Write();at.Write();
  TTree meta("Metadata","study metadata");std::string input=in,pidDef="ReachedHGC OR ReachedAGC";
  std::string warning="AGCNPE/HGCNPE are event-level; target-track presence is correlation, not NPE attribution";
  long long read=entries;meta.Branch("InputFile",&input);meta.Branch("TargetPDG",&targetPDG);
  meta.Branch("PIDRegionDefinition",&pidDef);meta.Branch("NPEInterpretation",&warning);
  meta.Branch("TrackEntriesRead",&read);meta.Branch("MissingPAEvents",&missingPA);
  meta.Branch("MalformedPrimaryEvents",&badPrimary);meta.Branch("DuplicateEventTrackIDs",&duplicates);
  meta.Fill();meta.Write();fo.Close();
  std::cout<<"TIMING write_ROOT_s "<<phase.RealTime()<<"\nTOTAL_s "<<total.RealTime()<<"\nOUTPUT "<<out<<"\n";
  return (missingPA||badPrimary||duplicates)?3:0;
}
