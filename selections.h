// See selections.C file for details

#ifndef selections_h
#define selections_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class selections {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TTree          *nRooTracker;
   TTree          *outChain;
   Int_t           fCurrent; //!current Tree number in a TChain
   Int_t           fCurrent2; //!current Tree number in a TChain
   Int_t           fCurrent3; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           evt;
   Int_t           nSubEvts;
   Int_t           nClusters;
   Int_t           recoNRings[20]; //[nClusters]
   Int_t           subevt[20];   //[nSubevents]
   Int_t           ring[20];   //[nSubevents]
   Bool_t          isHighE[20];   //[nSubevents]
   Int_t           ringPEs[20];   //[nSubevents]
   Double_t        trueVtxX;
   Double_t        trueVtxY;
   Double_t        trueVtxZ;
   Double_t        trueTime;
   Double_t        trueDirX;
   Double_t        trueDirY;
   Double_t        trueDirZ;
   Double_t        trueKE;
   Double_t        diffKE;
   Double_t        diffVtxX;
   Double_t        diffVtxY;
   Double_t        diffVtxZ;
   Double_t        diffTime;
   Double_t        diffDirX;
   Double_t        diffDirY;
   Double_t        diffDirZ;
   Double_t        diffVtxAbs;
   Double_t        diffDirAbs;
   Int_t           mode;
   Double_t        neutrinoE;
   Double_t        neutrinoDirX;
   Double_t        neutrinoDirY;
   Double_t        neutrinoDirZ;
   Int_t           neutrinoPID;
   Int_t           nNeutrons;
   Int_t           neutronCount;
   Int_t           nCaptures;
   Double_t        trueToWall;
   Double_t        trueDWall;
   Double_t        recoToWall;
   Double_t        recoDWall;
   Int_t           recoCaptures;

   Double_t        recoVtxXLowE[20];   //[nClusters]
   Double_t        recoVtxYLowE[20];   //[nClusters]
   Double_t        recoVtxZLowE[20];   //[nClusters]
   Double_t        recoTimeLowE[20];   //[nClusters]
   Double_t        recoDirXLowE[20];   //[nClusters]
   Double_t        recoDirYLowE[20];   //[nClusters]
   Double_t        recoDirZLowE[20];   //[nClusters]
   Double_t        recoChkvAngleLowE[20];   //[nClusters]
   Double_t        recoEnergyLowE[20];   //[nClusters]

   Double_t        recoVtxXHighEMuon[20];   //[nSubevents]
   Double_t        recoVtxYHighEMuon[20];   //[nSubevents]
   Double_t        recoVtxZHighEMuon[20];   //[nSubevents]
   Double_t        recoTimeHighEMuon[20];   //[nSubevents]
   Double_t        recoDirXHighEMuon[20];   //[nSubevents]
   Double_t        recoDirYHighEMuon[20];   //[nSubevents]
   Double_t        recoDirZHighEMuon[20];   //[nSubevents]
   Double_t        recoChkvAngleHighEMuon[20];   //[nSubevents]
   Double_t        recoEnergyHighEMuon[20];   //[nSubevents]
   Double_t        recoLnLHighEMuon[20];   //[nSubevents]

   Double_t        recoVtxXHighEElectron[20];   //[nSubevents]
   Double_t        recoVtxYHighEElectron[20];   //[nSubevents]
   Double_t        recoVtxZHighEElectron[20];   //[nSubevents]
   Double_t        recoTimeHighEElectron[20];   //[nSubevents]
   Double_t        recoDirXHighEElectron[20];   //[nSubevents]
   Double_t        recoDirYHighEElectron[20];   //[nSubevents]
   Double_t        recoDirZHighEElectron[20];   //[nSubevents]
   Double_t        recoChkvAngleHighEElectron[20];   //[nSubevents]
   Double_t        recoEnergyHighEElectron[20];   //[nSubevents]
   Double_t        recoLnLHighEElectron[20];   //[nSubevents]

   Double_t        recoVtxX[20];   //[nSubevents]
   Double_t        recoVtxY[20];   //[nSubevents]
   Double_t        recoVtxZ[20];   //[nSubevents]
   Double_t        recoTime[20];   //[nSubevents]
   Double_t        recoDirX[20];   //[nSubevents]
   Double_t        recoDirY[20];   //[nSubevents]
   Double_t        recoDirZ[20];   //[nSubevents]
   Double_t        recoChkvAngle[20];   //[nSubevents]
   Double_t        recoEnergy[20];   //[nSubevents]
   Int_t           recoPID[20];   //[nSubevents]

   Int_t           NEneutmode;
   Int_t           NEnvc;
   Int_t           NEipvc[100];   //[NEnvc]
   Float_t         NEpvc[100][3];   //[NEnvc]
   Int_t           NEiorgvc[100];   //[NEnvc]
   Int_t           NEiflgvc[100];   //[NEnvc]
   Int_t           NEicrnvc[100];   //[NEnvc]
   Float_t         NEcrsx;
   Float_t         NEcrsy;
   Float_t         NEcrsz;
   Float_t         NEcrsphi;
   Int_t           NEnvert;
   Float_t         NEposvert[100][3];   //[NEnvert]
   Int_t           NEiflgvert[100];   //[NEnvert]
   Int_t           NEnvcvert;
   Float_t         NEdirvert[100][3];   //[NEnvcvert]
   Float_t         NEabspvert[100];   //[NEnvcvert]
   Float_t         NEabstpvert[100];   //[NEnvcvert]
   Int_t           NEipvert[100];   //[NEnvcvert]
   Int_t           NEiverti[100];   //[NEnvcvert]
   Int_t           NEivertf[100];   //[NEnvcvert]
   Int_t           StdHepPdg[24];   //[StdHepN]

   Int_t           nhits;
   Double_t        hit_time[1000000];   //[nhits]
   Double_t        hit_x[1000000];   //[nhits]
   Double_t        hit_y[1000000];   //[nhits]
   Double_t        hit_z[1000000];   //[nhits]
   Int_t           hit_PMTid[1000000];   //[nhits]
   Int_t           npart;
   Double_t        part_xStart[1000000];   //[npart]
   Double_t        part_yStart[1000000];   //[npart]
   Double_t        part_zStart[1000000];   //[npart]
   Double_t        part_tStart[1000000];   //[npart]
   Double_t        part_xEnd[1000000];   //[npart]
   Double_t        part_yEnd[1000000];   //[npart]
   Double_t        part_zEnd[1000000];   //[npart]
   Double_t        part_tEnd[1000000];   //[npart]
   Double_t        part_pxStart[1000000];   //[npart]
   Double_t        part_pyStart[1000000];   //[npart]
   Double_t        part_pzStart[1000000];   //[npart]
   Double_t        part_pxEnd[1000000];   //[npart]
   Double_t        part_pyEnd[1000000];   //[npart]
   Double_t        part_pzEnd[1000000];   //[npart]
   Double_t        part_KEstart[1000000];   //[npart]
   Double_t        part_KEend[1000000];   //[npart]
   Int_t           part_processStart[1000000];   //[npart]
   Int_t           part_processEnd[1000000];   //[npart]
   Int_t           part_parentid[1000000];   //[npart]
   Int_t           part_trackid[1000000];   //[npart]
   Int_t           part_pid[1000000];   //[npart]
   Int_t           ncapturecount;
   Int_t           neutroncount;
   Double_t        capt_x[500];   //[ncapturecount]
   Double_t        capt_y[500];   //[ncapturecount]
   Double_t        capt_z[500];   //[ncapturecount]
   Double_t        capt_t0[500];   //[ncapturecount]
   Double_t        capt_E[500];   //[ncapturecount]
   Int_t           capt_num[500];   //[ncapturecount]
   Int_t           capt_pid[500];   //[ncapturecount]
   Int_t           capt_nucleus[500];   //[ncapturecount]
   Int_t           capt_nphot[500];   //[ncapturecount]
   Int_t           capt_ngamma[500];   //[ncapturecount]
   Double_t        neutrino_E;
   Int_t           neutrino_id;
   Double_t        neutrino_px;
   Double_t        neutrino_py;
   Double_t        neutrino_pz;
   Int_t           ntrks;
   Int_t           nneutrons;
   Double_t        vtxx;
   Double_t        vtxy;
   Double_t        vtxz;
   Int_t           mpid[100];   //[ntrks]
   Double_t        px[100];   //[ntrks]
   Double_t        py[100];   //[ntrks]
   Double_t        pz[100];   //[ntrks]
   Double_t        KE[100];   //[ntrks]

   selections(TTree *tree=0, TTree * tracker=0, TTree * out=0);
   virtual ~selections();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, TTree *tracker, TTree * out);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

};
#endif

#ifdef selections_cxx
selections::selections(TTree *tree, TTree *tracker, TTree *out) : fChain(0), nRooTracker(0), outChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TChain *f=new TChain("Final_Reconstruction");
      TChain *le=new TChain("Low_E");
      TChain *hee=new TChain("High_E_Electron");
      TChain *hem=new TChain("High_E_Muon");
      TChain *d=new TChain("Debug");
      TChain *tr = new TChain("nRooTracker");
      TChain *o = new TChain("HitsTree");
      TString flavs[4] = {"numu","nue","antinumu","antinue"};
      TString hcs[4] = {"nu","antinu"};
      for(int k=0; k<1; k++){
         const char * h = hcs[k].Data();
         for (int j = 0; j < 4; j++) {
            const char *s = flavs[j].Data();
            for (int i = 1000; i < 1100; i++) {
               //Missing vector files for antinu mode:
               if (k == 1 && j == 0 && (i == 1096 || i==1037)) continue;
               if (k == 1 && j == 1 && (i == 1017 || i == 1053)) continue;
               if (k == 1 && j == 2 && i == 1078) continue;
               if (k == 1 && j == 3 && (i == 1018 || i == 1019 || i == 1053)) continue;
               char *file = Form("/data/hyperk/wchsandbox_reco/flav_%s/%s_%s_%i/%s_%s_%i_12in_reco_5.root", s, h, s, i, h, s,i);
               f->AddFile(file);
               le->AddFile(file);
               hee->AddFile(file);
               hem->AddFile(file);
               d->AddFile(file);
               tr->AddFile(Form("/data/hyperk/wchsandbox_reco/vectors/v00-01/flav_%s/genev_%s_cylinder_r551_z2200_Z_%s_1721827_%i.root",s, h, s, i), 1000);
               o->AddFile(Form("/data/hyperk/wchsandbox_reco/flav_%s/%s_%s_%i/%s_%s_%i_out_12in.root", s, h, s, i, h, s, i));
            }
         }
      }
      f->AddFriend(le);
      f->AddFriend(hem);
      f->AddFriend(hee);
      f->AddFriend(d);
      tree = f;
      tracker = tr;
      out = o;
   }
   tree->AddFriend(tracker);
   tree->AddFriend(out);
   Init(tree,tracker,out);
}

selections::~selections()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   if (!nRooTracker) return;
   delete nRooTracker->GetCurrentFile();
   if (!outChain) return;
   delete outChain->GetCurrentFile();
}

Int_t selections::GetEntry(Long64_t entry)
{
// Read contents of entry.
   int nb = 0;
   if (!fChain) return 0;
   nb += fChain->GetEntry(entry);
   if (!nRooTracker) return 0;
   nb += nRooTracker->GetEntry(entry);
   if (!outChain) return 0;
   nb += outChain->GetEntry(entry);
   return nb;
}
Long64_t selections::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   if (!nRooTracker) return -5;
   Long64_t centry2 = nRooTracker->LoadTree(entry);
   if (centry2 < 0) return centry2;
   if (nRooTracker->GetTreeNumber() != fCurrent2) {
      fCurrent2 = nRooTracker->GetTreeNumber();
      Notify();
   }
   if (!outChain) return -5;
   Long64_t centry3 = outChain->LoadTree(entry);
   if (centry3 < 0) return centry3;
   if (outChain->GetTreeNumber() != fCurrent3) {
      fCurrent3 = outChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void selections::Init(TTree *tree, TTree *tracker, TTree *out)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or fChain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   if (!tracker) return;
   if (!out) return;
   fChain = tree;
   nRooTracker = tracker;
   outChain = out;
   fCurrent = -1;
   fCurrent2 = -1;
   fCurrent3 = -1;
   fChain->SetMakeClass(1);
   nRooTracker->SetMakeClass(1);
   outChain->SetMakeClass(1);

   fChain->SetBranchAddress("evt", &evt);
   fChain->SetBranchAddress("nSubevents", &nSubEvts);
   fChain->SetBranchAddress("nClusters", &nClusters);
   fChain->SetBranchAddress("recoNRings", recoNRings);
   fChain->SetBranchAddress("cluster", subevt);
   fChain->SetBranchAddress("ring", ring);
   fChain->SetBranchAddress("isHighE", isHighE);
   fChain->SetBranchAddress("ringPEs", ringPEs);
   fChain->SetBranchAddress("trueVtxX", &trueVtxX);
   fChain->SetBranchAddress("trueVtxY", &trueVtxY);
   fChain->SetBranchAddress("trueVtxZ", &trueVtxZ);
   fChain->SetBranchAddress("trueTime", &trueTime);
   fChain->SetBranchAddress("trueDirX", &trueDirX);
   fChain->SetBranchAddress("trueDirY", &trueDirY);
   fChain->SetBranchAddress("trueDirZ", &trueDirZ);
   fChain->SetBranchAddress("trueKE", &trueKE);
   fChain->SetBranchAddress("diffKE", &diffKE);
   fChain->SetBranchAddress("diffVtxX", &diffVtxX);
   fChain->SetBranchAddress("diffVtxY", &diffVtxY);
   fChain->SetBranchAddress("diffVtxZ", &diffVtxZ);
   fChain->SetBranchAddress("diffTime", &diffTime);
   fChain->SetBranchAddress("diffDirX", &diffDirX);
   fChain->SetBranchAddress("diffDirY", &diffDirY);
   fChain->SetBranchAddress("diffDirZ", &diffDirZ);
   fChain->SetBranchAddress("diffVtxAbs", &diffVtxAbs);
   fChain->SetBranchAddress("diffDirAbs", &diffDirAbs);
   fChain->SetBranchAddress("mode", &mode);
   fChain->SetBranchAddress("neutrinoE", &neutrinoE);
   fChain->SetBranchAddress("neutrinoDirX", &neutrinoDirX);
   fChain->SetBranchAddress("neutrinoDirY", &neutrinoDirY);
   fChain->SetBranchAddress("neutrinoDirZ", &neutrinoDirZ);
   fChain->SetBranchAddress("neutrinoPID", &neutrinoPID);
   fChain->SetBranchAddress("nNeutrons", &nNeutrons);
   fChain->SetBranchAddress("neutronCount", &neutronCount);
   fChain->SetBranchAddress("nCaptures", &nCaptures);

//   fChain->SetBranchAddress("evt", &evt);
//   fChain->SetBranchAddress("nSubevents", &nSubevents);
//   fChain->SetBranchAddress("cluster", cluster);
//   fChain->SetBranchAddress("ring", ring);
   fChain->SetBranchAddress("recoVtxXLowE", recoVtxXLowE);
   fChain->SetBranchAddress("recoVtxYLowE", recoVtxYLowE);
   fChain->SetBranchAddress("recoVtxZLowE", recoVtxZLowE);
   fChain->SetBranchAddress("recoTimeLowE", recoTimeLowE);
   fChain->SetBranchAddress("recoDirXLowE", recoDirXLowE);
   fChain->SetBranchAddress("recoDirYLowE", recoDirYLowE);
   fChain->SetBranchAddress("recoDirZLowE", recoDirZLowE);
   fChain->SetBranchAddress("recoChkvAngleLowE", recoChkvAngleLowE);
   fChain->SetBranchAddress("recoEnergyLowE", recoEnergyLowE);

   fChain->SetBranchAddress("recoVtxXHighEMuon", recoVtxXHighEMuon);
   fChain->SetBranchAddress("recoVtxYHighEMuon", recoVtxYHighEMuon);
   fChain->SetBranchAddress("recoVtxZHighEMuon", recoVtxZHighEMuon);
   fChain->SetBranchAddress("recoTimeHighEMuon", recoTimeHighEMuon);
   fChain->SetBranchAddress("recoDirXHighEMuon", recoDirXHighEMuon);
   fChain->SetBranchAddress("recoDirYHighEMuon", recoDirYHighEMuon);
   fChain->SetBranchAddress("recoDirZHighEMuon", recoDirZHighEMuon);
   fChain->SetBranchAddress("recoChkvAngleHighEMuon", recoChkvAngleHighEMuon);
   fChain->SetBranchAddress("recoEnergyHighEMuon", recoEnergyHighEMuon);
   fChain->SetBranchAddress("recoLnLHighEMuon", recoLnLHighEMuon);

   fChain->SetBranchAddress("recoVtxXHighEElectron", recoVtxXHighEElectron);
   fChain->SetBranchAddress("recoVtxYHighEElectron", recoVtxYHighEElectron);
   fChain->SetBranchAddress("recoVtxZHighEElectron", recoVtxZHighEElectron);
   fChain->SetBranchAddress("recoTimeHighEElectron", recoTimeHighEElectron);
   fChain->SetBranchAddress("recoDirXHighEElectron", recoDirXHighEElectron);
   fChain->SetBranchAddress("recoDirYHighEElectron", recoDirYHighEElectron);
   fChain->SetBranchAddress("recoDirZHighEElectron", recoDirZHighEElectron);
   fChain->SetBranchAddress("recoChkvAngleHighEElectron", recoChkvAngleHighEElectron);
   fChain->SetBranchAddress("recoEnergyHighEElectron", recoEnergyHighEElectron);
   fChain->SetBranchAddress("recoLnLHighEElectron", recoLnLHighEElectron);

   fChain->SetBranchAddress("recoVtxX", recoVtxX);
   fChain->SetBranchAddress("recoVtxY", recoVtxY);
   fChain->SetBranchAddress("recoVtxZ", recoVtxZ);
   fChain->SetBranchAddress("recoTime", recoTime);
   fChain->SetBranchAddress("recoDirX", recoDirX);
   fChain->SetBranchAddress("recoDirY", recoDirY);
   fChain->SetBranchAddress("recoDirZ", recoDirZ);
   fChain->SetBranchAddress("recoChkvAngle", recoChkvAngle);
   fChain->SetBranchAddress("recoEnergy", recoEnergy);
   fChain->SetBranchAddress("recoPID", recoPID);
   fChain->SetBranchAddress("trueToWall", &trueToWall);
   fChain->SetBranchAddress("trueDWall", &trueDWall);
   fChain->SetBranchAddress("recoToWall", &recoToWall);
   fChain->SetBranchAddress("recoDWall", &recoDWall);
   fChain->SetBranchAddress("recoCaptures", &recoCaptures);



   nRooTracker->SetBranchAddress("NEneutmode", &NEneutmode);
   nRooTracker->SetBranchAddress("NEnvc", &NEnvc);
   nRooTracker->SetBranchAddress("NEipvc", NEipvc);
   nRooTracker->SetBranchAddress("NEpvc", NEpvc);
   nRooTracker->SetBranchAddress("NEiorgvc", NEiorgvc);
   nRooTracker->SetBranchAddress("NEiflgvc", NEiflgvc);
   nRooTracker->SetBranchAddress("NEicrnvc", NEicrnvc);
   nRooTracker->SetBranchAddress("NEcrsx", &NEcrsx);
   nRooTracker->SetBranchAddress("NEcrsy", &NEcrsy);
   nRooTracker->SetBranchAddress("NEcrsz", &NEcrsz);
   nRooTracker->SetBranchAddress("NEcrsphi", &NEcrsphi);
   nRooTracker->SetBranchAddress("NEnvert", &NEnvert);
   nRooTracker->SetBranchAddress("NEposvert", NEposvert);
   nRooTracker->SetBranchAddress("NEiflgvert", NEiflgvert);
   nRooTracker->SetBranchAddress("NEnvcvert", &NEnvcvert);
   nRooTracker->SetBranchAddress("NEdirvert", NEdirvert);
   nRooTracker->SetBranchAddress("NEabspvert", NEabspvert);
   nRooTracker->SetBranchAddress("NEabstpvert", NEabstpvert);
   nRooTracker->SetBranchAddress("NEipvert", NEipvert);
   nRooTracker->SetBranchAddress("NEiverti", NEiverti);
   nRooTracker->SetBranchAddress("NEivertf", NEivertf);
   nRooTracker->SetBranchAddress("StdHepPdg", StdHepPdg);

   //outChain->SetBranchAddress("evt", &evt);
   outChain->SetBranchAddress("nhits", &nhits);
   outChain->SetBranchAddress("hit_time", hit_time);
   outChain->SetBranchAddress("hit_x", hit_x);
   outChain->SetBranchAddress("hit_y", hit_y);
   outChain->SetBranchAddress("hit_z", hit_z);
   //outChain->SetBranchAddress("hit_PMTid", hit_PMTid);
   outChain->SetBranchAddress("npart", &npart);
   outChain->SetBranchAddress("part_xStart", &part_xStart);
   outChain->SetBranchAddress("part_yStart", &part_yStart);
   outChain->SetBranchAddress("part_zStart", &part_zStart);
   outChain->SetBranchAddress("part_tStart", &part_tStart);
   outChain->SetBranchAddress("part_xEnd", &part_xEnd);
   outChain->SetBranchAddress("part_yEnd", &part_yEnd);
   outChain->SetBranchAddress("part_zEnd", &part_zEnd);
   outChain->SetBranchAddress("part_tEnd", &part_tEnd);
   outChain->SetBranchAddress("part_pxStart", &part_pxStart);
   outChain->SetBranchAddress("part_pyStart", &part_pyStart);
   outChain->SetBranchAddress("part_pzStart", &part_pzStart);
   outChain->SetBranchAddress("part_pxEnd", &part_pxEnd);
   outChain->SetBranchAddress("part_pyEnd", &part_pyEnd);
   outChain->SetBranchAddress("part_pzEnd", &part_pzEnd);
   outChain->SetBranchAddress("part_KEstart", &part_KEstart);
   outChain->SetBranchAddress("part_KEend", &part_KEend);
   outChain->SetBranchAddress("part_processStart", &part_processStart);
   outChain->SetBranchAddress("part_processEnd", &part_processEnd);
   outChain->SetBranchAddress("part_parentid", &part_parentid);
   outChain->SetBranchAddress("part_trackid", &part_trackid);
   outChain->SetBranchAddress("part_pid", &part_pid);
   outChain->SetBranchAddress("ncapturecount", &ncapturecount);
   outChain->SetBranchAddress("neutroncount", &neutroncount);
   outChain->SetBranchAddress("capt_x", capt_x);
   outChain->SetBranchAddress("capt_y", capt_y);
   outChain->SetBranchAddress("capt_z", capt_z);
   outChain->SetBranchAddress("capt_t0", capt_t0);
   outChain->SetBranchAddress("capt_E", capt_E);
   outChain->SetBranchAddress("capt_num", capt_num);
   outChain->SetBranchAddress("capt_pid", capt_pid);
   outChain->SetBranchAddress("capt_nucleus", capt_nucleus);
   outChain->SetBranchAddress("capt_nphot", capt_nphot);
   outChain->SetBranchAddress("capt_ngamma", capt_ngamma);
   outChain->SetBranchAddress("mode", &mode);
   outChain->SetBranchAddress("neutrino_E", &neutrino_E);
   outChain->SetBranchAddress("neutrino_id", &neutrino_id);
   outChain->SetBranchAddress("neutrino_px", &neutrino_px);
   outChain->SetBranchAddress("neutrino_py", &neutrino_py);
   outChain->SetBranchAddress("neutrino_pz", &neutrino_pz);
   outChain->SetBranchAddress("ntrks", &ntrks);
   outChain->SetBranchAddress("nneutrons", &nneutrons);
   outChain->SetBranchAddress("vtxx", &vtxx);
   outChain->SetBranchAddress("vtxy", &vtxy);
   outChain->SetBranchAddress("vtxz", &vtxz);
   outChain->SetBranchAddress("mpid", mpid);
   outChain->SetBranchAddress("px", px);
   outChain->SetBranchAddress("py", py);
   outChain->SetBranchAddress("pz", pz);
   outChain->SetBranchAddress("KE", KE);
   
   Notify();
}

Bool_t selections::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void selections::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
   if (!nRooTracker) return;
   nRooTracker->Show(entry);
   if (!outChain) return;
   outChain->Show(entry);
}
Int_t selections::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef selections_cxx
