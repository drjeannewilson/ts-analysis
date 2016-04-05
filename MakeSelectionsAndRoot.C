//
// Makes the selection files in the format that Raj uses for his sensitivity studies and the root tree used by Jeanne at the same time
//
// Run in root:
//   .L MakeSelectionsAndRoot.C
//   MakeSelectionsAndRoot s(IsAntiNu);
//   s.Loop();
// or s.Loop(0); to omit outputting the full selection ntuple (big)
// IsAntiNu should be false for forward horn current (nu beam) or true for reverse horn current (antinu beam)
//

// Jeanne's modifications: 
// Added neutron counting
// Added decay electron cuts
// Added some ouput variables to the output ntuples
// modified the selection cuts (PID etc.) in line with tuning for optimised selection
// Added a weighting on 1/9 to NCpi0 events (just by neut code) justify that fitQun type algorithm gives 
// this level of improvement

#define MakeSelectionsAndRoot_cxx
#include "MakeSelectionsAndRoot.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>
#include <TRandom3.h>

void MakeSelectionsAndRoot::Loop(bool outputntuple)
{

   if (fChain == 0) return;

   float Abspvc[100];
   //int mode;
   int Numatom=1;
   double Erec, Etrue;
   int recoNeutrons;
   int ntrueNcap;
   int ntrueNeutrons;
   int ncapGd;
   int ncapD;
   int ncapOther;

   Float_t         pnu[100];   //[numnu]  GeV/c
   Float_t         dirnu[100][3];   //[numnu]
   Int_t Numbndn;
   Int_t Numbndp;
   Int_t Numfrep;
   Int_t Ibound;

   const char * s = IsAntiNu ? "RHC" : "FHC";
   TFile *out = new TFile(Form("selections_androot_tagged_%s.root",s),"RECREATE");
   //cout << "flav count " << flav_count[0] << " " << flav_count[1] << " " << flav_count[2] << " " << flav_count[3] << endl;
   TTree *numu_mulike_nontag = MakeNtuple("numu_mulike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[0],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *nue_mulike_nontag = MakeNtuple("nue_mulike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[1],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *numubar_mulike_nontag = MakeNtuple("numubar_mulike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[2],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *nuebar_mulike_nontag = MakeNtuple("nuebar_mulike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[3],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *numu_elike_nontag = MakeNtuple("numu_elike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[0],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *nue_elike_nontag = MakeNtuple("nue_elike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[1],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *numubar_elike_nontag = MakeNtuple("numubar_elike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[2],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *nuebar_elike_nontag = MakeNtuple("nuebar_elike_nontag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[3],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *numu_mulike_ntag = MakeNtuple("numu_mulike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[0],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *nue_mulike_ntag = MakeNtuple("nue_mulike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[1],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *numubar_mulike_ntag = MakeNtuple("numubar_mulike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[2],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *nuebar_mulike_ntag = MakeNtuple("nuebar_mulike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[3],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *numu_elike_ntag = MakeNtuple("numu_elike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[0],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *nue_elike_ntag = MakeNtuple("nue_elike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[1],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *numubar_elike_ntag = MakeNtuple("numubar_elike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[2],recoNeutrons,ntrueNcap, ntrueNeutrons);
   TTree *nuebar_elike_ntag = MakeNtuple("nuebar_elike_ntag", Abspvc, Numatom, Erec, Etrue, pnu, dirnu, Numbndn, Numbndp, Numfrep, Ibound, flav_count[3],recoNeutrons,ntrueNcap, ntrueNeutrons);

   // Output a root file of all events
   TTree * ntuple = new TTree("AllEvents","AllEvents");
   double lepton_KE;
   double recoPIDLikelihood;
   int npi0;
   int nneutron;
   int ring1PEs;
   int ring2PEs;
   bool isHighE0;
   double mu_reco_dwall;
   double mu_reco_towall;
   double mu_reco_nu_E;
   double mu_recoKE;
   double e_reco_dwall;
   double e_reco_towall;
   double e_reco_nu_E;
   double e_recoKE;
   int recoMichels;
  
   // paired down output ntuple using parameters for developing selections on (may want to add in more parameters from selections.C)
   ntuple->Branch("neutrino_id", &neutrino_id, "neutrino_id/I");
   ntuple->Branch("neutrino_E", &neutrino_E, "neutrino_E/D");
   ntuple->Branch("trueKE", &lepton_KE, "trueKE/D");
   ntuple->Branch("interaction_mode", &mode, "interaction_mode/I");
   ntuple->Branch("n_pi0", &npi0, "n_pi0/I");
   ntuple->Branch("nneutron", &nneutron, "nneutron/I");
   ntuple->Branch("ncapture", &ntrueNcap, "ncapture/I");
   ntuple->Branch("ncapGd", &ncapGd, "ncapGd/I");
   ntuple->Branch("ncapD", &ncapD, "ncapD/I");
   ntuple->Branch("ncapOther", &ncapOther, "ncapOther/I");
   ntuple->Branch("true_dwall", &trueDWall, "true_dwall/D");
   ntuple->Branch("true_towall", &trueToWall, "true_towall/D");
   ntuple->Branch("mu_reco_dwall", &mu_reco_dwall, "mu_reco_dwall/D");
   ntuple->Branch("mu_reco_towall", &mu_reco_towall, "mu_reco_towall/D");
   ntuple->Branch("mu_reco_nu_E", &mu_reco_nu_E, "mu_reco_nu_E/D");
   ntuple->Branch("mu_recoKE", &mu_recoKE, "mu_recoKE/D");
   ntuple->Branch("e_reco_dwall", &e_reco_dwall, "e_reco_dwall/D");
   ntuple->Branch("e_reco_towall", &e_reco_towall, "e_reco_towall/D");
   ntuple->Branch("e_reco_nu_E", &e_reco_nu_E, "e_reco_nu_E/D");
   ntuple->Branch("e_recoKE", &e_recoKE, "e_recoKE/D");
   ntuple->Branch("nClusters",&nClusters,"nClusters/I");
   ntuple->Branch("recoNRings",&recoNRings, "recoNRings[nClusters]/I");
   ntuple->Branch("ring1PEs", &ring1PEs, "ring1PEs/I");
   ntuple->Branch("ring2PEs", &ring2PEs, "ring2PEs/I");
   ntuple->Branch("recoPIDLikelihood", &recoPIDLikelihood, "recoPIDLikelihood/D");
   ntuple->Branch("isHighE", &isHighE0, "isHighE/O");
   ntuple->Branch("nSubevents",&nSubevents,"nSubevents/I");
   ntuple->Branch("recoTime",recoTime,"recoTime[nSubevents]/D");
   ntuple->Branch("recoEnergy",recoEnergy,"recoEnergy[nSubevents]/D");
   ntuple->Branch("recoCaptures",&recoCaptures,"recoCaptures/I");
   ntuple->Branch("recoMichels",&recoMichels,"recoMichel/I");

   Long64_t nentries = fChain->GetEntriesFast();

   TRandom3 rand(31415);
   TRandom3 myrand(987654);  

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = GetEntry(jentry);   nbytes += nb;
      // Need to extract parameters and fill ntuple before selections
      if (nClusters == 0 || !isHighE[0]) {
        ring1PEs = 0;
        ring2PEs = 0;
        isHighE0 = false;
        recoPIDLikelihood = 0;
        mu_reco_dwall = 0;
        mu_reco_towall = 0;  //
        mu_reco_nu_E = 0;
        mu_recoKE = 0;
        e_reco_dwall = 0;
        e_reco_towall = 0;
        e_reco_nu_E = 0;
        e_recoKE = 0;
      }else{
        ring1PEs = ringPEs[0];
        ring2PEs = (recoNRings[0] > 1 ? ringPEs[1] : 0);
        isHighE0 = isHighE[0];
        recoPIDLikelihood = recoLnLHighEMuon[0] - recoLnLHighEElectron[0];
     
        // parameters for reco (dwall, towall, KE, nu_e) depend on hypothesis
        // mu hypothesis --------------
        double vtxX = recoVtxXHighEMuon[0];
        double vtxY = recoVtxYHighEMuon[0];
        double vtxZ = recoVtxZHighEMuon[0];
        double recoVtxR2 = vtxX * vtxX + vtxY * vtxY;
        double dWallZ_mu = 1100-TMath::Abs(vtxZ);
        double dWallR_mu = 550-TMath::Sqrt(recoVtxR2);
        // take the smallest dwall
        mu_reco_dwall = dWallR_mu < dWallZ_mu ? dWallR_mu : dWallZ_mu;
        
        double mn = 939;
        double ml = 105;
        mu_recoKE= recoEnergyHighEMuon[0];
        double recoE = mu_recoKE + ml;
        double recoP = TMath::Sqrt(recoE*recoE- ml * ml);
        double dirZ = recoDirZHighEMuon[0];
        
        mu_reco_nu_E = (mn*recoE- ml * ml /2.)/(mn-recoE+recoP* dirZ);
        
        double a = 1- dirZ * dirZ;
        double dirX = recoDirXHighEMuon[0];
        double dirY = recoDirYHighEMuon[0];
        double b = vtxX * dirX + vtxY * dirY;
        double c = recoVtxR2-550*550;
        double recoToWallR = (TMath::Sqrt(b*b-a*c)-b)/a;
        double recoToWallZ = 1100 - vtxZ *TMath::Abs(dirZ);
        // take smallest towall
        mu_reco_towall = recoToWallR<recoToWallZ ? recoToWallR : recoToWallZ;
        
        // e hypothesis --------------
        vtxX = recoVtxXHighEElectron[0];
        vtxY = recoVtxYHighEElectron[0];
        vtxZ = recoVtxZHighEElectron[0];
        recoVtxR2 = vtxX * vtxX + vtxY * vtxY;
        double dWallZ_e = 1100-TMath::Abs(vtxZ);
        double dWallR_e = 550-TMath::Sqrt(recoVtxR2);
        // take the smallest dwall
        e_reco_dwall = dWallR_e < dWallZ_e ? dWallR_e : dWallZ_e;
        
        mn = 939;
        ml = 0.511;
        e_recoKE = recoEnergyHighEElectron[0];
        recoE = e_recoKE + ml;
        recoP = TMath::Sqrt(recoE*recoE- ml * ml);
        dirZ = recoDirZHighEElectron[0];
        
        e_reco_nu_E = (mn*recoE- ml * ml /2.)/(mn-recoE+recoP* dirZ);
        
        a = 1- dirZ * dirZ;
        dirX = recoDirXHighEElectron[0];
        dirY = recoDirYHighEElectron[0];
        b = vtxX * dirX + vtxY * dirY;
        c = recoVtxR2-550*550;
        recoToWallR = (TMath::Sqrt(b*b-a*c)-b)/a;
        recoToWallZ = 1100 - vtxZ *TMath::Abs(dirZ);
        // take smallest towall
        e_reco_towall = recoToWallR<recoToWallZ ? recoToWallR : recoToWallZ;
      }
      //  --------------
     
      // New method to count up the pi0 in the final state -> weight down events with pi0 as well as by neut code (see below)
      // not sure how to count pi0 so use 3 checks:
      int npi0_track = 0;
      int nneutron_track = 0;
      lepton_KE = 0;  // also store the lepton KE
      for (int itrk = 0; itrk < ntrks; itrk++) {
         // Find highest KE outgoing charged lepton
         if (TMath::Abs(mpid[itrk]) == 11 || (TMath::Abs(mpid[itrk]) == 13 && KE[itrk] > lepton_KE)) lepton_KE = KE[itrk];
         if (mpid[itrk] == 111) npi0_track++;
         if (mpid[itrk] == 2112) nneutron_track++;
      }
      // loop through SI
      int npi0_SI = 0;
      int nneutron_SI = 0;
      for(int ipart = 0; ipart < npart; ++ipart){
         if (part_pid[ipart] == 111) npi0_SI++;
         if (part_pid[ipart] == 2112) nneutron_SI++;
      }
      // loop through FI too
      int npi0_FS = 0;
      int nneutron_FS = 0;
       for(int ipart = 0; ipart< NEnvc; ++ipart){
         if (NEipvc[ipart] == 111) npi0_FS++;
         if (NEipvc[ipart] == 2112) nneutron_FS++;
        
      } 
      // if((npi0_track+npi0_SI+npi0_FS)>0)std::cout << " pi0: " << npi0_track << ", " << npi0_SI << " " << npi0_FS << std::endl;
      // need to decide which of these to use but use all pi0s present for the NC weight to be sure
      npi0 = npi0_FS;
      nneutron = nneutron_FS;
     
      Etrue = TMath::Sqrt(NEpvc[0][0]*NEpvc[0][0]+NEpvc[0][1]*NEpvc[0][1]+NEpvc[0][2]*NEpvc[0][2]);
      for(int i=0; i<NEnvc; i++) {
         Abspvc[i] = TMath::Sqrt(NEpvc[i][0]*NEpvc[i][0] + NEpvc[i][1]*NEpvc[i][1] + NEpvc[i][2]*NEpvc[i][2]);
         pnu[i] = Abspvc[i] / 1000.;
         for (int j = 0; j < 3; j++) {
            dirnu[i][j] = NEpvc[i][j] / Abspvc[i];
         }
      }
      if(StdHepPdg[1]==1000080160){ //Oxygen target
         Numatom = 16;
         Numbndn = 8;
         Numbndp = 8;
         Numfrep = 0;
         Ibound = 1;
      }
      else if(StdHepPdg[1]==1000010010){ //Hydrogen target
         Numatom = 1;
         Numbndn = 0;
         Numbndp = 0;
         Numfrep = 1;
         Ibound = 0;
      }
      else{
         Numatom = 0;
         Numbndn = 0;
         Numbndp = 0;
         Numfrep = 0;
         Ibound = 0;
      }

      bool neutron_tag = false;
      ntrueNcap = nCaptures;  // Get from original file
      ncapGd = 0;
      ncapD = 0;
      ncapOther = 0;
      for(int ic = 0 ; ic<nCaptures; ++ic){
      	// Check the capture nucleus
	if(capt_nucleus[ic]>=3336){
	  ncapGd++;
	}else if(capt_nucleus[ic]==3329){
	  ncapD++;
        }else{
	  ncapOther++;
	}
        //cout << ic << " " << capt_nucleus[ic] << ", ";
      }
      //cout << endl;
      ntrueNeutrons = nneutron;  // Counted from simulation FS info above
      recoMichels = 0;
      recoNeutrons = 0;
      for(int iSubevent = 0; iSubevent < nSubevents; iSubevent ++) {
            if(recoTime[iSubevent] > 135 && recoTime[iSubevent]<8000 && recoEnergy[iSubevent]>15.){
              recoMichels++;
            }
            if (recoEnergy[iSubevent] > 2 && recoEnergy[iSubevent] < 10 && recoTime[iSubevent] > 1000 && recoTime[iSubevent] < 100000) {
              neutron_tag = true;
              recoNeutrons++;
           }
      }
      // put all events in the ntuple
      if(outputntuple)ntuple->Fill();
     
      // Weight down the NCpi0 here:
      float myran = myrand.Rndm();
      if(((abs(NEneutmode)==31 || abs(NEneutmode)==32 || abs(NEneutmode)==36) || ((npi0_FS+npi0_SI+npi0_track)>0))&& myran<0.88889)continue;
     
// dwall
// towall
// 1ring
// PID
// michels
// energy
// NCpiweight
// N followers

      // Now we choose what to put in each tree for Valor
      bool eLike = e_reco_dwall>200 && e_reco_towall>200 && (isHighE0 && recoNRings[0]==1) && ((recoLnLHighEMuon[0] - recoLnLHighEElectron[0]) <-400) && recoMichels==0 && e_reco_nu_E>0&&e_reco_nu_E<2500&&e_recoKE>100&&e_recoKE<2500 ;

      bool muLike =  mu_reco_dwall>200 && mu_reco_towall>200 && (isHighE0 && recoNRings[0]==1) && ((recoLnLHighEMuon[0] - recoLnLHighEElectron[0]) > 0)&& recoMichels<=1  && mu_reco_nu_E>0&&mu_reco_nu_E<1250&&mu_recoKE>200&&mu_recoKE<2000;

      if(eLike && muLike) std::cout << " Both???" << std::endl; // This shouldn't happen!
//      if(!eLike && !muLike) std::cout << " Neither???" << std::endl;  // This can happen
     
      if(neutron_tag) {
         if (muLike) {
             Erec = mu_reco_nu_E;
            if (neutrino_id == 14)
               numu_mulike_ntag->Fill();
            else if (neutrino_id == -14)
               numubar_mulike_ntag->Fill();
            else if (neutrino_id == 12)
               nue_mulike_ntag->Fill();
            else if (neutrino_id == -12)
               nuebar_mulike_ntag->Fill();
         }
         /*else */if (eLike) {
            Erec = e_reco_nu_E;
            if (neutrino_id == 14)
               numu_elike_ntag->Fill();
            else if (neutrino_id == -14)
               numubar_elike_ntag->Fill();
            else if (neutrino_id == 12)
               nue_elike_ntag->Fill();
            else if (neutrino_id == -12)
               nuebar_elike_ntag->Fill();
         }
      }
      else {
         if (muLike) {
            Erec = mu_reco_nu_E;
            if (neutrino_id == 14)
               numu_mulike_nontag->Fill();
            else if (neutrino_id == -14)
               numubar_mulike_nontag->Fill();
            else if (neutrino_id == 12)
               nue_mulike_nontag->Fill();
            else if (neutrino_id == -12)
               nuebar_mulike_nontag->Fill();
         }
         /*else */if (eLike) {
            Erec = e_reco_nu_E;
            if (neutrino_id == 14)
               numu_elike_nontag->Fill();
            else if (neutrino_id == -14)
               numubar_elike_nontag->Fill();
            else if (neutrino_id == 12)
               nue_elike_nontag->Fill();
            else if (neutrino_id == -12)
               nuebar_elike_nontag->Fill();
         }
      }
   }// loop through entries, jentry
   numu_mulike_ntag->Write();
   numubar_mulike_ntag->Write();
   nue_mulike_ntag->Write();
   nuebar_mulike_ntag->Write();
   numu_elike_ntag->Write();
   numubar_elike_ntag->Write();
   nue_elike_ntag->Write();
   nuebar_elike_ntag->Write();
   numu_mulike_nontag->Write();
   numubar_mulike_nontag->Write();
   nue_mulike_nontag->Write();
   nuebar_mulike_nontag->Write();
   numu_elike_nontag->Write();
   numubar_elike_nontag->Write();
   nue_elike_nontag->Write();
   nuebar_elike_nontag->Write();
   if(outputntuple)ntuple->Write();
}

TTree *MakeSelectionsAndRoot::MakeNtuple(TString name, float Abspvc[100], int &Numatom, double &Erec, double &Etrue, Float_t pnu[100],
                                  Float_t dirnu[100][3], Int_t &Numbndn, Int_t &Numbndp, Int_t &Numfrep, Int_t &Ibound, Int_t &flavEvents, Int_t &nRecoNeutrons, Int_t &ntrueNcap, Int_t &ntrueNeutrons) {
   TTree * ntuple = new TTree(name, name);
   ntuple->Branch("evt",&evt,"evt/I");

   ntuple->Branch("Npvc",&NEnvc, "Npvc/I");
   ntuple->Branch("Ichvc", NEicrnvc, "Ichvc[Npvc]/I");
   ntuple->Branch("Ipvc", NEipvc, "Ipvc[Npvc]/I");
   ntuple->Branch("Abspvc",Abspvc, "Abspvc[Npvc]/F");
   ntuple->Branch("Pvc", NEpvc, "Pvc[Npvc][3]/F");
   ntuple->Branch("mode",&NEneutmode, "mode/I");

   ntuple->Branch("Crsx",&NEcrsx, "Crsx/F");
   ntuple->Branch("Crsy",&NEcrsy, "Crsy/F");
   ntuple->Branch("Crsz",&NEcrsz, "Crsz/F");
   ntuple->Branch("Crsphi",&NEcrsphi, "Crsphi/F");
   ntuple->Branch("Nvert",&NEnvert, "Nvert/I");
   ntuple->Branch("Iflgvert", NEiflgvert, "Iflgvert[Nvert]/I");
   ntuple->Branch("posvert", NEposvert, "posvert[Nvert][3]/F");
   ntuple->Branch("Nvcvert",&NEnvcvert, "Nvcvert/I");
   ntuple->Branch("Abspvert", NEabspvert, "Abspvert[Nvcvert]/F");
   ntuple->Branch("Abstpvert", NEabstpvert, "Abstpvert[Nvcvert]/F");
   ntuple->Branch("Ipvert", NEipvert, "Ipvert[Nvcvert]/I");
   ntuple->Branch("Iverti", NEiverti, "Iverti[Nvcvert]/I");
   ntuple->Branch("Ivertf", NEivertf, "Ivertf[Nvcvert]/I");
   ntuple->Branch("Dirvert", NEdirvert, "Dirvert[Nvcvert][3]/F");
   ntuple->Branch("numnu",&NEnvc, "numnu/I");
   ntuple->Branch("ipnu", NEipvc, "ipnu[numnu]/I");
   ntuple->Branch("pnu",pnu, "pnu[numnu][3]/F");
   ntuple->Branch("dirnu",dirnu, "dirnu[numnu][3]/F");
   ntuple->Branch("Numbndn",&Numbndn, "Numbndn/I");
   ntuple->Branch("Numbndp",&Numbndp, "Numbndp/I");
   ntuple->Branch("Numfrep",&Numfrep, "Numfrep/I");
   ntuple->Branch("Ibound",&Ibound, "Ibound/I");


   ntuple->Branch("Numatom",&Numatom, "Numatom/I");
   ntuple->Branch("Erec",&Erec,"Erec/D");
   ntuple->Branch("Etrue",&Etrue,"Etrue/D");

   ntuple->Branch("totFlavEvts",&flavEvents,"totFlavEvts/I");
   
   ntuple->Branch("nrecoNeutrons",&nRecoNeutrons,"nrecoNeutrons/I");
   ntuple->Branch("ntrueNcap",&ntrueNcap,"ntrueNcap/I");
   ntuple->Branch("ntrueNeutrons",&ntrueNeutrons,"ntrueNeutrons/I");
   
   ntuple->Branch("NEiorgvc", NEiorgvc, "NEiorgvc[Npvc]/I");
   ntuple->Branch("NEiflgvc", NEiflgvc, "NEiflgvc[Npvc]/I");
   ntuple->Branch("NEicrnvc", NEicrnvc, "NEicrnvc[Npvc]/I");
   ntuple->Branch("StdHepPdg", StdHepPdg, "StdHepPdg[24]/I");
   return ntuple;
}
