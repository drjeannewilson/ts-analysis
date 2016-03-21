//
// This script takes reconstruction output and puts it in a more useful format.
//
// To run in root:
//   .L selections.C+
//   selections s;
//   s.Loop();
//
// Branches in the output: (unless otherwise noted, reconstructed quantities are for the first ring found for the first cluster of hits)
//   evt                -  The number of the event in it's original file (will be repeated when combining multiple reconstruction files
//   neutrino_id        -  true neutrino flavour ID
//   neutrino_E         -  true neutrino energy
//   trueKE             -  true kinetic energy of charged lepton for CC events
//   interaction_mode   -  NEUT interaction mode
//   n_chkv_part        -  count of the true number of cherenkov producing particles
//   n_pi0              -  true number of neutral pions
//   n_e                -  true number of electrons
//   n_mu               -  true number of muons
//   n_pi_pm            -  true number of charged pions
//   n_k0               -  true number of neutral kaons
//   n_k_pm             -  true number of charged kaons
//   nneutron           -  true number of neutrons
//   nproton            -  true number of neutrons
//   ngamma             -  true number of gammas
//   nother             -  true number of particles not already accounted for
//   neutroncount       -  number of neutrons according to WChSandBox (sometimes wrong?)
//   n_captures         -  true number of neutron capture events according to WChSandBox
//   true_dwall         -  true distance between vertex and tank wall
//   true_towall        -  true distance from vertex to tank wall in direction of charged lepton
//   mu_reco_dwall      -  reconstructed distance from wall under muon hypothesis
//   mu_reco_towall     -  reconstructed distance to wall in track direction under muon hypothesis
//   mu_reco_nu_E       -  reconstructed neutrino energy (using CCQE formula) under muon hypothesis
//   mu_recoKE          -  reconstructed muon kinetic energy under muon hypothesis
//   mu_diffVtx         -  difference between reconstructed and true vertex under muon hypothesis
//   mu_diffDir         -  difference between reconstructed and true direction under muon hypothesis
//   mu_diffKE          -  difference between reconstructed and true muon KE under muon hypothesis
//   mu_diffVtxLong     -  difference between reconstructed and true vertex in direction of track under muon hypothesis
//   mu_diffVtxTrans    -  difference between reconstructed and true vertex perpundicular to track direction muon hypothesis
//   mu_diff_nu_E       -  difference between reconstructed and true neutrino energy under muon hypothesis
//   e_reco_dwall       -  reconstructed distance from wall under electron hypothesis
//   e_reco_towall      -  reconstructed distance to wall in track direction under electron hypothesis
//   e_reco_nu_E        -  reconstructed neutrino energy (using CCQE formula) under electron hypothesis
//   e_recoKE           -  reconstructed muon kinetic energy under electron hypothesis
//   e_diffVtx          -  difference between reconstructed and true vertex under electron hypothesis
//   e_diffDir          -  difference between reconstructed and true direction under electron hypothesis
//   e_diffKE           -  difference between reconstructed and true muon KE under electron hypothesis
//   e_diffVtxLong      -  difference between reconstructed and true vertex in direction of track under electron hypothesis
//   e_diffVtxTrans     -  difference between reconstructed and true vertex perpundicular to track direction electron hypothesis
//   e_diff_nu_E        -  difference between reconstructed and true neutrino energy under electron hypothesis
//   nClusters          -  number of clusters of hits found in reconstruction
//   ring1PEs           -  number of PEs assigned to the first ring for the first cluster of hits
//   ring2PEs           -  number of PEs assigned to the second ring for the first cluster of hits
//   recoPID            -  reconstructed charged lepton type
//   recoPIDLikelihood  -  difference in reconstructed likelihoods for muon minus that for electron
//   isHighE            -  whether reconstruction goes through high-E or stops at low-E
//   nSubevents         -  number of subevents (each ring of each cluster of hits counts as one subevent)
//   recoTime           -  reconstructed interaction time - array length nSubevents with element for each subevent
//   recoEnergy         -  reconstructed energy - array length nSubevents with element for each subevent
//   recoCaptures       -  reconstructed number of neutron captures, where one capture is one subevent with 2-10 MeV energy occurring after 1-100 microseconds

// Added by Jeanne
//   recoMichels - reconstructed number of decay electrons, where one decay electron is one subevent with E>15MeV and occuring after 135-8000 ns 
//-  recoNRings         - number of rings in each cluster
//-   E_pi_pm           - energy of any pi_pluses or pi_minuses created (n_pi_pm)
//-   isHighE_cluster    - did recon for this cluster go through highE chain?
//    mu_diffVtx_x
//    mu_diffVtx_y
//    mu_diffVtx_z
//    e_diffVtx_x
//    e_diffVtx_y
//    e_diffVtx_z


#define selections_cxx
#include "selections.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>

void selections::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L selections.C
//      Root > selections t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the fChain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;


   int nChkvPart;
//   double recoDWall;
//   double recoToWall;
//   double trueDWall;
//   double trueToWall;
   double recoNuE;
   double diffNuE;
   double diffVtxLong;
   double diffVtxTrans;
   double lepton_KE;
   double recoPIDLikelihood;
   double E_pi_pm[10];
   int ne;
   int nmu;
   int npi0;
   int npipm;
   int nk0;
   int nkpm;
   int nneutron;
   int nproton;
   int ngamma;
   int nother;
   int ring1PEs;
   int ring2PEs;
   bool isHighE0;
   int recoPID0;
   double mu_reco_dwall;
   double mu_reco_towall;
   double mu_reco_nu_E;
   double mu_recoKE;
   double mu_diffVtx;
   double mu_diffVtx_x;
   double mu_diffVtx_y;
   double mu_diffVtx_z;
   double e_diffVtx_x;
   double e_diffVtx_y;
   double e_diffVtx_z;
   double mu_diffDir;
   double mu_diffKE;
   double mu_diffVtxLong;
   double mu_diffVtxTrans;
   double mu_diff_nu_E;
   double e_reco_dwall;
   double e_reco_towall;
   double e_reco_nu_E;
   double e_recoKE;
   double e_diffVtx;
   double e_diffDir;
   double e_diffKE;
   double e_diffVtxLong;
   double e_diffVtxTrans;
   double e_diff_nu_E;
   int recoMichels;
   TFile *out = new TFile(Form("selection_test_%i_%i.root",TDatime().GetDate(),TDatime().GetTime()),"RECREATE");

   TTree * ntuple = new TTree("selection","selection");

   ntuple->Branch("evt",&evt,"evt/I");
   ntuple->Branch("neutrino_id", &neutrino_id, "neutrino_id/I");
   ntuple->Branch("neutrino_E", &neutrino_E, "neutrino_E/D");
   ntuple->Branch("trueKE", &lepton_KE, "trueKE/D");
   ntuple->Branch("interaction_mode", &mode, "interaction_mode/I");
   ntuple->Branch("n_chkv_part", &nChkvPart, "n_chkv_part/I");
   ntuple->Branch("n_pi0", &npi0, "n_pi0/I");
   ntuple->Branch("n_e", &ne, "n_e/I");
   ntuple->Branch("n_mu", &nmu, "n_mu/I");
   ntuple->Branch("n_pi_pm", &npipm, "n_pi_pm/I");
   ntuple->Branch("n_k0", &nk0, "n_k0/I");
   ntuple->Branch("n_k_pm", &nkpm, "n_k_pm/I");
   ntuple->Branch("E_pi_pm", &E_pi_pm, "E_pi_pm[n_pi_pm]/D");
   ntuple->Branch("nneutron", &nneutron, "nneutron/I");
   ntuple->Branch("nproton", &nproton, "nproton/I");
   ntuple->Branch("ngamma", &ngamma, "ngamma/I");
   ntuple->Branch("nother", &nother, "nother/I");
   ntuple->Branch("neutroncount", &neutroncount, "neutroncount/I");
   ntuple->Branch("n_captures", &ncapturecount, "n_captures/I");
   ntuple->Branch("true_dwall", &trueDWall, "true_dwall/D");
   ntuple->Branch("true_towall", &trueToWall, "true_towall/D");
   ntuple->Branch("mu_reco_dwall", &mu_reco_dwall, "mu_reco_dwall/D");
   ntuple->Branch("mu_reco_towall", &mu_reco_towall, "mu_reco_towall/D");
   ntuple->Branch("mu_reco_nu_E", &mu_reco_nu_E, "mu_reco_nu_E/D");
   ntuple->Branch("mu_recoKE", &mu_recoKE, "mu_recoKE/D");
   ntuple->Branch("mu_diffVtx", &mu_diffVtx, "mu_diffVtx/D");
   ntuple->Branch("mu_diffVtx_x", &mu_diffVtx_x, "mu_diffVtx_x/D");
   ntuple->Branch("mu_diffVtx_y", &mu_diffVtx_y, "mu_diffVtx_y/D");
   ntuple->Branch("mu_diffVtx_z", &mu_diffVtx_z, "mu_diffVtx_z/D");
   ntuple->Branch("mu_diffDir", &mu_diffDir, "mu_diffDir/D");
   ntuple->Branch("mu_diffKE", &mu_diffKE, "mu_diffKE/D");
   ntuple->Branch("mu_diffVtxLong", &mu_diffVtxLong, "mu_diffVtxLong/D");
   ntuple->Branch("mu_diffVtxTrans", &mu_diffVtxTrans, "mu_diffVtxTrans/D");
   ntuple->Branch("mu_diff_nu_E", &mu_diff_nu_E, "mu_diff_nu_E/D");
   ntuple->Branch("e_reco_dwall", &e_reco_dwall, "e_reco_dwall/D");
   ntuple->Branch("e_reco_towall", &e_reco_towall, "e_reco_towall/D");
   ntuple->Branch("e_reco_nu_E", &e_reco_nu_E, "e_reco_nu_E/D");
   ntuple->Branch("e_recoKE", &e_recoKE, "e_recoKE/D");
   ntuple->Branch("e_diffVtx", &e_diffVtx, "e_diffVtx/D");
   ntuple->Branch("e_diffVtx_x", &e_diffVtx_x, "e_diffVtx_x/D");
   ntuple->Branch("e_diffVtx_y", &e_diffVtx_y, "e_diffVtx_y/D");
   ntuple->Branch("e_diffVtx_z", &e_diffVtx_z, "e_diffVtx_z/D");
   ntuple->Branch("e_diffDir", &e_diffDir, "e_diffDir/D");
   ntuple->Branch("e_diffKE", &e_diffKE, "e_diffKE/D");
   ntuple->Branch("e_diffVtxLong", &e_diffVtxLong, "e_diffVtxLong/D");
   ntuple->Branch("e_diffVtxTrans", &e_diffVtxTrans, "e_diffVtxTrans/D");
   ntuple->Branch("e_diff_nu_E", &e_diff_nu_E, "e_diff_nu_E/D");
   ntuple->Branch("nClusters",&nClusters,"nClusters/I");
   ntuple->Branch("recoNRings",&recoNRings, "recoNRings[nClusters]/I");
   ntuple->Branch("isHighECluster",&isHighE, "isHighECluster[nClusters]/O");
   ntuple->Branch("ring1PEs", &ring1PEs, "ring1PEs/I");
   ntuple->Branch("ring2PEs", &ring2PEs, "ring2PEs/I");
   ntuple->Branch("recoPID", &recoPID0, "recoPID/I");
   ntuple->Branch("recoPIDLikelihood", &recoPIDLikelihood, "recoPIDLikelihood/D");
   ntuple->Branch("isHighE", &isHighE0, "isHighE/O");
   ntuple->Branch("nSubevents",&nSubEvts,"nSubevents/I");
   ntuple->Branch("recoTime",recoTime,"recoTime[nSubevents]/D");
   ntuple->Branch("recoEnergy",recoEnergy,"recoEnergy[nSubevents]/D");
   ntuple->Branch("recoCaptures",&recoCaptures,"recoCaptures/I");
   ntuple->Branch("recoMichels",&recoMichels,"recoMichel/I");

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = GetEntry(jentry);
      nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if ((jentry + 1) % 1000 == 0) std::cout << jentry + 1 << std::endl;

      lepton_KE = 0;
      nChkvPart = 0;
      ne = 0;
      nmu = 0;
      npi0 = 0;
      npipm = 0;
      nk0 = 0;
      nkpm = 0;
      nneutron = 0;
      nproton = 0;
      ngamma = 0;
      nother = 0;
      E_pi_pm = {0,0,0,0,0,0,0,0,0,0};
      for (int itrk = 0; itrk < ntrks; itrk++) {
         // Find highest KE outgoing charged lepton
         if (TMath::Abs(mpid[itrk]) == 11 || (TMath::Abs(mpid[itrk]) == 13 && KE[itrk] > lepton_KE)) lepton_KE = KE[itrk];

         //Count particle types
         if (TMath::Abs(mpid[itrk]) == 11) ne++;
         else if (TMath::Abs(mpid[itrk]) == 13) nmu++;
         else if (mpid[itrk] == 22) ngamma++;
         else if (mpid[itrk] == 111) npi0++;
         else if (TMath::Abs(mpid[itrk]) == 211){
             if(npipm<10)E_pi_pm[npipm] = KE[itrk];
             npipm++;
         }
         else if (mpid[itrk] == 130) nk0++;
         else if (mpid[itrk] == 310) nk0++;
         else if (TMath::Abs(mpid[itrk]) == 321)nkpm++;
         else if (mpid[itrk] == 2212) nproton++;
         else if (mpid[itrk] == 2112) nneutron++;
         else nother++;

         //Check for cherenkov producing particles

         //Ignore particles not producing rings
         if (TMath::Abs(mpid[itrk]) == 12 || TMath::Abs(mpid[itrk]) == 14 || TMath::Abs(mpid[itrk]) == 2112 ||
             TMath::Abs(mpid[itrk]) == 2212)
            continue;

         //Treat pi0, k0, etc as decaying giving 2 rings (since no info on decay products)
         if (TMath::Abs(mpid[itrk]) == 111 ||
             TMath::Abs(mpid[itrk]) == 130 ||
             TMath::Abs(mpid[itrk]) == 310 ||
             TMath::Abs(mpid[itrk]) == 221 ||
             TMath::Abs(mpid[itrk]) == 3122) {
            nChkvPart += 2;
            continue;
         }

         //Find mass of charged particles
         double m = (TMath::Abs(mpid[itrk]) == 321) ? 493.7 : // K+/-
                    (TMath::Abs(mpid[itrk]) == 211) ? 139.6 : // pi+/-
                    (TMath::Abs(mpid[itrk]) == 13) ? 105.7 : // mu+/-
                    0.511;                                   // e+/-
         const double n_ref = 1.34;
         double beta = TMath::Sqrt(1 - 1 / TMath::Power((KE[itrk] / m) + 1, 2));
         if (beta > 1. / n_ref) nChkvPart++;
      }
      if (nClusters == 0 || !isHighE[0]) {
         ring1PEs = 0;
         ring2PEs = 0;
         isHighE0 = false;
         recoPID0 = 0;
         recoPIDLikelihood = 0;
         mu_reco_dwall = 0;
         mu_reco_towall = 0;  //
         mu_reco_nu_E = 0;
         mu_recoKE = 0;
         mu_diffVtx = 0;
         mu_diffDir = 0;
         mu_diffKE = 0;
         mu_diffVtxLong = 0;  //
         mu_diffVtxTrans = 0; //
         mu_diff_nu_E = 0;
         e_reco_dwall = 0;
         e_reco_towall = 0;
         e_reco_nu_E = 0;
         e_recoKE = 0;
         e_diffVtx = 0;
         e_diffDir = 0;
         e_diffKE = 0;
         e_diffVtxLong = 0;
         e_diffVtxTrans = 0;
         e_diff_nu_E = 0;
      }
      else {
         ring1PEs = ringPEs[0];
         ring2PEs = (recoNRings[0] > 1 ? ringPEs[1] : 0);
         isHighE0 = isHighE[0];
         recoPID0 = recoPID[0];

         //muon hyp.
         mu_recoKE = recoEnergyHighEMuon[0];
         mu_diffKE = mu_recoKE-lepton_KE;

         double mn = 939;
         double ml = 105;
         double recoE = recoEnergyHighEMuon[0] + ml;
         double recoP = TMath::Sqrt(recoE * recoE - ml * ml);
         mu_reco_nu_E = (mn * recoE - ml * ml / 2.) / (mn - recoE + recoP * recoDirZHighEMuon[0]);
         mu_diff_nu_E = mu_reco_nu_E - neutrino_E;

         mu_diffVtx_x = recoVtxXHighEMuon[0]-trueVtxX;
         mu_diffVtx_y = recoVtxYHighEMuon[0]-trueVtxY;
         mu_diffVtx_z = recoVtxZHighEMuon[0]-trueVtxZ;
         mu_diffVtx = TMath::Sqrt(TMath::Power(recoVtxXHighEMuon[0]-trueVtxX,2)
                                  +TMath::Power(recoVtxYHighEMuon[0]-trueVtxY,2)
                                  +TMath::Power(recoVtxZHighEMuon[0]-trueVtxZ,2));
         mu_diffDir = TMath::ACos(recoDirXHighEMuon[0]*trueDirX
                                  +recoDirYHighEMuon[0]*trueDirY
                                  +recoDirZHighEMuon[0]*trueDirZ);

         double recoDWallZ = 1100 - TMath::Abs(recoVtxZHighEMuon[0]);
         double recoVtxR2 = recoVtxXHighEMuon[0] * recoVtxXHighEMuon[0] + recoVtxYHighEMuon[0] * recoVtxYHighEMuon[0];
         double recoDWallR = 550 - TMath::Sqrt(recoVtxR2);
         mu_reco_dwall = recoDWallR < recoDWallZ ? recoDWallR : recoDWallZ;
         double a = 1 - recoDirZHighEMuon[0] * recoDirZHighEMuon[0];
         double b = recoVtxXHighEMuon[0] * recoDirXHighEMuon[0] + recoVtxYHighEMuon[0] * recoDirYHighEMuon[0];
         double c = recoVtxR2 - 550 * 550;
         double recoToWallR = (TMath::Sqrt(b * b - a * c) - b) / a;
         double recoToWallZ = 1100 - recoVtxZHighEMuon[0] * TMath::Abs(recoDirZHighEMuon[0]);
         mu_reco_towall = recoToWallR < recoToWallZ ? recoToWallR : recoToWallZ;

         mu_diffVtxLong = (recoVtxXHighEMuon[0]-trueVtxX)*trueDirX
                          +(recoVtxYHighEMuon[0]-trueVtxY)*trueDirY
                          +(recoVtxZHighEMuon[0]-trueVtxZ)*trueDirZ;
         mu_diffVtxTrans = TMath::Sqrt(mu_diffVtx*mu_diffVtx - mu_diffVtxLong*mu_diffVtxLong);

         //electron hyp.
         e_recoKE = recoEnergyHighEElectron[0];
         e_diffKE = e_recoKE-lepton_KE;

         mn = 939;
         ml = 0.511;
         recoE = recoEnergyHighEElectron[0] + ml;
         recoP = TMath::Sqrt(recoE * recoE - ml * ml);
         e_reco_nu_E = (mn * recoE - ml * ml / 2.) / (mn - recoE + recoP * recoDirZHighEElectron[0]);
         e_diff_nu_E = e_reco_nu_E - neutrino_E;

         e_diffVtx_x = recoVtxXHighEElectron[0]-trueVtxX;
         e_diffVtx_y = recoVtxYHighEElectron[0]-trueVtxY;
         e_diffVtx_z = recoVtxZHighEElectron[0]-trueVtxZ;
         e_diffVtx = TMath::Sqrt(TMath::Power(recoVtxXHighEElectron[0]-trueVtxX,2)
                                  +TMath::Power(recoVtxYHighEElectron[0]-trueVtxY,2)
                                  +TMath::Power(recoVtxZHighEElectron[0]-trueVtxZ,2));
         e_diffDir = TMath::ACos(recoDirXHighEElectron[0]*trueDirX
                                  +recoDirYHighEElectron[0]*trueDirY
                                  +recoDirZHighEElectron[0]*trueDirZ);

         recoDWallZ = 1100 - TMath::Abs(recoVtxZHighEElectron[0]);
         recoVtxR2 = recoVtxXHighEElectron[0] * recoVtxXHighEElectron[0] + recoVtxYHighEElectron[0] * recoVtxYHighEElectron[0];
         recoDWallR = 550 - TMath::Sqrt(recoVtxR2);
         e_reco_dwall = recoDWallR < recoDWallZ ? recoDWallR : recoDWallZ;
         a = 1 - recoDirZHighEElectron[0] * recoDirZHighEElectron[0];
         b = recoVtxXHighEElectron[0] * recoDirXHighEElectron[0] + recoVtxYHighEElectron[0] * recoDirYHighEElectron[0];
         c = recoVtxR2 - 550 * 550;
         recoToWallR = (TMath::Sqrt(b * b - a * c) - b) / a;
         recoToWallZ = 1100 - recoVtxZHighEElectron[0] * TMath::Abs(recoDirZHighEElectron[0]);
         e_reco_towall = recoToWallR < recoToWallZ ? recoToWallR : recoToWallZ;

         e_diffVtxLong = (recoVtxXHighEElectron[0]-trueVtxX)*  trueDirX
                          +(recoVtxYHighEElectron[0]-trueVtxY)*trueDirY
                          +(recoVtxZHighEElectron[0]-trueVtxZ)*trueDirZ;
         e_diffVtxTrans = TMath::Sqrt(e_diffVtx*e_diffVtx - e_diffVtxLong*e_diffVtxLong);

         recoPIDLikelihood = recoLnLHighEMuon[0] - recoLnLHighEElectron[0];
         // Loop through the reconstruction and count the number of michels (should really do at reco stage but put here for now)
         recoMichels = 0;
         for(int is=0;is<nSubEvts;++is){
            if(recoTime[is] > 135 && recoTime[is]<8000 && recoEnergy[is]>15.){
              recoMichels++;
            }
         }

      }
      ntuple->Fill();
      // quick break for testing
 //     if(jentry>1e6)break;
   }
   ntuple->Write();
}
