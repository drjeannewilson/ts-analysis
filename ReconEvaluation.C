#include "TCut.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "iostream"
#include "TF1.h"
#include "TFile.h"
#include "TH3D.h"
#include "TTree.h"
#include "TROOT.h"

using namespace std;

void ReconEvaluation(bool fhc = true) {
  TFile * f= new TFile("selection_test_20160128_53815.root", "READ");;
  TTree *selection = (TTree *) f->Get("selection");

  TCut R1 = "isHighE&&(ring2PEs/ring1PEs)<0.09";
  TCut mulike = "recoPIDLikelihood>0";
  TCut elike = "recoPIDLikelihood<-400";
  TCut nodecay = "recoMichels==0";
  TCut onedecay = "recoMichels<=1";
  TCut muEnergyCut = "mu_reco_nu_E>0&&mu_reco_nu_E<1250&&mu_recoKE>200&&mu_recoKE<2000";
  TCut eEnergyCut = "e_reco_nu_E>0&&e_reco_nu_E<2500&&e_recoKE>100&&e_recoKE<2500";
  TCut mudwallcut = "mu_reco_dwall>200";
  TCut mutowallcut = "mu_reco_towall>200";
  TCut edwallcut = "e_reco_dwall>200";
  TCut etowallcut = "e_reco_towall>200";
  TCut nonfoll = "recoCaptures==0";
  TCut nfolls = "recoCaptures>0";
  
  // Is this the right way to evaluate the truth?
  TCut mufvtrue = "true_dwall>200&&true_towall>200";// &&neutrino_E>0&&neutrino_E<1250 &&trueKE>200&&trueKE<2000";
  TCut efvtrue = "true_dwall>200&&true_towall>200";// &&neutrino_E>0&&neutrino_E<2500 &&trueKE>100&&trueKE<2500";

  TCut mu1R = R1 && mulike && onedecay && muEnergyCut && mudwallcut && mutowallcut;
  TCut e1R = R1 && elike && nodecay && eEnergyCut && edwallcut && etowallcut;
  
  TCut mu1Rnofoll = mu1R && nonfoll;
  TCut mu1Rfolls  = mu1R && nfolls;
  TCut e1Rnofoll  = e1R && nonfoll;
  TCut e1Rfolls   = e1R && nfolls;
  
  TCut pi0 = "interaction_mode==31||interaction_mode==32||interaction_mode==36";
  TCut weightNCpi ="(abs(interaction_mode)==31||abs(interaction_mode)==32||abs(interaction_mode)==36)*(rndm>0.88889)||!(abs(interaction_mode)==31||abs(interaction_mode==32)||abs(interaction_mode)==36)";

  TCut truemu = "neutrino_id==14";
  TCut truee = "neutrino_id==12";
  TCut truecc = "interaction_mode==1";
  string tag = "FHC";

  TCut weightFHC = "(30.7419*(neutrino_id==14)+0.584919*(neutrino_id==12)+0.875846*(neutrino_id==-14)+0.0590713*(neutrino_id==-12))";
  TCut weightRHC = "(3.90061*(neutrino_id==14)+0.172118*(neutrino_id==12)+7.6609*(neutrino_id==-14)+0.147396*(neutrino_id==-12))";
  TCut weight = weightFHC;
  
  if(!fhc){
    truemu = "neutrino_id==-14";
    truee = "neutrino_id==-12";
    truecc = "interaction_mode==-1";
    tag = "FHC";
    weight = weightRHC;
  }
  
  // open an output root file
  TFile *o;
  o = new TFile(Form("recon_properties_%s.root",tag.c_str()), "RECREATE");


  // Firstly deal with the purities - here we evaluate purities within the FV only
  
  // count all the true e/mu like events with true vertices within the FV
  TH1F *muraw_neut = new TH1F("muraw_neut","mu all, neut codes",101,-50.5,50.5);
  TH1F *eraw_neut = new TH1F("eraw_neut","e all, neut codes",101,-50.5,50.5);
  
  // count all the events selected with true vertices within the FV (mu and e, with and without considering n followers)
  TH1F *mu1R_neut = new TH1F("mu1R_neut","mu1R selection, neut codes",101,-50.5,50.5);
  TH1F *mu1Rnofoll_neut = new TH1F("mu1Rnofoll_neut","mu1R selection, no n followers, neut codes",101,-50.5,50.5);
  TH1F *mu1Rfolls_neut = new TH1F("mu1Rfolls_neut","mu1R selection, with >=1 n, neut codes",101,-50.5,50.5);
  TH1F *e1R_neut = new TH1F("e1R_neut","e1R selection, neut codes",101,-50.5,50.5);
  TH1F *e1Rnofoll_neut = new TH1F("e1Rnofoll_neut","e1R selection, no n followers, neut codes",101,-50.5,50.5);
  TH1F *e1Rfolls_neut = new TH1F("e1Rfolls_neut","e1R selection, with >=1 n, neut codes",101,-50.5,50.5);

  mufvtrue = "1";
  efvtrue = "1";
  weightNCpi = "1";
  selection->Draw("interaction_mode>>muraw_neut",weight*(mufvtrue&&truemu));
  selection->Draw("interaction_mode>>eraw_neut",weight*(efvtrue&&truee));

  selection->Draw("interaction_mode>>mu1R_neut",weight*weightNCpi*(mu1R&&mufvtrue));
  selection->Draw("interaction_mode>>mu1Rnofoll_neut",weight*weightNCpi*(mu1Rnofoll&&mufvtrue));
  selection->Draw("interaction_mode>>mu1Rfolls_neut",weight*weightNCpi*(mu1Rfolls&&mufvtrue));

  selection->Draw("interaction_mode>>e1R_neut",weight*weightNCpi*(e1R&&efvtrue));
  selection->Draw("interaction_mode>>e1Rnofoll_neut",weight*weightNCpi*(e1Rnofoll&&efvtrue));
  selection->Draw("interaction_mode>>e1Rfolls_neut",weight*weightNCpi*(e1Rfolls&&efvtrue));

  // Write files out
  o->cd();
  muraw_neut->Write();
  eraw_neut->Write();
  
  mu1R_neut->Write();
  mu1Rnofoll_neut->Write();
  mu1Rfolls_neut->Write();
  e1R_neut->Write();
  e1Rnofoll_neut->Write();
  e1Rfolls_neut->Write();
  
  // Output the purities
  if(fhc){
    cout << tag << ": Purities: mu 1R all " << mu1R_neut->GetBinContent(52)/mu1R_neut->Integral() << " with no n followers = " << mu1Rnofoll_neut->GetBinContent(52)/mu1Rnofoll_neut->Integral() << endl;
    cout << tag << ": Purities: e 1R all " << e1R_neut->GetBinContent(52)/e1R_neut->Integral() << " with no n followers = " << e1Rnofoll_neut->GetBinContent(52)/e1Rnofoll_neut->Integral() << endl;
  }else{
     cout << tag << ": Purities: mu 1R all " << mu1R_neut->GetBinContent(50)/mu1R_neut->Integral() << " with n followers = " << mu1Rfolls_neut->GetBinContent(50)/mu1Rfolls_neut->Integral() << endl;
     cout << tag << ": Purities: e 1R all " << e1R_neut->GetBinContent(50)/e1R_neut->Integral() << " with n followers = " << e1Rfolls_neut->GetBinContent(50)/e1Rfolls_neut->Integral() << endl;
  }
  // Output the efficiencies
  if(fhc){
    cout << tag << ": Efficiency: mu 1R all " << mu1R_neut->GetBinContent(52)/muraw_neut->GetBinContent(52) << " with no n followers = " << mu1Rnofoll_neut->GetBinContent(52)/muraw_neut->GetBinContent(52)<< endl;
    cout << tag << ": Efficiency: e 1R all " << e1R_neut->GetBinContent(52)/eraw_neut->GetBinContent(52) << " with no n followers = " << e1Rnofoll_neut->GetBinContent(52)/eraw_neut->GetBinContent(52) << endl;
  }else{
     cout << tag << ": Efficiency: mu 1R all " << mu1R_neut->GetBinContent(50)/muraw_neut->GetBinContent(50) << " with n followers = " << mu1Rfolls_neut->GetBinContent(50)/muraw_neut->GetBinContent(50) << endl;
     cout << tag << ": Efficiency: e 1R all " << e1R_neut->GetBinContent(50)/eraw_neut->GetBinContent(50) << " with n followers = " << e1Rfolls_neut->GetBinContent(50)/eraw_neut->GetBinContent(50) << endl;
  }
/* --> results:
Within true FV
FHC: Purities: mu 1R all 0.80447 with no n followers = 0.881322
FHC: Purities: e 1R all 0.220577 with no n followers = 0.355524
FHC: Efficiency: mu 1R all 0.54313 with no n followers = 0.381125
FHC: Efficiency: e 1R all 0.368206 with no n followers = 0.280311

No FV cut
FHC: Purities: mu 1R all 0.803583 with no n followers = 0.880961
FHC: Purities: e 1R all 0.213834 with no n followers = 0.353961
FHC: Efficiency: mu 1R all 0.189269 with no n followers = 0.132795
FHC: Efficiency: e 1R all 0.133256 with no n followers = 0.101709

No NCpi weighting
FHC: Purities: mu 1R all 0.802075 with no n followers = 0.880414
FHC: Purities: e 1R all 0.11527 with no n followers = 0.159646
FHC: Efficiency: mu 1R all 0.189269 with no n followers = 0.132795
FHC: Efficiency: e 1R all 0.133256 with no n followers = 0.101709
*/
  TCanvas *CT = new TCanvas("CT");
  CT->cd();
  eraw_neut->Draw();
  e1R_neut->SetLineColor(2);
  e1R_neut->Draw("same");
  CT->Print("CT.C");
  
  // Now extract the energy resolutions for the selection
  // recon energies for true CC events, evaluate 68% quantile in MeV
  double q=0.68, r;
  
  TH1F *hEl_1Rmu_nofoll = new TH1F("hEl_1Rmu_nofoll","1Rmu selection, dE lepton, no n followers",100,-1000,1000);
  TH1F *hEl_1Re_nofoll = new TH1F("hEl_1Re_nofoll","1Re selection, dE lepton, no n followers",100,-1000,1000);
  TH1F *hEnu_1Rmu_nofoll = new TH1F("hEnu_1Rmu_nofoll","1Rmu selection, dE nu,  no n followers",100,-1000,1000);
  TH1F *hEnu_1Re_nofoll = new TH1F("hEnu_1Re_nofoll","1Re selection, dE nu, no n followers",100,-1000,1000);
  selection->Draw("mu_recoKE-trueKE>>hEl_1Rmu_nofoll",weight*(truemu&&truecc&&mu1Rnofoll));
  selection->Draw("e_recoKE-trueKE>>hEl_1Re_nofoll",weight*(e1Rnofoll&&truee&&truecc));
  selection->Draw("mu_reco_nu_E-neutrino_E>>hEnu_1Rmu_nofoll",weight*(mu1Rnofoll&&truemu&&truecc));
  selection->Draw("e_reco_nu_E-neutrino_E>>hEnu_1Re_nofoll",weight*(e1Rnofoll&&truee&&truecc));
  TCanvas *CrE = new TCanvas("CrE");
  CrE->Divide(2,2);
  CrE->cd(1);
  hEl_1Rmu_nofoll->Draw();
  r = hEl_1Rmu_nofoll->GetRMS();
  cout << "Muon lepton energy res: " << r << " MeV" << endl;
  CrE->cd(2);
  hEl_1Re_nofoll->Draw();
  r = hEl_1Re_nofoll->GetRMS();
  cout << "Electron lepton energy res: " << r << " MeV" << endl;
  CrE->cd(3);
  hEnu_1Rmu_nofoll->Draw();
  r = hEnu_1Rmu_nofoll->GetRMS();
  cout << "Muon neutrino energy res: " << r << " MeV" << endl;
  CrE->cd(4);
  hEnu_1Re_nofoll->Draw();
  r = hEnu_1Re_nofoll->GetRMS();
  cout << "Electron neutrino energy res: " << r << " MeV" << endl;
  CrE->Print("CrE.C");

  // now E in percent
  TH1F *hEl_1Rmu_pc_nofll = new TH1F("hEl_1Rmu_pc_nofll","1Rmu selection, dE lepton, no n followers",100,-100,100);
  TH1F *hEl_1Re_pc_nofll = new TH1F("hEl_1Re_pc_nofll","1Re selection, dE lepton, no n followers",100,-100,100);
  TH1F *hEnu_1Rmu_pc_nofll = new TH1F("hEnu_1Rmu_pc_nofll","1Rmu selection, dE nu,  no n followers",100,-100,100);
  TH1F *hEnu_1Re_pc_nofll = new TH1F("hEnu_1Re_pc_nofll","1Re selection, dE nu, no n followers",100,-100,100);
  selection->Draw("100*(mu_recoKE-trueKE)/trueKE>>hEl_1Rmu_pc_nofll",weight*(mu1Rnofoll&&truemu&&truecc));
  selection->Draw("100*(e_recoKE-trueKE)/trueKE>>hEl_1Re_pc_nofll",weight*(e1Rnofoll&&truee&&truecc));
  selection->Draw("100*(mu_reco_nu_E-neutrino_E)/neutrino_E>>hEnu_1Rmu_pc_nofll",weight*(mu1Rnofoll&&truemu&&truecc));
  selection->Draw("100*(e_reco_nu_E-neutrino_E)/neutrino_E>>hEnu_1Re_pc_nofll",weight*(e1Rnofoll&&truee&&truecc));
  TCanvas *CrEpc = new TCanvas("CrEpc");
  CrEpc->Divide(2,2);
  CrEpc->cd(1);
  hEl_1Rmu_pc_nofll->Draw();
  r = hEl_1Rmu_pc_nofll->GetRMS();
  cout << "Muon lepton energy res: " << r << " %" << endl;
  CrEpc->cd(2);
  hEl_1Re_pc_nofll->Draw();
  r = hEl_1Re_pc_nofll->GetRMS();
  cout << "Electron lepton energy res: " << r << " %" << endl;
  CrEpc->cd(3);
  hEnu_1Rmu_pc_nofll->Draw();
  r = hEnu_1Rmu_pc_nofll->GetRMS();
  cout << "Muon neutrino energy res: " << r << " %" << endl;
  CrEpc->cd(4);
  hEnu_1Re_pc_nofll->Draw();
  r = hEnu_1Re_pc_nofll->GetRMS();
  cout << "Electron neutrino energy res: " << r << " %" << endl;
  CrEpc->Print("CrEpcpc.C");
  
  // Position reconstruction
  TH1F *hvtx_1Rmu_nofll = new TH1F("hvtx_1Rmu_nofll","1Rmu selection, dVtx, no n followers",100,-200,200);
  TH1F *hvtx_1Re_nofll = new TH1F("hvtx_1Re_nofll","1Re selection, dVtx, no n followers",100,-200,200);
  selection->Draw("mu_diffVtx>>hvtx_1Rmu_nofll",weight*(mu1Rnofoll&&truemu&&truecc));
  selection->Draw("e_diffVtx>>hvtx_1Re_nofll",weight*(e1Rnofoll&&truee&&truecc));
  TCanvas *CVtx = new TCanvas("CVtx");
  CVtx->Divide(1,2);
  CVtx->cd(1);
  hvtx_1Rmu_nofll->Draw();
  hvtx_1Rmu_nofll->GetQuantiles(1,&r,&q);
  cout << "Muon position resolution: " << r << " cm " << endl;
  CVtx->cd(2);
  hvtx_1Re_nofll->Draw();
  hvtx_1Re_nofll->GetQuantiles(1,&r,&q);
  cout << "Electron position resolution: " << r << " cm " << endl;
  CVtx->Print("CVtx.C");

  // Position reconstruction
  TCanvas *CVtx_xyz = new TCanvas("CVtx_xyz");
  gStyle->SetOptFit(1);
  CVtx_xyz->Divide(2,3);
  CVtx_xyz->cd(1);
  TH1F *hvtx_x_1Rmu_nofll = new TH1F("hvtx_x_1Rmu_nofll","1Rmu selection, dVtx_x, no n followers",100,-200,200);
  TH1F *hvtx_x_1Re_nofll = new TH1F("hvtx_x_1Re_nofll","1Re selection, dVtx_x, no n followers",100,-200,200);
  selection->Draw("mu_diffVtx_x>>hvtx_x_1Rmu_nofll",weight*(mu1Rnofoll&&truemu&&truecc));
  selection->Draw("e_diffVtx_x>>hvtx_x_1Re_nofll",weight*(e1Rnofoll&&truee&&truecc));
  hvtx_x_1Rmu_nofll->Draw();hvtx_x_1Rmu_nofll->Fit("gaus","Q");
  //hvtx_x_1Rmu_nofll->GetQuantiles(1,&r,&q);
  float vtx_x_1Rmu_nofll_sigma =hvtx_x_1Rmu_nofll->GetFunction("gaus")->GetParameter(2);
  cout << "Muon _x position resolution: sigma =  " << vtx_x_1Rmu_nofll_sigma << " cm " << endl;
  CVtx_xyz->cd(2);
  hvtx_x_1Re_nofll->Draw();hvtx_x_1Re_nofll->Fit("gaus","Q");
  //hvtx_x_1Re_nofll->GetQuantiles(1,&r,&q);
  float vtx_x_1Re_nofll_sigma = hvtx_x_1Re_nofll->GetFunction("gaus")->GetParameter(2);
  cout << "Electron _x position resolution: sigma =  " << vtx_x_1Re_nofll_sigma << " cm " << endl;
  TH1F *hvtx_y_1Rmu_nofll = new TH1F("hvtx_y_1Rmu_nofll","1Rmu selection, dVtx_y, no n followers",100,-200,200);
  TH1F *hvtx_y_1Re_nofll = new TH1F("hvtx_y_1Re_nofll","1Re selection, dVtx_y, no n followers",100,-200,200);
  selection->Draw("mu_diffVtx_y>>hvtx_y_1Rmu_nofll",weight*(mu1Rnofoll&&truemu&&truecc));
  selection->Draw("e_diffVtx_y>>hvtx_y_1Re_nofll",weight*(e1Rnofoll&&truee&&truecc));
  CVtx_xyz->cd(3);
  hvtx_y_1Rmu_nofll->Draw();hvtx_y_1Rmu_nofll->Fit("gaus","Q");
  //hvtx_y_1Rmu_nofll->GetQuantiles(1,&r,&q);
  float vtx_y_1Rmu_nofll_sigma =hvtx_y_1Rmu_nofll->GetFunction("gaus")->GetParameter(2);
  cout << "Muon _y position resolution: sigma =  " << vtx_y_1Rmu_nofll_sigma << " cm  " << endl;
  CVtx_xyz->cd(4);
  hvtx_y_1Re_nofll->Draw();hvtx_y_1Re_nofll->Fit("gaus","Q");
  //hvtx_y_1Re_nofll->GetQuantiles(1,&r,&q);
  float vtx_y_1Re_nofll_sigma =hvtx_y_1Re_nofll->GetFunction("gaus")->GetParameter(2);
  cout << "Electron _y position resolution: sigma =  " << vtx_y_1Re_nofll_sigma << " cm  " << endl;
  TH1F *hvtx_z_1Rmu_nofll = new TH1F("hvtx_z_1Rmu_nofll","1Rmu selection, dVtx_z, no n followers",100,-200,200);
  TH1F *hvtx_z_1Re_nofll = new TH1F("hvtx_z_1Re_nofll","1Re selection, dVtx_z, no n followers",100,-200,200);
  selection->Draw("mu_diffVtx_z>>hvtx_z_1Rmu_nofll",weight*(mu1Rnofoll&&truemu&&truecc));
  selection->Draw("e_diffVtx_z>>hvtx_z_1Re_nofll",weight*(e1Rnofoll&&truee&&truecc));
  CVtx_xyz->cd(5);
  hvtx_z_1Rmu_nofll->Draw();hvtx_z_1Rmu_nofll->Fit("gaus","Q");
  //hvtx_z_1Rmu_nofll->GetQuantiles(1,&r,&q);
  float vtx_z_1Rmu_nofll_sigma =hvtx_z_1Rmu_nofll->GetFunction("gaus")->GetParameter(2);
  cout << "Muon _z position resolution: sigma =  " << vtx_z_1Rmu_nofll_sigma << " cm   " << endl;
  CVtx_xyz->cd(6);
  hvtx_z_1Re_nofll->Draw();hvtx_z_1Re_nofll->Fit("gaus","Q");
  //hvtx_z_1Re_nofll->GetQuantiles(1,&r,&q);
  float vtx_z_1Re_nofll_sigma =hvtx_z_1Re_nofll->GetFunction("gaus")->GetParameter(2);
  cout << "Electron _z position resolution: sigma =  " << vtx_z_1Re_nofll_sigma << " cm   " << endl;
  
  // combined resolution
  float vtx_sigma_1Rmu_nofll = pow((pow(vtx_x_1Rmu_nofll_sigma,2)+pow(vtx_y_1Rmu_nofll_sigma,2)+pow(vtx_z_1Rmu_nofll_sigma,2)),0.5);
  float vtx_sigma_1Re_nofll = pow((pow(vtx_x_1Re_nofll_sigma,2)+pow(vtx_y_1Re_nofll_sigma,2)+pow(vtx_z_1Re_nofll_sigma,2)),0.5);
  cout << "Muon vertex resolution = " <<vtx_sigma_1Rmu_nofll << " cm " << endl;
  cout << "Electron vertex resolution = " <<vtx_sigma_1Re_nofll << " cm " << endl;
  CVtx_xyz->Print("CVtx_xyz.C");

  // Direction reconstruction
  TH1F *hdir_1Rmu_nofll = new TH1F("hdir_1Rmu_nofll","1Rmu selection, dDir, no n followers",100,-100,100);
  TH1F *hdir_1Re_nofll = new TH1F("hdir_1Re_nofll","1Re selection, dDir, no n followers",100,-100,100);
  selection->Draw("mu_diffDir*180/TMath::Pi()>>hdir_1Rmu_nofll",weight*(mu1Rnofoll&&truemu&&truecc));
  selection->Draw("e_diffDir*180/TMath::Pi()>>hdir_1Re_nofll",weight*(e1Rnofoll&&truee&&truecc));
  TCanvas *CDir = new TCanvas("CDir");
  CDir->Divide(1,2);
  CDir->cd(1);
  hdir_1Rmu_nofll->Draw();
  hdir_1Rmu_nofll->GetQuantiles(1,&r,&q);
  cout << "Muon position resolution: " << r << " degrees " << endl;
  CDir->cd(2);
  hdir_1Re_nofll->Draw();
  hdir_1Re_nofll->GetQuantiles(1,&r,&q);
  cout << "Electron direction resolution: " << r << " degrees " << endl;
  CDir->Print("CDir.C");


  o->Close();
}
