// A script to optimise and examine the cut selections based on the output of the script selections.C/.h

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

void OptimiseSelection(bool fhc = true) {
  // This is the name of the file created by the selections.C/.h code
  TFile * f;
  if(fhc){
   f = new TFile("selections_androot_tagged_FHC.root", "READ");
  }esle{
   f = new TFile("selections_androot_tagged_RHC.root", "READ");
  }
  TTree *selection = (TTree *) f->Get("AllEvents");

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
  
  TCut pi0 = "abs(mu_recoKE)==31||abs(interaction_mode)==32||abs(interaction_mode)==36";
  TCut weightNCpi = "(((abs(interaction_mode)==31 || abs(interaction_mode)==32 || abs(interaction_mode)==36) || (n_pi0>0)) &&(rndm>0.88889)||!((abs(interaction_mode)==31 || abs(interaction_mode)==32 || abs(interaction_mode)==36)||n_pi0>0))";
  TCut truemu = "neutrino_id==14";
  TCut truee = "neutrino_id==12";
  TCut trueCCQE = "interaction_mode==1";
  TCut trueCCinc = "interaction_mode>=1&&interaction_mode<=26";
  string tag = "FHC";

  TCut weightFHC = "(30.7419*(neutrino_id==14)+0.584919*(neutrino_id==12)+0.875846*(neutrino_id==-14)+0.0590713*(neutrino_id==-12))";
  TCut weightRHC = "(3.90061*(neutrino_id==14)+0.172118*(neutrino_id==12)+7.6609*(neutrino_id==-14)+0.147396*(neutrino_id==-12))";
  TCut weight = weightFHC;
// Use this option for full weighting, but turn it off for cross check against Valor numbers
//  TCut weight = "1";
  
  if(!fhc){
    truemu = "neutrino_id==-14";
    truee = "neutrino_id==-12";
    trueCCQE = "interaction_mode==-1";
    trueCCinc = "interaction_mode<=-1&&interaction_mode>=-26";
    tag = "RHC";
    weight = weightRHC;
  }
 
  
//  Use these options in comparison to the MakeSelections code so that we can see the output for a particular neutrino type only
//  weightNCpi = "1";
//  weight = "1*neutrino_id==14+0*(neutrino_id!=14)";
//  weight = "1*neutrino_id==12+0*(neutrino_id!=12)";

  // Plot the neut IDs for the selections to check what is getting through
  TCanvas *Cneut = new TCanvas("Cneut");
  Cneut->Divide(1,2);
  Cneut->cd(1);
  gPad->SetLogy();
  selection->Draw("interaction_mode>>neut_esel",weight*weightNCpi*(R1&&elike&&nodecay&&eEnergyCut&&edwallcut&&etowallcut));
  Cneut->cd(2);
  gPad->SetLogy();
  selection->Draw("interaction_mode>>neut_musel",weight*weightNCpi*(R1&&mulike&&onedecay&&muEnergyCut&&mudwallcut&&mutowallcut));
  Cneut->Print("Cneut.C");
  TCanvas *Cneut2 = new TCanvas("Cneut2");
  Cneut2->Divide(1,2);
  Cneut2->cd(1);
  gPad->SetLogy();
  selection->Draw("interaction_mode>>neut_esel2",weight*weightNCpi*(R1&&elike&&nodecay&&eEnergyCut&&edwallcut&&etowallcut&&nonfoll));
  Cneut2->cd(2);
  gPad->SetLogy();
  selection->Draw("interaction_mode>>neut_musel2",weight*weightNCpi*(R1&&mulike&&onedecay&&muEnergyCut&&mudwallcut&&mutowallcut&&nonfoll));
  Cneut2->Print("Cneut2.C");
  
  float norm = neut_esel2->Integral(neut_esel2->FindBin(-50),neut_esel2->FindBin(50));
  for(int ii=-50;ii<50;++ii){
    int bin = neut_esel2->FindBin(ii);
    float num = neut_esel2->GetBinContent(bin);
    float frac = num/norm;
    cout << "mode " << ii << " contents = \t" << num << " events ( " << 100*frac << " % )" << endl;
  }
  
  // Make selections here - apply cuts one by one so we can see reduction effect
  // e1R
  // draw energy distributions - make a histogram for each to fill
  // show cut reduction
  const int ncuts = 9;
  float NInCuts_e1R_all[ncuts],NInCuts_e1R_all_int[ncuts];
  float NInCuts_e1R_CCQE[ncuts],NInCuts_e1R_CCQE_int[ncuts];
  float NInCuts_e1R_CCInc[ncuts],NInCuts_e1R_CCInc_int[ncuts];
  // note evaluating the number of events this way doesn't include the weights:
  NInCuts_e1R_all[0] = selection->Draw("e_recoKE>>InCuts_e_all0",weight);
  NInCuts_e1R_all[1] = selection->Draw("e_recoKE>>InCuts_e_all1",weight*(edwallcut));
  NInCuts_e1R_all[2] = selection->Draw("e_recoKE>>InCuts_e_all2",weight*(edwallcut&&etowallcut));
  NInCuts_e1R_all[3] = selection->Draw("e_recoKE>>InCuts_e_all3",weight*(edwallcut&&etowallcut&&R1));
  NInCuts_e1R_all[4] = selection->Draw("e_recoKE>>InCuts_e_all4",weight*(edwallcut&&etowallcut&&R1&&elike));
  NInCuts_e1R_all[5] = selection->Draw("e_recoKE>>InCuts_e_all5",weight*(edwallcut&&etowallcut&&R1&&elike&&nodecay));
  NInCuts_e1R_all[6] = selection->Draw("e_recoKE>>InCuts_e_all6",weight*(edwallcut&&etowallcut&&R1&&elike&&nodecay&&eEnergyCut));
  NInCuts_e1R_all[7] = selection->Draw("e_recoKE>>InCuts_e_all7",weight*weightNCpi*(R1&&elike&&nodecay&&eEnergyCut&&edwallcut&&etowallcut));
  NInCuts_e1R_all[8] = selection->Draw("e_recoKE>>InCuts_e_all8",weight*weightNCpi*(R1&&elike&&nodecay&&eEnergyCut&&edwallcut&&etowallcut&&nonfoll));
  // so use this method and integrate instead for weighted numbers
  NInCuts_e1R_all_int[0] = InCuts_e_all0->Integral();
  NInCuts_e1R_all_int[1] = InCuts_e_all1->Integral();
  NInCuts_e1R_all_int[2] = InCuts_e_all2->Integral();
  NInCuts_e1R_all_int[3] = InCuts_e_all3->Integral();
  NInCuts_e1R_all_int[4] = InCuts_e_all4->Integral();
  NInCuts_e1R_all_int[5] = InCuts_e_all5->Integral();
  NInCuts_e1R_all_int[6] = InCuts_e_all6->Integral();
  NInCuts_e1R_all_int[7] = InCuts_e_all7->Integral();
  NInCuts_e1R_all_int[8] = InCuts_e_all8->Integral();
  
  
  // Now do same for true CCQE
  NInCuts_e1R_CCQE[0] = selection->Draw("e_recoKE>>InCuts_e_CCQE0",weight*(truee&&trueCCQE));
  NInCuts_e1R_CCQE[1] = selection->Draw("e_recoKE>>InCuts_e_CCQE1",weight*(edwallcut&&truee&&trueCCQE));
  NInCuts_e1R_CCQE[2] = selection->Draw("e_recoKE>>InCuts_e_CCQE2",weight*(edwallcut&&etowallcut&&truee&&trueCCQE));
  NInCuts_e1R_CCQE[3] = selection->Draw("e_recoKE>>InCuts_e_CCQE3",weight*(edwallcut&&etowallcut&&R1&&(truee&&trueCCQE)));
  NInCuts_e1R_CCQE[4] = selection->Draw("e_recoKE>>InCuts_e_CCQE4",weight*(edwallcut&&etowallcut&&R1&&elike&&(truee&&trueCCQE)));
  NInCuts_e1R_CCQE[5] = selection->Draw("e_recoKE>>InCuts_e_CCQE5",weight*(edwallcut&&etowallcut&&R1&&elike&&nodecay&&(truee&&trueCCQE)));
  NInCuts_e1R_CCQE[6] = selection->Draw("e_recoKE>>InCuts_e_CCQE6",weight*(edwallcut&&etowallcut&&R1&&elike&&nodecay&&eEnergyCut&&(truee&&trueCCQE)));
  NInCuts_e1R_CCQE[7] = selection->Draw("e_recoKE>>InCuts_e_CCQE7",weight*weightNCpi*(R1&&elike&&nodecay&&eEnergyCut&&edwallcut&&etowallcut&&(truee&&trueCCQE)));
  NInCuts_e1R_CCQE[8] = selection->Draw("e_recoKE>>InCuts_e_CCQE8",weight*weightNCpi*(R1&&elike&&nodecay&&eEnergyCut&&edwallcut&&etowallcut&&nonfoll&&(truee&&trueCCQE)));
  NInCuts_e1R_CCQE_int[0] = InCuts_e_CCQE0->Integral();
  NInCuts_e1R_CCQE_int[1] = InCuts_e_CCQE1->Integral();
  NInCuts_e1R_CCQE_int[2] = InCuts_e_CCQE2->Integral();
  NInCuts_e1R_CCQE_int[3] = InCuts_e_CCQE3->Integral();
  NInCuts_e1R_CCQE_int[4] = InCuts_e_CCQE4->Integral();
  NInCuts_e1R_CCQE_int[5] = InCuts_e_CCQE5->Integral();
  NInCuts_e1R_CCQE_int[6] = InCuts_e_CCQE6->Integral();
  NInCuts_e1R_CCQE_int[7] = InCuts_e_CCQE7->Integral();
  NInCuts_e1R_CCQE_int[8] = InCuts_e_CCQE8->Integral();
  
  // And for true CCInc
  NInCuts_e1R_CCInc[0] = selection->Draw("e_recoKE>>InCuts_e_CCInc0",weight*(truee&&trueCCinc));
  NInCuts_e1R_CCInc[1] = selection->Draw("e_recoKE>>InCuts_e_CCInc1",weight*(edwallcut&&truee&&trueCCinc));
  NInCuts_e1R_CCInc[2] = selection->Draw("e_recoKE>>InCuts_e_CCInc2",weight*(edwallcut&&etowallcut&&truee&&trueCCinc));
  NInCuts_e1R_CCInc[3] = selection->Draw("e_recoKE>>InCuts_e_CCInc3",weight*(edwallcut&&etowallcut&&R1&&(truee&&trueCCinc)));
  NInCuts_e1R_CCInc[4] = selection->Draw("e_recoKE>>InCuts_e_CCInc4",weight*(edwallcut&&etowallcut&&R1&&elike&&(truee&&trueCCinc)));
  NInCuts_e1R_CCInc[5] = selection->Draw("e_recoKE>>InCuts_e_CCInc5",weight*(edwallcut&&etowallcut&&R1&&elike&&nodecay&&(truee&&trueCCinc)));
  NInCuts_e1R_CCInc[6] = selection->Draw("e_recoKE>>InCuts_e_CCInc6",weight*(edwallcut&&etowallcut&&R1&&elike&&nodecay&&eEnergyCut&&(truee&&trueCCinc)));
  NInCuts_e1R_CCInc[7] = selection->Draw("e_recoKE>>InCuts_e_CCInc7",weight*weightNCpi*(R1&&elike&&nodecay&&eEnergyCut&&edwallcut&&etowallcut&&(truee&&trueCCinc)));
  NInCuts_e1R_CCInc[8] = selection->Draw("e_recoKE>>InCuts_e_CCInc8",weight*weightNCpi*(R1&&elike&&nodecay&&eEnergyCut&&edwallcut&&etowallcut&&nonfoll&&(truee&&trueCCinc)));
  NInCuts_e1R_CCInc_int[0] = InCuts_e_CCInc0->Integral();
  NInCuts_e1R_CCInc_int[1] = InCuts_e_CCInc1->Integral();
  NInCuts_e1R_CCInc_int[2] = InCuts_e_CCInc2->Integral();
  NInCuts_e1R_CCInc_int[3] = InCuts_e_CCInc3->Integral();
  NInCuts_e1R_CCInc_int[4] = InCuts_e_CCInc4->Integral();
  NInCuts_e1R_CCInc_int[5] = InCuts_e_CCInc5->Integral();
  NInCuts_e1R_CCInc_int[6] = InCuts_e_CCInc6->Integral();
  NInCuts_e1R_CCInc_int[7] = InCuts_e_CCInc7->Integral();
  NInCuts_e1R_CCInc_int[8] = InCuts_e_CCInc8->Integral();
  
  // Draw cut reduction for the Nue selection
  TCanvas *C1 = new TCanvas("C1","Nu-e");
  C1->Divide(1,3);
  C1->cd(1);
  gPad->SetLogy(); gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  InCuts_e_all0->GetYaxis()->SetRangeUser(1,5e5); InCuts_e_all0->SetXTitle("Reconstructed E_{e} (MeV)");InCuts_e_all0->SetYTitle("Number Events Selected");
  InCuts_e_all0->Draw();
  InCuts_e_all1->SetLineColor(2); InCuts_e_all1->Draw("same");
  InCuts_e_all2->SetLineColor(3); InCuts_e_all2->Draw("same");
  InCuts_e_all3->SetLineColor(4); InCuts_e_all3->Draw("same");
  InCuts_e_all4->SetLineColor(5); InCuts_e_all4->Draw("same");
  InCuts_e_all5->SetLineColor(6); InCuts_e_all5->Draw("same");
  InCuts_e_all6->SetLineColor(7); InCuts_e_all6->Draw("same");
  InCuts_e_all7->SetLineColor(8); InCuts_e_all7->Draw("same");
  InCuts_e_all8->SetLineColor(9); InCuts_e_all8->Draw("same");
  
  TLegend *leg = new TLegend(0.8,0.6,1.0,1.0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(InCuts_e_all0,"All events","L");
  leg->AddEntry(InCuts_e_all1,"1Ring cut","L");
  leg->AddEntry(InCuts_e_all2,"+ e-like Cut","L");
  leg->AddEntry(InCuts_e_all3,"+ no decay e^{-}s","L");
  leg->AddEntry(InCuts_e_all4,"+ NC#pi^{0} reduction","L");
  leg->AddEntry(InCuts_e_all5,"+ energy cut","L");
  leg->AddEntry(InCuts_e_all6,"+ dwall cut","L");
  leg->AddEntry(InCuts_e_all7,"+ towall cut","L");
  leg->AddEntry(InCuts_e_all8,"+ no n followers","L");
  
  C1->cd(2);
  gPad->SetLogy(); gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  InCuts_e_CCQE0->GetYaxis()->SetRangeUser(1,5e4); InCuts_e_CCQE0->SetXTitle("Reconstructed E_{e} (MeV)");InCuts_e_CCQE0->SetYTitle("Number #nu_{e}CCQE Events Selected");
  InCuts_e_CCQE0->Draw();
  InCuts_e_CCQE1->SetLineColor(2); InCuts_e_CCQE1->Draw("same");
  InCuts_e_CCQE2->SetLineColor(3); InCuts_e_CCQE2->Draw("same");
  InCuts_e_CCQE3->SetLineColor(4); InCuts_e_CCQE3->Draw("same");
  InCuts_e_CCQE4->SetLineColor(5); InCuts_e_CCQE4->Draw("same");
  InCuts_e_CCQE5->SetLineColor(6); InCuts_e_CCQE5->Draw("same");
  InCuts_e_CCQE6->SetLineColor(7); InCuts_e_CCQE6->Draw("same");
  InCuts_e_CCQE7->SetLineColor(8); InCuts_e_CCQE7->Draw("same");
  InCuts_e_CCQE8->SetLineColor(9); InCuts_e_CCQE8->Draw("same");
  
  C1->cd(3);
  gPad->SetLogy(); gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  InCuts_e_CCInc0->GetYaxis()->SetRangeUser(1,5e4); InCuts_e_CCInc0->SetXTitle("Reconstructed E_{e} (MeV)");InCuts_e_CCInc0->SetYTitle("Number #nu_{e}CCInc Events Selected");
  InCuts_e_CCInc0->Draw();
  InCuts_e_CCInc1->SetLineColor(2); InCuts_e_CCInc1->Draw("same");
  InCuts_e_CCInc2->SetLineColor(3); InCuts_e_CCInc2->Draw("same");
  InCuts_e_CCInc3->SetLineColor(4); InCuts_e_CCInc3->Draw("same");
  InCuts_e_CCInc4->SetLineColor(5); InCuts_e_CCInc4->Draw("same");
  InCuts_e_CCInc5->SetLineColor(6); InCuts_e_CCInc5->Draw("same");
  InCuts_e_CCInc6->SetLineColor(7); InCuts_e_CCInc6->Draw("same");
  InCuts_e_CCInc7->SetLineColor(8); InCuts_e_CCInc7->Draw("same");
  InCuts_e_CCInc8->SetLineColor(9); InCuts_e_CCInc8->Draw("same");
  leg->Draw();
   cout << ncuts << endl;
  for(int i=0;i<ncuts;++i){
    cout << "After cut " << i << " \tIn 1Re selection =  " <<  NInCuts_e1R_all_int[i] << " \tNe CCQE = " << NInCuts_e1R_CCQE_int[i] << "  \tNe CCInc = " << NInCuts_e1R_CCInc_int[i]<< /*"\t \t" << NInCuts_e1R_all[i] << " \tNe CCQE = " << NInCuts_e1R_CCQE[i] << "  \tNe  = " << NInCuts_e1R_CCInc[i] << */endl;
  }
  C1->Print("ElectronSelectionReduction.pdf");
  
  //1Rmu
  TCanvas *C2 = new TCanvas("C2","Nu-mu");
  C2->Divide(1,3);
  C2->cd(1);
  // show cut reduction
  float NInCuts_mu1R_all[ncuts], NInCuts_mu1R_all_int[ncuts];
  float NInCuts_mu1R_CCQE[ncuts],NInCuts_mu1R_CCQE_int[ncuts];
  float NInCuts_mu1R_CCInc[ncuts],NInCuts_mu1R_CCInc_int[ncuts];
  NInCuts_mu1R_all[0] = selection->Draw("mu_recoKE>>InCuts_mu_all0",weight);
  NInCuts_mu1R_all[1] = selection->Draw("mu_recoKE>>InCuts_mu_all1",weight*(mudwallcut));
  NInCuts_mu1R_all[2] = selection->Draw("mu_recoKE>>InCuts_mu_all2",weight*(mudwallcut&&mutowallcut));
  NInCuts_mu1R_all[3] = selection->Draw("mu_recoKE>>InCuts_mu_all3",weight*(mudwallcut&&mutowallcut&&R1));
  NInCuts_mu1R_all[4] = selection->Draw("mu_recoKE>>InCuts_mu_all4",weight*(mudwallcut&&mutowallcut&&R1&&mulike));
  NInCuts_mu1R_all[5] = selection->Draw("mu_recoKE>>InCuts_mu_all5",weight*(mudwallcut&&mutowallcut&&R1&&mulike&&onedecay));
  NInCuts_mu1R_all[6] = selection->Draw("mu_recoKE>>InCuts_mu_all6",weight*(mudwallcut&&mutowallcut&&R1&&mulike&&onedecay&&muEnergyCut));
  NInCuts_mu1R_all[7] = selection->Draw("mu_recoKE>>InCuts_mu_all7",weight*weightNCpi*(R1&&mulike&&onedecay&&muEnergyCut&&mudwallcut&&mutowallcut));
  NInCuts_mu1R_all[8] = selection->Draw("mu_recoKE>>InCuts_mu_all8",weight*weightNCpi*(R1&&mulike&&onedecay&&muEnergyCut&&mudwallcut&&nonfoll&&mutowallcut));
  NInCuts_mu1R_all_int[0] = InCuts_mu_all0->Integral();
  NInCuts_mu1R_all_int[1] = InCuts_mu_all1->Integral();
  NInCuts_mu1R_all_int[2] = InCuts_mu_all2->Integral();
  NInCuts_mu1R_all_int[3] = InCuts_mu_all3->Integral();
  NInCuts_mu1R_all_int[4] = InCuts_mu_all4->Integral();
  NInCuts_mu1R_all_int[5] = InCuts_mu_all5->Integral();
  NInCuts_mu1R_all_int[6] = InCuts_mu_all6->Integral();
  NInCuts_mu1R_all_int[7] = InCuts_mu_all7->Integral();
  NInCuts_mu1R_all_int[8] = InCuts_mu_all8->Integral();
  
  NInCuts_mu1R_CCQE[0] = selection->Draw("mu_recoKE>>InCuts_mu_CCQE0",weight*(truemu&&trueCCQE));
  NInCuts_mu1R_CCQE[1] = selection->Draw("mu_recoKE>>InCuts_mu_CCQE1",weight*(mudwallcut&&truemu&&trueCCQE));
  NInCuts_mu1R_CCQE[2] = selection->Draw("mu_recoKE>>InCuts_mu_CCQE2",weight*(mudwallcut&&mutowallcut&&truemu&&trueCCQE));
  NInCuts_mu1R_CCQE[3] = selection->Draw("mu_recoKE>>InCuts_mu_CCQE3",weight*(mudwallcut&&mutowallcut&&R1&&(truemu&&trueCCQE)));
  NInCuts_mu1R_CCQE[4] = selection->Draw("mu_recoKE>>InCuts_mu_CCQE4",weight*(mudwallcut&&mutowallcut&&R1&&mulike&&(truemu&&trueCCQE)));
  NInCuts_mu1R_CCQE[5] = selection->Draw("mu_recoKE>>InCuts_mu_CCQE5",weight*(mudwallcut&&mutowallcut&&R1&&mulike&&onedecay&&(truemu&&trueCCQE)));
  NInCuts_mu1R_CCQE[6] = selection->Draw("mu_recoKE>>InCuts_mu_CCQE6",weight*(mudwallcut&&mutowallcut&&R1&&mulike&&onedecay&&muEnergyCut&&(truemu&&trueCCQE)));
  NInCuts_mu1R_CCQE[7] = selection->Draw("mu_recoKE>>InCuts_mu_CCQE7",weight*weightNCpi*(R1&&mulike&&onedecay&&muEnergyCut&&mudwallcut&&mutowallcut&&(truemu&&trueCCQE)));
  NInCuts_mu1R_CCQE[8] = selection->Draw("mu_recoKE>>InCuts_mu_CCQE8",weight*weightNCpi*(R1&&mulike&&onedecay&&muEnergyCut&&mudwallcut&&mutowallcut&&nonfoll&&(truemu&&trueCCQE)));
  NInCuts_mu1R_CCQE_int[0] = InCuts_mu_CCQE0->Integral();
  NInCuts_mu1R_CCQE_int[1] = InCuts_mu_CCQE1->Integral();
  NInCuts_mu1R_CCQE_int[2] = InCuts_mu_CCQE2->Integral();
  NInCuts_mu1R_CCQE_int[3] = InCuts_mu_CCQE3->Integral();
  NInCuts_mu1R_CCQE_int[4] = InCuts_mu_CCQE4->Integral();
  NInCuts_mu1R_CCQE_int[5] = InCuts_mu_CCQE5->Integral();
  NInCuts_mu1R_CCQE_int[6] = InCuts_mu_CCQE6->Integral();
  NInCuts_mu1R_CCQE_int[7] = InCuts_mu_CCQE7->Integral();
  NInCuts_mu1R_CCQE_int[8] = InCuts_mu_CCQE8->Integral();
  
  NInCuts_mu1R_CCInc[0] = selection->Draw("mu_recoKE>>InCuts_mu_CCInc0",weight*(truemu&&trueCCinc));
  NInCuts_mu1R_CCInc[1] = selection->Draw("mu_recoKE>>InCuts_mu_CCInc1",weight*(mudwallcut&&mutowallcut&&trueCCinc));
  NInCuts_mu1R_CCInc[2] = selection->Draw("mu_recoKE>>InCuts_mu_CCInc2",weight*(mudwallcut&&mutowallcut&&truemu&&trueCCinc));
  NInCuts_mu1R_CCInc[3] = selection->Draw("mu_recoKE>>InCuts_mu_CCInc3",weight*(mudwallcut&&mutowallcut&&R1&&(truemu&&trueCCinc)));
  NInCuts_mu1R_CCInc[4] = selection->Draw("mu_recoKE>>InCuts_mu_CCInc4",weight*(mudwallcut&&mutowallcut&&R1&&mulike&&(truemu&&trueCCinc)));
  NInCuts_mu1R_CCInc[5] = selection->Draw("mu_recoKE>>InCuts_mu_CCInc5",weight*(mudwallcut&&mutowallcut&&R1&&mulike&&onedecay&&(truemu&&trueCCinc)));
  NInCuts_mu1R_CCInc[6] = selection->Draw("mu_recoKE>>InCuts_mu_CCInc6",weight*(mudwallcut&&mutowallcut&&R1&&mulike&&onedecay&&muEnergyCut&&(truemu&&trueCCinc)));
  NInCuts_mu1R_CCInc[7] = selection->Draw("mu_recoKE>>InCuts_mu_CCInc7",weight*weightNCpi*(R1&&mulike&&onedecay&&muEnergyCut&&mudwallcut&&mutowallcut&&(truemu&&trueCCinc)));
  NInCuts_mu1R_CCInc[8] = selection->Draw("mu_recoKE>>InCuts_mu_CCInc8",weight*weightNCpi*(R1&&mulike&&onedecay&&muEnergyCut&&mudwallcut&&mutowallcut&&nonfoll&&(truemu&&trueCCinc)));
  NInCuts_mu1R_CCInc_int[0] = InCuts_mu_CCInc0->Integral();
  NInCuts_mu1R_CCInc_int[1] = InCuts_mu_CCInc1->Integral();
  NInCuts_mu1R_CCInc_int[2] = InCuts_mu_CCInc2->Integral();
  NInCuts_mu1R_CCInc_int[3] = InCuts_mu_CCInc3->Integral();
  NInCuts_mu1R_CCInc_int[4] = InCuts_mu_CCInc4->Integral();
  NInCuts_mu1R_CCInc_int[5] = InCuts_mu_CCInc5->Integral();
  NInCuts_mu1R_CCInc_int[6] = InCuts_mu_CCInc6->Integral();
  NInCuts_mu1R_CCInc_int[7] = InCuts_mu_CCInc7->Integral();
  NInCuts_mu1R_CCInc_int[8] = InCuts_mu_CCInc8->Integral();
  
  // Draw cut reduction for the Nue selection
  gPad->SetLogy(); gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  InCuts_mu_all0->GetYaxis()->SetRangeUser(1,5e5); InCuts_mu_all0->SetXTitle("Reconstructed E_{#mu} (MeV)");InCuts_mu_all0->SetYTitle("Number Events Selected");
  InCuts_mu_all0->Draw();
  InCuts_mu_all1->SetLineColor(2); InCuts_mu_all1->Draw("same");
  InCuts_mu_all2->SetLineColor(3); InCuts_mu_all2->Draw("same");
  InCuts_mu_all3->SetLineColor(4); InCuts_mu_all3->Draw("same");
  InCuts_mu_all4->SetLineColor(5); InCuts_mu_all4->Draw("same");
  InCuts_mu_all5->SetLineColor(6); InCuts_mu_all5->Draw("same");
  InCuts_mu_all6->SetLineColor(7); InCuts_mu_all6->Draw("same");
  InCuts_mu_all7->SetLineColor(8); InCuts_mu_all7->Draw("same");
  InCuts_mu_all8->SetLineColor(9); InCuts_mu_all7->Draw("same");
  
  TLegend *leg = new TLegend(0.8,0.6,1.0,1.0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(InCuts_mu_all0,"All events","L");
  leg->AddEntry(InCuts_mu_all1,"1Ring cut","L");
  leg->AddEntry(InCuts_mu_all2,"+ e-like Cut","L");
  leg->AddEntry(InCuts_mu_all3,"+ no decay e^{-}s","L");
  leg->AddEntry(InCuts_mu_all4,"+ NC#pi^{0} reduction","L");
  leg->AddEntry(InCuts_mu_all5,"+ energy cut","L");
  leg->AddEntry(InCuts_mu_all6,"+ dwall cut","L");
  leg->AddEntry(InCuts_mu_all7,"+ towall cut","L");
  leg->AddEntry(InCuts_mu_all8,"+ no n followers","L");
  
  C2->cd(2);
  gPad->SetLogy(); gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  InCuts_mu_CCQE0->GetYaxis()->SetRangeUser(1,5e5); InCuts_mu_CCQE0->SetXTitle("Reconstructed E_{#mu} (MeV)");InCuts_mu_CCQE0->SetYTitle("Number #nu_{#mu}CCQE Events Selected");
  InCuts_mu_CCQE0->Draw();
  InCuts_mu_CCQE1->SetLineColor(2); InCuts_mu_CCQE1->Draw("same");
  InCuts_mu_CCQE2->SetLineColor(3); InCuts_mu_CCQE2->Draw("same");
  InCuts_mu_CCQE3->SetLineColor(4); InCuts_mu_CCQE3->Draw("same");
  InCuts_mu_CCQE4->SetLineColor(5); InCuts_mu_CCQE4->Draw("same");
  InCuts_mu_CCQE5->SetLineColor(6); InCuts_mu_CCQE5->Draw("same");
  InCuts_mu_CCQE6->SetLineColor(7); InCuts_mu_CCQE6->Draw("same");
  InCuts_mu_CCQE7->SetLineColor(8); InCuts_mu_CCQE7->Draw("same");
  InCuts_mu_CCQE8->SetLineColor(9); InCuts_mu_CCQE8->Draw("same");
  
  C2->cd(3);
  gPad->SetLogy(); gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  InCuts_mu_CCInc0->GetYaxis()->SetRangeUser(1,5e5); InCuts_mu_CCInc0->SetXTitle("Reconstructed E_{#mu} (MeV)");InCuts_mu_CCInc0->SetYTitle("Number #nu_{#mu}CCInc Events Selected");
  InCuts_mu_CCInc0->Draw();
  InCuts_mu_CCInc1->SetLineColor(2); InCuts_mu_CCInc1->Draw("same");
  InCuts_mu_CCInc2->SetLineColor(3); InCuts_mu_CCInc2->Draw("same");
  InCuts_mu_CCInc3->SetLineColor(4); InCuts_mu_CCInc3->Draw("same");
  InCuts_mu_CCInc4->SetLineColor(5); InCuts_mu_CCInc4->Draw("same");
  InCuts_mu_CCInc5->SetLineColor(6); InCuts_mu_CCInc5->Draw("same");
  InCuts_mu_CCInc6->SetLineColor(7); InCuts_mu_CCInc6->Draw("same");
  InCuts_mu_CCInc7->SetLineColor(8); InCuts_mu_CCInc7->Draw("same");
  InCuts_mu_CCInc8->SetLineColor(9); InCuts_mu_CCInc8->Draw("same");
  leg->Draw();
   cout << ncuts << endl;
  for(int i=0;i<ncuts;++i){
    cout << "After cut " << i << " \tIn 1Rmu selection =  " <<  NInCuts_mu1R_all_int[i] << " \tNmu CCQE = " << NInCuts_mu1R_CCQE_int[i] << "  \tNmu CCInc = " << NInCuts_mu1R_CCInc_int[i]<< /*"\t \t" << NInCuts_mu1R_all[i]<< " \tNmu CCQE = " << NInCuts_mu1R_CCQE[i] << "  \tNmu  = " << NInCuts_mu1R_CCInc[i] << */ endl;
  }
  C2->Print("MuonSelectionReduction.pdf");

  TCanvas *C3 = new TCanvas("C3");
  // ....
  // Now get number of events selected
  
  float e1R_norm =InCuts_e_all8->Integral();      // 1Re stats
  // get fraction of events selected that are nueCCQE
  float ne1R_nueCCQE =InCuts_e_CCQE8->Integral(); // 1Re CCQE selected
  float frac_e1R_nueCCQE = ne1R_nueCCQE/e1R_norm; // 1Re CCQE purity
  // get fraction of events selected that are nueCC inclusive
  float ne1R_nueCCinc =InCuts_e_CCInc8->Integral();   // 1Re e selected
  float frac_e1R_nueCCinc = ne1R_nueCCinc/e1R_norm;   // 1Re e purity
  // get the total number of nueCCQE events
  float total_nueCCQE = InCuts_e_CCQE0->Integral();   // 1Re CCQE total
  float eff_nueCCQE = ne1R_nueCCQE/total_nueCCQE;     // 1RE CCQE efficiency
  float eff_nueCCQE_FV = ne1R_nueCCQE/InCuts_e_CCQE2->Integral();     // 1RE CCQE effiency in FV
  // get the total number of nueCCinclusive events
  float total_nueCCinc =InCuts_e_CCInc0->Integral();  // 1Re e total
  float eff_nueCCinc = ne1R_nueCCinc/total_nueCCinc;  // 1Re e efficiency
  float eff_nueCCinc_FV = ne1R_nueCCinc/InCuts_e_CCInc2->Integral();  // 1Re e efficiency in FV

  cout << " 1Re selection: stats " << e1R_norm << " CC QE stats " << ne1R_nueCCQE << " purity " << frac_e1R_nueCCQE << " e stats" << ne1R_nueCCinc << " purity" << frac_e1R_nueCCinc << " norm CCQE" << total_nueCCQE << " eff " << eff_nueCCQE << " eff FV " << eff_nueCCQE_FV << " norm e " <<   total_nueCCinc << " eff " << eff_nueCCinc << " eff FV " << eff_nueCCinc_FV << endl;

  //mu1R
  float mu1R_norm =InCuts_mu_all8->Integral();      // muon selection stats
  // get fraction of events selected that are nueCCQE
  float nmu1R_numuCCQE =InCuts_mu_CCQE8->Integral();    // muon selection true CCQE
  float frac_mu1R_numuCCQE = nmu1R_numuCCQE/mu1R_norm;  // purity CCQE
  // get fraction of events selected that are nueCC inclusive
  float nmu1R_numuCCinc =InCuts_mu_CCInc8->Integral();  // muon selection true mu
  float frac_mu1R_numuCCinc = nmu1R_numuCCinc/mu1R_norm;  // purity true mu
  // get the total number of nueCCQE events
  float total_numuCCQE = InCuts_mu_CCQE0->Integral();     // total true mu CCQE
  float eff_numuCCQE = nmu1R_numuCCQE/total_numuCCQE;     // efficiency true CQE
  float eff_numuCCQE_FV = nmu1R_numuCCQE/InCuts_mu_CCQE2->Integral(); // efficiency true CCQE in FV
  // get the total number of nueCCinclusive events
  float total_numuCCinc =InCuts_mu_CCInc0->Integral();    // total true mu
  float eff_numuCCinc = nmu1R_numuCCinc/total_numuCCinc;  // efficiency true mu
  float eff_numuCCinc_FV = nmu1R_numuCCinc/InCuts_mu_CCInc2->Integral();  // efficiency true mu in FV

  cout << " 1Rmu selection stats: " << mu1R_norm << " CC QE stats " << nmu1R_numuCCQE << " purity " << frac_mu1R_numuCCQE << " mu stats " << nmu1R_numuCCinc << " purity " << frac_mu1R_numuCCinc << " norm CCQE " << total_numuCCQE << " eff " << eff_numuCCQE << " eff FV " << eff_numuCCQE_FV << " norm mu " <<   total_numuCCinc << " eff " << eff_numuCCinc << " eff FV " << eff_numuCCinc_FV <<  endl;


  // Try a latex table
  cout << " Sample &  Purity & FV efficiency & Total Efficiency \\\\ \\hline " << endl;
  cout << setprecision(3)<< "1R$\\mu$ (CCQE) & " << 100*InCuts_mu_CCQE7->Integral()/InCuts_mu_all7->Integral() << " & " << 100*InCuts_mu_CCQE7->Integral()/InCuts_mu_CCQE2->Integral() << " & " << 100*InCuts_mu_CCQE7->Integral()/InCuts_mu_CCQE0->Integral() << " \\\\ \\hline" << endl;
  cout << setprecision(3)<< "1R$\\mu$ no n followers (CCQE) & " << 100*InCuts_mu_CCQE8->Integral()/InCuts_mu_all8->Integral() << " & " << 100*InCuts_mu_CCQE8->Integral()/InCuts_mu_CCQE2->Integral() << " & " << 100*InCuts_mu_CCQE8->Integral()/InCuts_mu_CCQE0->Integral() << " \\\\ \\hline" << endl;
  cout << setprecision(3)<< "1Re (CCQE) & " << 100*InCuts_e_CCQE6->Integral()/InCuts_e_all6->Integral() << " & " << 100*InCuts_e_CCQE6->Integral()/InCuts_e_CCQE2->Integral() << " & " << 100*InCuts_e_CCQE6->Integral()/InCuts_e_CCQE0->Integral() << " \\\\ \\hline" << endl;
  cout << setprecision(3)<< "1Re NC$\\pi^0$ reduction (CCQE) & " << 100*InCuts_e_CCQE7->Integral()/InCuts_e_all7->Integral() << " & " << 100*InCuts_e_CCQE7->Integral()/InCuts_e_CCQE2->Integral() << " & " << 100*InCuts_e_CCQE7->Integral()/InCuts_e_CCQE0->Integral() << " \\\\ \\hline" << endl;
  cout << setprecision(3)<< "1Re no n followers (CCQE) & " << 100*InCuts_e_CCQE8->Integral()/InCuts_e_all8->Integral() << " & " << 100*InCuts_e_CCQE8->Integral()/InCuts_e_CCQE2->Integral() << " & " << 100*InCuts_e_CCQE8->Integral()/InCuts_e_CCQE0->Integral() << " \\\\ \\hline" << endl;
  cout <<"\n"<< setprecision(3)<< "1R$\\mu$ n followers (CCQE) & " << 100*(InCuts_mu_CCQE7->Integral()-InCuts_mu_CCQE8->Integral())/(InCuts_mu_all7->Integral()-InCuts_mu_all8->Integral()) << " & " << 100*(InCuts_mu_CCQE7->Integral()-InCuts_mu_CCQE8->Integral())/InCuts_mu_CCQE2->Integral() << " & " << 100*(InCuts_mu_CCQE7->Integral()-InCuts_mu_CCQE8->Integral())/InCuts_mu_CCQE0->Integral() << " \\\\ \\hline" << endl;
  cout <<"\n"<< setprecision(3)<< "1Re n followers (CCQE) & " << 100*(InCuts_e_CCQE7->Integral()-InCuts_e_CCQE8->Integral())/(InCuts_e_all7->Integral()-InCuts_e_all8->Integral()) << " & " << 100*(InCuts_e_CCQE7->Integral()-InCuts_e_CCQE8->Integral())/InCuts_e_CCQE2->Integral() << " & " << 100*(InCuts_e_CCQE7->Integral()-InCuts_e_CCQE8->Integral())/InCuts_e_CCQE0->Integral() << " \\\\ \\hline" << endl;

  
  //NCpi0
  selection->Draw("mu_recoKE>>htot_NCpi0",weight*pi0);
  float total_NCpi0 = htot_NCpi0->Integral();
  cout << "total NCpi0 = " << total_NCpi0 << endl;

}
