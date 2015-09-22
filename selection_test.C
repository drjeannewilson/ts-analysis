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

void selection_test() {
  //bool b = gROOT->IsBatch();
  //gROOT->SetBatch(kTRUE);
  TFile * f = new TFile("selection_test_new.root", "READ");
  TTree *selection = (TTree *) f->Get("selection");

  TCut numuCut = "(neutrino_id)==14";
  TCut nueCut = "(neutrino_id)==12";
  TCut R1 = "isHighE&&ring2PEs/ring1PEs<0.09";
  TCut mulike = "recoPIDLikelihood>-200";
  TCut elike = "recoPIDLikelihood<-200";
  TCut mu1R = R1 && mulike;
  TCut e1R = R1 && elike;
  TCut numuCCQE = numuCut && "interaction_mode==1";//"n_mu==1&&n_chkv_part==1";
  TCut cce = nueCut && "interaction_mode==1";//&&"n_e==1&&n_mu==0";
  TCut muEnergyCut = "mu_reco_nu_E>0&&mu_reco_nu_E<1250&&mu_recoKE>200&&mu_recoKE<2000";
  TCut eEnergyCut = "e_reco_nu_E>0&&e_reco_nu_E<2500&&e_recoKE>100&&e_recoKE<2500";
  TCut nueCCQE = nueCut && "interaction_mode==1";
  TCut ccmu = numuCut && "interaction_mode==1";

  TFile *o = new TFile("selection_hists2.root", "RECREATE");


    cout << "Making initial hists" << endl;
    int nbins = 60;
    TH2D *hmu = new TH2D("hmu", "hmu", nbins, 0, 300, nbins, 0, 600);
    TH3D *hmuVtx = new TH3D("hmuVtx", "numu vertex", nbins, 0, 300, nbins, 0, 600, 200, 0, 200);
    TH3D *hmuDir = new TH3D("hmuDir", "numu direction", nbins, 0, 300, nbins, 0, 600, 200, 0, 20);
    TH3D *hmuEmu = new TH3D("hmuEmu", "numu Emu", nbins, 0, 300, nbins, 0, 600, 200, -500, 500);
    TH3D *hmuFEmu = new TH3D("hmuFEmu", "numu FEmu", nbins, 0, 300, nbins, 0, 600, 200, -1, 1.5);
    TH3D *hmuEnu = new TH3D("hmuEnu", "numu Enu", nbins, 0, 300, nbins, 0, 600, 200, -800, 800);
    TH3D *hmuFEnu = new TH3D("hmuFEnu", "numu FEnu", nbins, 0, 300, nbins, 0, 600, 200, -1, 2);
    selection->Draw("mu_reco_towall:mu_reco_dwall>>hmu", numuCCQE && mu1R, "goff", 400000, 00000);
    selection->Draw("mu_diffVtx:mu_reco_towall:mu_reco_dwall>>hmuVtx", numuCCQE && mu1R, "goff", 400000, 00000);
    selection->Draw("mu_diffDir*180/TMath::Pi():mu_reco_towall:mu_reco_dwall>>hmuDir", numuCCQE && mu1R, "goff", 400000,
                    00000);
    selection->Draw("mu_diffKE:mu_reco_towall:mu_reco_dwall>>hmuEmu", numuCCQE && mu1R, "goff", 400000, 00000);
    selection->Draw("mu_diffKE/trueKE:mu_reco_towall:mu_reco_dwall>>hmuFEmu", numuCCQE && mu1R, "goff", 400000, 00000);
    selection->Draw("mu_diff_nu_E:mu_reco_towall:mu_reco_dwall>>hmuEnu", numuCCQE && mu1R, "goff", 400000, 00000);
    selection->Draw("mu_diff_nu_E/neutrino_E:mu_reco_towall:mu_reco_dwall>>hmuFEnu", numuCCQE && mu1R, "goff", 400000,
                    00000);
    TH2D *hmuVtxRes = new TH2D("hmuVtxRes", "#nu_{#mu} vertex resolution", nbins + 1, -300 * 0.5 / nbins,
                               300 + 300 * 0.5 / nbins, nbins + 1, -600 * 0.5 / nbins, 600 + 600 * 0.5 / nbins);
    TH2D *hmuDirRes = new TH2D("hmuDirRes", "#nu_{#mu} direction resolution", nbins + 1, -300 * 0.5 / nbins,
                               300 + 300 * 0.5 / nbins, nbins + 1, -600 * 0.5 / nbins, 600 + 600 * 0.5 / nbins);
    TH2D *hmuEnuRes = new TH2D("hmuEnuRes", "#nu_{#mu} neutrino energy resolution [MeV]", nbins + 1, -300 * 0.5 / nbins,
                               300 + 300 * 0.5 / nbins, nbins + 1, -600 * 0.5 / nbins, 600 + 600 * 0.5 / nbins);
    TH2D *hmuFEnuRes = new TH2D("hmuFEnuRes", "#nu_{#mu} neutrino energy resolution [%]", nbins + 1, -300 * 0.5 / nbins,
                                300 + 300 * 0.5 / nbins, nbins + 1, -600 * 0.5 / nbins, 600 + 600 * 0.5 / nbins);
    TH2D *hmuEmuRes = new TH2D("hmuEmuRes", "#nu_{#mu} muon energy resolution [MeV]", nbins + 1, -300 * 0.5 / nbins,
                               300 + 300 * 0.5 / nbins, nbins + 1, -600 * 0.5 / nbins, 600 + 600 * 0.5 / nbins);
    TH2D *hmuFEmuRes = new TH2D("hmuFEmuRes", "#nu_{#mu} muon energy resolution [%]", nbins + 1, -300 * 0.5 / nbins,
                                300 + 300 * 0.5 / nbins, nbins + 1, -600 * 0.5 / nbins, 600 + 600 * 0.5 / nbins);
    TH1D *hmuEmuRes1 = new TH1D("hmuEmuRes1", "hmuEmuRes1", 200, -500, 500);
    TH1D *hmuFEmuRes1 = new TH1D("hmuFEmuRes1", "hmuFEmuRes1", 200, -1, 1.5);
    TH1D *hmuEnuRes1 = new TH1D("hmuEnuRes1", "hmuEnuRes1", 200, -800, 800);
    TH1D *hmuFEnuRes1 = new TH1D("hmuFEnuRes1", "hmuFEnuRes1", 200, -1, 2);
    TF1 *fit = new TF1("fit", "gaus");
    cout << "Making cumulative" << endl;
    double muTot = 0, muVtx = 0, muVtxNew = 0, muDir = 0, muDirNew = 0, muEnu = 0, muEmu = 0, muFEnu = 0, muFEmu = 0, res = 0;
    for (Int_t biny = nbins + 1; biny >= 1; --biny) {
      if (biny % 10 == 0) cout << biny / 10 << endl;
      for (Int_t binx = nbins + 1; binx >= 1; --binx) {
        muTot = hmu->GetBinContent(binx, biny);
        if (biny < nbins + 1) muTot += hmu->GetBinContent(binx, biny + 1);
        if (binx < nbins + 1) muTot += hmu->GetBinContent(binx + 1, biny);
        if (binx < nbins + 1 && biny < nbins + 1) muTot -= hmu->GetBinContent(binx + 1, biny + 1);
        hmu->SetBinContent(binx, biny, muTot);
        muVtx = 0;
        muDir = 0;
        for (int binz = 0; binz <= 201; binz++) {
          muVtxNew = hmuVtx->GetBinContent(binx, biny, binz);
          muDirNew = hmuDir->GetBinContent(binx, biny, binz);
          muEmu = hmuEmu->GetBinContent(binx, biny, binz);
          muFEmu = hmuFEmu->GetBinContent(binx, biny, binz);
          muEnu = hmuEnu->GetBinContent(binx, biny, binz);
          muFEnu = hmuFEnu->GetBinContent(binx, biny, binz);
          if (biny < nbins + 1) {
            muVtxNew += hmuVtx->GetBinContent(binx, biny + 1, binz);
            muDirNew += hmuDir->GetBinContent(binx, biny + 1, binz);
            muEmu += hmuEmu->GetBinContent(binx, biny + 1, binz);
            muEnu += hmuEnu->GetBinContent(binx, biny + 1, binz);
            muFEmu += hmuFEmu->GetBinContent(binx, biny + 1, binz);
            muFEnu += hmuFEnu->GetBinContent(binx, biny + 1, binz);
          }
          if (binx < nbins + 1) {
            muVtxNew += hmuVtx->GetBinContent(binx + 1, biny, binz);
            muDirNew += hmuDir->GetBinContent(binx + 1, biny, binz);
            muEmu += hmuEmu->GetBinContent(binx + 1, biny, binz);
            muEnu += hmuEnu->GetBinContent(binx + 1, biny, binz);
            muFEmu += hmuFEmu->GetBinContent(binx + 1, biny, binz);
            muFEnu += hmuFEnu->GetBinContent(binx + 1, biny, binz);
            if (biny < nbins + 1) {
              muVtxNew -= hmuVtx->GetBinContent(binx + 1, biny + 1, binz);
              muDirNew -= hmuDir->GetBinContent(binx + 1, biny + 1, binz);
              muEmu -= hmuEmu->GetBinContent(binx + 1, biny + 1, binz);
              muEnu -= hmuEnu->GetBinContent(binx + 1, biny + 1, binz);
              muFEmu -= hmuFEmu->GetBinContent(binx + 1, biny + 1, binz);
              muFEnu -= hmuFEnu->GetBinContent(binx + 1, biny + 1, binz);
            }
          }
          hmuVtx->SetBinContent(binx, biny, binz, muVtxNew);
          hmuDir->SetBinContent(binx, biny, binz, muDirNew);
          hmuEmu->SetBinContent(binx, biny, binz, muEmu);
          hmuEnu->SetBinContent(binx, biny, binz, muEnu);
          hmuFEmu->SetBinContent(binx, biny, binz, muFEmu);
          hmuFEnu->SetBinContent(binx, biny, binz, muFEnu);
          if (muVtx + muVtxNew > muTot * 0.68 && hmuVtxRes->GetBinContent(binx, biny) == 0) {
            res = hmuVtx->GetZaxis()->GetBinLowEdge(binz) +
                  hmuVtx->GetZaxis()->GetBinWidth(binz) * (muTot * 0.68 - muVtx) / (muVtxNew);
            //cout << binx << " " << biny << " " << binz << " " << muVtx << " " << muVtxNew << " " << muTot << " " <<res << endl;
            hmuVtxRes->SetBinContent(binx, biny, res);
          }
          muVtx += muVtxNew;
          if (muDir + muDirNew > muTot * 0.68 && hmuDirRes->GetBinContent(binx, biny) == 0) {
            res = hmuDir->GetZaxis()->GetBinLowEdge(binz) +
                  hmuDir->GetZaxis()->GetBinWidth(binz) * (muTot * 0.68 - muDir) / (muDirNew);
            hmuDirRes->SetBinContent(binx, biny, res);
          }
          muDir += muDirNew;
          hmuEmuRes1->SetBinContent(binz, muEmu);
          hmuEnuRes1->SetBinContent(binz, muEnu);
          hmuFEmuRes1->SetBinContent(binz, muFEmu);
          hmuFEnuRes1->SetBinContent(binz, muFEnu);
        }
        hmuEmuRes1->Fit(fit, "QN0");
        hmuEmuRes->SetBinContent(binx, biny, fit->GetParameter(2));
        hmuEnuRes1->Fit(fit, "QN0");
        hmuEnuRes->SetBinContent(binx, biny, fit->GetParameter(2));
        hmuFEmuRes1->Fit(fit, "QN0");
        hmuFEmuRes->SetBinContent(binx, biny, fit->GetParameter(2));
        hmuFEnuRes1->Fit(fit, "QN0");
        hmuFEnuRes->SetBinContent(binx, biny, fit->GetParameter(2));
      }
    }
    cout << hmuEnuRes->GetBinContent(hmuEnuRes->FindBin(100, 0)) << endl;

    o->cd();
    hmuVtxRes->GetXaxis()->SetTitle("dwall cut [cm]");
    hmuVtxRes->GetYaxis()->SetTitle("towall cut [cm]");
    hmuVtxRes->GetXaxis()->SetTitleSize(0.05);
    hmuVtxRes->GetYaxis()->SetTitleSize(0.05);
    hmuVtxRes->GetXaxis()->SetLabelSize(0.04);
    hmuVtxRes->GetYaxis()->SetLabelSize(0.04);
    hmuVtxRes->GetXaxis()->SetTitleOffset(0.85);
    hmuVtxRes->GetYaxis()->SetTitleOffset(0.9);
    hmuVtxRes->Write();
    hmuDirRes->GetXaxis()->SetTitle("dwall cut [cm]");
    hmuDirRes->GetYaxis()->SetTitle("towall cut [cm]");
    hmuDirRes->GetXaxis()->SetTitleSize(0.05);
    hmuDirRes->GetYaxis()->SetTitleSize(0.05);
    hmuDirRes->GetXaxis()->SetLabelSize(0.04);
    hmuDirRes->GetYaxis()->SetLabelSize(0.04);
    hmuDirRes->GetXaxis()->SetTitleOffset(0.85);
    hmuDirRes->GetYaxis()->SetTitleOffset(0.9);
    hmuDirRes->Write();
    hmuEmuRes->GetXaxis()->SetTitle("dwall cut [cm]");
    hmuEmuRes->GetYaxis()->SetTitle("towall cut [cm]");
    hmuEmuRes->GetXaxis()->SetTitleSize(0.05);
    hmuEmuRes->GetYaxis()->SetTitleSize(0.05);
    hmuEmuRes->GetXaxis()->SetLabelSize(0.04);
    hmuEmuRes->GetYaxis()->SetLabelSize(0.04);
    hmuEmuRes->GetXaxis()->SetTitleOffset(0.85);
    hmuEmuRes->GetYaxis()->SetTitleOffset(0.9);
    hmuEmuRes->Write();
    hmuEnuRes->GetXaxis()->SetTitle("dwall cut [cm]");
    hmuEnuRes->GetYaxis()->SetTitle("towall cut [cm]");
    hmuEnuRes->GetXaxis()->SetTitleSize(0.05);
    hmuEnuRes->GetYaxis()->SetTitleSize(0.05);
    hmuEnuRes->GetXaxis()->SetLabelSize(0.04);
    hmuEnuRes->GetYaxis()->SetLabelSize(0.04);
    hmuEnuRes->GetXaxis()->SetTitleOffset(0.85);
    hmuEnuRes->GetYaxis()->SetTitleOffset(0.9);
    hmuEnuRes->Write();
    hmuFEmuRes->GetXaxis()->SetTitle("dwall cut [cm]");
    hmuFEmuRes->GetYaxis()->SetTitle("towall cut [cm]");
    hmuFEmuRes->GetXaxis()->SetTitleSize(0.05);
    hmuFEmuRes->GetYaxis()->SetTitleSize(0.05);
    hmuFEmuRes->GetXaxis()->SetLabelSize(0.04);
    hmuFEmuRes->GetYaxis()->SetLabelSize(0.04);
    hmuFEmuRes->GetXaxis()->SetTitleOffset(0.85);
    hmuFEmuRes->GetYaxis()->SetTitleOffset(0.9);
    hmuFEmuRes->Write();
    hmuFEnuRes->GetXaxis()->SetTitle("dwall cut [cm]");
    hmuFEnuRes->GetYaxis()->SetTitle("towall cut [cm]");
    hmuFEnuRes->GetXaxis()->SetTitleSize(0.05);
    hmuFEnuRes->GetYaxis()->SetTitleSize(0.05);
    hmuFEnuRes->GetXaxis()->SetLabelSize(0.04);
    hmuFEnuRes->GetYaxis()->SetLabelSize(0.04);
    hmuFEnuRes->GetXaxis()->SetTitleOffset(0.85);
    hmuFEnuRes->GetYaxis()->SetTitleOffset(0.9);
    hmuFEnuRes->Write();


    cout << "Making initial hists" << endl;
    TH2D *he = new TH2D("he", "he", nbins, 0, 300, nbins, 0, 600);
    TH3D *heVtx = new TH3D("heVtx", "nue vertex", nbins, 0, 300, nbins, 0, 600, 200, 0, 200);
    TH3D *heDir = new TH3D("heDir", "nue direction", nbins, 0, 300, nbins, 0, 600, 200, 0, 20);
    TH3D *heEe = new TH3D("heEe", "nue Ee", nbins, 0, 300, nbins, 0, 600, 200, -500, 500);
    TH3D *heFEe = new TH3D("heFEe", "nue FEe", nbins, 0, 300, nbins, 0, 600, 200, -1, 1.5);
    TH3D *heEnu = new TH3D("heEnu", "nue Enu", nbins, 0, 300, nbins, 0, 600, 200, -800, 800);
    TH3D *heFEnu = new TH3D("heFEnu", "nue FEnu", nbins, 0, 300, nbins, 0, 600, 200, -1, 2);
    selection->Draw("e_reco_towall:e_reco_dwall>>he", nueCCQE && e1R, "goff", 400000, 00000);
    selection->Draw("e_diffVtx:e_reco_towall:e_reco_dwall>>heVtx", nueCCQE && e1R, "goff", 400000, 00000);
    selection->Draw("e_diffDir*180/TMath::Pi():e_reco_towall:e_reco_dwall>>heDir", nueCCQE && e1R, "goff", 400000,
                    00000);
    selection->Draw("e_diffKE:e_reco_towall:e_reco_dwall>>heEe", nueCCQE && e1R, "goff", 400000, 00000);
    selection->Draw("e_diffKE/trueKE:e_reco_towall:e_reco_dwall>>heFEe", nueCCQE && e1R, "goff", 400000, 00000);
    selection->Draw("e_diff_nu_E:e_reco_towall:e_reco_dwall>>heEnu", nueCCQE && e1R, "goff", 400000, 00000);
    selection->Draw("e_diff_nu_E/neutrino_E:e_reco_towall:e_reco_dwall>>heFEnu", nueCCQE && e1R, "goff", 400000, 00000);
    TH2D *heVtxRes = new TH2D("heVtxRes", "#nu_{e} vertex resolution", nbins + 1, -0.5, 300.5, nbins + 1, -1, 601);
    TH2D *heDirRes = new TH2D("heDirRes", "#nu_{e} direction resolution", nbins + 1, -0.5, 300.5, nbins + 1, -1, 601);
    TH2D *heEeRes = new TH2D("heEeRes", "#nu_{e} electron energy resolution [MeV]", nbins + 1, -0.5, 300.5, nbins + 1,
                             -1,
                             601);
    TH2D *heFEeRes = new TH2D("heFEeRes", "#nu_{e} electron energy resolution [%]", nbins + 1, -0.5, 300.5, nbins + 1,
                              -1,
                              601);
    TH2D *heEnuRes = new TH2D("heEnuRes", "#nu_{e} neutrino energy resolution [MeV]", nbins + 1, -0.5, 300.5, nbins + 1,
                              -1, 601);
    TH2D *heFEnuRes = new TH2D("heFEnuRes", "#nu_{e} neutrino energy resolution [%]", nbins + 1, -0.5, 300.5, nbins + 1,
                               -1, 601);
    TH1D *heEeRes1 = new TH1D("heEeRes1", "heEeRes1", 200, -500, 500);
    TH1D *heFEeRes1 = new TH1D("heFEeRes1", "heFEeRes1", 200, -1, 1.5);
    TH1D *heEnuRes1 = new TH1D("heEnuRes1", "heEnuRes1", 200, -800, 800);
    TH1D *heFEnuRes1 = new TH1D("heFEnuRes1", "heFEnuRes1", 200, -1, 2);
    fit = new TF1("fit", "gaus");
    cout << "Making cumulative" << endl;
    double eTot = 0, eVtx = 0, eVtxNew = 0, eDir = 0, eDirNew = 0, eEnu = 0, eEe = 0, eFEnu = 0, eFEe = 0;
    for (Int_t biny = nbins + 1; biny >= 1; --biny) {
      if (biny % 10 == 0) cout << biny / 10 << endl;
      for (Int_t binx = nbins + 1; binx >= 1; --binx) {
        eTot = he->GetBinContent(binx, biny);
        if (biny < nbins + 1) eTot += he->GetBinContent(binx, biny + 1);
        if (binx < nbins + 1) eTot += he->GetBinContent(binx + 1, biny);
        if (binx < nbins + 1 && biny < nbins + 1) eTot -= he->GetBinContent(binx + 1, biny + 1);
        he->SetBinContent(binx, biny, eTot);
        eVtx = 0;
        eDir = 0;
        for (int binz = 0; binz <= 201; binz++) {
          eVtxNew = heVtx->GetBinContent(binx, biny, binz);
          eDirNew = heDir->GetBinContent(binx, biny, binz);
          eEe = heEe->GetBinContent(binx, biny, binz);
          eFEe = heFEe->GetBinContent(binx, biny, binz);
          eEnu = heEnu->GetBinContent(binx, biny, binz);
          eFEnu = heFEnu->GetBinContent(binx, biny, binz);
          if (biny < nbins + 1) {
            eVtxNew += heVtx->GetBinContent(binx, biny + 1, binz);
            eDirNew += heDir->GetBinContent(binx, biny + 1, binz);
            eEe += heEe->GetBinContent(binx, biny + 1, binz);
            eEnu += heEnu->GetBinContent(binx, biny + 1, binz);
            eFEe += heFEe->GetBinContent(binx, biny + 1, binz);
            eFEnu += heFEnu->GetBinContent(binx, biny + 1, binz);
          }
          if (binx < nbins + 1) {
            eVtxNew += heVtx->GetBinContent(binx + 1, biny, binz);
            eDirNew += heDir->GetBinContent(binx + 1, biny, binz);
            eEe += heEe->GetBinContent(binx + 1, biny, binz);
            eEnu += heEnu->GetBinContent(binx + 1, biny, binz);
            eFEe += heFEe->GetBinContent(binx + 1, biny, binz);
            eFEnu += heFEnu->GetBinContent(binx + 1, biny, binz);
            if (biny < nbins + 1) {
              eVtxNew -= heVtx->GetBinContent(binx + 1, biny + 1, binz);
              eDirNew -= heDir->GetBinContent(binx + 1, biny + 1, binz);
              eEe -= heEe->GetBinContent(binx + 1, biny + 1, binz);
              eEnu -= heEnu->GetBinContent(binx + 1, biny + 1, binz);
              eFEe -= heFEe->GetBinContent(binx + 1, biny + 1, binz);
              eFEnu -= heFEnu->GetBinContent(binx + 1, biny + 1, binz);
            }
          }
          heVtx->SetBinContent(binx, biny, binz, eVtxNew);
          heDir->SetBinContent(binx, biny, binz, eDirNew);
          heEe->SetBinContent(binx, biny, binz, eEe);
          heEnu->SetBinContent(binx, biny, binz, eEnu);
          heFEe->SetBinContent(binx, biny, binz, eFEe);
          heFEnu->SetBinContent(binx, biny, binz, eFEnu);
          if (eVtx + eVtxNew > eTot * 0.68 && heVtxRes->GetBinContent(binx, biny) == 0) {
            res = heVtx->GetZaxis()->GetBinLowEdge(binz) +
                  heVtx->GetZaxis()->GetBinWidth(binz) * (eTot * 0.68 - eVtx) / (eVtxNew);
            //cout << binx << " " << biny << " " << binz << " " << eVtx << " " << eVtxNew << " " << eTot << " " <<res << endl;
            heVtxRes->SetBinContent(binx, biny, res);
          }
          eVtx += eVtxNew;
          if (eDir + eDirNew > eTot * 0.68 && heDirRes->GetBinContent(binx, biny) == 0) {
            res = heDir->GetZaxis()->GetBinLowEdge(binz) +
                  heDir->GetZaxis()->GetBinWidth(binz) * (eTot * 0.68 - eDir) / (eDirNew);
            heDirRes->SetBinContent(binx, biny, res);
          }
          eDir += eDirNew;
          heEeRes1->SetBinContent(binz, eEe);
          heEnuRes1->SetBinContent(binz, eEnu);
          heFEeRes1->SetBinContent(binz, eFEe);
          heFEnuRes1->SetBinContent(binz, eFEnu);
        }
        heEeRes1->Fit(fit, "QN0");
        heEeRes->SetBinContent(binx, biny, fit->GetParameter(2));
        heEnuRes1->Fit(fit, "QN0");
        heEnuRes->SetBinContent(binx, biny, fit->GetParameter(2));
        heFEeRes1->Fit(fit, "QN0");
        heFEeRes->SetBinContent(binx, biny, fit->GetParameter(2));
        heFEnuRes1->Fit(fit, "QN0");
        heFEnuRes->SetBinContent(binx, biny, fit->GetParameter(2));
      }
    }
    cout << heEnuRes->GetBinContent(hmuEnuRes->FindBin(100, 0)) << endl;

    o->cd();
    heVtxRes->GetXaxis()->SetTitle("dwall cut [cm]");
    heVtxRes->GetYaxis()->SetTitle("towall cut [cm]");
    heVtxRes->GetXaxis()->SetTitleSize(0.05);
    heVtxRes->GetYaxis()->SetTitleSize(0.05);
    heVtxRes->GetXaxis()->SetLabelSize(0.04);
    heVtxRes->GetYaxis()->SetLabelSize(0.04);
    heVtxRes->GetXaxis()->SetTitleOffset(0.85);
    heVtxRes->GetYaxis()->SetTitleOffset(0.9);
    heVtxRes->Write();
    heDirRes->GetXaxis()->SetTitle("dwall cut [cm]");
    heDirRes->GetYaxis()->SetTitle("towall cut [cm]");
    heDirRes->GetXaxis()->SetTitleSize(0.05);
    heDirRes->GetYaxis()->SetTitleSize(0.05);
    heDirRes->GetXaxis()->SetLabelSize(0.04);
    heDirRes->GetYaxis()->SetLabelSize(0.04);
    heDirRes->GetXaxis()->SetTitleOffset(0.85);
    heDirRes->GetYaxis()->SetTitleOffset(0.9);
    heDirRes->Write();
    heEeRes->GetXaxis()->SetTitle("dwall cut [cm]");
    heEeRes->GetYaxis()->SetTitle("towall cut [cm]");
    heEeRes->GetXaxis()->SetTitleSize(0.05);
    heEeRes->GetYaxis()->SetTitleSize(0.05);
    heEeRes->GetXaxis()->SetLabelSize(0.04);
    heEeRes->GetYaxis()->SetLabelSize(0.04);
    heEeRes->GetXaxis()->SetTitleOffset(0.85);
    heEeRes->GetYaxis()->SetTitleOffset(0.9);
    heEeRes->Write();
    heEnuRes->GetXaxis()->SetTitle("dwall cut [cm]");
    heEnuRes->GetYaxis()->SetTitle("towall cut [cm]");
    heEnuRes->GetXaxis()->SetTitleSize(0.05);
    heEnuRes->GetYaxis()->SetTitleSize(0.05);
    heEnuRes->GetXaxis()->SetLabelSize(0.04);
    heEnuRes->GetYaxis()->SetLabelSize(0.04);
    heEnuRes->GetXaxis()->SetTitleOffset(0.85);
    heEnuRes->GetYaxis()->SetTitleOffset(0.9);
    heEnuRes->Write();
    heFEeRes->GetXaxis()->SetTitle("dwall cut [cm]");
    heFEeRes->GetYaxis()->SetTitle("towall cut [cm]");
    heFEeRes->GetXaxis()->SetTitleSize(0.05);
    heFEeRes->GetYaxis()->SetTitleSize(0.05);
    heFEeRes->GetXaxis()->SetLabelSize(0.04);
    heFEeRes->GetYaxis()->SetLabelSize(0.04);
    heFEeRes->GetXaxis()->SetTitleOffset(0.85);
    heFEeRes->GetYaxis()->SetTitleOffset(0.9);
    heFEeRes->Write();
    heFEnuRes->GetXaxis()->SetTitle("dwall cut [cm]");
    heFEnuRes->GetYaxis()->SetTitle("towall cut [cm]");
    heFEnuRes->GetXaxis()->SetTitleSize(0.05);
    heFEnuRes->GetYaxis()->SetTitleSize(0.05);
    heFEnuRes->GetXaxis()->SetLabelSize(0.04);
    heFEnuRes->GetYaxis()->SetLabelSize(0.04);
    heFEnuRes->GetXaxis()->SetTitleOffset(0.85);
    heFEnuRes->GetYaxis()->SetTitleOffset(0.9);
    heFEnuRes->Write();

  
  
  /*
  selection->Draw("mu_diffVtx:reco_nu_E>>vtx_Enu",numuCCQE&&mu1R&&"reco_nu_E>0&&reco_nu_E<5000&&recoEnergy[0]<2999","goff",400000,00000);
  new TCanvas(); vtx_Enu->ProfileX()->Draw();
  selection->Draw("mu_diffVtx/neutrino_E:mu_recoKE>>vtx_Emu",numuCCQE&&mu1R&&"reco_nu_E>0&&reco_nu_E<5000&&recoEnergy[0]<2999","goff",400000,00000);
  new TCanvas(); vtx_Emu->ProfileX()->Draw();
  selection->Draw("e_diffVtx/neutrino_E:reco_nu_E>>vtx_Enu_e",nueCCQE&&e1R&&"reco_nu_E>0&&reco_nu_E<5000&&recoEnergy[0]<2999","goff",400000,00000);
  new TCanvas(); vtx_Enu_e->ProfileX()->Draw();
  selection->Draw("e_diffVtx/neutrino_E:e_recoKE>>vtx_Ee",nueCCQE&&e1R&&"reco_nu_E>0&&reco_nu_E<5000&&recoEnergy[0]<2999","goff",400000,00000);
  new TCanvas(); vtx_Ee->ProfileX()->Draw();
  */
/**/
  TH2D * dWallmuPIDmu = new TH2D("dWallmuPIDmu","dWallmuPIDmu",1400,-900,500,400,0,400);
  TH2D * dWallmuPIDe = new TH2D("dWallmuPIDe","dWallmuPIDe",1400,-900,500,400,0,400);
  double totalmu = selection->Draw("mu_reco_dwall:recoPIDLikelihood>>dWallmuPIDmu",numuCut&&R1,"goff",400000,00000);
  double totale = selection->Draw("mu_reco_dwall:recoPIDLikelihood>>dWallmuPIDe",nueCut&&R1,"goff",400000,00000);
  TH2D * dWallmuPIDpur   = new TH2D("dWallmuPIDpur",  "#nu_{#mu} sample purity",1401,-900.5,500.5,401,-0.5,400.5);
  TH2D * dWallmuPIDeff   = new TH2D("dWallmuPIDeff",  "#nu_{#mu} sample efficiency",1401,-900.5,500.5,401,-0.5,400.5);
  TH2D * dWallmuPIDmisID = new TH2D("dWallmuPIDmisID","#nu_{#mu} sample #nu_{e} mis-ID rate",1401,-900.5,500.5,401,-0.5,400.5);
  TH2D * dWallmuPID      = new TH2D("dWallmuPID",     "#nu_{#mu} sample purity*efficiency",1401,-900.5,500.5,401,-0.5,400.5);
  double cummu =0, cume =0;
  //First row
  for (Int_t binx = 1401, biny = 401; binx >= 1; --binx) {
    const Int_t bin = dWallmuPIDmu->GetBin(binx, biny);
    cummu += dWallmuPIDmu->GetBinContent(bin);
    cume += dWallmuPIDe->GetBinContent(bin);
    dWallmuPIDmu->SetBinContent(bin, cummu);
    dWallmuPIDe->SetBinContent(bin, cume);
    double eff = cummu/totalmu;
    double pur = cummu*4803420./(cummu*4803420.+cume*91393.6);
    double misID = cume/totale;
    dWallmuPID->SetBinContent(binx,biny,eff*pur);
    dWallmuPIDeff->SetBinContent(binx,biny,eff);
    dWallmuPIDpur->SetBinContent(binx,biny,pur);
    dWallmuPIDmisID->SetBinContent(binx,biny,misID);
  }
  //Subsequent rows
  for (Int_t biny = 400; biny >= 1; --biny) {
    cummu =0; cume =0;
    for (Int_t binx = 1401; binx >= 1; --binx) {
      const Int_t bin = dWallmuPIDmu->GetBin(binx, biny);
      cummu +=dWallmuPIDmu->GetBinContent(bin);
      cume +=dWallmuPIDe->GetBinContent(bin);
      Double_t mu = dWallmuPIDmu->GetBinContent(binx,biny+1)+ cummu;
      dWallmuPIDmu->SetBinContent(bin, mu);
      Double_t e = dWallmuPIDe->GetBinContent(binx,biny+1)+ cume;
      dWallmuPIDe->SetBinContent(bin, e);
      double eff = mu/totalmu;
      double pur = mu*4803420./(mu*4803420.+e*91393.6);
      double misID = e/totale;
      dWallmuPID->SetBinContent(binx,biny,eff*pur);
      dWallmuPIDeff->SetBinContent(binx,biny,eff);
      dWallmuPIDpur->SetBinContent(binx,biny,pur);
      dWallmuPIDmisID->SetBinContent(binx,biny,misID);
    }
  }
  TH2D * dWallePIDmu = new TH2D("dWallePIDmu","dWallePIDmu",1400,-900,500,400,0,400);
  TH2D * dWallePIDe  = new TH2D("dWallePIDe","dWallePIDe",1400,-900,500,400,0,400);
  totalmu = selection->Draw("e_reco_dwall:recoPIDLikelihood>>dWallePIDmu",numuCCQE&&R1,"goff",400000,00000);
  totale  = selection->Draw("e_reco_dwall:recoPIDLikelihood>>dWallePIDe" ,nueCCQE&&R1 ,"goff",400000,00000);
  TH2D * dWallePIDpur   = new TH2D("dWallePIDpur",  "#nu_{e} sample purity",1401,-900.5,500.5,401,-0.5,400.5);
  TH2D * dWallePIDeff   = new TH2D("dWallePIDeff",  "#nu_{e} sample efficiency",1401,-900.5,500.5,401,-0.5,400.5);
  TH2D * dWallePIDmisID = new TH2D("dWallePIDmisID","#nu_{e} sample #nu_{#mu} mis-ID rate",1401,-900.5,500.5,401,-0.5,400.5);
  TH2D * dWallePID      = new TH2D("dWallePID",     "#nu_{e} sample purity*efficiency",1401,-900.5,500.5,401,-0.5,400.5);
  cummu=0; cume=0;
  //First row
  for (Int_t binx = 0, biny = 401; binx <= 1400; ++binx) {
    const Int_t bin = dWallePIDmu->GetBin(binx, biny);
    cummu += dWallePIDmu->GetBinContent(bin);
    cume += dWallePIDe->GetBinContent(bin);
    dWallePIDmu->SetBinContent(bin, cummu);
    dWallePIDe->SetBinContent(bin, cume);
    double eff = cume/totale;
    double pur = cume*91393.6/(cummu*4803420.+cume*91393.6);
    double misID = cummu/totalmu;
    dWallePID->SetBinContent(binx+1,biny,eff*pur);
    dWallePIDeff->SetBinContent(binx,biny,eff);
    dWallePIDpur->SetBinContent(binx,biny,pur);
    dWallePIDmisID->SetBinContent(binx,biny, misID);
  }
  //Subsequent rows
  for (Int_t biny = 400; biny >= 1; --biny) {
    cummu =0; cume =0;
    for (Int_t binx = 0; binx <= 1400; ++binx) {
      const Int_t bin = dWallePIDmu->GetBin(binx, biny);
      cummu +=dWallePIDmu->GetBinContent(bin);
      cume +=dWallePIDe->GetBinContent(bin);
      Double_t mu = dWallePIDmu->GetBinContent(binx,biny+1)+ cummu;
      dWallePIDmu->SetBinContent(bin, mu);
      Double_t e = dWallePIDe->GetBinContent(binx,biny+1)+ cume;
      dWallePIDe->SetBinContent(bin,e);
      double eff = e/totale;
      double pur = e*91393.6/(mu*4803420.+e*91393.6);
      double misID = mu/totalmu;
      dWallePID->SetBinContent(binx+1,biny,eff*pur);
      dWallePIDeff->SetBinContent(binx,biny,eff);
      dWallePIDpur->SetBinContent(binx,biny,pur);
      dWallePIDmisID->SetBinContent(binx,biny, misID);
    }
  }
  TCanvas *cmu = new TCanvas();
  cmu->Divide(2,2);
  cmu->cd(1);
  dWallmuPIDeff->GetXaxis()->SetTitle("PID cut position");
  dWallmuPIDeff->GetYaxis()->SetTitle("Distance from wall cut [cm]");
  dWallmuPIDeff->Draw("colz");
  dWallmuPIDeff->GetXaxis()->SetTitleSize(0.05);
  dWallmuPIDeff->GetYaxis()->SetTitleSize(0.05);
  dWallmuPIDeff->GetXaxis()->SetLabelSize(0.04);
  dWallmuPIDeff->GetYaxis()->SetLabelSize(0.04);
  dWallmuPIDeff->GetXaxis()->SetTitleOffset(0.85);
  dWallmuPIDeff->GetYaxis()->SetTitleOffset(0.9);
  cmu->cd(2);
  dWallmuPIDpur->GetXaxis()->SetTitle("PID cut position");
  dWallmuPIDpur->GetYaxis()->SetTitle("Distance from wall cut [cm]");
  dWallmuPIDpur->GetXaxis()->SetTitleSize(0.05);
  dWallmuPIDpur->GetYaxis()->SetTitleSize(0.05);
  dWallmuPIDpur->GetXaxis()->SetLabelSize(0.04);
  dWallmuPIDpur->GetYaxis()->SetLabelSize(0.04);
  dWallmuPIDpur->GetXaxis()->SetTitleOffset(0.85);
  dWallmuPIDpur->GetYaxis()->SetTitleOffset(0.9);
  dWallmuPIDpur->Draw("colz");
  cmu->cd(3);
  dWallmuPIDmisID->GetXaxis()->SetTitle("PID cut position");
  dWallmuPIDmisID->GetYaxis()->SetTitle("Distance from wall cut [cm]");
  dWallmuPIDmisID->GetXaxis()->SetTitleSize(0.05);
  dWallmuPIDmisID->GetYaxis()->SetTitleSize(0.05);
  dWallmuPIDmisID->GetXaxis()->SetLabelSize(0.04);
  dWallmuPIDmisID->GetYaxis()->SetLabelSize(0.04);
  dWallmuPIDmisID->GetXaxis()->SetTitleOffset(0.85);
  dWallmuPIDmisID->GetYaxis()->SetTitleOffset(0.9);
  dWallmuPIDmisID->Draw("colz");
  cmu->SetLogz();
  cmu->cd(4);
  dWallmuPID->GetXaxis()->SetTitle("PID cut position");
  dWallmuPID->GetYaxis()->SetTitle("Distance from wall cut [cm]");
  dWallmuPID->GetXaxis()->SetTitleSize(0.05);
  dWallmuPID->GetYaxis()->SetTitleSize(0.05);
  dWallmuPID->GetXaxis()->SetLabelSize(0.04);
  dWallmuPID->GetYaxis()->SetLabelSize(0.04);
  dWallmuPID->GetXaxis()->SetTitleOffset(0.85);
  dWallmuPID->GetYaxis()->SetTitleOffset(0.9);
  dWallmuPID->Draw("colz");
  
  TCanvas *ce = new TCanvas();
  ce->Divide(2,2);
  ce->cd(1);
  dWallePIDeff->GetXaxis()->SetTitle("PID cut position");
  dWallePIDeff->GetYaxis()->SetTitle("Distance from wall cut [cm]");
  dWallePIDeff->GetXaxis()->SetTitleSize(0.05);
  dWallePIDeff->GetYaxis()->SetTitleSize(0.05);
  dWallePIDeff->GetXaxis()->SetLabelSize(0.04);
  dWallePIDeff->GetYaxis()->SetLabelSize(0.04);
  dWallePIDeff->GetXaxis()->SetTitleOffset(0.85);
  dWallePIDeff->GetYaxis()->SetTitleOffset(0.9);
  dWallePIDeff->Draw("colz");
  ce->cd(2);
  dWallePIDpur->GetXaxis()->SetTitle("PID cut position");
  dWallePIDpur->GetYaxis()->SetTitle("Distance from wall cut [cm]");
  dWallePIDpur->GetXaxis()->SetTitleSize(0.05);
  dWallePIDpur->GetYaxis()->SetTitleSize(0.05);
  dWallePIDpur->GetXaxis()->SetLabelSize(0.04);
  dWallePIDpur->GetYaxis()->SetLabelSize(0.04);
  dWallePIDpur->GetXaxis()->SetTitleOffset(0.85);
  dWallePIDpur->GetYaxis()->SetTitleOffset(0.9);
  dWallePIDpur->Draw("colz");
  ce->cd(3);
  dWallePIDmisID->GetXaxis()->SetTitle("PID cut position");
  dWallePIDmisID->GetYaxis()->SetTitle("Distance from wall cut [cm]");
  dWallePIDmisID->GetXaxis()->SetTitleSize(0.05);
  dWallePIDmisID->GetYaxis()->SetTitleSize(0.05);
  dWallePIDmisID->GetXaxis()->SetLabelSize(0.04);
  dWallePIDmisID->GetYaxis()->SetLabelSize(0.04);
  dWallePIDmisID->GetXaxis()->SetTitleOffset(0.85);
  dWallePIDmisID->GetYaxis()->SetTitleOffset(0.9);
  dWallePIDmisID->Draw("colz");
  ce->SetLogz();
  ce->cd(4);
  dWallePID->GetXaxis()->SetTitle("PID cut position");
  dWallePID->GetYaxis()->SetTitle("Distance from wall cut [cm]");
  dWallePID->GetXaxis()->SetTitleSize(0.05);
  dWallePID->GetYaxis()->SetTitleSize(0.05);
  dWallePID->GetXaxis()->SetLabelSize(0.04);
  dWallePID->GetYaxis()->SetLabelSize(0.04);
  dWallePID->GetXaxis()->SetTitleOffset(0.85);
  dWallePID->GetYaxis()->SetTitleOffset(0.9);
  dWallePID->Draw("colz");
  o->cd();
  dWallmuPIDeff->Write();
  dWallmuPIDpur->Write();
  dWallmuPIDmisID->Write();
  dWallmuPID->Write();
  dWallePIDeff->Write();
  dWallePIDpur->Write();
  dWallePIDmisID->Write();
  dWallePID->Write();
  cmu->Write();
  ce->Write();
  cmu->Close();
  ce->Close();
/*
  TH2D * toWallmuPIDmu = new TH2D("toWallmuPIDmu","toWallmuPIDmu",1000,-500,500,400,0,400);
  TH2D * toWallmuPIDe = new TH2D("toWallmuPIDe","toWallmuPIDe",1000,-500,500,400,0,400);
  totalmu = selection->Draw("mu_reco_towall:recoPIDLikelihood>>toWallmuPIDmu",numuCCQE&&R1&&muEnergyCut,"goff",400000,00000);
  totale = selection->Draw("mu_reco_towall:recoPIDLikelihood>>toWallmuPIDe",nueCCQE&&R1&&muEnergyCut,"goff",400000,00000);
  TH2D * toWallmuPIDpur   = new TH2D("toWallmuPIDpur",  "toWallmuPIDpur",1001,-500.5,500.5,401,-0.5,400.5);
  TH2D * toWallmuPIDeff   = new TH2D("toWallmuPIDeff",  "toWallmuPIDeff",1001,-500.5,500.5,401,-0.5,400.5);
  TH2D * toWallmuPIDmisID = new TH2D("toWallmuPIDmisID","toWallmuPIDmisID",1001,-500.5,500.5,401,-0.5,400.5);
  TH2D * toWallmuPID      = new TH2D("toWallmuPID",     "toWallmuPID",1001,-500.5,500.5,401,-0.5,400.5);
  cummu =0; cume =0;
  //First row
  for (Int_t binx = 1001, biny = 401; binx >= 1; --binx) {
    const Int_t bin = toWallmuPIDmu->GetBin(binx, biny);
    cummu += toWallmuPIDmu->GetBinContent(bin);
    cume += toWallmuPIDe->GetBinContent(bin);
    toWallmuPIDmu->SetBinContent(bin, cummu);
    toWallmuPIDe->SetBinContent(bin, cume);
    double eff = cummu/totalmu;
    double pur = cummu/(cummu+cume);
    double misID = cume/totale;
    toWallmuPID->SetBinContent(binx,biny,eff*pur);
    toWallmuPIDeff->SetBinContent(binx,biny,eff);
    toWallmuPIDpur->SetBinContent(binx,biny,pur);
    toWallmuPIDmisID->SetBinContent(binx,biny,misID);
  }
  //Subsequent rows
  for (Int_t biny = 400; biny >= 1; --biny) {
    cummu =0; cume =0;
    for (Int_t binx = 1001; binx >= 1; --binx) {
      const Int_t bin = toWallmuPIDmu->GetBin(binx, biny);
      cummu +=toWallmuPIDmu->GetBinContent(bin);
      cume +=toWallmuPIDe->GetBinContent(bin);
      Double_t mu = toWallmuPIDmu->GetBinContent(binx,biny+1)+ cummu;
      toWallmuPIDmu->SetBinContent(bin, mu);
      Double_t e = toWallmuPIDe->GetBinContent(binx,biny+1)+ cume;
      toWallmuPIDe->SetBinContent(bin, e);
      double eff = mu/totalmu;
      double pur = mu/(mu+e);
      double misID = e/totale;
      toWallmuPID->SetBinContent(binx,biny,eff*pur);
      toWallmuPIDeff->SetBinContent(binx,biny,eff);
      toWallmuPIDpur->SetBinContent(binx,biny,pur);
      toWallmuPIDmisID->SetBinContent(binx,biny,misID);
    }
  }
  TH2D * toWallePIDmu = new TH2D("toWallePIDmu","toWallePIDmu",1000,-500,500,400,0,400);
  TH2D * toWallePIDe  = new TH2D("toWallePIDe","toWallePIDe",1000,-500,500,400,0,400);
  totalmu = selection->Draw("e_reco_towall:recoPIDLikelihood>>toWallePIDmu",numuCCQE&&R1&&eEnergyCut,"goff",400000,00000);
  totale  = selection->Draw("e_reco_towall:recoPIDLikelihood>>toWallePIDe" ,nueCCQE&&R1&&eEnergyCut ,"goff",400000,00000);
  TH2D * toWallePIDpur   = new TH2D("toWallePIDpur",  "toWallePIDpur",1001,-500.5,500.5,401,-0.5,400.5);
  TH2D * toWallePIDeff   = new TH2D("toWallePIDeff",  "toWallePIDeff",1001,-500.5,500.5,401,-0.5,400.5);
  TH2D * toWallePIDmisID = new TH2D("toWallePIDmisID","toWallePIDmisID",1001,-500.5,500.5,401,-0.5,400.5);
  TH2D * toWallePID = new TH2D("toWallePID","toWallePID",1001,-500.5,500.5,401,-0.5,400.5);
  cummu=0; cume=0;
  //First row
  for (Int_t binx = 0, biny = 401; binx <= 1000; ++binx) {
    const Int_t bin = toWallePIDmu->GetBin(binx, biny);
    cummu += toWallePIDmu->GetBinContent(bin);
    cume += toWallePIDe->GetBinContent(bin);
    toWallePIDmu->SetBinContent(bin, cummu);
    toWallePIDe->SetBinContent(bin, cume);
    double eff = cume/totale;
    double pur = cume/(cummu+cume);
    double misID = cummu/totalmu;
    toWallePID->SetBinContent(binx+1,biny,eff*pur);
    toWallePIDeff->SetBinContent(binx,biny,eff);
    toWallePIDpur->SetBinContent(binx,biny,pur);
    toWallePIDmisID->SetBinContent(binx,biny, misID);
  }
  //Subsequent rows
  for (Int_t biny = 400; biny >= 1; --biny) {
    cummu =0; cume =0;
    for (Int_t binx = 0; binx <= 1000; ++binx) {
      const Int_t bin = toWallePIDmu->GetBin(binx, biny);
      cummu +=toWallePIDmu->GetBinContent(bin);
      cume +=toWallePIDe->GetBinContent(bin);
      Double_t mu = toWallePIDmu->GetBinContent(binx,biny+1)+ cummu;
      toWallePIDmu->SetBinContent(bin, mu);
      Double_t e = toWallePIDe->GetBinContent(binx,biny+1)+ cume;
      toWallePIDe->SetBinContent(bin,e);
      double eff = e/totale;
      double pur = e/(mu+e);
      double misID = mu/totalmu;
      toWallePID->SetBinContent(binx+1,biny,eff*pur);
      toWallePIDeff->SetBinContent(binx,biny,eff);
      toWallePIDpur->SetBinContent(binx,biny,pur);
      toWallePIDmisID->SetBinContent(binx,biny, misID);
    }
  }
  /**/
  //gROOT->SetBatch(b);

}
