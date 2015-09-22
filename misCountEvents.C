

void misCountEvents() {

    TFile sfile("selection_test_.root", "READ");
    TTree *selection = (TTree *) sfile.Get("selection");

    Int_t           evt;
    Int_t           neutrino_id;
    Double_t        neutrino_E;
    Double_t        trueKE;
    Int_t           interaction_mode;
    Int_t           n_chkv_part;
    Int_t           n_pi0;
    Int_t           n_e;
    Int_t           n_mu;
    Int_t           n_pi_pm;
    Int_t           n_k0;
    Int_t           n_k_pm;
    Int_t           nneutron;
    Int_t           nproton;
    Int_t           ngamma;
    Int_t           nother;
    Int_t           neutroncount;
    Int_t           n_captures;
    Double_t        true_dwall;
    Double_t        true_towall;
    Double_t        mu_reco_dwall;
    Double_t        mu_reco_towall;
    Double_t        mu_reco_nu_E;
    Double_t        mu_recoKE;
    Double_t        mu_diffVtx;
    Double_t        mu_diffDir;
    Double_t        mu_diffKE;
    Double_t        mu_diffVtxLong;
    Double_t        mu_diffVtxTrans;
    Double_t        mu_diff_nu_E;
    Double_t        e_reco_dwall;
    Double_t        e_reco_towall;
    Double_t        e_reco_nu_E;
    Double_t        e_recoKE;
    Double_t        e_diffVtx;
    Double_t        e_diffDir;
    Double_t        e_diffKE;
    Double_t        e_diffVtxLong;
    Double_t        e_diffVtxTrans;
    Double_t        e_diff_nu_E;
    Int_t           nClusters;
    Int_t           ring1PEs;
    Int_t           ring2PEs;
    Int_t           recoPID;
    Double_t        recoPIDLikelihood;
    Bool_t          isHighE;
    Int_t           nSubevents;
    Double_t        recoTime[20];   //[nSubevents]
    Double_t        recoEnergy[20];   //[nSubevents]
    Int_t           recoCaptures;

    selection->SetBranchAddress("evt", &evt);
    selection->SetBranchAddress("neutrino_id", &neutrino_id);
    selection->SetBranchAddress("neutrino_E", &neutrino_E);
    selection->SetBranchAddress("trueKE", &trueKE);
    selection->SetBranchAddress("interaction_mode", &interaction_mode);
    selection->SetBranchAddress("n_chkv_part", &n_chkv_part);
    selection->SetBranchAddress("n_pi0", &n_pi0);
    selection->SetBranchAddress("n_e", &n_e);
    selection->SetBranchAddress("n_mu", &n_mu);
    selection->SetBranchAddress("n_pi_pm", &n_pi_pm);
    selection->SetBranchAddress("n_k0", &n_k0);
    selection->SetBranchAddress("n_k_pm", &n_k_pm);
    selection->SetBranchAddress("nneutron", &nneutron);
    selection->SetBranchAddress("nproton", &nproton);
    selection->SetBranchAddress("ngamma", &ngamma);
    selection->SetBranchAddress("nother", &nother);
    selection->SetBranchAddress("neutroncount", &neutroncount);
    selection->SetBranchAddress("n_captures", &n_captures);
    selection->SetBranchAddress("true_dwall", &true_dwall);
    selection->SetBranchAddress("true_towall", &true_towall);
    selection->SetBranchAddress("mu_reco_dwall", &mu_reco_dwall);
    selection->SetBranchAddress("mu_reco_towall", &mu_reco_towall);
    selection->SetBranchAddress("mu_reco_nu_E", &mu_reco_nu_E);
    selection->SetBranchAddress("mu_recoKE", &mu_recoKE);
    selection->SetBranchAddress("mu_diffVtx", &mu_diffVtx);
    selection->SetBranchAddress("mu_diffDir", &mu_diffDir);
    selection->SetBranchAddress("mu_diffKE", &mu_diffKE);
    selection->SetBranchAddress("mu_diffVtxLong", &mu_diffVtxLong);
    selection->SetBranchAddress("mu_diffVtxTrans", &mu_diffVtxTrans);
    selection->SetBranchAddress("mu_diff_nu_E", &mu_diff_nu_E);
    selection->SetBranchAddress("e_reco_dwall", &e_reco_dwall);
    selection->SetBranchAddress("e_reco_towall", &e_reco_towall);
    selection->SetBranchAddress("e_reco_nu_E", &e_reco_nu_E);
    selection->SetBranchAddress("e_recoKE", &e_recoKE);
    selection->SetBranchAddress("e_diffVtx", &e_diffVtx);
    selection->SetBranchAddress("e_diffDir", &e_diffDir);
    selection->SetBranchAddress("e_diffKE", &e_diffKE);
    selection->SetBranchAddress("e_diffVtxLong", &e_diffVtxLong);
    selection->SetBranchAddress("e_diffVtxTrans", &e_diffVtxTrans);
    selection->SetBranchAddress("e_diff_nu_E", &e_diff_nu_E);
    selection->SetBranchAddress("nClusters", &nClusters);
    selection->SetBranchAddress("ring1PEs", &ring1PEs);
    selection->SetBranchAddress("ring2PEs", &ring2PEs);
    selection->SetBranchAddress("recoPID", &recoPID);
    selection->SetBranchAddress("recoPIDLikelihood", &recoPIDLikelihood);
    selection->SetBranchAddress("isHighE", &isHighE);
    selection->SetBranchAddress("nSubevents", &nSubevents);
    selection->SetBranchAddress("recoTime", recoTime);
    selection->SetBranchAddress("recoEnergy", recoEnergy);
    selection->SetBranchAddress("recoCaptures", &recoCaptures);

    TChain *f = new TChain("Final_Reconstruction");
    TChain *le = new TChain("Low_E");
    TChain *hee = new TChain("High_E_Electron");
    TChain *hem = new TChain("High_E_Muon");
    TChain *d = new TChain("Debug");
    TChain *h = new TChain("HitsTree");
    TChain *p = new TChain("PMTsTree");
    TChain *t=new TChain("tcardfile");
    TString flavs[4] = {"numu", "nue", "antinumu", "antinue"};
    TString hcs[4] = {"nu", "antinu"};
    for (int k = 0; k < 1; k++) {
        const char *hc = hcs[k].Data();
        for (int j = 0; j < 2; j++) {
            const char *s = flavs[j].Data();
            for (int i = 1000; i < 1100; i++) {
                //Missing vector files for antinu mode:
                if (k == 1 && j == 0 && (i == 1096 || i == 1037)) continue;
                if (k == 1 && j == 1 && (i == 1017 || i == 1053)) continue;
                if (k == 1 && j == 2 && i == 1078) continue;
                if (k == 1 && j == 3 && (i == 1018 || i == 1019 || i == 1053)) continue;
                char *file = Form("/data/hyperk/wchsandbox_reco/flav_%s/%s_%s_%i/%s_%s_%i_12in_reco_4.root", s, hc, s,
                                  i,
                                  hc, s, i);
                f->AddFile(file);
                le->AddFile(file);
                hee->AddFile(file);
                hem->AddFile(file);
                d->AddFile(file);
                h->AddFile(
                        Form("/data/hyperk/wchsandbox_reco/flav_%s/%s_%s_%i/%s_%s_%i_out_12in.root", s, hc, s, i, hc, s,
                             i));
                if (k == 0 && j == 0 && i == 1000)
                    p->AddFile(
                            Form("/data/hyperk/wchsandbox_reco/flav_%s/%s_%s_%i/%s_%s_%i_out_12in.root", s, hc, s, i,
                                 hc, s,
                                 i));
                t->AddFile(
                        Form("/data/hyperk/wchsandbox_reco/flav_%s/%s_%s_%i/%s_%s_%i_generatorcardfile.root", s, hc, s,
                             i, hc, s, i));
            }
        }
    }
    TFile gout("miscount_events_cardfile.root","RECREATE");
    TTree* t2 = (TTree*) t->CloneTree(0);
    TFile out("miscount_events_out.root","RECREATE");
    TTree* p2 = (TTree*) p->CloneTree();
    TTree* h2 = (TTree*) h->CloneTree(0);

    nentries = f->GetEntries();

    for(int i=0; i<100000; i++) {
        if(i%1000==0)cout << i << endl;
        h->GetEntry(i);
        selection->GetEntry(i);
        t->GetEntry(i);
        if((neutrino_id)==14
           &&interaction_mode>=11
           &&interaction_mode<=16
           &&n_chkv_part>1
           &&ring1PEs>10
           &&isHighE&&ring2PEs<0.09*ring1PEs
           &&recoPIDLikelihood<-200
           &&e_reco_dwall>100
           &&e_recoKE>100
           &&e_recoKE<2500
           &&e_reco_nu_E<2500){
            out.cd();
            h2->Fill();
            gout.cd();
            t2->Fill();
        }
    }
    gout.cd();
    t2->Write();
    gout.Close();
    out.cd();
    h2->Write();
    p2->Write();
    out.Close();
}