void recoAnalysis2(int n1 = -1,int n2 = -1)
{
    TString flavs[4] = {"numu", "nue", "antinumu", "antinue"};
    TString horns[2] = {"nu", "antinu"};
    TChain *f=new TChain("Final_Reconstruction");
    TChain *le=new TChain("Low_E");
    TChain *hee=new TChain("High_E_Electron");
    TChain *hem=new TChain("High_E_Muon");
    TChain *d=new TChain("Debug");
    TString suffix = "";
    if (n1!=-1) suffix = Form("_%i",n1);
    for (int iflav = 0; iflav < 2; iflav++) {
        const char *flav = flavs[iflav].Data();
        for (int ihorn = 0; ihorn < 1; ihorn++) {
            const char *horn = horns[ihorn].Data();
            for (int i = 1000; i < 1100; i++) {
                char *file = Form("/data/hyperk/wchsandbox_reco/flav_%s/%s_%s_%i/%s_%s_%i_12in_reco%s.root",flav,horn,flav,i,horn,flav,i,suffix.Data());
                f->AddFile(file);
                le->AddFile(file);
                hee->AddFile(file);
                hem->AddFile(file);
                d->AddFile(file);
            }
        }
    }
    f->AddFriend(le);
    f->AddFriend(hee);
    f->AddFriend(hem);
    f->AddFriend(d);
    f->SetLineStyle(1);

    TChain *f2=new TChain("Final_Reconstruction");
    TChain *le2=new TChain("Low_E");
    TChain *hee2=new TChain("High_E_Electron");
    TChain *hem2=new TChain("High_E_Muon");
    TChain *d2=new TChain("Debug");
    TString suffix2 = "";
    if (n2!=-1) suffix2 = Form("_%i",n2);
    for (int iflav = 0; iflav < 2; iflav++) {
        const char *flav = flavs[iflav].Data();
        for (int ihorn = 0; ihorn < 1; ihorn++) {
            const char *horn = horns[ihorn].Data();
            for (int i = 1000; i < 1100; i++) {
                char *file = Form("/data/hyperk/wchsandbox_reco/flav_%s/%s_%s_%i/%s_%s_%i_12in_reco%s.root",flav,horn,flav,i,horn,flav,i,suffix2.Data());
                f2->AddFile(file);
                le2->AddFile(file);
                hee2->AddFile(file);
                hem2->AddFile(file);
                d2->AddFile(file);
            }
        }
    }
    f2->AddFriend(le2);
    f2->AddFriend(hee2);
    f2->AddFriend(hem2);
    f2->AddFriend(d2);
    f2->SetLineStyle(2);

    gStyle->SetOptStat(0);

    double q=0.68, r;
    new TCanvas();
    TH1D * hMuVtx = new TH1D("hMuVtx","Muon vertex resolution",100,0,200);
    f->Draw("diffVtxAbs>>hMuVtx","isHighE[0]&&neutrinoPID==14&&recoPID[0]==13&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    hMuVtx->GetQuantiles(1,&r,&q);
    cout << "Muon vtx res 1:     " << r << endl;
    f2->Draw("diffVtxAbs>>hMuVtx2","isHighE[0]&&neutrinoPID==14&&recoPID[0]==13&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100","same");
    hMuVtx2->GetQuantiles(1,&r,&q);
    cout << "Muon vtx res 2:     " << r << endl;
    hMuVtx->GetXaxis()->SetTitle("Difference between true and reconstructed vertex [cm]");

    new TCanvas();
    TH1D * heVtx = new TH1D("heVtx","Electron vertex resolution",100,0,400);
    f->Draw("diffVtxAbs>>heVtx","isHighE[0]&&neutrinoPID==12&&recoPID[0]==11&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    heVtx->GetQuantiles(1,&r,&q);
    cout << "Electron vtx res 1: " << r << endl;
    f2->Draw("diffVtxAbs>>heVtx2","isHighE[0]&&neutrinoPID==12&&recoPID[0]==11&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100","same");
    heVtx2->GetQuantiles(1,&r,&q);
    cout << "Electron vtx res 2: " << r << endl;
    heVtx->GetXaxis()->SetTitle("Difference between true and reconstructed vertex [cm]");

    new TCanvas();
    TH1D * hMuDir = new TH1D("hMuDir","Muon direction resolution",100,0,20);
    f->Draw("diffDirAbs*180/TMath::Pi()>>hMuDir","isHighE[0]&&neutrinoPID==14&&recoPID[0]==13&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    hMuDir->GetQuantiles(1,&r,&q);
    cout << "Muon dir res 1:     " << r << endl;
    f2->Draw("diffDirAbs*180/TMath::Pi()>>hMuDir2","isHighE[0]&&neutrinoPID==14&&recoPID[0]==13&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100","same");
    hMuDir2->GetQuantiles(1,&r,&q);
    cout << "Muon dir res 2:     " << r << endl;
    hMuDir->GetXaxis()->SetTitle("Difference between true and reconstructed direction [deg]");

    new TCanvas();
    TH1D * heDir = new TH1D("heDir","Electron direction resolution",100,0,40);
    f->Draw("diffDirAbs*180/TMath::Pi()>>heDir","isHighE[0]&&neutrinoPID==12&&recoPID[0]==11&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    heDir->GetQuantiles(1,&r,&q);
    cout << "Electron dir res 1: " << r << endl;
    f2->Draw("diffDirAbs*180/TMath::Pi()>>heDir2","isHighE[0]&&neutrinoPID==12&&recoPID[0]==11&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100","same");
    heDir2->GetQuantiles(1,&r,&q);
    cout << "Electron dir res 2: " << r << endl;
    heDir->GetXaxis()->SetTitle("Difference between true and reconstructed direction [deg]");

    new TCanvas();
    f->Draw("diffKE>>hMuEne","isHighE[0]&&neutrinoPID==14&&recoPID[0]==13&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    f2->Draw("diffKE>>hMuEne2","isHighE[0]&&neutrinoPID==14&&recoPID[0]==13&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100","same");
    hMuEne->SetTitle("Muon energy resolution");
    hMuEne->GetXaxis()->SetTitle("E_{rec}-E_{true} [MeV]");

    new TCanvas();
    f->Draw("diffKE>>heEne","isHighE[0]&&neutrinoPID==12&&recoPID[0]==11&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    f2->Draw("diffKE>>heEne2","isHighE[0]&&neutrinoPID==12&&recoPID[0]==11&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100","same");
    heEne->SetTitle("Electron energy resolution");
    heEne->GetXaxis()->SetTitle("E_{rec}-E_{true} [MeV]");

    new TCanvas();
    f->Draw("diffKE/trueKE>>hMuEne3","isHighE[0]&&neutrinoPID==14&&recoPID[0]==13&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    f2->Draw("diffKE/trueKE>>hMuEne4","isHighE[0]&&neutrinoPID==14&&recoPID[0]==13&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100","same");
    hMuEne3->SetTitle("Muon energy resolution");
    hMuEne3->GetXaxis()->SetTitle("(E_{rec}-E_{true})/E_{true} [MeV]");

    new TCanvas();
    f->Draw("diffKE/trueKE>>heEne3","isHighE[0]&&neutrinoPID==12&&recoPID[0]==11&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    f2->Draw("diffKE/trueKE>>heEne4","isHighE[0]&&neutrinoPID==12&&recoPID[0]==11&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100","same");
    heEne3->SetTitle("Electron energy resolution");
    heEne3->GetXaxis()->SetTitle("(E_{rec}-E_{true})/E_{true} [MeV]");

    new TCanvas();
    TH1D * hMuPID = new TH1D("hMuPID","PID",100,-3000,3000);
    double muTot = f->Draw("recoLnLHighEMuon[0]-recoLnLHighEElectron[0]>>hMuPID","isHighE[0]&&neutrinoPID==14&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    hMuPID->SetLineColor(kRed);
    double eTot = f->Draw("recoLnLHighEMuon[0]-recoLnLHighEElectron[0]>>hePID","isHighE[0]&&neutrinoPID==12&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100","same");
    hePID->SetLineColor(kBlue);
    double muMisID = f->Draw("","isHighE[0]&&neutrinoPID==14&&recoLnLHighEMuon[0]-recoLnLHighEElectron[0]<0&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    cout << "Muon mis-ID percentage 1:     " << muMisID/muTot << endl;
    double eMisID = f->Draw("","isHighE[0]&&neutrinoPID==12&&recoLnLHighEMuon[0]-recoLnLHighEElectron[0]>0&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    cout << "Electron mis-ID percentage 1: " << eMisID/eTot << endl;
    TH1D * hMuPID2 = new TH1D("hMuPID2","PID",100,-3000,3000);
    double muTot2 = f2->Draw("recoLnLHighEMuon[0]-recoLnLHighEElectron[0]>>hMuPID2","isHighE[0]&&neutrinoPID==14&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100","same");
    hMuPID2->SetLineColor(kRed);
    hMuPID2->SetLineStyle(2);
    double eTot2 = f2->Draw("recoLnLHighEMuon[0]-recoLnLHighEElectron[0]>>hePID2","isHighE[0]&&neutrinoPID==12&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100","same");
    hePID2->SetLineColor(kBlue);
    hePID2->SetLineStyle(2);
    double muMisID2 = f2->Draw("","isHighE[0]&&neutrinoPID==14&&recoLnLHighEMuon[0]-recoLnLHighEElectron[0]<0&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    cout << "Muon mis-ID percentage 2:     " << muMisID2/muTot2 << endl;
    double eMisID2 = f2->Draw("","isHighE[0]&&neutrinoPID==12&&recoLnLHighEMuon[0]-recoLnLHighEElectron[0]>0&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    cout << "Electron mis-ID percentage 2: " << eMisID2/eTot2 << endl;
}
