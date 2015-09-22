void recoAnalysis(int n = -1)
{
    TString flavs[4] = {"numu", "nue", "antinumu", "antinue"};
    TString horns[2] = {"nu", "antinu"};
    TChain *f=new TChain("Final_Reconstruction");
    TChain *le=new TChain("Low_E");
    TChain *hee=new TChain("High_E_Electron");
    TChain *hem=new TChain("High_E_Muon");
    TChain *d=new TChain("Debug");
    TString suffix = "";
    if (n!=-1) suffix = Form("_%i",n);
    for (int iflav = 0; iflav < 4; iflav++) {
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
    double q=0.68, r;
    new TCanvas();
    f->Draw("diffVtxAbs>>hMuVtx","isHighE[0]&&neutrinoPID==14&&recoPID[0]==13&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    hMuVtx->GetQuantiles(1,&r,&q);
    cout << "Muon vtx res:     " << r << endl;
    new TCanvas();
    f->Draw("diffVtxAbs>>heVtx","isHighE[0]&&neutrinoPID==12&&recoPID[0]==11&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    heVtx->GetQuantiles(1,&r,&q);
    cout << "Electron vtx res: " << r << endl;
    new TCanvas();
    f->Draw("diffDirAbs*180/TMath::Pi()>>hMuDir","isHighE[0]&&neutrinoPID==14&&recoPID[0]==13&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    hMuDir->GetQuantiles(1,&r,&q);
    cout << "Muon dir res:     " << r << endl;
    new TCanvas();
    f->Draw("diffDirAbs*180/TMath::Pi()>>heDir","isHighE[0]&&neutrinoPID==12&&recoPID[0]==11&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    heDir->GetQuantiles(1,&r,&q);
    cout << "Electron dir res: " << r << endl;
    new TCanvas();
    f->Draw("diffKE>>hMuEne","isHighE[0]&&neutrinoPID==14&&recoPID[0]==13&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    new TCanvas();
    f->Draw("diffKE>>heEne","isHighE[0]&&neutrinoPID==12&&recoPID[0]==11&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    new TCanvas();
    TH1D * hMuPID = new TH1D("hMuPID","PID",100,-2000,2000);
    double muTot = f->Draw("recoLnLHighEMuon[0]-recoLnLHighEElectron[0]>>hMuPID","isHighE[0]&&neutrinoPID==14&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    hMuPID->SetLineColor(kRed);
    double eTot = f->Draw("recoLnLHighEMuon[0]-recoLnLHighEElectron[0]>>hePID","isHighE[0]&&neutrinoPID==12&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100","same");
    hePID->SetLineColor(kBlue);
    double muMisID = f->Draw("","isHighE[0]&&neutrinoPID==14&&recoLnLHighEMuon[0]-recoLnLHighEElectron[0]<0&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    cout << "Muon mis-ID percentage:     " << muMisID/muTot << endl;
    double eMisID = f->Draw("","isHighE[0]&&neutrinoPID==12&&recoLnLHighEMuon[0]-recoLnLHighEElectron[0]>0&&abs(mode)==1&&trueKE>400&&trueKE<1000&&trueDWall>100");
    cout << "Electron mis-ID percentage: " << eMisID/eTot << endl;
}
