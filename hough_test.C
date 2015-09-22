TGraph2D * hough_test(int i=0){
    TFile * tfile = new TFile("miscount_events_out.root");
    TTree* p = (TTree*) tfile->Get("PMTsTree");
    int nPMTs=p->GetEntries();
    double *pmt_x=new double[nPMTs];
    double *pmt_y=new double[nPMTs];
    double *pmt_z=new double[nPMTs];
    double *pmt_theta=new double[nPMTs];
    double *pmt_phi=new double[nPMTs];
    double *pmt_hits=new double[nPMTs];
    double pmtx, pmty, pmtz;
    p->SetBranchAddress("pmt_pos_x",&pmtx);
    p->SetBranchAddress("pmt_pos_y",&pmty);
    p->SetBranchAddress("pmt_pos_z",&pmtz);
    for(int iPMT=0; iPMT<nPMTs; iPMT++){
        p->GetEntry(iPMT);
        pmt_x[iPMT]=pmtx;
        pmt_y[iPMT]=pmty;
        pmt_z[iPMT]=pmtz;
        pmt_hits[iPMT]=0;
    }
    TTree* t = (TTree*) tfile->Get("HitsTree");
    int npart, nhits;
    double * part_pxStart = new double[100000];
    double * part_pyStart = new double[100000];
    double * part_pzStart = new double[100000];
    double * part_tStart = new double[100000];
    double * part_KEstart = new double[100000];
    int * part_pid = new int[100000];
    double * hit_x = new double[100000];
    double * hit_y = new double[100000];
    double * hit_z = new double[100000];
    double * hit_time = new double[100000];
    int * hit_PMTid = new int[100000];
    t->SetBranchAddress("npart",&npart);
    t->SetBranchAddress("part_pxStart",part_pxStart);
    t->SetBranchAddress("part_pyStart",part_pyStart);
    t->SetBranchAddress("part_pzStart",part_pzStart);
    t->SetBranchAddress("part_tStart",part_tStart);
    t->SetBranchAddress("part_KEstart",part_KEstart);
    t->SetBranchAddress("part_pid",part_pid);
    t->SetBranchAddress("nhits",&nhits);
    t->SetBranchAddress("hit_x",hit_x);
    t->SetBranchAddress("hit_y",hit_y);
    t->SetBranchAddress("hit_z",hit_z);
    t->SetBranchAddress("hit_time",hit_time);
    t->SetBranchAddress("hit_PMTid",hit_PMTid);
    TFile *file = new TFile("miscount_events_reco3.root");
    TTree* d = (TTree*) file->Get("Debug");
    TTree* l = (TTree*) file->Get("Low_E");
    TTree* f = (TTree*) file->Get("Final_Reconstruction");
    d->AddFriend(l);
    d->AddFriend(f);
    const int nBins = 2500;
    double* hough = new double[nBins];
    double * dirX = new double[20];
    double * dirY = new double[20];
    double * dirZ = new double[20];
    double * vtxX = new double[20];
    double * vtxY = new double[20];
    double * vtxZ = new double[20];
    int mode, evt;
    int * rings = new int[20];
    d->SetBranchAddress("houghSpace",hough);
    d->SetBranchAddress("mode",&mode);
    d->SetBranchAddress("evt",&evt);
    f->SetBranchAddress("recoNRings",rings);
    d->SetBranchAddress("pointVtxX",vtxX);
    d->SetBranchAddress("pointVtxY",vtxY);
    d->SetBranchAddress("pointVtxZ",vtxZ);
    double * theta = new double[nBins];
    double * phi = new double[nBins];
    double * theta2 = new double[nBins];
    double * phi2 = new double[nBins];
    //Distribute bins evenly around unit sphere using spiral method
    //Calculate shift to correct for end bins covering larger area
    double shift = 0.5;
    //Use Newton's method with initial guess 0.5 to solve equation:
    //  a*sin(a) - pi/N = 0
    //where
    //  a=sqrt(pi/(N-x))
    //  x is correct shift
    //only need one iteration for good enough estimate
    double a= TMath::Sqrt(TMath::Pi()/(nBins-shift));
    shift -= (a* TMath::Sin(a)- TMath::Pi()/nBins)*2.0*(nBins-shift)/(a*(TMath::Sin(a)+a* TMath::Cos(a)));
    double dz = 2.0/(nBins-shift); //Change in z=cos(theta) for consecutive bins
    double s = TMath::Sqrt(4* TMath::Pi()/(nBins-shift)); //Approx distance between bins
    double sOverdz = s/dz;
    //Calculate theta and phi position of each bin according to position on spiral
    d->GetEntry(i);
    f->GetEntry(i);
    t->GetEntry(i);
    cout << "evt: " << setw(3) << evt << "  mode: " << setw(2) << mode << "  Found rings: " << rings[0] << endl;
    double peak=0, peaktheta=0, peakphi=0;
    for(int iBin=0; iBin<nBins; iBin++){
        theta[iBin] = TMath::ACos(1-dz*(iBin+0.5*(1-shift)));
        phi[iBin] = theta[iBin]*sOverdz;
        theta2[iBin] = 90.-TMath::ACos(TMath::Sin(theta[iBin])*TMath::Cos(phi[iBin]))*180./TMath::Pi();
        phi2[iBin] = TMath::ATan2(TMath::Sin(theta[iBin])*TMath::Sin(phi[iBin]),TMath::Cos(theta[iBin]))*180./TMath::Pi();
        if(hough[iBin]>peak){
            peak=hough[iBin];
            peaktheta=theta2[iBin];
            peakphi=phi2[iBin];
        }
    }
    const Int_t NRGBs = 5;
    const Int_t NCont = 99;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    TGraph2D * g = new TGraph2D("Hough space","hough space",nBins,phi2,theta2, hough);
    TCanvas * c = new TCanvas();
    c->Divide(2,2);
    c->cd(1);
    g->Draw("aitoff");
    c->cd(2);
    //new TCanvas();
    g->Draw("col");
    TMarker * mPeak = new TMarker(peakphi,peaktheta,30);
    mPeak->Draw();
    for(int ipart=0; ipart<npart; ipart++){
        if(part_KEstart[ipart]<10 || part_tStart[ipart]>80) continue;
        if(abs(part_pid[ipart])!=11 && abs(part_pid[ipart])!=13 && part_pid[ipart]!=22 && abs(part_pid[ipart])!=211) continue;
        if(abs(part_pid[ipart])==13 && part_KEstart[ipart]<60) continue;
        if(abs(part_pid[ipart])==211 && part_KEstart[ipart]<70) continue;
        double part_theta = 90.-TMath::ACos(part_pxStart[ipart])*180./TMath::Pi();
        double part_phi = TMath::ATan2(part_pyStart[ipart],part_pzStart[ipart])*180./TMath::Pi();
        TMarker * m = new TMarker(part_phi,part_theta,2);
        if(part_pid[ipart]==13) m->SetMarkerStyle(26);
        else if(part_pid[ipart]==211) m->SetMarkerStyle(24);
        else if(part_pid[ipart]==11) m->SetMarkerStyle(5);
        m->Draw();
    }
    for(int iPMT=0; iPMT<nPMTs; iPMT++){
        pmt_hits[iPMT]=0;
        double x = pmt_x[iPMT]-vtxX[0];
        double y = pmt_y[iPMT]-vtxY[0];
        double z = pmt_z[iPMT]-vtxZ[0];
        pmt_theta[iPMT] = 90.-TMath::ACos((x/TMath::Sqrt(x*x+y*y+z*z)))*180./TMath::Pi();
        pmt_phi[iPMT] = TMath::ATan2(y,z)*180./TMath::Pi();
    }
    for(int iHit=0; iHit<nhits; iHit++){
        if(hit_time[iHit]<120)
            pmt_hits[hit_PMTid[iHit]]++;
    }
    TGraph2D * g2 = new TGraph2D("PMT hits","PMT hits",nPMTs, pmt_phi, pmt_theta, pmt_hits);
    c->cd(3);
    g2->Draw("aitoff");
    c->cd(4);
    //new TCanvas();
    g2->Draw("cont4");
    return g;
}
