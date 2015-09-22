#include "TCut.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>

using namespace std;

void selectioncounts()
{
    TFile f("/data/prouse/HK/hk-TITUSanalysis/WChSandBox/analysis/selection_test_20150902_231510.root");
    TTree* t = (TTree*) f.Get("selection");
    TCut numu = "neutrino_id==14";
    TCut nue = "neutrino_id==12";
    TCut anumu = "neutrino_id==-14";
    TCut anue = "neutrino_id==-12";
    TCut nc = "abs(interaction_mode)>30";
    TCut cc = "abs(interaction_mode)<31";
    TCut ccqe = "abs(interaction_mode)==1";
    TCut ccnqe = cc&&"abs(interaction_mode)!=1";
    TCut ccmec = "abs(interaction_mode)==2";
    TCut cc1pi = "abs(interaction_mode)==11||abs(interaction_mode)==12||abs(interaction_mode)==13||abs(interaction_mode)==16";
    TCut ccdis = "abs(interaction_mode)==26";
    TCut ccother = !(cc1pi||ccqe||ccdis||ccmec||nc);
    TCut ccnonqe = cc1pi||ccmec||ccdis||ccother;
    TCut mudwall = "mu_reco_dwall>100";
    TCut edwall = "e_reco_dwall>100";
    TCut fv = "true_dwall>100";
    TCut vol = mudwall&&"reco_towall>200";
    TCut r1 = "isHighE&&ring2PEs<0.09*ring1PEs";
    TCut E = "reco_nu_E>0&&reco_nu_E<2000&&recoEnergy[0]>200&&recoEnergy[0]<1000";
    TCut Emu = "mu_recoKE>200&&mu_recoKE<2000";
    TCut Ee = "e_recoKE>100&&e_recoKE<2500";
    TCut Enumu = "mu_reco_nu_E<1250";
    TCut Enue = "e_reco_nu_E<2500";
    TCut mu = "recoPIDLikelihood>-200";
    TCut e = "recoPIDLikelihood<-200";
    TCut n0 = "nneutron==0";
    TCut n1 = "nneutron==1";
    TCut n2 = "nneutron==2";
    TCut n3 = "nneutron==3";
    TCut n4 = "nneutron==4";
    TCut n5 = "nneutron==5";
    TCut nontag = "recoCaptures==0";
    TCut ntag = "recoCaptures>0";
//    TCut col[24]={numu&&ccqe, numu&&ccmec, numu&&cc1pi, numu&&ccdis, numu&&ccother, numu&&nc,
//                  nue&&ccqe,  nue&&ccmec,  nue&&cc1pi,  nue&&ccdis,  nue&&ccother,  nue&&nc,
//                  anumu&&ccqe,anumu&&ccmec,anumu&&cc1pi,anumu&&ccdis,anumu&&ccother,anumu&&nc,
//                  anue&&ccqe, anue&&ccmec, anue&&cc1pi, anue&&ccdis, anue&&ccother, anue&&nc};
//    TCut row[18]={"1==1",
//                  vol,
//                  vol&&r1,
//                  vol&&r1&&E,
//                  vol&&r1&&E&&mu,
//                  vol&&r1&&E&&mu&&n0,
//                  vol&&r1&&E&&mu&&n1,
//                  vol&&r1&&E&&mu&&n2,
//                  vol&&r1&&E&&mu&&n3,
//                  vol&&r1&&E&&mu&&n4,
//                  vol&&r1&&E&&mu&&n5,
//                  vol&&r1&&E&&e,
//                  vol&&r1&&E&&e&&n0,
//                  vol&&r1&&E&&e&&n1,
//                  vol&&r1&&E&&e&&n2,
//                  vol&&r1&&E&&e&&n3,
//                  vol&&r1&&E&&e&&n4,
//                  vol&&r1&&E&&e&&n5};
    TCut colmu[18]={numu&&ccqe,numu&&ccnqe,anumu&&cc,nue&&cc,anue&&cc,numu&&nc,anumu&&nc,nue&&nc,anue&&nc,
                    anumu&&ccqe,anumu&&ccnqe,numu&&cc,nue&&cc,anue&&cc,numu&&nc,anumu&&nc,nue&&nc,anue&&nc};
    TCut rowmu[9]={"",fv,
                   mudwall,
                   mudwall&&r1,
                   mudwall&&r1&&mu,
                   mudwall&&r1&&mu&&Emu,
                   mudwall&&r1&&mu&&Enumu,
                   mudwall&&r1&&mu&&Enumu&&nontag,
                   mudwall&&r1&&mu&&Enumu&&ntag};
    TCut cole[18]={nue&&ccqe,nue&&ccnqe,anue&&cc,numu&&cc,anumu&&cc,numu&&nc,anumu&&nc,nue&&nc,anue&&nc,
                    anue&&ccqe,anue&&ccnqe,nue&&cc,numu&&cc,anumu&&cc,numu&&nc,anumu&&nc,nue&&nc,anue&&nc};
    TCut rowe[9]={"",fv,
                   edwall,
                   edwall&&r1,
                   edwall&&r1&&e,
                   edwall&&r1&&e&&Ee,
                   edwall&&r1&&e&&Enue,
                   edwall&&r1&&e&&Enue&&nontag,
                   edwall&&r1&&e&&Enue&&ntag};

    double fhc_numu = 4803420./t->Draw("",numu,"goff",400000,0);
    double fhc_anumu = 136851./t->Draw("",anumu,"goff",400000,0);
    double fhc_nue = 91393.6/t->Draw("",nue,"goff",400000,0);
    double fhc_anue = 9229.89/t->Draw("",anue,"goff",400000,0);
    Long64_t nentries = t->GetEntries() -400000;
    double rhc_numu = 609470./t->Draw("",numu,"goff",nentries,400000);
    double rhc_anumu = 1197020./t->Draw("",anumu,"goff",nentries,400000);
    double rhc_nue = 26893.5/t->Draw("",nue,"goff",nentries,400000);
    double rhc_anue = 23030.6/t->Draw("",anue,"goff",nentries,400000);
    double counts[18][18], out[18][10];

    for (int i = 0; i < 18; i++) {
        for(int j =0; j <18; j++) {
            TCut colCut = i<9 ? colmu[j] : cole[j];
            TCut rowCut = i<9 ? rowmu[i] : rowe[i -9];
            Long64_t count = j<9 ? 400000 : nentries;
            int start = j<9 ? 0 : 400000;
            counts[i][j] = t->Draw("", colCut && rowCut, "goff", count, start);
            //cout << j << " " << i << " " << counts[i][j] << endl;
            if(i <9)
            {
                if(j<2 || j==5) counts[i][j] *= fhc_numu;
                else if(j==11 || j==14) counts[i][j] *= rhc_numu;
                else if(j==2 || j==6) counts[i][j] *= fhc_anumu;
                else if(j==9 || j==10 || j==15) counts[i][j] *= rhc_anumu;
                else if(j==3 || j==7) counts[i][j] *= fhc_nue;
                else if(j==12 || j==16) counts[i][j] *= rhc_nue;
                else if(j==4 || j== 8) counts[i][j] *= fhc_anue;
                else if(j==13 || j==17) counts[i][j] *= rhc_anue;
            }
            else{
                if(j==3 || j ==5) counts[i][j] *= fhc_numu;
                else if(j==12 || j==14) counts[i][j] *= rhc_numu;
                else if(j==4 || j==6)counts[i][j] *= fhc_anumu;
                else if(j==13 || j==15) counts[i][j] *= rhc_anumu;
                else if(j<2 || j==7)counts[i][j] *= fhc_nue;
                else if(j==11 || j==16) counts[i][j] *= rhc_nue;
                else if(j==2 || j==8) counts[i][j] *= fhc_anue;
                else if(j==9 || j==10 || j==17) counts[i][j] *= rhc_anue;
            }
        }
        out[i][0]=counts[i][0];
        out[i][1]=counts[i][1];
        out[i][2]=counts[i][2];
        out[i][3]=counts[i][3]+counts[i][4];
        out[i][4]=counts[i][5]+counts[i][6]+counts[i][7]+counts[i][8];
        out[i][5]=counts[i][9];
        out[i][6]=counts[i][10];
        out[i][7]=counts[i][11];
        out[i][8]=counts[i][12]+counts[i][13];
        out[i][9]=counts[i][14]+counts[i][15]+counts[i][16]+counts[i][17];
    }
    for(int i=0; i<18;i++){
        for(int j=0;j<10;j++) {
            cout << out[i][j];
            if(j<9) cout << " & ";
        }
        cout << endl;
    }


    TString text[7] = {"All       ","FV        ","dwall>100 ","1R        ","e-like    ","Ee<100    ","Enu<2.5   "};
    TCut col[10] = {nue&&cc,
                    nue&&ccqe,
                    numu&&cc,
                    numu&&ccqe,
                    numu&&"interaction_mode==2",
                    numu&&"interaction_mode>=11&&interaction_mode<=16",
                    numu&&"interaction_mode==21",
                    numu&&"interaction_mode==26",
                    numu&&"interaction_mode==17||(interaction_mode>=18&&interaction_mode<=20)||interaction_mode==22||interaction_mode==23",
                    ""};
    int jmax=9;
    cout << setw(24) << "     nue CC     "
         << setw(17) << "    nue CCQE    "
         << setw(17) << "    numu CC     "
         << setw(17) << "   numu CCQE    "
         << setw(17) << "  numu CC Coh   "
         << setw(17) << "  numu CC 1pi   "
         << setw(17) << "  numu CC 2+pi  "
         << setw(17) << "   numu CCDIS   "
         << setw(17) << " numu CC other  "
         << endl;
    for(int i=0;i<7;i++){
        cout<<text[i];
        for(int j=0;j<jmax;j++){
            counts[i][j]=t->Draw("",rowe[i]&&col[j],"goff",400000,0)*(j<2 ? fhc_nue : fhc_numu);
            cout << setiosflags(ios::fixed) << setprecision(0) << setw(7) << counts[i][j] << " "
                 << setw(7) <<  setprecision(1) << 100.*counts[i][j]/counts[0][j] << "  ";
        }
        cout<< endl;
    }
}