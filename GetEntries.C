// Macro to get the entries from the Valor input files
void GetEntries(){
  TFile g("selections_androot_tagged_FHC.root");
  cout << " FHC " << endl;

  int n_numu_mulike =numu_mulike_ntag->GetEntries()+numu_mulike_nontag->GetEntries();
  int n_nue_mulike =nue_mulike_ntag->GetEntries()+nue_mulike_nontag->GetEntries();
  int n_numu_elike =numu_elike_ntag->GetEntries()+numu_elike_nontag->GetEntries();
  int n_nue_elike =nue_elike_ntag->GetEntries()+nue_elike_nontag->GetEntries();
 
  int n_numubar_mulike =numubar_mulike_ntag->GetEntries()+numubar_mulike_nontag->GetEntries();
  int n_nuebar_mulike =nuebar_mulike_ntag->GetEntries()+nuebar_mulike_nontag->GetEntries();
  int n_numubar_elike =numubar_elike_ntag->GetEntries()+numubar_elike_nontag->GetEntries();
  int n_nuebar_elike =nuebar_elike_ntag->GetEntries()+nuebar_elike_nontag->GetEntries();
  
  cout << "True Nu \tSelection \tAll\tNtag \tNoNtag  " << endl;
  cout << "Numu \tMu-like \t" << n_numu_mulike << " \t" << numu_mulike_ntag->GetEntries() << " " << numu_mulike_nontag->GetEntries() << endl;
  cout << "Nue \tMu-like \t" << n_nue_mulike << " \t" << nue_mulike_ntag->GetEntries() << " " << nue_mulike_nontag->GetEntries() << endl;
  cout << "Numu \te-like \t" << n_numu_elike << " \t" << numu_elike_ntag->GetEntries() << " " << numu_elike_nontag->GetEntries() << endl;
  cout << "Nue \te-like \t" << n_nue_elike << " \t" << nue_elike_ntag->GetEntries() << " " << nue_elike_nontag->GetEntries() << endl;
  
  cout << "Numubar \tMu-like \t" << n_numubar_mulike << " \t" << numubar_mulike_ntag->GetEntries() << " " << numubar_mulike_nontag->GetEntries() << endl;
  cout << "Nuebar \tMu-like \t" << n_nuebar_mulike << " \t" << nuebar_mulike_ntag->GetEntries() << " " << nuebar_mulike_nontag->GetEntries() << endl;
  cout << "Numubar \te-like \t" << n_numubar_elike << " \t" << numubar_elike_ntag->GetEntries() << " " << numubar_elike_nontag->GetEntries() << endl;
  cout << "Nuebar \te-like \t" << n_nuebar_elike << " \t" << nuebar_elike_ntag->GetEntries() << " " << nuebar_elike_nontag->GetEntries() << endl;
  
}
