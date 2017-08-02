//INCLUDES
//#include <vector>
//#include "TH2F.h"

//DEFINITIONS
// PIDA vs energy for each plane
TH1F *h_MuonPida[3];
TH1F *h_PionPida[3];
TH1F *h_ProtonPida[3];

TH1F *h_MuonPida_cut[3];
TH1F *h_PionPida_cut[3];
TH1F *h_ProtonPida_cut[3];

// PID efficiency vs energy for each plane
TGraphAsymmErrors *h_EffEnergy[3][3];
TH1F *h_EffEnergy_pass[3][3]; // Count correct IDs
TH1F *h_EffEnergy_tot[3][3]; // Count all IDs, divide later

TGraphAsymmErrors *h_EffEnergy_cut[3][3];
TH1F *h_EffEnergy_pass_cut[3][3]; // Count correct IDs
TH1F *h_EffEnergy_tot_cut[3][3]; // Count all IDs, divide later

TH2F *h_PidEnergy[3][3];
TH2F *h_PidEnergy_cut[3][3];

double efficiency[3][3];
double total[3][3];
double efficiency_cut[3][3];
double total_cut[3][3];

// Vectors to store all histograms and graphs
std::vector<TH2*> v_AllTH2;
std::vector<TH1*> v_AllTH1;
std::vector<TGraphAsymmErrors*> v_AllTGraph;

int nContainedCut = 0;
int nMatchCut = 0;
int nEndPtCut = 0;
int nDaughtersCut = 0;

void calc_PIDAeff_Begin(TTree*){

  for (int i = 0; i < 3; i++){

    TString h_MuonPida_name = Form("h_MuonPida_plane%i",i); 
    h_MuonPida[i] = new TH1F(h_MuonPida_name,"",100,0,30);
    h_MuonPida[i]->GetXaxis()->SetTitle("PIDA Value");
    h_MuonPida[i]->GetYaxis()->SetTitle("Events (/0.3)");
    h_MuonPida[i]->GetYaxis()->SetTitleOffset(1.05);
    h_MuonPida[i]->SetLineColor(4);
    h_MuonPida[i]->SetLineWidth(3);
    v_AllTH1.push_back(h_MuonPida[i]);

    TString h_PionPida_name = Form("h_PionPida_plane%i",i); 
    h_PionPida[i] = new TH1F(h_PionPida_name,"",100,0,30);
    h_PionPida[i]->GetXaxis()->SetTitle("PIDA Value");
    h_PionPida[i]->GetYaxis()->SetTitle("Events (/0.3)");
    h_PionPida[i]->GetYaxis()->SetTitleOffset(1.05);
    h_PionPida[i]->SetLineColor(2);
    h_PionPida[i]->SetLineWidth(3);
    v_AllTH1.push_back(h_PionPida[i]);

    TString h_ProtonPida_name = Form("h_ProtonPida_plane%i",i); 
    h_ProtonPida[i] = new TH1F(h_ProtonPida_name,"",100,0,30);
    h_ProtonPida[i]->GetXaxis()->SetTitle("PIDA Value");
    h_ProtonPida[i]->GetYaxis()->SetTitle("Events (/0.3)");
    h_ProtonPida[i]->GetYaxis()->SetTitleOffset(1.05);
    h_ProtonPida[i]->SetLineColor(8);
    h_ProtonPida[i]->SetLineWidth(3);
    v_AllTH1.push_back(h_ProtonPida[i]);

    TString h_MuonPida_cut_name = Form("h_MuonPida_cut_plane%i",i); 
    h_MuonPida_cut[i] = new TH1F(h_MuonPida_cut_name,"",100,0,30);
    h_MuonPida_cut[i]->GetXaxis()->SetTitle("PIDA Value");
    h_MuonPida_cut[i]->GetYaxis()->SetTitle("Events (/0.3)");
    h_MuonPida_cut[i]->GetYaxis()->SetTitleOffset(1.05);
    h_MuonPida_cut[i]->SetLineColor(4);
    h_MuonPida_cut[i]->SetLineWidth(3);
    v_AllTH1.push_back(h_MuonPida_cut[i]);

    TString h_PionPida_cut_name = Form("h_PionPida_cut_plane%i",i); 
    h_PionPida_cut[i] = new TH1F(h_PionPida_cut_name,"",100,0,30);
    h_PionPida_cut[i]->GetXaxis()->SetTitle("PIDA Value");
    h_PionPida_cut[i]->GetYaxis()->SetTitle("Events (/0.3)");
    h_PionPida_cut[i]->GetYaxis()->SetTitleOffset(1.05);
    h_PionPida_cut[i]->SetLineColor(2);
    h_PionPida_cut[i]->SetLineWidth(3);
    v_AllTH1.push_back(h_PionPida_cut[i]);

    TString h_ProtonPida_cut_name = Form("h_ProtonPida_cut_plane%i",i); 
    h_ProtonPida_cut[i] = new TH1F(h_ProtonPida_cut_name,"",100,0,30);
    h_ProtonPida_cut[i]->GetXaxis()->SetTitle("PIDA Value");
    h_ProtonPida_cut[i]->GetYaxis()->SetTitle("Events (/0.3)");
    h_ProtonPida_cut[i]->GetYaxis()->SetTitleOffset(1.05);
    h_ProtonPida_cut[i]->SetLineColor(8);
    h_ProtonPida_cut[i]->SetLineWidth(3);
    v_AllTH1.push_back(h_ProtonPida_cut[i]);

    for (int j=0; j<3; j++){ //0 muon, 1 pion, 2 proton

      double max_energy = 2.3;
      if(j==0) max_energy = 1.3;

      TString h_EffEnergy_name = Form("h_EffEnergy_plane%i_%i",i,j);
      h_EffEnergy[i][j] = new TGraphAsymmErrors();
      h_EffEnergy[i][j]->SetName(h_EffEnergy_name);
      h_EffEnergy[i][j]->SetTitle(";Energy (GeV);PIDA Efficiency (%)");
      h_EffEnergy[i][j]->SetMarkerStyle(8);
      h_EffEnergy[i][j]->SetMarkerColor(kBlue);
      h_EffEnergy[i][j]->SetLineColor(kBlue);
      h_EffEnergy[i][j]->SetLineWidth(3);
      v_AllTGraph.push_back(h_EffEnergy[i][j]);

      TString h_EffEnergy_pass_name = Form("h_EffEnergy_pass_plane%i_%i",i,j);
      h_EffEnergy_pass[i][j] = new TH1F(h_EffEnergy_pass_name,"",10,0,max_energy);
      TString h_EffEnergy_tot_name = Form("h_EffEnergy_tot_plane%i_%i",i,j);
      h_EffEnergy_tot[i][j] = new TH1F(h_EffEnergy_tot_name,"",10,0,max_energy);

      TString h_EffEnergy_cut_name = Form("h_EffEnergy_cut_plane%i_%i",i,j);
      h_EffEnergy_cut[i][j] = new TGraphAsymmErrors();
      h_EffEnergy_cut[i][j]->SetName(h_EffEnergy_cut_name);
      h_EffEnergy_cut[i][j]->SetTitle(";Energy (GeV);PIDA Efficiency (%)");
      h_EffEnergy_cut[i][j]->SetMarkerStyle(8);
      h_EffEnergy_cut[i][j]->SetMarkerColor(kBlue);
      h_EffEnergy_cut[i][j]->SetLineColor(kBlue);
      h_EffEnergy_cut[i][j]->SetLineWidth(3);
      v_AllTGraph.push_back(h_EffEnergy_cut[i][j]);

      TString h_EffEnergy_pass_cut_name = Form("h_EffEnergy_pass_cut_plane%i_%i",i,j);
      h_EffEnergy_pass_cut[i][j] = new TH1F(h_EffEnergy_pass_cut_name,"",10,0,max_energy);
      TString h_EffEnergy_tot_cut_name = Form("h_EffEnergy_tot_cut_plane%i_%i",i,j);
      h_EffEnergy_tot_cut[i][j] = new TH1F(h_EffEnergy_tot_cut_name,"",10,0,max_energy);

      TString h_PidEnergy_name = Form("h_PidEnergy_plane%i_%i",i,j); 
      h_PidEnergy[i][j] = new TH2F(h_PidEnergy_name,"",10,0,max_energy,10,0,5);
      h_PidEnergy[i][j]->GetXaxis()->SetTitle("Energy (GeV)");
      h_PidEnergy[i][j]->GetYaxis()->SetTitle("#mu      #pi      K       p     ");
      h_PidEnergy[i][j]->GetYaxis()->SetTitleOffset(0.3);
      h_PidEnergy[i][j]->GetYaxis()->SetLabelSize(0.0);
      v_AllTH2.push_back(h_PidEnergy[i][j]);

      TString h_PidEnergy_cut_name = Form("h_PidEnergy_cut_plane%i_%i",i,j); 
      h_PidEnergy_cut[i][j] = new TH2F(h_PidEnergy_cut_name,"",10,0,max_energy,10,0,5);
      h_PidEnergy_cut[i][j]->GetXaxis()->SetTitle("Energy (GeV)");
      h_PidEnergy_cut[i][j]->GetYaxis()->SetTitle("#mu      #pi      K       p     ");
      h_PidEnergy_cut[i][j]->GetYaxis()->SetTitleOffset(0.3);
      h_PidEnergy_cut[i][j]->GetYaxis()->SetLabelSize(0.0);
      v_AllTH2.push_back(h_PidEnergy_cut[i][j]);

      efficiency[i][j]=0.0;
      total[i][j]=0.0;
      efficiency_cut[i][j]=0.0;
      total_cut[i][j]=0.0;

    }

  }

}

int calc_PIDAeff() {

/*
  std::cout<<endl<<"Event "<<event<<endl<<endl;

  std::cout<<"Primary particle = "<<cry_primaries_pdg[0]<<", Energy = "<<cry_Eng[0]<<" GeV"<<endl;
  std::cout<<"Number of daughters = "<<cry_ND[0]<<endl;
  std::cout<<"Number of reconstructed tracks = "<<ntracks_pmalgtrackmaker<<endl;

  for(int iTrk = 0; iTrk < ntracks_pmalgtrackmaker; iTrk++){
    std::cout<<"Track["<<iTrk<<"]:"<<endl;
    for(int iPln = 0; iPln < 3; iPln++){
      std::cout<<"Plane["<<iPln<<"]: Truth PDG = "<<trkpdgtruth_pmalgtrackmaker[iTrk][iPln]<<", PID PDG = "<<trkpidpdg_pmalgtrackmaker[iTrk][iPln]<<endl;
    }
    std::cout<<endl;
  }
*/

  double mu_max[3] = {6.3,8.3,8.7};
  double pi_max[3] = {10.5,15.0,14.0};

  //Here there is only one primary particle, in this case the primary particle is always top of the list
  for (int iPln=0; iPln<3; iPln++){
    double PIDA = trkpidpida_pmalgtrackmaker[0][iPln];
    if (PIDA>0){
      int pidcode = 0;
      if (PIDA < mu_max[iPln]) pidcode = 1; //Muon
      else if (PIDA >= mu_max[iPln] && PIDA <= pi_max[iPln]) pidcode = 2; //Pion
      else if (PIDA > pi_max[iPln]) pidcode = 4; //Proton
      if (pdg[0]==13) {
        h_MuonPida[iPln]->Fill(PIDA);
        h_PidEnergy[iPln][0]->Fill(cry_P[0],pidcode);
        total[iPln][0] += 1;
        if (pidcode == 1){ 
          h_EffEnergy_pass[iPln][0]->Fill(cry_P[0]);
          h_EffEnergy_tot[iPln][0]->Fill(cry_P[0]);
          efficiency[iPln][0] += 1;
        }
        else h_EffEnergy_tot[iPln][0]->Fill(cry_P[0]);
      }
      else if (abs(pdg[0])==211) {
        h_PionPida[iPln]->Fill(PIDA);
        h_PidEnergy[iPln][1]->Fill(cry_P[0],pidcode);
        total[iPln][1] += 1;
        if (pidcode == 2){
          h_EffEnergy_pass[iPln][1]->Fill(cry_P[0]);
          h_EffEnergy_tot[iPln][1]->Fill(cry_P[0]);
          efficiency[iPln][1] += 1;
        }
        else h_EffEnergy_tot[iPln][1]->Fill(cry_P[0]);
      }
      else if (pdg[0]==2212) { 
        h_ProtonPida[iPln]->Fill(PIDA);
        h_PidEnergy[iPln][2]->Fill(cry_P[0],pidcode);
        total[iPln][2] += 1;
        if (pidcode == 4){ 
          h_EffEnergy_pass[iPln][2]->Fill(cry_P[0]);
          h_EffEnergy_tot[iPln][2]->Fill(cry_P[0]);
          efficiency[iPln][2] += 1;
        }
        else h_EffEnergy_tot[iPln][2]->Fill(cry_P[0]);
      }
    }
  }

  //Check contained - for all particles
  //Contained means fully contained in TPC1 ONLY
  bool isContained = 0;
  //First particle in g4 array is always the particle that we want
  for (int i = 0; i < 1; i++){
    if (
        StartPointx[i] > 0 && StartPointx[i] < 200 
        &&  
        StartPointy[i] > -200 && StartPointy[i] < 200 
        &&  
        StartPointz[i] > 0 && StartPointz[i] < 500 
        &&  
        EndPointx[i] > 0 && EndPointx[i] < 200 
        &&  
        EndPointy[i] > -200 && EndPointy[i] < 200 
        &&  
        EndPointz[i] > 0 && EndPointz[i] < 500 
        ){  
      isContained = true;
    }
    else nContainedCut++;
  }

  //Now count how many reco tracks are associated to the particle we want - for all particles
  //First particle is always the particle we want
  int trueID = TrackId[0];
  int recoIndex = -99999;
  int NRecoTracksTruth = 0;
  //Loop over reco
  for (int i = 0; i < ntracks_pmalgtrackmaker; i++){
    //If the matched ID matches our particle then count it
    //Stuff is stored by plane so demand all 3 planes match
    if (trkidtruth_pmalgtrackmaker[i][0] == trueID && trkidtruth_pmalgtrackmaker[i][1] == trueID && trkidtruth_pmalgtrackmaker[i][2] == trueID){
      NRecoTracksTruth++;
      recoIndex = i; //yes it should equal i
    }   
  }
  if (NRecoTracksTruth!=1) nMatchCut++;

  //Apply selection
  if (!(isContained && NRecoTracksTruth==1)) return 0;
  //Demand track orientated the same way as the true particle
  TVector3 TrueStart(StartPointx[0],StartPointy[0],StartPointz[0]);
  TVector3 TrueEnd(EndPointx[0],EndPointy[0],EndPointz[0]);
  TVector3 RecoStart(trkstartx_pmalgtrackmaker[recoIndex],trkstarty_pmalgtrackmaker[recoIndex],trkstartz_pmalgtrackmaker[recoIndex]);
  TVector3 RecoEnd(trkendx_pmalgtrackmaker[recoIndex],trkendy_pmalgtrackmaker[recoIndex],trkendz_pmalgtrackmaker[recoIndex]);

  if (!((RecoStart-TrueStart).Mag() < (RecoEnd-TrueStart).Mag() && (RecoEnd-TrueEnd).Mag() < (RecoStart-TrueEnd).Mag())) return 0;

  if (abs(pdg[0])==211){
    if ((RecoEnd-TrueEnd).Mag() > 1.){nEndPtCut++; return 0;}
    if (NumberDaughters[0] > 3){nDaughtersCut++; return 0;}
  }
  else if (pdg[0]==2212){
    if ((RecoEnd-TrueEnd).Mag() > 0.3){nEndPtCut++; return 0;}
    if (NumberDaughters[0] > 3){nDaughtersCut++; return 0;}
  }
  else if (pdg[0]==13){
    if ((RecoEnd-TrueEnd).Mag() > 1.){nEndPtCut++; return 0;}

  }
  else {
    if ((RecoEnd-TrueEnd).Mag() > 1.){nEndPtCut++; return 0;}
  } 

  //Here there is only one primary particle, in this case the primary particle is always top of the list
  for (int iPln=0; iPln<3; iPln++){
    double PIDA = trkpidpida_pmalgtrackmaker[0][iPln];
    if (PIDA>0){
      int pidcode = 0;
      if (PIDA < mu_max[iPln]) pidcode = 1; //Muon
      else if (PIDA >= mu_max[iPln] && PIDA <= pi_max[iPln]) pidcode = 2; //Pion
      else if (PIDA > pi_max[iPln]) pidcode = 4;
      if (pdg[0]==13) {
        h_MuonPida_cut[iPln]->Fill(PIDA);
        h_PidEnergy_cut[iPln][0]->Fill(cry_P[0],pidcode);
        total_cut[iPln][0] += 1;
        if (pidcode == 1){ 
          h_EffEnergy_pass_cut[iPln][0]->Fill(cry_P[0]);
          h_EffEnergy_tot_cut[iPln][0]->Fill(cry_P[0]);
          efficiency_cut[iPln][0] += 1;
        }
        else h_EffEnergy_tot_cut[iPln][0]->Fill(cry_P[0]);
      }
      else if (abs(pdg[0])==211) {
        h_PionPida_cut[iPln]->Fill(PIDA);
        h_PidEnergy_cut[iPln][1]->Fill(cry_P[0],pidcode);
        total_cut[iPln][1] += 1;
        if (pidcode == 2){
          h_EffEnergy_pass_cut[iPln][1]->Fill(cry_P[0]);
          h_EffEnergy_tot_cut[iPln][1]->Fill(cry_P[0]);
          efficiency_cut[iPln][1] += 1;
        }
        else h_EffEnergy_tot_cut[iPln][1]->Fill(cry_P[0]);
      }
      else if (pdg[0]==2212) { 
        h_ProtonPida_cut[iPln]->Fill(PIDA);
        h_PidEnergy_cut[iPln][2]->Fill(cry_P[0],pidcode);
        total_cut[iPln][2] += 1;
        if (pidcode == 4){ 
          h_EffEnergy_pass_cut[iPln][2]->Fill(cry_P[0]);
          h_EffEnergy_tot_cut[iPln][2]->Fill(cry_P[0]);
          efficiency_cut[iPln][2] += 1;
        }
        else h_EffEnergy_tot_cut[iPln][2]->Fill(cry_P[0]);
      }
    }
  }

  return 0;
}

void calc_PIDAeff_Terminate(){

  TCanvas *c1 = new TCanvas("c1","",800,600);
  TLegend *leg = new TLegend( 0.586, 0.511, 0.812, 0.802 );
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  c1->SetBottomMargin(0.15);

  for (int i = 0; i<3; i++){
    c1->Clear();
    leg->Clear();
    Double_t muon_scale = 1./h_MuonPida[i]->Integral();
    h_MuonPida[i]->Scale(muon_scale);
    Double_t pion_scale = 1./h_PionPida[i]->Integral();
    h_PionPida[i]->Scale(pion_scale);
    Double_t proton_scale = 1./h_ProtonPida[i]->Integral();
    h_ProtonPida[i]->Scale(proton_scale);
    h_MuonPida[i]->Draw("HIST");
    h_PionPida[i]->Draw("HIST SAME");
    h_ProtonPida[i]->Draw("HIST SAME");
    leg->AddEntry(h_MuonPida[i],"#mu","l");
    leg->AddEntry(h_PionPida[i],"#pi","l");
    leg->AddEntry(h_ProtonPida[i],"p","l");
    leg->Draw();
    TString canv_name = Form("PIDA_plane%i.root",i);
    c1->SaveAs(canv_name);
  }

  for (int i = 0; i<3; i++){
    c1->Clear();
    leg->Clear();
    Double_t muon_scale = 1./h_MuonPida_cut[i]->Integral();
    h_MuonPida_cut[i]->Scale(muon_scale);
    Double_t pion_scale = 1./h_PionPida_cut[i]->Integral();
    h_PionPida_cut[i]->Scale(pion_scale);
    Double_t proton_scale = 1./h_ProtonPida_cut[i]->Integral();
    h_ProtonPida_cut[i]->Scale(proton_scale);
    h_MuonPida_cut[i]->Draw("HIST");
    h_PionPida_cut[i]->Draw("HIST SAME");
    h_ProtonPida_cut[i]->Draw("HIST SAME");
    leg->AddEntry(h_MuonPida_cut[i],"#mu","l");
    leg->AddEntry(h_PionPida_cut[i],"#pi","l");
    leg->AddEntry(h_ProtonPida_cut[i],"p","l");
    leg->Draw();
    TString canv_name = Form("PIDA_plane%i_cut.root",i);
    c1->SaveAs(canv_name);
  }

  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      h_EffEnergy[i][j]->BayesDivide(h_EffEnergy_pass[i][j],h_EffEnergy_tot[i][j]);
      h_EffEnergy_cut[i][j]->BayesDivide(h_EffEnergy_pass_cut[i][j],h_EffEnergy_tot_cut[i][j]);
      efficiency[i][j] /= total[i][j];
      efficiency_cut[i][j] /= total_cut[i][j];
      std::cout<<"["<<i<<"]["<<j<<"]: Before cuts = "<<efficiency[i][j]<<" After cuts = "<<efficiency_cut[i][j]<<std::endl;
    }
  }


  //Write all histograms and graphs to file
  for (unsigned int i = 0; i < v_AllTH2.size(); i++){
    v_AllTH2[i]->Write();
  }
  for (unsigned int i = 0; i < v_AllTH1.size(); i++){
    v_AllTH1[i]->Write();
  }
  for (unsigned int i = 0; i < v_AllTGraph.size(); i++){
    v_AllTGraph[i]->Write();
  }

  std::cout<<"Contained cuts:        "<<nContainedCut<<std::endl
           <<"Truth-reco match cuts: "<<nMatchCut<<std::endl
           <<"End point cuts:        "<<nEndPtCut<<std::endl
           <<"Num Daughters cuts:    "<<nDaughtersCut<<std::endl;

}
