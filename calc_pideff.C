//INCLUDES
//#include <vector>
//#include "TH2F.h"

//DEFINITIONS
TH2F *h_PidEnergy[3][3]; // PID vs energy for each plane and particle
/*
// Chi2 for each plane and particle
TH2F *h_Chi2p[3][3]; //for proton
TH2F *h_Chi2K[3][3]; //for kaon
TH2F *h_Chi2pi[3][3]; //for pion
TH2F *h_Chi2mu[3][3]; //for muon
*/
// PIDA vs energy for each plane and particle
//TH2F *h_PidaEnergy[3][3];
TH1F *h_Pida[3][3];

// PID efficiency vs energy for each plane and particle
TGraphAsymmErrors *h_EffEnergy[3][3];
TH1F *h_EffEnergy_pass[3][3]; // Count correct IDs
TH1F *h_EffEnergy_tot[3][3]; // Count all IDs, divide later

double efficiency[3][3];
double total[3][3];
double efficiency_cut[3][3];
double total_cut[3][3];

// Vectors to store all histograms and graphs
std::vector<TH2*> v_AllTH2;
std::vector<TH1*> v_AllTH1;
std::vector<TGraphAsymmErrors*> v_AllTGraph;

int nContainedCut[3];
int nMatchCut[3];
int nEndPtCut[3];
int nDaughtersCut[3];

TH1F *h_chi2null;
int nCut[3];
TH2F *h_dEdxVsResRg[3][3];
TH2F *h_dEdxVsResRg_cut[3][3];

void calc_pideff_Begin(TTree*){

  h_chi2null = new TH1F("hchi","",100,0,500);
  v_AllTH1.push_back(h_chi2null);

  //Loop over planes
  for (int i = 0; i < 3; i++){

    //Loop over particles, muon, pion, proton
    for (int j = 0; j < 3; j++){

      double max_energy = 0.0; std::string particle = "other"; int col = 1;
      if(j==0) {max_energy = 1.3; particle = "muon"; col = 4;} //Muons
      if(j==1) {max_energy = 2.3; particle = "pion"; col = 2;} //Pions
      if(j==2) {max_energy = 2.3; particle = "proton"; col = 8;} //Protons

      TString h_PidEnergy_name = Form("h_PidEnergy_plane%i_%s",i,particle.c_str()); 
      h_PidEnergy[i][j] = new TH2F(h_PidEnergy_name,"",10,0,max_energy,10,0,5);
      h_PidEnergy[i][j]->GetXaxis()->SetTitle("Momentum (GeV)");
      h_PidEnergy[i][j]->GetYaxis()->SetTitle("#mu      #pi      K       p     ");
      h_PidEnergy[i][j]->GetYaxis()->SetTitleOffset(0.3);
      h_PidEnergy[i][j]->GetYaxis()->SetLabelSize(0.0);
      v_AllTH2.push_back(h_PidEnergy[i][j]);
/*    
      TString h_PidaEnergy_name = Form("h_PidaEnergy_plane%i_%s",i,particle.c_str()); 
      h_PidaEnergy[i][j] = new TH2F(h_PidaEnergy_name,"",10,0,max_energy,100,0,30);
      h_PidaEnergy[i][j]->GetXaxis()->SetTitle("Momentum (GeV)");
      h_PidaEnergy[i][j]->GetYaxis()->SetTitle("PIDA Value");
      v_AllTH2.push_back(h_PidaEnergy[i][j]);
*/    
      TString h_Pida_name = Form("h_Pida_plane%i_%s",i,particle.c_str()); 
      h_Pida[i][j] = new TH1F(h_Pida_name,"",100,0,30);
      h_Pida[i][j]->GetXaxis()->SetTitle("PIDA Value");
      h_Pida[i][j]->GetYaxis()->SetTitle("Number Events");
      v_AllTH1.push_back(h_Pida[i][j]);
/*    
      TString h_Chi2p_name = Form("h_Chi2p_plane%i_%s",i,particle.c_str());
      h_Chi2p[i][j] = new TH2F(h_Chi2p_name,"",10,0,max_energy,50,0,500);
      h_Chi2p[i][j]->GetXaxis()->SetTitle("Momentum (GeV)");
      h_Chi2p[i][j]->GetYaxis()->SetTitle("#chi^{2} (p)");
      v_AllTH2.push_back(h_Chi2p[i][j]);
      TString h_Chi2K_name = Form("h_Chi2K_plane%i_%s",i,particle.c_str());
      h_Chi2K[i][j] = new TH2F(h_Chi2K_name,"",10,0,max_energy,50,0,500);
      h_Chi2K[i][j]->GetXaxis()->SetTitle("Momentum (GeV)");
      h_Chi2K[i][j]->GetYaxis()->SetTitle("#chi^{2} (K)");
      v_AllTH2.push_back(h_Chi2K[i][j]);
      TString h_Chi2pi_name = Form("h_Chi2pi_plane%i_%s",i,particle.c_str());
      h_Chi2pi[i][j] = new TH2F(h_Chi2pi_name,"",10,0,max_energy,50,0,500);
      h_Chi2pi[i][j]->GetXaxis()->SetTitle("Momentum (GeV)");
      h_Chi2pi[i][j]->GetYaxis()->SetTitle("#chi^{2} (#pi)");
      v_AllTH2.push_back(h_Chi2pi[i][j]);
      TString h_Chi2mu_name = Form("h_Chi2mu_plane%i_%s",i,particle.c_str());
      h_Chi2mu[i][j] = new TH2F(h_Chi2mu_name,"",10,0,max_energy,50,0,500);
      h_Chi2mu[i][j]->GetXaxis()->SetTitle("Momentum (GeV)");
      h_Chi2mu[i][j]->GetYaxis()->SetTitle("#chi^{2} (#mu)");
      v_AllTH2.push_back(h_Chi2mu[i][j]);
*/
      TString h_EffEnergy_name = Form("h_EffEnergy_plane%i_%s",i,particle.c_str());
      h_EffEnergy[i][j] = new TGraphAsymmErrors();
      h_EffEnergy[i][j]->SetName(h_EffEnergy_name);
      h_EffEnergy[i][j]->SetTitle(";Momentum (GeV);PID Efficiency (%)");
      h_EffEnergy[i][j]->SetMarkerStyle(8);
      h_EffEnergy[i][j]->SetMarkerColor(col);
      h_EffEnergy[i][j]->SetLineColor(col);
      h_EffEnergy[i][j]->SetLineWidth(3);
      h_EffEnergy[i][j]->GetYaxis()->SetRangeUser(0,1);
      h_EffEnergy[i][j]->Draw("P");
      v_AllTGraph.push_back(h_EffEnergy[i][j]);
     
      TString h_EffEnergy_pass_name = Form("h_EffEnergy_pass_plane%i_%s",i,particle.c_str());
      h_EffEnergy_pass[i][j] = new TH1F(h_EffEnergy_pass_name,"",10,0,max_energy);
      TString h_EffEnergy_tot_name = Form("h_EffEnergy_tot_plane%i_%s",i,particle.c_str());
      h_EffEnergy_tot[i][j] = new TH1F(h_EffEnergy_tot_name,"",10,0,max_energy);
     
      efficiency[i][j]=0.0;
      total[i][j]=0.0;
      efficiency_cut[i][j]=0.0;
      total_cut[i][j]=0.0;

      TString h_dEdxVsResRg_name = Form("h_dEdxVsResRg_plane%i_%i",i,j);
      h_dEdxVsResRg[i][j] = new TH2F(h_dEdxVsResRg_name,"",100,0,25,100,0,30);
      h_dEdxVsResRg[i][j]->GetXaxis()->SetTitle("Residual range (cm)");
      h_dEdxVsResRg[i][j]->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
      v_AllTH2.push_back(h_dEdxVsResRg[i][j]);

      TString h_dEdxVsResRg_cut_name = Form("h_dEdxVsResRg_cut_plane%i_%i",i,j);
      h_dEdxVsResRg_cut[i][j] = new TH2F(h_dEdxVsResRg_cut_name,"",100,0,25,100,0,30);
      h_dEdxVsResRg_cut[i][j]->GetXaxis()->SetTitle("Residual range (cm)");
      h_dEdxVsResRg_cut[i][j]->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
      v_AllTH2.push_back(h_dEdxVsResRg_cut[i][j]);

    }

    nContainedCut[i]=0;
    nMatchCut[i]=0;
    nEndPtCut[i]=0;
    nDaughtersCut[i]=0;

    nCut[i]=0;

  }

}

int calc_pideff() {

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

  //Check how many g4 tracks correspond to the incoming particle
  int nMerged = 0;
  for (int i = 0; i < geant_list_size; i++){
    if (MergedId[i]==MergedId[0]) nMerged++;
  }
  
  // Produce dE/dx vs res range curves and calculate efficiency before applying cuts
  for (int track_i = 0; track_i < ntracks_pmalgtrackmaker; track_i++){

    for (int plane_i = 0; plane_i < 3; plane_i++){

      //Uncomment to only look at primary track
      //if (trkidtruth_pmalgtrackmaker[track_i][plane_i]!=1) continue;
      int g4i = -1;
      for (int j = 0; j < geant_list_size; j++){
        if (trkidtruth_pmalgtrackmaker[track_i][plane_i]==TrackId[j]){ g4i =j; }
      }
      //Only primary particle
      if (g4i == -1 ) continue;
      if (MergedId[g4i]!=MergedId[0]) continue;

      int iPart = -1;
      if (trkpdgtruth_pmalgtrackmaker[track_i][plane_i]==13) iPart=0;
      else if (trkpdgtruth_pmalgtrackmaker[track_i][plane_i]==211) iPart=1;
      else if (trkpdgtruth_pmalgtrackmaker[track_i][plane_i]==2212) iPart=2;
      
      // dE/dx vs res range
      for (int hit_i = 0; hit_i < ntrkhits_pmalgtrackmaker[track_i][plane_i]; hit_i++){
        if(iPart != -1 && trkidtruth_pmalgtrackmaker[track_i][plane_i]!=-1 && trkresrg_pmalgtrackmaker[track_i][plane_i][hit_i]>0.1){
          h_dEdxVsResRg[plane_i][iPart]->Fill(trkresrg_pmalgtrackmaker[track_i][plane_i][hit_i],trkdedx_pmalgtrackmaker[track_i][plane_i][hit_i]);
        }
      }

      // Efficiency
      if( iPart!=-1 ){
        if (trkpidpdg_pmalgtrackmaker[track_i][plane_i]==trkpdgtruth_pmalgtrackmaker[track_i][plane_i]){
          efficiency[plane_i][iPart] += 1;
        } 
        total[plane_i][iPart] += 1;
      }
    }   
  }

  
/*
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
    else{
      if (pdg[0]==13) nContainedCut[0]++;
      if (pdg[0]==211) nContainedCut[1]++;
      if (pdg[0]==2212) nContainedCut[2]++;
    }
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
  if (NRecoTracksTruth!=1){
      if (pdg[0]==13) nMatchCut[0]++;
      if (pdg[0]==211) nMatchCut[1]++;
      if (pdg[0]==2212) nMatchCut[2]++;
  }

  //Apply selection
  if (!(isContained && NRecoTracksTruth==1)) return 0;
  //Demand track orientated the same way as the true particle
  TVector3 TrueStart(StartPointx[0],StartPointy[0],StartPointz[0]);
  TVector3 TrueEnd(EndPointx[0],EndPointy[0],EndPointz[0]);
  TVector3 RecoStart(trkstartx_pmalgtrackmaker[recoIndex],trkstarty_pmalgtrackmaker[recoIndex],trkstartz_pmalgtrackmaker[recoIndex]);
  TVector3 RecoEnd(trkendx_pmalgtrackmaker[recoIndex],trkendy_pmalgtrackmaker[recoIndex],trkendz_pmalgtrackmaker[recoIndex]);

  if (!((RecoStart-TrueStart).Mag() < (RecoEnd-TrueStart).Mag() && (RecoEnd-TrueEnd).Mag() < (RecoStart-TrueEnd).Mag())) return 0;

  if (abs(pdg[0])==211){
    if ((RecoEnd-TrueEnd).Mag() > 1.){nEndPtCut[1]++; return 0;}
    if (NumberDaughters[0] > 3){nDaughtersCut[1]++; return 0;}
  }
  else if (pdg[0]==2212){
    if ((RecoEnd-TrueEnd).Mag() > 0.3){nEndPtCut[2]++; return 0;}
    if (NumberDaughters[0] > 3){nDaughtersCut[2]++; return 0;}
  }
  else if (pdg[0]==13){
    if ((RecoEnd-TrueEnd).Mag() > 1.){nEndPtCut[0]++; return 0;}

  }
  else {
    if ((RecoEnd-TrueEnd).Mag() > 1.){return 0;}
  } 
*/

  // Calculate the chi2/ndf for a fit to a horizontal line to check if chi2 method is valid
  std::vector<std::vector<double>> chi2null;
  std::vector<std::vector<TGraph*>> v_dedxresrg;
  //Loop over reco tracks in event
  for (int track_i = 0; track_i < ntracks_pmalgtrackmaker; track_i++){
    //Loop over planes
    std::vector<double> chi2planes;
    std::vector<TGraph*> v_dedxresrg_pln;
    for (int plane_i = 0; plane_i < 3; plane_i++){
      //Loop over hits in tracks
      double chi2 = 0;
      int npt = 0;

      std::vector<double> dedx;
      std::vector<double> resrg;

      for (int hit_i = 0; hit_i < ntrkhits_pmalgtrackmaker[track_i][plane_i]; hit_i++){
        //Calculate the chi2 with a horizontal line at 3 MeV/cm
        if (trkdedx_pmalgtrackmaker[track_i][plane_i][hit_i]>1000) continue;
        if (trkresrg_pmalgtrackmaker[track_i][plane_i][hit_i]>25 || trkresrg_pmalgtrackmaker[track_i][plane_i][hit_i]<0.1) continue;
        // Needs to be done over the the same range as dedx_range_pro
        double err = 0.04231+0.0001783*trkdedx_pmalgtrackmaker[track_i][plane_i][hit_i]*trkdedx_pmalgtrackmaker[track_i][plane_i][hit_i];
        err *= trkdedx_pmalgtrackmaker[track_i][plane_i][hit_i];
        chi2 += pow((trkdedx_pmalgtrackmaker[track_i][plane_i][hit_i]-3)/std::sqrt(0.7+pow(err,2)),2); // Chi2 with horizontal line at 3+/-0.7 MeV/cm
        npt++;

        dedx.push_back(trkdedx_pmalgtrackmaker[track_i][plane_i][hit_i]);
        resrg.push_back(trkresrg_pmalgtrackmaker[track_i][plane_i][hit_i]);
  
      }   
      h_chi2null->Fill(chi2/npt);
      chi2planes.push_back(chi2/npt);

      //Create TGraphs of dEdx vs res range and fit with exponentials and 1st order polynomials
      TGraph *g_dedxresrg = new TGraph(dedx.size(),&(resrg[0]),&(dedx[0]));
      v_dedxresrg_pln.push_back(g_dedxresrg);
    }   
    chi2null.push_back(chi2planes);
    v_dedxresrg.push_back(v_dedxresrg_pln);
  }

  // Loop over all the tracks
  for (int i = 0; i<ntracks_pmalgtrackmaker; i++){
    //Here there is only one primary particle, in this case the primary particle is always top of the list
    for (int iPln=0; iPln<3; iPln++){
 
      //Check track corresponds to a true particle
      if (trkidtruth_pmalgtrackmaker[i][iPln]==-1) continue;
      //if (nMerged>1) continue;

      //Uncomment to only look at primary track
      //if (trkidtruth_pmalgtrackmaker[i][iPln]!=1) continue;
 
      //Uncomment to remove tracks which don't match the true length, to avoid joined tracks and delta rays
      double true_range = TMath::Sqrt(pow(StartPointx[0]-EndPointx[0],2)+pow(StartPointy[0]-EndPointy[0],2)+pow(StartPointz[0]-EndPointz[0],2));
      double reco_range = TMath::Sqrt(pow(trkstartx_pmalgtrackmaker[i]-trkendx_pmalgtrackmaker[i],2)+pow(trkstarty_pmalgtrackmaker[i]-trkendy_pmalgtrackmaker[i],2)+pow(trkstartz_pmalgtrackmaker[i]-trkendz_pmalgtrackmaker[i],2));
      //if (abs(true_range-reco_range)>10) continue;
 
      //If the smallest chi2 is larger than the chi2 of and invalid track skip the event
      if(chi2null[i][iPln]<trkpidchi_pmalgtrackmaker[i][iPln]){nCut[iPln]++; continue;}
      //TF1 *fexp = new TF1("fexp","expo",0,10);
      //TF1 *flin = new TF1("flin","pol1",0,10);
      //v_dedxresrg[i][iPln]->Fit("fexp","QRN");
      //v_dedxresrg[i][iPln]->Fit("flin","QRN");
      //if(fexp->GetChisquare()>flin->GetChisquare()){nCut[iPln]++; continue;}
 
      //Convert the true pdg code to an index
      int iPart = -1;
      if (abs(trkpdgtruth_pmalgtrackmaker[i][iPln])==13) iPart=0;
      else if (abs(trkpdgtruth_pmalgtrackmaker[i][iPln])==211) iPart=1;
      else if (abs(trkpdgtruth_pmalgtrackmaker[i][iPln])==2212) iPart=2;
 
      // Only analysing muons, pions and protons
      if( iPart!=-1 ){
 
        //Get the corresponding true momentum
        double mom = 0;
        int g4i = -1;
        for (int j = 0; j < geant_list_size; j++){
          if (trkidtruth_pmalgtrackmaker[i][iPln]==TrackId[j]){ mom = P[j]; g4i =j; }
        }
        //Only primary particle
        if (g4i == -1) continue;
        //if (MergedId[g4i]!=MergedId[0]) continue;
 
        int pdgcode = trkpidpdg_pmalgtrackmaker[i][iPln];
        if (pdgcode>0){
          //Define different codes for better looking plot (muon = 1, pion = 2, kaon = 3, proton = 4)
          int pidcode = 0;
          if (pdgcode == 13) pidcode = 1;
          else if (pdgcode == 211) pidcode = 2;
          else if (pdgcode == 321) pidcode = 3;
          else if (pdgcode == 2212) pidcode = 4;
          // Fill PID vs momentum histogram
          h_PidEnergy[iPln][iPart]->Fill(mom,pidcode);
        }
/* 
        // Fill chi2 vs momentum histograms for each particle comparison
        if (trkpidchipr_pmalgtrackmaker[i][iPln]>0) h_Chi2p[iPln][iPart]->Fill(mom,trkpidchipr_pmalgtrackmaker[i][iPln]);
        if (trkpidchika_pmalgtrackmaker[i][iPln]>0) h_Chi2K[iPln][iPart]->Fill(mom,trkpidchika_pmalgtrackmaker[i][iPln]);
        if (trkpidchipi_pmalgtrackmaker[i][iPln]>0) h_Chi2pi[iPln][iPart]->Fill(mom,trkpidchipi_pmalgtrackmaker[i][iPln]);
        if (trkpidchimu_pmalgtrackmaker[i][iPln]>0) h_Chi2mu[iPln][iPart]->Fill(mom,trkpidchimu_pmalgtrackmaker[i][iPln]);
*/
        // Calculate efficiency after cuts
        if (trkpidpdg_pmalgtrackmaker[i][iPln]==abs(trkpdgtruth_pmalgtrackmaker[i][iPln])){
          h_EffEnergy_pass[iPln][iPart]->Fill(mom);
          efficiency_cut[iPln][iPart] += 1;
        }
        h_EffEnergy_tot[iPln][iPart]->Fill(mom);
        total_cut[iPln][iPart] += 1;
        //Fill PIDA distributions
//        h_PidaEnergy[iPln][iPart]->Fill(mom,trkpidpida_pmalgtrackmaker[i][iPln]);
        h_Pida[iPln][iPart]->Fill(trkpidpida_pmalgtrackmaker[i][iPln]);
 
        // dE/dx vs res range curve for tracks which don't pass the cuts
        for (int hit_i = 0; hit_i < ntrkhits_pmalgtrackmaker[i][iPln]; hit_i++){
          if(trkresrg_pmalgtrackmaker[i][iPln][hit_i] > 0){
            h_dEdxVsResRg_cut[iPln][iPart]->Fill(trkresrg_pmalgtrackmaker[i][iPln][hit_i],trkdedx_pmalgtrackmaker[i][iPln][hit_i]);
          }
        }

      }
    }
  }

  return 0;
}

void calc_pideff_Terminate(){

  for (int i = 0; i<3; i++){
    for (int j = 0; j<3; j++){
      h_EffEnergy[i][j]->BayesDivide(h_EffEnergy_pass[i][j],h_EffEnergy_tot[i][j]);
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

  for (int i = 0; i < 3; i++){
    std::cout<<"["<<i<<"]:Contained cuts:        "<<nContainedCut[i]<<std::endl
             <<"Truth-reco match cuts: "<<nMatchCut[i]<<std::endl
             <<"End point cuts:        "<<nEndPtCut[i]<<std::endl
             <<"Num Daughters cuts:    "<<nDaughtersCut[i]<<std::endl
             <<"Invalid cuts:          "<<nCut[i]<<std::endl;
  }

}
