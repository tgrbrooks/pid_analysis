//Quick and dirty truth-reco PID comparison

//Use a flag to select which event to use
int evt_index = 0;
int exit_flag = 0;

//Called at start of event loop
void calc_evdMC_Begin(TTree*){

  evt_index = 0;
  exit_flag = 0;

}

int calc_evdMC() {

  //Look at first event where PID goes wrong
  if(/*pdg[0]!=trkpidpdg_pmalgtrackmaker[0][2] &&*/ evt_index<200 && !exit_flag){ 
    
    //Create a canvas and divide into 2, half for MC and half for reco
    TCanvas *c1 = new TCanvas("c1","evd",900,500);
    c1->Divide(2);
    c1->cd(1);
    // Creating a view for MC, size is size of one TPC
    Double_t rmin[3] = {0,-200,0};
    Double_t rmax[3] = {200,200,500};
    TView3D *view = new TView3D(1,rmin,rmax);
    view->ShowAxis();

    TPolyMarker3D *vertex = new TPolyMarker3D(1);
    vertex->SetPoint(0, StartPointx[0], StartPointy[0], StartPointz[0]);
    vertex->SetMarkerStyle(8);
    vertex->SetMarkerSize(1);
    vertex->SetMarkerColor(30);
    vertex->Draw();

    std::cout<<std::endl<<"Event"<<event<<" MC Truth:"<<std::endl;

    int nMerged = 0;
    
    //Loop over g4 particles
    for (int i = 0; i<geant_list_size; i++){
      //Get rid of any that don't have a defined start/end
      //if (abs(StartPointx[i])>1000||abs(EndPointx[i])>1000) continue;
      if (pdg[i]==-99999) continue;

      if (MergedId[i]==MergedId[0]) nMerged++;

      // Create a PolyLine3D for each particle, plot as straight line between start and end
      TPolyLine3D *part = new TPolyLine3D(2);
      part->SetPoint(0, StartPointx[i], StartPointy[i], StartPointz[i]);
      part->SetPoint(1, EndPointx[i], EndPointy[i], EndPointz[i]);

      //Make the primary particle stand out
      bool isPrimary = 0;
      //Count scattered particles as also primary
      if (MergedId[i]==MergedId[0]) isPrimary = 1; 
      if (!isPrimary) {part->SetLineWidth(1);}// part->SetLineStyle(5);}

      //Uncomment to get rid of mess of neutrons and nuclear fragments
      //if(pdg[i]==2112||pdg[i]>10000) continue;

      // Output info about each particle
      double range = TMath::Sqrt(pow(StartPointx[i]-EndPointx[i],2)+pow(StartPointy[i]-EndPointy[i],2)+pow(StartPointz[i]-EndPointz[i],2));

      //Don't show particles that don't cross two wires
      if (range < 0.3) continue;

      std::cout<<TrackId[i]<<": pdg = "<<pdg[i]<<" ["<<P[i]<<" GeV]   ("<<StartPointx[i]<<","<<StartPointy[i]<<","<<StartPointz[i]<<") -> ("<<EndPointx[i]<<","<<EndPointy[i]<<","<<EndPointz[i]<<") "<<range<<" cm, Mother = "<<Mother[i]<<", Merged ID = "<<MergedId[i]<<std::endl;

      // Colour code the particles by type
      if (abs(pdg[i])==13) part->SetLineColor(4); //muons are blue
      if (abs(pdg[i])==11) part->SetLineColor(6); //electrons are pink
      if (abs(pdg[i])==211) part->SetLineColor(2); //pions are red
      if (pdg[i]==2212) part->SetLineColor(8); //protons are green
      if (pdg[i]==2112) {part->SetLineStyle(kDashed); part->SetLineWidth(1); continue;} //neutrons are black and dashed
      if (pdg[i]==22) {part->SetLineStyle(kDashed); part->SetLineColor(kYellow+1);} //photons are yellow and dashed
      
      //Don't draw really long tracks
      if (range > 1000) continue;
      // Draw
      part->Draw();
    }

    //Do exactly the same but for reconstructed particles
    c1->cd(2);
    // Creating a view
    TView3D *view_reco = new TView3D(1,rmin,rmax);
    view_reco->ShowAxis();

    vertex->Draw();

    std::cout<<std::endl<<"Reco:"<<std::endl;

bool isValid = 0;
bool isWrong = 0;

    //Loop over reconstructed tracks
    for (int i = 0; i<ntracks_pmalgtrackmaker; i++){
      
      //Don't show particles with wrong start positions
      //if (abs(trkstartx_pmalgtrackmaker[i])>1000||abs(trkendx_pmalgtrackmaker[i])>1000) continue;

      // Create a first PolyLine3D
      TPolyLine3D *part = new TPolyLine3D(2);
      part->SetPoint(0, trkstartx_pmalgtrackmaker[i], trkstarty_pmalgtrackmaker[i], trkstartz_pmalgtrackmaker[i]);
      part->SetPoint(1, trkendx_pmalgtrackmaker[i], trkendy_pmalgtrackmaker[i], trkendz_pmalgtrackmaker[i]);
      //if (trkidtruth_pmalgtrackmaker[i][2]!=1) {part->SetLineWidth(2); part->SetLineStyle(5);}

      //Get rid of invalid tracks
      if (trkidtruth_pmalgtrackmaker[i][2]==-1) continue;

      int g4i = 0;
      for (int j = 0; j<geant_list_size; j++){
        if (trkidtruth_pmalgtrackmaker[i][2]==TrackId[j]){ g4i =j; }
      }   
      //Skip if track not in geant list
      if (g4i == -1 ) continue;

      //Loop over hits in tracks
      double chi2 = 0;
      int npt = 0;
      for (int hit_i = 0; hit_i < ntrkhits_pmalgtrackmaker[i][2]; hit_i++){
        //Calculate the chi2 with a horizontal line at 3 MeV/cm
        if (trkdedx_pmalgtrackmaker[i][2][hit_i]>1000) continue;
        if (trkresrg_pmalgtrackmaker[i][2][hit_i]>25) continue;
        // Needs to be done over the the same range as dedx_range_pro
        double err = 0.04231+0.0001783*trkdedx_pmalgtrackmaker[i][2][hit_i]*trkdedx_pmalgtrackmaker[i][2][hit_i];
        err *= trkdedx_pmalgtrackmaker[i][2][hit_i];
        chi2 += pow((trkdedx_pmalgtrackmaker[i][2][hit_i]-3)/std::sqrt(0.7+pow(err,2)),2); // Chi2 with horizontal line at 3+/-0.7 MeV/cm
        npt++;
      }
      if (chi2/npt>trkpidchi_pmalgtrackmaker[i][2]) std::cout<<"VAL: ";
      if (chi2/npt>trkpidchi_pmalgtrackmaker[i][2]&&MergedId[g4i]==MergedId[0]) isValid = 1;
      if (trkidtruth_pmalgtrackmaker[i][2]==1&&trkpidpdg_pmalgtrackmaker[i][2]!=trkpdgtruth_pmalgtrackmaker[i][2]) isWrong = 1;
      if (!isValid) part->SetLineWidth(1);

      //Use PID from V plane for no real reason, could change to best plane
      double range = TMath::Sqrt(pow(trkstartx_pmalgtrackmaker[i]-trkendx_pmalgtrackmaker[i],2)+pow(trkstarty_pmalgtrackmaker[i]-trkendy_pmalgtrackmaker[i],2)+pow(trkstartz_pmalgtrackmaker[i]-trkendz_pmalgtrackmaker[i],2));
      std::cout<<trkidtruth_pmalgtrackmaker[i][2]<<": pdg = "<<trkpidpdg_pmalgtrackmaker[i][2]<<" ("<<trkpdgtruth_pmalgtrackmaker[i][2]<<"):  ("<<trkstartx_pmalgtrackmaker[i]<<","<<trkstarty_pmalgtrackmaker[i]<<","<<trkstartz_pmalgtrackmaker[i]<<") -> ("<<trkendx_pmalgtrackmaker[i]<<","<<trkendy_pmalgtrackmaker[i]<<","<<trkendz_pmalgtrackmaker[i]<<") "<<range<<" cm"<<std::endl;

      //Colour particles by type
      if (trkpidpdg_pmalgtrackmaker[i][2]==13) part->SetLineColor(4); //muons are blue
      if (trkpidpdg_pmalgtrackmaker[i][2]==211) part->SetLineColor(2); //pions are red
      if (trkpidpdg_pmalgtrackmaker[i][2]==2212) part->SetLineColor(8); //protons are green
      if (trkpidpdg_pmalgtrackmaker[i][2]==321) part->SetLineColor(kOrange-3); //kaons are orange
      // Draw
      part->Draw();
    }
    // Save to a root file
    c1->SaveAs("evdMC_muon1.root");
    c1->Modified();
    c1->Update(); 
    // Wait for user input so canvas is displayed
    if (isValid&&isWrong) std::cin>>exit_flag;

    evt_index++;

  }

  return 0;
}

void calc_pideff_Terminate(){

}
