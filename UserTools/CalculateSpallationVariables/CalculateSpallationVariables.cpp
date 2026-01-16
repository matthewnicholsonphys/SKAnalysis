#include "CalculateSpallationVariables.h"

#include "MTreeReader.h"

#include "TH1D.h"
#include "geotnkC.h"  // for SK tank geometric constants

CalculateSpallationVariables::CalculateSpallationVariables():Tool(){}

bool CalculateSpallationVariables::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
  
  GetReaders();

  dt_bins = MakeLogBins(0.001, 60, nbins_log+1);
  dlt_bins = MakeLogBins(1, 5000, nbins_log+1);
  like_bins = MakeLogBins(0.01, 20, nbins_log_like+1);
  pre_dt_hist = TH1D("pre_dt_hist", "pre_dt_hist", nbins_log, dt_bins.data());
  post_dt_hist = TH1D("post_dt_hist", "post_dt_hist", nbins_log, dt_bins.data());

  CreateOutputFile();
  
  return true;

}

bool CalculateSpallationVariables::Execute(){

  GetRelicBranchValues(); // get relic branch values

  // "pre" means dt > 0, "post" is opposite - remember these are the times from the muon to the relic, so it's the opposite to the white paper.

  v_pairing_info.clear();
  v_dt.clear();
  v_dlt.clear();
  v_dll.clear();
  v_muqismsk.clear();
  v_resQ.clear();
  v_muontype.clear();
  v_bsenergy.clear();
  v_do_check.clear();


  //std::cout << "CalculateSpallationVariables: looping through " << MatchedTimeDiff_ptr->size() << " matched relics" << std::endl;

  // if size is zero, throw exception, we should always have lots of muons matched to a relic
  // if (MatchedTimeDiff_ptr->size() == 0){
  //   //throw std::runtime_error("CalculateSpallationVariables::Execute: relic not matched to a single muon! problem!");
  // }
  
  for (int muon_idx = 0; muon_idx < MatchedTimeDiff_ptr->size(); ++muon_idx){

    /* dt - time difference between muon and relic candidate */
    /* mu - relic */
    float dt = MatchedTimeDiff_ptr->at(muon_idx)/pow(10,9);
    //dt < 0 ? pre_dt_hist.Fill(abs(dt)) : post_dt_hist.Fill(abs(dt)); // change all signs .. (we are not fallible)
    
    muon_tree_ptr->GetEntry(MatchedOutEntryNums_ptr->at(muon_idx)); // get muons 
    GetMuonBranchValues();	// get muon values

    // do I need to return for badly reconstructed muons? yeah, probably - only when we use `all'
    if (MU_ptr->muboy_status == 0){

      v_pairing_info.push_back("");
      v_dt.push_back(0);
      v_dlt.push_back(0);
      v_dll.push_back(0);
      v_muqismsk.push_back(0);
      v_resQ.push_back(0);
      v_muontype.push_back(0);
      v_bsenergy.push_back(0);
      v_do_check.push_back(false);

      continue;
    }
    
    if (LOWE_ptr->bsenergy > 1000){
      //bad reconstruction, skipping

      v_pairing_info.push_back("");
      v_dt.push_back(0);
      v_dlt.push_back(0);
      v_dll.push_back(0);
      v_muqismsk.push_back(0);
      v_resQ.push_back(0);
      v_muontype.push_back(0);
      v_bsenergy.push_back(0);
      v_do_check.push_back(false);
      
      continue;
    }

    /* dlt - transverse distance between muon and relic candidate */
    float dlt = 0, appr = 0;

    float* muon_entrypoint = nullptr;
    float* muon_direction = nullptr;
    
    bool did_bff = MU_ptr->muinfo[6];
    double muon_tracklen = 0;
    
    if (did_bff){ // we never did do bff
      basic_array<float> bff_entrypoint(MU_ptr->mubff_entpos);
      basic_array<float> bff_dir = MU_ptr->mubff_dir;

      muon_entrypoint = const_cast<float*>(bff_entrypoint.data());
      muon_direction = const_cast<float*>(bff_dir.data());

      muon_tracklen = CalculateTrackLen(muon_entrypoint, muon_direction);
    } else {
      basic_array<float[10][4]> muboy_entrypoint(MU_ptr->muboy_entpos);
      basic_array<float> muboy_dir = MU_ptr->muboy_dir;

      int muboy_idx = MU_ptr->muinfo[7];
      muon_entrypoint = const_cast<float*>(muboy_entrypoint[muboy_idx].data());
      muon_direction = const_cast<float*>(muboy_dir.data());

      muon_tracklen = muboy_idx == 0 ? MU_ptr->muboy_length : CalculateTrackLen(muon_entrypoint, muon_direction);
      
    }
    
    float* relic_pos = const_cast<float*>(LOWE_ptr->bsvertex);
    
    getdl_(muon_direction,
	   &relic_pos[0],
	   &relic_pos[1],
	   &relic_pos[2],
	   muon_entrypoint,
	   &dlt,
	   &appr);

    if (dlt == 0.0){
      std::cout << "CalculateSpallationVariables::Execute - dlt = 0, dumping args of getdl_" << std::endl;
      std::cout<<"calling getdl_ with:\n"
	       <<"\trelic po: ("<<relic_pos[0]<<", "<<relic_pos[1]<<", "<<relic_pos[2]<<")\n"
	       <<"\tmuon entry point: ("<<muon_entrypoint[0]<<", "<<muon_entrypoint[1]<<", "<<muon_entrypoint[2]<<")\n"
	       <<"\tmuon entry dir: ("<<muon_direction[0]<<", "<<muon_direction[1]<<", "<<muon_direction[2]<<")"<<std::endl;
      throw std::runtime_error("CalculateSpallationVariables:: bad dlt");
    }
    
    // dt < 0 ? pre_dlt_hist.Fill(dlt) : post_dlt_hist.Fill(dlt);
    
    /* dll - longitudinal distance between the muon and relic candidate*/
    // basic_array<float> scott_dedx(MU_ptr->muboy_dedx);
    // float* muon_dedx = const_cast<float*>(scott_dedx.data());

    const float* muon_dedx = MU_ptr->muboy_dedx;

    double max_edep = 0;
    int max_edep_bin=0;
    for(int i=0;i<111;i++){
      double e_dep_in_window = 0 ;
      for(int j=0;j<9;j++){
	e_dep_in_window = e_dep_in_window + muon_dedx[i+j];
      }
      if(e_dep_in_window > max_edep){
	max_edep_bin = i+4;
	max_edep = e_dep_in_window;
      }
    }
    double max_energy_dep_pos = 50.*max_edep_bin;
    float dll = max_energy_dep_pos - appr;

    if (dll == 0.0){
      std::cout << "CalculateSpallationVariables::Execute - dll = 0, dumping args of getdl_" << std::endl;
      std::cout<<"calling getdl_ with:\n"
	       <<"\trelic po: ("<<relic_pos[0]<<", "<<relic_pos[1]<<", "<<relic_pos[2]<<")\n"
	       <<"\tmuon entry point: ("<<muon_entrypoint[0]<<", "<<muon_entrypoint[1]<<", "<<muon_entrypoint[2]<<")\n"
	       <<"\tmuon entry dir: ("<<muon_direction[0]<<", "<<muon_direction[1]<<", "<<muon_direction[2]<<")"<<")\n"
	       <<"\tmax_energy_dp_pos: "<< max_energy_dep_pos<< ", appr: "<<appr<<std::endl;
      throw std::runtime_error("CalculateSpallationVariables:: bad dll");
    }


    // dt < 0 ? pre_dll_hist.Fill(dll) : post_dll_hist.Fill(dll); // sign changes

    /* muqismsk - max charge deposited in the detector by the muon*/
    float muqismsk = (MU_ptr->muqismsk);
    // dt < 0 ? pre_muqismsk_hist.Fill(muqismsk) : post_muqismsk_hist.Fill(muqismsk);
    
    /* 
       resQ - residual charge deposited by the muon compared to the value expected from the min ionization. 
       reQ = muqismsk - q_MI * L, where q_MI is the number of photoelectrons per cm expected from the min ionization
       and L is the track length. L = sum_{i}(Li) for multiple tracks. Does that not make sense to you? yeah well get in line
    */

    double pe_per_coulomb = 30;     // this might come back to bite me but break in the sun till the sun breaks down old boy
    const float pe_per_cm = 26.78;
    double pe_from_muon = MU_ptr->muqismsk * (pe_per_cm / pe_per_coulomb);  // pe*cm^-1 / pe*C^-1 = C/cm
    double pe_from_MIP = muon_tracklen * pe_per_cm;
    float resQ = pe_from_muon - pe_from_MIP;
    // dt < 0 ? pre_resQ_hist.Fill(resQ) : post_resQ_hist.Fill(resQ);

    /* lastly, we need the type of muon event: misfit=0, single_through=1, single_stopping=2, multi=3,4, corner=5 */
    int muon_type = MU_ptr->muboy_status;

    /* oh we also need to the bsenergy: */
    float bse = LOWE_ptr->bsenergy;
    
    PairingInfo p = {dt, dlt, dll, muqismsk, resQ, muon_type, bse}; 
    std::string p_str = GetPairingString(p);
    //pairings[p_str].push_back(p);

    v_pairing_info.push_back(p_str);
    v_dt.push_back(dt);
    v_dlt.push_back(dlt);
    v_dll.push_back(dll);
    v_muqismsk.push_back(muqismsk);
    v_resQ.push_back(resQ);
    v_muontype.push_back(muon_type);
    v_bsenergy.push_back(bse);
    v_do_check.push_back(true); // check is false for badly reconstructuted muons and or relics, still save to preserve matching 

  }

  data_output_tree_ptr->Fill();
  
  return true;
}

bool CalculateSpallationVariables::Finalise(){

  // TTree* relic_tree_clone_ptr = relic_tree_ptr->GetTree()->CloneTree();
  
  data_output_file_ptr->cd();
  // relic_tree_clone_ptr->Write();
  data_output_tree_ptr->Write();

  data_output_file_ptr->Close();

  // std::string hist_outputfile_str = "";
  // m_variables.Get("hist_output_file_str", hist_outputfile_str);
  // if (hist_outputfile_str.empty()){throw std::runtime_error("CalculateSpallationVariables::Finalise(): no hist_output file specified!");}

  // TFile hist_output_file = TFile(hist_outputfile_str.c_str(), "RECREATE");
  // if (hist_output_file.IsZombie()){throw std::runtime_error("CalculateSpallationVariables::Finalise(): couldn't open hist_output file!");}
  // hist_output_file.cd();
  
  // for (int bin_idx = 1; bin_idx < pre_dt_hist.GetNbinsX()+1; ++bin_idx){
  //   pre_dt_hist.SetBinContent(bin_idx,
  // 			      pre_dt_hist.GetBinContent(bin_idx)/pre_dt_hist.GetBinWidth(bin_idx));
  //   post_dt_hist.SetBinContent(bin_idx,
  // 			      post_dt_hist.GetBinContent(bin_idx)/post_dt_hist.GetBinWidth(bin_idx));
  // }
  
  // pre_dt_hist.Write();
  // post_dt_hist.Write();
  // pre_dlt_hist.Write();
  // post_dlt_hist.Write();
  // pre_dll_hist.Write();
  // post_dll_hist.Write();
  // pre_muqismsk_hist.Write();
  // post_muqismsk_hist.Write();
  // pre_resQ_hist.Write();
  // post_resQ_hist.Write();

  //std::cout << "number of hists " << pairings.size() << std::endl;
  
  // for (const auto& [name, v_pairing] : pairings){
  //   std::cout << "name: " << name << std::endl;
  //   CreateLikelihood(name, v_pairing);
  // }
  
  // hist_output_file.Write();
  // hist_output_file.Close();
  
  return true;
}

void CalculateSpallationVariables::GetReaders(){
  std::string relic_reader_name = "";
  if (!m_variables.Get("relic_reader_name", relic_reader_name) || relic_reader_name.empty()){
    throw std::runtime_error("CalculateSpallationVariables::GetReader - no relic_reader_name specified!");
  }
  if (m_data->Trees.count(relic_reader_name) == 0){
    throw std::runtime_error("CalculateSpallationVariables::GetReader - relic reader not found!");
  }
  relic_tree_ptr = m_data->Trees.at(relic_reader_name);

  std::string muon_reader_name = "";
  if (!m_variables.Get("muon_reader_name", muon_reader_name) || muon_reader_name.empty()){
    throw std::runtime_error("CalculateSpallationVariables::GetReader - no muon_reader_name specified!");
  }
  if (m_data->Trees.count(muon_reader_name) == 0){
    throw std::runtime_error("CalculateSpallationVariables::GetReader - muon reader not found!");
  }
  muon_tree_ptr = m_data->Trees.at(muon_reader_name);

}

void CalculateSpallationVariables::GetRelicBranchValues(){
  bool ok = relic_tree_ptr->Get("MatchedTimeDiff", MatchedTimeDiff_ptr);
  if (!ok){throw std::runtime_error("CalculateSpallationVariables::GetRelicBranchValues: Couldn't retrieve MatchedTimeDiff!");}
  ok = relic_tree_ptr->Get("MatchedOutEntryNums", MatchedOutEntryNums_ptr);
  if (!ok){throw std::runtime_error("CalculateSpallationVariables::GetRelicBranchValues: Couldn't retrieve MatchedOutEntryNums!");}
  ok = relic_tree_ptr->Get("LOWE", LOWE_ptr);
  if (!ok){throw std::runtime_error("CalculateSpallationVariables::GetRelicBranchValues: Couldn't retrieve LOWE!!");}
  return;
}

void CalculateSpallationVariables::GetMuonBranchValues(){
  MU_ptr = nullptr;
  bool ok = muon_tree_ptr->Get("MU", MU_ptr);
  if (!ok || MU_ptr == nullptr){
    throw std::runtime_error("CalculateSpallationVariables::GetMuonBranchValues: can't get MU branch!");
  }
  return;
}

std::string CalculateSpallationVariables::GetPairingString(PairingInfo p) const {
  std::string hist_str = "";

  // return "all";
  
  const std::vector<std::string> type_strs = {"misfit", "singlethru", "singlestop", "multi", "multi", "corner"};
  hist_str+=type_strs.at(p.muon_type);

  //  return hist_str; //debug for now
  
  hist_str+="_";
  float bse = p.bse;
  if (bse > 8 && bse < 10){hist_str+="bse8-10";}
  else if (bse > 10 && bse < 12){hist_str+="bse10-12";}
  else if (bse > 12 && bse < 14){hist_str+="bse12-14";}
  else if (bse > 14 && bse < 16){hist_str+="bse14-16";}
  else if (bse > 16 && bse < 18){hist_str+="bse16-18";}
  else if (bse > 18 && bse < 20){hist_str+="bse18-20";}
  else if (bse > 20 && bse < 24){hist_str+="bse20-24";}

  hist_str+="_";
  const float dt = abs(p.dt);
  if (dt > 0 && dt <= 0.05){hist_str+="dt-short";
  } else if (dt > 0.05 && dt <= 0.5){hist_str+="dt-med";
  } else if (dt > 0.5/* && dt < 30*/){hist_str+="dt-long";}

  if (p.muon_type == 0){return hist_str;}

  const float dlt = p.dlt;
  if (dlt > 0 && dlt <= 300){hist_str+="_dlt-short";}
  else if (dlt > 300 && dlt <=  1000){hist_str+="_dlt-med";}
  else if (dlt > 1000/* && dlt < 1000*/){hist_str+="_dlt-long";}
  return hist_str;
}

void CalculateSpallationVariables::CreateLikelihood(const std::string& name, const std::vector<PairingInfo>& vp) const {
  // this isn't used but we keep it here if we ever want to do full splitting of muon/relic phase space.
  

  TH1D pre_dt("pre_dt", "pre_dt", nbins_log, dt_bins.data()), post_dt("post_dt", "post_dt", nbins_log, dt_bins.data());
  //TH1D pre_dlt("pre_dlt", "pre_dlt", nbins, 0, 5000), post_dlt("post_dlt", "post_dlt", nbins, 0, 5000);
  TH1D pre_dlt("pre_dlt", "pre_dlt", nbins_log, dlt_bins.data()), post_dlt("post_dlt", "post_dlt", nbins_log, dlt_bins.data());

  TH1D pre_dll("pre_dll", "pre_dll", nbins, -5000, 5000), post_dll("post_dll", "post_dll", nbins, -5000, 5000);
  TH1D pre_muqismsk("pre_muqismsk", "pre_muqismsk", nbins, 0, 250000), post_muqismsk("post_muqismsk", "post_muqismsk", nbins, 0, 250000);
  TH1D pre_resQ("pre_resQ", "pre_resQ", nbins, -100000, 1000000), post_resQ("post_resQ", "post_resQ", nbins, -100000, 1000000);

  
  for (const auto& p : vp){
    if (p.dt > 0){
      pre_dt.Fill(abs(p.dt));
      pre_dlt.Fill(p.dlt);
      pre_dll.Fill(p.dll);
      pre_muqismsk.Fill(p.muqismsk);
      pre_resQ.Fill(p.resQ);
    } else {
      post_dt.Fill(abs(p.dt));
      post_dlt.Fill(p.dlt);
      post_dll.Fill(p.dll);
      post_muqismsk.Fill(p.muqismsk);
      post_resQ.Fill(p.resQ);
    }
  }

  //come back to this when you get the damn likelihoods working
  // for (int bin_idx = 1; bin_idx < pre_dt.GetNbinsX()+1; ++bin_idx){
  //   pre_dt.SetBinContent(bin_idx,
  // 			      pre_dt.GetBinContent(bin_idx)/pre_dt.GetBinWidth(bin_idx));
  //   post_dt.SetBinContent(bin_idx,
  // 			      post_dt.GetBinContent(bin_idx)/post_dt.GetBinWidth(bin_idx));
  // }

  const double post_integral = post_dt.Integral();
  
  TH1D excess_dt("excess_dt", "excess;dt", nbins_log, dt_bins.data());
  for (int bin_idx = 1; bin_idx < excess_dt.GetNbinsX()+1; ++bin_idx){
    excess_dt.SetBinContent(bin_idx, (pre_dt.GetBinContent(bin_idx) / (post_integral * (post_dt.GetBinWidth(bin_idx)) / post_dt.GetXaxis()->GetXmax() - post_dt.GetXaxis()->GetXmin())) - 1); 
  }

  excess_dt.Write();
 
  TH1D excess_dlt("excess_dlt", "excess;dlt", nbins_log, dlt_bins.data());
  for (int bin_idx = 1; bin_idx < excess_dlt.GetNbinsX()+1; ++bin_idx){
    if (post_dlt.GetBinContent(bin_idx) != 0 ){
      excess_dlt.SetBinContent(bin_idx, (pre_dlt.GetBinContent(bin_idx) / post_dlt.GetBinContent(bin_idx)) - 1);
    } else {
      excess_dlt.SetBinContent(bin_idx, 0);
    }
  }
  excess_dlt.Write();

  TH1D excess_dll("excess_dll", "excess;dll", nbins, -5000, 5000);
  for (int bin_idx = 1; bin_idx < excess_dll.GetNbinsX()+1; ++bin_idx){
    if (post_dll.GetBinContent(bin_idx) != 0 ){
      excess_dll.SetBinContent(bin_idx, (pre_dll.GetBinContent(bin_idx) / post_dll.GetBinContent(bin_idx)) - 1); 
    } else {
      excess_dll.SetBinContent(bin_idx, 0);
    }

  }
  excess_dll.Write();
  
  TH1D excess_muqismsk("excess_muqismsk", "excess;muqismsk", nbins, 0, 250000);
  for (int bin_idx = 1; bin_idx < excess_muqismsk.GetNbinsX()+1; ++bin_idx){
    if (post_muqismsk.GetBinContent(bin_idx) != 0 ){
      excess_muqismsk.SetBinContent(bin_idx, (pre_muqismsk.GetBinContent(bin_idx) / post_muqismsk.GetBinContent(bin_idx)) - 1); 
    } else {
      excess_muqismsk.SetBinContent(bin_idx, 0);
    }
    
  }
  excess_muqismsk.Write();

  TH1D excess_resQ("excess_resQ", "excess;resQ", nbins, -100000, 100000);
  for (int bin_idx = 1; bin_idx < excess_resQ.GetNbinsX()+1; ++bin_idx){
    if (post_resQ.GetBinContent(bin_idx) != 0 ){
      excess_resQ.SetBinContent(bin_idx, (pre_resQ.GetBinContent(bin_idx) / post_resQ.GetBinContent(bin_idx)) - 1); 
    } else {
      excess_resQ.SetBinContent(bin_idx, 0);
    }
  }
  excess_resQ.Write();

  TH1D likepre(("likepre_"+name).c_str(), "likepre;L_{spall}", nbins_log_like, like_bins.data());
  TH1D likepost(("likepost_"+name).c_str(), "likepost;L_{spall}", nbins_log_like, like_bins.data());
  
  for (const auto& p : vp){

    const int dt_bin = excess_dt.FindBin(abs(p.dt));
    const int dlt_bin = excess_dlt.FindBin(abs(p.dlt));
    const int dll_bin = excess_dll.FindBin(abs(p.dll));
    const int muqismsk_bin = excess_muqismsk.FindBin(abs(p.muqismsk));
    const int resQ_bin = excess_resQ.FindBin(abs(p.resQ));
      
    //double likelihood = likeli_dt.GetBinContent(dt_bin);// *
    double likelihood = excess_dt.GetBinContent(dt_bin) *
      excess_dlt.GetBinContent(dlt_bin);// *
      //  excess_dll.GetBinContent(dll_bin) *
      //      excess_muqismsk.GetBinContent(muqismsk_bin) *
      //      excess_resQ.GetBinContent(resQ_bin);      
    
    likelihood = std::max(likelihood, 0.0);
    p.dt > 0 ? likepre.Fill(likelihood) : likepost.Fill(likelihood);
    
  }
  
  likepre.Write();
  likepost.Write();
  
}

double CalculateSpallationVariables::CalculateTrackLen(float* muon_entrypoint, float* muon_direction, double* exitpt){

  
  /*
    If you're reading this Marcus, Hello - you might be thinking that Matthew has copied this function from the RelicMuonPlots verbatim and you would indeed be correct. 
    You see, when I originally conceived of this tool I thought it wouldn't overlap quite so heavily with what you had done already - unforunately I was wrong. But nice thing for you, you won't need to check this when I push it because, hey ho, you wrote it didn't you. Feel free to bring up this sloppiness of mine in the form of a slack message but know that the only response you're going to get is: "that's terribly interesting, would you like a brew". 
  */
  
  // for reference HITKTK is the water volume height and DITKTK is its diameter,
  // these are #defined constants in geotnkC.h
	
  // sanity check muon is inward going or it's not going to go through the tank
  if( ((-muon_entrypoint[0]*muon_direction[0] + -muon_entrypoint[1]*muon_direction[1])<0) &&
      (muon_entrypoint[1]>=DITKTK/2.) ){
    // not radially inwards
    Log(m_unique_name+": Muon trajectory is not into tank!",v_error,m_verbose);
    return 0;
  }
  if( (muon_entrypoint[2] >= ( HITKTK/2.) && muon_direction[2]>0) ||
      (muon_entrypoint[2] <= (-HITKTK/2.) && muon_direction[2]<0) ){
    Log(m_unique_name+": Muon track points out of endcaps!",v_error,m_verbose);
    // pointing out of barrel
    return 0;
  }
	
  // calculate track length under assumption of a through-going muon,
  // based on its entry point and direction
  // first check for muons directed in the x-y plane
  if(std::abs(muon_direction[2]) > 0.1){
    double dist_to_endcap = (HITKTK/2.) - std::abs(muon_entrypoint[2]);
    if(std::signbit(muon_entrypoint[2]) != std::signbit(muon_direction[2])){
      dist_to_endcap = (HITKTK - dist_to_endcap);
    }
    // start with the simple case; project to the plane of the appropriate endcap
    double dxdz = muon_direction[0]/std::abs(muon_direction[2]);
    double dydz = muon_direction[1]/std::abs(muon_direction[2]);
    double proj_x = muon_entrypoint[0] + dxdz*dist_to_endcap;
    double proj_y = muon_entrypoint[1] + dydz*dist_to_endcap;
    if(((proj_x*proj_x)+(proj_y*proj_y)) < std::pow(DITKTK/2.,2.)){
      // if the projected point is within the tank radius, this is where it exits the tank
      double xtravel = proj_x-muon_entrypoint[0];
      double ytravel = proj_y-muon_entrypoint[1];
      double tracklen = std::sqrt(std::pow(xtravel,2)+
				  std::pow(ytravel,2.)+
				  std::pow(dist_to_endcap,2.));
			
      if(exitpt){
	exitpt[0] = proj_x;
	exitpt[1] = proj_y;
	exitpt[2] = muon_entrypoint[2] + dist_to_endcap*(muon_direction[2]>0 ? 1 : -1);
      }
      return tracklen;
    }
  }
  // if the muon is in the x-y plane, or the projected point is outside the tank,
  // then the track left through the barrel
	
  // get the angle of the track in the x-y plane
  double trackangle = std::atan2(muon_direction[1],muon_direction[0]);
	
  // get angle the track makes from entry point to the centre of the tank
  double entrypointpolarangle = std::atan2(-muon_entrypoint[1],-muon_entrypoint[0]);
	
  double chordlen;
  if(std::sqrt(std::pow(muon_entrypoint[0],2.)+std::pow(muon_entrypoint[1],2.))==DITKTK/2.){
    // with these we can calculate the angle subtended by the chord the muon track makes with:
    double angle_subtended = M_PI - 2.*(trackangle - entrypointpolarangle);
    // and from this we can work out the length of the chord:
    chordlen = std::abs(DITKTK * std::sin(angle_subtended/2.));
  } else {
    // ah, but if the muon entered through an endcap, the chord is truncated
    double a = std::sqrt(std::pow(muon_entrypoint[0],2.)+std::pow(muon_entrypoint[1],2.));
    // a is the distance from entry point to centre of endcap
    double C = entrypointpolarangle - trackangle;
    // C is the angle from the trajectory to tank centre
    // from a/SinA = c/SinC; -> sinA = a/c*SinC
    double A = std::asin(a/(DITKTK/2.) * std::sin(C));
    // A is the opening angle of the chord
    chordlen = a*std::cos(C) + (DITKTK/2.)*std::cos(A);
  }
	
  // amount of z travel is then chord length times z gradient
  double dzdr = muon_direction[2]/std::sqrt(std::pow(muon_direction[0],2.)+
					    std::pow(muon_direction[1],2.));
  double zdist = dzdr * chordlen;
  // sum is track length
  double tracklen = std::sqrt(std::pow(zdist,2.)+std::pow(chordlen,2.));
	
  if(exitpt){
    exitpt[0] = muon_entrypoint[0] + (chordlen * std::cos(trackangle));
    exitpt[1] = muon_entrypoint[1] + (chordlen * std::sin(trackangle));
    exitpt[2] = muon_entrypoint[2] + zdist;
  }
	
  return tracklen;
	
}

std::vector<double> CalculateSpallationVariables::MakeLogBins(double xmin, double xmax, int nbins){
        std::vector<double> binedges(nbins+1);
        double xxmin=log10(xmin);
        double xxmax = log10(xmax);
        for(int i=0; i<nbins; ++i){
                binedges[i] = pow(10,xxmin + (double(i)/(nbins-1.))*(xxmax-xxmin));
        }
        binedges[nbins+1] = xmax; // required                                                                                                                  
        return binedges;
}

void CalculateSpallationVariables::CreateOutputFile(){
  std::string data_output_file_str = "";
  m_variables.Get("data_output_file_str", data_output_file_str);
  if (data_output_file_str.empty()){
    throw std::runtime_error("CalculateSpallationVariables::CreateData_OutputFile: No data_output file specified!");
  }
  data_output_file_ptr = new TFile(data_output_file_str.c_str(), "RECREATE");
  data_output_tree_ptr = new TTree("spall", "spall");

  data_output_tree_ptr->Branch("pairing_string", &v_pairing_info);
  data_output_tree_ptr->Branch("dt", &v_dt);
  data_output_tree_ptr->Branch("dlt", &v_dlt);
  data_output_tree_ptr->Branch("dll", &v_dll);
  data_output_tree_ptr->Branch("muqismsk", &v_muqismsk);
  data_output_tree_ptr->Branch("resQ", &v_resQ);
  data_output_tree_ptr->Branch("muon_type", &v_muontype);
  data_output_tree_ptr->Branch("bsenergy", &v_bsenergy);
  data_output_tree_ptr->Branch("do_check", &v_do_check);
  
  return;
}


// change muon sided to relic sided, get rid of the execute hists, save pairing info and relic info, new toolchain for making likelihood hists, another one making cuts


