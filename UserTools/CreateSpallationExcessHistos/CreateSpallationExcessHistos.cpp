#include "CreateSpallationExcessHistos.h"

CreateSpallationExcessHistos::CreateSpallationExcessHistos():Tool(){}


bool CreateSpallationExcessHistos::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  GetReader();
  SetupHists();

  MakeOutputFile();
  
  return true;
}

bool CreateSpallationExcessHistos::Execute(){

  GetBranchValues();

  for (int muon_idx = 0; muon_idx < v_dt->size(); ++muon_idx){
    if (v_do_check->at(muon_idx) == false){continue;}
    if (v_dt->at(muon_idx) < 0){
      pre_dt.Fill(abs(v_dt->at(muon_idx)));
      pre_dlt.Fill(v_dlt->at(muon_idx));
      pre_dll.Fill(v_dll->at(muon_idx));
      pre_muqismsk.Fill(v_muqismsk->at(muon_idx));
      pre_resQ.Fill(v_resQ->at(muon_idx));		  
    } else {
      post_dt.Fill(v_dt->at(muon_idx));
      post_dlt.Fill(v_dlt->at(muon_idx));
      post_dll.Fill(v_dll->at(muon_idx));
      post_muqismsk.Fill(v_muqismsk->at(muon_idx));
      post_resQ.Fill(v_resQ->at(muon_idx));
    }
  } 

  return true;
}

bool CreateSpallationExcessHistos::Finalise(){

  output_file_ptr->cd();
  
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

  pre_dt.Write();
  post_dt.Write();
  
  output_file_ptr->Close();
  
  return true;
}

void CreateSpallationExcessHistos::GetReader(){
  std::string reader_name = "";
  if (!m_variables.Get("reader_name", reader_name) || reader_name.empty()){
    throw std::runtime_error("CalculateSpallationExcessHistos::GetReader - no reader_name specified!");
  }
  if (m_data->Trees.count(reader_name) == 0){
    throw std::runtime_error("CalculateSpallationExcessHistos::GetReader - reader not found!");
  }
  tree_ptr = m_data->Trees.at(reader_name);
  return;
}

void CreateSpallationExcessHistos::GetBranchValues(){
  bool ok = tree_ptr->Get("dt", v_dt);
  if (!ok){throw std::runtime_error("CreateSpallationExcessHistos::GetBranchValues: Couldn't get v_dt!");}
  ok = tree_ptr->Get("dlt", v_dlt);
  if (!ok){throw std::runtime_error("CreateSpallationExcessHistos::GetBranchValues: Couldn't get v_dlt!");}
  ok = tree_ptr->Get("dll", v_dll);
  if (!ok){throw std::runtime_error("CreateSpallationExcessHistos::GetBranchValues: Couldn't get v_dll!");}
  ok = tree_ptr->Get("muqismsk", v_muqismsk);
  if (!ok){throw std::runtime_error("CreateSpallationExcessHistos::GetBranchValues: Couldn't get v_muqismsk!");}
  ok = tree_ptr->Get("resQ", v_resQ);
  if (!ok){throw std::runtime_error("CreateSpallationExcessHistos::GetBranchValues: Couldn't get v_resQ!");}
  ok = tree_ptr->Get("do_check", v_do_check);
  if (!ok){throw std::runtime_error("CreateSpallationExcessHistos::GetBranchValues: Couldn't get v_do_check!");}
  return;
}

void CreateSpallationExcessHistos::SetupHists(){
  
  dt_bins = MakeLogBins(0.001, 60, nbins_log+1);
  dlt_bins = MakeLogBins(1, 5000, nbins_log+1);

  pre_dt = TH1D("pre_dt", "pre_dt", nbins_log, dt_bins.data()), post_dt = TH1D("post_dt", "post_dt", nbins_log, dt_bins.data());
  pre_dlt = TH1D("pre_dlt", "pre_dlt", nbins_log, dlt_bins.data()), post_dlt = TH1D("post_dlt", "post_dlt", nbins_log, dlt_bins.data());
  pre_dll = TH1D("pre_dll", "pre_dll", nbins, -5000, 5000), post_dll = TH1D("post_dll", "post_dll", nbins, -5000, 5000);
  pre_muqismsk = TH1D("pre_muqismsk", "pre_muqismsk", nbins, 0, 250000), post_muqismsk = TH1D("post_muqismsk", "post_muqismsk", nbins, 0, 250000);
  pre_resQ = TH1D("pre_resQ", "pre_resQ", nbins, -100000, 1000000), post_resQ = TH1D("post_resQ", "post_resQ", nbins, -100000, 1000000);

  return;
}

void CreateSpallationExcessHistos::MakeOutputFile(){

  std::string output_fname_str = "";
  m_variables.Get("output_fname", output_fname_str);
  if (output_fname_str.empty()){
    throw std::runtime_error("CreateSpallationExcessHistos::MakeOutputFile: no output file specified!");
  }
  output_file_ptr = TFile::Open(output_fname_str.c_str(), "RECREATE");
  if (output_file_ptr == nullptr || output_file_ptr->IsZombie()){
   throw std::runtime_error("CreateSpallationExcessHistos::MakeOutputFile: failed to open output file!"); 
  }
  
  return;
}

std::vector<double> CreateSpallationExcessHistos::MakeLogBins(const double& xmin, const double& xmax, const int& nbins) const {
        std::vector<double> binedges(nbins+1);
        double xxmin=log10(xmin);
        double xxmax = log10(xmax);
        for(int i=0; i<nbins; ++i){
                binedges[i] = pow(10,xxmin + (double(i)/(nbins-1.))*(xxmax-xxmin));
        }
        binedges[nbins] = xmax; // required    
        return binedges;
}


