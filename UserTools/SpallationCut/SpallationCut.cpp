#include "SpallationCut.h"

SpallationCut::SpallationCut():Tool(){}


bool SpallationCut::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  OpenSingleFile();
  GetReader();

  MakeOutputFile();
  
  return true;
}


bool SpallationCut::Execute(){
  
  GetBranchValues();
  //calculate likelihood

  v_likelihoods.clear();
  
  for (int i = 0; i < v_dt->size(); ++i){

    if (v_do_check->at(i) == false){continue;}
    
    const int dt_bin = excess_dt_ptr->FindBin(abs(v_dt->at(i)));
    const int dlt_bin = excess_dlt_ptr->FindBin(v_dlt->at(i));
    const int dll_bin = excess_dll_ptr->FindBin(v_dll->at(i));
    const int muqismsk_bin = excess_muqismsk_ptr->FindBin(v_muqismsk->at(i));
    const int resQ_bin = excess_resQ_ptr->FindBin(v_resQ->at(i));
  
    double likelihood_dt = std::max(excess_dt_ptr->GetBinContent(dt_bin), 0.0);
    double likelihood_dlt = std::max(excess_dlt_ptr->GetBinContent(dlt_bin), 0.0);
    double likelihood_dll = std::max(excess_dll_ptr->GetBinContent(dll_bin), 0.0);
    double likelihood_muqismsk = std::max(excess_muqismsk_ptr->GetBinContent(muqismsk_bin), 0.0);
    double likelihood_resQ = std::max(excess_resQ_ptr->GetBinContent(resQ_bin), 0.0);
    
    v_likelihoods.push_back({likelihood_dt, likelihood_dlt, likelihood_dll, likelihood_muqismsk, likelihood_resQ});
    
  }

  output_tree_ptr->Fill();
  return true;
}


bool SpallationCut::Finalise(){

  output_file_ptr->cd();
  output_tree_ptr->Write();
  output_file_ptr->Close();
  
  return true;
}

void SpallationCut::OpenSingleFile(){
  std::string fname = "";
  if (!m_variables.Get("single_file", fname) || fname.empty()){
    throw std::runtime_error("SpallationCut::OpenSingleFile: no histogram file!");
  }
  TFile* f_ptr = TFile::Open(fname.c_str(), "OPEN");
  if (f_ptr->IsZombie()){
    throw std::runtime_error("SpallationCut::OpenSingleFile: couldn't open file!");
  }

  excess_dt_ptr = (TH1D*) f_ptr->Get("excess_dt");
  if (!excess_dt_ptr){throw std::runtime_error("SpallationCut::OpenSingleFile: couldn't retrieve excess_dt!");}
  excess_dlt_ptr = (TH1D*) f_ptr->Get("excess_dlt");
  if (!excess_dlt_ptr){throw std::runtime_error("SpallationCut::OpenSingleFile: couldn't retrieve excess_dlt!");}
  excess_dll_ptr = (TH1D*) f_ptr->Get("excess_dll");
  if (!excess_dll_ptr){throw std::runtime_error("SpallationCut::OpenSingleFile: couldn't retrieve excess_dll!");}
  excess_muqismsk_ptr = (TH1D*) f_ptr->Get("excess_muqismsk");
  if (!excess_muqismsk_ptr){throw std::runtime_error("SpallationCut::OpenSingleFile: couldn't retrieve excess_muqismsk!");}
  excess_resQ_ptr = (TH1D*) f_ptr->Get("excess_resQ");
  if (!excess_resQ_ptr){throw std::runtime_error("SpallationCut::OpenSingleFile: couldn't retrieve excess_resQ!");}

  f_ptr->Close();
  
  return;
}

void SpallationCut::GetReader(){
  std::string reader_name = "";
  if (!m_variables.Get("reader_name", reader_name) || reader_name.empty()){
    throw std::runtime_error("SpallationCut::OpenSingleFile: no reader name!");
  }
  if (m_data->Trees.count(reader_name) == 0){
    throw std::runtime_error("SpallationCut::OpenSingleFile: reader not found!");
  }
  relic_tree_ptr = m_data->Trees.at(reader_name);
  return ;
}

void SpallationCut::GetBranchValues(){
 bool ok = relic_tree_ptr->Get("dt", v_dt);
  if (!ok){throw std::runtime_error("CreateSpallationExcessHistos::GetBranchValues: Couldn't get v_dt!");}
  ok = relic_tree_ptr->Get("dlt", v_dlt);
  if (!ok){throw std::runtime_error("CreateSpallationExcessHistos::GetBranchValues: Couldn't get v_dlt!");}
  ok = relic_tree_ptr->Get("dll", v_dll);
  if (!ok){throw std::runtime_error("CreateSpallationExcessHistos::GetBranchValues: Couldn't get v_dll!");}
  ok = relic_tree_ptr->Get("muqismsk", v_muqismsk);
  if (!ok){throw std::runtime_error("CreateSpallationExcessHistos::GetBranchValues: Couldn't get v_muqismsk!");}
  ok = relic_tree_ptr->Get("resQ", v_resQ);
  if (!ok){throw std::runtime_error("CreateSpallationExcessHistos::GetBranchValues: Couldn't get v_resQ!");}
  ok = relic_tree_ptr->Get("do_check", v_do_check);
  if (!ok){throw std::runtime_error("CreateSpallationExcessHistos::GetBranchValues: Couldn't get v_do_check!");}  
  return;
}

void SpallationCut::MakeOutputFile(){
  std::string output_fname_str = "";
  m_variables.Get("output_fname", output_fname_str);
  if (output_fname_str.empty()){
    throw std::runtime_error("CreateSpallationExcessHistos::MakeOutputFile: no output file specified!");
  }
  output_file_ptr = TFile::Open(output_fname_str.c_str(), "RECREATE");
  if (output_file_ptr == nullptr || output_file_ptr->IsZombie()){
   throw std::runtime_error("CreateSpallationExcessHistos::MakeOutputFile: failed to open output file!");
  }

  output_tree_ptr = new TTree("likelihoods", "likelihoods");

  output_tree_ptr->Branch("likelihoods", &v_likelihoods); 
  
}
