#include "CalculateNeutronCloudVertex.h"

#include "NeutronInfo.h"
#include "MTreeReader.h"

#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"

#include <algorithm>

CalculateNeutronCloudVertex::CalculateNeutronCloudVertex():Tool(){}

bool CalculateNeutronCloudVertex::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  GetTreeReader();

  CreateOutputFile();

  N_SLE_plot = TH1D("N_SLE_plot", "number of SLE triggers after muon", 20, 0, 20);
  mult_plot = TH1D("mult_plot", "multiplcity of neutron cloud;multiplcity", 20, 0, 20);
  dist_to_mu_plot = TH1D("dist_to_mu_plot", "distance to muon plot; distance [cm]", 100, 100, 100);
  
  return true;
}

bool CalculateNeutronCloudVertex::Execute(){

  std::vector<NeutronInfo> neutrons = {};
  m_data->CStore.Get("event_neutrons", neutrons);

  if (neutrons.empty()){
    mult = 0;
    std::fill(neutron_cloud_vertex.begin(), neutron_cloud_vertex.end(), 0);
    nvc_tree_ptr->Fill();
    return true;
  }
  
  for (const auto& neutron : neutrons){
    for (int dim = 0; dim < 3; ++dim){
      neutron_cloud_vertex.at(dim) += neutron.bs_vertex.at(dim) / neutrons.size();
    }
  }
  
  mult = neutrons.size();
  mult_plot.Fill(mult);

  dist_to_mu_plot.Fill(ClosestApproach(neutron_cloud_vertex));

  int N_SLE = 0;
  m_data->CStore.Get("N_SLE", N_SLE);
  N_SLE_plot.Fill(N_SLE);
  
  nvc_tree_ptr->Fill();

  if (MU_tree_reader->GetEntryNumber() % 1000 == 0){
    nvc_file_ptr->cd();
    nvc_tree_ptr->Write();
  }
  
  return true;
}

bool CalculateNeutronCloudVertex::Finalise(){

  std::string outfile_name = "";
  bool ok = m_variables.Get("outfile_name", outfile_name);
  if (!ok || outfile_name.empty()){ outfile_name = "calculateneutroncloudvertex_out.root";}

  TFile* outfile = TFile::Open(outfile_name.c_str(), "RECREATE");
  if (outfile == nullptr){
    throw std::runtime_error("CalculateNeutronCloudVertex::Finalise - Couldn't open output file");
  }

  outfile->cd();
  mult_plot.Write();
  N_SLE_plot.Write();
  dist_to_mu_plot.Write();

  nvc_file_ptr->cd();
  nvc_tree_ptr->Write();
  
  return true;
}

void CalculateNeutronCloudVertex::GetTreeReader(){
  std::string tree_reader_str = "";
  m_variables.Get("MU_TreeReader", tree_reader_str);
  if (m_data->Trees.count(tree_reader_str) == 0){
    throw std::runtime_error("CalculateNeutronCloudVertex::Execute - Failed to get treereader "+tree_reader_str+"!");
  }
  MU_tree_reader = m_data->Trees.at(tree_reader_str);
}


double CalculateNeutronCloudVertex::ClosestApproach(const std::vector<double>& vertex) {
  const std::vector<double> muon_ent = std::vector<double>(skroot_mu_.muentpoint, skroot_mu_.muentpoint + 3);
  muon_dir = std::vector<double>(skroot_mu_.mudir, skroot_mu_.mudir + 3);

  std::vector<double> diff = std::vector<double>(3, 0);
  for (int i = 0; i < 3; ++i){diff.at(i) = vertex.at(i) - muon_ent.at(i);}

  std::vector<double> proj = std::vector<double>(3, 0);
  std::vector<double> dist_vec = std::vector<double>(3, 0);

  const double diff_mudir_ip = std::inner_product(diff.begin(), diff.end(), muon_dir.begin(), 0);
  const double mudir_muent_ip = std::inner_product(muon_dir.begin(), muon_dir.end(), muon_dir.end(), 0);
  
  for (int i = 0; i < 3; ++i){
    proj.at(i) = (diff_mudir_ip / mudir_muent_ip) * muon_dir.at(i);
    dist_vec.at(i) = vertex.at(i) - proj.at(i) - muon_ent.at(i);
  }
  
  return sqrt(std::inner_product(dist_vec.begin(), dist_vec.begin(), dist_vec.end(), 0));
}

void CalculateNeutronCloudVertex::CreateOutputFile(){
  std::string nvc_file_str = "";
  m_variables.Get("nvc_file_str", nvc_file_str);
  if (nvc_file_str.empty()){
    throw std::runtime_error("CalculateNeutronCloudVertex::CreateOutputFile - no output file specified!");
  }
  nvc_file_ptr = TFile::Open(nvc_file_str.c_str(), "UPDATE");
  nvc_tree_ptr = new TTree("neutron_cloud_info", "neutron_cloud_info");

  nvc_tree_ptr->Branch("neutron_cloud_multiplicity", &mult);
  nvc_tree_ptr->Branch("neutron_cloud_vertex", &neutron_cloud_vertex);
  nvc_tree_ptr->Branch("muon_dir", &muon_dir);
  return;
}
