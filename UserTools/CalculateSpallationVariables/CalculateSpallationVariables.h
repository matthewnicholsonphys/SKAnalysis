#ifndef CalculateSpallationVariables_H
#define CalculateSpallationVariables_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"

#include "TH1D.h"

struct PairingInfo {
  float dt = 0;
  float dlt = 0;
  float dll = 0;
  float muqismsk = 0;
  float resQ = 0;
  int muon_type = 0;
  float bse = 0;
};

class CalculateSpallationVariables: public Tool {

public:

  CalculateSpallationVariables();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

private:

  void GetMuonBranchValues();
  void GetRelicBranchValues();
  std::string GetPairingString(PairingInfo) const;
  double CalculateTrackLen(float*, float*, double* exit=nullptr);
  void CreateLikelihood(const std::string&, const std::vector<PairingInfo>&) const;

  std::vector<double> MakeLogBins(double xmin, double xmax, int nbins);
  void CreateOutputFile();
  
  std::string run_type_str = "";
  
  const std::vector<float>* MatchedTimeDiff_ptr = nullptr;
  const std::vector<int>* MatchedOutEntryNums_ptr = nullptr;

  const LoweInfo* LOWE_ptr = nullptr;
  const MuInfo* MU_ptr = nullptr;
  
  MTreeReader* muon_tree_ptr = nullptr;
  MTreeReader* relic_tree_ptr = nullptr;

  int nbins = 200;
  int nbins_log = 25;
  int nbins_log_like = 100;
  std::vector<double> dt_bins = {};
  std::vector<double> dlt_bins = {};
  std::vector<double> like_bins = {};

  std::vector<std::string> v_pairing_info = {};
  std::vector<double> v_dt = {};
  std::vector<double> v_dlt = {};
  std::vector<double> v_dll = {};
  std::vector<double> v_muqismsk = {};
  std::vector<double> v_resQ = {};
  std::vector<double> v_muontype = {};
  std::vector<double> v_bsenergy = {};
  std::vector<bool> v_do_check = {};
    
  void GetReaders();

  std::map<std::string, std::vector<PairingInfo>> pairings = {};

  TH1D pre_dt_hist;
  TH1D post_dt_hist;
  TH1D pre_dlt_hist = TH1D("pre_dlt_hist", "pre_dlt_hist", nbins, 0, 5000);
  TH1D post_dlt_hist = TH1D("post_dlt_hist", "post_dlt_hist", nbins, 0, 5000);
  TH1D pre_dll_hist = TH1D("pre_dll_hist", "pre_dll_hist", nbins, -5000, 5000);
  TH1D post_dll_hist = TH1D("post_dll_hist", "post_dll_hist", nbins, -5000, 5000);
  TH1D pre_muqismsk_hist = TH1D("pre_muqismsk_hist", "pre_muqismsk_hist", nbins, 0, 250000);
  TH1D post_muqismsk_hist = TH1D("post_muqismsk_hist", "post_muqismsk_hist", nbins, 0, 250000);
  TH1D pre_resQ_hist = TH1D("pre_resQ_hist", "pre_resQ_hist", nbins, -100000, 1000000);
  TH1D post_resQ_hist = TH1D("post_resQ_hist", "post_resQ_hist", nbins, -100000, 1000000);

  TFile* data_output_file_ptr = nullptr;
  TTree* data_output_tree_ptr = nullptr;


  
};


#endif
