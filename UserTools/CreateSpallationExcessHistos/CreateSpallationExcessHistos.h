#ifndef CreateSpallationExcessHistos_H
#define CreateSpallationExcessHistos_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"
#include "TFile.h"
#include "TH1D.h"

class CreateSpallationExcessHistos: public Tool {

 public:

  CreateSpallationExcessHistos();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:

  MTreeReader* tree_ptr = nullptr;
  
  void GetReader();
  void GetBranchValues();
  void SetupHists();
  void MakeOutputFile();
  std::vector<double> MakeLogBins(const double&, const double&, const int&) const;
  
  TH1D pre_dt, post_dt;
  TH1D pre_dlt, post_dlt;
  TH1D pre_dll, post_dll;
  TH1D pre_muqismsk, post_muqismsk;
  TH1D pre_resQ, post_resQ;
  
  const std::vector<double>* v_dt = nullptr;
  const std::vector<double>* v_dlt = nullptr;
  const std::vector<double>* v_dll = nullptr;
  const std::vector<double>* v_muqismsk = nullptr;
  const std::vector<double>* v_resQ = nullptr;
  const std::vector<bool>* v_do_check = nullptr;

  std::vector<double> dt_bins = {};
  std::vector<double> dlt_bins = {};

  const int nbins_log = 25;
  const int nbins = 200;

  TFile* output_file_ptr = nullptr;
  
};

#endif
