#ifndef SpallationCut_H
#define SpallationCut_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "MTreeReader.h"

class SpallationCut: public Tool {

 public:

  SpallationCut();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

private:

  void OpenSingleFile();
  void GetReader();
  void MakeOutputFile();
  void GetBranchValues();

  TH1D* excess_dt_ptr = nullptr;
  TH1D* excess_dlt_ptr = nullptr;
  TH1D* excess_dll_ptr = nullptr;
  TH1D* excess_muqismsk_ptr = nullptr;
  TH1D* excess_resQ_ptr = nullptr;

  std::vector<double>* v_dt = nullptr;
  std::vector<double>* v_dlt = nullptr;
  std::vector<double>* v_dll = nullptr;
  std::vector<double>* v_muqismsk = nullptr;
  std::vector<double>* v_resQ = nullptr;
  std::vector<bool>* v_do_check = nullptr;
  
  MTreeReader* relic_tree_ptr = nullptr;

  TFile* output_file_ptr = nullptr;
  TTree* output_tree_ptr = nullptr;

  std::vector<std::array<double, 5>> v_likelihoods = {};
  
};

#endif
