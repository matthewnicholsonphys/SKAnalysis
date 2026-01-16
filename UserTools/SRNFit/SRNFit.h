#ifndef SRNFit_H
#define SRNFit_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"

#include "MTreeReader.h"

class SRNFit: public Tool {

public:

  SRNFit();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

private:

  void GetHists();
  void GetReader();
  void GetInputFiles();
  void CreateShapes();
  void UnbinnedFit();
  double FitFunction(double*, double*);

  MTreeReader* tree_ptr = nullptr;
  const LoweInfo* lowe_ptr = nullptr;

  std::vector<double> v_data = {};

  double bse_min = 8;
  double bse_max = 100;

  TH1D* srn_hist = nullptr;
  TH1D* li9_hist = nullptr;
  TH1D* reactor_hist = nullptr;
  TH1D* ncqe_hist = nullptr;
  TH1D* nonncqe_hist = nullptr;
  
  TF1 srn_shape;
  TF1 li9_shape;
  TF1 reactor_shape;
  TF1 ncqe_shape;
  TF1 nonncqe_shape;

  TF1 srn_fit;
  TH1D data_plot;

  TFile* ibd_hist_file_ptr = nullptr;
  TFile* atmos_hist_file_ptr = nullptr;
  TFile* data_file_ptr = nullptr;
};


#endif
