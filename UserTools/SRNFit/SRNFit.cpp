#include "SRNFit.h"

#include "MTreeReader.h"
#include "TROOT.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "Fit/Fitter.h"
#include "Fit/UnBinData.h"
#include "Fit/FitConfig.h"
#include "Math/WrappedMultiTF1.h"


SRNFit::SRNFit():Tool(){}

bool SRNFit::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;

  GetReader();
  GetInputFiles();
  GetHists();
  CreateShapes();

  // Set up of output file:
  std::string output_file_str = "";
  m_variables.Get("output_file_str", output_file_str);
  if (output_file_str.empty()){
    throw std::runtime_error("SRNFit::Finalise: No output file specified!");
  }
  data_file_ptr = new TFile(output_file_str.c_str(), "RECREATE");

  return true;
}


bool SRNFit::Execute(){

  // Get LOWE->bsenergy
  bool ok = tree_ptr->Get("LOWE", lowe_ptr);
  if (!ok){throw std::runtime_error("SRNFit::Execute: Couldn't get LOWE!");}

  v_data.push_back(lowe_ptr->bsenergy);
  
  return true;
}


bool SRNFit::Finalise(){

  // data_graph = TGraph("data", "data");
  // for (int i = 0; i < v_data.size(); ++i){
  //   data_graph.AddPoint(

  data_plot = TH1D("data", "data", 100, 0, 100);
  int c = 0;
  for (const auto& datum : v_data){
    // std::cout << datum;
    // c % 7 != 0 ? std::cout << "," : std::cout << "," << std::endl;
    // ++c;
    data_plot.Fill(datum);
  }
  
  UnbinnedFit();

  // Save results of shape fits
  data_file_ptr->cd();
  
  srn_hist->Write();
  srn_shape.Write();
  li9_hist->Write();
  li9_shape.Write();
  reactor_hist->Write();
  reactor_shape.Write();
  ncqe_hist->Write();
  ncqe_shape.Write();
  nonncqe_hist->Write();
  nonncqe_shape.Write();

  srn_fit.Write();

  data_plot.Write();
  
  data_file_ptr->Close();
  
  return true;
}

void SRNFit::GetInputFiles(){
  std::string filename_str = "";
  if (!m_variables.Get("ibd_hist_file", filename_str)){
    throw std::runtime_error("SRNFit::GetInputFiles: Couldn't get ibd_hist input file!");
  }
  ibd_hist_file_ptr = TFile::Open(filename_str.c_str(), "READ");
  if (ibd_hist_file_ptr->IsZombie()){
    throw std::runtime_error("SRNFit::GetInputFiles: Couldn't open ibd_hist file!");
  }
  filename_str.clear();
  if (!m_variables.Get("atmos_hist_file", filename_str)){
    throw std::runtime_error("SRNFit::GetInputFiles: Couldn't get atmos_hist input file!");
  }
  atmos_hist_file_ptr = TFile::Open(filename_str.c_str(), "READ");
  if (atmos_hist_file_ptr->IsZombie()){
    throw std::runtime_error("SRNFit::GetInputFiles: Couldn't open atmos_hist file!");
  }

  return;
}

void SRNFit::GetHists(){
  srn_hist = (TH1D*) ibd_hist_file_ptr->Get("weight_nakazato1_IBD_bsenergy");
  if (srn_hist == nullptr){throw std::runtime_error("SRNFit::GetHists: couldn't get srn hist!");}
  li9_hist = (TH1D*) ibd_hist_file_ptr->Get("weight_li9_IBD_bsenergy");
  if (li9_hist == nullptr){throw std::runtime_error("SRNFit::GetHists: couldn't get li9 hist!");}
  reactor_hist = (TH1D*) ibd_hist_file_ptr->Get("weight_reactor_IBD_bsenergy");
  if (reactor_hist == nullptr){throw std::runtime_error("SRNFit::GetHists: couldn't get reactor hist!");}
  ncqe_hist = (TH1D*) atmos_hist_file_ptr->Get("bsenergy");
  if (ncqe_hist == nullptr){throw std::runtime_error("SRNFit::GetHists: couldn't get ncqe hist!");}
  nonncqe_hist = (TH1D*) atmos_hist_file_ptr->Get("nonNCQE_bsenergy");
  if (nonncqe_hist == nullptr){throw std::runtime_error("SRNFit::GetHists: couldn't get nonncqe hist!");}
  return;
}

void SRNFit::CreateShapes(){
  // normalise hists
  srn_hist->Scale(1/srn_hist->Integral());
  li9_hist->Scale(1/li9_hist->Integral());
  reactor_hist->Scale(1/reactor_hist->Integral());
  ncqe_hist->Scale(1/ncqe_hist->Integral());
  nonncqe_hist->Scale(1/nonncqe_hist->Integral());

  // fit 6th order pols to each hist, we'll tweak them when I can see them
  std::cout << "fitting srn shape"<<std::endl;
  srn_shape = TF1("srn", "pol6", bse_min, 65);
  srn_hist->Fit(&srn_shape, "R");
  std::cout << "fitting li9 shape"<<std::endl;
  li9_shape = TF1("li9", "pol6", bse_min, 18);
  li9_hist->Fit(&li9_shape, "R");
  std::cout << "fitting reactor shape"<<std::endl;
  reactor_shape = TF1("reactor", "expo", 8, 15);
  reactor_hist->Fit(&reactor_shape, "R");
  std::cout << "fitting ncqe shape"<<std::endl;
  ncqe_shape = TF1("ncqe", "pol6", bse_min, bse_max);
  ncqe_hist->Fit(&ncqe_shape, "R");
  std::cout << "fitting nonncqe shape"<<std::endl;
  nonncqe_shape = TF1("nonncqe", "pol6", bse_min, bse_max);
  nonncqe_hist->Fit(&nonncqe_shape, "R");
  
  return;
}

void SRNFit::UnbinnedFit(){

  // create TF1 for fit, TF1 it's your lucky day, you're getting created
  const int fit_n = 5;
  const int fit_dim = 1;
  srn_fit = TF1("srn_fit", this, &SRNFit::FitFunction, bse_min, bse_max, fit_n, "SRNFit", "FitFunction");

  // set parameter names, can't name them after the beatles, too silly, not enough beatles:
  srn_fit.SetParName(0, "SRN");
  srn_fit.SetParName(1, "li9");
  srn_fit.SetParName(2, "reactor");
  srn_fit.SetParName(3, "NCQE");
  srn_fit.SetParName(4, "nonNQCE");

  // set parameter limits: TODO

  // convert TF1 to use in a unbinned fit:
  ROOT::Math::WrappedMultiTF1 srn_fit_func(srn_fit, srn_fit.GetNdim());
  
  // build a fitter, this is pretty exciting isn't it?
  ROOT::Fit::Fitter fittest;
  fittest.SetFunction(srn_fit_func);

  // convert data to UnBinData object, there's nothing else better to do is there?
  ROOT::Fit::UnBinData data(v_data.size());
  for (auto datum : v_data){data.Add(datum);}

  // do the fit, will it work? will it work?
  fittest.LikelihoodFit(data);

  // print the results:
  // std::cout "YOU HAVE FOUND SRN YOU BIG LUNK!" // not just yet
  ROOT::Fit::FitResult r=fittest.Result();
  r.Print(std::cout);

  std::cout<<"binned chi2 fit parameters = {";
  for(int i=0; i<(srn_fit.GetNpar()); i++){
    std::cout<<srn_fit.GetParameter(i);
    (i<(srn_fit.GetNpar()-1)) ? std::cout<<", " : std::cout<<"};" << std::endl;
  }

  
  //draw the results
  
  
  return;
}

double SRNFit::FitFunction(double* x, double* p){
  return p[0] * srn_shape.Eval(x[0]) +
    p[1] * li9_shape.Eval(x[0]) +
    p[2] * reactor_shape.Eval(x[0]) +
    p[3] * ncqe_shape.Eval(x[0]) +
    p[4] * nonncqe_shape.Eval(x[0]);  
}

void SRNFit::GetReader(){
  std::string reader_name = "";
  if (!m_variables.Get("reader_name", reader_name) || reader_name.empty()){
    throw std::runtime_error("SRNFit::GetReader - no reader_name specified!");
  }
  if (m_data->Trees.count(reader_name) == 0){
    throw std::runtime_error("SRNFit::GetReader - reader not found!");
  }
  tree_ptr = m_data->Trees.at(reader_name);
  return;
}
