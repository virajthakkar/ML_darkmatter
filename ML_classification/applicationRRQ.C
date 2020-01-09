/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides a simple example on how to use the trained classifiers
/// within an analysis module
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Exectuable: applicationRRQ
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker

 //   root -l ./applicationRRQ.C\(\"BDT\"\)


#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

void applicationRRQ( TString myMethodList = "" )
{

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // Cut optimisation
   Use["Cuts"]            = 1;
   Use["CutsD"]           = 1;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   //
   // 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 1;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 1;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 1;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 1; // k-nearest neighbour method
   //
   // Linear Discriminant Analysis
   Use["LD"]              = 1; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // Function Discriminant analysis
   Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   Use["DNN_CPU"] = 0;         // CUDA-accelerated DNN training.
   Use["DNN_GPU"] = 0;         // Multi-core accelerated DNN.
   //
   // Support Vector Machine
   Use["SVM"]             = 1;
   //
   // Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
   //
   // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 1;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start applicationRRQ" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "Color:!Silent" );

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   /*Float_t var1, var2;
   Float_t var3, var4;
   reader->AddVariable( "myvar1 := var1+var2", &var1 );
   reader->AddVariable( "myvar2 := var1-var2", &var2 );
   reader->AddVariable( "var3",                &var3 );
   reader->AddVariable( "var4",                &var4 );*/
   
  float precoiltNF,ytNF, qzpartOF,qrpart1OF,pzpartOF,prpartOF,qsummaxOF;
float PTNFchisq,ptNF,PTOFchisq,PTglitch1OFchisq,PTlfnoise1OFchisq,
QS1OFchisq,qsum1OF,QS2OFchisq,qsum2OF,qo1OF,qo2OF,qi1OF,qi2OF,qimean;

  
  reader->AddVariable( "precoiltNF",&precoiltNF);
  
   reader->AddVariable( "ytNF", &ytNF );
  reader->AddVariable( "qzpartOF", &qzpartOF  );
  reader->AddVariable( "qrpart1OF", &qrpart1OF );
  reader->AddVariable( "pzpartOF", &pzpartOF ); 
  reader->AddVariable( "prpartOF", &prpartOF );
    reader->AddVariable( "qsummaxOF", &qsummaxOF );
   
   //reader->AddVariable( "DCAdiff:=abs(dcad1-dcad2)", &dcad1 );
 
   // Spectator variables declared in the training have to be added to the reader, too
 /*  Float_t spec1,spec2;
   reader->AddSpectator( "spec1 := var1*2",   &spec1 );
   reader->AddSpectator( "spec2 := var1*3",   &spec2 );  */
   
   
   
   double cutvalue=-0.0528; //like
   cout<<"THE CUT VALUE APPLIED IS========================================="<<cutvalue<<endl;
   

  /* Float_t Category_cat1, Category_cat2, Category_cat3;
   if (Use["Category"]){
      // Add artificial spectators for distinguishing categories
      reader->AddSpectator( "Category_cat1 := var3<=0",             &Category_cat1 );
      reader->AddSpectator( "Category_cat2 := (var3>0)&&(var4<0)",  &Category_cat2 );
      reader->AddSpectator( "Category_cat3 := (var3>0)&&(var4>=0)", &Category_cat3 );
   }
   */

   // Book the MVA methods
//  Path : /home/viraj/Programs/root/tutorials/tmva/dataset/weights/Kstarweights/4featuresbdt    //weight file name:   Kstarrun_BDT.weights.xml
   //TString dir    = "/home/viraj/Programs/root/tutorials/tmva/";
   TString dir    = "dataset/weights/";  //"/home/viraj/Programs/root/tutorials/tmva/";
   
   TString prefix = "RRQtrain";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile );
      }
   }

   // Book output histograms
   UInt_t nbin = 100;
   TH1F *histLk(0);
   TH1F *histLkD(0);
   TH1F *histLkPCA(0);
   TH1F *histLkKDE(0);
   TH1F *histLkMIX(0);
   TH1F *histPD(0);
   TH1F *histPDD(0);
   TH1F *histPDPCA(0);
   TH1F *histPDEFoam(0);
   TH1F *histPDEFoamErr(0);
   TH1F *histPDEFoamSig(0);
   TH1F *histKNN(0);
   TH1F *histHm(0);
   TH1F *histFi(0);
   TH1F *histFiG(0);
   TH1F *histFiB(0);
   TH1F *histLD(0);
   TH1F *histNn(0);
   TH1F *histNnbfgs(0);
   TH1F *histNnbnn(0);
   TH1F *histNnC(0);
   TH1F *histNnT(0);
   TH1F *histBdt(0);
   TH1F *histBdtG(0);
   TH1F *histBdtB(0);
   TH1F *histBdtD(0);
   TH1F *histBdtF(0);
   TH1F *histRf(0);
   TH1F *histSVMG(0);
   TH1F *histSVMP(0);
   TH1F *histSVML(0);
   TH1F *histFDAMT(0);
   TH1F *histFDAGA(0);
   TH1F *histCat(0);
   TH1F *histPBdt(0);
   TH1F *histDnnGpu(0);
   TH1F *histDnnCpu(0);
   
   TCanvas *can;
   // ************************Mass histogram
   
   
double precoiltNFrange_1=0;
double precoiltNFrange_2=170; 
int precoiltNFbinsize=450; 
int qbinsize=500; // Prev: 500, 160
double q_y_range_1=0; //Prev: -10, -2
double q_y_range_2=160; //Prev: 80, 14

   TH2D *hERNR = new TH2D("hERNR","ER NR", precoiltNFbinsize, precoiltNFrange_1, precoiltNFrange_2, qbinsize, q_y_range_1, q_y_range_2);
hERNR->SetXTitle("precoiltNF [keV]");
hERNR->SetYTitle("qimean");
hERNR->GetYaxis()->SetTitleOffset(1.5);
hERNR->GetYaxis()->SetLabelSize(0.03);
hERNR->GetYaxis()->SetTitleSize(0.03);
//hERNR->SetLineColor(kGray+3);
//hERNR->SetFillColor(kGray+3);
//hERNR->SetMarkerColor(kGray+3);
   


   TH2D *hERNR_ML = new TH2D("hERNR_ML","ER NR ML-classification", precoiltNFbinsize, precoiltNFrange_1, precoiltNFrange_2, qbinsize, q_y_range_1, q_y_range_2);
hERNR_ML->SetXTitle("precoiltNF [keV]");
hERNR_ML->SetYTitle("qimean");
hERNR_ML->GetYaxis()->SetTitleOffset(1.5);
hERNR_ML->GetYaxis()->SetLabelSize(0.03);
hERNR_ML->GetYaxis()->SetTitleSize(0.03);
//hERNR_ML->SetLineColor(kGray+3);
//hERNR_ML->SetFillColor(kGray+3);
//hERNR_ML->SetMarkerColor(kGray+3);
   

TH2D *hyieldE = new TH2D("hyieldE","Yield vs Recoil Energy", precoiltNFbinsize, precoiltNFrange_1, precoiltNFrange_2, 100, 0, 2);
hyieldE->SetXTitle("precoiltNF [keV]");
hyieldE->SetYTitle("Yield");
hyieldE->GetYaxis()->SetTitleOffset(1.5);
hyieldE->GetYaxis()->SetLabelSize(0.03);
hyieldE->GetYaxis()->SetTitleSize(0.03);


TH2D *hyieldE_ML = new TH2D("hyieldE_ML","Yield vs Recoil Energy", precoiltNFbinsize, precoiltNFrange_1, precoiltNFrange_2, 100, 0, 2);
hyieldE_ML->SetXTitle("precoiltNF [keV]");
hyieldE_ML->SetYTitle("Yield");
hyieldE_ML->GetYaxis()->SetTitleOffset(1.5);
hyieldE_ML->GetYaxis()->SetLabelSize(0.03);
hyieldE_ML->GetYaxis()->SetTitleSize(0.03);


   
   
   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["DNN_GPU"]) histDnnGpu = new TH1F("MVA_DNN_GPU", "MVA_DNN_GPU", nbin, -0.1, 1.1);
   if (Use["DNN_CPU"]) histDnnCpu = new TH1F("MVA_DNN_CPU", "MVA_DNN_CPU", nbin, -0.1, 1.1);
   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
   if (Use["BDTB"])          histBdtB    = new TH1F( "MVA_BDTB",          "MVA_BDTB",          nbin, -1.0, 1.0 );
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTF"])          histBdtF    = new TH1F( "MVA_BDTF",          "MVA_BDTF",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TFile *input(0);
   
   //  **********************************************************CHANGE INPUT FILE NAME********************************************************************
   TString fname = "./RRQ_WS.root";
  // TString fname = "./RRQ_NR_cf.root";
  
   //TString fname = "./likesignplusplusfeatures.root";
 // TString fname = "./likesignplusplusfeatures_0.8pt1.2.root";
 
   
   if (!gSystem->AccessPathName( fname )) {
      input = TFile::Open( fname ); // check if file in local directory exists
   }
   else {
      TFile::SetCacheFileDir(".");
      input = TFile::Open("http://root.cern.ch/files/tmva_class_example.root", "CACHEREAD"); // if not: download from ROOT server
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

   // Event loop

   // Prepare the event tree
   // - Here the variable names have to corresponds to your tree
   // - You can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input->Get("tree_bkg");
  // TTree* theTree = (TTree*)input->Get("tree_bkg");
  // Float_t userVar1, userVar2;

  /* theTree->SetBranchAddress( "var1", &userVar1 );
   theTree->SetBranchAddress( "var2", &userVar2 );
   theTree->SetBranchAddress( "var3", &var3 );
   theTree->SetBranchAddress( "var4", &var4 );*/

   //  Int_t Event,mult_s;
 
 theTree->SetBranchAddress("PTNFchisq",&PTNFchisq);
theTree->SetBranchAddress("ptNF",&ptNF);
theTree->SetBranchAddress("PTOFchisq",&PTOFchisq);
theTree->SetBranchAddress("PTglitch1OFchisq",&PTglitch1OFchisq);
theTree->SetBranchAddress("PTlfnoise1OFchisq",&PTlfnoise1OFchisq);
theTree->SetBranchAddress("QS1OFchisq",&QS1OFchisq);
theTree->SetBranchAddress("qsum1OF",&qsum1OF);
theTree->SetBranchAddress("QS2OFchisq",&QS2OFchisq);
theTree->SetBranchAddress("qsum2OF",&qsum2OF);
theTree->SetBranchAddress("qo1OF",&qo1OF);
theTree->SetBranchAddress("qo2OF",&qo2OF);
theTree->SetBranchAddress("qi1OF",&qi1OF);
theTree->SetBranchAddress("qi2OF",&qi2OF);
theTree->SetBranchAddress("qimean",&qimean);


theTree->SetBranchAddress("precoiltNF",&precoiltNF);
theTree->SetBranchAddress("ytNF", &ytNF);
theTree->SetBranchAddress("qzpartOF", &qzpartOF);
theTree->SetBranchAddress("qrpart1OF", &qrpart1OF);
theTree->SetBranchAddress("pzpartOF", &pzpartOF);
theTree->SetBranchAddress("prpartOF", &prpartOF);
theTree->SetBranchAddress("qsummaxOF", &qsummaxOF);

  // theTree->SetBranchAddress( "abs(dcad1-dcad2)", &userVar1 );
  cout<<"The number of entries in theTree======================================"<<theTree->GetEntries()<<endl; 
  
  // Efficiency calculator for cut method
  Int_t    nSelCutsGA = 0;
  Double_t effS       = 0.7;
  
  std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests
  
  std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
  TStopwatch sw;
  sw.Start();
  
  //*****************************************************************EVENT LOOP ***********************************************************************************
   
  
   
  
 float precoiltNF_ER,ytNF_ER, qzpartOF_ER,qrpart1OF_ER,pzpartOF_ER,prpartOF_ER,qsummaxOF_ER;
float PTNFchisq_ER,ptNF_ER,PTOFchisq_ER,PTglitch1OFchisq_ER,PTlfnoise1OFchisq_ER,
QS1OFchisq_ER,qsum1OF_ER,QS2OFchisq_ER,qsum2OF_ER,qo1OF_ER,qo2OF_ER,qi1OF_ER,qi2OF_ER,qimean_ER;


 //TFile *fout = new TFile("A.root","RECREATE");
 TFile *fout = new TFile("A_WS.root","RECREATE");
  TTree *tree_ER = new TTree("tree_ER","tree_ER");
  tree_ER->Branch("precoiltNF",&precoiltNF_ER,"precoiltNF_ER/F");
  tree_ER->Branch("ytNF",&ytNF_ER,"ytNF_ER/F");
  tree_ER->Branch("qzpartOF",&qzpartOF_ER,"qzpartOF_ER/F");
  tree_ER->Branch("qrpart1OF",&qrpart1OF_ER,"qrpart1OF_ER/F");
  tree_ER->Branch("pzpartOF",&pzpartOF_ER,"pzpartOF_ER/F");
  tree_ER->Branch("prpartOF",&prpartOF_ER,"prpartOF_ER/F");
  tree_ER->Branch("qsummaxOF",&qsummaxOF_ER,"qsummaxOF_ER/F");
  
  tree_ER->Branch("PTNFchisq",&PTNFchisq_ER,"PTNFchisq_ER/F");
  tree_ER->Branch("ptNF",&ptNF_ER,"ptNF_ER/F");
  tree_ER->Branch("PTOFchisq",&PTOFchisq_ER,"PTOFchisq_ER/F");
  tree_ER->Branch("PTglitch1OFchisq",&PTglitch1OFchisq_ER,"PTglitch1OFchisq_ER/F");
  tree_ER->Branch("PTlfnoise1OFchisq",&PTlfnoise1OFchisq_ER,"PTlfnoise1OFchisq_ER/F");
  tree_ER->Branch("QS1OFchisq",&QS1OFchisq_ER,"QS1OFchisq_ER/F");
  tree_ER->Branch("qsum1OF",&qsum1OF_ER,"qsum1OF_ER/F");
  tree_ER->Branch("QS2OFchisq",&QS2OFchisq_ER,"QS2OFchisq_ER/F");
  tree_ER->Branch("qsum2OF",&qsum2OF_ER,"qsum2OF_ER/F");
  tree_ER->Branch("qo1OF",&qo1OF_ER,"qo1OF_ER/F");
  tree_ER->Branch("qo2OF",&qo2OF_ER,"qo2OF_ER/F");
  tree_ER->Branch("qi1OF",&qi1OF_ER,"qi1OF_ER/F");
  tree_ER->Branch("qi2OF",&qi2OF_ER,"qi2OF_ER/F");
  tree_ER->Branch("qimean",&qimean_ER,"qimean_ER/F");
  
  float precoiltNF_NR,ytNF_NR, qzpartOF_NR,qrpart1OF_NR,pzpartOF_NR,prpartOF_NR,qsummaxOF_NR;
float PTNFchisq_NR,ptNF_NR,PTOFchisq_NR,PTglitch1OFchisq_NR,PTlfnoise1OFchisq_NR,
QS1OFchisq_NR,qsum1OF_NR,QS2OFchisq_NR,qsum2OF_NR,qo1OF_NR,qo2OF_NR,qi1OF_NR,qi2OF_NR,qimean_NR;


  TTree *tree_NR = new TTree("tree_NR","tree_NR");
  tree_NR->Branch("precoiltNF",&precoiltNF_NR,"precoiltNF_NR/F");
  tree_NR->Branch("ytNF",&ytNF_NR,"ytNF_NR/F");
  tree_NR->Branch("qzpartOF",&qzpartOF_NR,"qzpartOF_NR/F");
  tree_NR->Branch("qrpart1OF",&qrpart1OF_NR,"qrpart1OF_NR/F");
  tree_NR->Branch("pzpartOF",&pzpartOF_NR,"pzpartOF_NR/F");
  tree_NR->Branch("prpartOF",&prpartOF_NR,"prpartOF_NR/F");
  tree_NR->Branch("qsummaxOF",&qsummaxOF_NR,"qsummaxOF_NR/F");
  
  tree_NR->Branch("PTNFchisq",&PTNFchisq_NR,"PTNFchisq_NR/F");
  tree_NR->Branch("ptNF",&ptNF_NR,"ptNF_NR/F");
  tree_NR->Branch("PTOFchisq",&PTOFchisq_NR,"PTOFchisq_NR/F");
  tree_NR->Branch("PTglitch1OFchisq",&PTglitch1OFchisq_NR,"PTglitch1OFchisq_NR/F");
  tree_NR->Branch("PTlfnoise1OFchisq",&PTlfnoise1OFchisq_NR,"PTlfnoise1OFchisq_NR/F");
  tree_NR->Branch("QS1OFchisq",&QS1OFchisq_NR,"QS1OFchisq_NR/F");
  tree_NR->Branch("qsum1OF",&qsum1OF_NR,"qsum1OF_NR/F");
  tree_NR->Branch("QS2OFchisq",&QS2OFchisq_NR,"QS2OFchisq_NR/F");
  tree_NR->Branch("qsum2OF",&qsum2OF_NR,"qsum2OF_NR/F");
  tree_NR->Branch("qo1OF",&qo1OF_NR,"qo1OF_NR/F");
  tree_NR->Branch("qo2OF",&qo2OF_NR,"qo2OF_NR/F");
  tree_NR->Branch("qi1OF",&qi1OF_NR,"qi1OF_NR/F");
  tree_NR->Branch("qi2OF",&qi2OF_NR,"qi2OF_NR/F");
  tree_NR->Branch("qimean",&qimean_NR,"qimean_NR/F");
  

  
  
  
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) 
   {
     // cout<<ievt<<endl;
     if (ievt%50000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     theTree->GetEntry(ievt);
  
	 //	 theTree->GetEntry(k1);
	 //var1 = userVar1 + userVar2;
	 //var2 = userVar1 - userVar2;
	 
      // Return the MVA outputs and fill into histograms

      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
      if (Use["DNN_GPU"]) histDnnGpu->Fill(reader->EvaluateMVA("DNN_GPU method"));
      if (Use["DNN_CPU"]) histDnnCpu->Fill(reader->EvaluateMVA("DNN_CPU method"));
      
      
      if (Use["BDT"          ])    // *************************************************BDT**************************************************************
      {
          if(qimean<1) continue;
          
          histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
          if(reader->EvaluateMVA( "BDT method")>cutvalue)
	    {
	      hERNR_ML->Fill(precoiltNF,qimean);
          hyieldE_ML->Fill(precoiltNF,ytNF);
          
     precoiltNF_NR=precoiltNF;
       ytNF_NR=ytNF;
       qzpartOF_NR=qzpartOF;
       qrpart1OF_NR=qrpart1OF;
       pzpartOF_NR=pzpartOF;
       prpartOF_NR=prpartOF;
       qsummaxOF_NR=qsummaxOF;
       
       PTNFchisq_NR=PTNFchisq;
       ptNF_NR=ptNF;
       PTOFchisq_NR=PTOFchisq;
       PTglitch1OFchisq_NR=PTglitch1OFchisq;
       PTlfnoise1OFchisq_NR=PTlfnoise1OFchisq;
       QS1OFchisq_NR=QS1OFchisq;
       qsum1OF_NR=qsum1OF;
       QS2OFchisq_NR=QS2OFchisq;
       qsum2OF_NR=qsum2OF;
       qo1OF_NR=qo1OF;
       qo2OF_NR=qo2OF;
       qi1OF_NR=qi1OF;
       qi2OF_NR=qi2OF;
       qimean_NR=0.5*(qi1OF+qi2OF);

           tree_NR->Fill();
          
	      
	    }
	  //  hERNR->Fill(precoiltNF,qimean);
	     if(reader->EvaluateMVA( "BDT method")<cutvalue)
	    {
         hERNR->Fill(precoiltNF,qimean);
         hyieldE->Fill(precoiltNF,ytNF);
         
         
              precoiltNF_ER=precoiltNF;
       ytNF_ER=ytNF;
       qzpartOF_ER=qzpartOF;
       qrpart1OF_ER=qrpart1OF;
       pzpartOF_ER=pzpartOF;
       prpartOF_ER=prpartOF;
       qsummaxOF_ER=qsummaxOF;
       
       PTNFchisq_ER=PTNFchisq;
       ptNF_ER=ptNF;
       PTOFchisq_ER=PTOFchisq;
       PTglitch1OFchisq_ER=PTglitch1OFchisq;
       PTlfnoise1OFchisq_ER=PTlfnoise1OFchisq;
       QS1OFchisq_ER=QS1OFchisq;
       qsum1OF_ER=qsum1OF;
       QS2OFchisq_ER=QS2OFchisq;
       qsum2OF_ER=qsum2OF;
       qo1OF_ER=qo1OF;
       qo2OF_ER=qo2OF;
       qi1OF_ER=qi1OF;
       qi2OF_ER=qi2OF;
       qimean_ER=0.5*(qi1OF+qi2OF);
       
        tree_ER->Fill();

         
         
        }
         //  cout<< reader->EvaluateMVA( "BDT method"           )<<"  -------------------------------------"<<invmass<<endl;
          
      }
      
      
      
      
      
      
      
      
      
      if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) );
      if (Use["BDTB"         ])   histBdtB   ->Fill( reader->EvaluateMVA( "BDTB method"          ) );
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTF"         ])   histBdtF   ->Fill( reader->EvaluateMVA( "BDTF method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

      
      //cout<< reader->EvaluateMVA( "BDT method"           )<<"  -------------------------------------"<<costheta<<endl;
      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
         probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
         rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
      }
  
   }//-----------------------------------------------------------------------------------EVENT LOOP ENDS -------------------------------------------------------------------//
      can = new TCanvas("can","",10,10,600,600);
      can->SetLeftMargin(0.2);
      can->SetRightMargin(0.05);
      can->SetTopMargin(0.05);
      can->SetBottomMargin(0.13);
      can->cd();
      gStyle->SetOptStat(0);   
      hERNR_ML->GetXaxis()->SetTitle("precoiltNF (keV)");
      hERNR_ML->GetXaxis()->SetTitleOffset(1.2);
      hERNR_ML->GetYaxis()->SetTitleOffset(1.7);
      hERNR_ML->GetYaxis()->SetTitle("qimean");
      hERNR_ML->SetMarkerColor(3);  //3=Green
      hERNR_ML->SetMarkerStyle(15);
      hERNR_ML->SetMarkerSize(1.0);
      hERNR_ML->Draw();
      
      hERNR->GetXaxis()->SetTitle("precoiltNF (keV)");
      hERNR->GetXaxis()->SetTitleOffset(1.2);
      hERNR->GetYaxis()->SetTitleOffset(1.7);
      hERNR->GetYaxis()->SetTitle("qimean");
      hERNR->SetMarkerColor(2);  //2=Red
      hERNR->SetMarkerStyle(15);
      hERNR->SetMarkerSize(1.0);
      hERNR->Draw("same");
      TLine *line1 = new TLine(10,7,170,95);
  //TLine *line1 = new TLine((3-c)/m,3,80,ylinefun(80));
 //   line1->SetLineWidth(3);
   // line1->SetLineColor(2);
 //  line1->Draw("same");
 
          
    TLegend *l1 = new TLegend(.5365772,0.7801394,.6708054,0.8797909);
   l1->SetTextFont(42);
   l1->SetBorderSize(0);
   l1->SetFillStyle(0);
   //l1->SetFillColor(0);
   l1->SetMargin(0.25);
   l1->SetTextSize(0.03);
   l1->SetEntrySeparation(0.5);
   l1->AddEntry(hERNR_ML,"NR Signal by ML ","lpf");
   l1->AddEntry(hERNR,"ER Background","lpf");
   l1->Draw("same");

      
      can->SaveAs("h2d.jpg");
      
      TCanvas *can2;
      can2 = new TCanvas("can2","",10,10,600,600);
      can2->SetLeftMargin(0.2);
      can2->SetRightMargin(0.05);
      can2->SetTopMargin(0.05);
      can2->SetBottomMargin(0.13);
      can2->cd();
      
       hyieldE_ML->GetXaxis()->SetTitle("precoiltNF (keV)");
      hyieldE_ML->GetXaxis()->SetTitleOffset(1.2);
      hyieldE_ML->GetYaxis()->SetTitleOffset(1.7);
      hyieldE_ML->GetYaxis()->SetTitle("yield");
      hyieldE_ML->SetMarkerColor(3);  //3=Green
      hyieldE_ML->SetMarkerStyle(15);
      hyieldE_ML->SetMarkerSize(1.0);
      hyieldE_ML->Draw();
      
      hyieldE->GetXaxis()->SetTitle("precoiltNF (keV)");
      hyieldE->GetXaxis()->SetTitleOffset(1.2);
      hyieldE->GetYaxis()->SetTitleOffset(1.7);
      hyieldE->GetYaxis()->SetTitle("qimean");
      hyieldE->SetMarkerColor(96);  //2=Red
      hyieldE->SetMarkerStyle(15);
      hyieldE->SetMarkerSize(1.0);
      hyieldE->Draw("same");
      
   
     TLegend *l2 = new TLegend(.5365772,0.7801394,.6708054,0.8797909);
   l2->SetTextFont(42);
   l2->SetBorderSize(0);
   l2->SetFillStyle(0);
   //l1->SetFillColor(0);
   l2->SetMargin(0.25);
   l2->SetTextSize(0.03);
   l2->SetEntrySeparation(0.5);
   l2->AddEntry(hyieldE_ML,"NR Band by ML ","p");
   l2->AddEntry(hyieldE,"ER Band by ML","p");
   l2->Draw("same");

    

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
   
   // Get efficiency for cuts classifier
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {
     
     // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer
     TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;
     
      if (mcuts) {
	std::vector<Double_t> cutsMin;
	std::vector<Double_t> cutsMax;
	mcuts->GetCuts( 0.7, cutsMin, cutsMax );
	std::cout << "--- -------------------------------------------------------------" << std::endl;
	std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
	for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
	  std::cout << "... Cut: "
		    << cutsMin[ivar]
		    << " < \""
		    << mcuts->GetInputVar(ivar)
		    << "\" <= "
		    << cutsMax[ivar] << std::endl;
	}
	std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }
   
   // Write histograms
   
   TFile *target  = new TFile( "TMVApp.root","RECREATE" );
   if (Use["Likelihood"   ])   histLk     ->Write();
   if (Use["LikelihoodD"  ])   histLkD    ->Write();
   if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
   if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
   if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
   if (Use["PDERS"        ])   histPD     ->Write();
   if (Use["PDERSD"       ])   histPDD    ->Write();
   if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
   if (Use["KNN"          ])   histKNN    ->Write();
   if (Use["HMatrix"      ])   histHm     ->Write();
   if (Use["Fisher"       ])   histFi     ->Write();
   if (Use["FisherG"      ])   histFiG    ->Write();
   if (Use["BoostedFisher"])   histFiB    ->Write();
   if (Use["LD"           ])   histLD     ->Write();
   if (Use["MLP"          ])   histNn     ->Write();
   if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
   if (Use["MLPBNN"       ])   histNnbnn  ->Write();
   if (Use["CFMlpANN"     ])   histNnC    ->Write();
   if (Use["TMlpANN"      ])   histNnT    ->Write();
   if (Use["DNN_GPU"]) histDnnGpu->Write();
   if (Use["DNN_CPU"]) histDnnCpu->Write();
  
   if (Use["BDT"          ])   {
       
       
       histBdt    ->Write();
       hERNR->Write();
       hERNR_ML->Write();
       hyieldE->Write();
       hyieldE_ML->Write();
       
}
   
   
   if (Use["BDTG"         ])   histBdtG   ->Write();
   if (Use["BDTB"         ])   histBdtB   ->Write();
   if (Use["BDTD"         ])   histBdtD   ->Write();
   if (Use["BDTF"         ])   histBdtF   ->Write();
   if (Use["RuleFit"      ])   histRf     ->Write();
   if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
   if (Use["SVM_Poly"     ])   histSVMP   ->Write();
   if (Use["SVM_Lin"      ])   histSVML   ->Write();
   if (Use["FDA_MT"       ])   histFDAMT  ->Write();
   if (Use["FDA_GA"       ])   histFDAGA  ->Write();
   if (Use["Category"     ])   histCat    ->Write();
   if (Use["Plugin"       ])   histPBdt   ->Write();

   // Write also error and significance histos
   if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }

   // Write also probability hists
   if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }
   target->Close();

   std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;

   delete reader;

   std::cout << "==> applicationRRQ is done!" << std::endl << std::endl;
   
   cout<<"The total number of entries-------------------------------------------------------------"<<theTree->GetEntries()<<endl;
   
   cout<<" Number of entries in Background Invariant Mass histogram"<<setw(10)<<hERNR_ML->Integral(1,60)<<endl;
   cout<<" Number of entries in Invariant Mass histogram           "<<setw(10)<<hERNR->Integral(1,60)<<endl;
   cout<<" Cut Value =  "<<cutvalue<<endl;
   
   
   
   fout->cd();
  tree_ER->Write();
  tree_NR->Write();
  fout->Close();
   
}


int main( int argc, char** argv )
{
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   applicationRRQ(methodList);
   return 0;
}
