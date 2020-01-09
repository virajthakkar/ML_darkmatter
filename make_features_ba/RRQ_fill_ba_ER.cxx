/*
THIS CODE IS TO RRQ plots USING Ba DATA.
data is in SLAC
*/
#include "example_loadAnalysisFilesnewba.C" //Load Prodv5-6-3
#include "BasicCuts.h"                    //Load only Basic Cuts
#include <iostream>
#include <fstream>
#include <iomanip>
#include <TMath.h>
#include <iostream>
#include <fstream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>


void RRQ_fill_ba_ER(string weekno="1", Int_t zip=4)
{
TGaxis::SetMaxDigits(3);
gStyle->SetLabelSize(0.025, "Z");

//Define Filename variable, Not used anywhere, if unnecesary then discard this line later
string filename,fileaddress;
filename="ba";
fileaddress="ba";

//Convert weekno from string to int for wherever it might be convenient
Int_t wno= atoi(Form("%s", weekno.c_str()));

//Define source for proper naming of output file
string source, sourcefullname;


source="Ba";
sourcefullname="Barium";
//Initialize zlite
Int_t zlite;
if(zip==4)
{zlite=14;}
else
{
 if(zip==14)
 {zlite=4;}
 else
 {
  if(zip==5  || weekno=="Allz4" )
  {
      if(wno>=1 && wno<7)
      zlite=14; cout<<"zlite: "<< zlite<< endl;
      
  }
  else
  {
   if(zip==5  || weekno=="Allz14")
   { 
       if(wno>6 && wno<=12)
       zlite=4; cout<<"zlite: "<< zlite<< endl;
   }
   else
   {
    cout<<"Incorrect zip or week number. Please select zip4, zip5 or zip14 and week numbers between 1-12, or All or Allz4 or Allz14"<<endl;
    
   }
  }
 }
}

//Initialize weekno2
string weekno2=weekno;
if(weekno=="All" || weekno=="Allz4" || weekno=="Allz14")
{
weekno="*";
}

if(weekno2=="Allz4" && zip==4) //Added for convenience in running shell script
{
cout<<"Incorrect zip or week number. Please select zip4, zip5 or zip14 and week numbers between 1-12, or All or Allz4 or Allz14"<<endl;
    
}

if(weekno2=="Allz14" && zip==14)  //Added for convenience in running shell script
{
cout<<"Incorrect zip or week number. Please select zip4, zip5 or zip14 and week numbers between 1-12, or All or Allz4 or Allz14"<<endl;
    
}

//Initialize det
string det;
if(zip==4)
{det="T2Z1";}
else
{
 if(zip==14)
  {det="T5Z2";}
   else
    {
      if(zip==5)
        {det="T2Z2";}
          else
            {cout<<"Please enter either zip 4, 5 or 14";
             
            }
        }
    }
//Initialize weektitle, Not used in this script but maybe useful in others
string weektitle;
if(weekno2=="All" && det=="T2Z1")
{
weektitle="Weeks1-6";
}
else
{

 if(weekno2=="All" && det=="T5Z2")
 {
  weektitle="Weeks7-12";
 }
 else
 {
  if(det=="T2Z2" && weekno2=="Allz4")
  {weektitle="Weeks5-6";}
  else
  {
   if(det=="T2Z2" && weekno2=="Allz14")
   {weektitle="Weeks7-12";}
   else
   {
    if(source=="Sb" && weekno2=="All")
    {//For Sb 
     weektitle="Weeks1-9";
    }
    else
    {weektitle="Week";
     weektitle=weektitle+weekno2;
    }
   }
  }
 }
}



//--------------------- Output File Address -------------------------
//string opfileaddress=Form("/home/viraj/Documents/DarkMatter/HT_Analysis/mf%s", det.c_str());
string opfileaddress="/home/u1/rik/Viraj/ba_bkg";

//Setting up detector indices
Int_t DetIndex = zip-1;
Int_t tid = (DetIndex/3)+1;
Int_t zid = (DetIndex%3)+1;
cout<<"\nWe are looking at: T"<<tid<<"Z"<<zid<<endl;

//Confirming BiasFlashtime
//cout<<"BiasFlashTime= "<<flashtime()<<endl;

// ----------------------- Location of Data --------------------------
// --- Loading Prodv5-6-3 data ---
//string dataDir = "/home/viraj/Documents/DarkMatter/HT_Analysis/";
string dataDir = "/data/R134/dataReleases/Prodv5-3-5/merged/all/";

string fileaddress2;
fileaddress2=fileaddress;
dataDir = dataDir+fileaddress2;
string seriesweek="week";
seriesweek=seriesweek+weekno;

//Call function that loads data
TChain *chain2 = chainDataAllSpecial(zip, dataDir, weekno, fileaddress, zlite);
//chain2->Scan();
cout<<"Number of events: "<<chain2->GetEntries()<<endl;
cout<<"Name of Tree: "<<chain2->GetName()<<endl<<endl;

//Define Variables
double ptNFrange_1=-10; // Prev: -10, 0
double ptNFrange_2=170; // Prev: 170, 30
double ptNFbinsize=450; // Prev: 450, 300
double qbinsize=500; // Prev: 500, 160
double q_y_range_1=-10; //Prev: -10, -2
double q_y_range_2= 80; //Prev: 80, 14

const double binw= (ptNFrange_2-ptNFrange_1)/ptNFbinsize;


//Define some Quality cuts
TCut cLEChisq("PTNFchisq < (0.015*ptNF*ptNF+4500) && PTNFchisq >3600");
TCut cDeltaChisqGlitch("(PTOFchisq-PTglitch1OFchisq)<((-2.3*(ptNF*ptNF))+25)");
TCut cDeltaChiSqLFNLine("(PTOFchisq-PTlfnoise1OFchisq)< 14.2");
TCut cDeltaChiSqLFNParabola("(PTOFchisq-PTlfnoise1OFchisq)<((-13.5152*(ptNF*ptNF))-(2.18824*ptNF)+27.692)");
TCut cDeltaChiSqLFN=cDeltaChiSqLFNLine+cDeltaChiSqLFNParabola;
TCut cChargeChiSqS1("QS1OFchisq <= ((( 0.00578106)*qsum1OF*qsum1OF*qsum1OF)-( 0.427417*qsum1OF*qsum1OF)+( 9.88895*qsum1OF)+ 5448.68)");
TCut cChargeChiSqS2("QS2OFchisq <= ((( 0.00327275)*qsum2OF*qsum2OF*qsum2OF)-( 0.25359*qsum2OF*qsum2OF)+( 7.05146*qsum2OF)+ 6241.38)");
TCut cChargeChiSq=cChargeChiSqS1+cChargeChiSqS2;

TCut cRadialCutS1("qo1OF<3.0 || qi1OF>3.0");
TCut cRadialCutS2("qo2OF<1.8 || qi2OF>3.0");
TCut cRadialCut=cRadialCutS1+cRadialCutS2;

TCut cFVPolyUp("qi2OF<(-0.00233946*qi1OF*qi1OF)+1.11534*qi1OF+2.06405");
TCut cFVPolyLow("qi2OF>(-0.0035754*qi1OF*qi1OF)+1.04275*qi1OF-2.85197");
TCut cFVPoly=cFVPolyUp+cFVPolyLow;
TCut cHorizontalLine("qi2OF<3.5"); //mu+5.5*sigma of qi2OF
TCut cVerticalLine("qi1OF<3.5"); // mu+13*sigma of qi1OF
TCut cFlatLines=cHorizontalLine+cVerticalLine;
TCut cSymmetricFV=(cFVPoly)||(cFlatLines);

TCut cZeroCharge("((qi1OF+qi2OF)/2)<1 && ((qi1OF+qi2OF)/2)>-1");
TCut cPhononRadial("prpartOF<(TMath::Exp(-2.42429-0.0990673*ptNF)+0.241877)"); // 2 sigma
TCut cNotPhononRadial("prpartOF>(TMath::Exp(-2.42429-0.0990673*ptNF)+0.241877)"); // 2 sigma

//TCut threesigmacut("( (qi1OF+qi2OF)/2.0 ) > (70/170)*ptNF -10");  // simulation for 3 sigma NR band cut
TCut threesigmacut("( (qi1OF+qi2OF)/2.0 ) > 0.95*precoiltNF -10 && ((qi1OF+qi2OF)/2.0) >3 ");  // simulation for 3 sigma NR band cut with ptNF on x axis

TCut cQualityCuts=cLEChisq+cDeltaChisqGlitch+cDeltaChiSqLFN+cChargeChiSq+cSymmetricFV;



 cout<<"\nNumber of events: "<<chain2->GetEntries()<<endl;
cout<<"Name of Tree: "<<chain2->GetName()<<endl;

//chain2->Scan("SeriesNumber","SeriesNumber==11509070319", "colsize=15 precision=13");



double livetime;
double Empty, EventCategory, volt;
double livetimesec;
double netlivetime = 0.0;
//chain2->SetBranchStatus("*",0);

//reactivate only the branches we need before setting branch adress
//chain2->SetBranchStatus("LiveTime",1);
chain2->SetBranchAddress("LiveTime", &livetime);


double precoiltNF_local,ytNF_local, qzpartOF_local,qrpart1OF_local,pzpartOF_local,prpartOF_local,qsummaxOF_local;
  
double PTNFchisq_local,
ptNF_local,
PTOFchisq_local,
PTglitch1OFchisq_local,
PTlfnoise1OFchisq_local,
QS1OFchisq_local,
qsum1OF_local,
QS2OFchisq_local,
qsum2OF_local,
qo1OF_local,
qo2OF_local,
qi1OF_local,
qi2OF_local;


chain2->SetBranchAddress("PTNFchisq",&PTNFchisq_local);
chain2->SetBranchAddress("ptNF",&ptNF_local);
chain2->SetBranchAddress("PTOFchisq",&PTOFchisq_local);
chain2->SetBranchAddress("PTglitch1OFchisq",&PTglitch1OFchisq_local);
chain2->SetBranchAddress("PTlfnoise1OFchisq",&PTlfnoise1OFchisq_local);
chain2->SetBranchAddress("QS1OFchisq",&QS1OFchisq_local);
chain2->SetBranchAddress("qsum1OF",&qsum1OF_local);
chain2->SetBranchAddress("QS2OFchisq",&QS2OFchisq_local);
chain2->SetBranchAddress("qsum2OF",&qsum2OF_local);
chain2->SetBranchAddress("qo1OF",&qo1OF_local);
chain2->SetBranchAddress("qo2OF",&qo2OF_local);
chain2->SetBranchAddress("qi1OF",&qi1OF_local);
chain2->SetBranchAddress("qi2OF",&qi2OF_local);



chain2->SetBranchAddress("precoiltNF",&precoiltNF_local);
chain2->SetBranchAddress("ytNF", &ytNF_local);
chain2->SetBranchAddress("qzpartOF", &qzpartOF_local);
chain2->SetBranchAddress("qrpart1OF", &qrpart1OF_local);
chain2->SetBranchAddress("pzpartOF", &pzpartOF_local);
chain2->SetBranchAddress("prpartOF", &prpartOF_local);
chain2->SetBranchAddress("qsummaxOF", &qsummaxOF_local);
chain2->SetBranchAddress("EventCategory", &EventCategory);


int ctr = 0;
double flashtime;
//chain2->SetBranchStatus("BiasFlashTime",1);
chain2->SetBranchAddress("BiasFlashTime", &flashtime);

float precoiltNF_b,ytNF_b, qzpartOF_b,qrpart1OF_b,pzpartOF_b,prpartOF_b,qsummaxOF_b;
float PTNFchisq_b,ptNF_b,PTOFchisq_b,PTglitch1OFchisq_b,PTlfnoise1OFchisq_b,
QS1OFchisq_b,qsum1OF_b,QS2OFchisq_b,qsum2OF_b,qo1OF_b,qo2OF_b,qi1OF_b,qi2OF_b,qimean_b;


 TFile *fout = new TFile("RRQ_ER.root","RECREATE");
  TTree *tree_bkg = new TTree("tree_bkg","Bkg tree");
  tree_bkg->Branch("precoiltNF",&precoiltNF_b,"precoiltNF_b/F");
  tree_bkg->Branch("ytNF",&ytNF_b,"ytNF_b/F");
  tree_bkg->Branch("qzpartOF",&qzpartOF_b,"qzpartOF_b/F");
  tree_bkg->Branch("qrpart1OF",&qrpart1OF_b,"qrpart1OF_b/F");
  tree_bkg->Branch("pzpartOF",&pzpartOF_b,"pzpartOF_b/F");
  tree_bkg->Branch("prpartOF",&prpartOF_b,"prpartOF_b/F");
  tree_bkg->Branch("qsummaxOF",&qsummaxOF_b,"qsummaxOF_b/F");
  
  tree_bkg->Branch("PTNFchisq",&PTNFchisq_b,"PTNFchisq_b/F");
  tree_bkg->Branch("ptNF",&ptNF_b,"ptNF_b/F");
  tree_bkg->Branch("PTOFchisq",&PTOFchisq_b,"PTOFchisq_b/F");
  tree_bkg->Branch("PTglitch1OFchisq",&PTglitch1OFchisq_b,"PTglitch1OFchisq_b/F");
  tree_bkg->Branch("PTlfnoise1OFchisq",&PTlfnoise1OFchisq_b,"PTlfnoise1OFchisq_b/F");
  tree_bkg->Branch("QS1OFchisq",&QS1OFchisq_b,"QS1OFchisq_b/F");
  tree_bkg->Branch("qsum1OF",&qsum1OF_b,"qsum1OF_b/F");
  tree_bkg->Branch("QS2OFchisq",&QS2OFchisq_b,"QS2OFchisq_b/F");
  tree_bkg->Branch("qsum2OF",&qsum2OF_b,"qsum2OF_b/F");
  tree_bkg->Branch("qo1OF",&qo1OF_b,"qo1OF_b/F");
  tree_bkg->Branch("qo2OF",&qo2OF_b,"qo2OF_b/F");
  tree_bkg->Branch("qi1OF",&qi1OF_b,"qi1OF_b/F");
  tree_bkg->Branch("qi2OF",&qi2OF_b,"qi2OF_b/F");
  tree_bkg->Branch("qimean",&qimean_b,"qimean_b/F");
  
  
  
 
 int filledcount;
 filledcount=0;
  
  
double precoiltNFrange_1=10; 
double precoiltNFrange_2=170;
double precoiltNFbinsize=8; 


double qzpartOFrange_1=-0.3; 
double qzpartOFrange_2=0.3; 


double qrpart1OFrange_1=-0.2; 
double qrpart1OFrange_2=0.6; 


double ytNFrange_1=0.6; 
double ytNFrange_2=1.4; 


double prpartOFrange_1=0.1; 
double prpartOFrange_2=0.3; 


double pzpartOFrange_1=-0.1; 
double pzpartOFrange_2=0.2;

double qsummaxOFrange_1=0; 
double qsummaxOFrange_2=150; 

/*TH1D *h_precoiltNF = new TH1D("h_precoiltNF", Form("T%dZ%d: %s-%s, %s,  precoiltNF", tid, zid, sourcefullname.c_str(), filename.c_str(), weektitle.c_str()), precoiltNFbinsize, precoiltNFrange_1, precoiltNFrange_2);
h_precoiltNF->SetXTitle("precoiltNF [keV]");
h_precoiltNF->SetYTitle("Counts");
h_precoiltNF->GetYaxis()->SetTitleOffset(1.5);
h_precoiltNF->GetYaxis()->SetLabelSize(0.03);
h_precoiltNF->GetYaxis()->SetTitleSize(0.03);
*/

/*double ptNFrange_1=-10; // Prev: -10, 0
double ptNFrange_2=170; // Prev: 170, 30
double ptNFbinsize=450; // Prev: 450, 300
double qbinsize=500; // Prev: 500, 160
double q_y_range_1=-10; //Prev: -10, -2
double q_y_range_2= 80; //Prev: 80, 14

const double binw= (ptNFrange_2-ptNFrange_1)/ptNFbinsize;
double qimean=0;

TCanvas* qptNF= new TCanvas("qptNF", "qimean vs precoiltNF (RRQ_fill script)", 600, 600);
TH2D *h_ERNRallcuts = new TH2D("h_ERNRallcuts", Form("T%dZ%d: %s-%s, %s, Ionization vs precoiltNF", tid, zid, sourcefullname.c_str(), filename.c_str(), weektitle.c_str()), ptNFbinsize, ptNFrange_1, ptNFrange_2, qbinsize, q_y_range_1, q_y_range_2);
h_ERNRallcuts->SetXTitle("precoiltNF [keV]");
h_ERNRallcuts->SetYTitle("qimean");
h_ERNRallcuts->GetYaxis()->SetTitleOffset(1.5);
h_ERNRallcuts->GetYaxis()->SetLabelSize(0.03);
h_ERNRallcuts->GetYaxis()->SetTitleSize(0.03);
*/
//TCanvas* c_precoiltNF= new TCanvas("c_precoiltNF", "c_precoiltNF (All Cuts) 1", 600, 600);
cout<<"start while loop"<<endl;
//while(chain2->GetEntry(ctr))
while(chain2->GetEntry(ctr))
{
   // cout<<"outside if condition "<<"precoiltNF: "<<precoiltNF_local<<endl;
if(ctr%50000==0)cout<<"Number of events processed =========== "<<ctr<<endl;
   // if(/*flashtime < 3300.00 && EventCategory!=1 && */ precoiltNF_local>10)//T2Z1: 3300 , T5Z2: 10800
      
   //  if(PTNFchisq_local < (0.015*ptNF_local*ptNF_local+4500) && PTNFchisq_local >3600
    //   && (PTOFchisq_local-PTglitch1OFchisq_local)<((-2.3*(ptNF_local*ptNF_local))+25)
     //  && (PTOFchisq_local-PTlfnoise1OFchisq_local)< 14.2 
     //  && (PTOFchisq_local-PTlfnoise1OFchisq_local)<((-13.5152*(ptNF_local*ptNF_local))-(2.18824*ptNF_local)+27.692)
     //  && QS1OFchisq_local <= ((( 0.00578106)*qsum1OF_local*qsum1OF_local*qsum1OF_local)-( 0.427417*qsum1OF_local*qsum1OF_local)+( 9.88895*qsum1OF_local)+ 5448.68)
     //  && QS2OFchisq_local <= ((( 0.00327275)*qsum2OF_local*qsum2OF_local*qsum2OF_local)-( 0.25359*qsum2OF_local*qsum2OF_local)+( 7.05146*qsum2OF_local)+ 6241.38)
     //  && (qo1OF_local<3.0 || qi1OF_local>3.0) && (qo2OF_local<1.8 || qi2OF_local>3.0)
  //&& ( (/*parallel*/ qi2OF_local<(-0.00233946*qi1OF_local*qi1OF_local)+1.11534*qi1OF_local+2.06405 && qi2OF_local>(-0.0035754*qi1OF_local*qi1OF_local)+1.04275*qi1OF_local-2.85197 )

//|| (qi2OF_local<3.5 &&  qi1OF_local<3.5) /*cSymmetricFV*/) 
       
// && (  (  ((qi1OF_local+qi2OF_local)/2.0 ) > 0.95*precoiltNF_b -10 ) && ((qi1OF_local+qi2OF_local)/2.0) >3  )  
     // ) //if cuts ends
    // {
       if(precoiltNF_local>precoiltNFrange_1 && precoiltNF_local<precoiltNFrange_2 && ((qi1OF_local+qi2OF_local)/2.0 ) < 160 
            && ytNF_local>0 && ytNF_local<2 
        && qzpartOF_local>qzpartOFrange_1 && qzpartOF_local< qzpartOFrange_2 && qrpart1OF_local> qrpart1OFrange_1 && qrpart1OF_local < qrpart1OFrange_2
        && prpartOF_local> prpartOFrange_1 && prpartOF_local < prpartOFrange_2 && pzpartOF_local>pzpartOFrange_1 && pzpartOF_local<pzpartOFrange_2
        //&& qsummaxOF_local > qsummaxOFrange_1 && qsummaxOF_local < qsummaxOFrange_2
         )// for range of features T2Z1: 3300 , T5Z2: 10800
        {
        filledcount++;
      //  cout<<"livetime : "<<livetime<<endl;
    //   cout<<" precoiltNF : "<<precoiltNF_local<<endl;
      //  h_precoiltNF->Fill(precoiltNF_local);
        
       precoiltNF_b=precoiltNF_local;
       ytNF_b=ytNF_local;
       qzpartOF_b=qzpartOF_local;
       qrpart1OF_b=qrpart1OF_local;
       pzpartOF_b=pzpartOF_local;
       prpartOF_b=prpartOF_local;
       qsummaxOF_b=qsummaxOF_local;
       
       PTNFchisq_b=PTNFchisq_local;
       ptNF_b=ptNF_local;
       PTOFchisq_b=PTOFchisq_local;
       PTglitch1OFchisq_b=PTglitch1OFchisq_local;
       PTlfnoise1OFchisq_b=PTlfnoise1OFchisq_local;
       QS1OFchisq_b=QS1OFchisq_local;
       qsum1OF_b=qsum1OF_local;
       QS2OFchisq_b=QS2OFchisq_local;
       qsum2OF_b=qsum2OF_local;
       qo1OF_b=qo1OF_local;
       qo2OF_b=qo2OF_local;
       qi1OF_b=qi1OF_local;
       qi2OF_b=qi2OF_local;
       qimean_b=0.5*(qi1OF_local+qi2OF_local);

       
    //qimean = 0.5*( qi1OF_local+ qi2OF_local  );
       //h_ERNRallcuts->Fill( precoiltNF_local ,  qimean );
        tree_bkg->Fill();
        }  // range of variables IF
        
  //  } // TCut IF
 
    ctr++;
   // if(ctr>1000000) break;
    
} // while loop ends
//break;
//h_precoiltNF->Draw();
//cout<<"h integral :"<<h_precoiltNF->Integral()<<endl;
/*TLegend *leg=new TLegend(0.55,0.57,0.90,0.77);
leg->SetBorderSize(0);
leg->AddEntry(h_precoiltNF,"after Quality+PRC cuts+ 3sigma NR cut","lpf");
leg->Draw("same");

c_precoiltNF->Update();
*/

cout<<"ctr (counts) "<<ctr<<endl;
cout<<"filled (counts) "<<filledcount<<endl;

  fout->cd();
  tree_bkg->Scan();
  tree_bkg->Write();
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("precoiltNF"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("ytNF"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("qzpartOF"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("qrpart1OF"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("pzpartOF"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("prpartOF"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("qsummaxOF"));

  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("PTNFchisq"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("ptNF"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("PTOFchisq"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("PTglitch1OFchisq"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("PTlfnoise1OFchisq"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("QS1OFchisq"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("qsum1OF"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("QS2OFchisq"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("qsum2OF"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("qo1OF"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("qo2OF"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("qi1OF"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("qi2OF"));
  tree_bkg->ResetBranchAddress(tree_bkg->GetBranch("qimean"));
 
  
  tree_bkg->Scan();
  fout->Close();

/*
Double_t ylinefun(double x) 
{ 
    return 0.95*x + 10 ; }
 const double m= 0.95; //slope 
const double c=-10; // y intecept

  //TLine *line1 = new TLine(0,-10,170,60);
  //TLine *line1 = new TLine((3-c)/m,3,80,ylinefun(80));
   // line1->SetLineWidth(4);
   // line1->SetLineColor(2);
  //line1->Draw("same");
//TLine *line2 = new TLine(0,3,(3-c)/m,3); // y=3 line
  //  line2->SetLineWidth(4);
    //line2->SetLineColor(2);
  //line2->Draw("same");
h_ERNRallcuts->Print();


qptNF->Update();

// Save the output plot in .gif form
qptNF->SaveAs(Form("qiprecoiltNF_RRQfill_%s_%s_%s.gif", fileaddress.c_str(),det.c_str(),weektitle.c_str()));
*/

   
}
