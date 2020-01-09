/*
THIS CODE for charge radial cut side 2.
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

Double_t Normal(Double_t *x, Double_t *par)
{
  //return (par[0]*par[1]*0.01*1)/(2*3.14159)/((x[0]-par[2])*(x[0]-par[2])+(par[1]*par[1])/4.);
  Double_t y;
  y=par[0]*exp(  -0.5*pow( (x-par[1])/par[2] ,2 ) )      ;
 return y;
}



void cradialcut_ba2(string weekno="1", Int_t zip=4)
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
string opfileaddress="/home/u1/rik/Viraj/FV/chargeradialcut";

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

TCut virajchargeradialcut("qo1OF<1.05087 && qo2OF<1.23711");
TCut virajPhononRadial("prpartOF<(TMath::Exp(-0.424675-3.53768*ptNF)+0.295494)"); // 2 sigma

 cout<<"\nNumber of events: "<<chain2->GetEntries()<<endl;
cout<<"Name of Tree: "<<chain2->GetName()<<endl;

//chain2->Scan("SeriesNumber","SeriesNumber==11509070319", "colsize=15 precision=13");

double qi1OF_range_1=-2;
double qi1OF_range_2=10;
double qi2OF_range_1=-4; 
double qi2OF_range_2=10;
double qo1OF_range_1=-2;
double qo1OF_range_2=5;
double qo2OF_range_1=-2;
double qo2OF_range_2=10;
double qsum1OF_range_1=-4;
double qsum1OF_range_2=15;
double qsum2OF_range_1=-5;
double qsum2OF_range_2=15;
double qi1OF_bins=100;
double qo1OF_bins=100;
double qi2OF_bins=100;
double qo2OF_bins=100;
double qsum1OF_bins=100;
double qsum2OF_bins=100;

TH2D *h_charge = new TH2D("h_charge", Form("Z%d: Barium-%s, Side 1: Outer Charge vs Inner Charge",zid, filename.c_str()), qi1OF_bins, qi1OF_range_1, qi1OF_range_2, qo1OF_bins, qo1OF_range_1, qo1OF_range_2);
h_charge->SetXTitle("Inner Charge [keV]");
h_charge->SetYTitle("Outer Charge [keV]");
h_charge->GetYaxis()->SetTitleOffset(1.5);
h_charge->GetYaxis()->SetLabelSize(0.03);
h_charge->GetYaxis()->SetTitleSize(0.03);
h_charge->SetMarkerColor(kRed);
h_charge->SetFillColor(kRed);
h_charge->SetStats(kFALSE);

TCanvas *ch2d= new TCanvas("ch2d", "cTH2D outer charge vs inner charge", 600, 600);
chain2->Draw("qo2OF:qi2OF>>h_charge",cBasicCuts);
//chain2->Draw("qo2OF:qi2OF>>h_charge");
//chain2->Draw("qo1OF:qi1OF>>h_charge", cNotEmpty+cRandom+cGoodFlash+cPstd+cBaseTemp+cVoltageBias, "same");
//ch2d->SaveAs(Form("%s/OuterCharge vs Innercharge_side1_%s_%s_%s.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));



TH1D *S1proj_inb, *S1proj_onb, *S1proj_inevents;

TCanvas* CanQP2= new TCanvas("CanQP2", "S2: Outer Noise Events Projection", 600, 600);
S1proj_onb=h_charge->ProjectionY("S1proj_inevents", h_charge->GetXaxis()->FindBin(-2), h_charge->GetXaxis()->FindBin(2));
S1proj_onb->SetTitle(Form("Z%d: Barium-%s,%s, Side 2 Outer Noise Spectrum Projection",zid, filename.c_str(), weektitle.c_str()));
S1proj_onb->SetYTitle("counts");
S1proj_onb->GetYaxis()->SetTitleOffset(1.5);
S1proj_onb->GetYaxis()->SetLabelSize(0.03);
S1proj_onb->GetYaxis()->SetTitleSize(0.03);

//S1proj_onb->GetXaxis()->SetRange(-1,1);
//S1proj_onb->SetAxisRange(-0.8, 1,"X");

//S1proj_onb->SetStats(kFALSE);
S1proj_onb->Draw();
cout<<" Bin 0: "<<h_charge->GetXaxis()->FindBin(-2)<<"  Bin 1: "<< h_charge->GetXaxis()->FindBin(2)<<endl;

/*
double fitlow=-1.0;
double fithi=1.0;
 TF1 *fitFunNormal = new TF1("fitFunNormal",Normal,fitlow,fithi,3);
fitFunNormal->SetParameter(0,250);
fitFunNormal->SetParameter(1,0.5);  // mean set
fitFunNormal->SetParameter(2,0.5);  // sigma set
  //fitFunNormal->FixParameter(1,0.048);

  //fitFunNormal->SetParLimits(0,100,1000);
  //fitFunNormal->SetParLimits(1,-1,1);
  //fitFunNormal->SetParLimits(2,0,1);
  //fitFunNormal->SetLineColor(3);
  fitFunNormal->SetParNames("Amplitude","Mean","Sigma");

  S1proj_onb->Fit(fitFunNormal,"IER+");

  float A=fitFunNormal->GetParameter(0);
  float mu=fitFunNormal->GetParameter(1);
  float sigma=fitFunNormal->GetParameter(2);

*/

S1proj_onb->Fit("gaus");
gaus->SetLineColor(3);

float A=gaus->GetParameter(0);
  float mu=gaus->GetParameter(1);
  float sigma=gaus->GetParameter(2);

  cout<<" mu "<< mu<<" 3 X sigma= "<<3*sigma<<endl;
  cout<<" mu "<< mu<<"mu + 3sigma= "<<mu+3*sigma<<endl;
  
   TLegend *leg = new TLegend(0.2,0.7,0.3,0.8,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->AddEntry("S1proj_onb","Projection Outer Charge ","lpf");
   leg->AddEntry("gaus","Gaussian Fit","lpf");
  
   leg->Draw("same");
//cout<<"Outer Noise Mean: "<<S1proj_onb->GetMean()<<endl;
//cout<<"Outer Noise StdDev: "<<S1proj_onb->GetStdDev()<<endl;
//cout<<"Outer Noise Mean+5*StdDev(D)= "<<S1proj_onb->GetMean()+5*S1proj_onb->GetStdDev()<<endl;
//CanQP2->SaveAs(Form("%s/canvas_S2OuterNoiseProj_%s_%s_%s.C", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));


TCanvas *ch2dcut_side2= new TCanvas("ch2dcut_side2", "cTH2D outer charge vs inner charge", 600, 600);
ch2dcut_side2->cd();
chain2->Draw("qo2OF:qi2OF>>h_charge",cBasicCuts);
TLine *line1 = new TLine(-2,mu+3*sigma,10,mu+3*sigma);
    line1->SetLineWidth(4);
    line1->SetLineColor(3);
  line1->Draw("same");
  TLegend *legcut = new TLegend(0.5365772,0.7801394,0.6708054,0.8797909,NULL,"brNDC");
   legcut->SetBorderSize(0);
   legcut->SetTextSize(0.03);
   legcut->SetLineColor(1);
   legcut->SetLineStyle(1);
   legcut->SetLineWidth(1);
   legcut->SetFillColor(0);
   legcut->SetFillStyle(0);
   TLegendEntry *entryc=legcut->AddEntry("line1","Charge radial cut ","p");
   legcut->Draw("same");
   entryc->SetTextFont(42);
//chain2->Draw("qo1OF:qi1OF>>h_charge");
//chain2->Draw("qo1OF:qi1OF>>h_charge", cNotEmpty+cRandom+cGoodFlash+cPstd+cBaseTemp+cVoltageBias, "same");
//ch2dcut->SaveAs(Form("%s/canvas_chargeradialcut_OuterCharge vs Innercharge_side2_%s_%s_%s.C", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));
//ch2dcut->SaveAs(Form("%s/canvas_cut_OuterCharge vs Innercharge_side2_%s_%s_%s.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));


   
}
