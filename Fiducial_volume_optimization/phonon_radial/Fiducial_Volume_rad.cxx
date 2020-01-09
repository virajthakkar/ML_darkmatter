#include "example_loadAnalysisFilesnew.C" //Load Prodv5-6-3
#include "BasicCuts.h"                    //Load only Basic Cuts
#include "QualityCuts.h"                  //Load Quality Cuts
#include "DeltaChiSq.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <TMath.h>


void Fiducial_Volume_rad(string fileaddress="ybe", string weekno="All", Int_t zip=4)
{

TGaxis::SetMaxDigits(3);
gStyle->SetLabelSize(0.025, "Z");

//Define Filename variable, Not used anywhere, if unnecesary then discard this line later
string filename;
if(fileaddress=="ybe")
{
filename="Beryllium";
}

if(fileaddress=="yblank")
{
filename="Blank";
}

if(fileaddress!="yblank" && fileaddress!="ybe")
{
cout<<"ERROR ! Please enter either ybe or yblank as first argument."<<endl;
break;
}


//Convert weekno from string to int for wherever it might be convenient
Int_t wno= atoi(Form("%s", weekno.c_str()));

//Output file address
string opfileaddress=Form("/home/vijay/Analysis/Photoneutrons/Ebook/190916_VJ/Fiducial_Volume");
 
//Define source for proper naming of output file
string source="Y";

//Initialize zlite
Int_t zlite;
if(zip==4)
{zlite=14;}
if(zip==14)
{zlite=4;}

/*
if(zip==5 && (wno>=1 && wno<7) || weekno=="Allz4" )
{zlite=14; cout<<"zlite: "<< zlite<< endl;}

if(zip==5 && (wno>6 && wno<=12) || weekno=="Allz14")
{zlite=4; cout<<"zlite: "<< zlite<< endl;}
*/



//Initialize weekno2
string weekno2=weekno;
if(weekno=="All" || weekno=="Allz4" || weekno=="Allz14")
{
weekno="*";
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
             break;
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
                                                                              
//Initialize fileaddress2
string fileaddress2;
if(fileaddress=="ybe")
{fileaddress2="YBe";}
else
{
 if(fileaddress=="yblank")
 {fileaddress2="YBlank";}
 else
  {
   if(fileaddress=="sbbe")
   {fileaddress2="SbBe";}
   else
   {
    if(fileaddress=="sbblank")
    {fileaddress2="SbBlank";}
    else
    {
    cout<<"Please enter either ybe, yblank, sbbe or sbblank"<<endl;
    break;
    }
   }
  }
}

//Set Bias flashtime per detector
if(zip==4)
{x=3300;}
else
{
 if(zip==14)
 {x=10800;}
}

// Set Pstd Cut
if(zip==4)
{pstd=0;}
else
{pstd=1;}

//Setting up detectot indices
Int_t DetIndex = zip-1;
Int_t tid = (DetIndex/3)+1;
Int_t zid = (DetIndex%3)+1;
cout<<"\nWe are looking at: T"<<tid<<"Z"<<zid<<endl;

//Confirming BiasFlashtime
cout<<"BiasFlashTime= "<<flashtime()<<endl;

if(returnpstd()==0)
{cout<<"Using Pstd cuts defined for T2Z1"<<endl;}
else
{cout<<"Using Pstd cuts defined for T5Z2"<<endl;}

// --- Loading Prdv5-6-3 data ---
string dataDir = "/home/vijay/Analysis/Photoneutrons/Data/merged/all/";

dataDir = dataDir+fileaddress2;
string seriesweek="week";
seriesweek=seriesweek+weekno;

// --- Loading Prodv5-6-3 data ---
string dataDir2 = "/home/vijay/Analysis/Photoneutrons/Data/merged/all/YBlank";
string fileaddress3="yblank";

//Call function that loads data
TChain *chain2 = chainDataAllSpecial(zip, dataDir, seriesweek, fileaddress, zlite);
cout<<"Number of events: "<<chain2->GetEntries()<<endl;
cout<<"Name of Tree: "<<chain2->GetName()<<endl;

TChain *chain1 = chainDataAllSpecial(zip, dataDir2, seriesweek, fileaddress3, zlite);
cout<<"Number of events: "<<chain1->GetEntries()<<endl;
cout<<"Name of Tree: "<<chain1->GetName()<<endl;

//Define Variables
double qi1OF_range_1=-2;
double qi1OF_range_2=10;
double qi2OF_range_1=-4; 
double qi2OF_range_2=10;
double qo1OF_range_1=-2;
double qo1OF_range_2=10;
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

//Variables to calculate Yield
Double_t k=0.157; // For Ge: Z=32, A=72.64
Double_t Ernr=8.1;
Double_t g, epsilon, Y;
epsilon=11.5*Ernr*(3.076e-4);
g=3*TMath::Power(epsilon,0.15)+0.7*TMath::Power(epsilon,0.6)+epsilon;
Y=k*(double)g/(1+(k*g));
cout<<endl<<"Using the following parameters to calculate total phonon energy:"<<endl;
cout<<"E(r,nr): "<<Ernr<<"\t"<<"Y: "<<Y<<"\t"<<"g: "<<g<<"\t"<<"epsilon: "<<epsilon<<endl;

Double_t Et, Eree, V, epsl;
V=4;
epsl=3;
Et=Ernr*(1+(V*Y/epsl));
cout<<"The total phonon energy(Et) calculated using the following parameters:"<<endl;
cout<<"Et: "<<Et<<"\t"<<"V: "<<V<<"\t"<<"epsl: "<<epsl<<"\t"<<"Y: "<<Y<<endl;

Eree=Et/(1+(V/epsl));
cout<<"Single scatter recoil energy from Y/Be neutrons in Electron Equivalent energy scale:"<<endl;
cout<<"E(r,ee): "<<Eree<<endl<<endl;

//break;

// Define variables for TGraph
const Int_t n=9;
const Int_t n2=11;
Double_t Qo[n], Qi[n];
int i1;
for (Int_t i=0;i<n;i++)
{
i1=i-2;
Qi[i]=i1;
Qo[i]=Eree-Qi[i];
//cout<<"Qo = "<<Qo[i]<<"\t Qi = "<<Qi[i]<<"\t Qo+Qi = "<<Qo[i]+Qi[i]<<endl;
}
//cout<<endl;
Double_t Qo2[n2], Qi2[n2];
int i2;
for (Int_t i=0;i<n2;i++)
{
i2=i-4;
Qi2[i]=i2;
Qo2[i]=Eree-Qi2[i];
//cout<<"Qo2 = "<<Qo2[i]<<"\t Qi2 = "<<Qi2[i]<<"\t Qo2+Qi2 = "<<Qo2[i]+Qi2[i]<<endl;
}
//break;

//Define Histograms
TH2D *h_charge = new TH2D("h_charge", Form("T%dZ%d: Yttrium-%s, Side 1: Outer Charge vs Inner Charge", tid, zid, filename.c_str()), qi1OF_bins, qi1OF_range_1, qi1OF_range_2, qo1OF_bins, qo1OF_range_1, qo1OF_range_2);
h_charge->SetXTitle("Inner Charge [keV]");
h_charge->SetYTitle("Outer Charge [keV]");
h_charge->GetYaxis()->SetTitleOffset(1.5);
h_charge->GetYaxis()->SetLabelSize(0.03);
h_charge->GetYaxis()->SetTitleSize(0.03);
h_charge->SetMarkerColor(kRed);
h_charge->SetFillColor(kRed);
h_charge->SetStats(kFALSE);

TH2D *h_charge2 = new TH2D("h_charge2", Form("T%dZ%d: Yttrium-%s, Side 1: Outer Charge vs Inner Charge", tid, zid, filename.c_str()), qi1OF_bins, qi1OF_range_1, qi1OF_range_2, qo1OF_bins, qo1OF_range_1, qo1OF_range_2);
h_charge2->SetXTitle("Inner Charge [keV]");
h_charge2->SetYTitle("Outer Charge [keV]");
h_charge2->GetYaxis()->SetTitleOffset(1.5);
h_charge2->GetYaxis()->SetLabelSize(0.03);
h_charge2->GetYaxis()->SetTitleSize(0.03);
h_charge2->SetMarkerColor(kRed+3);
h_charge2->SetFillColor(kRed+3);
//h_charge2->SetMarkerStyle(29);
//h_charge2->SetMarkerSize(1);
h_charge2->SetStats(kFALSE);

//TH2D *h_charge3 = new TH2D("h_charge3", Form("T%dZ%d: Yttrium-%s, Side 2: Outer Charge vs Inner Charge", tid, zid, filename.c_str()), qi2OF_bins, qi2OF_range_1, qi2OF_range_2, qo2OF_bins, qo2OF_range_1, qo2OF_range_2);
TH2D *h_charge3 = new TH2D("h_charge3", Form("T%dZ%d: Yttrium-%s, Side 2: Outer Charge vs Inner Charge", tid, zid, filename.c_str()), qi2OF_bins, -3, 3, qo2OF_bins, -1, 1);
h_charge3->SetXTitle("Inner Charge [keV]");
h_charge3->SetYTitle("Outer Charge [keV]");
h_charge3->GetYaxis()->SetTitleOffset(1.5);
h_charge3->GetYaxis()->SetLabelSize(0.03);
h_charge3->GetYaxis()->SetTitleSize(0.03);
h_charge3->SetMarkerColor(kRed);
h_charge3->SetFillColor(kRed);
h_charge3->SetStats(kFALSE);

TH2D *h_charge4 = new TH2D("h_charge4", Form("T%dZ%d: Yttrium-%s, Side 2: Outer Charge vs Inner Charge", tid, zid, filename.c_str()), qi2OF_bins, qi2OF_range_1, qi2OF_range_2, qo2OF_bins, qo2OF_range_1, qo2OF_range_2);
h_charge4->SetXTitle("Inner Charge [keV]");
h_charge4->SetYTitle("Outer Charge [keV]");
h_charge4->GetYaxis()->SetTitleOffset(1.5);
h_charge4->GetYaxis()->SetLabelSize(0.03);
h_charge4->GetYaxis()->SetTitleSize(0.03);
h_charge4->SetMarkerColor(kRed+3);
h_charge4->SetFillColor(kRed+3);
h_charge4->SetStats(kFALSE);


TH2D *h_charge5 = new TH2D("h_charge5", Form("T%dZ%d: Yttrium-%s, Side 1 Charge vs Side 2 Charge", tid, zid, filename.c_str()), qi2OF_bins, qsum1OF_range_1, qsum1OF_range_2, qo2OF_bins, qsum2OF_range_1, qsum2OF_range_2);
h_charge5->SetXTitle("Side 1 Charge [keV]");
h_charge5->SetYTitle("Side 2 Charge [keV]");
h_charge5->GetYaxis()->SetTitleOffset(1.5);
h_charge5->GetYaxis()->SetLabelSize(0.03);
h_charge5->GetYaxis()->SetTitleSize(0.03);
h_charge5->SetMarkerColor(kRed);
h_charge5->SetFillColor(kRed);
h_charge5->SetStats(kFALSE);


TH2D *h_charge6 = new TH2D("h_charge6", Form("T%dZ%d: Yttrium-%s, Side 1 Charge vs Side 2 Charge", tid, zid, filename.c_str()), qi2OF_bins, qsum1OF_range_1, qsum1OF_range_2, qo2OF_bins, qsum2OF_range_1, qsum2OF_range_2);
h_charge6->SetXTitle("Side 1 Charge [keV]");
h_charge6->SetYTitle("Side 2 Charge [keV]");
h_charge6->GetYaxis()->SetTitleOffset(1.5);
h_charge6->GetYaxis()->SetLabelSize(0.03);
h_charge6->GetYaxis()->SetTitleSize(0.03);
h_charge6->SetMarkerColor(kRed+3);
h_charge6->SetFillColor(kRed+3);
h_charge6->SetStats(kFALSE);

// Define Histograms that will have Radial Cut
//TH2D *h_charge12, *h_charge34, *h_charge56;

TH2D *h_charge12 = new TH2D("h_charge12", Form("T%dZ%d: Yttrium-%s, Side 1: Outer Charge vs Inner Charge", tid, zid, filename.c_str()), qi1OF_bins, qi1OF_range_1, qi1OF_range_2, qo1OF_bins, qo1OF_range_1, qo1OF_range_2);
h_charge12->SetXTitle("Inner Charge [keV]");
h_charge12->SetYTitle("Outer Charge [keV]");
h_charge12->GetYaxis()->SetTitleOffset(1.5);
h_charge12->GetYaxis()->SetLabelSize(0.03);
h_charge12->GetYaxis()->SetTitleSize(0.03);
h_charge12->SetMarkerColor(kBlack);
h_charge12->SetFillColor(kBlack);
//h_charge12->SetMarkerStyle(29);
//h_charge12->SetMarkerSize(1);
h_charge12->SetStats(kFALSE);

TH2D *h_charge34 = new TH2D("h_charge34", Form("T%dZ%d: Yttrium-%s, Side 2: Outer Charge vs Inner Charge", tid, zid, filename.c_str()), qi2OF_bins, qi2OF_range_1, qi2OF_range_2, qo2OF_bins, qo2OF_range_1, qo2OF_range_2);
h_charge34->SetXTitle("Inner Charge [keV]");
h_charge34->SetYTitle("Outer Charge [keV]");
h_charge34->GetYaxis()->SetTitleOffset(1.5);
h_charge34->GetYaxis()->SetLabelSize(0.03);
h_charge34->GetYaxis()->SetTitleSize(0.03);
h_charge34->SetMarkerColor(kBlack);
h_charge34->SetFillColor(kBlack);
//h_charge34->SetMarkerStyle(29);
//h_charge34->SetMarkerSize(1);
h_charge34->SetStats(kFALSE);

TH2D *h_charge7 = new TH2D("h_charge7", Form("T%dZ%d: Yttrium-%s, Side 1: Outer Charge vs Inner Charge", tid, zid, filename.c_str()), qi1OF_bins, qi1OF_range_1, qi1OF_range_2, qo1OF_bins, qo1OF_range_1, qo1OF_range_2);
h_charge7->SetXTitle("Inner Charge [keV]");
h_charge7->SetYTitle("Outer Charge [keV]");
h_charge7->GetYaxis()->SetTitleOffset(1.5);
h_charge7->GetYaxis()->SetLabelSize(0.03);
h_charge7->GetYaxis()->SetTitleSize(0.03);
h_charge7->SetMarkerColor(kRed);
h_charge7->SetFillColor(kRed);
//h_charge7->SetMarkerStyle(29);
//h_charge7->SetMarkerSize(1);
h_charge7->SetStats(kFALSE);

TH2D *h_charge8 = new TH2D("h_charge8", Form("T%dZ%d: Yttrium-%s, Side 2: Outer Charge vs Inner Charge", tid, zid, filename.c_str()), qi2OF_bins, qi2OF_range_1, qi2OF_range_2, qo2OF_bins, qo2OF_range_1, qo2OF_range_2);
h_charge8->SetXTitle("Inner Charge [keV]");
h_charge8->SetYTitle("Outer Charge [keV]");
h_charge8->GetYaxis()->SetTitleOffset(1.5);
h_charge8->GetYaxis()->SetLabelSize(0.03);
h_charge8->GetYaxis()->SetTitleSize(0.03);
h_charge8->SetMarkerColor(kRed);
h_charge8->SetFillColor(kRed);
h_charge8->SetStats(kFALSE);

//Plot Histogram
TCanvas* CanQ= new TCanvas("CanQ", "S1 outer charge Vs inner charge", 600, 600);

chain2->Draw("qo1OF:qi1OF>>h_charge2",cBasicCuts);
chain1->Draw("qo1OF:qi1OF>>+h_charge2",cBasicCuts);
chain2->Draw("qo1OF:qi1OF>>h_charge", cNotEmpty+cRandom+cGoodFlash+cPstd+cBaseTemp+cVoltageBias, "same");
chain1->Draw("qo1OF:qi1OF>>+h_charge", cNotEmpty+cRandom+cGoodFlash+cPstd+cBaseTemp+cVoltageBias, "same");
//cout<<"Bin Number range of Inner Events: "<<h_charge2->GetXaxis()->FindBin(4)<<" to "<<h_charge2->GetXaxis()->FindBin(6)<<endl;
//cout<<"Bin Number range of Inner Noise Events: "<<h_charge->GetXaxis()->FindBin(-1.5)<<" to "<<h_charge->GetXaxis()->FindBin(1.5)<<endl;
//cout<<"Bin Number range of Outer Noise Events: "<<h_charge->GetYaxis()->FindBin(-1.5)<<" to "<<h_charge->GetYaxis()->FindBin(1.5)<<endl;
//break;
h_charge12->Add(h_charge, h_charge2, 1, 1);
//h_charge->Print();
//h_charge2->Print();
//h_charge12->Print();

leg1=new TLegend(0.55,0.6,0.89,0.8);
leg1->SetBorderSize(0);
leg1->AddEntry(h_charge,"Noise spectrum","f");
leg1->AddEntry(h_charge2,"Events spectrum","f");
leg1->Draw();

CanQ->Update();

CanQ->SaveAs(Form("%s/RadialFV/S1FlatRadial_%s_%s_%s.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));

//break;
/*
// Projection of Noise blob and Events Side 1
TH1D *S1proj_inb, *S1proj_onb, *S1proj_inevents;

TCanvas* CanQP1= new TCanvas("CanQP1", "S1 Inner Events Projection", 600, 600);
S1proj_inevents=h_charge2->ProjectionX("S1proj_inevents", h_charge2->GetXaxis()->FindBin(4), h_charge2->GetXaxis()->FindBin(6));
S1proj_inevents->SetTitle(Form("T%dZ%d: Yttrium-%s,%s, Side 1 Inner Event Spectrum Projection", tid, zid, filename.c_str(), weektitle.c_str()));
S1proj_inevents->SetYTitle("counts");
S1proj_inevents->GetYaxis()->SetTitleOffset(1.5);
S1proj_inevents->GetYaxis()->SetLabelSize(0.03);
S1proj_inevents->GetYaxis()->SetTitleSize(0.03);
S1proj_inevents->SetStats(kFALSE);
S1proj_inevents->Draw();
cout<<"------------------ Side 1-------------------"<<endl;
//cout<<"Inner Events Mean: "<<S1proj_inevents->GetMean()<<endl;
//cout<<"Inner Events StdDev: "<<S1proj_inevents->GetStdDev()<<endl;
cout<<"Inner Events Mean+StdDev(A)= "<<S1proj_inevents->GetMean()+S1proj_inevents->GetStdDev()<<endl;
//CanQP1->SaveAs(Form("%s/RadialFV/S1InnerEventsProj_%s_%s_%s.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));

TCanvas* CanQP2= new TCanvas("CanQP2", "S1 Outer Noise Events Projection", 600, 600);
S1proj_onb=h_charge->ProjectionY("S1proj_inevents", h_charge->GetYaxis()->FindBin(-1.5), h_charge->GetYaxis()->FindBin(1.5));
S1proj_onb->SetTitle(Form("T%dZ%d: Yttrium-%s,%s, Side 1 Outer Noise Spectrum Projection", tid, zid, filename.c_str(), weektitle.c_str()));
S1proj_onb->SetYTitle("counts");
S1proj_onb->GetYaxis()->SetTitleOffset(1.5);
S1proj_onb->GetYaxis()->SetLabelSize(0.03);
S1proj_onb->GetYaxis()->SetTitleSize(0.03);
S1proj_onb->SetStats(kFALSE);
S1proj_onb->Draw();
//cout<<"Outer Noise Mean: "<<S1proj_onb->GetMean()<<endl;
//cout<<"Outer Noise StdDev: "<<S1proj_onb->GetStdDev()<<endl;
cout<<"Outer Noise Mean+5*StdDev(D)= "<<S1proj_onb->GetMean()+5*S1proj_onb->GetStdDev()<<endl;
//CanQP2->SaveAs(Form("%s/RadialFV/S1OuterNoiseProj_%s_%s_%s.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));

TCanvas* CanQP3= new TCanvas("CanQP3", "S1 Inner Noise Events Projection", 600, 600);
S1proj_inb=h_charge->ProjectionX("S1proj_inevents", h_charge->GetXaxis()->FindBin(-1.5), h_charge->GetXaxis()->FindBin(1.5));
S1proj_inb->SetTitle(Form("T%dZ%d: Yttrium-%s,%s, Side 1 Outer Noise Spectrum Projection", tid, zid, filename.c_str(), weektitle.c_str()));
//S1proj_onb->SetXTitle("Outer Charge [keV]");
S1proj_inb->SetYTitle("counts");
S1proj_inb->GetYaxis()->SetTitleOffset(1.5);
S1proj_inb->GetYaxis()->SetLabelSize(0.03);
S1proj_inb->GetYaxis()->SetTitleSize(0.03);
S1proj_inb->SetStats(kFALSE);
S1proj_inb->Draw();
//cout<<"Inner Noise Mean: "<<S1proj_inb->GetMean()<<endl;
//cout<<"Inner Noise StdDev: "<<S1proj_inb->GetStdDev()<<endl;
cout<<"Inner Noise Mean+6.5*StdDev(C)= "<<S1proj_inb->GetMean()+6.5*S1proj_inb->GetStdDev()<<endl;
//CanQP3->SaveAs(Form("%s/RadialFV/S1InnerNoiseProj_%s_%s_%s.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));
//break;
*/

TCanvas* CanQ2= new TCanvas("CanQ2", "S2 outer charge Vs inner charge", 600, 600);

chain2->Draw("qo2OF:qi2OF>>h_charge4",cBasicCuts);
chain1->Draw("qo2OF:qi2OF>>+h_charge4",cBasicCuts);
chain2->Draw("qo2OF:qi2OF>>h_charge3", cNotEmpty+cRandom+cGoodFlash+cPstd+cBaseTemp+cVoltageBias, "same");
chain1->Draw("qo2OF:qi2OF>>+h_charge3", cNotEmpty+cRandom+cGoodFlash+cPstd+cBaseTemp+cVoltageBias, "same");
//h_charge34->Add(h_charge3, h_charge4, 1, 1);
//h_charge3->Print();
//h_charge4->Print();

leg2=new TLegend(0.55,0.6,0.89,0.8);
leg2->SetBorderSize(0);
leg2->AddEntry(h_charge3,"Noise spectrum","f");
leg2->AddEntry(h_charge4,"Events spectrum","f");
leg2->Draw();

CanQ2->Update();

CanQ2->SaveAs(Form("%s/RadialFV/S2FlatRadial_%s_%s_%s.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));
//break;
/*
// Projection of Noise blob and Events Side 2
TH1D *S2proj_inb, *S2proj_onb, *S2proj_inevents;

TCanvas* CanQP4= new TCanvas("CanQP4", "S2 Inner Events Projection", 600, 600);
S2proj_inevents=h_charge4->ProjectionX("S2proj_inevents", h_charge2->GetXaxis()->FindBin(4), h_charge2->GetXaxis()->FindBin(6));
S2proj_inevents->SetTitle(Form("T%dZ%d: Yttrium-%s,%s, Side 2 Inner Event Spectrum Projection", tid, zid, filename.c_str(), weektitle.c_str()));
S2proj_inevents->SetYTitle("counts");
S2proj_inevents->GetYaxis()->SetTitleOffset(1.5);
S2proj_inevents->GetYaxis()->SetLabelSize(0.03);
S2proj_inevents->GetYaxis()->SetTitleSize(0.03);
S2proj_inevents->SetStats(kFALSE);
S2proj_inevents->Draw();
cout<<"------------------ Side 2-------------------"<<endl;
//cout<<"Inner Events Mean: "<<S2proj_inevents->GetMean()<<endl;
//cout<<"Inner Events StdDev: "<<S2proj_inevents->GetStdDev()<<endl;
cout<<"Inner Events Mean+StdDev (A)= "<<S2proj_inevents->GetMean()+S2proj_inevents->GetStdDev()<<endl;
//CanQP4->SaveAs(Form("%s/RadialFV/S2InnerEventsProj_%s_%s_%s.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));

TCanvas* CanQP5= new TCanvas("CanQP5", "S2 Outer Noise Events Projection", 600, 600);
S2proj_onb=h_charge3->ProjectionY("S2proj_inevents", h_charge->GetYaxis()->FindBin(-1), h_charge->GetYaxis()->FindBin(1));
S2proj_onb->SetTitle(Form("T%dZ%d: Yttrium-%s,%s, Side 2 Outer Noise Spectrum Projection", tid, zid, filename.c_str(), weektitle.c_str()));
S2proj_onb->SetYTitle("counts");
S2proj_onb->GetYaxis()->SetTitleOffset(1.5);
S2proj_onb->GetYaxis()->SetLabelSize(0.03);
S2proj_onb->GetYaxis()->SetTitleSize(0.03);
S2proj_onb->SetStats(kFALSE);
S2proj_onb->Draw();
//cout<<"Outer Noise Mean: "<<S2proj_onb->GetMean()<<endl;
//cout<<"Outer Noise StdDev: "<<S2proj_onb->GetStdDev()<<endl;
cout<<"Outer Noise Mean+4.5*StdDev (D)= "<<S2proj_onb->GetMean()+4.5*S2proj_onb->GetStdDev()<<endl;
//CanQP5->SaveAs(Form("%s/RadialFV/S2OuterNoiseProj_%s_%s_%s.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));

TCanvas* CanQP6= new TCanvas("CanQP6", "S2 Inner Noise Events Projection", 600, 600);
S2proj_inb=h_charge3->ProjectionX("S2proj_inevents", h_charge->GetXaxis()->FindBin(-2.5), h_charge->GetXaxis()->FindBin(2.5));
S2proj_inb->SetTitle(Form("T%dZ%d: Yttrium-%s,%s, Side 2 Outer Noise Spectrum Projection", tid, zid, filename.c_str(), weektitle.c_str()));
S2proj_inb->SetYTitle("counts");
S2proj_inb->GetYaxis()->SetTitleOffset(1.5);
S2proj_inb->GetYaxis()->SetLabelSize(0.03);
S2proj_inb->GetYaxis()->SetTitleSize(0.03);
S2proj_inb->SetStats(kFALSE);
S2proj_inb->Draw();
//cout<<"Inner Noise Mean: "<<S2proj_inb->GetMean()<<endl;
//cout<<"Inner Noise StdDev: "<<S2proj_inb->GetStdDev()<<endl;
cout<<"Inner Noise Mean+2.5*StdDev (C)= "<<S2proj_inb->GetMean()+2.5*S2proj_inb->GetStdDev()<<endl;
//CanQP6->SaveAs(Form("%s/RadialFV/S2InnerNoiseProj_%s_%s_%s.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));
//break;
*/

// ------------------------Effect of Radial Cut---------------------------------

//For T2Z1 Side 1 YBE: -A=inevnt_mu+1.5*inevnt_sigma , B= 4.2, C=inb_mu+7*inb_sig, D=onb_mu+5*onb_sig
//TCut cRadialCutS1("qo1OF<(-0.8/(1+TMath::Exp(-4.2*(qi1OF-2.6))) + 1.6)");
TCut cRadialCutS1("qo1OF<1.6");
//For T2Z1 Side 2 YBE: -A=inevnt_mu+0.5*inevnt_sigma , B= 4.2, C=inb_mu+2.5*inb_sig, D=onb_mu+4.5*onb_sig
//TCut cRadialCutS2("qo2OF<(-0.7/(1+TMath::Exp(-4.2*(qi2OF-2.6))) + 1.3)");
TCut cRadialCutS2("qo2OF<1.3");
//For T5Z2 Side 1 YBE: -A=inevnt_mu+0.3*inevnt_sigma , B= 4.2, C=inb_mu+4*inb_sig, D=onb_mu+4*onb_sig
//TCut cRadialCutS1("qo1OF<(-1.8/(1+TMath::Exp(-4.2*(qi1OF-2.9))) + 3)");
//TCut cRadialCutS1("qo1OF<3")
//For T5Z2 Side 2 YBE: -A=inevnt_mu+0.5*inevnt_sigma , B= 4.2, C=inb_mu+4*inb_sig, D=onb_mu+4*onb_sig
//TCut cRadialCutS2("qo2OF<(-0.9/(1+TMath::Exp(-4.2*(qi2OF-1.9))) + 1.8)");
//TCut cRadialCutS1("qo2OF<1.8")



TCanvas* CanQ4= new TCanvas("CanQ4", "S1 outer charge Vs inner charge", 600, 600);

//h_charge12->Draw();
//chain2->Draw("qo1OF:qi1OF>>h_charge7", cBasicCuts+cRadialCutS1, "same");
//chain2->Draw("qo1OF:qi1OF>>+h_charge7", cNotEmpty+cRandom+cGoodFlash+cPstd+cBaseTemp+cVoltageBias+cRadialCutS1, "same");
//h_charge7->Print();
h_charge2->Draw();
h_charge->Draw("same");
leg4=new TLegend(0.55,0.6,0.89,0.85);
leg4->SetBorderSize(0);
//leg4->AddEntry(h_charge12,"without radial cut","f");
//leg4->AddEntry(h_charge7,"with radial cut","f");
//leg4->Draw("same");

//Draw Radial Cut
line1=new TLine(qi1OF_range_1,1.6,qi1OF_range_2,1.6);
//line1=new TLine(qi1OF_range_1,1.6,2.6,1.6);
line1->SetLineWidth(2);
line1->SetLineColor(kGreen);
line1->Draw("same");
/*
line1a=new TLine(2.6,1.6,2.6,qo1OF_range_2);
line1a->SetLineWidth(2);
line1a->SetLineColor(kGreen);
line1a->Draw("same");
*/

//Draw TGraph
gr1 = new TGraph(n,Qi,Qo);
gr1->SetName("gr1");
gr1->SetLineColor(kBlue);
gr1->SetLineWidth(2);
gr1->Draw("C same");

leg4->AddEntry(h_charge,"Noise spectrum","f");
leg4->AddEntry(h_charge2,"Events spectrum","f");
leg4->AddEntry(line1,"Radial Cut","l");
leg4->AddEntry(gr1,"Neutron SS edge","l");
leg4->Draw("same");


CanQ4->Update();

CanQ4->SaveAs(Form("%s/RadialFV/S1LooseRadialCut_%s_%s_%s.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));
//break;

TCanvas* CanQ5= new TCanvas("CanQ5", "S2 outer charge Vs inner charge", 600, 600);

//h_charge34->Draw();
//chain2->Draw("qo2OF:qi2OF>>h_charge8", cBasicCuts+cRadialCutS2, "same");
//chain2->Draw("qo2OF:qi2OF>>+h_charge8", cNotEmpty+cRandom+cGoodFlash+cPstd+cBaseTemp+cVoltageBias+cRadialCutS2, "same");
//h_charge8->Print();
h_charge4->Draw();
h_charge3->Draw("same");
leg5=new TLegend(0.55,0.6,0.89,0.8);
leg5->SetBorderSize(0);
//leg5->AddEntry(h_charge34,"without radial cut","f");
//leg5->AddEntry(h_charge8,"with radial cut","f");
//leg5->Draw("same");

//Draw Radial Cut
line2=new TLine(qi2OF_range_1,1.3,qi2OF_range_2,1.3);
//line2=new TLine(qi2OF_range_1,1.3,2.6,1.3);
line2->SetLineWidth(2);
line2->SetLineColor(kGreen);
line2->Draw("same");
/*
line2a=new TLine(2.6,1.3,2.6,qo2OF_range_2);
line2a->SetLineWidth(2);
line2a->SetLineColor(kGreen);
line2a->Draw("same");
*/

//Draw TGraph
gr2 = new TGraph(n2,Qi2,Qo2);
gr2->SetName("gr2");
gr2->SetLineColor(kBlue);
gr2->SetLineWidth(2);
gr2->Draw("C same");

leg5->AddEntry(h_charge3,"Noise spectrum","f");
leg5->AddEntry(h_charge4,"Events spectrum","f");
leg5->AddEntry(line2,"Radial Cut","l");
leg5->AddEntry(gr2,"Neutron SS edge","l");
leg5->Draw("same");


CanQ5->Update();

CanQ5->SaveAs(Form("%s/RadialFV/S2LooseRadialCut_%s_%s_%s.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));
//break;
}
