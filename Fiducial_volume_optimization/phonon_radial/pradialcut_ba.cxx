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



void pradialcut_ba(string weekno="1", Int_t zip=4)
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
string opfileaddress="/home/u1/rik/Viraj/FV/phononradial";

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

TCut virajPhononRadial("prpartOF<(TMath::Exp(-0.424675-3.53768*ptNF)+0.295494)"); // 2 sigma
TCut virajchargeradialcut("qo1OF<1.05087 && qo2OF<1.23711");

 cout<<"\nNumber of events: "<<chain2->GetEntries()<<endl;
cout<<"Name of Tree: "<<chain2->GetName()<<endl;

//chain2->Scan("SeriesNumber","SeriesNumber==11509070319", "colsize=15 precision=13");


double ptNFbinsize=150;
double ptNFrange_1a=0;
double ptNFrange_2a=15;

double prpartOF_bins=200;
double prpartOF_range1=0.1;
double prpartOF_range2=0.5;

gStyle->SetOptFit(1);


TH2D *h_phonon = new TH2D("h_phonon", Form("Z%d: Barium-%s,prpart vs ptNF",zid, filename.c_str()), ptNFbinsize, ptNFrange_1a, ptNFrange_2a, 
			  prpartOF_bins, prpartOF_range1, prpartOF_range2);
h_phonon->SetXTitle("ptNF [keV]");
h_phonon->SetYTitle("prpartOF");
h_phonon->GetYaxis()->SetTitleOffset(1.5);
h_phonon->GetYaxis()->SetLabelSize(0.03);
h_phonon->GetYaxis()->SetTitleSize(0.03);
h_phonon->GetYaxis()->CenterTitle(true);
h_phonon->GetXaxis()->CenterTitle(true);
h_phonon->SetMarkerColor(kRed);
h_phonon->SetFillColor(kRed);
//h_phonon->SetStats(kFALSE);

TCanvas* canvas_th2dphonon= new TCanvas("canvas_th2dphonon", "Prpart Vs ptNF canvas", 600, 600);
chain2->Draw("prpartOF:ptNF>>h_phonon",cBasicCuts);

cout<<" Bin 1st: "<<h_phonon->GetXaxis()->FindBin(0)<<"  Bin last: "<< h_phonon->GetXaxis()->FindBin(1)<<endl;
//Projection of Y 
TF1 *gausfit=new TF1("gausfit","gaus", prpartOF_range1, prpartOF_range2); 

TF1 *polyfit=new TF1("polyfit","TMath::Exp(-[0]-[1]*x)+[2]", ptNFrange_1a, ptNFrange_2a);
polyfit->SetLineColor(1); // 3 sigma
polyfit->SetParameters(-0.8,5,0.2); // For 3 and 2.5 sigma fits



double binw=(ptNFrange_2a-ptNFrange_1a)/ptNFbinsize;
//cout<<endl<<"binw of ptNF is: "<<binw<<endl;
//TCanvas* CanprojY= new TCanvas("CanprojY", "Projection Y", 600, 600);

/*
TH1D *hist_projY;
//hist_projY=h_phonon->ProjectionY("hist_projY", h_phonon->GetXaxis()->FindBin(0), h_phonon->GetXaxis()->FindBin(1));
hist_projY=h_phonon->ProjectionY("hist_projY", h_phonon->GetXaxis()->FindBin(1), h_phonon->GetXaxis()->FindBin(1));
hist_projY->SetYTitle("counts");
hist_projY->GetYaxis()->SetTitleOffset(1.5);
hist_projY->Rebin(5);
hist_projY->Fit("gausfit");
TF1 *f = hist_projY->GetListOfFunctions()->FindObject("gausfit");
if(f){
    f->SetLineColor(3); //3=green 4 = "blue" 
    f->SetLineWidth(4.5); //f->SetLineStyle(1); // 2 = "- - -" 
    }


//hist_projY->SetAxisRange(0.1, 0.2,"X");
hist_projY->Draw();
*/
//Define histogram for mean of psumo vs ptnf side 1
TH2D *h_psumptnfmeanS1 = new TH2D("h_psumptnfmeanS1", Form("T%dZ%d: Barium, %s, prpartOF vs ptNF", tid, zid, weektitle.c_str()), ptNFbinsize, ptNFrange_1a, ptNFrange_2a, prpartOF_bins, prpartOF_range1, prpartOF_range2);
h_psumptnfmeanS1->SetXTitle("ptNF [keV]");
h_psumptnfmeanS1->SetYTitle("prpartOF");
h_psumptnfmeanS1->GetYaxis()->SetTitleOffset(1.5);
h_psumptnfmeanS1->GetYaxis()->SetLabelSize(0.03);
h_psumptnfmeanS1->GetYaxis()->SetTitleSize(0.03);
h_psumptnfmeanS1->SetMarkerColor(5);
h_psumptnfmeanS1->SetMarkerSize(1);
h_psumptnfmeanS1->SetMarkerStyle(20);
//h_psumptnfmeanS1->SetStats(kFALSE);

//Define histogram for 1 sigma of psumo vs ptnf side 1
TH2D *h_psumptnfsigmaS1 = new TH2D("h_psumptnfsigmaS1", Form("T%dZ%d: Barium, %s, prpartOF vs ptNF", tid, zid, weektitle.c_str()), ptNFbinsize, ptNFrange_1a, ptNFrange_2a, prpartOF_bins, prpartOF_range1, prpartOF_range2);
h_psumptnfsigmaS1->SetXTitle("ptNF [keV]");
h_psumptnfsigmaS1->SetYTitle("prpartOF");
h_psumptnfsigmaS1->GetYaxis()->SetTitleOffset(1.5);
h_psumptnfsigmaS1->GetYaxis()->SetLabelSize(0.03);
h_psumptnfsigmaS1->GetYaxis()->SetTitleSize(0.03);
h_psumptnfsigmaS1->SetMarkerColor(4);
h_psumptnfsigmaS1->SetMarkerSize(3);
h_psumptnfsigmaS1->SetMarkerStyle(3);
//h_psumptnfsigmaS1->SetStats(kFALSE);

const int nbins=20; //Previously 222

TH1D *py[nbins];
TCanvas* poptS1[nbins];
TLatex *t3[nbins];


double var1_ptNF, var2_ptNF; 
for(int i=0; i<nbins; i++)
{cout<<"----------------- Bin No: "<<i<<"------------------"<<endl;
poptS1[i]= new TCanvas(Form("poptS1_%d", i), Form("prpartOF bin%d",i), 600, 600);
py[i]=h_phonon->ProjectionY(Form("py%d", i), i, i);
var1_ptNF=i*binw; var2_ptNF=var1_ptNF+binw; cout<<"var1_ptNF = "<<var1_ptNF<<"var2_ptNF = "<<var2_ptNF<<endl;
py[i]->Rebin(5);
py[i]->SetTitle(Form("T%dZ%d: Barium, %s, prpartOF, Bin No. %d, %.2f <ptNF < %.2f", tid, zid, weektitle.c_str(), i,var1_ptNF, var2_ptNF ));
py[i]->SetYTitle("counts");
py[i]->GetYaxis()->SetTitleOffset(1.5);
py[i]->GetYaxis()->SetLabelSize(0.03);
py[i]->GetYaxis()->SetTitleSize(0.03);
//py[i]->SetStats(kFALSE);
poptS1[i]->cd();
py[i]->Draw();

py[i]->Fit("gausfit");
TF1 *f = py[i]->GetListOfFunctions()->FindObject("gausfit");
if(f){
    f->SetLineColor(3); //3=green 4 = "blue" 
    f->SetLineWidth(4.5); //f->SetLineStyle(1); // 2 = "- - -" 
    }
//t3[i]=new TLatex(.55,.8,Form("#scale[0.7]{%.2f #leq ptNF #leq %.2f}", (nbins-(nbins-(i-1)))*binw, i*binw));
//t3[i]->SetNDC(kTRUE);
//t3[i]->Draw("same");
//poptS1[i]->SaveAs(Form("%s/New/Bins/prptbin%d_%s_%s_%s.gif", opfileaddress.c_str(), i, fileaddress.c_str(),det.c_str(),weektitle.c_str()));

h_psumptnfmeanS1->Fill(((i*binw)-(binw/2)), gausfit->GetParameter(1));
h_psumptnfsigmaS1->Fill(((i*binw)-(binw/2)), (gausfit->GetParameter(1)+2*gausfit->GetParameter(2)));
//h_psumptnf3sigmaS1->Fill(((i*binw)-(binw/2)-10), (gausfit->GetParameter(1)+(2*gausfit->GetParameter(2))));

}

TCanvas* canvas_mean_sigma= new TCanvas("canvas_mean_sigma", "Side 1 All Means and sigmas histogram", 600, 600);
h_phonon->Draw();
//h_psumptnfmeanS1->Draw("same");
h_psumptnfsigmaS1->Draw("same");
h_psumptnfsigmaS1->Fit("polyfit", "R");

   TLegend *leg = new TLegend(0.35,0.7,0.45,0.8,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->AddEntry("h_phonon","Data","lpf");
   leg->AddEntry("h_psumptnfsigmaS1","mean+sigma of ProjectionY","lpf");
   leg->AddEntry("polyfit","Exponential Fit","lpf");
   leg->Draw("same");
   
  // TCanvas* canvas_aftercut= new TCanvas("canvas_after cut", "After phonon radial cut", 600, 600);
//chain2->Draw("prpartOF:ptNF",cBasicCuts && virajPhononRadial && "ptNF<15" && "ptNF>0" && "prpartOF>0.1" && "prpartOF<0.5");
chain2->Draw("0.5*(qi1OF+qi2OF):precoiltNF",cBasicCuts && virajPhononRadial && "qi1OF+qi2OF>0" && "0.5*(qi1OF+qi2OF)<160 " && "precoiltNF>0" && "precoiltNF<160");

   
}
