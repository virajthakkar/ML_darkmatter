/*
THIS CODE IS TO plot RRQ plots USING Ba DATA.
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

void ERNR_cf_PRCba2(string weekno="1", Int_t zip=4)
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
string opfileaddress=Form("/home/u1/rik/Viraj/ba_bkg/mf%s", det.c_str());

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

Double_t ylinefun(double x) 
{ m1=0.95; c2=-10;
    return m1*x + c2 ; }
//TCut threesigmacut("( (qi1OF+qi2OF)/2.0 ) > (70/170)*precoiltNF -10 && ((qi1OF+qi2OF)/2.0) >3 ");  // simulation for 3 sigma NR band cut with ptNF on x axis
TCut threesigmacut("( (qi1OF+qi2OF)/2.0 ) > 0.95*precoiltNF -10 && ((qi1OF+qi2OF)/2.0) >3 ");  // simulation for 3 sigma NR band cut with ptNF on x axis
TCut cQualityCuts=cLEChisq+cDeltaChisqGlitch+cDeltaChiSqLFN+cChargeChiSq+cSymmetricFV;





//Define Histogram for plotting Avg Charge Energy vs ptNF with Basic Cuts
/*TH2D *h_ERNRallcuts = new TH2D("h_ERNRallcuts", Form("T%dZ%d: %s-%s, %s, Ionization vs precoiltNF", tid, zid, sourcefullname.c_str(), filename.c_str(), weektitle.c_str()), ptNFbinsize, ptNFrange_1, ptNFrange_2, qbinsize, q_y_range_1, q_y_range_2);
h_ERNRallcuts->SetXTitle("precoiltNF [keV]");
h_ERNRallcuts->SetYTitle("qimean ");
h_ERNRallcuts->GetYaxis()->SetTitleOffset(1.5);
h_ERNRallcuts->GetYaxis()->SetLabelSize(0.03);
h_ERNRallcuts->GetYaxis()->SetTitleSize(0.03);
h_ERNRallcuts->SetLineColor(kGray+3);
h_ERNRallcuts->SetFillColor(kGray+3);
h_ERNRallcuts->SetMarkerColor(kGray+3);
//h_ERNRallcuts->SetMaximum(90);
//h_ERNRallcuts->SetMinimum(1);
h_ERNRallcuts->SetStats(kFALSE);
*/


double precoiltNFrange_1=10; 
double precoiltNFrange_2=170; 
double precoiltNFbinsize=8; 
TH1D *h_precoiltNF = new TH1D("h_precoiltNF", Form("T%dZ%d: %s-%s, %s,  precoiltNF", tid, zid, sourcefullname.c_str(), filename.c_str(), weektitle.c_str()), precoiltNFbinsize, precoiltNFrange_1, precoiltNFrange_2);
h_precoiltNF->SetXTitle("precoiltNF [keV]");
h_precoiltNF->SetYTitle("Counts");
h_precoiltNF->GetYaxis()->SetTitleOffset(1.5);
h_precoiltNF->GetYaxis()->SetLabelSize(0.03);
h_precoiltNF->GetYaxis()->SetTitleSize(0.03);
//h_precoiltNF->SetLineColor(kGray+3);
//h_precoiltNF->SetFillColor(kGray+3);
//h_precoiltNF->SetMarkerColor(kGray+3);
//h_precoiltNF->SetMaximum(90);
//h_precoiltNF->SetMinimum(1);
//h_precoiltNF->SetStats(kFALSE);
TCanvas* c_precoiltNF= new TCanvas("c_precoiltNF", "c_precoiltNF (All Cuts) 1", 600, 600);
chain2->Draw("precoiltNF>>h_precoiltNF", cBasicCuts+cQualityCuts+threesigmacut);
h_precoiltNF->Print();
//cout<<"h_precoiltNF Integral :"<<h_precoiltNF->Integral()<<endl;

TLegend *leg=new TLegend(0.55,0.57,0.90,0.77);
leg->SetBorderSize(0);
leg->AddEntry(h_precoiltNF,"after Quality+PRC cuts+ 3sigma NR cut","lpf");
leg->Draw("same");

c_precoiltNF->Update();
c_precoiltNF->SaveAs(Form("%s/c_precoiltNF_%s_%s_%s_zout.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));


double qzpartOFrange_1=-0.3; 
double qzpartOFrange_2=0.3; 
double qzpartOFbinsize=6; 
TH1D *h_qzpartOF = new TH1D("h_qzpartOF", Form("T%dZ%d: %s-%s, %s,  qzpartOF", tid, zid, sourcefullname.c_str(), filename.c_str(), weektitle.c_str()), qzpartOFbinsize, qzpartOFrange_1, qzpartOFrange_2);
h_qzpartOF->SetXTitle("qzpartOF");
h_qzpartOF->SetYTitle("Counts");
h_qzpartOF->GetYaxis()->SetTitleOffset(1.5);
h_qzpartOF->GetYaxis()->SetLabelSize(0.03);
h_qzpartOF->GetYaxis()->SetTitleSize(0.03);

TCanvas* c_qzpartOF= new TCanvas("c_qzpartOF", "c_qzpartOF (All Cuts) 1", 600, 600);
chain2->Draw("qzpartOF>>h_qzpartOF", cBasicCuts+cQualityCuts+threesigmacut);
h_qzpartOF->Print();

TLegend *leg2=new TLegend(0.55,0.57,0.90,0.77);
leg2->SetBorderSize(0);
leg2->AddEntry(h_qzpartOF,"after Quality+PRC cuts+ 3sigma NR cut","lpf");
leg2->Draw("same");

c_qzpartOF->Update();
c_qzpartOF->SaveAs(Form("%s/c_qzpartOF_%s_%s_%s_zout.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));


double qrpart1OFrange_1=-0.2; 
double qrpart1OFrange_2=0.6; 
double qrpart1OFbinsize=8; 
TH1D *h_qrpart1OF = new TH1D("h_qrpart1OF", Form("T%dZ%d: %s-%s, %s,  qrpart1OF", tid, zid, sourcefullname.c_str(), filename.c_str(), weektitle.c_str()), qrpart1OFbinsize, qrpart1OFrange_1, qrpart1OFrange_2);
h_qrpart1OF->SetXTitle("qrpart1OF");
h_qrpart1OF->SetYTitle("Counts");
h_qrpart1OF->GetYaxis()->SetTitleOffset(1.5);
h_qrpart1OF->GetYaxis()->SetLabelSize(0.03);
h_qrpart1OF->GetYaxis()->SetTitleSize(0.03);
//h_qrpart1OF->SetLineColor(kGray+3);
//h_qrpart1OF->SetFillColor(kGray+3);
//h_qrpart1OF->SetMarkerColor(kGray+3);
//h_qrpart1OF->SetMaximum(90);
//h_qrpart1OF->SetMinimum(1);
//h_qrpart1OF->SetStats(kFALSE);
TCanvas* c_qrpart1OF= new TCanvas("c_qrpart1OF", "c_qrpart1OF (All Cuts) 1", 600, 600);
chain2->Draw("qrpart1OF>>h_qrpart1OF", cBasicCuts+cQualityCuts+threesigmacut);
h_qrpart1OF->Print();

TLegend *leg3=new TLegend(0.55,0.57,0.90,0.77);
leg3->SetBorderSize(0);
leg3->AddEntry(h_qrpart1OF,"after Quality+PRC cuts+ 3sigma NR cut","lpf");
leg3->Draw("same");

c_qrpart1OF->Update();
c_qrpart1OF->SaveAs(Form("%s/c_qrpart1OF_%s_%s_%s_zout.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));

double ytNFrange_1=0.6; 
double ytNFrange_2=1.4; 
double ytNFbinsize=8; 
TH1D *h_ytNF = new TH1D("h_ytNF", Form("T%dZ%d: %s-%s, %s,  ytNF", tid, zid, sourcefullname.c_str(), filename.c_str(), weektitle.c_str()), ytNFbinsize, ytNFrange_1, ytNFrange_2);
h_ytNF->SetXTitle("ytNF");
h_ytNF->SetYTitle("Counts");
h_ytNF->GetYaxis()->SetTitleOffset(1.5);
h_ytNF->GetYaxis()->SetLabelSize(0.03);
h_ytNF->GetYaxis()->SetTitleSize(0.03);

TCanvas* c_ytNF= new TCanvas("c_ytNF", "c_ytNF (All Cuts) 1", 600, 600);
chain2->Draw("ytNF>>h_ytNF", cBasicCuts+cQualityCuts+threesigmacut);
h_ytNF->Print();

TLegend *leg4=new TLegend(0.55,0.57,0.90,0.77);
leg4->SetBorderSize(0);
leg4->AddEntry(h_ytNF,"after Quality+PRC cuts+ 3sigma NR cut","lpf");
leg4->Draw("same");

c_ytNF->Update();
c_ytNF->SaveAs(Form("%s/c_ytNF_%s_%s_%s_zout.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));


double prpartOFrange_1=0.1; 
double prpartOFrange_2=0.3; 
double prpartOFbinsize=10; 
TH1D *h_prpartOF = new TH1D("h_prpartOF", Form("T%dZ%d: %s-%s, %s,  prpartOF", tid, zid, sourcefullname.c_str(), filename.c_str(), weektitle.c_str()), prpartOFbinsize, prpartOFrange_1, prpartOFrange_2);
h_prpartOF->SetXTitle("prpartOF");
h_prpartOF->SetYTitle("Counts");
h_prpartOF->GetYaxis()->SetTitleOffset(1.5);
h_prpartOF->GetYaxis()->SetLabelSize(0.03);
h_prpartOF->GetYaxis()->SetTitleSize(0.03);

TCanvas* c_prpartOF= new TCanvas("c_prpartOF", "c_prpartOF (All Cuts) 1", 600, 600);
chain2->Draw("prpartOF>>h_prpartOF", cBasicCuts+cQualityCuts+threesigmacut);
h_prpartOF->Print();

TLegend *leg5=new TLegend(0.55,0.57,0.90,0.77);
leg5->SetBorderSize(0);
leg5->AddEntry(h_prpartOF,"after Quality+PRC cuts+ 3sigma NR cut","lpf");
leg5->Draw("same");

c_prpartOF->Update();
c_prpartOF->SaveAs(Form("%s/c_prpartOF_%s_%s_%s_zout.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));

double pzpartOFrange_1=-0.1; 
double pzpartOFrange_2=0.2; 
double pzpartOFbinsize=15; 
TH1D *h_pzpartOF = new TH1D("h_pzpartOF", Form("T%dZ%d: %s-%s, %s,  pzpartOF", tid, zid, sourcefullname.c_str(), filename.c_str(), weektitle.c_str()), pzpartOFbinsize, pzpartOFrange_1, pzpartOFrange_2);
h_pzpartOF->SetXTitle("pzpartOF");
h_pzpartOF->SetYTitle("Counts");
h_pzpartOF->GetYaxis()->SetTitleOffset(1.5);
h_pzpartOF->GetYaxis()->SetLabelSize(0.03);
h_pzpartOF->GetYaxis()->SetTitleSize(0.03);

TCanvas* c_pzpartOF= new TCanvas("c_pzpartOF", "c_pzpartOF (All Cuts) 1", 600, 600);
chain2->Draw("pzpartOF>>h_pzpartOF", cBasicCuts+cQualityCuts+threesigmacut);
h_pzpartOF->Print();

TLegend *leg6=new TLegend(0.55,0.57,0.90,0.77);
leg6->SetBorderSize(0);
leg6->AddEntry(h_pzpartOF,"after Quality+PRC cuts+ 3sigma NR cut","lpf");
leg6->Draw("same");

c_pzpartOF->Update();
c_pzpartOF->SaveAs(Form("%s/c_pzpartOF_%s_%s_%s_zout.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));


double qsummaxOFrange_1=0; 
double qsummaxOFrange_2=150; 
double qsummaxOFbinsize=30; 
TH1D *h_qsummaxOF = new TH1D("h_qsummaxOF", Form("T%dZ%d: %s-%s, %s,  qsummaxOF", tid, zid, sourcefullname.c_str(), filename.c_str(), weektitle.c_str()), qsummaxOFbinsize, qsummaxOFrange_1, qsummaxOFrange_2);
h_qsummaxOF->SetXTitle("qsummaxOF");
h_qsummaxOF->SetYTitle("Counts");
h_qsummaxOF->GetYaxis()->SetTitleOffset(1.5);
h_qsummaxOF->GetYaxis()->SetLabelSize(0.03);
h_qsummaxOF->GetYaxis()->SetTitleSize(0.03);

TCanvas* c_qsummaxOF= new TCanvas("c_qsummaxOF", "c_qsummaxOF (All Cuts) 1", 600, 600);
chain2->Draw("qsummaxOF>>h_qsummaxOF", cBasicCuts+cQualityCuts+threesigmacut);
h_qsummaxOF->Print();

TLegend *leg7=new TLegend(0.55,0.27,0.90,0.37);
leg7->SetBorderSize(0);
leg7->AddEntry(h_qsummaxOF,"after Quality+PRC cuts+ 3sigma NR cut","lpf");
leg7->Draw("same");

c_qsummaxOF->Update();
c_qsummaxOF->SaveAs(Form("%s/c_qsummaxOF_%s_%s_%s_zout.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));



TFile *outfile = new TFile(Form("%s/Final_plots_%s_%s_%s.root", opfileaddress.c_str(), det.c_str(),weektitle.c_str(),filename.c_str()), "recreate");
h_precoiltNF->Write();
h_qzpartOF->Write();
h_qrpart1OF->Write();
h_ytNF->Write();
h_prpartOF->Write();
h_pzpartOF->Write();
h_qsummaxOF->Write();
outfile->Close();


/*TCanvas* qptNF= new TCanvas("qptNF", "Ionization vs precoiltNF (All Cuts) 1", 600, 600);
chain2->Draw("((qi1OF+qi2OF)/2):precoiltNF>>h_ERNRallcuts", cBasicCuts+cQualityCuts, "colz");
const double m= 0.95; //slope 
const double c=-10; // y intecept

  //TLine *line1 = new TLine(0,-10,170,60);
  TLine *line1 = new TLine((3-c)/m,3,80,ylinefun(80));
    line1->SetLineWidth(4);
    line1->SetLineColor(2);
  line1->Draw("same");
TLine *line2 = new TLine(0,3,(3-c)/m,3); // y=3 line
    line2->SetLineWidth(4);
    line2->SetLineColor(2);
  line2->Draw("same");
  //chain2->Draw("((qi1OF+qi2OF)/2):ptNF>>h_ERNRallcuts2", cBasicCuts+cQualityCuts+cPhononRadial+cRadialCut, "colz");
h_ERNRallcuts->Print();
//h_ERNRallcuts2->Print();

TLegend *leg2=new TLegend(0.15,0.60,0.50,0.80);
leg2->SetBorderSize(0);
leg2->AddEntry(h_ERNRallcuts,"after Quality+PRC cuts","f");
//leg->AddEntry(h_ERNRallcuts2,"Failing PRC","f");
//leg->Draw("same");

qptNF->Update();

// Save the output plot in .gif form
qptNF->SaveAs(Form("%s/qiprecoiltNF_%s_%s_%s.gif", opfileaddress.c_str(), fileaddress.c_str(),det.c_str(),weektitle.c_str()));

*/

   
}
