#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TMinuit.h>
#include <TFile.h>
#include <TLegend.h>
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include <TRandom2.h>
#include <TRandom3.h>
#include <TRandom.h>
#include <TFitResultPtr.h>
#include <RConfig.h>
#include <TStopwatch.h>


#define a1 2.00e-2 //Linear conversion factors from time to ADC counts. Produced by timecal.C
#define da1 2.00e-4
#define of1 -3.2e-1
#define dof1 2e-2
#define a2 2.28e-2
#define da2 3.00e-4
#define of2 -3.5e-1
#define dof2 2e-2
#define a3 -1.44e-1 //Linear conversion factors from position to time. Produced by positions.C
#define da3 2e-3
#define of3 4.89e1
#define dof3 2e-1
#define c_b 14
#define sp_res 0.85
#define Zu 175
#define Dcx 140
#define Dcy 2

TRandom3* Rx = new TRandom3(23);
//Velocità nella luce nel vuoto
Double_t C = 29.9792458 ;
//Dimensioni fisichedel sistema
//Velocità di propagazione della luce nella sbarra
Double_t beta_s = 0.4697 * C;

//Lunghezze in metri
//Lunghezze lastra sopra
Double_t Ux = 279;
Double_t Uy = 4;

//Lunghezze lastra sotto
Double_t Dx = 14;
Double_t Dy = 23;


//Qulache variabile per la simulazione
//Massa del muone in GeV
Double_t M_mu = 0.105;

//Intervallo di energia della simulazione in GeV
Double_t Emin = M_mu+ 0.04;
Double_t Emax = 10 ; 

//Prametri per il fit
Double_t BIN_SIZE = 1e-1;
Double_t TIN = 0e;
Double_t TFIN = 20e;

Double_t pdf(Double_t* x, Double_t *par)
{
Double_t y=0;
Double_t  a=par[0];
Double_t  s1=par[1];
Double_t  s2=par[2];
Double_t  offset=par[3];
Double_t  Norm=par[4];

Double_t s1_eff=sqrt(s1*s1+s1*s1+s2*s2);
Double_t s3_eff=s2;
Double_t s2_eff=sqrt(s1*s1+s2*s2);

Double_t g1=1/(sqrt(2*TMath::Pi())*s1_eff)*exp(-(pow((x[0]-offset)/(2*s1_eff),2.0)));
Double_t g2=1/(sqrt(2*TMath::Pi())*s2_eff)*exp(-(pow((x[0]-offset)/(2*s2_eff),2.0)));
Double_t g3=1/(sqrt(2*TMath::Pi())*s3_eff)*exp(-(pow((x[0]-offset)/(2*s3_eff),2.0)));
y= Norm*BIN_SIZE*(a*a*g1+ 2*a*(1-a)*g2 +(1-a)*(1-a)*g3 );
	return y;
}


TString ToString(int num){
	ostringstream start;
	start<<num;
	TString start1=start.str();
	return start1;	
}


int Advanced_cal(){

	double trash;
	double temp1;
	double temp2;
	double tof;
	double dist;
	
	unsigned i=0, N = 1e7;
	//Counter vari
	unsigned Missed =0;
	

	//Variabili dell'evento
	Double_t Xu=0, Yu=0, Xd=0, Yd=0, C_Theta=0, Phi=0;
	//Variabile di un evento andato a segno
	Double_t E=0, beta_mu=0, Ts=0, Td=0, Vs=0, Vd=0;
	
	//Variabile di ricostruzione
	Double_t Tsr=0, Tdr =0, Xur=0, beta_mur=0, Vd2=0, beta_mur2=0, Td2r=0, Xur2=0, Td2=0;
	
	//TString filename = "data/distribution/1V6_barcal3.dat";
	//TString filename = "../data/2G7_beta1_mis1.dat";
	//TString filename = "../data/3L1_beta1_mis2.dat";
	//TString filename = "../data/2V2_beta0_mis2.dat";
	//TString filename = "../data/resolutions/2M9_fine_tcal23_50_100.dat";
	//TString filename = "../data/resolutions/2G2_fine_tcal23_200_100.dat";
	//TString filename = "../data/resolutions/2G1_fine_tcal23_150_100.dat";
	//TString filename = "../data/resolutions/2L4_fine_tcal12_200_100.dat";
	

	TH1D * hist_tof = new TH1D("hist_tof", "hist_tof", (TFIN - TIN)/BIN_SIZE, TIN, TFIN);
	TH1D * hist_tof_cen = new TH1D("hist_tof_cen", "hist_tof_cen", 1400, -20, 50);
	TH1D * hist_dist = new TH1D("hist_dist", "hist_dis", 100, -100, 400);
	TH1D * hist_beta = new TH1D("hist_beta", "hist_beta", 300, 0, 10);
	TH2D * hist_dist_tof = new TH2D("hist_dist_tof", "hist_dist_tof", 120, -30, 330, 50, -5, 20);	
	
	std::ifstream myfile;
	myfile.open(filename);

   while ( true ) {
     myfile >> trash;
     myfile >> temp1;
     myfile >> temp2;
     if(temp2 != 0){
		tof = (temp1/a1 - temp2/a2)/2 ;
		dist = -(temp2/a2 + temp1/a1)/a3;
		hist_dist -> Fill(dist);
		hist_dist_tof -> Fill(dist, tof);
		if ( dist > 250 && dist < 500) hist_tof_cen -> Fill(tof);
		dist = TMath::Sqrt((dist-140)*(dist-140)+Zu*Zu);
		hist_tof -> Fill(tof); 
		hist_beta -> Fill(dist/tof/30);
	 }
     if( myfile.eof() ) break;
   };
   	fitl = new TF1("fitl", pdf,TIN, TFIN,5); 
	fitl -> SetParameter(0,0.03);
	fitl -> SetParameter(1,hist_tof->GetStdDev()*8);
	fitl -> SetParameter(2,hist_tof->GetStdDev()/10);
	fitl -> SetParameter(3,hist_tof->GetMean()*1.01);
	fitl -> SetParameter(4,hist_tof->GetEntries()*1.1);

   
   TCanvas * c1 = new TCanvas("c1", "c1", 1);
   c1 -> cd();
   hist_tof -> GetXaxis() -> SetTitle("Time of flight [ns]");
   hist_tof -> GetYaxis() -> SetTitle("Number of events");
   hist_tof -> Draw();
   	hist_tof -> Fit(fitl, "L", "R", TIN, TFIN);
	fitl -> Draw("same");
   cout << hist_tof -> GetMean() << endl;
   



return 0;
}