#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>  
#include <TMinuit.h>
#include <TFile.h>
#include "TFile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TString.h"
#include <fstream>
#include <string>
#include <cmath>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TRandom.h>
#include <TFitResultPtr.h>
#include <RConfig.h>
#include <TStopwatch.h>

TRandom3* Rx = new TRandom3(23);
//Velocità nella luce nel vuoto
Double_t C = 299792458.0 ;
//Dimensioni fisichedel sistema
//Velocità di propagazione della luce nella sbarra
Double_t beta_s = 0.4697 * C;

//Lunghezze in metri
//Lunghezze lastra sopra
Double_t Ux = 2.79;
Double_t Uy = 0.04;
//Altezza della lastra alta rispetto alla lastra bassa
Double_t Zu = 0.05;

//Lunghezze lastra sotto
Double_t Dx = 0.14;
Double_t Dy = 0.23;


//Posizione del centro della lastra sotto rispetto a quella sopra
Double_t Dcx = 1.40;
Double_t Dcy = 0.02;

//Prametri per il fit
Double_t BIN_SIZE = 1e-10;
Double_t TIN = -100e-9;
Double_t TFIN = 100e-9;
Double_t S_small =2e-9;
Double_t S_big = 10e-9;

Double_t pdf(Double_t* x, Double_t *par)
{
Double_t y=0;
Double_t  a=par[0];
Double_t  s1=par[1];
Double_t  s2=par[2];
Double_t  offset=par[3];
Double_t  Norm=par[4];


Double_t s1_eff=s1*sqrt(2.0);
Double_t s3_eff=s2*sqrt(2.0);
Double_t s2_eff=sqrt(s1*s1+s2*s2);

Double_t g1=1/(sqrt(2*TMath::Pi())*s1_eff)*exp(-(pow((x[0]-offset)/(2*s1_eff),2.0)));
Double_t g2=1/(sqrt(2*TMath::Pi())*s2_eff)*exp(-(pow((x[0]-offset)/(2*s2_eff),2.0)));
Double_t g3=1/(sqrt(2*TMath::Pi())*s3_eff)*exp(-(pow((x[0]-offset)/(2*s3_eff),2.0)));
y= Norm*BIN_SIZE*(a*a*g1+ 2*a*(1-a)*g2 +(1-a)*(1-a)*g3 );
	return y;
}
	

void Tiny_x_MC(){
	//Counter e numero eventi da generare
	unsigned i=0, N = 1e6;
	//Counter vari
	unsigned Missed =0;
	
	
	
	
	//Variabili dell'evento
	Double_t Xu=0, Yu=0, Xd=0, Yd=0, C_Theta=0, Phi=0;
	//Variabile di un evento andato a segno
	Double_t E=0, beta_mu=0, Ts=0, Td=0, Vs=0, Vd=0, temp1=0,temp2=0, ratio=0.1,Tcen=0;
	
	TH1F* histo = new TH1F ("histo", "histo", (TFIN - TIN)/BIN_SIZE, TIN, TFIN);
	TH1F* histo2 = new TH1F ("histo2", "histo2", (TFIN - TIN)/BIN_SIZE, TIN, TFIN);
	TH1F* histo3 = new TH1F ("histo2", "histo2", (TFIN - TIN)/BIN_SIZE, TIN, TFIN);
	Tcen = (Ux - 2 * Dcx)/beta_s ;
	//Ciclo di generazione
	for (i=0; i<N; i++){
		
		//Estraggo posizione e direzione iniziale
		Xu = ((Rx -> Rndm())-0.5) * Dx + Dcx;
		Yu = ((Rx -> Rndm())-0.5) * Dy + Dcy;
		C_Theta =TMath::Power(Rx -> Rndm(),1.0/3.0);
		Phi = (Rx -> Rndm()) * 2 * (TMath::Pi());
		
		//Calcolo posizione finale
		Xd = Xu - Zu * TMath::Cos(Phi) * sqrt(1 - pow(C_Theta, 2.0))/C_Theta;
		Yd = Yu - Zu * TMath::Sin(Phi) * sqrt(1 - pow(C_Theta, 2.0))/C_Theta;
		if(0 <= Xd && Xd<= Ux && 0 <= Yd && Yd <= Uy){
			
			//Genero i tempi in lettura nella barra
			Ts = (Ux - 2 * Xd)/beta_s ;
			if(Rx -> Rndm()>ratio){ temp1=Rx->Gaus(0,S_small);}
			else{temp1=Rx->Gaus(0,S_big);}
			if(Rx -> Rndm()>ratio){ temp2=Rx->Gaus(0,S_small);}
			else{temp2=Rx->Gaus(0,S_big);}
			histo2 -> Fill(Tcen+temp1+temp2);
			histo -> Fill (Ts);
			histo3 -> Fill(temp1+temp2+Ts);
		}
		else {Missed ++;}
	}
	fitl = new TF1("fitl", pdf,TIN, TFIN,5); 
	fitl -> SetParameter(0,ratio);
	fitl -> SetParameter(1,S_big);
	fitl -> SetParameter(2,S_small);
	fitl -> SetParameter(3,Tcen);
	fitl -> SetParameter(4,1);
	histo3 -> Scale(1/histo3->GetEntries());
	histo3 -> Fit(fitl, "L", "R", TIN, TFIN);
	
	
	
	
	//histo -> Draw();
	//histo2 -> Draw();
	//histo2-> SetLineColor(2);
	histo3 -> Draw();
	//histo3-> SetLineColor(3);
	fitl -> Draw("same");
}
