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
//Velocit� nella luce nel vuoto
Double_t C = 299792458.0 ;
//Dimensioni fisichedel sistema
//Velocit� di propagazione della luce nella sbarra
Double_t beta_s = 0.4697 * C;

//Lunghezze in metri
//Lunghezze lastra sopra
Double_t Ux = 2.79;
Double_t Uy = 0.04;
//Altezza della lastra alta rispetto alla lastra bassa
Double_t Zu = 1.75;

//Lunghezze lastra sotto
Double_t Dx = 0.14;
Double_t Dy = 0.23;


//Posizione del centro della lastra sotto rispetto a quella sopra
Double_t Dcx = 1.40;
Double_t Dcy = 0.02;

//Qulache variabile per la simulazione
//Massa del muone in GeV
Double_t M_mu = 0.105;

//Intervallo di energia della simulazione in GeV
Double_t Emin = M_mu;
Double_t Emax = M_mu + 0.05 ; 


//Variabile per simulare la lettura del primo TAC
Double_t a2=0.0228*1e9, c2=-0.35, s2=0.05, delay2= 30.5 * 1e-9;

//Variabile per simulare la lettura del secondo TAC
Double_t a1=0.0200*1e9, c1=-0.322, s1=0.05, delay1= 30.5 * 1e-9;


void Montecarlo_confronto(){
	//Counter e numero eventi da generare
	unsigned i=0, N = 1e8;
	//Counter vari
	unsigned Missed =0;
	
	
	TH1F* histo_beta = new TH1F("histo_beta","histo_beta", 1000,0,8*1e8);
	TH1F* histo_TOF = new TH1F("histoTOF","histoTOF", 1000,0,2*1e-8);
	TH1F* histo_x = new TH1F("histox","histox", 280,-0.005,2.795);
	
	TH1F* histo_beta2 = new TH1F("histo_beta2","histo_beta2", 1000,0,8*1e8);
	TH1F* histo_TOF2 = new TH1F("histoTOF2","histoTOF2", 1000,0,2*1e-8);
	TH1F* histo_x2 = new TH1F("histox2","histox2", 280,-0.005,2.795);
	//Variabili dell'evento
	Double_t Xu=0, Yu=0, Xd=0, Yd=0, C_Theta=0, Phi=0;
	//Variabile di un evento andato a segno
	Double_t E=0, beta_mu=0, Ts=0, Td=0, Vs=0, Vd=0, Td2=0, Vd2=0;
	
	//Variabile di ricostruzione
	Double_t Tsr=0, Tdr =0,Td2r =0, Xur=0, Xur2 =0, beta_mur=0, beta_mur2=0;
	
	
	//Apro il file su cui salvare
	/**char* filename= "Montecarlo.txt";
	ofstream fout(filename);
	if ( ! fout ) { cerr << " can't open input - " << filename <<endl; return 1; }
	fout.precision(3);**/
	
	//Ciclo di generazione
	for (i=0; i<N; i++){
		
		//Estraggo posizione e direzione iniziale
		Xu = (Rx -> Rndm()) * Ux;
		Yu = (Rx -> Rndm()) * Uy;
		C_Theta = TMath::Power(Rx -> Rndm(),1.0/3.0);
		Phi = (Rx -> Rndm()) * 2 * (TMath::Pi());
		
		//Calcolo posizione finale
		Xd = Xu - Zu * TMath::Cos(Phi) * sqrt(1 - pow(C_Theta, 2.0))/C_Theta;
		Yd = Yu - Zu * TMath::Sin(Phi) * sqrt(1 - pow(C_Theta, 2.0))/C_Theta;
		if(Dcx - Dx/2 <= Xd && Xd<= Dcx + Dx/2 && Dcy - Dy/2 <= Yd && Yd <= Dcy + Dy/2){
			//Genero l'energia e la velocit�
			E = 111;//(Rx -> Rndm()) *( Emax - Emin) + Emin;
			beta_mu = sqrt (1 - pow(M_mu/E, 2.0)) * C;
			
			//Genero i tempi in lettura nella barra
			Ts = (Ux - 2 * Xu)/beta_s + delay1;
			Vs = Rx -> Gaus(a1 * Ts + c1,s1);
			
			//Genero i TOF
			Td = sqrt(pow(Zu,2.0) + pow( Xu-Xd, 2.0) + pow( Yu-Yd, 2.0))/beta_mu + delay2 - Xu/beta_s;
			Vd = Rx -> Gaus(a2 * Td + c2,s2);
			
			Td2 = sqrt(pow(Zu,2.0) + pow( Xu-Xd, 2.0) + pow( Yu-Yd, 2.0))/beta_mu + delay1 + Xu/beta_s;
			Vd2 = Rx -> Gaus(a1 * Td2 + c1,s1);
			//Ricostruisco
			Tdr = (Vd-c2)/a2;
			Tsr = (Vs-c1)/a1;
			Td2r = (Vd2 - c1)/a1;

			Xur = (Ux-(Tsr - delay1)*beta_s)/2;
			Xur2 = (Td2r - Tdr - delay1 +delay2)*beta_s/2;
			beta_mur = sqrt(pow(Zu,2.0) + pow( Xur-Dcx, 2.0))/(Tdr -delay2 + Xur/beta_s);
			beta_mur2 = sqrt(pow(Zu,2.0) + pow( Xur2-Dcx, 2.0))/(Tdr -delay2 + Xur2/beta_s);
			
			histo_beta -> Fill(beta_mur);
			histo_TOF -> Fill(Tdr -delay2 + Xur/beta_s);
			histo_x -> Fill(Xur);
			histo_beta2 -> Fill(beta_mur2);
			histo_TOF2 -> Fill(Tdr -delay2 + Xur2/beta_s);
			histo_x2 -> Fill(Xur2);
			// Stampo le variabili generate e la velcit� reale
			//fout  << Vs << '\t'  << Vd << '\t'  << beta_mu << endl;
		}
		else {Missed ++;}
	}
	//fout.close();
	//Disegno il tutto
	TCanvas *c_beta = new TCanvas("c_beta","c_beta",1);  	 
	c_beta -> cd();
	histo_beta -> Draw();
	histo_beta2 -> Draw("same");
	histo_beta2 -> SetLineColor(2);
	
	TCanvas *c_TOF = new TCanvas("c_TOF","c_TOF",2);  	 
	c_TOF -> cd();
	histo_TOF -> Draw();
	histo_TOF2 -> Draw("same");
	histo_TOF2 -> SetLineColor(2);
	
	TCanvas *c_x = new TCanvas("c_x","c_x",3);  	 
	c_x -> cd();
	histo_x -> Draw();
	histo_x2 -> Draw("same");
	histo_x2 -> SetLineColor(2);
	
	cout << histo_beta2 -> GetMean() << '\t' << histo_beta2 -> GetStdDev()<< endl;
}
