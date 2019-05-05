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
Double_t beta_s = 0.5 * C;

//Lunghezze in metri
//Lunghezze lastra sopra
Double_t Ux = 2.79;
Double_t Uy = 0.04;
//Altezza rispetto alla lastra bassa
Double_t Zu = 1.75;

//Lunghezze lastra sotto
Double_t Dx = 0.23;
Double_t Dy = 0.14;
//Altezza rispetto alla lastra bassa
Double_t Zd = 0;

//Posizione del centro della lastra sotto rispetto a quella sopra
Double_t Dcx = 1.40;
Double_t Dcy = 0.02;

//Qulache variabile per la simulazione
//Massa del muone in GeV
Double_t M_mu = 0.105;

//Intervallo di energia della simulazione in GeV
Double_t Emin = M_mu;
Double_t Emax = 1 ; 


//Variabile per simulare la lettura del primo TAC
Double_t a1=0.0228*1e9, c1=-0.35, s1=0.05, delay1= 30.5 * 1e-9;

//Variabile per simulare la lettura del secondo TAC
Double_t a2=0.0228*1e9, c2=-0.35, s2=0.05, delay2= 30.5 * 1e-9;


void Montecarlo_general(){
	//Counter e numero eventi da generare
	unsigned i=0, N = 1e7;
	//Counter vari
	unsigned Missed =0;
	
	
	histo = new TH1F("histo", "histo", 2500, 0, 2.5);
	
	
	//Variabili dell'evento
	Double_t Xu=0, Yu=0, Xd=0, Yd=0, C_Theta=0, Phi=0;
	//Variabile di un evento andato a segno
	Double_t E=0, beta_mu=0, Ts=0, Td=0, Vs=0, Vd=0;
	
	
	//Ciclo di generazione
	for (i=0; i<N; i++){
		
		//Estraggo posizione e direzione iniziale
		Xu = (Rx -> Rndm()) * Ux;
		Yu = (Rx -> Rndm()) * Uy;
		C_Theta = TMath::Power(Rx -> Rndm(),1.0/3.0);
		Phi = (Rx -> Rndm()) * 2 * (TMath::Pi());
		
		//Calcolo posizione finale
		Xd = Xu - Zu * TMath::Cos(Phi) * sqrt(1 - pow(C_Theta, 2.0)/C_Theta);
		Yd = Yu - Zu * TMath::Sin(Phi) * sqrt(1 - pow(C_Theta, 2.0)/C_Theta);
		if(Dcx - Dx/2 <= Xd && Xd<= Dcx + Dx/2 && Dcy - Dy/2 <= Yd && Yd <= Dcy + Dy/2){
			//Genero l'energia e la velocità
			E = (Rx -> Rndm()) *( Emax - Emin) + Emin;
			beta_mu = sqrt (1 - pow(M_mu/E, 2.0)) * C;
			
			//Genero i tempi in lettura nella barra
			Ts = (Ux - 2 * Xu)/beta_s + delay1;
			Vs = Rx -> Gaus(a1 * Ts + c1,s1);
			
			//Genero i TOF
			Td = sqrt(pow(Zu,2.0) + pow( Xu-Xd, 2.0) + pow( Yu-Yd, 2.0))/beta_mu + delay2;
			
			Vd = Rx -> Gaus(a2 * Td + c2,s2);
			
		}
		else {Missed ++;}
	}
}












