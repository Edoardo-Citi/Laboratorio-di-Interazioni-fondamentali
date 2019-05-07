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




void Tiny_x_MC(){
	//Counter e numero eventi da generare
	unsigned i=0, N = 1e8;
	//Counter vari
	unsigned Missed =0;
	
	
	
	
	//Variabili dell'evento
	Double_t Xu=0, Yu=0, Xd=0, Yd=0, C_Theta=0, Phi=0;
	//Variabile di un evento andato a segno
	Double_t E=0, beta_mu=0, Ts=0, Td=0, Vs=0, Vd=0;
	
	TH1F* histo = new TH1F ("histo", "histo", 10000, 25 1e-9, 35* 1e-9);
	
	//Ciclo di generazione
	for (i=0; i<N; i++){
		
		//Estraggo posizione e direzione iniziale
		Xu = (Rx -> Rndm()) * Dx + Dcx;
		Yu = (Rx -> Rndm()) * Dy + Dcy;
		C_Theta =TMath::Power(Rx -> Rndm(),1.0/3.0);
		Phi = (Rx -> Rndm()) * 2 * (TMath::Pi());
		
		//Calcolo posizione finale
		Xd = Xu - Zu * TMath::Cos(Phi) * sqrt(1 - pow(C_Theta, 2.0))/C_Theta;
		Yd = Yu - Zu * TMath::Sin(Phi) * sqrt(1 - pow(C_Theta, 2.0))/C_Theta;
		if(0 <= Xd && Xd<= Ux && 0 <= Yd && Yd <= Uy){
			
			//Genero i tempi in lettura nella barra
			Ts = (Ux - 2 * Xd)/beta_s +30 * 1e-9;
			histo -> Fill (Ts);

		}
		else {Missed ++;}
	}
histo -> Draw();
}
