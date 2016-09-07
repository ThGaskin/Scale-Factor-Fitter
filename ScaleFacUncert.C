#include <iostream>
#include "Analyse.h"

void ScaleFacUncert_v2() {
	// initialise variables and open files
	TFile f("ele_scale_factors_v3.root");
	TFile *g = TFile::Open("qqZZ.root");
	TFile h("SF_qqZZ.root", "RECREATE");
	TH2F *ScaleFac = (TH2F *)f.Get("ele_scale_factors");
	TH2F *ScaleUncert = (TH2F *)f.Get("ele_scale_factors_uncertainties");
	//create Histograms
	p_Scalefactor = new TProfile("p_Scalefactor", ";ZZ mass; scale facor", 1000,0, 1000, 0.8, 1.2);
	h_2D_Scalefactor = new TH2F("h_2D_Scalefactor", "", 1000,0, 800, 100, 0.8, 1.2);
	p_Uncertainty_Corr = new TProfile ("p_Uncertainty_Corr", "Correlated SF Uncertainties", 1000,0,1000,0,1);
	p_Uncertainty_Uncorr = new TProfile ("p_Uncertainty_Uncorr", "Uncorrelated SF Uncertainties", 1000,0,1000,0,1);
	p_Uncertainty_Anticorr = new TProfile ("p_Uncertainty_Anticorr", "Anticorrelated SF Uncertainties", 1000,0,1000,0,1);
	p_Uncertainty_Semicorr = new TProfile ("p_Uncertainty_Semicorr", "Half-correlated SF Uncertainties (Correlation: 0.5)", 1000,0,1000,0,1);
	//create Tree Reader and its data containers
	TTreeReader myReader("ZZTree/candTree", g);
	TTreeReaderValue<short> Z1Flav(myReader, "Z1Flav");
	TTreeReaderValue<short> Z2Flav(myReader, "Z2Flav");
	TTreeReaderValue<float> ZZMass(myReader, "ZZMass");
	TTreeReaderArray<float> LepPt(myReader, "LepPt");
	TTreeReaderArray<float> LepEta(myReader, "LepEta");
	//loop over all events, selecting only electrons
	while (myReader.Next()) {
		if (*Z1Flav * *Z2Flav ==121*121) {
			std::vector<float> ScaleFactor;
			std::vector<float> Uncertainty;
			//initialise correlation matrices
			float CorrMatrix [4][4];
			float AnticorrMatrix [4][4];
			float SemicorrMatrix [4][4];
			//Fill correlation matrices
			for (int i=0; i<4; ++i){
				for (int j=0; j<4; ++j){
					if (i==j) {
						CorrMatrix[i][j]=1;
						AnticorrMatrix[i][j]=1;
						SemicorrMatrix[i][j]=1;
					}
					else {
						CorrMatrix[i][j]=1;
						AnticorrMatrix[i][j]=(-1);
						SemicorrMatrix[i][j]=0.5;
					}					
				}
		
			}	
			float s=1;		
			float delta_s_Corr;
			float delta_s_Uncorr;
			float delta_s_Anticorr;
			float delta_s_Semicorr;
			//get scale factors and uncerts for each Lepton, and add to vector
			for (int i=0; i<4; ++i){
				float lookupPt=LepPt[i];
				if (lookupPt>200) lookupPt = 200;
				ScaleFactor.push_back(ScaleFac->GetBinContent(ScaleFac->GetXaxis()->FindBin(abs(LepEta[i])), ScaleFac->GetYaxis()->FindBin(lookupPt)));	
				Uncertainty.push_back(ScaleUncert->GetBinContent(ScaleUncert->GetXaxis()->FindBin(abs(LepEta[i])), ScaleUncert->GetYaxis()->FindBin(lookupPt)));
		
			}
			float delta_s_Uncorr_sq=0;
			float delta_s_Corr_sq=0;
			float delta_s_Anticorr_sq=0;
			float delta_s_Semicorr_sq=0;
			//calculate s and delta_s for ensemble:Gaussian
			for (int i=0; i<4; ++i){
				s=s*ScaleFactor.at(i);
				delta_s_Uncorr_sq=(Uncertainty.at(i)*Uncertainty.at(i))/(ScaleFactor.at(i)*ScaleFactor.at(i));
			}
			delta_s_Uncorr=sqrt(delta_s_Uncorr_sq);
			//calculate s and delta_s for ensemble:non-Gaussian
			//First: Matrixmultiplication
			std::vector<float> Multipliedvector_Corr={0,0,0,0};
			std::vector<float> Multipliedvector_Anticorr={0,0,0,0};
			std::vector<float> Multipliedvector_Semicorr={0,0,0,0};
			for (int j=0; j<4; ++j){
				for (int i=0; i<4; ++i){
					Multipliedvector_Corr.at(j)=(Multipliedvector_Corr.at(j)+(Uncertainty.at(i)/(ScaleFactor.at(i))*CorrMatrix[i][j]));
					Multipliedvector_Anticorr.at(j)=(Multipliedvector_Anticorr.at(j)+(Uncertainty.at(i)/(ScaleFactor.at(i))*AnticorrMatrix[i][j]));
					Multipliedvector_Semicorr.at(j)=(Multipliedvector_Semicorr.at(j)+(Uncertainty.at(i)/(ScaleFactor.at(i))*SemicorrMatrix[i][j]));
					}
			}
			//Second: Scalar product
			for (int k=0; k<4; ++k){
				delta_s_Corr_sq=(delta_s_Corr_sq+(Multipliedvector_Corr.at(k)*Uncertainty.at(k)/(ScaleFactor.at(k))));
				delta_s_Anticorr_sq=(delta_s_Anticorr_sq+(Multipliedvector_Anticorr.at(k)*Uncertainty.at(k)/(ScaleFactor.at(k))));
				delta_s_Semicorr_sq=(delta_s_Semicorr_sq+(Multipliedvector_Semicorr.at(k)*Uncertainty.at(k)/(ScaleFactor.at(k))));
			}
			delta_s_Corr=sqrt(delta_s_Corr_sq);
			delta_s_Anticorr=sqrt(delta_s_Anticorr_sq);
			delta_s_Semicorr=sqrt(delta_s_Semicorr_sq);
			/*//check Lepton Eta
			int EtaCounter=0;
			for (int i=0; i<4; ++i){
				if (abs(LepEta[i])>1.1479) ++EtaCounter;
			}
			if (EtaCounter==0) h_2D_4Barrel->Fill(s);
			if (EtaCounter==1) h_2D_3Barrel->Fill(s);
			if (EtaCounter==2) h_2D_2Barrel->Fill(s);	
			if (EtaCounter==3) h_2D_1Barrel->Fill(s);*/
			//Fill Histograms
			p_Scalefactor->Fill(*ZZMass, s);
			h_2D_Scalefactor->Fill(*ZZMass, s);
			p_Uncertainty_Corr->Fill(*ZZMass, delta_s_Corr);
			p_Uncertainty_Uncorr->Fill(*ZZMass, delta_s_Uncorr);
			p_Uncertainty_Anticorr->Fill(*ZZMass, delta_s_Anticorr);
			p_Uncertainty_Semicorr->Fill(*ZZMass, delta_s_Semicorr);
		}
	
	}
	h_2D_Scalefactor->Write();
	p_Scalefactor->Write();
	p_Uncertainty_Corr->Write();
	p_Uncertainty_Uncorr->Write();
	p_Uncertainty_Anticorr->Write();
	p_Uncertainty_Semicorr->Write();
}



