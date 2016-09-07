#include <iostream>
#include "Analyse.h"

void relerror() {
	// initialise variables and open files
	TFile f("ele_scale_factors_v3.root");
	TFile *g = TFile::Open("qqZZ.root");
	TFile h("ScaleFactors2.root", "RECREATE");
	TH1F *h_ZZMass;
	TH1F *h_ZZMass_weights;
	TH2F *ScaleFac = (TH2F *)f.Get("ele_scale_factors");
	TH2F *ScaleUncert = (TH2F *)f.Get("ele_scale_factors_uncertainties");
	
	//create Histograms
	h_ZZMass = new TH1F("ZZMass", "Unweighted ZZMass", 1000, 0, 800);
	h_ZZMass_weights = new TH1F("ZZMass_weights", "Weighted ZZ Mass", 1000, 0, 800);
	p_Scalefactor = new TProfile("p_Scalefactor", ";ZZ mass; scale facor", 100,0, 1000, 0.8, 1.2);
	p_Uncertainty = new TProfile("p_Uncertainty","ZZ mass scale factor uncertainty",100,0,1000,0,1);
	p_Uncert_Corr = new TProfile("p_Uncert_Corr", "Correlated uncertainties", 100,0,1000,0,1);
	h_2D_Scalefactor = new TH2F("h_2D_Scalefactor", "", 1000,0, 800, 100, 0.8, 1.2);
	h_2D_1Barrel = new TH1F("h_2D_1Barrel", "Scalefactor for one electron in barrel", 50,0.8, 1.2);
	h_2D_2Barrel = new TH1F("h_2D_2Barrel", "Scalefactor for two electrons in barrel", 50,0.8, 1.2);
	h_2D_3Barrel = new TH1F("h_2D_3Barrel", "Scalefactor for three electrons in barrel", 50,0.8, 1.2);
	h_2D_4Barrel = new TH1F("h_2D_4Barrel", "Scalefactor for four electrons in barrel", 50,0.8, 1.2);

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
			float s=1;
			float delta_s;

			//get scale factors and uncerts for each Lepton, and add to vector
			for (int i=0; i<4; ++i){
				float lookupPt=LepPt[i];
				if (lookupPt>200) lookupPt = 200;
				ScaleFactor.push_back(ScaleFac->GetBinContent(ScaleFac->GetXaxis()->FindBin(abs(LepEta[i])), ScaleFac->GetYaxis()->FindBin(lookupPt)));	
				Uncertainty.push_back(ScaleUncert->GetBinContent(ScaleUncert->GetXaxis()->FindBin(abs(LepEta[i])), ScaleUncert->GetYaxis()->FindBin(lookupPt)));
		
			}
			float delta_s_sq=0;
			//calculate s and delta_s for ensemble:Gaussian
			for (int i=0; i<4; ++i){
				s=s*ScaleFactor.at(i);
				delta_s_sq=(Uncertainty.at(i)*Uncertainty.at(i))/(ScaleFactor.at(i)*ScaleFactor.at(i));
			}
			delta_s=sqrt(delta_s_sq);
			//check Lepton Eta
			int EtaCounter=0;
			for (int i=0; i<4; ++i){
				if (abs(LepEta[i])>1.1479) ++EtaCounter;
			}
			if (EtaCounter==0) h_2D_4Barrel->Fill(s);
			if (EtaCounter==1) h_2D_3Barrel->Fill(s);
			if (EtaCounter==2) h_2D_2Barrel->Fill(s);	
			if (EtaCounter==3) h_2D_1Barrel->Fill(s);
			//Fill Histograms
			p_Scalefactor->Fill(*ZZMass, s);
			h_2D_Scalefactor->Fill(*ZZMass, s);
			p_Uncertainty->Fill(*ZZMass, delta_s);
			h_ZZMass_weights->Fill(*ZZMass, s);
			h_ZZMass->Fill(*ZZMass);	
		}
	}
	h_ZZMass->Write();
	h_ZZMass_weights->Write();
	p_Scalefactor->Write();
	p_Uncertainty->Write();
	h_2D_Scalefactor->Write();
	h_2D_1Barrel->Write();
	h_2D_2Barrel->Write();
	h_2D_3Barrel->Write();
	h_2D_4Barrel->Write();
}






/*
Float_t relerror (Float_t lepPt, Float_t Eta){ 
	//cout<<"Enter electron momentum"<<endl;
	//cin >> lepPt;
	//cout<<"Enter electron eta"<<endl;
	//cin >> Eta;
	if (lepPt>199) {
		cout<<"Lepton momentum too high"<<endl;
		return 0;
	}
	Float_t s=1;
	Float_t deltas=2;
	TFile f("ele_scale_factors_v3.root");
	TH2F *ScaleFac = (TH2F *)f.Get("ele_scale_factors");
	TH2F *ScaleUncert = (TH2F *)f.Get("ele_scale_factors_uncertainties");
	s =ScaleFac->GetBinContent(ScaleFac->GetXaxis()->FindBin(Eta), ScaleFac->GetYaxis()->FindBin(lepPt));
	deltas =ScaleUncert->GetBinContent(ScaleUncert->GetXaxis()->FindBin(Eta), ScaleFac->GetYaxis()->FindBin(lepPt));
	cout<<"deltas:"<<deltas<<endl;

	cout<<"s:"<<s<<endl;
	return (deltas*deltas)/(s*s);
	
}
int main(){
	cout<<(relerror(157, 1.5))<<endl;
}
*/

