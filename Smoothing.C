#include <iostream>

void Smoothing() {
	// initialise variables and open files
	TFile f("ele_scale_factors_v3.root");
	TFile h("Smoothing.root", "RECREATE");
	TH2F *ScaleFac = (TH2F *)f.Get("ele_scale_factors");
	TH2F *ScaleFacUncert = (TH2F *)f.Get("ele_scale_factors_uncertainties");
	//project Histograms
	TH1D *h_1D_EtaBin1 = ScaleFac->ProjectionY("SF for 0<Eta<0.8", 1, 1);
	TH1D *h_1D_EtaBin2 = ScaleFac->ProjectionY("SF for 0.8<Eta<1.479", 2, 2);
	TH1D *h_1D_EtaBin3 = ScaleFac->ProjectionY("SF for 1.479<Eta<2", 3, 3);
	TH1D *h_1D_EtaBin4 = ScaleFac->ProjectionY("SF for 2<Eta<2.5", 4, 4);
	//set errors
	int Nmax = h_1D_EtaBin1->GetNbinsX();
	for (int i=0; i<Nmax; ++i){	
		TH1D *h_1D_EtaBin1_Uncert = ScaleFacUncert->ProjectionY("SFUncert1", 1, 1);
		TH1D *h_1D_EtaBin2_Uncert = ScaleFacUncert->ProjectionY("SFUncert2", 2, 2);
		TH1D *h_1D_EtaBin3_Uncert = ScaleFacUncert->ProjectionY("SFUncert3", 3, 3);
		TH1D *h_1D_EtaBin4_Uncert = ScaleFacUncert->ProjectionY("SFUncert4", 4, 4);
		double Error1 = h_1D_EtaBin1_Uncert->GetBinContent(i);
		h_1D_EtaBin1->SetBinError(i, Error1);
		double Error2 = h_1D_EtaBin2_Uncert->GetBinContent(i);
		h_1D_EtaBin2->SetBinError(i, Error2);
		double Error3 = h_1D_EtaBin3_Uncert->GetBinContent(i);
		h_1D_EtaBin3->SetBinError(i, Error3);
		double Error4 = h_1D_EtaBin4_Uncert->GetBinContent(i);
		h_1D_EtaBin4->SetBinError(i, Error4);
	}
	h_1D_EtaBin1->SetTitle("SF for 0<Eta<0.8");
	h_1D_Eta
	h_1D_EtaBin1->Write();
	h_1D_EtaBin2->SetTitle("SF for 0.8<Eta<1.479");
	h_1D_EtaBin2->Write();
	h_1D_EtaBin3->SetTitle("SF for 1.479<Eta<2");
	h_1D_EtaBin3->Write();
	h_1D_EtaBin4->SetTitle("SF for 2<Eta<2.5");
	h_1D_EtaBin4->Write();
}

