//#include "Analyse.h"
//#include "Analyse.C"

void AnalyzeTree()
{
	 // open the file
   	TFile *f = TFile::Open("Higgs.root");
  	 if (f == 0) {
  		 // if cannot open file, print error and return
     		 printf("Error: cannot open file\n");
      		return;
   	}

TFile *MyPlots = new TFile("MyPlots.root", "RECREATE");

// Variables used to store the data
  
TH1F *lep1;  // leading electron
TH1F *lep2;  // next to leading
TH1F *lep3;  // etc.
TH1F *lep4;  
TH1F *ZZMass;
TH1F *ZZMassWeights;
TH1F *leptotal;
TH1F *leptotalWeights;

// create the TH1F histogram
lep1 = new TH1F("lep1", "Leading Electron", 1000, 0, 400);
lep2 = new TH1F("lep2", "2nd Electron", 1000, 0, 400);
lep3 = new TH1F("lep3", "3rd Electron", 1000, 0, 400);
lep4 = new TH1F("lep4", "4th Electron", 1000, 0, 400);
leptotal = new TH1F("leptotal", "Lepton Momentum", 1000, 0, 400);
leptotalWeights = new TH1F("leptotalWeights", "Lepton Momentum with weights", 1000, 0, 400);
ZZMass = new TH1F("ZZMass", "ZZ Mass", 1000, 0, 400);
ZZMassWeights = new TH1F("ZZMassWeights", "ZZ Mass with weights", 1000, 0, 400);

// Create the tree reader and its data containers
TTreeReader myReader("ZZTree/candTree", f);
// First, ascertain lepton flavour
TTreeReaderValue<Short_t> Z1Flav(myReader, "Z1Flav");
TTreeReaderValue<Short_t> Z2Flav(myReader, "Z2Flav");
TTreeReaderValue<Float_t> Mass(myReader, "ZZMass");
TTreeReaderValue<Float_t> Weight(myReader, "dataMCWeight");
// Second, read all LepPt vectors
TTreeReaderArray<Float_t> LepPt(myReader, "LepPt");
   

// loop over all events and select electrons
while (myReader.Next()) {
	if (*Z1Flav * *Z2Flav==169*169) {
   		std::vector<Float_t> newlepPt;
	  	for(int i = 0; i< 4; ++i){	
			newlepPt.push_back(LepPt[i]);
			leptotal->Fill(LepPt[i]);
  			leptotalWeights->Fill(LepPt[i], *Weight);
			cout<<*Weight<<endl;
  		}
		//newlepPt.insert(newlepPt.end(), &LepPt[0], &LepPt[4]);

	   	// sort newlpPt
   		std::sort(newlepPt.begin(), newlepPt.end());
 	        lep1->Fill(newlepPt[0]);
		lep2->Fill(newlepPt[1]);
		lep3->Fill(newlepPt[2]);
		lep4->Fill(newlepPt[3]);
		ZZMass->Fill(*Mass);
   		ZZMassWeights->Fill(*Mass, *Weight);
	}
	
    }

lep1->SetLineColor(kRed);
lep1->Draw();

lep2->SetLineColor(kGreen);
lep2->Draw("SAME");

lep3->SetLineColor(kYellow);
lep3->Draw("SAME");

lep4->SetLineColor(kBlue);
lep4->Draw("SAME");

leptotalWeights->SetLineColor(kRed);

ZZMassWeights->SetLineColor(kRed);

lep1->Write();
lep2->Write();
lep3->Write();
lep4->Write();
leptotal->Write();
leptotalWeights->Write();
ZZMass->Write();
ZZMassWeights->Write();

MyPlots->Close();
f->Close();
}
