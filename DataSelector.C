#include "Analyse.h"
#include "Analyse.C"
#include "TreeAnalyser.C"
#include "AnalyseTree.C"

class FullEventDataSelector : public TSelector {
public :

   TH1   *ZZMass;  // ZZMass hist

   // Variables used to access and store the data
   TTreeReader fReader;                            // The tree reader
   TTreeReaderArray<Float_t> ZZMass;      // ZZ Mass
   TTreeReaderArray<Float_t> ZZPt;  // ZZ momentum

   FullEventDataSelector(TTree * = 0): ZZMass(0), fMass(fReader, "b_ZZMass"),
                                       fMomentum(fReader, "b_ZZPt") { }
   virtual ~FullEventDataSelector() { }
   virtual void    Init(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual Bool_t  Process(Long64_t entry);
   virtual void    Terminate();
   virtual Int_t   Version() const { return 2; }

   ClassDef(FullEventDataSelector,0);
};

void FullEventDataSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Associate the reader and the tree 
   fReader.SetTree(tree);
}

void FullEventDataSelector::SlaveBegin(TTree *tree)
{
   // SlaveBegin() is a good place to create histograms. 
   // For PROOF, this is called for each worker.
   // The TTree* is there for backward compatibility; e.g. PROOF passes 0.

   hist2 = new TH1F("ZZMass", "ZZ Mass", 10, 0, 100);
   // Add to output list (needed for PROOF)
   GetOutputList()->Add(fPosX);
}

Bool_t FullEventDataSelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree to be 
   // processed. The entry argument specifies which entry in the currently
   // loaded tree is to be processed.
   // It can be passed to either EventSelector::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the TTree.
   //
   // This function should contain the "body" of the analysis: select relevant
   // tree entries, run algorithms on the tree entry and typically fill histograms.

   // *** 1. *** Tell the reader to load the data for this entry:
   fReader.SetEntry(entry);

   // *** 2. *** Do the actual analysis
   for (unsigned int iParticle = 0; iParticle < fMass.GetSize(); ++iParticle) {
      if (fMomentum[iParticle] > 40.0)
         hist2->Fill(fMass[iParticle]);
   }

   return kTRUE;
}

void FullEventDataSelector::Terminate()
{

   hist2->Draw();
}
