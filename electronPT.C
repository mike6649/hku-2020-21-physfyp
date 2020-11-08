#include "TChain.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include <iostream>

using namespace std;

void electronPT (char* fileName) {
  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(fileName);
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");

  // Book histograms
  TH1 *histElectronPT = new TH1F("P_{T}", "P_{T}", 50, 0.0, 100.0);
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    if (branchElectron->GetEntries() > 0)
    {

      Electron *e1 = (Electron*) branchElectron->At(0);
      

      // Plot electron PT
      histElectronPT->Fill(e1->PT);

      // Print electron PT
      cout << e1->PT << endl;
    }

  }

  // Show resulting histograms
  histElectronPT->Draw();
}