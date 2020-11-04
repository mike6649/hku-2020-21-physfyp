#include "TChain.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include <iostream>
using namespace std;

void test (char* fileName) {
  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(fileName);
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");

  // Book histograms
  TH1 *histElectronPT = new TH1F("electron pt", "electron P_{T}", 50, 0.0, 100.0);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    // If event contains at least 1 electron
    if (branchElectron->GetEntries() > 0)
    {
      // Take first electron
      Electron *electron = (Electron*) branchElectron->At(0);
      
      // Plot electron transverse momentum
      histElectronPT->Fill(electron->PT);
      
      // Print electron transverse momentum
      cout << electron->PT << endl;
    }

  }

  // Show resulting histograms
  histElectronPT->Draw();
}