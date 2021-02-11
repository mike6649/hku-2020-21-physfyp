#include "TChain.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include <iostream>

using namespace std;

void diElectronMass (char* fileName) {
  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(fileName);
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  
  // Book histograms
  TH1 *histDiElectronMass = new TH1F("M", "M_{inv}(e_{1},e_{2})", 100, 40, 140);
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    if (branchElectron->GetEntries() > 1)
    {

      Electron *e1 = (Electron*) branchElectron->At(0);
      Electron *e2 = (Electron*) branchElectron->At(1);
      
      Double_t invMass = (e1->P4()+e2->P4()).M();

      // Plot electron mass
      histDiElectronMass->Fill(invMass);

      // Print electron mass
      cout << invMass << endl;
    }

  }

  // Show resulting histograms
  histDiElectronMass->Draw();
}