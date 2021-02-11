#include "TChain.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include <iostream>

using namespace std;

void misIdCount (char* fileName) {
  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(fileName);
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  cout << "** Chain contains " << numberOfEntries << " events" << endl;
  // Loop over all events

  Long64_t count = 0;
  Long64_t noElectrons = 0;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // cout << "** Branch contains " << branchElectron->GetEntriesFast() << " items" << endl;
    for(Int_t i = 0; i < branchElectron->GetEntriesFast(); ++i)
    {
      noElectrons++;
      
      Electron* electron = (Electron*) branchElectron->At(i);
      GenParticle* particle = (GenParticle*) electron->Particle.GetObject();

      // cout << particle->Charge << " "  << electron->Charge << endl;

      if (particle->Charge  != electron->Charge ) { // charge misid
        count++;
      } 

    }

  }

  // Show resulting histograms
  // histElectronPT->Draw();
  cout << "MisID " << count << " out of " << noElectrons << " electrons" << endl;
}