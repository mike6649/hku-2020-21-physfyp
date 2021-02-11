// This code is based on the paper:
// Search for chargino and neutralino production in final states with two same-sign leptons,
// jets and missing transverse momentum at sqrt(s) = 13 TeV with the ATLAS detector
// Cheuk Yee LO, 2019

// Usage:
// g++ -o likelihood.exe likelihood.C `root-config --cflags --libs`
// ./likelihood.exe
// https://root.cern.ch/root/html534/tutorials/fit/NumericalMinimization.C.html
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

// represent one bin representing the pT and eta
class Bin {

public:
  // Bin(i,j) is the bin with
  // pT_bins[i-1] <= pT < pT_bins[i-1] , (in MeV)
  // eta_bins[j-1] <= eta < eta_bins[j-1]
  static constexpr Float_t pT_bins[] = {25000,  60000,  90000,
                                        130000, 150000, 1000000};
  static constexpr Float_t eta_bins[] = {0.50, 1, 1.37, 1.52, 1.8, 2, 2.47};
  int pT_bin_index;
  int eta_bin_index;
  Bin(int pT_bin_index, int eta_bin_index) {
    this->pT_bin_index = pT_bin_index;
    this->eta_bin_index = eta_bin_index;
  }
  Bin(Float_t pT, Float_t eta) {
    bool done = 0;
    for (int i = 0; i < size(this->pT_bins); i++) {
      if (pT < this->pT_bins[i]) {
        this->pT_bin_index = i;
        done = true;
        break;
      }
    }
    if (!done)
      this->pT_bin_index = size(this->pT_bins);

    done = 0;
    for (int i = 0; i < size(this->eta_bins); i++) {
      if (eta < this->eta_bins[i]) {
        this->eta_bin_index = i;
        done = true;
        break;
      }
    }
    if (!done)
      this->eta_bin_index = size(this->eta_bins);
  }

  friend ostream &operator<<(ostream &os, const Bin &b) {
    os << "(" << b.pT_bin_index << "," << b.eta_bin_index << ")";
    return os;
  }

  // compare pT before eta
  // no physical signifance, for storage in Map only
  friend bool operator<(const Bin &a, const Bin &b) {
    if (a.pT_bin_index < b.pT_bin_index) {
      return true;
    } else if (a.pT_bin_index > b.pT_bin_index) {
      return false;
    } else {
      if (a.eta_bin_index < b.eta_bin_index) {
        return true;
      }
    }
    return false;
  }
};

// Store number of events binned by 2 pairs of pt,eta bins
// Represents N_{ij}
class Events {
public:
  Events() {}

  std::map<Bin, std::map<Bin, Long64_t>> myMap;   // first key is the leading lepton bin, second is the subleading lepton bin

  Long64_t get(const Bin &leadingLepton, const Bin &subLeadingLepton) {
    if (!this->myMap.count(leadingLepton) ||
        !this->myMap[leadingLepton].count(subLeadingLepton)) {
      // cout << "EVENTS no mapping at " << leadingLepton << " " <<
      // subLeadingLepton;
      return 0;
    } else {
      return this->myMap[leadingLepton][subLeadingLepton];
    }
  }

  void put(const Bin& leadingLepton, const Bin& subLeadingLepton) {
    if (!this->myMap.count(leadingLepton)) {
      std::map<Bin, Long64_t> newMap;
      newMap[subLeadingLepton] = 0;
      this->myMap[leadingLepton] = newMap;
    }
    if (!this->myMap[leadingLepton].count(subLeadingLepton)) {
      this->myMap[leadingLepton][subLeadingLepton] = 0;
    }
    this->myMap[leadingLepton][subLeadingLepton] += 1;
    // cout << leadingLepton << " " << subLeadingLepton <<
    // myMap[leadingLepton][subLeadingLepton] << endl;
  }

  friend ostream &operator<<(ostream &os, Events &events) {
    for (int i = 0; i < size(Bin::pT_bins); i++) {
      for (int j = 0; j < size(Bin::eta_bins); j++) {
        Bin a(i, j);
        os << "=======\n Leading Electron in " << a << endl;
        for (int k = 0; k < size(Bin::pT_bins); k++) {
          for (int l = 0; l < size(Bin::eta_bins); l++) {
            os << events.get(a, Bin(k, l)) << "\t";
          }
          os << "\n";
        }
      }
      os << endl;
    }
    return os;
  }
};

Float_t likelihood_function(const Double_t &totalEvents, const Double_t &ssEvents,
                       Float_t e_i, Float_t e_j) {
  return - (ssEvents * TMath::Log(totalEvents*(e_i*(1-e_j)+(1-e_i)*e_j)) - totalEvents*(e_i*(1-e_j)+(1-e_i)*e_j));
}

void likelihood(char *fileName) {

  Events allEvents;

  TFile *myFile = TFile::Open(fileName);
  // Create chain of root trees
  // TTreeReader myReader("nominal_Loose;724", myFile);
  TTree *t1 = (TTree *)myFile->Get("nominal_Loose;724");

  vector<int> *el_charge = 0;
  TBranch *b_el_charge = 0;

  vector<Float_t> *el_pt = 0;
  TBranch *b_el_pt = 0;

  vector<Float_t> *el_eta = 0;
  TBranch *b_el_eta = 0;
  // t1->Print();
  t1->ResetBranchAddresses();
  t1->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
  t1->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
  t1->SetBranchAddress("el_eta", &el_eta, &b_el_eta);

  Long64_t nentries = t1->GetEntries();

  for (Long64_t i = 0; i < nentries; i++) {
    Long64_t tentry = t1->LoadTree(i);
    b_el_charge->GetEntry(tentry);
    b_el_pt->GetEntry(tentry);
    b_el_eta->GetEntry(tentry);

    if ((*el_pt)[0] > (*el_pt)[1]) { // first element is leading lepton
      allEvents.put(Bin((*el_pt)[0], (*el_eta)[0]),
                    Bin((*el_pt)[1], (*el_eta)[1]));
    } else {
      allEvents.put(Bin((*el_pt)[1], (*el_eta)[1]),
                    Bin((*el_pt)[0], (*el_eta)[0]));
    }
  }

  cout << allEvents;

  t1->ResetBranchAddresses();

  myFile->Close();
}

int main(int argc, char **argv) {
  likelihood(argv[1]);
  return 0;
}