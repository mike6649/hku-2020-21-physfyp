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
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// #include "Math/KahanSum.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
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
  static constexpr Double_t pT_bins[] = {25000,  60000,  90000,
                                        130000, 150000, 1000000};
  static constexpr Double_t eta_bins[] = {0, 0.50, 1, 1.8, 2, 2.47};

  static constexpr int numPTBins = size(Bin::pT_bins) - 1;
  static constexpr int numEtaBins = size(Bin::eta_bins) - 1;

  int pT_bin_index;
  int eta_bin_index;
  Bin(int pT_bin_index, int eta_bin_index) {
    this->pT_bin_index = pT_bin_index;
    this->eta_bin_index = eta_bin_index;
  }
  Bin(Float_t pT, Float_t eta) {
    this->pT_bin_index = 0;
    this->eta_bin_index = 0;

    for (int i = 0; i < size(this->pT_bins) - 1; i++) {
      if (pT > this->pT_bins[i] && pT < this->pT_bins[i+1]) {
        this->pT_bin_index = i;
        break;
      }
    }

    for (int i = 0; i < size(this->eta_bins) - 1; i++) {
      if (eta > this->eta_bins[i] && eta < this->eta_bins[i+1]) {
        this->eta_bin_index = i;
        break;
      }
    }

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
    for (int i = 0; i < Bin::numPTBins; i++) {
      for (int j = 0; j < Bin::numEtaBins ; j++) {
        Bin a(i, j);
        os << "=======\n Leading Electron in " << a << endl;
        for (int k = 0; k < Bin::numPTBins; k++) {
          for (int l = 0; l < Bin::numEtaBins ; l++) {
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

Events ssEvents;
Events allEvents;

double likelihood_single(const Long64_t &numTotal, const Long64_t &numSS,
                       const double &e_i, const double &e_j) {
  return - (numSS * TMath::Log(numTotal*(e_i*(1-e_j)+(1-e_i)*e_j)) - numTotal*(e_i*(1-e_j)+(1-e_i)*e_j));
}

double likelihood_sum(const double* e) {
  ROOT::Math::KahanSum<double,4> k;
    k.Add(1000000.);
  std::vector<double> numbers;
  // double sum = 0;
    for (int i = 0; i < Bin::numPTBins ; i++) {
      for (int j = 0; j < Bin::numEtaBins ; j++) {
        Bin a(i, j);
        for (int k = 0; k < Bin::numPTBins; k++) {
          for (int l = 0; l < Bin::numEtaBins; l++) {
            Bin b(k,l);
            if (allEvents.get(a,b) == 0 ) {  // avoid divide by zero
              continue;
            }
            numbers.push_back(likelihood_single(allEvents.get(a,b),ssEvents.get(a,b),e[i*Bin::numEtaBins+j],e[k*Bin::numEtaBins+l]));
          }
        }
      }
    }
  k.Add(numbers);
  return k.Sum();
}

inline double getXValueAt(const double* xs, const int& pt_index, const int& eta_index ) {
  return xs[pt_index*Bin::numEtaBins + eta_index];
}

void minimize(){

  const char* minName = "Minuit2";
  const char* algoName = "";
  int randomSeed = -1;
  const int numDimensions = Bin::numPTBins * Bin::numEtaBins;   // 6 pt, 5 eta

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName,algoName);

  // set tolerance, etc...
  min->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2 
  min->SetStrategy(2);
   min->SetMaxIterations(10000000);  // for GSL 
   min->SetTolerance(1);
   min->SetPrintLevel(2);

   // create funciton wrapper for minmizer
   // a IMultiGenFunction type 
   ROOT::Math::Functor f(&likelihood_sum,numDimensions);

  // step
  double step[numDimensions];
  // starting point
  double variable[numDimensions];
   for (int i=0;i<numDimensions;i++){
     step[i] = 0.01;
     variable[i] = 0.002;
   }
 
   min->SetFunction(f);
   // Set the free variables to be minimized!
   for (int i=0;i<numDimensions;i++){
   min->SetVariable(0,"e"+std::to_string(i/Bin::numEtaBins) + std::to_string(i%Bin::numEtaBins),variable[i], step[i]);
   }
  //  for (int i=0;i<5;i++){
  //    min->FixVariable(i);
  //  }
 
   // do the minimization
   min->Minimize();

   const double *xs = min->X();
   for (int i=0;i<Bin::numPTBins;i++){
     for (int j=0;j<Bin::numEtaBins;j++){
      std::cout << xs[i*Bin::numEtaBins + j] << "\t";
     }
     std::cout << std::endl;
   }
   std::cout << min->MinValue() << std::endl;

   auto c = new TCanvas("c","c",600,600);

   double* pts = new double[Bin::numPTBins];

  TMultiGraph *mg = new TMultiGraph();
  TGraph g[Bin::numEtaBins];
   for (int i=0;i<Bin::numEtaBins;i++){
     // construct array of X's for eta_index = i
      for (int j=0;j<Bin::numPTBins;j++){
        pts[j] = getXValueAt(xs,j,i);
      }
    g[i] = TGraph(Bin::numPTBins,Bin::pT_bins , pts);
    mg->Add(g);
   }
   mg->Draw("AL*");
   c->BuildLegend();
}


void likelihood(char *fileName) {

  auto graph = new TH1F("dielectron mass","#",100,0,120000);
  auto c = new TCanvas();

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

  vector<Float_t> *el_e = 0;
  TBranch *b_el_e = 0;

  vector<Float_t> *el_phi = 0;
  TBranch *b_el_phi = 0;

  t1->Print();
  t1->ResetBranchAddresses();
  t1->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
  t1->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
  t1->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
  t1->SetBranchAddress("el_e", &el_e, &b_el_e);
  t1->SetBranchAddress("el_phi", &el_phi, &b_el_phi);

  // MCTruth
  Events mc_total;   // not use first bins
  Events mc_flip;
  vector<int> *el_true_pdg = 0; 
  TBranch *b_el_true_pdg = 0;
  t1->SetBranchAddress("el_true_pdg", &el_true_pdg, &b_el_true_pdg);

  Long64_t nentries = t1->GetEntries();

  Float_t pt_0;
  Float_t pt_1;
  Float_t eta_0;
  Float_t eta_1;
  Float_t phi_0,phi_1;
  Float_t px,py,pz, energy, mass;

  Bin dummy (0,0);

  for (Long64_t i = 0; i < nentries; i++) {
    Long64_t tentry = t1->LoadTree(i);
    b_el_charge->GetEntry(tentry);
    b_el_pt->GetEntry(tentry);
    b_el_eta->GetEntry(tentry);
    b_el_e->GetEntry(tentry);
    b_el_phi->GetEntry(tentry);

    pt_0 = (*el_pt)[0];
    pt_1 = (*el_pt)[1];
    eta_0 = (*el_eta)[0];
    eta_1 = (*el_eta)[1];
    phi_0 = (*el_phi)[0];
    phi_1 = (*el_phi)[1];

    // apply preselection
    px = pt_0 * TMath::Cos(phi_0) + pt_1 * TMath::Cos(phi_1);
    py = pt_0 * TMath::Sin(phi_0) + pt_1 * TMath::Sin(phi_1);
    pz = pt_0 * TMath::SinH(eta_0) + pt_1 * TMath::SinH(eta_1);
    energy = (*el_e)[0] + (*el_e)[1];
    mass = TMath::Power(energy,2) - px*px - py*py - pz*pz;
    mass = TMath::Sqrt(mass);
    // cout << TMath::Sqrt(mass) << endl;
    if (mass < 80000 || mass > 100000) continue;  // mass between 80-100 MeV

    eta_0 = TMath::Abs(eta_0);
    eta_1 = TMath::Abs(eta_1);

    if (pt_0 > pt_1) { // first element is leading lepton
      allEvents.put(Bin(pt_0, eta_0),
                    Bin(pt_1, eta_1));
      if ((*el_charge)[0] == (*el_charge)[1]){  // same sign
      ssEvents.put(Bin(pt_0, eta_0),
                    Bin(pt_1, eta_1));
      }
    } else {
      allEvents.put(Bin(pt_1, eta_1),
                    Bin(pt_0, eta_0));
      if ((*el_charge)[0] == (*el_charge)[1]){
      ssEvents.put(Bin(pt_1, eta_1),
                    Bin(pt_0, eta_0));
      }
    }

    // MCTruth
    b_el_true_pdg->GetEntry(tentry);
    float true_charge = 0;

    for (Long64_t j = 0; j < 2; j++)
    {
      pt_0 = (*el_pt)[j];
      eta_0 = (*el_eta)[j];
      if (el_true_pdg->at(j)!=11 && el_true_pdg->at(j)!=-11 )
      {
        continue;
      }
      true_charge = el_true_pdg->at(j) == 11 ? -1 : 1;
      mc_total.put(dummy,Bin(pt_0,eta_0));
      if (true_charge != el_charge->at(j))
      {
        mc_flip.put(dummy,Bin(pt_0,eta_0));
      }
    }
  }

  cout << ssEvents << endl;

  t1->ResetBranchAddresses();

  myFile->Close();

  minimize();

  // output truth
  cout << "MC TRUTH ----------\n";
  for (int i = 0; i < Bin::numPTBins ; i++) {
    for (int j = 0; j < Bin::numEtaBins ; j++) {
      Bin b(i,j);
      if (mc_total.get(dummy,b) == 0) {
        cout << "0\t";
      } else {
      cout << (float)mc_flip.get(dummy,b)/ (float)mc_total.get(dummy,b) << "\t";
      }
    }
    cout << endl;
  }

}

int main(int argc, char **argv) {
  likelihood(argv[1]);
  return 0;
}