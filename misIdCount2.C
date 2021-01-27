// Usage:
// g++ -o misIdCount2.exe misIdCount2.C `root-config --cflags --libs`
// ./misIdCount2.exe
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <iostream>
#include <vector>
#include <string>
#ifdef __MAKECINT__
#pragma link C++ class vector < float> + ;
#endif

using namespace std;

void misIdCount2(char *fileName)
{

  Long64_t totalE = 0;
  Long64_t totalAnti = 0;
  Long64_t tEToAnti = 0;
  Long64_t tAntiToE = 0;

  TFile *myFile = TFile::Open(fileName);
  // Create chain of root trees
  // TTreeReader myReader("nominal_Loose;724", myFile);
  TTree *t1 = (TTree *)myFile->Get("nominal_Loose;724");

  // TTreeReader reader("nominal_Loose;724", myFile);
  // t1->Print();

  vector<int> *el_true_pdg = 0;
  TBranch *b_el_true_pdg = 0;
  vector<int> *el_charge = 0;
  TBranch *b_el_charge = 0;
  t1->ResetBranchAddresses();
  t1->SetBranchAddress("el_true_pdg", &el_true_pdg, &b_el_true_pdg);
  t1->SetBranchAddress("el_charge", &el_charge, &b_el_charge);

  Long64_t nentries = t1->GetEntries();
  for (Long64_t i = 0; i < 10000; i++)
  {
    Long64_t tentry = t1->LoadTree(i);
    b_el_true_pdg->GetEntry(tentry);
    b_el_charge->GetEntry(tentry);

    float true_charge = 0;
    if (el_true_pdg->size()!=2)
    {
      cout << el_true_pdg->size() << endl;
    }
    for (Long64_t j = 0; j < el_true_pdg->size(); j++)
    {
      if (el_true_pdg->at(j)!=11 && el_true_pdg->at(j)!=-11 )
      {
        // cout << i << "\t" << el_true_pdg->at(j) << endl;
        continue;
      }
      true_charge = el_true_pdg->at(j) == 11 ? -1 : 1;
      if (true_charge == 1)
      {
        totalAnti++;
      }
      else
      {
        totalE++;
      }
      if (true_charge != el_charge->at(j))
      {
        // cout << el_true_pdg->at(j) << "\t" << el_charge->at(j) << endl;
        if (true_charge == 1)
        {
          tAntiToE++;
        }
        else
        {
          tEToAnti++;
        }
      }
    }
  }
  cout << totalE << "\t" << totalAnti << "\t" << tAntiToE << "\t" << tEToAnti;
  t1->ResetBranchAddresses();
}

int main(int argc, char **argv)
{
  misIdCount2(argv[1]);
  return 0;
}