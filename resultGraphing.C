// g++ -o resultGraphing.exe resultGraphing.C `root-config --cflags --libs`
// ./resultGraphing.exe
#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

static constexpr Float_t pT_bins[] = {1, 60, 90, 130, 1000};
static constexpr Float_t pT_bin_center[] = {30, 75, 110, 230};
static constexpr Float_t pT_bin_width[] = {30, 15, 20, 100};
static constexpr Float_t eta_bins[] = {0, 0.6, 1.1, 1.52, 1.7, 2.3, 2.5};

static constexpr int numPTBins = std::size(pT_bins) - 1;
static constexpr int numEtaBins = std::size(eta_bins) - 1;

static constexpr int canvas_nx = 2;
static constexpr int canvas_ny = 3;

Float_t** getArray(std::string fileName) {
  std::ifstream myfile;
  myfile.open(fileName);
  Float_t** arr;
  arr = new Float_t*[numEtaBins];
  for (int i = 0; i < numEtaBins; i++) {
    arr[i] = new Float_t[numPTBins];
  }
  std::string line;
  for (int i = 0; i < numPTBins; i++) {
    if (myfile.eof()) {
      std::cout << i;
      return 0;
    }

    std::getline(myfile, line);
    std::istringstream iss(line);
    Float_t f;
    for (int j = 0; j < numEtaBins; j++) {
      if (iss >> f) {
        arr[j][i] = f;
      } else {
        std::cout << i << " " << j;
        return 0;
      }
    }
  }
  return arr;
}

void outputArray(Float_t** arr) {
  for (int i = 0; i < numEtaBins; i++) {
    for (int j = 0; j < numPTBins; j++) {
      std::cout << arr[i][j] << "\t";
    }
    std::cout << "\n";
  }
}

void draw(TCanvas* canvas, Float_t** arr, char* title, Float_t** err) {
  canvas->SetLogx();  // log pt

  TH2* hist = new TH2F("hist", "Histogram Title", numPTBins, pT_bins,
                       numEtaBins, eta_bins);

  hist->SetTitle(title);
  hist->GetXaxis()->SetTitle("p_{T} (GeV)");
  hist->GetYaxis()->SetTitle("|#eta|");
  for (int i = 0; i < numEtaBins; i++) {
    for (int j = 0; j < numPTBins; j++) {
      hist->SetBinContent(j + 1, i + 1, arr[i][j] * 1000);
      if (err != 0) hist->SetBinError(j + 1, i + 1, err[i][j] * 1000);
    }
  }
  // auto *st = (TPaveStats*)h->FindObject("stats");

  hist->Draw("COLZ TEXT45");
  hist->GetXaxis()->CenterTitle(true);
  hist->GetYaxis()->CenterTitle(true);
  canvas->Update();
  // canvas->Modified();
}

void resultGraphing() {
  // gStyle->SetTitleFontSize(1);
  Float_t** mctruth = getArray("mctruth.tsv");
  Float_t** likelihood = getArray("likelihood.tsv");
  Float_t** errors = getArray("errors.tsv");
  if (mctruth == 0 || likelihood == 0) {
    return;
  }
  std::cout << "MC TRUTH\n";
  outputArray(mctruth);
  std::cout << "LIKELIHOOD\n";
  outputArray(likelihood);
  std::cout << "ERRORS\n";
  outputArray(errors);
  auto* canvas =
      new TCanvas("mc_likelihood", "Charge Mis-identification Rates", 700, 900);
  canvas->Divide(canvas_nx, canvas_ny, 0.02, 0.02);

  for (int i = 0; i < numEtaBins; i++) {
    canvas->cd(i + 1);
    canvas->SetLogx();
    TMultiGraph* mg = new TMultiGraph();
    auto gr1 = new TGraph(numPTBins, pT_bin_center,
                          mctruth[i]);  // gr1->SetLineColor(kBlue);
    auto gr2 =
        new TGraphErrors(numPTBins, pT_bin_center, likelihood[i], pT_bin_width,
                         errors[i]);  // gr2->SetLineColor(kRed);
    gr1->SetMarkerColor(kRed);
    gr2->SetMarkerColor(kBlue);
    mg->Add(gr1);
    mg->Add(gr2);
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << eta_bins[i] << " #leq |#eta| < "
       << eta_bins[i + 1];
    std::string str = ss.str();
    mg->SetTitle(str.c_str());
    mg->GetXaxis()->SetTitle("P_{T} (GeV)");
    mg->GetYaxis()->SetTitle("Mis-id Rate");
    mg->GetXaxis()->SetRangeUser(0, 250);
    mg->GetYaxis()->SetMaxDigits(2);
    mg->Draw("A*");
    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    auto legend = new TLegend(0.15, 0.7, 0.4, 0.85);
    legend->AddEntry(gr1, "MC Truth", "P");
    legend->AddEntry(gr2, "Likelihood", "P");
    legend->Draw();
    canvas->Update();
    canvas->Modified();
  }
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.18);
  gStyle->SetPaintTextFormat("4.2f");
  auto* canvas_mc = new TCanvas("mc_truth", "MC Truth", 800, 800);
  draw(canvas_mc, mctruth, "Charge flip rate from MC Truth (x10^{-3})", 0);

  auto* canvas_l = new TCanvas("likelihood", "Likelihood", 850, 800);
  draw(canvas_l, likelihood, "Charge flip rate from Likelihood (x10^{-3})",
       errors);
  return;
}

int main(int argc, char** argv) {
  resultGraphing();
  return 0;
}