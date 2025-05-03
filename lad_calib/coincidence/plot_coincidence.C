#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TPad.h>
#include <TString.h>
#include <TTree.h>
#include <iostream>
#include <vector>

const int TDC_NBins = 150;
const int TDC_TMIN  = -3000;
const int TDC_TMAX  = 1000;

const int ADC_NBins = 50;
const int ADC_TMIN  = -250;
const int ADC_TMAX  = 250;

static const double TDC2NS = 0.09766; // TDC to ns conversion factor
static const double ADC2NS = 0.0625;  // TDC to ADC conversion factor

void plot_coincidence(const char *inputFileName, const char *outputPdfName, int plotsPerPage = 4, bool use_raw = false) {
  // Open the ROOT file
  TFile *file = TFile::Open(inputFileName);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open file " << inputFileName << std::endl;
    return;
  }

  // Get the TTree
  TTree *tree = (TTree *)file->Get("T"); // Replace "T" with the actual tree name if different
  if (!tree) {
    std::cerr << "Error: Cannot find tree in file " << inputFileName << std::endl;
    file->Close();
    return;
  }

  // Branches and labels
  std::vector<TString> branches;
  std::vector<TString> types  = {"TdcTime", "AdcPulseTime"};
  if (use_raw) {
    for (TString &type : types) {
      type = type + "Raw";
    }
  }
  std::vector<TString> planes = {"000", "001", "100", "101", "200"};
  std::vector<TString> sides  = {"Top", "Btm"};

  for (const auto &type : types) {
    for (const auto &plane : planes) {
      for (const auto &suffix : sides) {
        branches.push_back(Form("H.ladhod.%s.%s%s", plane.Data(), suffix.Data(), type.Data()));
      }
    }
  }

  TString hmsLabel  = "HMS";
  TString shmsLabel = "SHMS";

  // Create a canvas
  TCanvas *canvas = new TCanvas("canvas", "Plots", 800, 600);
  int nRows       = std::sqrt(plotsPerPage);
  int nCols       = std::ceil((double)plotsPerPage / nRows);
  canvas->Divide(nCols, nRows);

  // Open the PDF
  TString pdfName = outputPdfName;
  canvas->Print((pdfName + "[").Data()); // Open the PDF

  int plotCount = 0;

  // Loop over branches
  for (const auto &branch : branches) {
    TString hmsBranch  = branch;
    TString shmsBranch = branch;
    hmsBranch.ReplaceAll("H.", "H.");
    shmsBranch.ReplaceAll("H.", "P.");

    // Loop over paddles (0 to 10)
    for (int paddle = 0; paddle < 11; ++paddle) {
      TString hmsBranchPaddle  = Form("%s[%d]", hmsBranch.Data(), paddle);
      TString shmsBranchPaddle = Form("%s[%d]", shmsBranch.Data(), paddle);

      // Create histograms
      int nBins = branch.Contains("TdcTime") ? TDC_NBins : ADC_NBins;
      int tMin  = branch.Contains("TdcTime") ? TDC_TMIN : ADC_TMIN;
      int tMax  = branch.Contains("TdcTime") ? TDC_TMAX : ADC_TMAX;

      if (use_raw) {
        tMax = tMax - tMin;
        tMin = 0;
      }

      TH1F *hHMS  = new TH1F("hHMS", hmsBranchPaddle, nBins, tMin, tMax);
      TH1F *hSHMS = new TH1F("hSHMS", shmsBranchPaddle, nBins, tMin, tMax);

      // Draw histograms
      if (use_raw || branch.Contains("TdcTime")) {
        TString hmsBranchScaled = Form("(%s)*%f", hmsBranchPaddle.Data(), branch.Contains("TdcTime") ? TDC2NS : ADC2NS);
        TString shmsBranchScaled = Form("(%s)*%f", shmsBranchPaddle.Data(), branch.Contains("TdcTime") ? TDC2NS : ADC2NS);
        tree->Draw(hmsBranchScaled + ">>hHMS", "", "goff");
        tree->Draw(shmsBranchScaled + ">>hSHMS", "", "goff");
      } else {
        tree->Draw(hmsBranchPaddle + ">>hHMS", "", "goff");
        tree->Draw(shmsBranchPaddle + ">>hSHMS", "", "goff");
      }

      // Determine the maximum y value for both histograms
      double maxHMS  = hHMS->GetMaximum();
      double maxSHMS = hSHMS->GetMaximum();
      double maxY    = std::max(maxHMS, maxSHMS);

      // Set the same y-axis range for both histograms
      hHMS->SetMaximum(maxY * 1.1); // Add some padding
      hSHMS->SetMaximum(maxY * 1.1);

      // Style histograms
      hHMS->SetLineColor(kRed);
      hHMS->SetLineWidth(2);
      hSHMS->SetLineColor(kBlue);
      hSHMS->SetLineWidth(2);

      // Create a legend
      TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
      legend->AddEntry(hHMS, hmsLabel, "l");
      legend->AddEntry(hSHMS, shmsLabel, "l");

      // Check if histograms have entries
      if (hHMS->GetEntries() == 0 && hSHMS->GetEntries() == 0) {
        delete hHMS;
        delete hSHMS;
        delete legend;
        continue; // Skip this paddle if both histograms are empty
      }

      // Draw on canvas
      canvas->cd((plotCount % plotsPerPage) + 1);
      hHMS->Draw();
      hSHMS->Draw("SAME");
      legend->Draw();

      TString xAxisLabel;
      if (use_raw) {
        xAxisLabel = Form("Raw Time (ns) - Paddle %d", paddle);
      } else {
        xAxisLabel = Form("Reftime Subtracted Time (ns) - Paddle %d", paddle);
      }
      hHMS->GetXaxis()->SetTitle(xAxisLabel);
      hSHMS->GetXaxis()->SetTitle(xAxisLabel);

      plotCount++;

      // Save page if full
      if (plotCount % plotsPerPage == 0) {
        canvas->Print(pdfName.Data());
        canvas->Clear();              // Clear only after saving the page
        canvas->Divide(nCols, nRows); // Re-divide the canvas
      }

      // Clean up
      delete hHMS;
      delete hSHMS;
      delete legend;
    }
  }

  // Save remaining plots
  if (plotCount % plotsPerPage != 0) {
    canvas->Print(pdfName.Data());
  }

  // Close the PDF
  canvas->Print((pdfName + "]").Data());

  // Clean up
  delete canvas;
  file->Close();
}