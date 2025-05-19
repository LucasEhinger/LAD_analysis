#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>

const int samp2ns    = 24;
const int offset     = 4572;
const int gem_window = 6 * samp2ns;
const int rebin      = 2; // Rebin factor for histograms

void make_gem_window_plots() {
  gROOT->SetBatch(kTRUE); // Enable batch mode to suppress graphical output
  std::map<int, int> runLatencyMap = {
      {22562, 119}, {22563, 119}, {22564, 115}, {22565, 121}, {22566, 113}, {22567, 123},
      // Add more run numbers and latencies as needed
  };

  // Vector of histogram names
  std::vector<std::string> histogramNames = {"time_wTrack/h_all_time_top",
                                             "time_wTrack/h_all_time_btm",
                                             "all_times/h_all_times_top",
                                             "all_times/h_all_times_btm",
                                             "time_ratio/h_time_ratio_top_all_planes",
                                             "time_ratio/h_time_ratio_btm_all_planes"};

  // Create a single PDF file to save all pages
  std::string pdfFileName = "gem_window_summary.pdf";
  bool firstPage          = true;

  // Loop through the map of run numbers to latencies
  for (const auto &entry : runLatencyMap) {
    int runNumber = entry.first;
    int latency   = entry.second;

    // Construct the ROOT file name based on the run number
    std::string fileName = "gem_window_" + std::to_string(runNumber) + "_-1_P.root";

    // Open the ROOT file
    TFile *file = TFile::Open(fileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open file " << fileName << std::endl;
      continue;
    }

    // Create a blank page with header information
    TCanvas headerCanvas("headerCanvas", "Header Page", 1200, 800);
    TPaveText header(0.1, 0.4, 0.9, 0.6, "NDC");
    header.AddText(("Run: " + std::to_string(runNumber) + ", Latency: " + std::to_string(latency)).c_str());
    header.SetFillColor(0);
    header.SetTextAlign(22); // Center align
    header.Draw();
    headerCanvas.SaveAs((pdfFileName + (firstPage ? "(" : "")).c_str());
    firstPage = false;

    // Create a canvas for drawing histograms
    TCanvas canvas("canvas", "GEM Window Plots", 1200, 800);
    canvas.Divide(1, 2); // Divide canvas into a 1x2 grid for histograms

    TLine line1[2], line2[2];

    int padIndex = 1;

    // Loop through the histogram names
    for (const auto &histName : histogramNames) {
      // Retrieve the histogram from the file
      TH1 *hist = dynamic_cast<TH1 *>(file->Get(histName.c_str()));
      if (!hist) {
        std::cerr << "Warning: Histogram " << histName << " not found in file " << fileName << std::endl;
        continue;
      }

      // Draw the histogram on the corresponding pad
      canvas.cd(padIndex);
      hist->Rebin(rebin);
      hist->Draw();

      // Calculate the positions for the vertical lines
      double line1Pos = offset - latency * samp2ns;
      double line2Pos = offset - latency * samp2ns + gem_window;

      // Draw the first vertical line
      line1[padIndex - 1] = TLine(line1Pos, hist->GetMinimum(), line1Pos, hist->GetMaximum());
      line1[padIndex - 1].SetLineColor(kRed);
      line1[padIndex - 1].SetLineStyle(2); // Dashed line
      line1[padIndex - 1].Draw("same");

      // Draw the second vertical line
      line2[padIndex - 1] = TLine(line2Pos, hist->GetMinimum(), line2Pos, hist->GetMaximum());
      line2[padIndex - 1].SetLineColor(kRed);
      line2[padIndex - 1].SetLineStyle(2); // Dashed line
      line2[padIndex - 1].Draw("same");

      padIndex++;

      // Save the current canvas page to the PDF if both histograms are drawn
      if (padIndex > 2 || &histName == &histogramNames.back()) {
        canvas.SaveAs(pdfFileName.c_str());

        // Reset for the next page if needed
        if (&histName != &histogramNames.back()) {
          canvas.Clear();
          canvas.Divide(1, 2);
          padIndex = 1;
        }
      }
    }

    // Close the ROOT file
    file->Close();
    delete file;
  }

  // Close the PDF file
  TCanvas dummyCanvas;
  dummyCanvas.SaveAs((pdfFileName + ")").c_str());
}