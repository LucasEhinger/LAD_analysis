#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TROOT.h>
#include <iostream>
#include <string>
#include <vector>

void get_gem_separation() {
  // Hardcoded mapping of GEM separation to run numbers
  std::vector<std::pair<int, double>> run_to_separation = {{1, 20},    {2, 10},    {3, 15},    {4, 18},   {5, 17},
                                                           {6, 18},    {7, 19},    {7, 19},    {8, 22},   {9, 24},
                                                           {10, 18.3}, {11, 18.6}, {12, 17.7}, {13, 17.4}};

  // Set ROOT to batch mode
  gROOT->SetBatch(kTRUE);

  // Adjustable threshold percentage
  const double threshold_percent = 40.0; // 20% of the maximum height

  // Vectors to store results
  std::vector<double> separations;
  std::vector<double> peak_widths;

  for (const auto &pair : run_to_separation) {
    int run_number    = pair.first;
    double separation = pair.second;

    // Construct file name
    std::string file_name =
        "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_296_-1_setting" +
        std::to_string(run_number) + "_filtered.root";

    // Open the file
    TFile *file = TFile::Open(file_name.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Error opening file: " << file_name << std::endl;
      continue;
    }

    // Get the histogram
    TH1 *hist = dynamic_cast<TH1 *>(file->Get("hist_delta_pos_trans"));
    if (!hist) {
      std::cerr << "Error: Histogram 'hits_delta_pos_trans' not found in " << file_name << std::endl;
      file->Close();
      continue;
    }

    // Clone the histogram to avoid modifying the original
    TH1 *hist_clone = dynamic_cast<TH1 *>(hist->Clone());
    hist_clone->SetDirectory(nullptr);

    // Calculate the threshold as a percentage of the maximum height
    double max_height = hist_clone->GetMaximum();
    double threshold  = (threshold_percent / 100.0) * max_height;

    // Find the largest continuous peak above the threshold
    int max_start_bin = -1;
    int max_end_bin   = -1;
    double max_width  = 0.0;

    int current_start_bin = -1;
    for (int bin = 1; bin <= hist_clone->GetNbinsX(); ++bin) {
      if (hist_clone->GetBinContent(bin) >= threshold) {
        if (current_start_bin == -1) {
          current_start_bin = bin; // Start of a new peak
        }
      } else {
        if (current_start_bin != -1) {
          // End of the current peak
          int current_end_bin = bin - 1;
          double current_width =
              hist_clone->GetBinCenter(current_end_bin) - hist_clone->GetBinCenter(current_start_bin);
          if (current_width > max_width) {
            max_width     = current_width;
            max_start_bin = current_start_bin;
            max_end_bin   = current_end_bin;
          }
          current_start_bin = -1; // Reset for the next peak
        }
      }
    }

    // Handle the case where the last bin is part of the largest peak
    if (current_start_bin != -1) {
      int current_end_bin  = hist_clone->GetNbinsX();
      double current_width = hist_clone->GetBinCenter(current_end_bin) - hist_clone->GetBinCenter(current_start_bin);
      if (current_width > max_width) {
        max_width     = current_width;
        max_start_bin = current_start_bin;
        max_end_bin   = current_end_bin;
      }
    }

    // Store results
    separations.push_back(separation);
    peak_widths.push_back(max_width);

    // Clean up
    delete hist_clone;
    file->Close();
  }

  // Create a graph for peak width vs GEM separation
  TCanvas *canvas = new TCanvas("canvas", "GEM Separation Analysis", 800, 600);

  TGraph *graph_width = new TGraph(separations.size(), separations.data(), peak_widths.data());
  graph_width->SetTitle("Peak Width vs GEM Separation;GEM Separation;Peak Width");
  graph_width->SetMarkerStyle(20);

  // Draw the graph
  graph_width->Draw("AP");

  // Save the canvas
  canvas->SaveAs("gem_separation_peak_width_analysis.pdf");
}