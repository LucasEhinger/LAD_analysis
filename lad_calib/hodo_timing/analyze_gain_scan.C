#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// Define a pair called window to hold the lower and upper bounds
std::pair<double, double> window        = {1700, 1850};
std::pair<double, double> peak          = {1760, 1790};
std::pair<double, double> bkd_window    = {1700, 1750};
std::pair<double, double> bkd_no_window = {1850, 1950};

void analyze_gain_scan() {
  // Set ROOT to batch mode to suppress GUI
  gROOT->SetBatch(kTRUE);
  // Vector of run numbers
  // Map of run numbers to GEM gains
  std::map<int, double> run_gain_map = {
      {22609, 2946}, {22610, 3027}, {22611, 3109}, {22613, 3154},
      {22614, 3190},
      // {22615, 3232}
  };

  char spec_prefix = 'H';
  TString file_prefix =
      "/home/ehingerl/hallc/analysis/lad_calib/hodo_timing/files/hodo_timing_plots/hodo_timing_plots_%d_%c.root";

  // Name of the histogram to extract from each file
  std::string hist_name =
      "KIN/LADKIN_Time/trans_50_long_100_d0_10/h_LADKIN_Time_trans_50_long_100_d0_10_final_sum"; // <-- Edit as needed

  // Prepare vectors for plotting
  std::map<int, double> peak_int;
  std::map<int, double> bkd_window_int;
  std::map<int, double> window_int;
  std::map<int, double> bkd_no_window_int;

  for (const auto &run_gain : run_gain_map) {
    int run = run_gain.first;
    // Construct file name
    std::string filename = Form(file_prefix, run, spec_prefix);

    // Open ROOT file
    TFile *f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "Could not open file: " << filename << std::endl;
      continue;
    }

    // Get histogram
    TH1 *h = dynamic_cast<TH1 *>(f->Get(hist_name.c_str()));
    if (!h) {
      std::cerr << "Could not find histogram " << hist_name << " in " << filename << std::endl;
      f->Close();
      continue;
    }

    // Calculate integrals
    double integral_peak          = h->Integral(h->FindBin(peak.first), h->FindBin(peak.second));
    double integral_window        = h->Integral(h->FindBin(window.first), h->FindBin(window.second));
    double integral_bkd_window    = h->Integral(h->FindBin(bkd_window.first), h->FindBin(bkd_window.second));
    double integral_bkd_no_window = h->Integral(h->FindBin(bkd_no_window.first), h->FindBin(bkd_no_window.second));
    // Store results in maps
    peak_int[run]          = integral_peak;
    window_int[run]        = integral_window;
    bkd_window_int[run]    = integral_bkd_window;
    bkd_no_window_int[run] = integral_bkd_no_window;

    f->Close();
  }

  // Make graph
  TGraph *g_peak           = new TGraph();
  TGraph *g_window         = new TGraph();
  TGraph *g_bkd_window     = new TGraph();
  TGraph *g_bkd_no_window  = new TGraph();
  TGraph *g_peak_bkd_sub   = new TGraph();
  TGraph *g_window_bkd_sub = new TGraph();

  int i = 0;
  for (const auto &run_gain : run_gain_map) {
    int run     = run_gain.first;
    double gain = run_gain.second;

    g_peak->SetPoint(i, gain, peak_int[run]);
    g_window->SetPoint(i, gain, window_int[run]);
    g_bkd_window->SetPoint(i, gain, bkd_window_int[run]);
    g_bkd_no_window->SetPoint(i, gain, bkd_no_window_int[run]);
    // Calculate normalization factors for background windows
    double bkd_nowindow_width = bkd_no_window.second - bkd_no_window.first;
    double peak_width         = peak.second - peak.first;
    double window_width       = window.second - window.first;

    // Normalize background integrals to the peak width
    double bkd_no_window_peak_norm       = bkd_no_window_int[run] * (peak_width / bkd_nowindow_width);
    double bkd_no_window_for_window_norm = bkd_no_window_int[run] * (window_width / bkd_nowindow_width);

    g_peak_bkd_sub->SetPoint(i, gain, peak_int[run] - bkd_no_window_peak_norm);
    g_window_bkd_sub->SetPoint(i, gain, bkd_window_int[run] - bkd_no_window_for_window_norm);
    i++;
  }
  // Create canvases for each graph
  TCanvas *c_peak = new TCanvas("c_peak", "Peak Integral vs Gain", 800, 600);
  g_peak->SetMarkerStyle(20);
  g_peak->SetMarkerColor(kRed);
  g_peak->SetLineColor(kRed);
  g_peak->SetTitle("Peak Integral vs Gain; Gain; Integral");
  g_peak->Draw("AP");

  TCanvas *c_window = new TCanvas("c_window", "Window Integral vs Gain", 800, 600);
  g_window->SetMarkerStyle(21);
  g_window->SetMarkerColor(kBlue);
  g_window->SetLineColor(kBlue);
  g_window->SetTitle("Window Integral vs Gain; Gain; Integral");
  g_window->Draw("AP");

  TCanvas *c_bkd_window = new TCanvas("c_bkd_window", "Bkd Window Integral vs Gain", 800, 600);
  g_bkd_window->SetMarkerStyle(22);
  g_bkd_window->SetMarkerColor(kGreen);
  g_bkd_window->SetLineColor(kGreen);
  g_bkd_window->SetTitle("Bkd Window Integral vs Gain; Gain; Integral");
  g_bkd_window->Draw("AP");

  TCanvas *c_bkd_no_window = new TCanvas("c_bkd_no_window", "Bkd No Window Integral vs Gain", 800, 600);
  g_bkd_no_window->SetMarkerStyle(23);
  g_bkd_no_window->SetMarkerColor(kMagenta);
  g_bkd_no_window->SetLineColor(kMagenta);
  g_bkd_no_window->SetTitle("Bkd No Window Integral vs Gain; Gain; Integral");
  g_bkd_no_window->Draw("AP");

  TCanvas *c_peak_bkd_sub = new TCanvas("c_peak_bkd_sub", "Peak - Bkd No Window vs Gain", 800, 600);
  g_peak_bkd_sub->SetMarkerStyle(24);
  g_peak_bkd_sub->SetMarkerColor(kCyan);
  g_peak_bkd_sub->SetLineColor(kCyan);
  g_peak_bkd_sub->SetTitle("Peak - Bkd No Window vs Gain; Gain; Integral");
  g_peak_bkd_sub->Draw("AP");

  TCanvas *c_bkd_sub = new TCanvas("c_bkd_sub", "Bkd Window - Bkd No Window vs Gain", 800, 600);
  g_window_bkd_sub->SetMarkerStyle(25);
  g_window_bkd_sub->SetMarkerColor(kOrange);
  g_window_bkd_sub->SetLineColor(kOrange);
  g_window_bkd_sub->SetTitle("Bkd Window - Bkd No Window vs Gain; Gain; Integral");
  g_window_bkd_sub->Draw("AP");

  // Save all canvases and graphs to a ROOT file
  TFile *outFile = new TFile("files/gain_scan_ana/gain_scan_graphs.root", "RECREATE");
  g_peak->Write("g_peak");
  g_window->Write("g_window");
  g_bkd_window->Write("g_bkd_window");
  g_bkd_no_window->Write("g_bkd_no_window");
  g_peak_bkd_sub->Write("g_peak_bkd_sub");
  g_window_bkd_sub->Write("g_window_bkd_sub");
  c_peak->Write();
  c_window->Write();
  c_bkd_window->Write();
  c_bkd_no_window->Write();
  c_peak_bkd_sub->Write();
  c_bkd_sub->Write();
  outFile->Close();
}