#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

const int nPlanes                 = 5; // Number of GEM planes
const string plane_names[nPlanes] = {"000", "001", "100", "101", "200"};
const int nCuts                   = 5; // Number of cuts to analyze
const int nHists                  = 3; // Number of histograms to analyze (peak, window, bkd_no_window)
const string hist_names[nHists]   = {"time/%s/cut%d/Times/h_hodo_time_GEM_all_x_punchthrough_cut%d_%s",
                                     "time/%s/cut%d/Times/h_hodo_time_GEM0_x_punchthrough_cut%d_%s",
                                     "time/%s/cut%d/Times/h_hodo_time_GEM1_x_punchthrough_cut%d_%s"};

// Define a pair called window to hold the lower and upper bounds
std::pair<double, double> window        = {1625, 1775};
std::pair<double, double> peak          = {1750, 1780};
std::pair<double, double> bkd_no_window = {1850, 1950};

void analyze_gain_scan_by_clust() {
  // Set ROOT to batch mode to suppress GUI
  gROOT->SetBatch(kTRUE);
  // Vector of run numbers
  // Map of run numbers to GEM gains
  std::map<int, double> run_gain_map = {{22609, 2946}, {22610, 3027}, {22611, 3109},
                                        {22613, 3154}, {22614, 3190}, {22615, 3232}};

  char spec_prefix = 'H';
  TString file_prefix =
      "/home/ehingerl/hallc/analysis/lad_calib/gem_tracking/files/tracking_gem/tracking_gem_%d_-1_%c.root";

  // Prepare arrays of maps for plotting: [nPlanes][nCuts]
  std::vector<std::vector<std::vector<std::map<int, double>>>> peak_int(
      nPlanes, std::vector<std::vector<std::map<int, double>>>(nCuts, std::vector<std::map<int, double>>(nHists)));
  std::vector<std::vector<std::vector<std::map<int, double>>>> window_int(
      nPlanes, std::vector<std::vector<std::map<int, double>>>(nCuts, std::vector<std::map<int, double>>(nHists)));
  std::vector<std::vector<std::vector<std::map<int, double>>>> bkd_no_window_int(
      nPlanes, std::vector<std::vector<std::map<int, double>>>(nCuts, std::vector<std::map<int, double>>(nHists)));

  for (const auto &run_gain : run_gain_map) {
    for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
      for (int i_cut = 0; i_cut < nCuts; ++i_cut) {
        for (int i_hist = 0; i_hist < nHists; ++i_hist) {

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
          TH1 *h = dynamic_cast<TH1 *>(f->Get(Form(hist_names[i_hist].c_str(), plane_names[i_plane].c_str(), i_cut,
                                                   i_cut, plane_names[i_plane].c_str())));
          if (!h) {
            std::cerr << "Could not find histogram " << hist_names[i_hist] << " in " << filename << std::endl;
            f->Close();
            continue;
          }

          // Calculate integrals
          double integral_peak   = h->Integral(h->FindBin(peak.first), h->FindBin(peak.second));
          double integral_window = h->Integral(h->FindBin(window.first), h->FindBin(window.second));
          double integral_bkd_no_window =
              h->Integral(h->FindBin(bkd_no_window.first), h->FindBin(bkd_no_window.second));
          // Store results in maps
          peak_int[i_plane][i_cut][i_hist][run]          = integral_peak;
          window_int[i_plane][i_cut][i_hist][run]        = integral_window;
          bkd_no_window_int[i_plane][i_cut][i_hist][run] = integral_bkd_no_window;

          f->Close();
        }
      }
    }
  }

  // Create output file
  cout << "Creating output file..." << endl;
  TFile *outFile = new TFile("files/gain_scan_clust/gain_scan_graphs.root", "RECREATE");

  // Loop over planes, cuts, and hists
  for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
    for (int i_cut = 0; i_cut < nCuts; ++i_cut) {
      for (int i_hist = 0; i_hist < nHists; ++i_hist) {
        // Create directory for plane/cut/hist
        TString dir_name = Form("plane_%s/cut_%d/hist_%d", plane_names[i_plane].c_str(), i_cut, i_hist);
        outFile->mkdir(dir_name);
        outFile->cd(dir_name);

        TGraph *g_peak           = new TGraph();
        TGraph *g_peak_bkd_sub   = new TGraph();
        TGraph *g_window         = new TGraph();
        TGraph *g_window_bkd_sub = new TGraph();
        TGraph *g_bkd_no_window  = new TGraph();

        int i = 0;
        for (const auto &run_gain : run_gain_map) {
          int run     = run_gain.first;
          double gain = run_gain.second;

          double peak_val          = peak_int[i_plane][i_cut][i_hist][run];
          double window_val        = window_int[i_plane][i_cut][i_hist][run];
          double bkd_no_window_val = bkd_no_window_int[i_plane][i_cut][i_hist][run];

          double bkd_nowindow_width = bkd_no_window.second - bkd_no_window.first;
          double peak_width         = peak.second - peak.first;
          double window_width       = window.second - window.first;

          double bkd_no_window_peak_norm       = bkd_no_window_val * (peak_width / bkd_nowindow_width);
          double bkd_no_window_for_window_norm = bkd_no_window_val * (window_width / bkd_nowindow_width);

          g_peak->SetPoint(i, gain, peak_val);
          g_window->SetPoint(i, gain, window_val);
          g_bkd_no_window->SetPoint(i, gain, bkd_no_window_val);
          g_peak_bkd_sub->SetPoint(i, gain, peak_val - bkd_no_window_peak_norm);
          g_window_bkd_sub->SetPoint(i, gain, window_val - bkd_no_window_for_window_norm);
          i++;
        }

        g_peak->SetMarkerStyle(20);
        g_window->SetMarkerStyle(20);
        g_bkd_no_window->SetMarkerStyle(20);
        g_peak_bkd_sub->SetMarkerStyle(20);
        g_window_bkd_sub->SetMarkerStyle(20);

        g_peak->Write("g_peak");
        g_window->Write("g_window");
        g_bkd_no_window->Write("g_bkd_no_window");
        g_peak_bkd_sub->Write("g_peak_bkd_sub");
        g_window_bkd_sub->Write("g_window_bkd_sub");

        // Clean up
        delete g_peak;
        delete g_window;
        delete g_bkd_no_window;
        delete g_peak_bkd_sub;
        delete g_window_bkd_sub;
      }
    }
  }
  outFile->Close();
}