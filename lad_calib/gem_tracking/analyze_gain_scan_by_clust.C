#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

enum class HistType { Integral = 0, Average = 1, Average_y = 2 };

struct hist_param {
  string name;  // Histogram name
  string yaxis; // Y-axis label
  HistType op;  // Whether to use integral or avg
  double xmin;  // Peak value
  double xmax;  // Window value
  double ymin;  // Minimum Y value
  double ymax;  // Maximum Y value
  int shared_canvas;
  hist_param(const string &n, const string &y, HistType operation, double x_min, double x_max)
      : name(n), yaxis(y), op(operation), xmin(x_min), xmax(x_max), ymin(0), ymax(1e36), shared_canvas(-1) {}
  hist_param(const string &n, const string &y, HistType operation, double x_min, double x_max, double y_min,
             double y_max)
      : name(n), yaxis(y), op(operation), xmin(x_min), xmax(x_max), ymin(y_min), ymax(y_max), shared_canvas(-1) {}
  hist_param(const string &n, const string &y, HistType operation, double x_min, double x_max, double y_min,
             double y_max, int shared)
      : name(n), yaxis(y), op(operation), xmin(x_min), xmax(x_max), ymin(y_min), ymax(y_max), shared_canvas(shared) {}
};

const double window_start         = 1625; // Start of the window
const double window_end           = 1775; // End of the window
const double GEM_dim              = 70;
const double ymin_avg             = 0.1; // Minimum Y value for average
const double ymax_avg             = 1.5; // Maximum Y value for average
const int nPlanes                 = 5;   // Number of GEM planes
const string plane_names[nPlanes] = {"000", "001", "100", "101", "200"};
const int nCuts                   = 5;                // Number of cuts to analyze
const int cut_values[nCuts]       = {1, 2, 4, 6, 10}; // dxy cuts in cm

vector<hist_param> hist_params = {
    // {"time/%s/cut_dx%d_dy%d/Times/h_hodo_time_GEM_all_x_punchthrough_cut_dx%d_dy%d_%s", "Signal Integral",
    //  HistType::Integral, window_start, window_end},
    // {"time/%s/cut_dx%d_dy%d/Times/h_hodo_time_GEM0_x_punchthrough_cut_dx%d_dy%d_%s", "Signal Integral",
    //  HistType::Integral, window_start, window_end},
    // {"time/%s/cut_dx%d_dy%d/Times/h_hodo_time_GEM1_x_punchthrough_cut_dx%d_dy%d_%s", "Signal Integral",
    //  HistType::Integral, window_start, window_end},
    // // Efficiency Bkd Sub
    {"GEM_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM0_x_eff_punchthrough_cut_dx%d_dy%d_%s",
     "GEM0 x Punchthrough Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg, 0},
    {"GEM_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM1_x_eff_punchthrough_cut_dx%d_dy%d_%s",
     "GEM1 x Punchthrough Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg, 1},
    {"GEM_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM0_y_eff_punchthrough_cut_dx%d_dy%d_%s",
     "GEM1 y Punchthrough Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg, 0},
    {"GEM_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM1_y_eff_punchthrough_cut_dx%d_dy%d_%s",
     "GEM1 y Punchthrough Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg, 1},
    {"GEM_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM0_x_eff_proton_cut_dx%d_dy%d_%s", "GEM0 x Proton Efficiency",
     HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg, 0},
    {"GEM_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM1_x_eff_proton_cut_dx%d_dy%d_%s", "GEM1 x Proton Efficiency",
     HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg, 1},
    {"GEM_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM0_y_eff_proton_cut_dx%d_dy%d_%s", "GEM1 y Proton Efficiency",
     HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg, 0},
    {"GEM_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM1_y_eff_proton_cut_dx%d_dy%d_%s", "GEM1 y Proton Efficiency",
     HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg, 1},

    // Efficiency No Bkd Sub
    // {"GEM_Efficiency/%s/cut_dx%d_dy%d/h_GEM0_x_eff_punchthrough_cut_dx%d_dy%d_%s", "GEM0 Efficiency",
    //  HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg},
    // {"GEM_Efficiency/%s/cut_dx%d_dy%d/h_GEM1_x_eff_punchthrough_cut_dx%d_dy%d_%s", "GEM1 Efficiency",
    //  HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg},
    // // Conditional Efficiency
    {"GEM_cond_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM0_x_cond_eff_punchthrough_cut_dx%d_dy%d_%s",
     "GEM0 Conditional Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg},
    {"GEM_cond_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM1_x_cond_eff_punchthrough_cut_dx%d_dy%d_%s",
     "GEM1 Conditional Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg},
    {"GEM_cond_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM0_y_cond_eff_punchthrough_cut_dx%d_dy%d_%s",
     "GEM1 y Conditional Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg},
    {"GEM_cond_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM1_y_cond_eff_punchthrough_cut_dx%d_dy%d_%s",
     "GEM1 y Conditional Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg},
    {"GEM_cond_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM0_x_cond_eff_proton_cut_dx%d_dy%d_%s",
     "GEM0 x Conditional Proton Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg},
    {"GEM_cond_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM1_x_cond_eff_proton_cut_dx%d_dy%d_%s",
     "GEM1 x Conditional Proton Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg},
    {"GEM_cond_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM0_y_cond_eff_proton_cut_dx%d_dy%d_%s",
     "GEM1 y Conditional Proton Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg},
    {"GEM_cond_Efficiency_BkdSub/%s/cut_dx%d_dy%d/h_GEM1_y_cond_eff_proton_cut_dx%d_dy%d_%s",
     "GEM1 y Conditional Proton Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg},

    // // Conditional Efficiency No Bkd Sub
    // {"GEM_cond_Efficiency/%s/cut_dx%d_dy%d/h_GEM0_x_cond_eff_punchthrough_cut_dx%d_dy%d_%s",
    //  "GEM0 Conditional Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg},
    // {"GEM_cond_Efficiency/%s/cut_dx%d_dy%d/h_GEM1_x_cond_eff_punchthrough_cut_dx%d_dy%d_%s",
    //  "GEM1 Conditional Efficiency", HistType::Average_y, -GEM_dim, GEM_dim, ymin_avg, ymax_avg},
    // // ADC/Edep Avg
    {"time/%s/cut_dx%d_dy%d/ADC_Amp_BkdSub/h_hodo_adc_amp_GEM_all_x_punchthrough_cut_dx%d_dy%d_%s",
     "Average Hodo ADC (mV)", HistType::Average, 0, 1e36},
    {"time/%s/cut_dx%d_dy%d/ADC_Amp_BkdSub/h_hodo_adc_amp_GEM0_x_punchthrough_cut_dx%d_dy%d_%s",
     "Average Hodo ADC (mV)", HistType::Average, 0, 1e36},
    {"time/%s/cut_dx%d_dy%d/ADC_Amp_BkdSub/h_hodo_adc_amp_GEM1_x_punchthrough_cut_dx%d_dy%d_%s",
     "Average Hodo ADC (mV)", HistType::Average, 0, 1e36},
    // Adc/Edep No Bkd Sub
    {"time/%s/cut_dx%d_dy%d/ADC_Amp/h_hodo_adc_amp_GEM_all_x_punchthrough_cut_dx%d_dy%d_%s", "Average Hodo ADC (mV)",
     HistType::Average, 0, 1e36},
    {"time/%s/cut_dx%d_dy%d/ADC_Amp/h_hodo_adc_amp_GEM0_x_punchthrough_cut_dx%d_dy%d_%s", "Average Hodo ADC (mV)",
     HistType::Average, 0, 1e36},
    {"time/%s/cut_dx%d_dy%d/ADC_Amp/h_hodo_adc_amp_GEM1_x_punchthrough_cut_dx%d_dy%d_%s", "Average Hodo ADC (mV)",
     HistType::Average, 0, 1e36},
    // // Time Difference Avg
    // {"time/%s/cut_dx%d_dy%d/Diff_Times_BkdSub/h_hodo_time_diff_GEM_all_x_punchthrough_cut_dx%d_dy%d_%s",
    //  "Average Time Difference", HistType::Average, window_start, window_end},
    // {"time/%s/cut_dx%d_dy%d/Diff_Times_BkdSub/h_hodo_time_diff_GEM0_x_punchthrough_cut_dx%d_dy%d_%s",
    //  "Average Time Difference", HistType::Average, window_start, window_end},
    // {"time/%s/cut_dx%d_dy%d/Diff_Times_BkdSub/h_hodo_time_diff_GEM1_x_punchthrough_cut_dx%d_dy%d_%s",
    //  "Average Time Difference", HistType::Average, window_start, window_end},
};

// Define a pair called window to hold the lower and upper bounds
pair<double, double> window        = {1625, 1775};
pair<double, double> peak          = {1750, 1780};
pair<double, double> bkd_no_window = {1850, 1950};

void analyze_gain_scan_by_clust() {
  // Set ROOT to batch mode to suppress GUI
  gROOT->SetBatch(kTRUE);
  // Vector of run numbers
  // Map of run numbers to GEM gains
  // map<int, double> run_gain_map = {{22609, 2946}, {22610, 3027}, {22611, 3109},
  //                                  {22613, 3154}, {22614, 3190}, {22615, 3232}};

  // map<int, double> run_gain_map = {{23105, 3190}, {23106, 3232}, {23107, 3272}, {23108, 3314}};

  // map<int, double> run_gain_map = {{22609, 2946}, {22610, 3027}, {22611, 3109}, {22613, 3154},
  //                                  {22614, 3190}, {22615, 3232}, {23105, 3190}, {23106, 3232},
  //                                  {23107, 3272}, {23108, 3314}, {23109, 3354}};

  map<int, double> run_gain_map = {{23105, 3190}, {23106, 3232}, {23107, 3272}, {23108, 3314},{23109, 3354}};

  char spec_prefix = 'P';
  TString file_prefix =
      "/home/ehingerl/hallc/analysis/lad_calib/gem_tracking/files/tracking_gem/tracking_gem_%d_-1_%c.root";

  // Prepare arrays of maps for plotting: [nPlanes][nCuts]
  const int nHists = hist_params.size();
  const int nRuns  = run_gain_map.size();

  // 4D vectors: [nPlanes][nCuts][nHists][run]
  vector<vector<vector<map<int, double>>>> summary_vals(
      nPlanes, vector<vector<map<int, double>>>(nCuts, vector<map<int, double>>(nHists)));

  for (const auto &run_gain : run_gain_map) {
    for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
      for (int i_cut = 0; i_cut < nCuts; ++i_cut) {
        for (int i_hist = 0; i_hist < nHists; ++i_hist) {

          int run = run_gain.first;
          // Construct file name
          string filename = Form(file_prefix, run, spec_prefix);

          // Open ROOT file
          TFile *f = TFile::Open(filename.c_str(), "READ");
          if (!f || f->IsZombie()) {
            cerr << "Could not open file: " << filename << endl;
            continue;
          }

          // Get histogram
          TH1 *h = dynamic_cast<TH1 *>(
              f->Get(Form(hist_params[i_hist].name.c_str(), plane_names[i_plane].c_str(), cut_values[i_cut],
                          cut_values[i_cut], cut_values[i_cut], cut_values[i_cut], plane_names[i_plane].c_str())));

          double int_or_avg = 0;
          // Calculate integrals
          if (hist_params[i_hist].op == HistType::Integral) {
            int_or_avg = h->Integral(h->FindBin(hist_params[i_hist].xmin), h->FindBin(hist_params[i_hist].xmax));
          } else if (hist_params[i_hist].op == HistType::Average) {
            double sum    = 0;
            double weight = 0;
            int bin_min   = h->FindBin(hist_params[i_hist].xmin);
            int bin_max   = h->FindBin(hist_params[i_hist].xmax);
            for (int bin = bin_min; bin <= bin_max; ++bin) {
              double content = h->GetBinContent(bin);
              double center  = h->GetBinCenter(bin);
              if (content > hist_params[i_hist].ymax || content < hist_params[i_hist].ymin)
                continue;
              sum += content * center; // Multiply by bin center for average
              weight += content;       // Use content as weight
            }
            if (weight > 0)
              int_or_avg = sum / weight;
            else
              int_or_avg = 0;
          }

          else if (hist_params[i_hist].op == HistType::Average_y) {
            double sum    = 0;
            double weight = 0;
            int bin_min   = h->FindBin(hist_params[i_hist].xmin);
            int bin_max   = h->FindBin(hist_params[i_hist].xmax);
            for (int bin = bin_min; bin <= bin_max; ++bin) {
              double content = h->GetBinContent(bin);
              double center  = h->GetBinCenter(bin);
              if (content > hist_params[i_hist].ymax || content < hist_params[i_hist].ymin)
                continue;
              sum += content;
              weight += 1;
            }
            if (weight > 0)
              int_or_avg = sum / weight;
            else
              int_or_avg = 0;
          }

          // Store results in summary_vals
          summary_vals[i_plane][i_cut][i_hist][run] = int_or_avg;

          f->Close();
        }
      }
    }
  }

  // Create output file
  cout << "Creating output file..." << endl;
  TFile *outFile = new TFile("files/gain_scan_clust/gain_scan_graphs_all_P.root", "RECREATE");

  for (int i_hist = 0; i_hist < nHists; ++i_hist) {
    for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
      for (int i_cut = 0; i_cut < nCuts; ++i_cut) {
        vector<double> x_vals_1, y_vals_1;
        vector<double> x_vals_2, y_vals_2;
        bool is_scan_1 = true;
        // is_scan_1=false;
        for (const auto &run_gain : run_gain_map) {
          int run     = run_gain.first;
          double gain = run_gain.second;
          auto it     = summary_vals[i_plane][i_cut][i_hist].find(run);
          if (it != summary_vals[i_plane][i_cut][i_hist].end()) {
            if (is_scan_1) {
              x_vals_1.push_back(gain);
              y_vals_1.push_back(it->second);
            } else {
              x_vals_2.push_back(gain);
              y_vals_2.push_back(it->second);
            }
            // if (run_gain.second == 3232) {
            //   is_scan_1 = false; // Skip runs outside the range
            // }
          }
        }
        if (x_vals_1.empty())
          continue;

        TGraph *gr = new TGraph(x_vals_1.size(), &x_vals_1[0], &y_vals_1[0]);
        gr->SetTitle(
            Form("%s Plane %s Cut %d", hist_params[i_hist].yaxis.c_str(), plane_names[i_plane].c_str(), i_cut));
        gr->GetXaxis()->SetTitle("GEM Gain (V)");
        gr->GetYaxis()->SetTitle(hist_params[i_hist].yaxis.c_str());
        TGraph *gr2 = nullptr;
        if (!x_vals_2.empty()) {
          gr2 = new TGraph(x_vals_2.size(), &x_vals_2[0], &y_vals_2[0]);
          gr2->SetMarkerStyle(kOpenCircle);
          gr2->SetMarkerSize(1.2);
          gr2->SetLineColor(kRed);
          gr2->SetLineWidth(2);
        }

        // TString dirName = Form("%s/%s/cut%d", hist_params[i_hist].name.c_str(), plane_names[i_plane].c_str(), i_cut);
        TString dirName = Form(hist_params[i_hist].name.c_str(), plane_names[i_plane].c_str(), cut_values[i_cut],
                               cut_values[i_cut], cut_values[i_cut], cut_values[i_cut], plane_names[i_plane].c_str());

        TString dirNameStr = dirName;
        Ssiz_t lastSlash   = dirNameStr.Last('/');
        TString histName;
        if (lastSlash != kNPOS) {
          histName = dirNameStr(lastSlash + 1, dirNameStr.Length() - lastSlash - 1);
          dirNameStr.Remove(lastSlash); // Remove last directory
        } else {
          histName   = dirNameStr;
          dirNameStr = "";
        }
        if (!outFile->GetDirectory(dirNameStr)) {
          outFile->mkdir(dirNameStr);
        }
        outFile->cd(dirNameStr);
        gr->SetTitle(histName);

        // Set graph style
        gr->SetMarkerStyle(kFullCircle);
        gr->SetMarkerSize(1.2);
        gr->SetLineWidth(2);
        // Draw graph to a canvas and save the canvas to the ROOT file
        TCanvas *c = new TCanvas(histName, histName, 800, 600);
        gr->Draw("ALP");
        if (gr2) {
          gr2->Draw("LPsame");
          // gr2->SetTitle(histName);
          // gr2->SetMarkerStyle(kOpenCircle);
          // gr2->SetMarkerSize(1.2);
          // gr2->SetLineColor(kRed);
          // gr2->SetLineWidth(2);
          // Add some spacing to the min and max values for better plot appearance
          double xmin = std::min(*std::min_element(x_vals_1.begin(), x_vals_1.end()),
                                 x_vals_2.empty() ? 1e36 : *std::min_element(x_vals_2.begin(), x_vals_2.end()));
          double xmax = std::max(*std::max_element(x_vals_1.begin(), x_vals_1.end()),
                                 x_vals_2.empty() ? -1e36 : *std::max_element(x_vals_2.begin(), x_vals_2.end()));
          double ymin = std::min(*std::min_element(y_vals_1.begin(), y_vals_1.end()),
                                 y_vals_2.empty() ? 1e36 : *std::min_element(y_vals_2.begin(), y_vals_2.end()));
          double ymax = std::max(*std::max_element(y_vals_1.begin(), y_vals_1.end()),
                                 y_vals_2.empty() ? -1e36 : *std::max_element(y_vals_2.begin(), y_vals_2.end()));

          // Add 5% margin to each side
          double xmargin = 0.05 * (xmax - xmin);
          double ymargin = 0.05 * (ymax - ymin);
          xmin -= xmargin;
          xmax += xmargin;
          ymin -= ymargin;
          ymax += ymargin;
          gr->GetXaxis()->SetLimits(xmin, xmax);
          gr->GetYaxis()->SetRangeUser(ymin, ymax);
          if (gr2) {
            gr2->GetXaxis()->SetLimits(xmin, xmax);
            gr2->GetYaxis()->SetRangeUser(ymin, ymax);
          }
        }
        c->Write(histName);
        delete c;
        outFile->cd();
      }
    }
  }

  outFile->Close();
}