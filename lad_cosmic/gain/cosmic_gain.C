#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TTree.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace ROOT;
using namespace std;

bool is_laser                     = true;
const static int N_BARS           = 11;
const static int N_PLANES         = 6;
const static string plane_names[] = {"000", "001", "100", "101", "200", "REFBAR"};
const static string sub_dir       = "cosmic2"; // laser_run2
const static string fit_type      = "exponential"; // "linear" or "exponential"

const static int AMP_MAX_LOWER = 10;
const static int AMP_MAX_UPPER = 800;

const static bool use_median              = true;
const static double TARGET_ADC            = 600; // was 100 for cosmic
const static int N_POINTS_FIT             = 100;
const static int MIN_N_ENTRIES_HISTO_MODE = 100;
const static int FONT_SIZE                = -1;

const static string voltage_file_name = "pmt_voltages.csv";

int paddle_ID_to_idx(int paddle_idx) {
  switch (paddle_idx / 100) {
  case 000:
    return paddle_idx % 100;
  case 001:
    return N_BARS + paddle_idx % 100;
  case 100:
    return 2 * N_BARS + paddle_idx % 100;
  case 101:
    return 3 * N_BARS + paddle_idx % 100;
  case 200:
    return 4 * N_BARS + paddle_idx % 100;
  case 999:
    return 5 * N_BARS + paddle_idx % 100;
  default:
    return -1;
  }
}
void getPeakPlane(TFile *file, int plane_idx, int run_idx, vector<vector<double>> &amplitudes_btm,
                  vector<vector<double>> &amplitudes_top) {

  TDirectory *subdir1 = file->GetDirectory("Pulse_Amp");
  TDirectory *subdir2 = subdir1->GetDirectory("Event_Amp");
  TDirectory *subdir3 = subdir2->GetDirectory(plane_names[plane_idx].c_str());

  string hodo_name_top = "Top_Evt_Amp";
  string hodo_name_btm = "Btm_Evt_Amp";

  for (int i = 0; i < N_BARS; i++) {
    TH1F *h1_amp_top =
        dynamic_cast<TH1F *>(subdir3->Get(Form("%s_%s_%d", hodo_name_top.c_str(), plane_names[plane_idx].c_str(), i)));
    TH1F *h1_amp_btm =
        dynamic_cast<TH1F *>(subdir3->Get(Form("%s_%s_%d", hodo_name_btm.c_str(), plane_names[plane_idx].c_str(), i)));
    if (use_median) {
      double median, q;
      q = 0.5;
      h1_amp_top->GetQuantiles(1, &median, &q);
      if (median < AMP_MAX_LOWER)
        median = AMP_MAX_LOWER;
      if (median > AMP_MAX_UPPER)
        median = AMP_MAX_UPPER;
      amplitudes_top[N_BARS * plane_idx + i][run_idx] = median;

      h1_amp_btm->GetQuantiles(1, &median, &q);
      if (median < AMP_MAX_LOWER)
        median = AMP_MAX_LOWER;
      if (median > AMP_MAX_UPPER)
        median = AMP_MAX_UPPER;
      amplitudes_btm[N_BARS * plane_idx + i][run_idx] = median;
    } else { // Get maximum
      int max_bin    = h1_amp_top->GetMaximumBin();
      double max_amp = h1_amp_top->GetBinCenter(max_bin);
      if (max_amp < AMP_MAX_LOWER)
        max_amp = AMP_MAX_LOWER;
      if (max_amp > AMP_MAX_UPPER)
        max_amp = AMP_MAX_UPPER;
      amplitudes_top[N_BARS * plane_idx + i][run_idx] = max_amp;

      max_bin = h1_amp_btm->GetMaximumBin();
      max_amp = h1_amp_btm->GetBinCenter(max_bin);
      if (max_amp < AMP_MAX_LOWER)
        max_amp = AMP_MAX_LOWER;
      if (max_amp > AMP_MAX_UPPER)
        max_amp = AMP_MAX_UPPER;
      amplitudes_btm[N_BARS * plane_idx + i][run_idx] = max_amp;
    }
    if (h1_amp_top->GetEntries() < MIN_N_ENTRIES_HISTO_MODE) {
      amplitudes_top[N_BARS * plane_idx + i][run_idx] = -1;
    }
    if (h1_amp_btm->GetEntries() < MIN_N_ENTRIES_HISTO_MODE) {
      amplitudes_btm[N_BARS * plane_idx + i][run_idx] = -1;
    }
  }
}

double doLinFit(double *x, double *y, int n, double *xout, double *yout, pair<double, double> &fit_results) {
  double sum_x  = 0;
  double sum_y  = 0;
  double sum_x2 = 0;
  double sum_xy = 0;
  for (int i = 0; i < n; i++) {
    sum_x += x[i];
    sum_y += y[i];
    sum_x2 += x[i] * x[i];
    sum_xy += x[i] * y[i];
  }
  double a = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
  double b = (sum_y - a * sum_x) / n;

  // Store the fit results
  fit_results.first  = a;
  fit_results.second = b;
  // Calculate the fitted values
  for (int i = 0; i < N_POINTS_FIT; i++) {
    yout[i] = a * xout[i] + b;
  }

  double target = (TARGET_ADC - b) / a;
  if (target < 0 || isnan(target)) {
    return -1;
  }
  return target;
}

double doExpFit(double *x, double *y, int n, double *xout, double *yout, pair<double, double> &fit_results) {
  double sum_x      = 0;
  double sum_y      = 0;
  double sum_x2     = 0;
  double sum_xy     = 0;
  double sum_ln_y   = 0;
  double sum_x_ln_y = 0;

  for (int i = 0; i < n; i++) {
    sum_x += x[i];
    sum_y += y[i];
    sum_x2 += x[i] * x[i];
    sum_xy += x[i] * y[i];
    sum_ln_y += log(y[i]);
    sum_x_ln_y += x[i] * log(y[i]);
  }

  double denom = n * sum_x2 - sum_x * sum_x;
  double a     = (sum_x2 * sum_ln_y - sum_x * sum_x_ln_y) / denom;
  double b     = (n * sum_x_ln_y - sum_x * sum_ln_y) / denom;

  a = exp(a);

  // Store the fit results
  fit_results.first  = a;
  fit_results.second = b;
  // Calculate the fitted values
  for (int i = 0; i < N_POINTS_FIT; i++) {
    yout[i] = a * exp(b * xout[i]);
  }

  double target = log(TARGET_ADC / a) / b;
  if (target < 0 || isnan(target)) {
    return -1;
  }
  return target;
}

double doFit(double *x, double *y, int n, double *xout, double *yout, pair<double, double> &fit_results) {
  double target = -1;
  if (fit_type == "linear") {
    target = doLinFit(x, y, n, xout, yout, fit_results);
  } else if (fit_type == "exponential") {
    target = doExpFit(x, y, n, xout, yout, fit_results);
  } else {
    cerr << "Error: Invalid fit type. Use 'linear' or 'exponential'." << endl;
  }

  return target;
}
void cosmic_gain() {
  map<string, pair<double, double>> pmt_voltage_limits;

  // Initialize the map with min and max voltage values for each PMT
  for (int i = 0; i < N_PLANES; i++) {
    for (int j = 0; j < N_BARS; j++) {
      string pmt_top = plane_names[i] + to_string(j) + "U";
      string pmt_btm = plane_names[i] + to_string(j) + "D";

      pmt_voltage_limits[pmt_top] = make_pair(0, 2500);
      pmt_voltage_limits[pmt_btm] = make_pair(0, 2500);
    }
  }
  // Update the upper PMT voltage limits for specific bars
  pmt_voltage_limits["2007D"].second = 2100;
  pmt_voltage_limits["2008D"].second = 2100;
  pmt_voltage_limits["2003U"].second = 2100;
  
  gROOT->SetBatch(kTRUE); // Prevent ROOT from opening a GUI window
  string path = "/work/hallc/c-lad/ehingerl/analysis/lad_cosmic/general/histos/";

  ifstream voltage_file(path + "../../gain/files/" + sub_dir.c_str() + "/" + voltage_file_name);
  if (!voltage_file.is_open()) {
    cerr << "Error opening file: " << voltage_file_name << endl;
    return;
  }

  string line;
  // Read the first line to determine the number of runs (columns - 1) and get the run numbers
  getline(voltage_file, line);
  stringstream ss(line);
  string cell;
  vector<int> run_numbers;
  int n_runs = 0;
  while (getline(ss, cell, ',')) {
    if (n_runs > 0) { // Skip the first column (paddle number)
      run_numbers.push_back(stoi(cell));
    }
    n_runs++;
  }
  n_runs--; // Subtract 1 for the paddle number column

  vector<vector<double>> voltages_top(N_BARS * N_PLANES, vector<double>(n_runs, 0));
  vector<vector<double>> voltages_btm(N_BARS * N_PLANES, vector<double>(n_runs, 0));
  vector<vector<double>> amplitudes_top(N_BARS * N_PLANES, vector<double>(n_runs, 0));
  vector<vector<double>> amplitudes_btm(N_BARS * N_PLANES, vector<double>(n_runs, 0));
  vector<double> voltage_target_top(N_BARS * N_PLANES, 0);
  vector<double> voltage_target_btm(N_BARS * N_PLANES, 0);
  vector<pair<double, double>> fit_params(N_BARS * N_PLANES, make_pair(0.0, 0.0));

  int paddle_num;
  bool top;

  for (int i = 0; i < N_BARS * N_PLANES * 2; ++i) {

    getline(voltage_file, line);
    stringstream ss2(line);
    string cell2;

    // Read paddle number
    getline(ss2, cell2, ',');
    if (cell2.back() == 'U') {
      top = true;
    } else if (cell2.back() == 'D') {
      top = false;
    } else {
      cerr << "Error: Invalid paddle identifier " << cell2 << endl;
      continue;
    }
    int paddle_idx = paddle_ID_to_idx(stoi(cell2.substr(0, cell2.size() - 1)));

    // Read voltages for each run
    for (int j = 0; j < n_runs; ++j) {
      getline(ss2, cell2, ',');
      if (top) {
        voltages_top[paddle_idx][j] = stod(cell2);
      } else {
        voltages_btm[paddle_idx][j] = stod(cell2);
      }
    }
  }

  voltage_file.close();

  for (int run_idx = 0; run_idx < n_runs; run_idx++) {
    string run_file = path + "cosmic_histos_wREF_" + to_string(run_numbers[run_idx]) + "_output.root";
    TFile *file     = new TFile(run_file.c_str());
    for (int plane_idx = 0; plane_idx < N_PLANES; plane_idx++) {
      getPeakPlane(file, plane_idx, run_idx, amplitudes_btm, amplitudes_top);
    }
    file->Close();
    delete file;
  }

  TCanvas *c3 = new TCanvas("c3", "Fit Parameters", 800, 600);
  c3->Divide(2, 3); // Divide canvas into 1 column and 2 rows
  for (int plane_idx = 0; plane_idx < N_PLANES; plane_idx++) {
    TCanvas *c1 = new TCanvas(Form("c1_plane_%d", plane_idx),
                              Form("Peak Amplitude Top Plane %s", plane_names[plane_idx].c_str()), 800, 600);
    TGraph *graphs_top[N_BARS][n_runs];
    TGraph *fit_graphs_top[N_BARS];

    c1->Divide(2, 6);                                      // Divide canvas into 2 columns and 5 rows
    TLegend *legend_top = new TLegend(0.1, 0.3, 0.2, 0.9); // Create a legend
    legend_top->SetHeader("Legend", "C");                  // Set legend header
    for (int bar_num = 0; bar_num < N_BARS; bar_num++) {
      c1->cd(bar_num + 1); // Move to the next pad

      // Plot the data
      for (int i = 0; i < n_runs; i++) {
        graphs_top[bar_num][i] = new TGraph(1);
        if (amplitudes_top[N_BARS * plane_idx + bar_num][i] > 0)
          graphs_top[bar_num][i]->SetPoint(0, voltages_top[N_BARS * plane_idx + bar_num][i],
                                           amplitudes_top[N_BARS * plane_idx + bar_num][i]);

        graphs_top[bar_num][i]->SetMarkerColor(i + 1);
        graphs_top[bar_num][i]->SetMarkerStyle(20);

        graphs_top[bar_num][i]->SetTitle(Form("Bar %d; Voltage [mV]; Peak Amplitude [mV]", bar_num));
        double min_voltage = *min_element(voltages_top[N_BARS * plane_idx + bar_num].begin(),
                                          voltages_top[N_BARS * plane_idx + bar_num].end()) -
                             50;
        double max_voltage = *max_element(voltages_top[N_BARS * plane_idx + bar_num].begin(),
                                          voltages_top[N_BARS * plane_idx + bar_num].end()) +
                             50;
        double min_amplitude = min(*min_element(amplitudes_top[N_BARS * plane_idx + bar_num].begin(),
                                                amplitudes_top[N_BARS * plane_idx + bar_num].end()) -
                                       10,
                                   TARGET_ADC - 10);
        double max_amplitude = max(*max_element(amplitudes_top[N_BARS * plane_idx + bar_num].begin(),
                                                amplitudes_top[N_BARS * plane_idx + bar_num].end()) +
                                       10,
                                   TARGET_ADC + 10);
        graphs_top[bar_num][i]->GetXaxis()->SetLimits(min_voltage, max_voltage);
        graphs_top[bar_num][i]->GetYaxis()->SetRangeUser(min_amplitude, max_amplitude);

        if (FONT_SIZE > 0) {
          graphs_top[bar_num][i]->GetXaxis()->SetLabelSize(FONT_SIZE);
          graphs_top[bar_num][i]->GetYaxis()->SetLabelSize(FONT_SIZE);
        }

        if (i == 0) {
          graphs_top[bar_num][i]->Draw("AP");
        } else {
          graphs_top[bar_num][i]->Draw("P SAME");
        }
        if (bar_num == 0) {
          // Add the graph to the legend
          legend_top->AddEntry(graphs_top[bar_num][i], Form("Run %d", run_numbers[i]), "P");
        }
      }

      // Fit the data
      double x[n_runs];
      double y[n_runs];
      int n_data = 0;
      for (int i = 0; i < n_runs; i++) {
        if (amplitudes_top[N_BARS * plane_idx + bar_num][i] > AMP_MAX_LOWER &&
        amplitudes_top[N_BARS * plane_idx + bar_num][i] < AMP_MAX_UPPER &&
        voltages_top[N_BARS * plane_idx + bar_num][i] > pmt_voltage_limits[plane_names[plane_idx] + to_string(bar_num) + "U"].first &&
        voltages_top[N_BARS * plane_idx + bar_num][i] < pmt_voltage_limits[plane_names[plane_idx] + to_string(bar_num) + "U"].second) {
          x[n_data] = voltages_top[N_BARS * plane_idx + bar_num][i];
          y[n_data] = amplitudes_top[N_BARS * plane_idx + bar_num][i];
          n_data++;
        }
      }
      fit_graphs_top[bar_num] = new TGraph(N_POINTS_FIT, x, y);
      if (n_data > 1) {
        double xmin = *min_element(voltages_top[N_BARS * plane_idx + bar_num].begin(),
                                   voltages_top[N_BARS * plane_idx + bar_num].end()) -
                      50;
        double xmax = *max_element(voltages_top[N_BARS * plane_idx + bar_num].begin(),
                                   voltages_top[N_BARS * plane_idx + bar_num].end()) +
                      50;
        double x_vals[N_POINTS_FIT] = {0};
        double step                 = (xmax - xmin) / (N_POINTS_FIT - 1);
        for (int i = 0; i < N_POINTS_FIT; i++) {
          x_vals[i] = xmin + i * step;
        }
        double *y_vals = new double[N_POINTS_FIT];
        voltage_target_top[N_BARS * plane_idx + bar_num] =
            doFit(x, y, n_data, x_vals, y_vals, fit_params[N_BARS * plane_idx + bar_num]);
        if (voltage_target_top[N_BARS * plane_idx + bar_num] > 0) {
          fit_graphs_top[bar_num]->SetLineColor(kBlack);
          fit_graphs_top[bar_num]->SetLineWidth(2);
          fit_graphs_top[bar_num]->SetLineStyle(2); // Set line style to dashed
          for (int i = 0; i < N_POINTS_FIT; i++) {
            fit_graphs_top[bar_num]->SetPoint(i, x_vals[i], y_vals[i]);
          }
          fit_graphs_top[bar_num]->Draw("L SAME");
          if (bar_num == 0)
            legend_top->AddEntry(fit_graphs_top[bar_num], "Fit", "L");

          TLine *line = new TLine(xmin, TARGET_ADC, xmax, TARGET_ADC);
          line->SetLineColor(kRed);
          line->SetLineWidth(2);
          line->SetLineStyle(2); // Set line style to dashed
          line->Draw("SAME");
        }
        delete[] y_vals;
      }
    }
    // Draw the legend on the first pad
    c1->cd(1);
    legend_top->Draw();
    c1->Update();
    c1->cd();
    TPaveText *title_top = new TPaveText(0.3, 0.95, 0.7, 1, "brNDC");
    title_top->AddText(Form("Peak Amplitude Top Plane %s", plane_names[plane_idx].c_str()));
    title_top->SetFillStyle(0);   // Remove the fill
    title_top->SetBorderSize(0);  // Remove the border
    title_top->SetShadowColor(0); // Remove the shadow
    title_top->SetTextAlign(22);  // Center align the text
    title_top->SetTextSize(0.04); // Adjust text size to make it more prominent
    title_top->Draw();

    if (use_median)
      c1->SaveAs(Form("files/%s/peak_amp_top_median_plane_%s.pdf", sub_dir.c_str(), plane_names[plane_idx].c_str()));
    else
      c1->SaveAs(Form("files/%s/peak_amp_top_mode_plane_%s.pdf", sub_dir.c_str(), plane_names[plane_idx].c_str()));

    TCanvas *c2 = new TCanvas(Form("c2_plane_%d", plane_idx),
                              Form("Peak Amplitude Bottom Plane %s", plane_names[plane_idx].c_str()), 800, 600);
    TGraph *graphs_btm[N_BARS][n_runs];
    TGraph *fit_graphs_btm[N_BARS];

    c2->Divide(2, 6);                                      // Divide canvas into 2 columns and 5 rows
    TLegend *legend_btm = new TLegend(0.1, 0.3, 0.2, 0.9); // Create a legend
    legend_btm->SetHeader("Legend", "C");                  // Set legend header
    for (int bar_num = 0; bar_num < N_BARS; bar_num++) {
      c2->cd(bar_num + 1); // Move to the next pad

      for (int i = 0; i < n_runs; i++) {
        graphs_btm[bar_num][i] = new TGraph(1);
        if (amplitudes_btm[N_BARS * plane_idx + bar_num][i] > 0)
          graphs_btm[bar_num][i]->SetPoint(0, voltages_btm[N_BARS * plane_idx + bar_num][i],
                                           amplitudes_btm[N_BARS * plane_idx + bar_num][i]);

        graphs_btm[bar_num][i]->SetMarkerColor(i + 1);
        graphs_btm[bar_num][i]->SetMarkerStyle(20);

        graphs_btm[bar_num][i]->SetTitle(Form("Bar %d; Voltage [mV]; Peak Amplitude [mV]", bar_num));
        double min_voltage = *min_element(voltages_btm[N_BARS * plane_idx + bar_num].begin(),
                                          voltages_btm[N_BARS * plane_idx + bar_num].end()) -
                             50;
        double max_voltage = *max_element(voltages_btm[N_BARS * plane_idx + bar_num].begin(),
                                          voltages_btm[N_BARS * plane_idx + bar_num].end()) +
                             50;
        double min_amplitude = min(*min_element(amplitudes_btm[N_BARS * plane_idx + bar_num].begin(),
                                                amplitudes_btm[N_BARS * plane_idx + bar_num].end()) -
                                       10,
                                   TARGET_ADC - 10);
        double max_amplitude = max(*max_element(amplitudes_btm[N_BARS * plane_idx + bar_num].begin(),
                                                amplitudes_btm[N_BARS * plane_idx + bar_num].end()) +
                                       10,
                                   TARGET_ADC + 10);
        graphs_btm[bar_num][i]->GetXaxis()->SetLimits(min_voltage, max_voltage);
        graphs_btm[bar_num][i]->GetYaxis()->SetRangeUser(min_amplitude, max_amplitude);

        if (FONT_SIZE > 0) {
          graphs_btm[bar_num][i]->GetXaxis()->SetLabelSize(FONT_SIZE);
          graphs_btm[bar_num][i]->GetYaxis()->SetLabelSize(FONT_SIZE);
        }

        if (i == 0) {
          graphs_btm[bar_num][i]->Draw("AP");
        } else {
          graphs_btm[bar_num][i]->Draw("P SAME");
        }
        if (bar_num == 0) {
          // Add the graph to the legend
          legend_btm->AddEntry(graphs_btm[bar_num][i], Form("Run %d", run_numbers[i]), "P");
        }
      }

      double x[n_runs];
      double y[n_runs];
      int n_data = 0;
      for (int i = 0; i < n_runs; i++) {
        if (amplitudes_btm[N_BARS * plane_idx + bar_num][i] > AMP_MAX_LOWER &&
        amplitudes_btm[N_BARS * plane_idx + bar_num][i] < AMP_MAX_UPPER &&
        voltages_btm[N_BARS * plane_idx + bar_num][i] > pmt_voltage_limits[plane_names[plane_idx] + to_string(bar_num) + "D"].first &&
        voltages_btm[N_BARS * plane_idx + bar_num][i] < pmt_voltage_limits[plane_names[plane_idx] + to_string(bar_num) + "D"].second) {
          x[n_data] = voltages_btm[N_BARS * plane_idx + bar_num][i];
          y[n_data] = amplitudes_btm[N_BARS * plane_idx + bar_num][i];
          n_data++;
        }
      }

      fit_graphs_btm[bar_num] = new TGraph(N_POINTS_FIT, x, y);
      if (n_data > 1) {
        double xmin = *min_element(voltages_btm[N_BARS * plane_idx + bar_num].begin(),
                                   voltages_btm[N_BARS * plane_idx + bar_num].end()) -
                      50;
        double xmax = *max_element(voltages_btm[N_BARS * plane_idx + bar_num].begin(),
                                   voltages_btm[N_BARS * plane_idx + bar_num].end()) +
                      50;
        double x_vals[N_POINTS_FIT] = {0};
        double step                 = (xmax - xmin) / (N_POINTS_FIT - 1);
        for (int i = 0; i < N_POINTS_FIT; i++) {
          x_vals[i] = xmin + i * step;
        }
        double *y_vals = new double[N_POINTS_FIT];
        voltage_target_btm[N_BARS * plane_idx + bar_num] =
            doFit(x, y, n_data, x_vals, y_vals, fit_params[N_BARS * plane_idx + bar_num]);
        if (voltage_target_btm[N_BARS * plane_idx + bar_num] > 0) {
          fit_graphs_btm[bar_num]->SetLineColor(kBlack);
          fit_graphs_btm[bar_num]->SetLineWidth(2);
          fit_graphs_btm[bar_num]->SetLineStyle(2); // Set line style to dashed
          for (int i = 0; i < N_POINTS_FIT; i++) {
            fit_graphs_btm[bar_num]->SetPoint(i, x_vals[i], y_vals[i]);
          }
          fit_graphs_btm[bar_num]->Draw("L SAME");
          if (bar_num == 0)
            legend_btm->AddEntry(fit_graphs_btm[bar_num], "Fit", "L");

          TLine *line = new TLine(xmin, TARGET_ADC, xmax, TARGET_ADC);
          line->SetLineColor(kRed);
          line->SetLineWidth(2);
          line->SetLineStyle(2); // Set line style to dashed
          line->Draw("SAME");
        }
        delete[] y_vals;
      }
    }
    // Draw the legend on the first pad
    c2->cd(1);
    legend_btm->Draw();
    c2->Update();
    c2->cd();
    TPaveText *title_btm = new TPaveText(0.3, 0.95, 0.7, 1, "brNDC");
    title_btm->AddText(Form("Peak Amplitude Bottom Plane %s", plane_names[plane_idx].c_str()));
    title_btm->SetFillStyle(0);   // Remove the fill
    title_btm->SetBorderSize(0);  // Remove the border
    title_btm->SetShadowColor(0); // Remove the shadow
    title_btm->SetTextAlign(22);  // Center align the text
    title_btm->SetTextSize(0.04); // Adjust text size to make it more prominent
    title_btm->Draw();

    if (use_median)
      c2->SaveAs(Form("files/%s/peak_amp_btm_median_plane_%s.pdf", sub_dir.c_str(), plane_names[plane_idx].c_str()));
    else
      c2->SaveAs(Form("files/%s/peak_amp_btm_mode_plane_%s.pdf", sub_dir.c_str(), plane_names[plane_idx].c_str()));

    c3->cd();
    c3->cd(plane_idx + 1);       // Move to the correct pad for the fit parameters
    TGraph *fit_param_graphs[2]; // 0 for top, 1 for bottom
    fit_param_graphs[0] = new TGraph(); // Top PMT
    fit_param_graphs[1] = new TGraph(); // Bottom PMT

    int point_idx_top = 0;
    int point_idx_btm = 0;

    for (int bar_num = 0; bar_num < N_BARS; bar_num++) {
      if (fit_params[N_BARS * plane_idx + bar_num].first != 0) {
      fit_param_graphs[0]->SetPoint(point_idx_top, bar_num, fit_params[N_BARS * plane_idx + bar_num].first);
      point_idx_top++;
      }
      if (fit_params[N_BARS * plane_idx + bar_num].second != 0) {
      fit_param_graphs[1]->SetPoint(point_idx_btm, bar_num, fit_params[N_BARS * plane_idx + bar_num].second);
      point_idx_btm++;
      }
    }

    fit_param_graphs[0]->SetMarkerColor(kBlue);
    fit_param_graphs[0]->SetMarkerStyle(20);
    fit_param_graphs[0]->SetMarkerSize(1.5);
    fit_param_graphs[0]->SetLineColor(kBlue);
    fit_param_graphs[0]->SetLineWidth(2);
    fit_param_graphs[0]->SetLineStyle(1); // Set line style to solid
    fit_param_graphs[0]->SetTitle(
      Form("Fit Parameters Plane %s; Bar Number; Fit Parameter Value", plane_names[plane_idx].c_str()));

    fit_param_graphs[1]->SetMarkerColor(kRed);
    fit_param_graphs[1]->SetMarkerStyle(20);
    fit_param_graphs[1]->SetMarkerSize(1.5);
    fit_param_graphs[1]->SetLineColor(kRed);
    fit_param_graphs[1]->SetLineWidth(2);
    fit_param_graphs[1]->SetLineStyle(1); // Set line style to solid

    double min_y = numeric_limits<double>::max();
    double max_y = numeric_limits<double>::lowest();
    for (int i = 0; i < fit_param_graphs[0]->GetN(); i++) {
      double x, y;
      fit_param_graphs[0]->GetPoint(i, x, y);
      if (y < min_y) min_y = y;
      if (y > max_y) max_y = y;
    }
    for (int i = 0; i < fit_param_graphs[1]->GetN(); i++) {
      double x, y;
      fit_param_graphs[1]->GetPoint(i, x, y);
      if (y < min_y) min_y = y;
      if (y > max_y) max_y = y;
    }

    fit_param_graphs[0]->GetXaxis()->SetLimits(0, N_BARS - 1);
    fit_param_graphs[0]->GetYaxis()->SetRangeUser(min_y, max_y);

    fit_param_graphs[0]->Draw("AP");
    fit_param_graphs[1]->Draw("P SAME");

    TLegend *legend_fit_params = new TLegend(0.1, 0.8, 0.3, 0.9); // Create a legend
    legend_fit_params->SetHeader("Legend", "C");                  // Set legend header
    legend_fit_params->AddEntry(fit_param_graphs[0], "Fit Parameter Top", "P");
    legend_fit_params->AddEntry(fit_param_graphs[1], "Fit Parameter Bottom", "P");
    legend_fit_params->Draw();

    // TPaveText *title_fit_params = new TPaveText(0.3, 0.95, 0.7, 1, "brNDC");
    // title_fit_params->AddText(Form("Fit Parameters Plane %s", plane_names[plane_idx].c_str()));
    // title_fit_params->SetFillStyle(0);   // Remove the fill
    // title_fit_params->SetBorderSize(0);  // Remove the border
    // title_fit_params->SetShadowColor(0); // Remove the shadow
    // title_fit_params->SetTextAlign(22);  // Center align the text
    // title_fit_params->SetTextSize(0.04); // Adjust text size to make it more prominent
    // title_fit_params->Draw();

    c3->Update();

    // delete fit_param_graphs[0];
    // delete fit_param_graphs[1];

    delete c1;
    delete c2;
    delete title_top;
    delete title_btm;
    for (int bar_num = 0; bar_num < N_BARS; bar_num++) {
      delete fit_graphs_top[bar_num];
      // delete fit_graphs_btm[bar_num];
      for (int i = 0; i < n_runs; i++) {
        delete graphs_top[bar_num][i];
        delete graphs_btm[bar_num][i];
      }
    }
    delete legend_top;
    delete legend_btm;
  }

  c3->cd();
  // TPaveText *title_fit_params = new TPaveText(0.3, 0.95, 0.7, 1, "brNDC");
  // title_fit_params->AddText("Fit Parameters");
  // title_fit_params->SetFillStyle(0);   // Remove the fill
  // title_fit_params->SetBorderSize(0);  // Remove the border
  // title_fit_params->SetShadowColor(0); // Remove the shadow
  // title_fit_params->SetTextAlign(22);  // Center align the text
  // title_fit_params->SetTextSize(0.04); // Adjust text size to make it more prominent
  // title_fit_params->Draw();
  c3->SaveAs(Form("files/%s/fit_params.pdf", sub_dir.c_str()));
  delete c3;
  // delete title_fit_params;

  cout << "Delete Worked" << endl;

  ofstream csv_file(Form("files/%s/target_voltages_corrections.csv", sub_dir.c_str()));
  csv_file << "Bar, Voltage\n";
  for (int i = 0; i < N_PLANES; i++) {
    for (int j = 0; j < N_BARS; j++) {
      csv_file << plane_names[i] << setw(2) << setfill('0') << j << "U, " << round(voltage_target_top[N_BARS * i + j])
               << "\n";
      csv_file << plane_names[i] << setw(2) << setfill('0') << j << "D, " << round(voltage_target_btm[N_BARS * i + j])
               << "\n";
    }
  }
  csv_file.close();

  ofstream fit_file(Form("files/%s/fit_params.csv", sub_dir.c_str()));
  fit_file << "Bar, a, b\n";

  for (int i = 0; i < N_PLANES; i++) {
    for (int j = 0; j < N_BARS; j++) {
      fit_file << plane_names[i] << setw(2) << setfill('0') << j << "U, " << fit_params[N_BARS * i + j].first << ", "
               << fit_params[N_BARS * i + j].second << "\n";
      fit_file << plane_names[i] << setw(2) << setfill('0') << j << "D, " << fit_params[N_BARS * i + j].first << ", "
               << fit_params[N_BARS * i + j].second << "\n";
    }
  }
}
