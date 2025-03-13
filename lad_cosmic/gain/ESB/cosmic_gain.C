#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TTree.h>
#include <fstream>
#include <iostream>

using namespace ROOT;
using namespace std;

const static int N_BARS = 10;

const static int N_DATA_MAX = 100;
const static int AMP_N_BINS = 300;
const static int AMP_MIN    = 0;
const static int AMP_MAX    = 3000;

const static int AMP_MAX_LOWER = 50;
const static int AMP_MAX_UPPER = 1000;

const static bool use_median   = true;
const static double TARGET_ADC = 100;

double **getModes(TFile *file) {

  TDirectory *subdir1 = file->GetDirectory("Pulse_Amp");
  TDirectory *subdir2 = subdir1->GetDirectory("Event_Amp");

  string hodo_name_top = "Top_Evt_Amp";
  string hodo_name_btm = "Btm_Evt_Amp";

  int n_data_btm = 0;
  int n_data_top = 0;

  double **max_peaks = new double *[2];
  for (int i = 0; i < 2; ++i) {
    max_peaks[i] = new double[N_BARS - 2];
  }

  for (int i = 0; i < N_BARS - 2; i++) {
    TH1F *h1_amp_top = dynamic_cast<TH1F *>(subdir2->Get(Form("%s_%d", hodo_name_top.c_str(), i)));
    TH1F *h1_amp_btm = dynamic_cast<TH1F *>(subdir2->Get(Form("%s_%d", hodo_name_btm.c_str(), i)));
    if (use_median) {
      double median, q;
      q = 0.5;
      h1_amp_top->GetQuantiles(1, &median, &q);
      if (median < AMP_MAX_LOWER)
        median = AMP_MAX_LOWER;
      if (median > AMP_MAX_UPPER)
        median = AMP_MAX_UPPER;
      max_peaks[0][i] = median;

      h1_amp_btm->GetQuantiles(1, &median, &q);
      if (median < AMP_MAX_LOWER)
        median = AMP_MAX_LOWER;
      if (median > AMP_MAX_UPPER)
        median = AMP_MAX_UPPER;
      max_peaks[1][i] = median;
    } else {
      int max_bin    = h1_amp_top->GetMaximumBin();
      double max_amp = h1_amp_top->GetBinCenter(max_bin);
      if (max_amp < AMP_MAX_LOWER)
        max_amp = AMP_MAX_LOWER;
      if (max_amp > AMP_MAX_UPPER)
        max_amp = AMP_MAX_UPPER;
      max_peaks[0][i] = max_amp;

      max_bin = h1_amp_btm->GetMaximumBin();
      max_amp = h1_amp_btm->GetBinCenter(max_bin);
      if (max_amp < AMP_MAX_LOWER)
        max_amp = AMP_MAX_LOWER;
      if (max_amp > AMP_MAX_UPPER)
        max_amp = AMP_MAX_UPPER;
      max_peaks[1][i] = max_amp;
    }
  }

  // cout << max_peaks[0][0] << " " << max_peaks[0][1] << " " << max_peaks[0][2] << " " << max_peaks[0][3] << " "
  //      << max_peaks[0][4] << " " << max_peaks[0][5] << " " << max_peaks[0][6] << " " << max_peaks[0][7] << " " <<
  //      endl;
  // cout << max_peaks[0][0] << " " << max_peaks[0][1] << " " << max_peaks[0][2] << " " << max_peaks[0][3] << " "
  //      << max_peaks[0][4] << " " << max_peaks[0][5] << " " << max_peaks[0][6] << " " << max_peaks[0][7] << " "
  //      << max_peaks[0][8] << " " << max_peaks[0][9] << endl;
  return max_peaks;
}

double *getLinFit(double *x, double *y, int n) {
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
  double a    = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
  double b    = (sum_y - a * sum_x) / n;
  double *fit = new double[2];
  fit[0]      = a;
  fit[1]      = b;
  return fit;
}
void cosmic_gain() {
  const int n_runs                      = 8;
  vector<int> run_numbers               = {22, 65, 67, 92, 93, 94, 95, 96};
  double v_init                         = 0;
  vector<double> voltage_offsets        = {0, 0, -100, -100, -300, -200, -200, -200};
  string path                           = "/work/hallc/c-lad/ehingerl/analysis/lad_cosmic/";
  double peak_amp_top[N_BARS][n_runs]   = {0};
  double peak_amp_btm[N_BARS][n_runs]   = {0};
  double voltage_correction_top[N_BARS] = {0};
  double voltage_correction_btm[N_BARS] = {0};

  vector<map<int, int>> map_vector(n_runs);
  for (int i = 0; i < n_runs; i++) {
    for (int j = 0; j < N_BARS - 2; j++) {
      map_vector[i][j] = j + 1;
    }
  }

  for (int i = 5; i < n_runs; i++) {
    map_vector[i][0]          = 0;
    map_vector[i][N_BARS - 3] = N_BARS;
  }

  for (int i = 0; i < n_runs; i++) {
    int run_number  = run_numbers[i];
    string run_file = path + "cosmic_histos_" + to_string(run_number) + "_output.root";
    TFile *file     = new TFile(run_file.c_str());

    double **max_peaks = getModes(file);
    for (int j = 0; j < N_BARS - 2; j++) {
      peak_amp_top[map_vector.at(i)[j]][i] = max_peaks[0][j];
      peak_amp_btm[map_vector.at(i)[j]][i] = max_peaks[1][j];
    }
    cout << "Processed run " << run_number << endl;
  }

  TCanvas *c1 = new TCanvas("c1", "Peak Amplitude Top", 800, 600);
  TGraph *graphs[N_BARS][n_runs];

  c1->Divide(2, 5);                                  // Divide canvas into 2 columns and 5 rows
  TLegend *legend = new TLegend(0.1, 0.3, 0.2, 0.9); // Create a legend
  legend->SetHeader("Legend", "C");                  // Set legend header
  for (int bar_num = 0; bar_num < N_BARS; bar_num++) {
    c1->cd(bar_num + 1); // Move to the next pad

    for (int i = 0; i < n_runs; i++) {
      graphs[bar_num][i] = new TGraph(1);
      if (peak_amp_top[bar_num][i] > 0)
        graphs[bar_num][i]->SetPoint(0, v_init + voltage_offsets[i], peak_amp_top[bar_num][i]);
      // cout << "Bar " << bar_num + 1 << " Run " << i << " " << peak_amp_top[bar_num][i] << endl;

      graphs[bar_num][i]->SetMarkerColor(i + 1);
      // graphs[bar_num][i]->SetMarkerColor(kBlack);
      graphs[bar_num][i]->SetMarkerStyle(20);

      graphs[bar_num][i]->SetTitle(Form("Bar %d; Voltage [mV]; Peak Amplitude [mV]", bar_num));
      double min_voltage = v_init + *min_element(voltage_offsets.begin(), voltage_offsets.end()) - 50;
      double max_voltage = v_init + *max_element(voltage_offsets.begin(), voltage_offsets.end()) + 50;
      graphs[bar_num][i]->GetXaxis()->SetLimits(min_voltage, max_voltage);
      graphs[bar_num][i]->GetYaxis()->SetRangeUser(AMP_MAX_LOWER - 10, AMP_MAX_UPPER + 10);

      if (i == 0) {
        graphs[bar_num][i]->Draw("AP");
      } else {
        graphs[bar_num][i]->Draw("P SAME");
      }
      if (bar_num == 0) {
        // Add the graph to the legend
        legend->AddEntry(graphs[bar_num][i], Form("Run %d", run_numbers[i]), "P");
      }
    }

    double x[n_runs];
    double y[n_runs];
    int n_data = 0;
    for (int i = 0; i < n_runs; i++) {
      if (peak_amp_top[bar_num][i] > AMP_MAX_LOWER && peak_amp_top[bar_num][i] < AMP_MAX_UPPER) {
        // if (peak_amp_top[bar_num][i] < AMP_MAX_UPPER) {
        x[n_data] = v_init + voltage_offsets[i];
        y[n_data] = peak_amp_top[bar_num][i];
        n_data++;
      }
    }
    if (n_data > 1) {
      double *fit = getLinFit(x, y, n_data);
      if (isfinite(fit[0]) && isfinite(fit[1])) {
        TGraph *fit_graph = new TGraph(2, x, y);
        fit_graph->SetLineColor(kBlack);
        fit_graph->SetLineWidth(2);
        fit_graph->SetLineStyle(2); // Set line style to dashed
        double xmin = v_init + *min_element(voltage_offsets.begin(), voltage_offsets.end()) - 50;
        double xmax = v_init + *max_element(voltage_offsets.begin(), voltage_offsets.end()) + 50;
        fit_graph->SetPoint(0, xmin, fit[0] * xmin + fit[1]);
        fit_graph->SetPoint(1, xmax, fit[0] * xmax + fit[1]);
        fit_graph->Draw("L SAME");
        legend->AddEntry(fit_graph, "Fit", "L");

        voltage_correction_top[bar_num] = (TARGET_ADC - fit[1]) / fit[0];
      }
    }
  }
  // Draw the legend on the first pad
  c1->cd(1);
  legend->Draw();
  c1->Update();
  c1->cd();
  TPaveText *title = new TPaveText(0.3, 0.95, 0.7, 1, "brNDC");
  title->AddText("Peak Amplitude Top");
  title->SetFillStyle(0); // Remove the fill
  // title->SetFillColor(0);
  title->SetBorderSize(0);  // Remove the border
  title->SetShadowColor(0); // Remove the shadow
  title->SetTextAlign(22);  // Center align the text
  title->SetTextSize(0.04); // Adjust text size to make it more prominent
  title->Draw();

  if (use_median)
    c1->SaveAs("peak_amp_top_median.pdf");
  else
    c1->SaveAs("peak_amp_top_mode.pdf");

  TCanvas *c2 = new TCanvas("c2", "Peak Amplitude Bottom", 800, 600);
  TGraph *graphs_btm[N_BARS][n_runs];

  c2->Divide(2, 5);                                      // Divide canvas into 2 columns and 5 rows
  TLegend *legend_btm = new TLegend(0.1, 0.3, 0.2, 0.9); // Create a legend
  legend_btm->SetHeader("Legend", "C");                  // Set legend header
  for (int bar_num = 0; bar_num < N_BARS; bar_num++) {
    c2->cd(bar_num + 1); // Move to the next pad

    for (int i = 0; i < n_runs; i++) {
      graphs_btm[bar_num][i] = new TGraph(1);
      if (peak_amp_btm[bar_num][i] > 0)
        graphs_btm[bar_num][i]->SetPoint(0, v_init + voltage_offsets[i], peak_amp_btm[bar_num][i]);
      // cout << "Bar " << bar_num + 1 << " Run " << i << " " << peak_amp_btm[bar_num][i] << endl;

      graphs_btm[bar_num][i]->SetMarkerColor(i + 1);
      // graphs_btm[bar_num][i]->SetMarkerColor(kBlack);
      graphs_btm[bar_num][i]->SetMarkerStyle(20);

      graphs_btm[bar_num][i]->SetTitle(Form("Bar %d; Voltage [V]; Peak Amplitude [mV]", bar_num));
      double min_voltage = v_init + *min_element(voltage_offsets.begin(), voltage_offsets.end()) - 50;
      double max_voltage = v_init + *max_element(voltage_offsets.begin(), voltage_offsets.end()) + 50;
      graphs_btm[bar_num][i]->GetXaxis()->SetLimits(min_voltage, max_voltage);
      graphs_btm[bar_num][i]->GetYaxis()->SetRangeUser(AMP_MAX_LOWER - 10, AMP_MAX_UPPER + 10);

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
      if (peak_amp_btm[bar_num][i] > AMP_MAX_LOWER && peak_amp_btm[bar_num][i] < AMP_MAX_UPPER) {
        x[n_data] = v_init + voltage_offsets[i];
        y[n_data] = peak_amp_btm[bar_num][i];
        n_data++;
      }
    }
    if (n_data > 1) {
      double *fit = getLinFit(x, y, n_data);
      if (isfinite(fit[0]) && isfinite(fit[1])) {
        TGraph *fit_graph = new TGraph(2, x, y);
        fit_graph->SetLineColor(kBlack);
        fit_graph->SetLineWidth(2);
        fit_graph->SetLineStyle(2); // Set line style to dashed
        double xmin = v_init + *min_element(voltage_offsets.begin(), voltage_offsets.end()) - 50;
        double xmax = v_init + *max_element(voltage_offsets.begin(), voltage_offsets.end()) + 50;
        fit_graph->SetPoint(0, xmin, fit[0] * xmin + fit[1]);
        fit_graph->SetPoint(1, xmax, fit[0] * xmax + fit[1]);
        fit_graph->Draw("L SAME");
        legend_btm->AddEntry(fit_graph, "Fit", "L");

        voltage_correction_btm[bar_num] = (TARGET_ADC - fit[1]) / fit[0];
      }
    }
  }
  // Draw the legend on the first pad
  c2->cd(1);
  c2->Update();

  c2->cd();
  TPaveText *title_btm = new TPaveText(0.3, 0.95, 0.7, 1, "brNDC");
  title_btm->AddText("Peak Amplitude Bottom");
  title_btm->SetFillStyle(0); // Remove the fill
  // title_btm->SetFillColor(0);
  title_btm->SetBorderSize(0);  // Remove the border
  title_btm->SetShadowColor(0); // Remove the shadow
  title_btm->SetTextAlign(22);  // Center align the text
  title_btm->SetTextSize(0.04); // Adjust text size to make it more prominent
  title_btm->Draw();

  if (use_median)
    c2->SaveAs("peak_amp_btm_median.pdf");
  else
    c2->SaveAs("peak_amp_btm_mode.pdf");

  ofstream csv_file("voltage_corrections.csv");
  csv_file << "Bar,Top Voltage Correction,Bottom Voltage Correction\n";
  for (int bar_num = 0; bar_num < N_BARS; bar_num++) {
    csv_file << bar_num + 1 << "," << voltage_correction_top[bar_num] << "," << voltage_correction_btm[bar_num] << "\n";
  }
  csv_file.close();

}