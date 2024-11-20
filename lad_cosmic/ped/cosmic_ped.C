#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TTree.h>
#include <iostream>

using namespace ROOT;
using namespace std;

const static int N_BARS = 10;

const static int N_DATA_MAX = 100;
const static int AMP_N_BINS = 100;
const static int AMP_MIN    = 0;
const static int AMP_MAX    = 1000;

const static int PED_LOWER = 30;
const static int PED_UPPER = 50;

double **getModes(TTree *tree, int bar_num, bool isTop) {

  string ped_branch_name_top = "L.hod.000.TopAdcPed";
  string ped_branch_name_btm = "L.hod.000.BtmAdcPed";

  string bar_branch_name_top = "L.hod.000.TopAdcCounter";
  string bar_branch_name_btm = "L.hod.000.BtmAdcCounter";

  double peak_peds_top[N_DATA_MAX] = {0};
  double peak_peds_btm[N_DATA_MAX] = {0};
  double hit_bars_top[N_DATA_MAX]  = {0};
  double hit_bars_btm[N_DATA_MAX]  = {0};

  int n_data_btm = 0;
  int n_data_top = 0;

  tree->SetBranchAddress("Ndata.L.hod.000.TopAdcCounter", &n_data_top);
  tree->SetBranchAddress("Ndata.L.hod.000.BtmAdcCounter", &n_data_btm);
  tree->SetBranchAddress(ped_branch_name_top.c_str(), peak_peds_top);
  tree->SetBranchAddress(ped_branch_name_btm.c_str(), peak_peds_btm);

  tree->SetBranchAddress(bar_branch_name_top.c_str(), hit_bars_top);
  tree->SetBranchAddress(bar_branch_name_btm.c_str(), hit_bars_btm);

  TH1F h1_amp_top[N_BARS - 2];
  TH1F h1_amp_btm[N_BARS - 2];
  for (int i = 0; i < N_BARS - 2; i++) {
    h1_amp_top[i] = TH1F(Form("Pulse_Ped_Top_%d", i), Form("Pulse_Ped_%d; Pulse Ped [mV]; Counts", i), AMP_N_BINS,
                         AMP_MIN, AMP_MAX);
    h1_amp_btm[i] = TH1F(Form("Pulse_Ped_Btm_%d", i), Form("Pulse_Ped_%d; Pulse Ped [mV]; Counts", i), AMP_N_BINS,
                         AMP_MIN, AMP_MAX);
  }

  for (int i_evt = 0; i_evt < tree->GetEntries(); i_evt++) {
    tree->GetEntry(i_evt);

    
    for (int i = 0; i < n_data_top; i++) {
      h1_amp_top[int(hit_bars_top[i])-1].Fill(peak_peds_top[i]);
    }

    for (int i = 0; i < n_data_btm; i++) {
      h1_amp_btm[int(hit_bars_btm[i])-1].Fill(peak_peds_btm[i]);
    }
  }

  double **ped_medians = new double *[2];
  for (int i = 0; i < 2; ++i) {
    ped_medians[i] = new double[N_BARS - 2];
  }

  for (int i = 0; i < N_BARS - 2; i++) {
    double max_amp = h1_amp_top[i].GetMean();
    ped_medians[0][i] = max_amp;

    max_amp = h1_amp_btm[i].GetMean();
    ped_medians[1][i] = max_amp;
  }

  // cout << max_peaks[0][0] << " " << max_peaks[0][1] << " " << max_peaks[0][2] << " " << max_peaks[0][3] << " "
  //      << max_peaks[0][4] << " " << max_peaks[0][5] << " " << max_peaks[0][6] << " " << max_peaks[0][7] << " " <<
  //      endl;
  // cout << ped_medians[0][0] << " " << ped_medians[0][1] << " " << ped_medians[0][2] << " " << ped_medians[0][3] << " "
  //      << ped_medians[0][4] << " " << ped_medians[0][5] << " " << ped_medians[0][6] << " " << ped_medians[0][7] << " "
  //      << ped_medians[0][8] << " " << ped_medians[0][9] << endl;
  return ped_medians;
}

void cosmic_ped() {
  const int n_runs                    = 5;
  vector<int> run_numbers             = {22, 65, 67, 92, 93};
  double v_init                       = 1000;
  vector<double> voltage_offsets      = {0, 0, -100, -100, -300};
  string path                         = "/volatile/hallc/c-lad/ehingerl/ROOTfiles/COSMICS/";
  double peak_amp_top[N_BARS][n_runs] = {0};
  double peak_amp_btm[N_BARS][n_runs] = {0};

  vector<map<int, int>> map_vector(n_runs);
  for (int i = 0; i < n_runs; i++) {
    for (int j = 0; j < N_BARS - 2; j++) {
      map_vector[i][j] = j + 1;
    }
  }

  for (int i = 0; i < n_runs; i++) {
    int run_number  = run_numbers[i];
    string run_file = path + "LAD_cosmic_" + to_string(run_number) + "_-1.root";
    TFile *file     = new TFile(run_file.c_str());
    TTree *tree     = dynamic_cast<TTree *>(file->Get("T"));

    double **max_peaks = getModes(tree, i, true);
    for (int j = 0; j < N_BARS - 2; j++) {
      peak_amp_top[map_vector.at(i)[j]][i] = max_peaks[0][j];
      peak_amp_btm[map_vector.at(i)[j]][i] = max_peaks[1][j];
    }
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
      graphs[bar_num][i]->SetPoint(0, v_init + voltage_offsets[i], peak_amp_top[bar_num][i]);
      // cout << "Bar " << bar_num + 1 << " Run " << i << " " << peak_amp_top[bar_num][i] << endl;

      graphs[bar_num][i]->SetMarkerColor(i + 1);
      // graphs[bar_num][i]->SetMarkerColor(kBlack);
      graphs[bar_num][i]->SetMarkerStyle(20);

      graphs[bar_num][i]->SetTitle(Form("Bar %d; Voltage [mV]; Peak Amplitude [mV]", bar_num));
      double min_voltage = v_init + *min_element(voltage_offsets.begin(), voltage_offsets.end()) - 50;
      double max_voltage = v_init + *max_element(voltage_offsets.begin(), voltage_offsets.end()) + 50;
      graphs[bar_num][i]->GetXaxis()->SetLimits(min_voltage, max_voltage);
      graphs[bar_num][i]->GetYaxis()->SetRangeUser(PED_LOWER, PED_UPPER);

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
  title->SetBorderSize(0); // Remove the border
  title->SetShadowColor(0); // Remove the shadow
  title->SetTextAlign(22); // Center align the text
  title->SetTextSize(0.04); // Adjust text size to make it more prominent
  title->Draw();
  

  c1->SaveAs("peak_ped_top.pdf");

  TCanvas *c2 = new TCanvas("c2", "Peak Amplitude Bottom", 800, 600);
  TGraph *graphs_btm[N_BARS][n_runs];

  c2->Divide(2, 5);                                  // Divide canvas into 2 columns and 5 rows
  TLegend *legend_btm = new TLegend(0.1, 0.3, 0.2, 0.9); // Create a legend
  legend_btm->SetHeader("Legend", "C");                  // Set legend header 
  for (int bar_num = 0; bar_num < N_BARS; bar_num++) {
    c2->cd(bar_num + 1); // Move to the next pad

    for (int i = 0; i < n_runs; i++) {
      graphs_btm[bar_num][i] = new TGraph(1);
      graphs_btm[bar_num][i]->SetPoint(0, v_init + voltage_offsets[i], peak_amp_btm[bar_num][i]);
      // cout << "Bar " << bar_num + 1 << " Run " << i << " " << peak_amp_btm[bar_num][i] << endl;

      graphs_btm[bar_num][i]->SetMarkerColor(i + 1);
      // graphs_btm[bar_num][i]->SetMarkerColor(kBlack);
      graphs_btm[bar_num][i]->SetMarkerStyle(20);

      graphs_btm[bar_num][i]->SetTitle(Form("Bar %d; Voltage [V]; Peak Amplitude [mV]", bar_num));
      double min_voltage = v_init + *min_element(voltage_offsets.begin(), voltage_offsets.end()) - 50;
      double max_voltage = v_init + *max_element(voltage_offsets.begin(), voltage_offsets.end()) + 50;
      graphs_btm[bar_num][i]->GetXaxis()->SetLimits(min_voltage, max_voltage);
      graphs_btm[bar_num][i]->GetYaxis()->SetRangeUser(PED_LOWER, PED_UPPER);

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
  }
  // Draw the legend on the first pad
  c2->cd(1);
  c2->Update();

  c2->cd();
  TPaveText *title_btm = new TPaveText(0.3, 0.95, 0.7, 1, "brNDC");
  title_btm->AddText("Peak Amplitude Bottom");
  title_btm->SetFillStyle(0); // Remove the fill
  // title_btm->SetFillColor(0);
  title_btm->SetBorderSize(0); // Remove the border
  title_btm->SetShadowColor(0); // Remove the shadow
  title_btm->SetTextAlign(22); // Center align the text
  title_btm->SetTextSize(0.04); // Adjust text size to make it more prominent
  title_btm->Draw();

  c2->SaveAs("peak_ped_btm.pdf");
}