#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
double epics_changes[5]       = {325, 917, 1961, 2991, 3519};
std::string setting_labels[6] = {"Original", "+1 Y", "-1 Y", "-1 X", "+1 X", "Original"};

// Helper: Convert "YYYY-MM-DD HH:MM:SS" to seconds since epoch
double parse_datetime(const std::string &datetime) {
  struct tm tm = {};
  strptime(datetime.c_str(), "%Y-%m-%d %H:%M:%S", &tm);
  return mktime(&tm);
}

void steering_study_ion_chamber() {
  std::ifstream infile("ion_chamber_rates.dat");
  if (!infile) {
    std::cerr << "Could not open file.\n";
    return;
  }

  std::string line;
  // Skip header
  std::getline(infile, line);

  std::vector<double> v_time, v_IIC3H04Pk, v_IIC3H07Pk, v_IICSHMSPk, v_XPOS, v_YPOS;

  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    std::string date, time;
    double IIC3H04Pk, IIC3H07Pk, IICSHMSPk, XPOS, YPOS;
    if (!(iss >> date >> time >> IIC3H04Pk >> IIC3H07Pk >> IICSHMSPk >> XPOS >> YPOS))
      continue;
    std::string datetime = date + " " + time;
    double t             = parse_datetime(datetime);
    v_time.push_back(t);
    v_IIC3H04Pk.push_back(IIC3H04Pk);
    v_IIC3H07Pk.push_back(IIC3H07Pk);
    v_IICSHMSPk.push_back(IICSHMSPk);
    v_XPOS.push_back(XPOS);
    v_YPOS.push_back(YPOS);
  }

  // Helper: Compute n-point moving average
  auto moving_average = [](const std::vector<double> &data, int n) {
    std::vector<double> avg;
    int N = data.size();
    avg.reserve(N);
    double sum = 0;
    for (int i = 0; i < N; ++i) {
      sum += data[i];
      if (i >= n)
        sum -= data[i - n];
      if (i >= n - 1)
        avg.push_back(sum / n);
      else
        avg.push_back(data[i]); // For first (n-1) points, just copy original
    }
    return avg;
  };

  int n = 10; // Set your moving average window size here

  std::vector<double> v_IIC3H04Pk_avg = moving_average(v_IIC3H04Pk, n);
  std::vector<double> v_IIC3H07Pk_avg = moving_average(v_IIC3H07Pk, n);
  std::vector<double> v_IICSHMSPk_avg = moving_average(v_IICSHMSPk, n);
  std::vector<double> v_XPOS_avg      = moving_average(v_XPOS, n);
  std::vector<double> v_YPOS_avg      = moving_average(v_YPOS, n);

  TFile *fout         = new TFile("ion_chamber_graphs.root", "RECREATE");
  TGraph *g_IIC3H04Pk = new TGraph(v_time.size(), v_time.data(), v_IIC3H04Pk.data());
  TGraph *g_IIC3H07Pk = new TGraph(v_time.size(), v_time.data(), v_IIC3H07Pk.data());
  TGraph *g_IICSHMSPk = new TGraph(v_time.size(), v_time.data(), v_IICSHMSPk.data());
  TGraph *g_XPOS      = new TGraph(v_time.size(), v_time.data(), v_XPOS.data());
  TGraph *g_YPOS      = new TGraph(v_time.size(), v_time.data(), v_YPOS.data());

  // Moving average graphs
  TGraph *g_IIC3H04Pk_avg = new TGraph(v_time.size(), v_time.data(), v_IIC3H04Pk_avg.data());
  TGraph *g_IIC3H07Pk_avg = new TGraph(v_time.size(), v_time.data(), v_IIC3H07Pk_avg.data());
  TGraph *g_IICSHMSPk_avg = new TGraph(v_time.size(), v_time.data(), v_IICSHMSPk_avg.data());

  g_IIC3H04Pk_avg->SetName("IIC3H04Pk_vs_time_moving_avg");
  g_IIC3H04Pk_avg->SetTitle("IIC3H04Pk Moving Average;Time;IIC3H04Pk");

  g_IIC3H07Pk_avg->SetName("IIC3H07Pk_vs_time_moving_avg");
  g_IIC3H07Pk_avg->SetTitle("IIC3H07Pk Moving Average;Time;IIC3H07Pk");

  g_IICSHMSPk_avg->SetName("IICSHMSPk_vs_time_moving_avg");
  g_IICSHMSPk_avg->SetTitle("IICSHMSPk Moving Average;Time;IICSHMSPk");

  g_IIC3H04Pk->SetName("IIC3H04Pk_vs_time");
  g_IIC3H04Pk->SetTitle("IIC3H04Pk vs Time;Time;IIC3H04Pk");

  g_IIC3H07Pk->SetName("IIC3H07Pk_vs_time");
  g_IIC3H07Pk->SetTitle("IIC3H07Pk vs Time;Time;IIC3H07Pk");

  g_IICSHMSPk->SetName("IICSHMSPk_vs_time");
  g_IICSHMSPk->SetTitle("IICSHMSPk vs Time;Time;IICSHMSPk");

  g_XPOS->SetName("XPOS_vs_time");
  g_XPOS->SetTitle("XPOS vs Time;Time;XPOS");

  g_YPOS->SetName("YPOS_vs_time");
  g_YPOS->SetTitle("YPOS vs Time;Time;YPOS");

  // Draw vertical lines at epics_changes and label regions with setting_labels
  auto draw_lines_and_labels = [&](TGraph *g, const char *yaxis = "Y") {
    double ymin = g->GetYaxis()->GetXmin();
    double ymax = g->GetYaxis()->GetXmax();
    double xmin = g->GetXaxis()->GetXmin();
    double xmax = g->GetXaxis()->GetXmax();

    // Draw vertical lines
    for (int i = 0; i < 5; ++i) {
      double x = v_time[static_cast<int>(epics_changes[i])];
      TLine *line = new TLine(x, ymin, x, ymax);
      line->SetLineColor(kRed + 2);
      line->SetLineStyle(2);
      line->Draw("same");
    }

    // Draw labels inside the graph, at the bottom
    double y_label = ymin + 0.05 * (ymax - ymin); // 5% above ymin, inside the plot
    for (int i = 0; i < 6; ++i) {
      double x_left = (i == 0) ? xmin : v_time[static_cast<int>(epics_changes[i - 1])];
      double x_right = (i == 5) ? xmax : v_time[static_cast<int>(epics_changes[i])];
      double x_label = 0.5 * (x_left + x_right);
      TLatex *latex = new TLatex(x_label, y_label, setting_labels[i].c_str());
      latex->SetTextAlign(23);
      latex->SetTextSize(0.03);
      latex->Draw("same");
    }
  };

  // Draw on a canvas for each graph, draw the graph, then draw lines/labels, then save
  TCanvas *c1 = new TCanvas("c1", "IIC3H04Pk Moving Avg", 800, 600);
  g_IIC3H04Pk_avg->Draw("AL");
  draw_lines_and_labels(g_IIC3H04Pk_avg);
  c1->Write();

  TCanvas *c2 = new TCanvas("c2", "IIC3H07Pk Moving Avg", 800, 600);
  g_IIC3H07Pk_avg->Draw("AL");
  draw_lines_and_labels(g_IIC3H07Pk_avg);
  c2->Write();

  TCanvas *c3 = new TCanvas("c3", "IICSHMSPk Moving Avg", 800, 600);
  g_IICSHMSPk_avg->Draw("AL");
  draw_lines_and_labels(g_IICSHMSPk_avg);
  c3->Write();

  TCanvas *c4 = new TCanvas("c4", "IIC3H04Pk", 800, 600);
  g_IIC3H04Pk->Draw("AL");
  draw_lines_and_labels(g_IIC3H04Pk);
  c4->Write();

  TCanvas *c5 = new TCanvas("c5", "IIC3H07Pk", 800, 600);
  g_IIC3H07Pk->Draw("AL");
  draw_lines_and_labels(g_IIC3H07Pk);
  c5->Write();

  TCanvas *c6 = new TCanvas("c6", "IICSHMSPk", 800, 600);
  g_IICSHMSPk->Draw("AL");
  draw_lines_and_labels(g_IICSHMSPk);
  c6->Write();

  TCanvas *c7 = new TCanvas("c7", "XPOS", 800, 600);
  g_XPOS->Draw("AL");
  draw_lines_and_labels(g_XPOS);
  c7->Write();

  TCanvas *c8 = new TCanvas("c8", "YPOS", 800, 600);
  g_YPOS->Draw("AL");
  draw_lines_and_labels(g_YPOS);
  c8->Write();


  fout->Close();
  std::cout << "Graphs saved to ion_chamber_graphs.root\n";
  return;
}