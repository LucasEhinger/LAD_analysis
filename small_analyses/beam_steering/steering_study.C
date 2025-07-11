#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TTree.h>
#include <string>

double ereal_axis_max         = 500.0;
double ereal_axis_min         = 0.0;
double three_four_axis_max    = 1e3;
double three_four_axis_min    = 0.0;
double s1x_axis_max           = 1e7;
double s1x_axis_min           = 0.0;
double ibcm_axis_max          = 2.0;
double ibcm_axis_min          = 1.8;
double scaler_changes[4]      = {290, 1163, 2011, 2441};
double epics_changes[4]       = {107, 432, 737, 902};
std::string setting_labels[5] = {"+1 Y", "-1 Y", "-1 X", "+1 X", "Original"};

void DrawRawAndAvgGraphs(TFile *fout, const char *basename, TGraph *g1, TGraph *g2, const char *g1_label,
                         const char *g2_label, const double *scaler_changes, int n_changes,
                         const std::string *setting_labels, int n_labels, int nAvg = 10) {
  // Find min/max for both graphs
  double y1_min = *std::min_element(g1->GetY(), g1->GetY() + g1->GetN());
  double y1_max = *std::max_element(g1->GetY(), g1->GetY() + g1->GetN());
  double y2_min = *std::min_element(g2->GetY(), g2->GetY() + g2->GetN());
  double y2_max = *std::max_element(g2->GetY(), g2->GetY() + g2->GetN());

  // Use calculated min/max for left axis, and scale right axis accordingly
  double left_min = y1_min, left_max = y1_max;
  double right_min = y2_min, right_max = y2_max;
  // Avoid division by zero
  double scale = (right_max - right_min) != 0 ? (left_max - left_min) / (right_max - right_min) : 1.0;
  double shift = left_min - right_min * scale;

  // Draw raw graphs with left/right axes
  std::string cname_raw  = std::string("c_") + basename + "_raw";
  std::string ctitle_raw = std::string(basename) + " Raw";
  TCanvas *c_raw         = new TCanvas(cname_raw.c_str(), ctitle_raw.c_str(), 800, 600);

  // Draw g1 on left axis
  g1->GetYaxis()->SetRangeUser(left_min, left_max);
  g1->Draw("AP");

  // Draw g2 on right axis (scaled)
  std::vector<double> x2(g2->GetN()), y2_scaled(g2->GetN());
  for (int i = 0; i < g2->GetN(); ++i) {
    x2[i]        = g2->GetX()[i];
    y2_scaled[i] = g2->GetY()[i] * scale + shift;
  }
  TGraph *g2_scaled = new TGraph(g2->GetN(), x2.data(), y2_scaled.data());
  g2_scaled->SetLineColor(g2->GetLineColor());
  g2_scaled->SetMarkerColor(g2->GetMarkerColor());
  g2_scaled->SetMarkerStyle(g2->GetMarkerStyle());
  g2_scaled->Draw("P SAME");

  // Draw right axis
  double x_min       = g1->GetXaxis()->GetXmin();
  double x_max       = g1->GetXaxis()->GetXmax();
  TGaxis *axis_right = new TGaxis(x_max, left_min, x_max, left_max, right_min, right_max, 510, "+L");
  axis_right->SetLineColor(g2->GetLineColor());
  axis_right->SetLabelColor(g2->GetLineColor());
  axis_right->SetTitle(g2_label);
  axis_right->SetTitleColor(g2->GetLineColor());
  axis_right->SetTitleOffset(1.2);
  axis_right->Draw();

  // Draw vertical lines and labels
  for (int i = 0; i < n_changes; ++i) {
    double x = scaler_changes[i];
    TLine *l = new TLine(x, left_min, x, left_max);
    l->SetLineColor(kBlack);
    l->SetLineStyle(2);
    l->Draw();
  }
  for (int i = 0; i < n_labels; ++i) {
    double x_left  = (i == 0) ? 0 : scaler_changes[i - 1];
    double x_right = (i == n_changes) ? g1->GetXaxis()->GetXmax() : scaler_changes[i];
    double x_label = 0.5 * (x_left + x_right);
    double y_label = left_min + 0.05 * (left_max - left_min);
    TLatex *latex  = new TLatex(x_label, y_label, setting_labels[i].c_str());
    latex->SetTextAlign(22);
    latex->SetTextSize(0.03);
    latex->Draw();
  }
  TLegend *leg = new TLegend(0.15, 0.75, 0.35, 0.88);
  leg->AddEntry(g1, g1_label, "lp");
  leg->AddEntry(g2_scaled, g2_label, "lp");
  leg->Draw();
  c_raw->Write();

  // Compute moving averages
  std::vector<double> x1, y1, x2_orig, y2;
  for (int i = 0; i < g1->GetN(); ++i) {
    x1.push_back(g1->GetX()[i]);
    y1.push_back(g1->GetY()[i]);
  }
  for (int i = 0; i < g2->GetN(); ++i) {
    x2_orig.push_back(g2->GetX()[i]);
    y2.push_back(g2->GetY()[i]);
  }
  std::vector<double> x1_avg, y1_avg, x2_avg, y2_avg;
  for (int i = 0; i <= (int)x1.size() - nAvg; ++i) {
    double sumx = 0, sumy = 0;
    for (int j = 0; j < nAvg; ++j) {
      sumx += x1[i + j];
      sumy += y1[i + j];
    }
    x1_avg.push_back(sumx / nAvg);
    y1_avg.push_back(sumy / nAvg);
  }
  for (int i = 0; i <= (int)x2_orig.size() - nAvg; ++i) {
    double sumx = 0, sumy = 0;
    for (int j = 0; j < nAvg; ++j) {
      sumx += x2_orig[i + j];
      sumy += y2[i + j];
    }
    x2_avg.push_back(sumx / nAvg);
    y2_avg.push_back(sumy / nAvg);
  }
  // Scale averaged g2 to left axis
  std::vector<double> y2_avg_scaled(y2_avg.size());
  for (size_t i = 0; i < y2_avg.size(); ++i) {
    y2_avg_scaled[i] = y2_avg[i] * scale + shift;
  }
  TGraph *g1_avg        = new TGraph(x1_avg.size(), x1_avg.data(), y1_avg.data());
  TGraph *g2_avg_scaled = new TGraph(x2_avg.size(), x2_avg.data(), y2_avg_scaled.data());
  g1_avg->SetName((std::string(g1->GetName()) + "_avg").c_str());
  g1_avg->SetTitle((std::string(g1->GetTitle()) + " (Averaged)").c_str());
  g1_avg->SetLineColor(g1->GetLineColor() + 2);
  g1_avg->SetMarkerColor(g1->GetMarkerColor() + 2);
  g1_avg->SetMarkerStyle(24);
  g1_avg->SetLineWidth(2);
  g2_avg_scaled->SetName((std::string(g2->GetName()) + "_avg").c_str());
  g2_avg_scaled->SetTitle((std::string(g2->GetTitle()) + " (Averaged)").c_str());
  g1_avg->GetXaxis()->SetTitle(g1->GetXaxis()->GetTitle());
  g1_avg->GetYaxis()->SetTitle(g1->GetYaxis()->GetTitle());
  g2_avg_scaled->GetXaxis()->SetTitle(g2->GetXaxis()->GetTitle());
  g2_avg_scaled->GetYaxis()->SetTitle(g2->GetYaxis()->GetTitle());
  g2_avg_scaled->SetLineColor(g2->GetLineColor() + 2);
  g2_avg_scaled->SetMarkerColor(g2->GetMarkerColor() + 2);
  g2_avg_scaled->SetMarkerStyle(25);
  g2_avg_scaled->SetLineWidth(2);

  // Draw averaged graphs with left/right axes
  std::string cname_avg  = std::string("c_") + basename + "_avg";
  std::string ctitle_avg = std::string(basename) + " (Averaged)";
  TCanvas *c_avg         = new TCanvas(cname_avg.c_str(), ctitle_avg.c_str(), 800, 600);
  g1_avg->GetYaxis()->SetRangeUser(left_min, left_max);
  g1_avg->Draw("AP");
  g2_avg_scaled->Draw("P SAME");
  // Draw right axis for avg
  TGaxis *axis_right_avg = new TGaxis(x_max, left_min, x_max, left_max, right_min, right_max, 510, "+L");
  axis_right_avg->SetLineColor(g2->GetLineColor() + 2);
  axis_right_avg->SetLabelColor(g2->GetLineColor() + 2);
  axis_right_avg->SetTitle((std::string(g2_label) + " (avg)").c_str());
  axis_right_avg->SetTitleColor(g2->GetLineColor() + 2);
  axis_right_avg->SetTitleOffset(1.2);
  axis_right_avg->Draw();

  for (int i = 0; i < n_changes; ++i) {
    double x = scaler_changes[i];
    TLine *l = new TLine(x, left_min, x, left_max);
    l->SetLineColor(kBlack);
    l->SetLineStyle(2);
    l->Draw();
  }
  for (int i = 0; i < n_labels; ++i) {
    double x_left  = (i == 0) ? 0 : scaler_changes[i - 1];
    double x_right = (i == n_changes) ? g1_avg->GetXaxis()->GetXmax() : scaler_changes[i];
    double x_label = 0.5 * (x_left + x_right);
    double y_label = left_min + 0.05 * (left_max - left_min);
    TLatex *latex  = new TLatex(x_label, y_label, setting_labels[i].c_str());
    latex->SetTextAlign(22);
    latex->SetTextSize(0.03);
    latex->Draw();
  }
  TLegend *leg_avg = new TLegend(0.15, 0.75, 0.45, 0.88);
  leg_avg->AddEntry(g1_avg, (std::string(g1_label) + " (avg)").c_str(), "lp");
  leg_avg->AddEntry(g2_avg_scaled, (std::string(g2_label) + " (avg)").c_str(), "lp");
  leg_avg->Draw();
  c_avg->Write();
}

void steering_study() {
  const char *infile =
      "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/SCALAR/SCALAR_LAD_COIN_23738_0_10_-1.root";
  const char *outfile = "output.root";

  // Open input ROOT file
  TFile *fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    printf("Cannot open input file %s\n", infile);
    return;
  }

  // Get trees
  TTree *TSP = (TTree *)fin->Get("TSP");
  TTree *TSH = (TTree *)fin->Get("TSH");
  TTree *E   = (TTree *)fin->Get("E");
  if (!TSP || !TSH || !E) {
    printf("Cannot find trees TSP or TSH\n");
    fin->Close();
    return;
  }

  // H.S1X.scalerRate P.hTRIG1.scalerRate P.hTRIG3.scalerRate

  // Variables for TSP
  Double_t tsp_time, tsp_rate, tsp_34_rate, tsp_s1x_rate;
  TSP->SetBranchAddress("P.1MHz.scalerTime", &tsp_time);
  TSP->SetBranchAddress("P.pEL_REAL.scalerRate", &tsp_rate);
  TSP->SetBranchAddress("P.pTRIG1.scalerRate", &tsp_34_rate);
  TSP->SetBranchAddress("P.S1X.scalerRate", &tsp_s1x_rate);

  // Variables for TSH
  Double_t tsh_time, tsh_rate, tsh_34_rate, tsh_s1x_rate;
  TSH->SetBranchAddress("H.1MHz.scalerTime", &tsh_time);
  TSH->SetBranchAddress("H.hEL_REAL.scalerRate", &tsh_rate);
  TSH->SetBranchAddress("H.hTRIG3.scalerRate", &tsh_34_rate);
  TSH->SetBranchAddress("H.S1X.scalerRate", &tsh_s1x_rate);

  // Variables for E
  Double_t e_time, IPM3HO7A_X, IPM3HO7A_Y, IPM3HO7B_X, IPM3HO7B_Y, IPM3HO7C_X, IPM3HO7C_Y, IBCM3H04B;
  E->SetBranchAddress("IPM3H07A.XPOS", &IPM3HO7A_X);
  E->SetBranchAddress("IPM3H07A.YPOS", &IPM3HO7A_Y);
  E->SetBranchAddress("IPM3H07B.XPOS", &IPM3HO7B_X);
  E->SetBranchAddress("IPM3H07B.YPOS", &IPM3HO7B_Y);
  E->SetBranchAddress("IPM3H07C.XPOS", &IPM3HO7C_X);
  E->SetBranchAddress("IPM3H07C.YPOS", &IPM3HO7C_Y);
  E->SetBranchAddress("ibcm3H04B", &IBCM3H04B);

  // Fill TGraph for TSP
  int nTSP = TSP->GetEntries();
  std::vector<double> tsp_times, tsp_rates, tsp_34_rates, tsp_s1x_rates;
  for (int i = 0; i < nTSP; ++i) {
    TSP->GetEntry(i);
    tsp_times.push_back(tsp_time);
    tsp_rates.push_back(tsp_rate);
    tsp_34_rates.push_back(tsp_34_rate);
    tsp_s1x_rates.push_back(tsp_s1x_rate);
  }
  TGraph *gTSP = new TGraph(nTSP, tsp_times.data(), tsp_rates.data());
  gTSP->SetLineColor(kBlue);
  gTSP->SetMarkerColor(kBlue);
  gTSP->SetMarkerStyle(20);
  gTSP->SetTitle("SHMS EREAL Rate;Scaler Time (s);P.hEL_REAL Rate (Hz)");

  TGraph *gTSP_34 = new TGraph(nTSP, tsp_times.data(), tsp_34_rates.data());
  gTSP_34->SetLineColor(kBlue);
  gTSP_34->SetMarkerColor(kBlue);
  gTSP_34->SetMarkerStyle(21);
  gTSP_34->SetTitle("SHMS 3/4 Rate;Scaler Time (s);P.pTRIG1 Rate (Hz)");

  TGraph *gTSP_S1X = new TGraph(nTSP, tsp_times.data(), tsp_s1x_rates.data());
  gTSP_S1X->SetLineColor(kBlue);
  gTSP_S1X->SetMarkerColor(kBlue);
  gTSP_S1X->SetMarkerStyle(22);
  gTSP_S1X->SetTitle("SHMS S1X Rate;Scaler Time (s);P.S1X Rate (Hz)");

  // Fill TGraph for TSH
  int nTSH = TSH->GetEntries();
  std::vector<double> tsh_times, tsh_rates, tsh_34_rates, tsh_s1x_rates;
  for (int i = 0; i < nTSH; ++i) {
    TSH->GetEntry(i);
    tsh_times.push_back(tsh_time);
    tsh_rates.push_back(tsh_rate);
    tsh_34_rates.push_back(tsh_34_rate);
    tsh_s1x_rates.push_back(tsh_s1x_rate);
  }
  TGraph *gTSH = new TGraph(nTSH, tsh_times.data(), tsh_rates.data());
  gTSH->SetLineColor(kRed);
  gTSH->SetMarkerColor(kRed);
  gTSH->SetMarkerStyle(21);
  gTSH->SetTitle("HMS EREAL Rate ;Scaler Time (s);H.hEL_REAL Rate (Hz)");

  TGraph *gTSH_34 = new TGraph(nTSH, tsh_times.data(), tsh_34_rates.data());
  gTSH_34->SetLineColor(kRed);
  gTSH_34->SetMarkerColor(kRed);
  gTSH_34->SetMarkerStyle(22);
  gTSH_34->SetTitle("HMS 3/4 Rate");
  TGraph *gTSH_S1X = new TGraph(nTSH, tsh_times.data(), tsh_s1x_rates.data());
  gTSH_S1X->SetLineColor(kRed);
  gTSH_S1X->SetMarkerColor(kRed);
  gTSH_S1X->SetMarkerStyle(23);
  gTSH_S1X->SetTitle("HMS S1X Rate;Scaler Time (s);H.S1X Rate (Hz)");

  int nE = E->GetEntries();
  std::vector<double> e_times(nE), a_x(nE), a_y(nE), b_x(nE), b_y(nE), c_x(nE), c_y(nE), ibcm3h04b(nE);

  for (int i = 0; i < nE; ++i) {
    E->GetEntry(i);
    e_times[i]   = i;
    a_x[i]       = IPM3HO7A_X;
    a_y[i]       = IPM3HO7A_Y;
    b_x[i]       = IPM3HO7B_X;
    b_y[i]       = IPM3HO7B_Y;
    c_x[i]       = IPM3HO7C_X;
    c_y[i]       = IPM3HO7C_Y;
    ibcm3h04b[i] = IBCM3H04B; // Store IBCM3H04B if needed
  }

  // Create TGraphs for each branch in E
  TGraph *gA_X = new TGraph(nE, e_times.data(), a_x.data());
  gA_X->SetTitle("IPM3HO7A_X;EPICS time (index);X position (mm)");
  gA_X->GetXaxis()->SetTitle("EPICS time (index)");
  gA_X->GetYaxis()->SetTitle("X position (mm)");
  gA_X->SetLineColor(kBlue);
  gA_X->SetMarkerColor(kBlue);
  gA_X->SetMarkerStyle(22);

  TGraph *gA_Y = new TGraph(nE, e_times.data(), a_y.data());
  gA_Y->SetTitle("IPM3HO7A_Y;EPICS time (index);Y position (mm)");
  gA_Y->GetXaxis()->SetTitle("EPICS time (index)");
  gA_Y->GetYaxis()->SetTitle("Y position (mm)");
  gA_Y->SetLineColor(kRed);
  gA_Y->SetMarkerColor(kRed);
  gA_Y->SetMarkerStyle(23);

  TGraph *gB_X = new TGraph(nE, e_times.data(), b_x.data());
  gB_X->SetTitle("IPM3HO7B_X;EPICS time (index);X position (mm)");
  gB_X->GetXaxis()->SetTitle("EPICS time (index)");
  gB_X->GetYaxis()->SetTitle("X position (mm)");
  gB_X->SetLineColor(kBlue);
  gB_X->SetMarkerColor(kBlue);
  gB_X->SetMarkerStyle(24);

  TGraph *gB_Y = new TGraph(nE, e_times.data(), b_y.data());
  gB_Y->SetTitle("IPM3HO7B_Y;EPICS time (index);Y position (mm)");
  gB_Y->GetXaxis()->SetTitle("EPICS time (index)");
  gB_Y->GetYaxis()->SetTitle("Y position (mm)");
  gB_Y->SetLineColor(kRed);
  gB_Y->SetMarkerColor(kRed);
  gB_Y->SetMarkerStyle(25);

  TGraph *gC_X = new TGraph(nE, e_times.data(), c_x.data());
  gC_X->SetTitle("IPM3HO7C_X;EPICS time (index);X position (mm)");
  gC_X->GetXaxis()->SetTitle("EPICS time (index)");
  gC_X->GetYaxis()->SetTitle("X position (mm)");
  gC_X->SetLineColor(kBlue);
  gC_X->SetMarkerColor(kBlue);
  gC_X->SetMarkerStyle(26);

  TGraph *gC_Y = new TGraph(nE, e_times.data(), c_y.data());
  gC_Y->SetTitle("IPM3HO7C_Y;EPICS time (index);Y position (mm)");
  gC_Y->GetXaxis()->SetTitle("EPICS time (index)");
  gC_Y->GetYaxis()->SetTitle("Y position (mm)");
  gC_Y->SetLineColor(kRed);
  gC_Y->SetMarkerColor(kRed);
  gC_Y->SetMarkerStyle(32);

  // Optionally, write these graphs to the output file later

  // Helper function to draw raw and moving average graphs for a pair of TGraphs
  // Open output file before calling the function
  TFile *fout = TFile::Open(outfile, "RECREATE");

  // Draw scaler rate
  DrawRawAndAvgGraphs(fout, "scalerRate", gTSP, gTSH, "SHMS EREAL Rate", "HMS EREAL Rate", scaler_changes, 4,
                      setting_labels, 5);

  // Draw 3/4 rate
  DrawRawAndAvgGraphs(fout, "34Rate", gTSP_34, gTSH_34, "SHMS 3/4 Rate", "HMS 3/4 Rate", scaler_changes, 4,
                      setting_labels, 5);

  // Draw S1X rate
  DrawRawAndAvgGraphs(fout, "S1XRate", gTSP_S1X, gTSH_S1X, "SHMS S1X Rate", "HMS S1X Rate", scaler_changes, 4,
                      setting_labels, 5);

  // Create a canvas for E branch graphs
  // Draw IPM3HO7A on its own canvas
  // Find min/max for A
  double a_ymin = std::min(*std::min_element(a_x.begin(), a_x.end()), *std::min_element(a_y.begin(), a_y.end()));
  double a_ymax = std::max(*std::max_element(a_x.begin(), a_x.end()), *std::max_element(a_y.begin(), a_y.end()));

  // Draw IPM3HO7A on its own canvas
  TCanvas *cA = new TCanvas("cA", "IPM3HO7A", 800, 600);
  gA_X->GetYaxis()->SetRangeUser(a_ymin, a_ymax);
  gA_X->Draw("AP");
  gA_Y->Draw("P SAME");
  // Draw vertical lines at epics_changes
  for (int i = 0; i < 4; ++i) {
    double x = epics_changes[i];
    TLine *l = new TLine(x, a_ymin, x, a_ymax);
    l->SetLineColor(kBlack);
    l->SetLineStyle(2);
    l->Draw();
  }
  for (int i = 0; i < 5; ++i) {
    double x_left  = (i == 0) ? 0 : epics_changes[i - 1];
    double x_right = (i == 4) ? gA_X->GetXaxis()->GetXmax() : epics_changes[i];
    double x_label = 0.5 * (x_left + x_right);
    double y_label = a_ymin + 0.05 * (a_ymax - a_ymin); // 5% from bottom
    TLatex *latex  = new TLatex(x_label, y_label, setting_labels[i].c_str());
    latex->SetTextAlign(22); // center
    latex->SetTextSize(0.03);
    latex->Draw();
  }
  TLegend *legA = new TLegend(0.15, 0.75, 0.35, 0.88);
  legA->AddEntry(gA_X, "IPM3HO7A_X", "lp");
  legA->AddEntry(gA_Y, "IPM3HO7A_Y", "lp");
  legA->Draw();
  cA->Write();

  // Find min/max for B
  double b_ymin = std::min(*std::min_element(b_x.begin(), b_x.end()), *std::min_element(b_y.begin(), b_y.end()));
  double b_ymax = std::max(*std::max_element(b_x.begin(), b_x.end()), *std::max_element(b_y.begin(), b_y.end()));

  // Draw IPM3HO7B on its own canvas
  TCanvas *cB = new TCanvas("cB", "IPM3HO7B", 800, 600);
  gB_X->GetYaxis()->SetRangeUser(b_ymin, b_ymax);
  gB_X->Draw("AP");
  gB_Y->Draw("P SAME");
  for (int i = 0; i < 4; ++i) {
    double x = epics_changes[i];
    TLine *l = new TLine(x, b_ymin, x, b_ymax);
    l->SetLineColor(kBlack);
    l->SetLineStyle(2);
    l->Draw();
  }
  for (int i = 0; i < 5; ++i) {
    double x_left  = (i == 0) ? 0 : epics_changes[i - 1];
    double x_right = (i == 4) ? gB_X->GetXaxis()->GetXmax() : epics_changes[i];
    double x_label = 0.5 * (x_left + x_right);
    double y_label = b_ymin + 0.05 * (b_ymax - b_ymin); // 5% from bottom
    TLatex *latex  = new TLatex(x_label, y_label, setting_labels[i].c_str());
    latex->SetTextAlign(22); // center
    latex->SetTextSize(0.03);
    latex->Draw();
  }
  TLegend *legB = new TLegend(0.15, 0.75, 0.35, 0.88);
  legB->AddEntry(gB_X, "IPM3HO7B_X", "lp");
  legB->AddEntry(gB_Y, "IPM3HO7B_Y", "lp");
  legB->Draw();
  cB->Write();

  // Find min/max for C
  double c_ymin = std::min(*std::min_element(c_x.begin(), c_x.end()), *std::min_element(c_y.begin(), c_y.end()));
  double c_ymax = std::max(*std::max_element(c_x.begin(), c_x.end()), *std::max_element(c_y.begin(), c_y.end()));

  // Draw IPM3HO7C on its own canvas
  TCanvas *cC = new TCanvas("cC", "IPM3HO7C", 800, 600);
  gC_X->GetYaxis()->SetRangeUser(c_ymin, c_ymax);
  gC_X->Draw("AP");
  gC_Y->Draw("P SAME");
  for (int i = 0; i < 4; ++i) {
    double x = epics_changes[i];
    TLine *l = new TLine(x, c_ymin, x, c_ymax);
    l->SetLineColor(kBlack);
    l->SetLineStyle(2);
    l->Draw();
  }
  for (int i = 0; i < 5; ++i) {
    double x_left  = (i == 0) ? 0 : epics_changes[i - 1];
    double x_right = (i == 4) ? gC_X->GetXaxis()->GetXmax() : epics_changes[i];
    double x_label = 0.5 * (x_left + x_right);
    double y_label = c_ymin + 0.05 * (c_ymax - c_ymin); // 5% from bottom
    TLatex *latex  = new TLatex(x_label, y_label, setting_labels[i].c_str());
    latex->SetTextAlign(22); // center
    latex->SetTextSize(0.03);
    latex->Draw();
  }
  TLegend *legC = new TLegend(0.15, 0.75, 0.35, 0.88);
  legC->AddEntry(gC_X, "IPM3HO7C_X", "lp");
  legC->AddEntry(gC_Y, "IPM3HO7C_Y", "lp");
  legC->Draw();
  cC->Write();

  // Draw IBCM3H04B on its own canvas
  TCanvas *cIBCM = new TCanvas("cIBCM", "IBCM3H04B", 800, 600);
  TGraph *gIBCM  = new TGraph(nE, e_times.data(), ibcm3h04b.data());
  gIBCM->SetTitle("IBCM3H04B;EPICS time;IBCM3H04B (uA)");
  gIBCM->SetLineColor(kGreen);
  gIBCM->SetMarkerColor(kGreen);
  gIBCM->SetMarkerStyle(20);
  gIBCM->GetYaxis()->SetRangeUser(ibcm_axis_min, ibcm_axis_max);
  gIBCM->Draw("AP");
  // Draw vertical lines at epics_changes
  for (int i = 0; i < 4; ++i) {
    double x = epics_changes[i];
    TLine *l = new TLine(x, ibcm_axis_min, x, ibcm_axis_max);
    l->SetLineColor(kBlack);
    l->SetLineStyle(2);
    l->Draw();
  }
  for (int i = 0; i < 5; ++i) {
    double x_left  = (i == 0) ? 0 : epics_changes[i - 1];
    double x_right = (i == 4) ? gIBCM->GetXaxis()->GetXmax() : epics_changes[i];
    double x_label = 0.5 * (x_left + x_right);
    double y_label = ibcm_axis_min + 0.05 * (ibcm_axis_max - ibcm_axis_min); // 5% from bottom
    TLatex *latex  = new TLatex(x_label, y_label, setting_labels[i].c_str());
    latex->SetTextAlign(22); // center
    latex->SetTextSize(0.03);
    latex->Draw();
  }
  cIBCM->Write();
  // Clean up
  fout->Close();
  fin->Close();
}

/*
IIC3H04Pk
IIC3H07Pk
IICSHMSPk
*/