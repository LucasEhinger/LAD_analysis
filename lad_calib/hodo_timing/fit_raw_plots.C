#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>
#include <vector>

using namespace std;

const int N_PLANES  = 5;
const int N_PADDLES = 11;

const string plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};

void fit_raw_plots() {
  gROOT->SetBatch(kTRUE);
  // Path to your ROOT file
  TString filename = "files/raw_hodo_timing_plots/raw_hodo_timing_plots_C3_H.root";
  TFile *file      = TFile::Open(filename, "READ");
  if (!file || file->IsZombie()) {
    cerr << "Error opening file: " << filename << endl;
    return;
  }

  // List of histogram names
  map<string, TString> canvasNames = {
      {"ToF", "KIN/ToF/Anti-Matching_Hit_Tol_0/c_ToF_Anti-Matching_Hit_Tol_0_plane_%s"},
      {"Time_Hodo_Sub", "KIN/ToF_Trigger/Anti-Matching_Hit_Tol_0/c_ToF_Trigger_Anti-Matching_Hit_Tol_0_plane_%s"},
      {"Time_Avg", "KIN/Time_Avg/Anti-Matching_Hit_Tol_0/c_Time_Avg_Anti-Matching_Hit_Tol_0_plane_%s"}
      // Add more short_name to canvasName mappings as needed
  };

  vector<vector<vector<double>>> peak_width(canvasNames.size(),
                                            vector<vector<double>>(N_PLANES, vector<double>(N_PADDLES)));
  vector<vector<vector<double>>> peak_center(canvasNames.size(),
                                             vector<vector<double>>(N_PLANES, vector<double>(N_PADDLES)));
  size_t canvasIdx = 0;
  for (const auto &name : canvasNames) {
    for (int plane = 0; plane < N_PLANES; ++plane) {

      TString fullName = Form(name.second.Data(), plane_names[plane].c_str());
      TCanvas *canvas  = dynamic_cast<TCanvas *>(file->Get(fullName));
      if (!canvas) {
        cerr << "Canvas " << name.second << " not found in file." << endl;
        continue;
      }
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        if (paddle >= canvas->GetListOfPrimitives()->GetSize())
          continue;
        // Get the TPad corresponding to this paddle
        TObject *obj = canvas->GetListOfPrimitives()->At(paddle);
        TPad *pad    = dynamic_cast<TPad *>(obj);
        if (!pad) {
          cerr << "Primitive at paddle " << paddle << " is not a TPad." << endl;
          continue;
        }
        // Assume the first primitive in the pad is the histogram
        TObject *hobj = pad->GetListOfPrimitives()->First();
        TH1 *hist     = dynamic_cast<TH1 *>(hobj);
        if (hist) {
          // Find the maximum bin and its content
          int maxBin  = hist->GetMaximumBin();
          double maxY = hist->GetBinContent(maxBin);

          // Define threshold as a fraction of the peak height (e.g., 50%)
          double sumY = 0;
          int nBins   = hist->GetNbinsX();
          for (int i = 1; i <= nBins; ++i) {
            sumY += hist->GetBinContent(i);
          }
          double averageY  = sumY / nBins;
          double threshold = averageY + (maxY - averageY) * 0.5;

          // Find left crossing
          int leftBin = maxBin;
          while (leftBin > 1 && hist->GetBinContent(leftBin) > threshold)
            --leftBin;
          double xLeft = hist->GetBinCenter(leftBin);

          // Find right crossing
          int rightBin = maxBin;
          while (rightBin < hist->GetNbinsX() && hist->GetBinContent(rightBin) > threshold)
            ++rightBin;
          double xRight = hist->GetBinCenter(rightBin);

          // Calculate width at threshold
          double width = xRight - xLeft;
          // Store the peak width and center
          peak_width[canvasIdx][plane][paddle]  = width;
          peak_center[canvasIdx][plane][paddle] = hist->GetBinCenter(maxBin);
        } else {
          cerr << "No histogram found in pad for paddle " << paddle << "." << endl;
        }
      }
    }
    ++canvasIdx;
  }
  file->Close();

  // Create output ROOT file
  TFile *outfile = TFile::Open("output_summary.root", "RECREATE");

  // Colors for different canvasNames
  vector<int> colors = {kRed, kBlue, kGreen + 2, kMagenta, kCyan + 2, kOrange + 7};

  // Store the keys (short names) of canvasNames in a vector for indexed access
  vector<string> canvasKeys;
  for (const auto &kv : canvasNames) {
    canvasKeys.push_back(kv.first);
  }

  for (int plane = 0; plane < N_PLANES; ++plane) {
    // Peak width plot
    TCanvas *c_width = new TCanvas(Form("c_width_plane_%s", plane_names[plane].c_str()),
                                   Form("Peak Widths Plane %s", plane_names[plane].c_str()), 800, 600);
    c_width->cd();
    vector<TH1D *> hists_width;
    // Find min and max for y-axis range
    double minY_width = std::numeric_limits<double>::max();
    double maxY_width = -std::numeric_limits<double>::max();
    for (size_t i = 0; i < canvasKeys.size(); ++i) {
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        double val = peak_width[i][plane][paddle];
        if (val < minY_width)
          minY_width = val;
        if (val > maxY_width)
          maxY_width = val;
      }
    }
    // Add some margin
    double margin_width = 0.05 * (maxY_width - minY_width);
    minY_width -= margin_width;
    maxY_width += margin_width;

    for (size_t i = 0; i < canvasKeys.size(); ++i) {
      TString hname = Form("width_%s_plane_%s", canvasKeys[i].c_str(), plane_names[plane].c_str());
      TH1D *h       = new TH1D(
          hname, Form("Peak Widths %s Plane %s;Paddle;Width (ns)", canvasKeys[i].c_str(), plane_names[plane].c_str()),
          N_PADDLES, 0 - 0.5, N_PADDLES - 0.5);
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        h->SetBinContent(paddle + 1, peak_width[i][plane][paddle]);
      }
      h->SetLineColor(colors[i % colors.size()]);
      h->SetMarkerColor(colors[i % colors.size()]);
      h->SetMarkerStyle(20 + i);
      h->SetMinimum(minY_width);
      h->SetMaximum(maxY_width);
      hists_width.push_back(h);
      if (i == 0)
        h->Draw("PL");
      else
        h->Draw("PL SAME");
    }
    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (size_t i = 0; i < canvasKeys.size(); ++i) {
      legend->AddEntry(hists_width[i], canvasKeys[i].c_str(), "pl");
    }
    legend->Draw();
    c_width->Write();

    // Peak center plot
    TCanvas *c_center = new TCanvas(Form("c_center_plane_%s", plane_names[plane].c_str()),
                                    Form("Peak Centers Plane %s", plane_names[plane].c_str()), 800, 600);
    c_center->cd();
    vector<TH1D *> hists_center;
    // Find min and max for y-axis range
    double minY = std::numeric_limits<double>::max();
    double maxY = -std::numeric_limits<double>::max();
    for (size_t i = 0; i < canvasKeys.size(); ++i) {
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        double val = peak_center[i][plane][paddle];
        if (val < minY)
          minY = val;
        if (val > maxY)
          maxY = val;
      }
    }
    // Add some margin
    double margin = 0.05 * (maxY - minY);
    minY -= margin;
    maxY += margin;

    for (size_t i = 0; i < canvasKeys.size(); ++i) {
      TString hname = Form("center_%s_plane_%s", canvasKeys[i].c_str(), plane_names[plane].c_str());
      TH1D *h       = new TH1D(
          hname, Form("Peak Centers %s Plane %s;Paddle;Center", canvasKeys[i].c_str(), plane_names[plane].c_str()),
          N_PADDLES, 0 - 0.5, N_PADDLES - 0.5);
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        h->SetBinContent(paddle + 1, peak_center[i][plane][paddle]);
      }
      h->SetLineColor(colors[i % colors.size()]);
      h->SetMarkerColor(colors[i % colors.size()]);
      h->SetMarkerStyle(20 + i);
      h->SetMinimum(minY);
      h->SetMaximum(maxY);
      hists_center.push_back(h);
      if (i == 0)
        h->Draw("PL");
      else
        h->Draw("PL SAME");
    }
    auto legend_center = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (size_t i = 0; i < canvasKeys.size(); ++i) {
      legend_center->AddEntry(hists_center[i], canvasKeys[i].c_str(), "pl");
    }
    legend_center->Draw();
    c_center->Write();
  }

  outfile->Close();
}