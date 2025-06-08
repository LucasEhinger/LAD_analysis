#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TTree.h>
#include <iostream>
#include <string>

using namespace std;

const int MAX_DATA                = 1000; // Maximum number of data points
const int nplanes                 = 5;    // Number of planes (SHMS and HMS)
const string plane_names[nplanes] = {"000", "001", "100", "101", "200"};
const string spec_names[2]        = {"P", "H"};
const double HMS_SCALE            = 3.1;
const double TDC2NS               = 0.09766; // ns per TDC channel
const double ADC_OFFSET[2]        = {-1850, -1797};

int compare_SHMS_HMS() {
  gROOT->SetBatch(kTRUE);
  // Open the ROOT file
  TFile *file = TFile::Open(
      "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/bad_timing/LAD_COIN_22844_0_4_-1.root");
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file!" << std::endl;
    return 1;
  }

  // Get the tree named "T"
  TTree *tree = (TTree *)file->Get("T");
  if (!tree) {
    std::cerr << "Tree 'T' not found!" << std::endl;
    file->Close();
    return 1;
  }

  // Get the number of branches
  Double_t hodo_hit_time[2][nplanes][MAX_DATA];
  Int_t ndata_hodo_hit_time[2][nplanes];
  Double_t TopTdcTimeRaw[2][nplanes][MAX_DATA];
  Int_t ndata_TopTdcTimeRaw[2][nplanes];
  Double_t BtmTdcTimeRaw[2][nplanes][MAX_DATA];
  Int_t ndata_BtmTdcTimeRaw[2][nplanes];
  Double_t TopAdcTime[2][nplanes][MAX_DATA];
  Int_t ndata_TopAdcTimeRaw[2][nplanes];
  Double_t BtmAdcTimeRaw[2][nplanes][MAX_DATA];
  Int_t ndata_BtmAdcTimeRaw[2][nplanes];
  Double_t TdcRefTime[2];

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < nplanes; ++j) {
      tree->SetBranchAddress(Form("%s.ladhod.%s.HodoHitTime", spec_names[i].c_str(), plane_names[j].c_str()),
                             hodo_hit_time[i][j]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.TopTdcTimeRaw", spec_names[i].c_str(), plane_names[j].c_str()),
                             TopTdcTimeRaw[i][j]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.BtmTdcTimeRaw", spec_names[i].c_str(), plane_names[j].c_str()),
                             BtmTdcTimeRaw[i][j]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.TopAdcPulseTime", spec_names[i].c_str(), plane_names[j].c_str()),
                             TopAdcTime[i][j]);
      tree->SetBranchAddress(Form("%s.ladhod.%s.BtmAdcPulseTime", spec_names[i].c_str(), plane_names[j].c_str()),
                             BtmAdcTimeRaw[i][j]);

      tree->SetBranchAddress(Form("Ndata.%s.ladhod.%s.HodoHitTime", spec_names[i].c_str(), plane_names[j].c_str()),
                             &ndata_hodo_hit_time[i][j]);
      tree->SetBranchAddress(Form("Ndata.%s.ladhod.%s.TopTdcTimeRaw", spec_names[i].c_str(), plane_names[j].c_str()),
                             &ndata_TopTdcTimeRaw[i][j]);
      tree->SetBranchAddress(Form("Ndata.%s.ladhod.%s.BtmTdcTimeRaw", spec_names[i].c_str(), plane_names[j].c_str()),
                             &ndata_BtmTdcTimeRaw[i][j]);
      tree->SetBranchAddress(Form("Ndata.%s.ladhod.%s.TopAdcPulseTime", spec_names[i].c_str(), plane_names[j].c_str()),
                             &ndata_TopAdcTimeRaw[i][j]);
      tree->SetBranchAddress(Form("Ndata.%s.ladhod.%s.BtmAdcPulseTime", spec_names[i].c_str(), plane_names[j].c_str()),
                             &ndata_BtmAdcTimeRaw[i][j]);
    }
  }

  tree->SetBranchAddress("T.shms.shmsTrigLAD_tdcTimeRaw", &TdcRefTime[0]);
  tree->SetBranchAddress("T.hms.hmsTrigLAD_tdcTimeRaw", &TdcRefTime[1]);

  TH1F *hodo_hit_time_hist[nplanes][2];
  TH1F *TopTdcTimeRaw_hist[nplanes][2];
  TH1F *BtmTdcTimeRaw_hist[nplanes][2];
  TH1F *TopAdcTime_hist[nplanes][2];
  TH1F *BtmAdcTime_hist[nplanes][2];
  TH1F *RefTimeRaw_hist[2];

  for (int i = 0; i < nplanes; ++i) {
    for (int j = 0; j < 2; ++j) {
      hodo_hit_time_hist[i][j] = new TH1F(
          Form("hodo_hit_time_%s_%s", plane_names[i].c_str(), spec_names[j].c_str()),
          Form("Hodo Hit Time %s %s;Time (ns)", plane_names[i].c_str(), spec_names[j].c_str()), 100, 1500, 2100);
      TopTdcTimeRaw_hist[i][j] = new TH1F(
          Form("TopTdcTimeRaw_%s_%s", plane_names[i].c_str(), spec_names[j].c_str()),
          Form("Top TDC Time Raw %s %s;Time (ns)", plane_names[i].c_str(), spec_names[j].c_str()), 100, 0, 400);
      BtmTdcTimeRaw_hist[i][j] = new TH1F(
          Form("BtmTdcTimeRaw_%s_%s", plane_names[i].c_str(), spec_names[j].c_str()),
          Form("Bottom TDC Time Raw %s %s;Time (ns)", plane_names[i].c_str(), spec_names[j].c_str()), 100, 0, 400);
      TopAdcTime_hist[i][j] =
          new TH1F(Form("TopAdcTime_%s_%s", plane_names[i].c_str(), spec_names[j].c_str()),
                   Form("Top ADC Time %s %s;Time (ns)", plane_names[i].c_str(), spec_names[j].c_str()), 100, -300, 300);
      BtmAdcTime_hist[i][j] = new TH1F(
          Form("BtmAdcTime_%s_%s", plane_names[i].c_str(), spec_names[j].c_str()),
          Form("Bottom ADC Time %s %s;Time (ns)", plane_names[i].c_str(), spec_names[j].c_str()), 100, -300, 300);
    }
  }
  RefTimeRaw_hist[0] = new TH1F("RefTimeRaw_SHMS", "SHMS Reference Time Raw;Time (ns)", 100, 1480, 1600);
  RefTimeRaw_hist[1] = new TH1F("RefTimeRaw_HMS", "HMS Reference Time Raw;Time (ns)", 100, 1480, 1600);

  // Loop over the entries in the tree
  Long64_t nentries = tree->GetEntries();
  // nentries          = 5000;
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);

    for (int j = 0; j < nplanes; ++j) {
      for (int k = 0; k < 2; ++k) {
        for (int l = 0; l < ndata_hodo_hit_time[k][j]; ++l) {
          hodo_hit_time_hist[j][k]->Fill(hodo_hit_time[k][j][l]);
        }
        for (int l = 0; l < ndata_TopTdcTimeRaw[k][j]; ++l) {
          TopTdcTimeRaw_hist[j][k]->Fill(TopTdcTimeRaw[k][j][l] * TDC2NS);
        }
        for (int l = 0; l < ndata_BtmTdcTimeRaw[k][j]; ++l) {
          BtmTdcTimeRaw_hist[j][k]->Fill(BtmTdcTimeRaw[k][j][l] * TDC2NS);
        }
        for (int l = 0; l < ndata_TopAdcTimeRaw[k][j]; ++l) {
          TopAdcTime_hist[j][k]->Fill(TopAdcTime[k][j][l] + ADC_OFFSET[k]);
        }
        for (int l = 0; l < ndata_BtmAdcTimeRaw[k][j]; ++l) {
          BtmAdcTime_hist[j][k]->Fill(BtmAdcTimeRaw[k][j][l] + ADC_OFFSET[k]);
        }
      }
    }
    // Fill the reference time histogram

    RefTimeRaw_hist[0]->Fill(TdcRefTime[0] * TDC2NS);
    RefTimeRaw_hist[1]->Fill(TdcRefTime[1] * TDC2NS);

    if (i % 10000 == 0 || i == nentries - 1) {
      double progress = 100.0 * (i + 1) / nentries;
      printf("\rProgress: %.2f%%", progress);
      fflush(stdout);
    }
  }

  TFile *outfile = new TFile("HodoTimingComparison.root", "RECREATE");
  // Save histograms to the file
  outfile->cd();

  // Create a canvas and draw SHMS (P) and HMS (H) histograms for each plane
  TCanvas *c = new TCanvas("c", "SHMS vs HMS Hodo Timing", 1200, 800);

  for (int i = 0; i < nplanes; ++i) {
    c->Clear();
    c->Divide(2, 2);

    // Draw and save for each variable
    const char *varnames[4] = {"TopTdcTimeRaw", "BtmTdcTimeRaw", "TopAdcTime", "BtmAdcTime"};
    TH1F *(*hists[4])[2]    = {TopTdcTimeRaw_hist, BtmTdcTimeRaw_hist, TopAdcTime_hist, BtmAdcTime_hist};

    for (int v = 0; v < 4; ++v) {
      c->cd(v + 1);
      // Clone and scale SHMS (P) histogram
      TH1F *shms = (TH1F *)((hists[v])[i][0]->Clone());
      shms->Scale(HMS_SCALE);
      shms->SetLineColor(kRed);
      shms->SetTitle(Form("%s %s: SHMS (red, scaled) vs HMS (black)", varnames[v], plane_names[i].c_str()));

      // Find max y for both histograms
      double max_shms = shms->GetMaximum();
      double max_hms  = hists[v][i][1]->GetMaximum();
      double max_y    = std::max(max_shms, max_hms) * 1.1; // add 10% headroom

      shms->SetMaximum(max_y);
      hists[v][i][1]->SetMaximum(max_y);

      shms->Draw("hist");
      // Draw HMS (H) histogram
      hists[v][i][1]->SetLineColor(kBlack);
      hists[v][i][1]->Draw("hist same");
      shms->Draw("hist same");
      // Add legend
      TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
      leg->AddEntry(shms, "SHMS (scaled)", "l");
      leg->AddEntry(hists[v][i][1], "HMS", "l");
      leg->Draw();
    }

    // c->SaveAs(Form("HodoTimingComparison/Plane_%s_Comparison.png", plane_names[i].c_str()));
    c->Write(Form("Plane_%s_Comparison", plane_names[i].c_str()));
  }

  // Draw and save reference time histograms
  c->Clear();
  // Scale SHMS (RefTimeRaw_hist[0]) by HMS_SCALE and draw both on the same axis
  TH1F *shms_ref_scaled = (TH1F *)RefTimeRaw_hist[0]->Clone("shms_ref_scaled");
  shms_ref_scaled->Scale(HMS_SCALE);
  shms_ref_scaled->SetLineColor(kRed);
  shms_ref_scaled->SetTitle("Reference Time Raw: SHMS (red, scaled) vs HMS (black)");
  shms_ref_scaled->Draw("hist");

  RefTimeRaw_hist[1]->SetLineColor(kBlack);
  RefTimeRaw_hist[1]->Draw("hist same");

  // Add legend
  TLegend *leg_ref = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg_ref->AddEntry(shms_ref_scaled, "SHMS (scaled)", "l");
  leg_ref->AddEntry(RefTimeRaw_hist[1], "HMS", "l");
  leg_ref->Draw();

  c->Write("ReferenceTimeRaw_Comparison");

  file->Close();
  return 0;
}