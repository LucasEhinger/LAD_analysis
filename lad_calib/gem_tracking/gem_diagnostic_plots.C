#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TTree.h>
#include <iostream>
#include <vector>

const int MAX_DATA        = 10000;
const int d0_NBINS        = 100;
const double d0_MIN       = 0.0;
const double d0_MAX       = 100.0;
const int projz_NBINS     = 50;
const double projz_MIN    = -20.0;
const double projz_MAX    = 20.0;
const int projy_NBINS     = 100;
const double projy_MIN    = -20.0;
const double projy_MAX    = 20.0;
const int deltapos_NBINS  = 100;
const double deltapos_MIN = -500.0;
const double deltapos_MAX = 500.0;
const int time_NBINS      = 100;
const double time_MIN     = 1700.0;
const double time_MAX     = 2000.0;
const double theta_MIN    = 0.0;
const double theta_MAX    = 180.0;
const int theta_NBINS     = 180;

const double phi_MIN = -180.0;
const double phi_MAX = 180.0;
const int phi_NBINS  = 360;

const double paddle_MIN = 0.0;
const double paddle_MAX = 100.0;
const int paddle_NBINS  = 100;

const double hit_position_MIN = -500.0;
const double hit_position_MAX = 500.0;
const int hit_position_NBINS  = 100;

const double edep_MIN = 0.0;
const double edep_MAX = 10.0;
const int edep_NBINS  = 100;

const int nTransCuts               = 5;
const double transCuts[nTransCuts] = {500, 300, 200, 100, 50};
const int nLongCuts                = 5;
const double longCuts[nLongCuts]   = {200, 150, 100, 50};
const int nD0Cuts                  = 5;
const double d0Cuts[nD0Cuts]       = {100, 50, 30, 20, 10};

const double delta_pos_trans_sig = 80.0;
const double delta_pos_long_sig  = 80.0;

const int maxTracks = 30;

void gem_diagnostic_plots(int runnum = 0) {
  // Set batch mode to suppress graphical output
  gROOT->SetBatch(kTRUE);
  // Open the ROOT file
  TString fileName = Form("/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/"
                          "LAD_COIN_%d_0_0_-1.root",
                          runnum);
  //  "LAD_COIN_22073_-1.root";
  //  "LAD_COIN_22282_0_0_500002.root";
  TString outputFileName = Form("files/gem_diagnostic_plots/gem_diagnostic_plots_%d_-1.root", runnum);
  // Open the ROOT file
  TFile *file = TFile::Open(fileName);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open the ROOT file!" << std::endl;
    return;
  }

  // Get the TTree
  TTree *T = dynamic_cast<TTree *>(file->Get("T"));
  if (!T) {
    std::cerr << "Error: Cannot find the TTree named 'T'!" << std::endl;
    file->Close();
    return;
  }

  // Define arrays to hold the data
  Double_t trk_d0[MAX_DATA];
  Double_t trk_d0_good[MAX_DATA];
  Double_t trk_projz[MAX_DATA], trk_projy[MAX_DATA];
  Double_t trk_time[MAX_DATA], trk_dt[MAX_DATA];
  Double_t kin_trackID_0[MAX_DATA], kin_trackID_1[MAX_DATA];
  Double_t kin_plane_0[MAX_DATA], kin_plane_1[MAX_DATA];
  Double_t kin_paddle_0[MAX_DATA], kin_paddle_1[MAX_DATA];
  Double_t kin_hittime_0[MAX_DATA], kin_hittime_1[MAX_DATA];
  Double_t kin_hittheta_0[MAX_DATA], kin_hittheta_1[MAX_DATA];
  Double_t kin_hitphi_0[MAX_DATA], kin_hitphi_1[MAX_DATA];
  Double_t kin_hitedep_0[MAX_DATA], kin_hitedep_1[MAX_DATA];
  Double_t kin_deltapostrans_0[MAX_DATA], kin_deltapostrans_1[MAX_DATA];
  Double_t kin_deltaposlong_0[MAX_DATA], kin_deltaposlong_1[MAX_DATA];
  Int_t nTracks, nGoodHits;
  Double_t reactz;

  char spec_prefix = 'H'; // Define the character variable for the prefix

  T->SetBranchAddress(Form("Ndata.%c.gem.trk.d0", spec_prefix), &nTracks);
  T->SetBranchAddress(Form("Ndata.%c.ladkin.goodhit_trackid_0", spec_prefix), &nGoodHits);
  T->SetBranchAddress(Form("%c.gem.trk.d0", spec_prefix), &trk_d0);
  T->SetBranchAddress(Form("%c.gem.trk.d0_good", spec_prefix), &trk_d0_good);
  T->SetBranchAddress(Form("%c.gem.trk.projz", spec_prefix), &trk_projz);
  T->SetBranchAddress(Form("%c.gem.trk.projy", spec_prefix), &trk_projy);
  T->SetBranchAddress(Form("%c.gem.trk.t", spec_prefix), &trk_time);
  T->SetBranchAddress(Form("%c.gem.trk.dt", spec_prefix), &trk_dt);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_trackid_0", spec_prefix), &kin_trackID_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_trackid_1", spec_prefix), &kin_trackID_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_plane_0", spec_prefix), &kin_plane_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_plane_1", spec_prefix), &kin_plane_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_paddle_0", spec_prefix), &kin_paddle_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_paddle_1", spec_prefix), &kin_paddle_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hittime_0", spec_prefix), &kin_hittime_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hittime_1", spec_prefix), &kin_hittime_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hittheta_0", spec_prefix), &kin_hittheta_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hittheta_1", spec_prefix), &kin_hittheta_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hitphi_0", spec_prefix), &kin_hitphi_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hitphi_1", spec_prefix), &kin_hitphi_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hitedep_0", spec_prefix), &kin_hitedep_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hitedep_1", spec_prefix), &kin_hitedep_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_dTrkHoriz_0", spec_prefix), &kin_deltapostrans_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_dTrkHoriz_1", spec_prefix), &kin_deltapostrans_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_dTrkVert_0", spec_prefix), &kin_deltaposlong_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_dTrkVert_1", spec_prefix), &kin_deltaposlong_1);
  T->SetBranchAddress(Form("%c.react.z", spec_prefix), &reactz);

  // Create a histogram to store the data
  TH1F *h_d0_all =
      new TH1F("h_d0", "d0 Distribution;d0: #Delta track proj to target (cm);Counts", d0_NBINS, d0_MIN, d0_MAX);
  TH1F *h_d0_good = new TH1F("h_d0_good", "d0 Good Distribution;d0: #Delta track proj to target (cm);Counts", d0_NBINS,
                             d0_MIN, d0_MAX);
  TH1F *h_d0_bad =
      new TH1F("h_d0_bad", "d0 Bad Distribution;d0: #Delta track proj to target (cm);Counts", d0_NBINS, d0_MIN, d0_MAX);

  // Create histograms for projz
  TH1F *h_projz_all  = new TH1F("h_projz", "projz Distribution;projz: track proj to z-axis (cm);Counts", projz_NBINS,
                                projz_MIN, projz_MAX);
  TH1F *h_projz_good = new TH1F("h_projz_good", "projz Good Distribution;projz: track proj to z-axis (cm);Counts",
                                projz_NBINS, projz_MIN, projz_MAX);
  TH1F *h_projz_bad  = new TH1F("h_projz_bad", "projz Bad Distribution;projz: track proj to z-axis (cm);Counts",
                                projz_NBINS, projz_MIN, projz_MAX);
  // Create histograms for projy
  TH1F *h_projy_all  = new TH1F("h_projy", "projy Distribution;projy: track proj to y-axis (cm);Counts", projy_NBINS,
                                projy_MIN, projy_MAX);
  TH1F *h_projy_good = new TH1F("h_projy_good", "projy Good Distribution;projy: track proj to y-axis (cm);Counts",
                                projy_NBINS, projy_MIN, projy_MAX);
  TH1F *h_projy_bad  = new TH1F("h_projy_bad", "projy Bad Distribution;projy: track proj to y-axis (cm);Counts",
                                projy_NBINS, projy_MIN, projy_MAX);

  TH1F *h_delta_pos_long =
      new TH1F("h_delta_pos_long", "Delta Pos Long Distribution;#Delta track projection to LAD long (cm);Counts",
               deltapos_NBINS, deltapos_MIN, deltapos_MAX);
  TH1F *h_delta_pos_trans =
      new TH1F("h_delta_pos_trans", "Delta Pos Trans Distribution;#Delta track projection to LAD trans (cm);Counts",
               deltapos_NBINS, deltapos_MIN, deltapos_MAX);

  TH1F *h_hittime = new TH1F("h_hittime", "Hit Time Distribution;Hit Time (ns);Counts", time_NBINS, time_MIN, time_MAX);

  TH2F *h_d0_vs_deltapostrans = new TH2F(
      "h_d0_vs_deltapostrans",
      "Track d0 vs Delta Pos Trans;d0: #Delta track proj to target (cm);#Delta track projection to LAD trans (cm)",
      d0_NBINS, d0_MIN, d0_MAX, deltapos_NBINS, deltapos_MIN, deltapos_MAX);

  TH2F *h_projz_vs_deltapostrans =
      new TH2F("h_projz_vs_deltapostrans",
               "Projz vs Delta Pos Trans;projz: #Delta track proj to target (cm);#Delta track projection to LAD (cm)",
               projz_NBINS, projz_MIN, projz_MAX, deltapos_NBINS, deltapos_MIN, deltapos_MAX);

  TH2F *h_projy_vs_deltapostrans = new TH2F(
      "h_projy_vs_deltapostrans",
      "Track projy vs Delta Pos Trans;projy: track proj to y-axis (cm);#Delta track projection to LAD trans (cm)",
      projy_NBINS, projy_MIN, projy_MAX, deltapos_NBINS, deltapos_MIN, deltapos_MAX);
  // Create 2D histograms for d0 vs delta pos long and delta pos long vs delta pos trans
  TH2F *h_d0_vs_deltaposlong = new TH2F(
      "h_d0_vs_deltaposlong",
      "Track d0 vs Delta Pos Long;d0: #Delta track proj to target (cm);#Delta track projection to LAD long (cm)",
      d0_NBINS, d0_MIN, d0_MAX, deltapos_NBINS, deltapos_MIN, deltapos_MAX);

  TH2F *h_deltaposlong_vs_deltapostrans =
      new TH2F("h_deltaposlong_vs_deltapostrans",
               "Delta Pos Long vs Delta Pos Trans;#Delta track projection to LAD trans (cm);#Delta track projection to "
               "LAD long (cm)",
               deltapos_NBINS, deltapos_MIN, deltapos_MAX, deltapos_NBINS, deltapos_MIN, deltapos_MAX);

  // Create histograms for each trans cut
  std::vector<TH1F *> h_d0_all_trans(nTransCuts);
  std::vector<TH1F *> h_d0_good_trans(nTransCuts);
  std::vector<TH1F *> h_d0_bad_trans(nTransCuts);
  std::vector<TH1F *> h_projz_all_trans(nTransCuts);
  std::vector<TH1F *> h_projz_good_trans(nTransCuts);
  std::vector<TH1F *> h_projz_bad_trans(nTransCuts);
  std::vector<TH1F *> h_projy_all_trans(nTransCuts);
  std::vector<TH1F *> h_projy_good_trans(nTransCuts);
  std::vector<TH1F *> h_projy_bad_trans(nTransCuts);

  std::vector<TH1F *> h_d0_all_long(nLongCuts);
  std::vector<TH1F *> h_d0_good_long(nLongCuts);
  std::vector<TH1F *> h_d0_bad_long(nLongCuts);
  std::vector<TH1F *> h_projz_all_long(nLongCuts);
  std::vector<TH1F *> h_projz_good_long(nLongCuts);
  std::vector<TH1F *> h_projz_bad_long(nLongCuts);
  std::vector<TH1F *> h_projy_all_long(nLongCuts);
  std::vector<TH1F *> h_projy_good_long(nLongCuts);
  std::vector<TH1F *> h_projy_bad_long(nLongCuts);

  std::vector<TH1F *> h_delta_pos_long_d0cut(nD0Cuts);
  std::vector<TH1F *> h_delta_pos_trans_d0cut(nD0Cuts);

  for (int i = 0; i < nTransCuts; ++i) {
    h_d0_all_trans[i] =
        new TH1F(Form("h_d0_all_trans_%d", i),
                 Form("d0 All Distribution (Trans Cut %.1f);d0: #Delta track proj to target (cm);Counts", transCuts[i]),
                 d0_NBINS, d0_MIN, d0_MAX);
    h_d0_good_trans[i] = new TH1F(
        Form("h_d0_good_trans_%d", i),
        Form("d0 Good Distribution (Trans Cut %.1f);d0: #Delta track proj to target (cm);Counts", transCuts[i]),
        d0_NBINS, d0_MIN, d0_MAX);
    h_d0_bad_trans[i] =
        new TH1F(Form("h_d0_bad_trans_%d", i),
                 Form("d0 Bad Distribution (Trans Cut %.1f);d0: #Delta track proj to target (cm);Counts", transCuts[i]),
                 d0_NBINS, d0_MIN, d0_MAX);

    h_projz_all_trans[i] =
        new TH1F(Form("h_projz_all_trans_%d", i),
                 Form("projz All Distribution (Trans Cut %.1f);projz: track proj to z-axis (cm);Counts", transCuts[i]),
                 projz_NBINS, projz_MIN, projz_MAX);
    h_projz_good_trans[i] =
        new TH1F(Form("h_projz_good_trans_%d", i),
                 Form("projz Good Distribution (Trans Cut %.1f);projz: track proj to z-axis (cm);Counts", transCuts[i]),
                 projz_NBINS, projz_MIN, projz_MAX);
    h_projz_bad_trans[i] =
        new TH1F(Form("h_projz_bad_trans_%d", i),
                 Form("projz Bad Distribution (Trans Cut %.1f);projz: track proj to z-axis (cm);Counts", transCuts[i]),
                 projz_NBINS, projz_MIN, projz_MAX);

    h_projy_all_trans[i] =
        new TH1F(Form("h_projy_all_trans_%d", i),
                 Form("projy All Distribution (Trans Cut %.1f);projy: track proj to y-axis (cm);Counts", transCuts[i]),
                 projy_NBINS, projy_MIN, projy_MAX);
    h_projy_good_trans[i] =
        new TH1F(Form("h_projy_good_trans_%d", i),
                 Form("projy Good Distribution (Trans Cut %.1f);projy: track proj to y-axis (cm);Counts", transCuts[i]),
                 projy_NBINS, projy_MIN, projy_MAX);
    h_projy_bad_trans[i] =
        new TH1F(Form("h_projy_bad_trans_%d", i),
                 Form("projy Bad Distribution (Trans Cut %.1f);projy: track proj to y-axis (cm);Counts", transCuts[i]),
                 projy_NBINS, projy_MIN, projy_MAX);
  }

  for (int i = 0; i < nLongCuts; ++i) {
    h_d0_all_long[i] =
        new TH1F(Form("h_d0_all_long_%d", i),
                 Form("d0 All Distribution (Long Cut %.1f);d0: #Delta track proj to target (cm);Counts", longCuts[i]),
                 d0_NBINS, d0_MIN, d0_MAX);
    h_d0_good_long[i] =
        new TH1F(Form("h_d0_good_long_%d", i),
                 Form("d0 Good Distribution (Long Cut %.1f);d0: #Delta track proj to target (cm);Counts", longCuts[i]),
                 d0_NBINS, d0_MIN, d0_MAX);
    h_d0_bad_long[i] =
        new TH1F(Form("h_d0_bad_long_%d", i),
                 Form("d0 Bad Distribution (Long Cut %.1f);d0: #Delta track proj to target (cm);Counts", longCuts[i]),
                 d0_NBINS, d0_MIN, d0_MAX);

    h_projz_all_long[i] =
        new TH1F(Form("h_projz_all_long_%d", i),
                 Form("projz All Distribution (Long Cut %.1f);projz: track proj to z-axis (cm);Counts", longCuts[i]),
                 projz_NBINS, projz_MIN, projz_MAX);
    h_projz_good_long[i] =
        new TH1F(Form("h_projz_good_long_%d", i),
                 Form("projz Good Distribution (Long Cut %.1f);projz: track proj to z-axis (cm);Counts", longCuts[i]),
                 projz_NBINS, projz_MIN, projz_MAX);
    h_projz_bad_long[i] =
        new TH1F(Form("h_projz_bad_long_%d", i),
                 Form("projz Bad Distribution (Long Cut %.1f);projz: track proj to z-axis (cm);Counts", longCuts[i]),
                 projz_NBINS, projz_MIN, projz_MAX);
    h_projy_all_long[i] =
        new TH1F(Form("h_projy_all_long_%d", i),
                 Form("projy All Distribution (Long Cut %.1f);projy: track proj to y-axis (cm);Counts", longCuts[i]),
                 projy_NBINS, projy_MIN, projy_MAX);
    h_projy_good_long[i] =
        new TH1F(Form("h_projy_good_long_%d", i),
                 Form("projy Good Distribution (Long Cut %.1f);projy: track proj to y-axis (cm);Counts", longCuts[i]),
                 projy_NBINS, projy_MIN, projy_MAX);
    h_projy_bad_long[i] =
        new TH1F(Form("h_projy_bad_long_%d", i),
                 Form("projy Bad Distribution (Long Cut %.1f);projy: track proj to y-axis (cm);Counts", longCuts[i]),
                 projy_NBINS, projy_MIN, projy_MAX);
  }

  for (int i = 0; i < nD0Cuts; ++i) {
    h_delta_pos_long_d0cut[i] = new TH1F(
        Form("h_delta_pos_long_d0cut_%d", i),
        Form("Delta Pos Long Distribution (d0 Cut %.1f);#Delta track projection to LAD long (cm);Counts", d0Cuts[i]),
        deltapos_NBINS, deltapos_MIN, deltapos_MAX);

    h_delta_pos_trans_d0cut[i] = new TH1F(
        Form("h_delta_pos_trans_d0cut_%d", i),
        Form("Delta Pos Trans Distribution (d0 Cut %.1f);#Delta track projection to LAD trans (cm);Counts", d0Cuts[i]),
        deltapos_NBINS, deltapos_MIN, deltapos_MAX);
  }

  // Create histograms for box cut
  TH1F *h_d0_all_box_cut  = new TH1F("h_d0_all_box_cut", "d0 All Box Cut;d0: #Delta track proj to target (cm);Counts",
                                     d0_NBINS, d0_MIN, d0_MAX);
  TH1F *h_d0_good_box_cut = new TH1F("h_d0_good_box_cut", "d0 Good Box Cut;d0: #Delta track proj to target (cm);Counts",
                                     d0_NBINS, d0_MIN, d0_MAX);
  TH1F *h_d0_bad_box_cut  = new TH1F("h_d0_bad_box_cut", "d0 Bad Box Cut;d0: #Delta track proj to target (cm);Counts",
                                     d0_NBINS, d0_MIN, d0_MAX);
  TH1F *h_projz_all_box_cut =
      new TH1F("h_projz_all_box_cut", "projz All Box Cut;projz: track proj to z-axis (cm);Counts", projz_NBINS,
               projz_MIN, projz_MAX);
  TH1F *h_projz_good_box_cut =
      new TH1F("h_projz_good_box_cut", "projz Good Box Cut;projz: track proj to z-axis (cm);Counts", projz_NBINS,
               projz_MIN, projz_MAX);
  TH1F *h_projz_bad_box_cut =
      new TH1F("h_projz_bad_box_cut", "projz Bad Box Cut;projz: track proj to z-axis (cm);Counts", projz_NBINS,
               projz_MIN, projz_MAX);
  TH1F *h_projy_all_box_cut =
      new TH1F("h_projy_all_box_cut", "projy All Box Cut;projy: track proj to y-axis (cm);Counts", projy_NBINS,
               projy_MIN, projy_MAX);
  TH1F *h_projy_good_box_cut =
      new TH1F("h_projy_good_box_cut", "projy Good Box Cut;projy: track proj to y-axis (cm);Counts", projy_NBINS,
               projy_MIN, projy_MAX);
  TH1F *h_projy_bad_box_cut =
      new TH1F("h_projy_bad_box_cut", "projy Bad Box Cut;projy: track proj to y-axis (cm);Counts", projy_NBINS,
               projy_MIN, projy_MAX);

  TH1F *h_hittime_box_cut =
      new TH1F("h_hittime_box_cut", "Hit Time Box Cut;Hit Time (ns);Counts", time_NBINS, time_MIN, time_MAX);

  TH2D *h_projz_vs_reactz_sig =
      new TH2D("h_projz_vs_reactz_sig",
               "Projz vs reactz;projz: GEM track proj to z-axis (cm) sig;Electron z: track proj to z-axis (cm)",
               projz_NBINS, projz_MIN, projz_MAX, projz_NBINS, projz_MIN, projz_MAX);
  TH2D *h_projz_vs_reactz_bdk =
      new TH2D("h_projz_vs_reactz_bdk",
               "Projz vs reactz;projz: GEM track proj to z-axis (cm) bdk;Electron z: track proj to z-axis (cm)",
               projz_NBINS, projz_MIN, projz_MAX, projz_NBINS, projz_MIN, projz_MAX);
  // Create histograms for box cut track dt
  TH1F *h_gem_dt_all_box_cut =
      new TH1F("h_dt_all_box_cut", "Track dt All Box Cut;Track dt (ns);Counts", time_NBINS, -50, 50);
  TH1F *h_gem_dt_good_box_cut =
      new TH1F("h_dt_good_box_cut", "Track dt Good Box Cut;Track dt (ns);Counts", time_NBINS, -50, 50);
  TH1F *h_gem_dt_bad_box_cut =
      new TH1F("h_dt_bad_box_cut", "Track dt Bad Box Cut;Track dt (ns);Counts", time_NBINS, -50, 50);

  // Create time histograms for the box cut
  TH1F *h_gem_time_all_box_cut =
      new TH1F("h_time_all_box_cut", "Time All Box Cut;Time (ns);Counts", time_NBINS, 0, 150);
  TH1F *h_gem_time_good_box_cut =
      new TH1F("h_time_good_box_cut", "Time Good Box Cut;Time (ns);Counts", time_NBINS, 0, 150);
  TH1F *h_gem_time_bad_box_cut =
      new TH1F("h_time_bad_box_cut", "Time Bad Box Cut;Time (ns);Counts", time_NBINS, 0, 150);

  // Create histograms for final cuts
  TH1F *h_zprojection = new TH1F("h_zprojection", "Z Projection Distribution;Z Projection (cm);Counts", projz_NBINS,
                                 projz_MIN, projz_MAX);
  TH1F *h_yprojection = new TH1F("h_yprojection", "Y Projection Distribution;Y Projection (cm);Counts", projy_NBINS,
                                 projy_MIN, projy_MAX);
  TH1F *h_theta       = new TH1F("h_theta", "Theta Distribution;Theta (degrees);Counts", 180, 0, 180);
  TH1F *h_phi         = new TH1F("h_phi", "Phi Distribution;Phi (degrees);Counts", 360, -180, 180);
  TH1F *h_paddle      = new TH1F("h_paddle", "Paddle Distribution;Paddle ID;Counts", 100, 0, 100);
  TH1F *h_hit_position =
      new TH1F("h_hit_position", "Hit Position Distribution;Hit Position (cm);Counts", 100, -500, 500);
  TH1F *h_edep = new TH1F("h_edep", "Energy Deposition Distribution;Energy Deposition (MeV);Counts", 100, 0, 10);
  // Loop over entries and fill arrays
  Long64_t nEntries = T->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    T->GetEntry(i);
    if (nTracks > maxTracks) {
      // std::cerr << "Warning: Number of tracks exceeds maxTracks. Skipping entry " << i << std::endl;
      continue;
    }
    for (int j = 0; j < nTracks; ++j) {
      h_d0_all->Fill(trk_d0[j]);
      h_projz_all->Fill(trk_projz[j]);
      h_projy_all->Fill(trk_projy[j]);

      if (trk_d0_good[j] > 0) {
        h_d0_good->Fill(trk_d0[j]);
        h_projz_good->Fill(trk_projz[j]);
        h_projy_good->Fill(trk_projy[j]);
      } else {
        h_d0_bad->Fill(trk_d0[j]);
        h_projz_bad->Fill(trk_projz[j]);
        h_projy_bad->Fill(trk_projy[j]);
      }
    }
    for (int k = 0; k < nGoodHits; ++k) {
      h_delta_pos_trans->Fill(kin_deltapostrans_0[k]);
      h_delta_pos_long->Fill(kin_deltaposlong_0[k]);
      h_hittime->Fill(kin_hittime_0[k]);

      for (int cutIdx = 0; cutIdx < nD0Cuts; ++cutIdx) {
        int trackID = static_cast<int>(kin_trackID_0[k]);
        if (trk_d0[trackID] < d0Cuts[cutIdx]) {
          h_delta_pos_trans_d0cut[cutIdx]->Fill(kin_deltapostrans_0[k]);
          h_delta_pos_long_d0cut[cutIdx]->Fill(kin_deltaposlong_0[k]);
        }
      }
      for (int cutIdx = 0; cutIdx < nTransCuts; ++cutIdx) {
        if (kin_deltapostrans_0[k] < transCuts[cutIdx]) {
          int trackID = static_cast<int>(kin_trackID_0[k]);
          if (trackID >= 0 && trackID < nTracks) {
            h_d0_all_trans[cutIdx]->Fill(trk_d0[trackID]);
            h_projz_all_trans[cutIdx]->Fill(trk_projz[trackID]);
            h_projy_all_trans[cutIdx]->Fill(trk_projy[trackID]);

            if (trk_d0_good[trackID] > 0) {
              h_d0_good_trans[cutIdx]->Fill(trk_d0[trackID]);
              h_projz_good_trans[cutIdx]->Fill(trk_projz[trackID]);
              h_projy_good_trans[cutIdx]->Fill(trk_projy[trackID]);
            } else {
              h_d0_bad_trans[cutIdx]->Fill(trk_d0[trackID]);
              h_projz_bad_trans[cutIdx]->Fill(trk_projz[trackID]);
              h_projy_bad_trans[cutIdx]->Fill(trk_projy[trackID]);
            }
          }
        }
      }

      for (int cutIdx = 0; cutIdx < nLongCuts; ++cutIdx) {
        if (kin_deltaposlong_0[k] < longCuts[cutIdx]) {
          int trackID = static_cast<int>(kin_trackID_0[k]);
          if (trackID >= 0 && trackID < nTracks) {
            h_d0_all_long[cutIdx]->Fill(trk_d0[trackID]);
            h_projz_all_long[cutIdx]->Fill(trk_projz[trackID]);
            h_projy_all_long[cutIdx]->Fill(trk_projy[trackID]);

            if (trk_d0_good[trackID] > 0) {
              h_d0_good_long[cutIdx]->Fill(trk_d0[trackID]);
              h_projz_good_long[cutIdx]->Fill(trk_projz[trackID]);
              h_projy_good_long[cutIdx]->Fill(trk_projy[trackID]);
            } else {
              h_d0_bad_long[cutIdx]->Fill(trk_d0[trackID]);
              h_projz_bad_long[cutIdx]->Fill(trk_projz[trackID]);
              h_projy_bad_long[cutIdx]->Fill(trk_projy[trackID]);
            }
          }
        }
      }
      if (abs(kin_deltaposlong_0[k]) < delta_pos_long_sig && abs(kin_deltapostrans_0[k]) < delta_pos_trans_sig) {
        int trackID = static_cast<int>(kin_trackID_0[k]);
        h_d0_all_box_cut->Fill(trk_d0[trackID]);
        h_projz_all_box_cut->Fill(trk_projz[trackID]);
        h_projz_vs_reactz_sig->Fill(trk_projz[trackID], reactz); // FIXME
        h_projy_all_box_cut->Fill(trk_projy[trackID]);
        h_hittime_box_cut->Fill(kin_hittime_0[k]);
        h_gem_time_all_box_cut->Fill(trk_time[trackID]);
        h_gem_dt_all_box_cut->Fill(trk_dt[trackID]);

        if (trk_d0_good[trackID] > 0) {
          h_d0_good_box_cut->Fill(trk_d0[trackID]);
          h_projz_good_box_cut->Fill(trk_projz[trackID]);
          h_projy_good_box_cut->Fill(trk_projy[trackID]);
          h_gem_time_good_box_cut->Fill(trk_time[trackID]);
          h_gem_dt_good_box_cut->Fill(trk_dt[trackID]);
        } else {
          h_d0_bad_box_cut->Fill(trk_d0[trackID]);
          h_projz_bad_box_cut->Fill(trk_projz[trackID]);
          h_projy_bad_box_cut->Fill(trk_projy[trackID]);
          h_gem_time_bad_box_cut->Fill(trk_time[trackID]);
          h_gem_dt_bad_box_cut->Fill(trk_dt[trackID]);
        }
      }
      // Fill the 2D histograms
      if (kin_hitedep_0[k] > 2) {
        // if (true) {
        int trackID = static_cast<int>(kin_trackID_0[k]);
        h_d0_vs_deltapostrans->Fill(trk_d0[trackID], kin_deltapostrans_0[k]);
        h_projz_vs_deltapostrans->Fill(trk_projz[k], kin_deltapostrans_0[k]);
        h_projy_vs_deltapostrans->Fill(trk_projy[k], kin_deltapostrans_0[k]);
        h_d0_vs_deltaposlong->Fill(trk_d0[trackID], kin_deltaposlong_0[k]);
        h_deltaposlong_vs_deltapostrans->Fill(kin_deltapostrans_0[k], kin_deltaposlong_0[k]);
      }
      // h_d0_vs_deltapostrans->Fill(trk_d0[k], kin_deltapostrans_0[k]);
      // h_projz_vs_deltapostrans->Fill(trk_projz[k], kin_deltapostrans_0[k]);
    }
    // Print the status as a percentage
    int progress = static_cast<int>((100.0 * i) / nEntries);
    std::cout << "\rProcessing: " << progress << "% completed" << std::flush;
  }
  std::cout << "\nProcessing completed." << std::endl;
  // End Event Loop

  // Create output file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
    return;
  }

  // Write histograms to file
  outputFile->mkdir("d0");
  outputFile->cd("d0");
  h_d0_all->Write();
  h_d0_good->Write();
  h_d0_bad->Write();

  outputFile->mkdir("projz");
  outputFile->cd("projz");
  h_projz_all->Write();
  h_projz_good->Write();
  h_projz_bad->Write();

  outputFile->mkdir("projy");
  outputFile->cd("projy");
  h_projy_all->Write();
  h_projy_good->Write();
  h_projy_bad->Write();

  outputFile->mkdir("lad_hits");
  outputFile->cd("lad_hits");
  h_delta_pos_long->Write();
  h_delta_pos_trans->Write();
  h_hittime->Write();

  for (int i = 0; i < nTransCuts; ++i) {
    outputFile->mkdir(Form("trans_cut_%.1f", transCuts[i]));
    outputFile->cd(Form("trans_cut_%.1f", transCuts[i]));
    h_d0_all_trans[i]->Write();
    h_d0_good_trans[i]->Write();
    h_d0_bad_trans[i]->Write();
    h_projz_all_trans[i]->Write();
    h_projz_good_trans[i]->Write();
    h_projz_bad_trans[i]->Write();
    h_projy_all_trans[i]->Write();
    h_projy_good_trans[i]->Write();
    h_projy_bad_trans[i]->Write();
  }
  for (int i = 0; i < nLongCuts; ++i) {
    outputFile->mkdir(Form("long_cut_%.1f", longCuts[i]));
    outputFile->cd(Form("long_cut_%.1f", longCuts[i]));
    h_d0_all_long[i]->Write();
    h_d0_good_long[i]->Write();
    h_d0_bad_long[i]->Write();
    h_projz_all_long[i]->Write();
    h_projz_good_long[i]->Write();
    h_projz_bad_long[i]->Write();
    h_projy_all_long[i]->Write();
    h_projy_good_long[i]->Write();
    h_projy_bad_long[i]->Write();
  }
  for (int i = 0; i < nD0Cuts; ++i) {
    outputFile->mkdir(Form("d0_cut_%.1f", d0Cuts[i]));
    outputFile->cd(Form("d0_cut_%.1f", d0Cuts[i]));
    h_delta_pos_long_d0cut[i]->Write();
    h_delta_pos_trans_d0cut[i]->Write();
  }
  outputFile->mkdir("box_cut");
  outputFile->cd("box_cut");
  h_d0_all_box_cut->Write();
  h_d0_good_box_cut->Write();
  h_d0_bad_box_cut->Write();
  // Create a canvas for projz and projy (All, Good, Bad)
  TCanvas *c_projz_projy_all = new TCanvas("c_projz_projy_all", "Projz and Projy All Box Cut", 1200, 600);
  c_projz_projy_all->Divide(2, 1);
  c_projz_projy_all->cd(1);
  h_projz_all_box_cut->Draw("HIST");
  c_projz_projy_all->cd(2);
  h_projy_all_box_cut->Draw("HIST");
  c_projz_projy_all->Write();

  TCanvas *c_projz_projy_good = new TCanvas("c_projz_projy_good", "Projz and Projy Good Box Cut", 1200, 600);
  c_projz_projy_good->Divide(2, 1);
  c_projz_projy_good->cd(1);
  h_projz_good_box_cut->Draw("HIST");
  c_projz_projy_good->cd(2);
  h_projy_good_box_cut->Draw("HIST");
  c_projz_projy_good->Write();

  TCanvas *c_projz_projy_bad = new TCanvas("c_projz_projy_bad", "Projz and Projy Bad Box Cut", 1200, 600);
  c_projz_projy_bad->Divide(2, 1);
  c_projz_projy_bad->cd(1);
  h_projz_bad_box_cut->Draw("HIST");
  c_projz_projy_bad->cd(2);
  h_projy_bad_box_cut->Draw("HIST");
  c_projz_projy_bad->Write();
  h_hittime_box_cut->Write();
  h_gem_time_all_box_cut->Write();
  h_gem_time_good_box_cut->Write();
  h_gem_time_bad_box_cut->Write();
  h_gem_dt_all_box_cut->Write();
  h_gem_dt_good_box_cut->Write();
  h_gem_dt_bad_box_cut->Write();

  // Create a canvas with 1x2 pads
  TCanvas *c_projz_vs_elz = new TCanvas("c_projz_vs_elz", "Projz vs ELZ with Boxes", 1200, 600);
  c_projz_vs_elz->Divide(2, 1);

  // Draw the 2D histogram on the left pad
  c_projz_vs_elz->cd(1);
  h_projz_vs_reactz_sig->Draw("COLZ");

  // Draw lines at x = ±5 and y = ±5
  TLine *line_x1 = new TLine(-5, projz_MIN, -5, projz_MAX);
  TLine *line_x2 = new TLine(5, projz_MIN, 5, projz_MAX);
  TLine *line_y1 = new TLine(projz_MIN, -5, projz_MAX, -5);
  TLine *line_y2 = new TLine(projz_MIN, 5, projz_MAX, 5);
  line_x1->SetLineColor(kRed);
  line_x2->SetLineColor(kRed);
  line_y1->SetLineColor(kRed);
  line_y2->SetLineColor(kRed);
  line_x1->Draw();
  line_x2->Draw();
  line_y1->Draw();
  line_y2->Draw();

  // Calculate the sum of events in each 3x3 box
  c_projz_vs_elz->cd(2);
  TH2D *h_box_sums =
      new TH2D("h_box_sums", "Sum of Events in each foil combo;ELZ (cm);Projz (cm)", 3, -15, 15, 3, -15, 15);

  for (int x_bin = 1; x_bin <= 3; ++x_bin) {
    for (int y_bin = 1; y_bin <= 3; ++y_bin) {
      double x_min = h_box_sums->GetXaxis()->GetBinLowEdge(x_bin);
      double x_max = h_box_sums->GetXaxis()->GetBinUpEdge(x_bin);
      double y_min = h_box_sums->GetYaxis()->GetBinLowEdge(y_bin);
      double y_max = h_box_sums->GetYaxis()->GetBinUpEdge(y_bin);

      int bin_x_min = h_projz_vs_reactz_sig->GetXaxis()->FindBin(x_min);
      int bin_x_max = h_projz_vs_reactz_sig->GetXaxis()->FindBin(x_max);
      int bin_y_min = h_projz_vs_reactz_sig->GetYaxis()->FindBin(y_min);
      int bin_y_max = h_projz_vs_reactz_sig->GetYaxis()->FindBin(y_max);

      double sum = h_projz_vs_reactz_sig->Integral(bin_x_min, bin_x_max, bin_y_min, bin_y_max);
      h_box_sums->SetBinContent(x_bin, y_bin, sum);
    }
  }

  h_box_sums->Draw("COLZ TEXT");

  // Write the canvas to the output file
  c_projz_vs_elz->Write();

  // Write canvases
  // Create a new directory for canvases
  outputFile->mkdir("canvases");
  outputFile->cd("canvases");

  // Create canvases and overlay plots for d0
  TCanvas *c_d0_all = new TCanvas("c_d0_all", "d0 All Distributions with Cuts", 800, 600);
  h_d0_all->SetLineColor(kBlack);
  h_d0_all->Draw("HIST");
  for (int i = 0; i < nTransCuts; ++i) {
    h_d0_all_trans[i]->SetLineColor(i + 2); // Different colors for each cut
    h_d0_all_trans[i]->Draw("HIST SAME");
  }
  TLegend *leg_d0_all = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_d0_all->AddEntry(h_d0_all, "No Cut", "l");
  for (int i = 0; i < nTransCuts; ++i) {
    leg_d0_all->AddEntry(h_d0_all_trans[i], Form("Cut %.1f", transCuts[i]), "l");
  }
  leg_d0_all->Draw();
  c_d0_all->Write();

  TCanvas *c_d0_good = new TCanvas("c_d0_good", "d0 Good Distributions with Cuts", 800, 600);
  h_d0_good->SetLineColor(kBlack);
  h_d0_good->Draw("HIST");
  for (int i = 0; i < nTransCuts; ++i) {
    h_d0_good_trans[i]->SetLineColor(i + 2);
    h_d0_good_trans[i]->Draw("HIST SAME");
  }
  TLegend *leg_d0_good = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_d0_good->AddEntry(h_d0_good, "No Cut", "l");
  for (int i = 0; i < nTransCuts; ++i) {
    leg_d0_good->AddEntry(h_d0_good_trans[i], Form("Cut %.1f", transCuts[i]), "l");
  }
  leg_d0_good->Draw();
  c_d0_good->Write();

  TCanvas *c_d0_bad = new TCanvas("c_d0_bad", "d0 Bad Distributions with Cuts", 800, 600);
  h_d0_bad->SetLineColor(kBlack);
  h_d0_bad->Draw("HIST");
  for (int i = 0; i < nTransCuts; ++i) {
    h_d0_bad_trans[i]->SetLineColor(i + 2);
    h_d0_bad_trans[i]->Draw("HIST SAME");
  }
  TLegend *leg_d0_bad = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_d0_bad->AddEntry(h_d0_bad, "No Cut", "l");
  for (int i = 0; i < nTransCuts; ++i) {
    leg_d0_bad->AddEntry(h_d0_bad_trans[i], Form("Cut %.1f", transCuts[i]), "l");
  }
  leg_d0_bad->Draw();
  c_d0_bad->Write();

  // Create canvases and overlay plots for projz
  TCanvas *c_projz_all = new TCanvas("c_projz_all", "projz All Distributions with Cuts", 800, 600);
  h_projz_all->SetLineColor(kBlack);
  h_projz_all->Draw("HIST");
  for (int i = 0; i < nTransCuts; ++i) {
    h_projz_all_trans[i]->SetLineColor(i + 2);
    h_projz_all_trans[i]->Draw("HIST SAME");
  }
  TLegend *leg_projz_all = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_projz_all->AddEntry(h_projz_all, "No Cut", "l");
  for (int i = 0; i < nTransCuts; ++i) {
    leg_projz_all->AddEntry(h_projz_all_trans[i], Form("Cut %.1f", transCuts[i]), "l");
  }
  leg_projz_all->Draw();
  c_projz_all->Write();

  TCanvas *c_projz_good = new TCanvas("c_projz_good", "projz Good Distributions with Cuts", 800, 600);
  h_projz_good->SetLineColor(kBlack);
  h_projz_good->Draw("HIST");
  for (int i = 0; i < nTransCuts; ++i) {
    h_projz_good_trans[i]->SetLineColor(i + 2);
    h_projz_good_trans[i]->Draw("HIST SAME");
  }
  TLegend *leg_projz_good = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_projz_good->AddEntry(h_projz_good, "No Cut", "l");
  for (int i = 0; i < nTransCuts; ++i) {
    leg_projz_good->AddEntry(h_projz_good_trans[i], Form("Cut %.1f", transCuts[i]), "l");
  }
  leg_projz_good->Draw();
  c_projz_good->Write();

  TCanvas *c_projz_bad = new TCanvas("c_projz_bad", "projz Bad Distributions with Cuts", 800, 600);
  h_projz_bad->SetLineColor(kBlack);
  h_projz_bad->Draw("HIST");
  for (int i = 0; i < nTransCuts; ++i) {
    h_projz_bad_trans[i]->SetLineColor(i + 2);
    h_projz_bad_trans[i]->Draw("HIST SAME");
  }
  TLegend *leg_projz_bad = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_projz_bad->AddEntry(h_projz_bad, "No Cut", "l");
  for (int i = 0; i < nTransCuts; ++i) {
    leg_projz_bad->AddEntry(h_projz_bad_trans[i], Form("Cut %.1f", transCuts[i]), "l");
  }
  leg_projz_bad->Draw();
  c_projz_bad->Write();

  // Create canvases and overlay plots for projy
  TCanvas *c_projy_all = new TCanvas("c_projy_all", "projy All Distributions with Cuts", 800, 600);
  h_projy_all->SetLineColor(kBlack);
  h_projy_all->Draw("HIST");
  for (int i = 0; i < nTransCuts; ++i) {
    h_projy_all_trans[i]->SetLineColor(i + 2);
    h_projy_all_trans[i]->Draw("HIST SAME");
  }
  TLegend *leg_projy_all = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_projy_all->AddEntry(h_projy_all, "No Cut", "l");
  for (int i = 0; i < nTransCuts; ++i) {
    leg_projy_all->AddEntry(h_projy_all_trans[i], Form("Cut %.1f", transCuts[i]), "l");
  }
  leg_projy_all->Draw();
  c_projy_all->Write();

  TCanvas *c_projy_good = new TCanvas("c_projy_good", "projy Good Distributions with Cuts", 800, 600);
  h_projy_good->SetLineColor(kBlack);
  h_projy_good->Draw("HIST");
  for (int i = 0; i < nTransCuts; ++i) {
    h_projy_good_trans[i]->SetLineColor(i + 2);
    h_projy_good_trans[i]->Draw("HIST SAME");
  }
  TLegend *leg_projy_good = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_projy_good->AddEntry(h_projy_good, "No Cut", "l");
  for (int i = 0; i < nTransCuts; ++i) {
    leg_projy_good->AddEntry(h_projy_good_trans[i], Form("Cut %.1f", transCuts[i]), "l");
  }
  leg_projy_good->Draw();
  c_projy_good->Write();

  TCanvas *c_projy_bad = new TCanvas("c_projy_bad", "projy Bad Distributions with Cuts", 800, 600);
  h_projy_bad->SetLineColor(kBlack);
  h_projy_bad->Draw("HIST");
  for (int i = 0; i < nTransCuts; ++i) {
    h_projy_bad_trans[i]->SetLineColor(i + 2);
    h_projy_bad_trans[i]->Draw("HIST SAME");
  }
  TLegend *leg_projy_bad = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_projy_bad->AddEntry(h_projy_bad, "No Cut", "l");
  for (int i = 0; i < nTransCuts; ++i) {
    leg_projy_bad->AddEntry(h_projy_bad_trans[i], Form("Cut %.1f", transCuts[i]), "l");
  }
  leg_projy_bad->Draw();
  c_projy_bad->Write();

  // Create directories and write 2D histograms
  outputFile->mkdir("2D_histograms");
  outputFile->cd("2D_histograms");
  h_d0_vs_deltapostrans->Write();
  h_projz_vs_deltapostrans->Write();
  h_projy_vs_deltapostrans->Write();
  h_d0_vs_deltaposlong->Write();
  h_deltaposlong_vs_deltapostrans->Write();

  // Close the output file
  outputFile->Close();
  delete outputFile;
  // Clean up
  file->Close();
  delete file;
}