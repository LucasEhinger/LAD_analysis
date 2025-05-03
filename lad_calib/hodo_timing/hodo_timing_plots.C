// Lucas Ehinger
// General LAD hodo plotting script
// Multithreaded, to speed up the process
// It's absolutely terrible when it comes to memory management, but it works
// Multi-threading with ROOT which is inherently not thread-safe is hard, and I wasn't smart enough to do it elegantly, but this works.
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include </usr/lib/gcc/x86_64-redhat-linux/11/include/omp.h>

using namespace std;

const int MAX_DATA     = 1000;
const int MAX_DATA_GEM = 1000;

const int d0_NBINS        = 100;
const double d0_MIN       = 0.0;
const double d0_MAX       = 100.0;
const int projz_NBINS     = 100;
const double projz_MIN    = -20.0;
const double projz_MAX    = 20.0;
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
const double edep_MAX = 100.0;
const int edep_NBINS  = 1000;

const int MINT_EVTS_PER_THREAD = 10000;

struct track_cut {
  double delta_pos_long_cut;
  double delta_pos_trans_cut;
  double d0_cut;
  string cut_name;
};

const int nTrackCuts             = 4;
track_cut track_cuts[nTrackCuts] = {{2000, 2000.0, 200.0, "loose"},
                                    {100.0, 2000.0, 200.0, "trans_50"},
                                    {100.0, 100.0, 100.0, "trans_50_long_100"},
                                    {100.0, 100.0, 30.0, "trans_50_long_100_d0_30"}};
// const double delta_pos_long_cut  = 200.0;
// const double delta_pos_trans_cut = 800.0;
// const double d0_cut              = 200.0;

const int maxTracks = 30;

const int N_PLANES                 = 5;
const int N_PADDLES                = 11;
const string plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const int N_SIDES                  = 2;
const string side_names[N_SIDES]   = {"Top", "Btm"};

template <typename HistType>
void write_to_canvas(HistType *hist_arr_bar[N_PLANES][N_SIDES][N_PADDLES], TFile *file, TString dir, TString var_name,
                     bool single_side = false) {
  // Create and navigate to the directory
  file->mkdir(dir);
  file->cd(dir);

  // Create a canvas for each plane and side
  HistType *hist_plane_sum[N_PLANES][N_SIDES] = {nullptr};
  for (int plane = 0; plane < N_PLANES; ++plane) {
    for (int side = 0; side < N_SIDES; ++side) {
      if (single_side && side == 1)
        continue; // Skip the second side if single_side is true
      // Create canvas for bars
      TCanvas *c1 = new TCanvas(Form("c_%s_plane_%d_%s", var_name.Data(), plane, side_names[side].c_str()),
                                Form("%s Plane %d %s", var_name.Data(), plane, side_names[side].c_str()), 800, 600);
      c1->Divide(4, 3); // Divide canvas into subpads
      for (int bar = 0; bar < N_PADDLES; ++bar) {
        c1->cd(bar + 1);
        hist_arr_bar[plane][side][bar]->Draw();
      }
      c1->Write();
      delete c1;

      // Sum all bars for this plane and side into a single histogram
      if (!hist_plane_sum[plane][side]) {
        hist_plane_sum[plane][side] = (HistType *)hist_arr_bar[plane][side][0]->Clone(
            Form("h_%s_plane_%d_%s_all_bars", var_name.Data(), plane, side_names[side].c_str()));
        hist_plane_sum[plane][side]->Reset();
        hist_plane_sum[plane][side]->SetTitle(
            Form("%s Plane %d %s All Bars", var_name.Data(), plane, side_names[side].c_str()));
      }
      for (int bar = 0; bar < N_PADDLES; ++bar) {
        hist_plane_sum[plane][side]->Add(hist_arr_bar[plane][side][bar]);
      }

      // Write the summed histogram for this plane and side
      // hist_plane_sum[plane][side]->Write();
    }
  }

  // Create a canvas to overlay all sides for each plane
  for (int plane = 0; plane < N_PLANES; ++plane) {
    TCanvas *c_all_sides = new TCanvas(Form("c_%s_plane_%d_all_sides", var_name.Data(), plane),
                                       Form("%s Plane %d All Sides", var_name.Data(), plane), 1200, 800);
    c_all_sides->Divide(2, 1); // Divide canvas into subpads for sides
    for (int side = 0; side < N_SIDES; ++side) {
      if (single_side && side == 1)
        continue; // Skip the second side if single_side is true
      c_all_sides->cd(side + 1);
      hist_plane_sum[plane][side]->SetLineColor(side + 1); // Assign different colors for each side
      hist_plane_sum[plane][side]->Draw("HIST");
    }
    c_all_sides->Write();
    delete c_all_sides;
  }

  // Sum all plane and side histograms into one final histogram
  HistType *hist_final_sum = (HistType *)hist_plane_sum[0][0]->Clone(Form("h_%s_final_sum", var_name.Data()));
  hist_final_sum->Reset();
  hist_final_sum->SetTitle(Form("%s All Planes", var_name.Data()));
  for (int plane = 0; plane < N_PLANES; ++plane) {
    for (int side = 0; side < N_SIDES; ++side) {
      if (single_side && side == 1)
        continue; // Skip the second side if single_side is true
      hist_final_sum->Add(hist_plane_sum[plane][side]);
    }
  }

  // Write the final summed histogram
  hist_final_sum->Write();
  delete hist_final_sum;
}

template <typename HistType>
void write_to_canvas_plane(HistType *hist_arr[N_PLANES][N_PADDLES], TFile *file, TString dir, TString var_name) {
  // Create and navigate to the directory
  file->mkdir(dir);
  file->cd(dir);
  // Create a canvas for each plane
  TH1F *hist_plane_sum[N_PLANES] = {nullptr};
  for (int plane = 0; plane < N_PLANES; ++plane) {
    // Create canvas for top bars
    TCanvas *c1 = new TCanvas(Form("c_%s_plane_%d", var_name.Data(), plane),
                              Form("%s Plane %d", var_name.Data(), plane), 800, 600);
    c1->Divide(4, 3); // Divide canvas into subpads
    for (int bar = 0; bar < N_PADDLES; ++bar) {
      c1->cd(bar + 1);
      hist_arr[plane][bar]->Draw();
    }
    c1->Write();
    delete c1;

    // Sum all bars for this plane into a single histogram
    if (!hist_plane_sum[plane]) {
      hist_plane_sum[plane] = (TH1F *)hist_arr[plane][0]->Clone(Form("h_%s_plane_%d_all_bars", var_name.Data(), plane));
      hist_plane_sum[plane]->Reset();
      hist_plane_sum[plane]->SetTitle(Form("%s Plane %d All Bars", var_name.Data(), plane));
    }
    for (int bar = 0; bar < N_PADDLES; ++bar) {
      hist_plane_sum[plane]->Add(hist_arr[plane][bar]);
    }

    // Write the summed histogram for this plane
    // hist_plane_sum[plane]->Write();
  }

  // Create a canvas to overlay all planes
  TCanvas *c_all_planes =
      new TCanvas(Form("c_%s_all_planes", var_name.Data()), Form("%s All Planes", var_name.Data()), 1200, 800);
  c_all_planes->Divide(2, 3); // Divide canvas into subpads
  for (int plane = 0; plane < N_PLANES; ++plane) {
    c_all_planes->cd(plane + 1);
    hist_plane_sum[plane]->SetLineColor(plane + 1); // Assign different colors for each plane
    hist_plane_sum[plane]->Draw("HIST");
  }
  c_all_planes->Write();
  delete c_all_planes;

  // Sum all plane histograms into one final histogram
  TH1F *hist_final_sum = (TH1F *)hist_plane_sum[0]->Clone(Form("h_%s_final_sum", var_name.Data()));
  hist_final_sum->Reset();
  hist_final_sum->SetTitle(Form("%s All Planes", var_name.Data()));
  for (int plane = 0; plane < N_PLANES; ++plane) {
    hist_final_sum->Add(hist_plane_sum[plane]);
  }

  // Write the final summed histogram
  hist_final_sum->Write();
  delete hist_final_sum;
}

void process_chunk(int i_thread, int start, int end, TString fileName,
                   map<string, TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]> &hist_map,
                   map<string, TH1 *(*)[nTrackCuts][N_PLANES][N_PADDLES]> &hist_map_trk_cuts) {

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
  Double_t trk_d0[MAX_DATA_GEM];
  Double_t trk_d0_good[MAX_DATA_GEM];
  Double_t trk_projz[MAX_DATA_GEM];
  Double_t kin_trackID_0[MAX_DATA_GEM], kin_trackID_1[MAX_DATA_GEM];
  Double_t kin_plane_0[MAX_DATA_GEM], kin_plane_1[MAX_DATA_GEM];
  Double_t kin_paddle_0[MAX_DATA_GEM], kin_paddle_1[MAX_DATA_GEM];
  Double_t kin_hittime_0[MAX_DATA_GEM], kin_hittime_1[MAX_DATA_GEM];
  Double_t kin_hittheta_0[MAX_DATA_GEM], kin_hittheta_1[MAX_DATA_GEM];
  Double_t kin_hitphi_0[MAX_DATA_GEM], kin_hitphi_1[MAX_DATA_GEM];
  Double_t kin_hitedep_0[MAX_DATA_GEM], kin_hitedep_1[MAX_DATA_GEM];
  Double_t kin_deltapostrans_0[MAX_DATA_GEM], kin_deltapostrans_1[MAX_DATA_GEM];
  Double_t kin_deltaposlong_0[MAX_DATA_GEM], kin_deltaposlong_1[MAX_DATA_GEM];
  Int_t nTracks, nGoodHits;

  // Hodo-level data
  Double_t tdc_time[N_PLANES][N_SIDES][MAX_DATA];
  Double_t tdc_counter[N_PLANES][N_SIDES][MAX_DATA];
  Double_t tdc_time_UnCorr[N_PLANES][N_SIDES][MAX_DATA], tdc_time_TWCorr[N_PLANES][N_SIDES][MAX_DATA],
      tdc_time_Corr[N_PLANES][N_SIDES][MAX_DATA];
  Double_t fullhit_time_avg[N_PLANES][MAX_DATA], fullhit_adc_avg[N_PLANES][MAX_DATA],
      fullhit_paddle[N_PLANES][MAX_DATA];
  Int_t fullhit_n[N_PLANES];

  Double_t adc_time[N_PLANES][N_SIDES][MAX_DATA];
  Double_t adc_counter[N_PLANES][N_SIDES][MAX_DATA];
  Double_t adc_time_good[N_PLANES][N_SIDES][MAX_DATA];
  Double_t adc_amp[N_PLANES][N_SIDES][MAX_DATA];
  Double_t adc_amp_good[N_PLANES][N_SIDES][MAX_DATA];
  Double_t adc_int[N_PLANES][N_SIDES][MAX_DATA];
  Double_t adc_int_good[N_PLANES][N_SIDES][MAX_DATA];

  Int_t nData_adc[N_PLANES][N_SIDES];
  Int_t nData_tdc[N_PLANES][N_SIDES];
  Double_t hodo_start_time, pTRIG1, pTRIG2, pTRIG3, pTRIG4;

  T->SetBranchAddress("Ndata.H.gem.trk.d0", &nTracks);
  T->SetBranchAddress("Ndata.H.ladkin.goodhit_trackid_0", &nGoodHits);
  T->SetBranchAddress("H.gem.trk.d0", &trk_d0);
  T->SetBranchAddress("H.gem.trk.d0_good", &trk_d0_good);
  T->SetBranchAddress("H.gem.trk.projz", &trk_projz);
  T->SetBranchAddress("H.ladkin.goodhit_trackid_0", &kin_trackID_0);
  T->SetBranchAddress("H.ladkin.goodhit_trackid_1", &kin_trackID_1);
  T->SetBranchAddress("H.ladkin.goodhit_plane_0", &kin_plane_0);
  T->SetBranchAddress("H.ladkin.goodhit_plane_1", &kin_plane_1);
  T->SetBranchAddress("H.ladkin.goodhit_paddle_0", &kin_paddle_0);
  T->SetBranchAddress("H.ladkin.goodhit_paddle_1", &kin_paddle_1);
  T->SetBranchAddress("H.ladkin.goodhit_hittime_0", &kin_hittime_0);
  T->SetBranchAddress("H.ladkin.goodhit_hittime_1", &kin_hittime_1);
  T->SetBranchAddress("H.ladkin.goodhit_hittheta_0", &kin_hittheta_0);
  T->SetBranchAddress("H.ladkin.goodhit_hittheta_1", &kin_hittheta_1);
  T->SetBranchAddress("H.ladkin.goodhit_hitphi_0", &kin_hitphi_0);
  T->SetBranchAddress("H.ladkin.goodhit_hitphi_1", &kin_hitphi_1);
  T->SetBranchAddress("H.ladkin.goodhit_hitedep_0", &kin_hitedep_0);
  T->SetBranchAddress("H.ladkin.goodhit_hitedep_1", &kin_hitedep_1);
  T->SetBranchAddress("H.ladkin.goodhit_deltapostrans_0", &kin_deltapostrans_0);
  T->SetBranchAddress("H.ladkin.goodhit_deltapostrans_1", &kin_deltapostrans_1);
  T->SetBranchAddress("H.ladkin.goodhit_deltaposlong_0", &kin_deltaposlong_0);
  T->SetBranchAddress("H.ladkin.goodhit_deltaposlong_1", &kin_deltaposlong_1);

  T->SetBranchAddress("H.hod.starttime", &hodo_start_time);
  T->SetBranchAddress("T.hms.pTRIG1_tdcTime", &pTRIG1);
  T->SetBranchAddress("T.hms.pTRIG2_tdcTime", &pTRIG2);
  T->SetBranchAddress("T.hms.pTRIG3_tdcTime", &pTRIG3);
  T->SetBranchAddress("T.hms.pTRIG4_tdcTime", &pTRIG4);

  for (int plane = 0; plane < N_PLANES; ++plane) {
    for (int side = 0; side < N_SIDES; ++side) {
      T->SetBranchAddress(Form("H.ladhod.%s.%sTdcTime", plane_names[plane].c_str(), side_names[side].c_str()),
                          &tdc_time[plane][side]);
      T->SetBranchAddress(Form("H.ladhod.%s.Good%sTdcTimeUnCorr", plane_names[plane].c_str(), side_names[side].c_str()),
                          &tdc_time_UnCorr[plane][side]);
      T->SetBranchAddress(
          Form("H.ladhod.%s.Good%sTdcTimeWalkCorr", plane_names[plane].c_str(), side_names[side].c_str()),
          &tdc_time_TWCorr[plane][side]);
      T->SetBranchAddress(Form("H.ladhod.%s.Good%sTdcTimeCorr", plane_names[plane].c_str(), side_names[side].c_str()),
                          &tdc_time_Corr[plane][side]);
      T->SetBranchAddress(Form("H.ladhod.%s.%sAdcPulseTime", plane_names[plane].c_str(), side_names[side].c_str()),
                          &adc_time[plane][side]);
      T->SetBranchAddress(Form("H.ladhod.%s.Good%sAdcPulseTime", plane_names[plane].c_str(), side_names[side].c_str()),
                          &adc_time_good[plane][side]);
      T->SetBranchAddress(Form("H.ladhod.%s.%sAdcPulseAmp", plane_names[plane].c_str(), side_names[side].c_str()),
                          &adc_amp[plane][side]);
      T->SetBranchAddress(Form("H.ladhod.%s.Good%sAdcPulseAmp", plane_names[plane].c_str(), side_names[side].c_str()),
                          &adc_amp_good[plane][side]);
      T->SetBranchAddress(Form("H.ladhod.%s.%sAdcPulseInt", plane_names[plane].c_str(), side_names[side].c_str()),
                          &adc_int[plane][side]);
      T->SetBranchAddress(Form("H.ladhod.%s.Good%sAdcPulseInt", plane_names[plane].c_str(), side_names[side].c_str()),
                          &adc_int_good[plane][side]);
      T->SetBranchAddress(Form("H.ladhod.%s.%sTdcCounter", plane_names[plane].c_str(), side_names[side].c_str()),
                          &tdc_counter[plane][side]);
      T->SetBranchAddress(Form("H.ladhod.%s.%sAdcCounter", plane_names[plane].c_str(), side_names[side].c_str()),
                          &adc_counter[plane][side]);
      T->SetBranchAddress(Form("Ndata.H.ladhod.%s.%sAdcCounter", plane_names[plane].c_str(), side_names[side].c_str()),
                          &nData_adc[plane][side]);
      T->SetBranchAddress(Form("Ndata.H.ladhod.%s.%sTdcCounter", plane_names[plane].c_str(), side_names[side].c_str()),
                          &nData_tdc[plane][side]);
    }
    T->SetBranchAddress(Form("H.ladhod.%s.HodoHitTime", plane_names[plane].c_str()), &fullhit_time_avg[plane]);
    T->SetBranchAddress(Form("H.ladhod.%s.HodoHitEdep", plane_names[plane].c_str()), &fullhit_adc_avg[plane]);
    T->SetBranchAddress(Form("H.ladhod.%s.HodoHitPaddleNum", plane_names[plane].c_str()), &fullhit_paddle[plane]);
    T->SetBranchAddress(Form("Ndata.H.ladhod.%s.HodoHitTime", plane_names[plane].c_str()), &fullhit_n[plane]);
  }
  ////////////////////////////////////////////////////
  // Start Event Loop

  for (int i = start; i < end; ++i) {
    T->GetEntry(i);

    for (int plane = 0; plane < N_PLANES; ++plane) {
      for (int side = 0; side < N_SIDES; ++side) {

        // Fill Raw Values
        for (int j = 0; j < nData_tdc[plane][side]; ++j) {
          hist_map["h_tdc_bar"][0][plane][side][int(tdc_counter[plane][side][j] - 1)]->Fill(tdc_time[plane][side][j]);
        }
        for (int j = 0; j < nData_adc[plane][side]; ++j) {
          hist_map["h_adc_bar"][0][plane][side][int(adc_counter[plane][side][j] - 1)]->Fill(adc_time[plane][side][j]);
          hist_map["h_adc_amp_bar"][0][plane][side][int(adc_counter[plane][side][j] - 1)]->Fill(
              adc_amp[plane][side][j]);
          hist_map["h_adc_int_bar"][0][plane][side][int(adc_counter[plane][side][j] - 1)]->Fill(
              adc_int[plane][side][j]);
        }

        // Fill Good Values
        for (int j = 0; j < N_PADDLES; ++j) {

          hist_map["h_adc_good_bar"][0][plane][side][j]->Fill(adc_time_good[plane][side][j]);
          hist_map["h_adc_amp_good_bar"][0][plane][side][j]->Fill(adc_amp_good[plane][side][j]);
          hist_map["h_adc_int_bar"][0][plane][side][j]->Fill(adc_int[plane][side][j]);
          hist_map["h_adc_int_good_bar"][0][plane][side][j]->Fill(adc_int_good[plane][side][j]);
          hist_map["h_tdc_UnCorr_bar"][0][plane][side][j]->Fill(tdc_time_UnCorr[plane][side][j]);
          hist_map["h_tdc_TWCorr_bar"][0][plane][side][j]->Fill(tdc_time_TWCorr[plane][side][j]);
          hist_map["h_tdc_Corr_bar"][0][plane][side][j]->Fill(tdc_time_Corr[plane][side][j]);
        }
      }
      for (int j = 0; j < fullhit_n[plane]; j++) {
        hist_map["h_tdc_avg_bar"][0][plane][0][int(fullhit_paddle[plane][j] - 1)]->Fill(fullhit_time_avg[plane][j]);
        hist_map["h_edep_avg_bar"][0][plane][0][int(fullhit_paddle[plane][j] - 1)]->Fill(fullhit_adc_avg[plane][j]);
        hist_map["h_tdc_avg_vs_edep_bar"][0][plane][0][int(fullhit_paddle[plane][j] - 1)]->Fill(
            fullhit_time_avg[plane][j], fullhit_adc_avg[plane][j]);
      }
      for (int i_goodhit = 0; i_goodhit < nGoodHits; i_goodhit++) {
        for (int i_track_cut = 0; i_track_cut < nTrackCuts; i_track_cut++) {
          if (abs(kin_deltaposlong_0[i_goodhit]) < track_cuts[i_track_cut].delta_pos_long_cut &&
              abs(kin_deltapostrans_0[i_goodhit]) < track_cuts[i_track_cut].delta_pos_trans_cut &&
              trk_d0[int(kin_trackID_0[i_goodhit])] < track_cuts[i_track_cut].d0_cut) {
            hist_map_trk_cuts["h_ladkin_time_bar"][0][i_track_cut][int(kin_plane_0[i_goodhit])]
                             [int(kin_paddle_0[i_goodhit])]
                                 ->Fill(kin_hittime_0[i_goodhit]);
            hist_map_trk_cuts["h_ladkin_edep_bar"][0][i_track_cut][int(kin_plane_0[i_goodhit])]
                             [int(kin_paddle_0[i_goodhit])]
                                 ->Fill(kin_hitedep_0[i_goodhit]);
            hist_map_trk_cuts["h_ladkin_theta_bar"][0][i_track_cut][int(kin_plane_0[i_goodhit])]
                             [int(kin_paddle_0[i_goodhit])]
                                 ->Fill(kin_hittheta_0[i_goodhit]);
            hist_map_trk_cuts["h_ladkin_phi_bar"][0][i_track_cut][int(kin_plane_0[i_goodhit])]
                             [int(kin_paddle_0[i_goodhit])]
                                 ->Fill(kin_hitphi_0[i_goodhit]);
            hist_map_trk_cuts["h_ladkin_deltaposlong_bar"][0][i_track_cut][int(kin_plane_0[i_goodhit])]
                             [int(kin_paddle_0[i_goodhit])]
                                 ->Fill(kin_deltaposlong_0[i_goodhit]);
            hist_map_trk_cuts["h_ladkin_deltapostrans_bar"][0][i_track_cut][int(kin_plane_0[i_goodhit])]
                             [int(kin_paddle_0[i_goodhit])]
                                 ->Fill(kin_deltapostrans_0[i_goodhit]);
            hist_map_trk_cuts["h_ladkin_d0_bar"][0][i_track_cut][int(kin_plane_0[i_goodhit])]
                             [int(kin_paddle_0[i_goodhit])]
                                 ->Fill(trk_d0_good[int(kin_trackID_0[i_goodhit])]);
          }
          if (abs(kin_deltaposlong_1[i_goodhit]) < track_cuts[i_track_cut].delta_pos_long_cut &&
              abs(kin_deltapostrans_1[i_goodhit]) < track_cuts[i_track_cut].delta_pos_trans_cut &&
              trk_d0[int(kin_trackID_1[i_goodhit])] < track_cuts[i_track_cut].d0_cut) {
            hist_map_trk_cuts["h_ladkin_time_bar"][0][i_track_cut][int(kin_plane_1[i_goodhit])]
                             [int(kin_paddle_1[i_goodhit])]
                                 ->Fill(kin_hittime_1[i_goodhit]);
            hist_map_trk_cuts["h_ladkin_edep_bar"][0][i_track_cut][int(kin_plane_1[i_goodhit])]
                             [int(kin_paddle_1[i_goodhit])]
                                 ->Fill(kin_hitedep_1[i_goodhit]);
            hist_map_trk_cuts["h_ladkin_theta_bar"][0][i_track_cut][int(kin_plane_1[i_goodhit])]
                             [int(kin_paddle_1[i_goodhit])]
                                 ->Fill(kin_hittheta_1[i_goodhit]);
            hist_map_trk_cuts["h_ladkin_phi_bar"][0][i_track_cut][int(kin_plane_1[i_goodhit])]
                             [int(kin_paddle_1[i_goodhit])]
                                 ->Fill(kin_hitphi_1[i_goodhit]);
            hist_map_trk_cuts["h_ladkin_deltaposlong_bar"][0][i_track_cut][int(kin_plane_1[i_goodhit])]
                             [int(kin_paddle_1[i_goodhit])]
                                 ->Fill(kin_deltaposlong_1[i_goodhit]);
            hist_map_trk_cuts["h_ladkin_deltapostrans_bar"][0][i_track_cut][int(kin_plane_1[i_goodhit])]
                             [int(kin_paddle_1[i_goodhit])]
                                 ->Fill(kin_deltapostrans_1[i_goodhit]);
            hist_map_trk_cuts["h_ladkin_d0_bar"][0][i_track_cut][int(kin_plane_1[i_goodhit])]
                             [int(kin_paddle_1[i_goodhit])]
                                 ->Fill(trk_d0_good[int(kin_trackID_1[i_goodhit])]);
          }
        }
      }
    } // End Plane Loop

    // Print the status as a percentage
    if (i % ((end - start) / 100) == 0 && i_thread == 0) {
      std::cout << "\rProcessing: " << int((i - start) * 100.0 / (end - start)) << "% completed." << std::flush;
    }
  } // End Event Loop

  return;
}

void hodo_timing_plots() {
  // Set batch mode to suppress graphical output
  gROOT->SetBatch(kTRUE);
  ROOT::EnableThreadSafety();

  // Open the ROOT file
  TString fileName       = "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/"
                          //  "LAD_COIN_22282_500000.root";
                          "LAD_COIN_22382_0_5_500000.root";
  TString outputFileName = "hodo_timing_plots_22382.root";
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

  // Number of entries in the TTree
  int nEntries = T->GetEntries();

  // Number of threads to use
  int numThreads = std::thread::hardware_concurrency();
  int chunkSize  = nEntries / numThreads;

  // Adjust the number of threads if the chunk size is too small
  if (chunkSize < MINT_EVTS_PER_THREAD) {
    numThreads = std::max(1, nEntries / MINT_EVTS_PER_THREAD);
    chunkSize  = nEntries / numThreads;
  }

  // map<string, TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]> hist_map;
  // map<string, TH1 *(*)[nTrackCuts][N_PLANES][N_PADDLES]> hist_map_trk_cuts;

  std::vector<map<string, TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>> hist_map_vec(numThreads);
  std::vector<map<string, TH1 *(*)[nTrackCuts][N_PLANES][N_PADDLES]>> hist_map_trk_cuts_vec(numThreads);

  // Create histograms for each bar in each plane

  TH1F *h_tdc_bar[numThreads][N_PLANES][N_SIDES][N_PADDLES];
  TH1F *h_adc_bar[numThreads][N_PLANES][N_SIDES][N_PADDLES];
  TH1F *h_adc_good_bar[numThreads][N_PLANES][N_SIDES][N_PADDLES];
  TH1F *h_adc_amp_bar[numThreads][N_PLANES][N_SIDES][N_PADDLES];
  TH1F *h_adc_amp_good_bar[numThreads][N_PLANES][N_SIDES][N_PADDLES];
  TH1F *h_adc_int_bar[numThreads][N_PLANES][N_SIDES][N_PADDLES];
  TH1F *h_adc_int_good_bar[numThreads][N_PLANES][N_SIDES][N_PADDLES];
  TH1F *h_tdc_UnCorr_bar[numThreads][N_PLANES][N_SIDES][N_PADDLES];
  TH1F *h_tdc_TWCorr_bar[numThreads][N_PLANES][N_SIDES][N_PADDLES];
  TH1F *h_tdc_Corr_bar[numThreads][N_PLANES][N_SIDES][N_PADDLES];

  TH1F *h_tdc_avg_bar[numThreads][N_PLANES][N_SIDES][N_PADDLES]; // N_SIDES not used. Just use 0 index
  TH1F *h_edep_avg_bar[numThreads][N_PLANES][N_SIDES][N_PADDLES];
  TH2F *h_tdc_avg_vs_edep_bar[numThreads][N_PLANES][N_SIDES][N_PADDLES];

  // LADKIN Histograms
  TH1F *h_ladkin_time_bar[numThreads][nTrackCuts][N_PLANES][N_PADDLES];
  TH1F *h_ladkin_edep_bar[numThreads][nTrackCuts][N_PLANES][N_PADDLES];
  TH1F *h_ladkin_theta_bar[numThreads][nTrackCuts][N_PLANES][N_PADDLES];
  TH1F *h_ladkin_phi_bar[numThreads][nTrackCuts][N_PLANES][N_PADDLES];
  TH1F *h_ladkin_deltaposlong_bar[numThreads][nTrackCuts][N_PLANES][N_PADDLES];
  TH1F *h_ladkin_deltapostrans_bar[numThreads][nTrackCuts][N_PLANES][N_PADDLES];
  TH1F *h_ladkin_d0_bar[numThreads][nTrackCuts][N_PLANES][N_PADDLES];
  for (int thread = 0; thread < numThreads; ++thread) {
    for (int plane = 0; plane < N_PLANES; ++plane) {
      for (int bar = 0; bar < N_PADDLES; ++bar) {
        for (int side = 0; side < N_SIDES; ++side) {
          h_tdc_bar[thread][plane][side][bar] =
              new TH1F(Form("h_tdc_%s_plane_%d_bar_%d_thread_%d", side_names[side].c_str(), plane, bar, thread),
                       Form("TDC %s Plane %d Bar %d Thread %d", side_names[side].c_str(), plane, bar, thread),
                       time_NBINS, time_MIN, time_MAX);
          h_adc_bar[thread][plane][side][bar] =
              new TH1F(Form("h_adc_%s_plane_%d_bar_%d_thread_%d", side_names[side].c_str(), plane, bar, thread),
                       Form("ADC %s Plane %d Bar %d Thread %d", side_names[side].c_str(), plane, bar, thread),
                       time_NBINS, time_MIN, time_MAX);
          h_adc_good_bar[thread][plane][side][bar] =
              new TH1F(Form("h_adc_good_%s_plane_%d_bar_%d_thread_%d", side_names[side].c_str(), plane, bar, thread),
                       Form("ADC Good %s Plane %d Bar %d Thread %d", side_names[side].c_str(), plane, bar, thread),
                       time_NBINS, time_MIN, time_MAX);
          h_adc_amp_bar[thread][plane][side][bar] =
              new TH1F(Form("h_adc_amp_%s_plane_%d_bar_%d_thread_%d", side_names[side].c_str(), plane, bar, thread),
                       Form("ADC Amp %s Plane %d Bar %d Thread %d", side_names[side].c_str(), plane, bar, thread),
                       edep_NBINS, edep_MIN, edep_MAX);
          h_adc_amp_good_bar[thread][plane][side][bar] = new TH1F(
              Form("h_adc_amp_good_%s_plane_%d_bar_%d_thread_%d", side_names[side].c_str(), plane, bar, thread),
              Form("ADC Amp Good %s Plane %d Bar %d Thread %d", side_names[side].c_str(), plane, bar, thread),
              edep_NBINS, edep_MIN, edep_MAX);
          h_adc_int_bar[thread][plane][side][bar] =
              new TH1F(Form("h_adc_int_%s_plane_%d_bar_%d_thread_%d", side_names[side].c_str(), plane, bar, thread),
                       Form("ADC Int %s Plane %d Bar %d Thread %d", side_names[side].c_str(), plane, bar, thread),
                       edep_NBINS, edep_MIN, edep_MAX);
          h_adc_int_good_bar[thread][plane][side][bar] = new TH1F(
              Form("h_adc_int_good_%s_plane_%d_bar_%d_thread_%d", side_names[side].c_str(), plane, bar, thread),
              Form("ADC Int Good %s Plane %d Bar %d Thread %d", side_names[side].c_str(), plane, bar, thread),
              edep_NBINS, edep_MIN, edep_MAX);
          h_tdc_UnCorr_bar[thread][plane][side][bar] =
              new TH1F(Form("h_tdc_UnCorr_%s_plane_%d_bar_%d_thread_%d", side_names[side].c_str(), plane, bar, thread),
                       Form("TDC UnCorr %s Plane %d Bar %d Thread %d", side_names[side].c_str(), plane, bar, thread),
                       time_NBINS, time_MIN, time_MAX);
          h_tdc_TWCorr_bar[thread][plane][side][bar] =
              new TH1F(Form("h_tdc_TWCorr_%s_plane_%d_bar_%d_thread_%d", side_names[side].c_str(), plane, bar, thread),
                       Form("TDC TWCorr %s Plane %d Bar %d Thread %d", side_names[side].c_str(), plane, bar, thread),
                       time_NBINS, time_MIN, time_MAX);
          h_tdc_Corr_bar[thread][plane][side][bar] =
              new TH1F(Form("h_tdc_Corr_%s_plane_%d_bar_%d_thread_%d", side_names[side].c_str(), plane, bar, thread),
                       Form("TDC Corr %s Plane %d Bar %d Thread %d", side_names[side].c_str(), plane, bar, thread),
                       time_NBINS, time_MIN, time_MAX);

          // Initializing side 0 and 1 for following histograms to prevent segfaults. Only side 0 is used (data values
          // are avg of top and btm)
          if (side == 1) {
            h_tdc_avg_bar[thread][plane][side][bar] =
                new TH1F(Form("h_tdc_avg_dummy_plane_%d_bar_%d_thread_%d", plane, bar, thread),
                         Form("Dummy TDC Avg Plane %d Bar %d Thread %d", plane, bar, thread), 1, 0, 1);
            h_edep_avg_bar[thread][plane][side][bar] =
                new TH1F(Form("h_edep_avg_dummy_plane_%d_bar_%d_thread_%d", plane, bar, thread),
                         Form("Dummy EDEP Avg Plane %d Bar %d Thread %d", plane, bar, thread), 1, 0, 1);
            h_tdc_avg_vs_edep_bar[thread][plane][side][bar] =
                new TH2F(Form("h_tdc_avg_vs_edep_dummy_plane_%d_bar_%d_thread_%d", plane, bar, thread),
                         Form("Dummy TDC Avg vs EDEP Plane %d Bar %d Thread %d", plane, bar, thread), 1, 0, 1, 1, 0, 1);
          } else {
            h_tdc_avg_bar[thread][plane][side][bar] =
                new TH1F(Form("h_tdc_avg_plane_%d_bar_%d_thread_%d", plane, bar, thread),
                         Form("TDC Avg Plane %d Bar %d Thread %d", plane, bar, thread), time_NBINS, time_MIN, time_MAX);
            h_edep_avg_bar[thread][plane][side][bar] = new TH1F(
                Form("h_edep_avg_plane_%d_bar_%d_thread_%d", plane, bar, thread),
                Form("EDEP Avg Plane %d Bar %d Thread %d", plane, bar, thread), edep_NBINS, edep_MIN, edep_MAX);
            h_tdc_avg_vs_edep_bar[thread][plane][side][bar] =
                new TH2F(Form("h_tdc_avg_vs_edep_plane_%d_bar_%d_thread_%d", plane, bar, thread),
                         Form("TDC Avg vs EDEP Plane %d Bar %d Thread %d", plane, bar, thread), time_NBINS, time_MIN,
                         time_MAX, edep_NBINS, edep_MIN, edep_MAX);
          }
        }
        for (int i_trk_cut = 0; i_trk_cut < nTrackCuts; i_trk_cut++) {
          h_ladkin_time_bar[thread][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_time_plane_%d_bar_%d_%s_thread_%d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN Time Plane %d Bar %d %s Thread %d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       time_NBINS, time_MIN, time_MAX);
          h_ladkin_edep_bar[thread][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_edep_plane_%d_bar_%d_%s_thread_%d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN EDEP Plane %d Bar %d %s Thread %d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       edep_NBINS, edep_MIN, edep_MAX);
          h_ladkin_theta_bar[thread][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_theta_plane_%d_bar_%d_%s_thread_%d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN Theta Plane %d Bar %d %s Thread %d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       theta_NBINS, theta_MIN, theta_MAX);
          h_ladkin_phi_bar[thread][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_phi_plane_%d_bar_%d_%s_thread_%d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN Phi Plane %d Bar %d %s Thread %d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       phi_NBINS, phi_MIN, phi_MAX);
          h_ladkin_deltaposlong_bar[thread][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_deltaposlong_plane_%d_bar_%d_%s_thread_%d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN Delta Pos Long Plane %d Bar %d %s Thread %d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       deltapos_NBINS, deltapos_MIN, deltapos_MAX);
          h_ladkin_deltapostrans_bar[thread][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_deltapostrans_plane_%d_bar_%d_%s_thread_%d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN Delta Pos Trans Plane %d Bar %d %s Thread %d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       deltapos_NBINS, deltapos_MIN, deltapos_MAX);
          h_ladkin_d0_bar[thread][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_d0_plane_%d_bar_%d_%s_thread_%d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN D0 Plane %d Bar %d %s Thread %d", plane, bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       d0_NBINS, d0_MIN, d0_MAX);
        }
      }
    }
  }

  for (int thread = 0; thread < numThreads; ++thread) {
    hist_map_vec[thread]["h_tdc_bar"] = reinterpret_cast<TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>(h_tdc_bar[thread]);
    hist_map_vec[thread]["h_adc_bar"] = reinterpret_cast<TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>(h_adc_bar[thread]);
    hist_map_vec[thread]["h_adc_good_bar"] =
        reinterpret_cast<TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>(h_adc_good_bar[thread]);
    hist_map_vec[thread]["h_adc_amp_bar"] =
        reinterpret_cast<TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>(h_adc_amp_bar[thread]);
    hist_map_vec[thread]["h_adc_amp_good_bar"] =
        reinterpret_cast<TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>(h_adc_amp_good_bar[thread]);
    hist_map_vec[thread]["h_adc_int_bar"] =
        reinterpret_cast<TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>(h_adc_int_bar[thread]);
    hist_map_vec[thread]["h_adc_int_good_bar"] =
        reinterpret_cast<TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>(h_adc_int_good_bar[thread]);
    hist_map_vec[thread]["h_tdc_UnCorr_bar"] =
        reinterpret_cast<TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>(h_tdc_UnCorr_bar[thread]);
    hist_map_vec[thread]["h_tdc_TWCorr_bar"] =
        reinterpret_cast<TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>(h_tdc_TWCorr_bar[thread]);
    hist_map_vec[thread]["h_tdc_Corr_bar"] =
        reinterpret_cast<TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>(h_tdc_Corr_bar[thread]);
    hist_map_vec[thread]["h_tdc_avg_bar"] =
        reinterpret_cast<TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>(h_tdc_avg_bar[thread]);
    hist_map_vec[thread]["h_edep_avg_bar"] =
        reinterpret_cast<TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>(h_edep_avg_bar[thread]);
    hist_map_vec[thread]["h_tdc_avg_vs_edep_bar"] =
        reinterpret_cast<TH1 *(*)[N_PLANES][N_SIDES][N_PADDLES]>(h_tdc_avg_vs_edep_bar[thread]);
  }

  for (int thread = 0; thread < numThreads; ++thread) {
    hist_map_trk_cuts_vec[thread]["h_ladkin_time_bar"] =
        reinterpret_cast<TH1 *(*)[nTrackCuts][N_PLANES][N_PADDLES]>(h_ladkin_time_bar[thread]);
    hist_map_trk_cuts_vec[thread]["h_ladkin_edep_bar"] =
        reinterpret_cast<TH1 *(*)[nTrackCuts][N_PLANES][N_PADDLES]>(h_ladkin_edep_bar[thread]);
    hist_map_trk_cuts_vec[thread]["h_ladkin_theta_bar"] =
        reinterpret_cast<TH1 *(*)[nTrackCuts][N_PLANES][N_PADDLES]>(h_ladkin_theta_bar[thread]);
    hist_map_trk_cuts_vec[thread]["h_ladkin_phi_bar"] =
        reinterpret_cast<TH1 *(*)[nTrackCuts][N_PLANES][N_PADDLES]>(h_ladkin_phi_bar[thread]);
    hist_map_trk_cuts_vec[thread]["h_ladkin_deltaposlong_bar"] =
        reinterpret_cast<TH1 *(*)[nTrackCuts][N_PLANES][N_PADDLES]>(h_ladkin_deltaposlong_bar[thread]);
    hist_map_trk_cuts_vec[thread]["h_ladkin_deltapostrans_bar"] =
        reinterpret_cast<TH1 *(*)[nTrackCuts][N_PLANES][N_PADDLES]>(h_ladkin_deltapostrans_bar[thread]);
    hist_map_trk_cuts_vec[thread]["h_ladkin_d0_bar"] =
        reinterpret_cast<TH1 *(*)[nTrackCuts][N_PLANES][N_PADDLES]>(h_ladkin_d0_bar[thread]);
  }

  // Create an output ROOT file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
    return;
  }

  // start threads
  cout << "Starting " << numThreads << " threads..." << endl;
  std::vector<std::thread> threads;
  for (int i_thread = 0; i_thread < numThreads; ++i_thread) {
    int start = i_thread * chunkSize;
    int end   = (i_thread == numThreads - 1) ? nEntries : start + chunkSize;
    threads.emplace_back(process_chunk, i_thread, start, end, fileName, ref(hist_map_vec[i_thread]),
                         ref(hist_map_trk_cuts_vec[i_thread]));
  }
  // Wait for all threads to finish
  for (auto &thread : threads) {
    thread.join();
  }
  std::cout << "\rProcessing: 100% completed. \nMerging histograms" << std::endl;
  // Merge histograms from all threads by adding to the first thread
  for (int i_thread = 0; i_thread < numThreads; ++i_thread) {
    for (const auto &hist_pair : hist_map_vec[i_thread]) {
      for (int plane = 0; plane < N_PLANES; ++plane) {
        for (int side = 0; side < N_SIDES; ++side) {
          for (int bar = 0; bar < N_PADDLES; ++bar) {
            if (i_thread == 0) {
              (hist_map_vec[0][hist_pair.first][0][plane][side][bar])
                  ->SetTitle(TString(hist_map_vec[0][hist_pair.first][0][plane][side][bar]->GetTitle())
                                 .ReplaceAll(Form(" Thread %d", i_thread), ""));
              (hist_map_vec[0][hist_pair.first][0][plane][side][bar])
                  ->SetName(TString(hist_map_vec[0][hist_pair.first][0][plane][side][bar]->GetName())
                                .ReplaceAll(Form("_thread_%d", i_thread), ""));
            } else {

              (hist_map_vec[0][hist_pair.first][0][plane][side][bar])->Add(hist_pair.second[0][plane][side][bar]);
            }
          }
        }
      }
    }
    for (const auto &hist_pair : hist_map_trk_cuts_vec[i_thread]) {
      for (int i_trk_cut = 0; i_trk_cut < nTrackCuts; ++i_trk_cut) {
        for (int plane = 0; plane < N_PLANES; ++plane) {
          for (int bar = 0; bar < N_PADDLES; ++bar) {
            if (i_thread == 0) {
              (hist_map_trk_cuts_vec[0][hist_pair.first][0][i_trk_cut][plane][bar])
                  ->SetTitle(TString(hist_map_trk_cuts_vec[0][hist_pair.first][0][i_trk_cut][plane][bar]->GetTitle())
                                 .ReplaceAll(Form(" Thread %d", i_thread), ""));
              (hist_map_trk_cuts_vec[0][hist_pair.first][0][i_trk_cut][plane][bar])
                  ->SetName(TString(hist_map_trk_cuts_vec[0][hist_pair.first][0][i_trk_cut][plane][bar]->GetName())
                                .ReplaceAll(Form("_thread_%d", i_thread), ""));
            } else {
              (hist_map_trk_cuts_vec[0][hist_pair.first][0][i_trk_cut][plane][bar])
                  ->Add(hist_pair.second[0][i_trk_cut][plane][bar]);
            }
          }
        }
      }
    }
  }

  cout << "Writing histograms to file..." << endl;
  // Write histograms to the output file
  outputFile->cd();
  write_to_canvas(hist_map_vec[0]["h_tdc_bar"][0], outputFile, "Raw/TDC_Time", "TDC");
  write_to_canvas(hist_map_vec[0]["h_adc_bar"][0], outputFile, "Raw/ADC_Time", "ADC");
  write_to_canvas(hist_map_vec[0]["h_adc_amp_bar"][0], outputFile, "Raw/Amp", "ADC_Amp");
  write_to_canvas(hist_map_vec[0]["h_adc_int_bar"][0], outputFile, "Raw/Int", "ADC_Int");
  write_to_canvas(hist_map_vec[0]["h_adc_good_bar"][0], outputFile, "Good/ADC_Time", "ADC_Good");
  write_to_canvas(hist_map_vec[0]["h_adc_amp_good_bar"][0], outputFile, "Good/Amp", "ADC_Amp_Good");
  write_to_canvas(hist_map_vec[0]["h_adc_int_good_bar"][0], outputFile, "Good/Int", "ADC_Int_Good");
  write_to_canvas(hist_map_vec[0]["h_tdc_UnCorr_bar"][0], outputFile, "Good/TDC_Time/UnCorr", "TDC_UnCorr");
  write_to_canvas(hist_map_vec[0]["h_tdc_TWCorr_bar"][0], outputFile, "Good/TDC_Time/TWCorr", "TDC_TWCorr");
  write_to_canvas(hist_map_vec[0]["h_tdc_Corr_bar"][0], outputFile, "Good/TDC_Time/Corr", "TDC_Corr");
  write_to_canvas(hist_map_vec[0]["h_tdc_avg_bar"][0], outputFile, "Good/TDC_Time/Average", "TDC_Avg", true);
  write_to_canvas(hist_map_vec[0]["h_edep_avg_bar"][0], outputFile, "Good/edep/", "EDEP_Avg", true);
  write_to_canvas(hist_map_vec[0]["h_tdc_avg_vs_edep_bar"][0], outputFile, "Good/TDC_vs_Edep", "TDC_vs_Edep", true);

  for (int i_trk_cut = 0; i_trk_cut < nTrackCuts; ++i_trk_cut) {
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_time_bar"][0][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_Time/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_Time_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_edep_bar"][0][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_Edep/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_Edep_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_theta_bar"][0][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_Theta/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_Theta_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_phi_bar"][0][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_Phi/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_Phi_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_deltaposlong_bar"][0][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_DeltaPosLong/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_DeltaPosLong_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_deltapostrans_bar"][0][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_DeltaPosTrans/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_DeltaPosTrans_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_d0_bar"][0][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_D0/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_D0_%s", track_cuts[i_trk_cut].cut_name.c_str()));
  }

  // outputFile->mkdir("TIME_COMP");
  // outputFile->cd("TIME_COMP");
  // // Create a canvas to overlay TDC time, ADC time, Average time, and all hittimes
  // for (int plane = 0; plane < N_PLANES; ++plane) {
  //   TCanvas *c_overlay_plane =
  //       new TCanvas(Form("c_overlay_plane_%d", plane), Form("Overlay Plane %d", plane), 1200, 800);
  //   c_overlay_plane->Divide(4, 3); // Divide canvas into subpads

  //   for (int bar = 0; bar < N_PADDLES; ++bar) {
  //     c_overlay_plane->cd(bar + 1);

  //     // Adjust the y-axis range to fit all histograms
  //     double max_y = 0;
  //     max_y        = std::max(max_y, h_tdc_bar[plane][0][bar]->GetMaximum());
  //     max_y        = std::max(max_y, h_tdc_bar[plane][1][bar]->GetMaximum());
  //     max_y        = std::max(max_y, h_adc_bar[plane][0][bar]->GetMaximum());
  //     max_y        = std::max(max_y, h_adc_bar[plane][1][bar]->GetMaximum());
  //     max_y        = std::max(max_y, h_tdc_avg_bar[plane][bar]->GetMaximum());
  //     for (int i_trk_cut = 0; i_trk_cut < nTrackCuts; ++i_trk_cut) {
  //       max_y = std::max(max_y, h_ladkin_time_bar[i_trk_cut][plane][bar]->GetMaximum());
  //     }
  //     max_y *= 1.1; // Add some padding to the y-axis

  //     h_tdc_bar[plane][0][bar]->SetMaximum(max_y);

  //     // Draw TDC time for both sides
  //     h_tdc_bar[plane][0][bar]->SetLineColor(kRed);
  //     h_tdc_bar[plane][0][bar]->SetLineWidth(2);
  //     h_tdc_bar[plane][0][bar]->Draw("HIST");

  //     h_tdc_bar[plane][1][bar]->SetLineColor(kRed + 2);
  //     h_tdc_bar[plane][1][bar]->SetLineWidth(2);
  //     h_tdc_bar[plane][1][bar]->Draw("HIST SAME");

  //     // Draw ADC time for both sides
  //     h_adc_bar[plane][0][bar]->SetLineColor(kBlue);
  //     h_adc_bar[plane][0][bar]->SetLineWidth(2);
  //     h_adc_bar[plane][0][bar]->Draw("HIST SAME");

  //     h_adc_bar[plane][1][bar]->SetLineColor(kBlue + 2);
  //     h_adc_bar[plane][1][bar]->SetLineWidth(2);
  //     h_adc_bar[plane][1][bar]->Draw("HIST SAME");

  //     // Draw Average time
  //     h_tdc_avg_bar[plane][bar]->SetLineColor(kGreen);
  //     h_tdc_avg_bar[plane][bar]->SetLineWidth(2);
  //     h_tdc_avg_bar[plane][bar]->Draw("HIST SAME");

  //     // Draw hittimes for all track cuts
  //     for (int i_trk_cut = 0; i_trk_cut < nTrackCuts; ++i_trk_cut) {
  //       h_ladkin_time_bar[i_trk_cut][plane][bar]->SetLineColor(kMagenta + i_trk_cut);
  //       h_ladkin_time_bar[i_trk_cut][plane][bar]->SetLineWidth(2);
  //       h_ladkin_time_bar[i_trk_cut][plane][bar]->Draw("HIST SAME");
  //     }
  //   }

  //   // Add a legend only once for the entire plane
  //   TLegend *legend = new TLegend(0.7, 0.6, 0.9, 0.9);
  //   legend->AddEntry(h_tdc_bar[plane][0][0], "TDC Time Top", "l");
  //   legend->AddEntry(h_tdc_bar[plane][1][0], "TDC Time Bottom", "l");
  //   legend->AddEntry(h_adc_bar[plane][0][0], "ADC Time Top", "l");
  //   legend->AddEntry(h_adc_bar[plane][1][0], "ADC Time Bottom", "l");
  //   legend->AddEntry(h_tdc_avg_bar[plane][0], "Average Time", "l");
  //   for (int i_trk_cut = 0; i_trk_cut < nTrackCuts; ++i_trk_cut) {
  //     legend->AddEntry(h_ladkin_time_bar[i_trk_cut][plane][0],
  //                      Form("Hittime %s", track_cuts[i_trk_cut].cut_name.c_str()), "l");
  //   }
  //   c_overlay_plane->cd(12); // Use the last pad for the legend
  //   legend->Draw();

  //   // Write the canvas to the output file
  //   c_overlay_plane->Write();
  //   delete c_overlay_plane;

  //   // Create a canvas for the entire plane as a whole
  //   // TCanvas *c_plane_whole = new TCanvas(Form("c_plane_whole_%d", plane), Form("Plane %d Whole", plane), 800,
  //   600);

  //   // // Adjust the y-axis range to fit all histograms
  //   // double max_y_plane = 0;
  //   // max_y_plane        = std::max(max_y_plane, h_tdc_plane[plane][0]->GetMaximum());
  //   // max_y_plane        = std::max(max_y_plane, h_tdc_plane[plane][1]->GetMaximum());
  //   // max_y_plane        = std::max(max_y_plane, h_adc_plane[plane][0]->GetMaximum());
  //   // max_y_plane        = std::max(max_y_plane, h_adc_plane[plane][1]->GetMaximum());
  //   // max_y_plane        = std::max(max_y_plane, h_tdc_avg_plane[plane]->GetMaximum());
  //   // for (int i_trk_cut = 0; i_trk_cut < nTrackCuts; ++i_trk_cut) {
  //   //   max_y_plane = std::max(max_y_plane, h_ladkin_time[i_trk_cut][plane]->GetMaximum());
  //   // }
  //   // max_y_plane *= 1.1; // Add some padding to the y-axis

  //   // h_tdc_plane[plane][0]->SetMaximum(max_y_plane);

  //   // Draw TDC time for both sides
  //   // h_tdc_plane[plane][0]->SetLineColor(kRed);
  //   // h_tdc_plane[plane][0]->SetLineWidth(2);
  //   // h_tdc_plane[plane][0]->Draw("HIST");

  //   // h_tdc_plane[plane][1]->SetLineColor(kRed + 2);
  //   // h_tdc_plane[plane][1]->SetLineWidth(2);
  //   // h_tdc_plane[plane][1]->Draw("HIST SAME");

  //   // // Draw ADC time for both sides
  //   // h_adc_plane[plane][0]->SetLineColor(kBlue);
  //   // h_adc_plane[plane][0]->SetLineWidth(2);
  //   // h_adc_plane[plane][0]->Draw("HIST SAME");

  //   // h_adc_plane[plane][1]->SetLineColor(kBlue + 2);
  //   // h_adc_plane[plane][1]->SetLineWidth(2);
  //   // h_adc_plane[plane][1]->Draw("HIST SAME");

  //   // Draw Average time
  //   // h_tdc_avg_plane[plane]->SetMaximum(max_y_plane);
  //   // h_tdc_avg_plane[plane]->SetLineColor(kGreen);
  //   // h_tdc_avg_plane[plane]->SetLineWidth(2);
  //   // h_tdc_avg_plane[plane]->Draw("HIST");

  //   // // Draw hittimes for all track cuts
  //   // for (int i_trk_cut = 0; i_trk_cut < nTrackCuts; ++i_trk_cut) {
  //   //   h_ladkin_time[i_trk_cut][plane]->SetLineColor(kMagenta + i_trk_cut);
  //   //   h_ladkin_time[i_trk_cut][plane]->SetLineWidth(2);
  //   //   h_ladkin_time[i_trk_cut][plane]->Draw("HIST SAME");
  //   // }

  //   // // Add a legend for the whole plane
  //   // TLegend *legend_plane = new TLegend(0.7, 0.6, 0.9, 0.9);
  //   // legend_plane->AddEntry(h_tdc_plane[plane][0], "TDC Time Top", "l");
  //   // legend_plane->AddEntry(h_tdc_plane[plane][1], "TDC Time Bottom", "l");
  //   // legend_plane->AddEntry(h_adc_plane[plane][0], "ADC Time Top", "l");
  //   // legend_plane->AddEntry(h_adc_plane[plane][1], "ADC Time Bottom", "l");
  //   // legend_plane->AddEntry(h_tdc_avg_plane[plane], "Average Time", "l");
  //   // for (int i_trk_cut = 0; i_trk_cut < nTrackCuts; ++i_trk_cut) {
  //   //   legend_plane->AddEntry(h_ladkin_time[i_trk_cut][plane],
  //   //                          Form("Hittime %s", track_cuts[i_trk_cut].cut_name.c_str()), "l");
  //   // }
  //   // legend_plane->Draw();

  //   // // Write the canvas to the output file
  //   // c_plane_whole->Write();
  //   // delete c_plane_whole;
  // }
  outputFile->Close();
  delete outputFile;

  std::cout << "Done! All histograms written to file. Cleaning up..." << std::endl;
  // Clean up
  numThreads = std::thread::hardware_concurrency();
  std::vector<std::thread> cleanup_threads;
  for (int thread = 0; thread < numThreads; ++thread) {
    cleanup_threads.emplace_back([&, thread]() {
      for (int plane = 0; plane < N_PLANES; ++plane) {
        for (int bar = 0; bar < N_PADDLES; ++bar) {
          for (int side = 0; side < N_SIDES; ++side) {
            delete h_tdc_bar[thread][plane][side][bar];
            delete h_adc_bar[thread][plane][side][bar];
            delete h_adc_good_bar[thread][plane][side][bar];
            delete h_adc_amp_bar[thread][plane][side][bar];
            delete h_adc_amp_good_bar[thread][plane][side][bar];
            delete h_adc_int_bar[thread][plane][side][bar];
            delete h_adc_int_good_bar[thread][plane][side][bar];
            delete h_tdc_UnCorr_bar[thread][plane][side][bar];
            delete h_tdc_TWCorr_bar[thread][plane][side][bar];
            delete h_tdc_Corr_bar[thread][plane][side][bar];

            // Initializing side 0 and 1 for following histograms to prevent segfaults. Only side 0 is used (data values
            // are avg of top and btm)
            if (side == 1) {
              delete h_tdc_avg_bar[thread][plane][side][bar];
              delete h_edep_avg_bar[thread][plane][side][bar];
              delete h_tdc_avg_vs_edep_bar[thread][plane][side][bar];
            }
          }
          for (int i_trk_cut = 0; i_trk_cut < nTrackCuts; ++i_trk_cut) {
            delete h_ladkin_time_bar[thread][i_trk_cut][plane][bar];
            delete h_ladkin_edep_bar[thread][i_trk_cut][plane][bar];
            delete h_ladkin_theta_bar[thread][i_trk_cut][plane][bar];
            delete h_ladkin_phi_bar[thread][i_trk_cut][plane][bar];
            delete h_ladkin_deltaposlong_bar[thread][i_trk_cut][plane][bar];
            delete h_ladkin_deltapostrans_bar[thread][i_trk_cut][plane][bar];
            delete h_ladkin_d0_bar[thread][i_trk_cut][plane][bar];
          }
        }
      }
    });
  }

  // Wait for all cleanup threads to finish
  for (auto &cleanup_thread : cleanup_threads) {
    cleanup_thread.join();
  }
   

  cout << "Done!" << endl;
  // Close the output file
  

  return;
}
