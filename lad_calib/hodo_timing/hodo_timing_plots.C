// Lucas Ehinger
// General LAD hodo plotting script
// Multithreaded, to speed up the process
// It's absolutely terrible when it comes to memory management, but it works
// Multi-threading with ROOT which is inherently not thread-safe is hard, and I wasn't smart enough to do it elegantly,
// but this works.
#include </usr/lib/gcc/x86_64-redhat-linux/11/include/omp.h>
#include <TCanvas.h>
#include <TChain.h>
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

using namespace std;

const int MAX_DATA     = 500;
const int MAX_DATA_GEM = 3000;

struct hist_params {
  int NBINS;
  double MIN;
  double MAX;
};

const hist_params d0_params           = {100, 0.0, 100.0};
const hist_params projz_params        = {100, -20.0, 20.0};
const hist_params gem_X_params        = {100, -70.0, 70.0};
const hist_params gem_Y_params        = {100, -40.0, 40.0};
const hist_params dTrk_params         = {100, -500.0, 500.0};
const hist_params time_params         = {100, 1600.0, 2000.0};
const hist_params tof_params          = {100, -300, 500};
const hist_params raw_tdc_time_params = {100, 0, 4000};
const hist_params raw_adc_time_params = {100, 0, 500};
const hist_params theta_params        = {180, 0.0, 180.0};
const hist_params phi_params          = {360, -180.0, 180.0};
const hist_params paddle_params       = {100, 0.0, 100.0};
const hist_params hit_position_params = {100, -500.0, 500.0};
const hist_params edep_params         = {1000, 0.0, 100.0};
const hist_params amp_params          = {1000, 0.0, 1000.0};
const hist_params int_params          = {1000, 0.0, 2.0};
const double TDC2NS                   = 0.09766;
;                             // TDC to ns conversion factor
const double ADC2NS = 0.0625; // ADC to ns conversion factor
// const double ADC2PC = 0.1;   // ADC to percent conversion factor
// const double ADC2mV = 0.1;  // ADC to mV conversion factor

const int MINT_EVTS_PER_THREAD = 20000;
const  char spec_prefix = 'H'; // Default value
struct track_cut {
  double dTrkVert_cut;
  double dTrkHoriz_cut;
  double d0_cut;
  string cut_name;
};

const int nTrackCuts             = 5;
track_cut track_cuts[nTrackCuts] = {{2000, 2000.0, 200.0, "loose"},
                                    {100.0, 2000.0, 200.0, "trans_50"},
                                    {100.0, 100.0, 100.0, "trans_50_long_100"},
                                    {100.0, 100.0, 25.0, "trans_50_long_100_d0_25"},
                                    {100.0, 100.0, 10.0, "trans_50_long_100_d0_10"}};
// const int nTrackCuts             = 1;
// track_cut track_cuts[nTrackCuts] = {{2000, 2000.0, 200.0, "loose"}};
// const double dTrkVert_cut  = 200.0;
// const double dTrkHoriz_cut = 800.0;
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
      TString side_label = single_side ? "Top&Btm" : side_names[side].c_str();
      TCanvas *c1 =
          new TCanvas(Form("c_%s_plane_%s_%s", var_name.Data(), plane_names[plane].c_str(), side_label.Data()),
                      Form("%s Plane %s %s", var_name.Data(), plane_names[plane].c_str(), side_label.Data()), 800, 600);

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
            Form("h_%s_plane_%s_%s_all_bars", var_name.Data(), plane_names[plane].c_str(), side_label.Data()));
        hist_plane_sum[plane][side]->Reset();
        hist_plane_sum[plane][side]->SetTitle(
            Form("%s Plane %s %s All Bars", var_name.Data(), plane_names[plane].c_str(), side_label.Data()));
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
    TCanvas *c_all_sides =
        new TCanvas(Form("c_%s_plane_%s_all_sides", var_name.Data(), plane_names[plane].c_str()),
                    Form("%s Plane %s All Sides", var_name.Data(), plane_names[plane].c_str()), 1200, 800);
    if (!single_side)
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
    TCanvas *c1 = new TCanvas(Form("c_%s_plane_%s", var_name.Data(), plane_names[plane].c_str()),
                              Form("%s Plane %s", var_name.Data(), plane_names[plane].c_str()), 800, 600);
    c1->Divide(4, 3); // Divide canvas into subpads
    for (int bar = 0; bar < N_PADDLES; ++bar) {
      c1->cd(bar + 1);
      hist_arr[plane][bar]->Draw();
    }
    c1->Write();
    delete c1;

    // Sum all bars for this plane into a single histogram
    if (!hist_plane_sum[plane]) {
      hist_plane_sum[plane] = (TH1F *)hist_arr[plane][0]->Clone(
          Form("h_%s_plane_%s_all_bars", var_name.Data(), plane_names[plane].c_str()));
      hist_plane_sum[plane]->Reset();
      hist_plane_sum[plane]->SetTitle(Form("%s Plane %s All Bars", var_name.Data(), plane_names[plane].c_str()));
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

void process_chunk(int i_thread, int start, int end, std::vector<TString> &fileNames,
                   map<string, TH1 *[N_PLANES][N_SIDES][N_PADDLES]> &hist_map,
                   map<string, TH1 *[nTrackCuts][N_PLANES][N_PADDLES]> &hist_map_trk_cuts) {

  TChain *T = new TChain("T");
  for (const auto &fileName : fileNames) {
    T->Add(fileName);
  }

  if (T->GetEntries() == 0) {
    std::cerr << "Error: Cannot open the ROOT files or no entries in the TChain!" << std::endl;
    return;
  }

  // Define arrays to hold the data
  Double_t trk_d0[MAX_DATA_GEM];
  Double_t trk_d0_good[MAX_DATA_GEM];
  Double_t trk_projz[MAX_DATA_GEM];
  Double_t trk_xloc[2][MAX_DATA_GEM], trk_yloc[2][MAX_DATA_GEM];
  Double_t kin_trackID_0[MAX_DATA_GEM], kin_trackID_1[MAX_DATA_GEM];
  Double_t kin_plane_0[MAX_DATA_GEM], kin_plane_1[MAX_DATA_GEM];
  Double_t kin_paddle_0[MAX_DATA_GEM], kin_paddle_1[MAX_DATA_GEM];
  Double_t kin_hittime_0[MAX_DATA_GEM], kin_hittime_1[MAX_DATA_GEM];
  Double_t kin_hit_tof_0[MAX_DATA_GEM], kin_hit_tof_1[MAX_DATA_GEM];
  Double_t kin_hittheta_0[MAX_DATA_GEM], kin_hittheta_1[MAX_DATA_GEM];
  Double_t kin_hitphi_0[MAX_DATA_GEM], kin_hitphi_1[MAX_DATA_GEM];
  Double_t kin_hitedep_0[MAX_DATA_GEM], kin_hitedep_1[MAX_DATA_GEM];
  Double_t kin_dTrkHoriz_0[MAX_DATA_GEM], kin_dTrkHoriz_1[MAX_DATA_GEM];
  Double_t kin_dTrkVert_0[MAX_DATA_GEM], kin_dTrkVert_1[MAX_DATA_GEM];
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
  Double_t hodo_start_time, pTRIG1, pTRIG2, pTRIG3, pTRIG4, evtyp;

  T->SetBranchAddress(Form("Ndata.%c.gem.trk.d0", spec_prefix), &nTracks);
  T->SetBranchAddress(Form("Ndata.%c.ladkin.goodhit_trackid_0", spec_prefix), &nGoodHits);
  T->SetBranchAddress(Form("%c.gem.trk.d0", spec_prefix), &trk_d0);
  T->SetBranchAddress(Form("%c.gem.trk.d0_good", spec_prefix), &trk_d0_good);
  T->SetBranchAddress(Form("%c.gem.trk.projz", spec_prefix), &trk_projz);
  T->SetBranchAddress(Form("%c.gem.trk.x1_local", spec_prefix), &trk_xloc[0]);
  T->SetBranchAddress(Form("%c.gem.trk.x2_local", spec_prefix), &trk_xloc[1]);
  T->SetBranchAddress(Form("%c.gem.trk.y1_local", spec_prefix), &trk_yloc[0]);
  T->SetBranchAddress(Form("%c.gem.trk.y2_local", spec_prefix), &trk_yloc[1]);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_trackid_0", spec_prefix), &kin_trackID_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_trackid_1", spec_prefix), &kin_trackID_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_plane_0", spec_prefix), &kin_plane_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_plane_1", spec_prefix), &kin_plane_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_paddle_0", spec_prefix), &kin_paddle_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_paddle_1", spec_prefix), &kin_paddle_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hittime_0", spec_prefix), &kin_hittime_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hittime_1", spec_prefix), &kin_hittime_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_tof_0", spec_prefix), &kin_hit_tof_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_tof_1", spec_prefix), &kin_hit_tof_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hittheta_0", spec_prefix), &kin_hittheta_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hittheta_1", spec_prefix), &kin_hittheta_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hitphi_0", spec_prefix), &kin_hitphi_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hitphi_1", spec_prefix), &kin_hitphi_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hitedep_0", spec_prefix), &kin_hitedep_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hitedep_1", spec_prefix), &kin_hitedep_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_dTrkHoriz_0", spec_prefix), &kin_dTrkHoriz_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_dTrkHoriz_1", spec_prefix), &kin_dTrkHoriz_1);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_dTrkVert_0", spec_prefix), &kin_dTrkVert_0);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_dTrkVert_1", spec_prefix), &kin_dTrkVert_1);

  T->SetBranchAddress(Form("%c.hod.starttime", spec_prefix), &hodo_start_time);
  T->SetBranchAddress("T.hms.pTRIG1_tdcTime", &pTRIG1);
  T->SetBranchAddress("T.hms.pTRIG2_tdcTime", &pTRIG2);
  T->SetBranchAddress("T.hms.pTRIG3_tdcTime", &pTRIG3);
  T->SetBranchAddress("T.hms.pTRIG4_tdcTime", &pTRIG4);
  T->SetBranchAddress("g.evtyp", &evtyp);

  for (int plane = 0; plane < N_PLANES; ++plane) {
    for (int side = 0; side < N_SIDES; ++side) {
      T->SetBranchAddress(
          Form("%c.ladhod.%s.%sTdcTimeRaw", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &tdc_time[plane][side]);
      T->SetBranchAddress(
          Form("%c.ladhod.%s.Good%sTdcTimeUnCorr", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &tdc_time_UnCorr[plane][side]);
      T->SetBranchAddress(
          Form("%c.ladhod.%s.Good%sTdcTimeWalkCorr", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &tdc_time_TWCorr[plane][side]);
      T->SetBranchAddress(
          Form("%c.ladhod.%s.Good%sTdcTimeCorr", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &tdc_time_Corr[plane][side]);
      T->SetBranchAddress(
          Form("%c.ladhod.%s.%sAdcPulseTimeRaw", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &adc_time[plane][side]);
      T->SetBranchAddress(
          Form("%c.ladhod.%s.Good%sAdcPulseTime", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &adc_time_good[plane][side]);
      T->SetBranchAddress(
          Form("%c.ladhod.%s.%sAdcPulseAmp", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &adc_amp[plane][side]);
      T->SetBranchAddress(
          Form("%c.ladhod.%s.Good%sAdcPulseAmp", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &adc_amp_good[plane][side]);
      T->SetBranchAddress(
          Form("%c.ladhod.%s.%sAdcPulseInt", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &adc_int[plane][side]);
      T->SetBranchAddress(
          Form("%c.ladhod.%s.Good%sAdcPulseInt", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &adc_int_good[plane][side]);
      T->SetBranchAddress(
          Form("%c.ladhod.%s.%sTdcCounter", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &tdc_counter[plane][side]);
      T->SetBranchAddress(
          Form("%c.ladhod.%s.%sAdcCounter", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &adc_counter[plane][side]);
      T->SetBranchAddress(
          Form("Ndata.%c.ladhod.%s.%sAdcCounter", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &nData_adc[plane][side]);
      T->SetBranchAddress(
          Form("Ndata.%c.ladhod.%s.%sTdcCounter", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &nData_tdc[plane][side]);
    }
    T->SetBranchAddress(Form("%c.ladhod.%s.HodoHitTime", spec_prefix, plane_names[plane].c_str()),
                        &fullhit_time_avg[plane]);
    T->SetBranchAddress(Form("%c.ladhod.%s.HodoHitEdep", spec_prefix, plane_names[plane].c_str()),
                        &fullhit_adc_avg[plane]);
    T->SetBranchAddress(Form("%c.ladhod.%s.HodoHitPaddleNum", spec_prefix, plane_names[plane].c_str()),
                        &fullhit_paddle[plane]);
    T->SetBranchAddress(Form("Ndata.%c.ladhod.%s.HodoHitTime", spec_prefix, plane_names[plane].c_str()),
                        &fullhit_n[plane]);
  }
  ////////////////////////////////////////////////////
  // Start Event Loop

  for (int i = start; i < end; ++i) {
    T->GetEntry(i);
    if ((spec_prefix == 'H' && int(evtyp) != 2) || (spec_prefix == 'P' && int(evtyp) != 1))
      continue; // Only process events with the correct evtyp based on spec_prefix

    for (int plane = 0; plane < N_PLANES; ++plane) {
      for (int side = 0; side < N_SIDES; ++side) {

        // Fill Raw Values
        for (int j = 0; j < nData_tdc[plane][side]; ++j) {
          hist_map["h_tdc_bar"][plane][side][int(tdc_counter[plane][side][j] - 1)]->Fill(tdc_time[plane][side][j] *
                                                                                         TDC2NS);
        }
        for (int j = 0; j < nData_adc[plane][side]; ++j) {
          hist_map["h_adc_bar"][plane][side][int(adc_counter[plane][side][j] - 1)]->Fill(adc_time[plane][side][j] *
                                                                                         ADC2NS);
          hist_map["h_adc_amp_bar"][plane][side][int(adc_counter[plane][side][j] - 1)]->Fill(adc_amp[plane][side][j]);
          hist_map["h_adc_int_bar"][plane][side][int(adc_counter[plane][side][j] - 1)]->Fill(adc_int[plane][side][j]);
        }

        // Fill Good Values
        for (int j = 0; j < N_PADDLES; ++j) {

          hist_map["h_adc_good_bar"][plane][side][j]->Fill(adc_time_good[plane][side][j]);
          hist_map["h_adc_amp_good_bar"][plane][side][j]->Fill(adc_amp_good[plane][side][j]);
          hist_map["h_adc_int_bar"][plane][side][j]->Fill(adc_int[plane][side][j]);
          hist_map["h_adc_int_good_bar"][plane][side][j]->Fill(adc_int_good[plane][side][j]);
          hist_map["h_tdc_UnCorr_bar"][plane][side][j]->Fill(tdc_time_UnCorr[plane][side][j]);
          hist_map["h_tdc_TWCorr_bar"][plane][side][j]->Fill(tdc_time_TWCorr[plane][side][j]);
          hist_map["h_tdc_Corr_bar"][plane][side][j]->Fill(tdc_time_Corr[plane][side][j]);
        }
      }
      for (int j = 0; j < fullhit_n[plane]; j++) {
        hist_map["h_tdc_avg_bar"][plane][0][int(fullhit_paddle[plane][j] - 1)]->Fill(fullhit_time_avg[plane][j]);
        hist_map["h_edep_avg_bar"][plane][0][int(fullhit_paddle[plane][j] - 1)]->Fill(fullhit_adc_avg[plane][j]);
        hist_map["h_tdc_avg_vs_edep_bar"][plane][0][int(fullhit_paddle[plane][j] - 1)]->Fill(fullhit_time_avg[plane][j],
                                                                                             fullhit_adc_avg[plane][j]);
      }
    } // End Plane Loop

    for (int i_goodhit = 0; i_goodhit < nGoodHits; i_goodhit++) {
      for (int i_track_cut = 0; i_track_cut < nTrackCuts; i_track_cut++) {
        if (abs(kin_dTrkVert_0[i_goodhit]) < track_cuts[i_track_cut].dTrkVert_cut &&
            abs(kin_dTrkHoriz_0[i_goodhit]) < track_cuts[i_track_cut].dTrkHoriz_cut &&
            trk_d0[int(kin_trackID_0[i_goodhit])] < track_cuts[i_track_cut].d0_cut) {

          hist_map_trk_cuts["h_ladkin_time_bar"][i_track_cut][int(kin_plane_0[i_goodhit])][int(kin_paddle_0[i_goodhit])]
              ->Fill(kin_hittime_0[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_tof_bar"][i_track_cut][int(kin_plane_0[i_goodhit])][int(kin_paddle_0[i_goodhit])]
            ->Fill(kin_hit_tof_0[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_edep_bar"][i_track_cut][int(kin_plane_0[i_goodhit])][int(kin_paddle_0[i_goodhit])]
              ->Fill(kin_hitedep_0[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_theta_bar"][i_track_cut][int(kin_plane_0[i_goodhit])]
                           [int(kin_paddle_0[i_goodhit])]
                               ->Fill(kin_hittheta_0[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_phi_bar"][i_track_cut][int(kin_plane_0[i_goodhit])][int(kin_paddle_0[i_goodhit])]
              ->Fill(kin_hitphi_0[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_dTrkVert_bar"][i_track_cut][int(kin_plane_0[i_goodhit])]
                           [int(kin_paddle_0[i_goodhit])]
                               ->Fill(kin_dTrkVert_0[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_dTrkHoriz_bar"][i_track_cut][int(kin_plane_0[i_goodhit])]
                           [int(kin_paddle_0[i_goodhit])]
                               ->Fill(kin_dTrkHoriz_0[i_goodhit]);
          hist_map_trk_cuts["h_edep_vs_time_bar"][i_track_cut][int(kin_plane_0[i_goodhit])]
                           [int(kin_paddle_0[i_goodhit])]
                               ->Fill(kin_hittime_0[i_goodhit], kin_hitedep_0[i_goodhit]);
          hist_map_trk_cuts["h_edep_vs_tof_bar"][i_track_cut][int(kin_plane_0[i_goodhit])][int(kin_paddle_0[i_goodhit])]
              ->Fill(kin_hit_tof_0[i_goodhit], kin_hitedep_0[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_d0_bar"][i_track_cut][int(kin_plane_0[i_goodhit])][int(kin_paddle_0[i_goodhit])]
              ->Fill(trk_d0[int(kin_trackID_0[i_goodhit])]);
          hist_map_trk_cuts["h_projz_bar"][i_track_cut][int(kin_plane_0[i_goodhit])][int(kin_paddle_0[i_goodhit])]
              ->Fill(trk_projz[int(kin_trackID_0[i_goodhit])]);
          hist_map_trk_cuts["h_XY_GEM0_bar"][i_track_cut][int(kin_plane_0[i_goodhit])][int(kin_paddle_0[i_goodhit])]
              ->Fill(trk_xloc[0][int(kin_trackID_0[i_goodhit])], trk_yloc[0][int(kin_trackID_0[i_goodhit])]);
          hist_map_trk_cuts["h_XY_GEM1_bar"][i_track_cut][int(kin_plane_0[i_goodhit])][int(kin_paddle_0[i_goodhit])]
              ->Fill(trk_xloc[1][int(kin_trackID_0[i_goodhit])], trk_yloc[1][int(kin_trackID_0[i_goodhit])]);
        }
        if (abs(kin_dTrkVert_1[i_goodhit]) < track_cuts[i_track_cut].dTrkVert_cut &&
            abs(kin_dTrkHoriz_1[i_goodhit]) < track_cuts[i_track_cut].dTrkHoriz_cut &&
            trk_d0[int(kin_trackID_1[i_goodhit])] < track_cuts[i_track_cut].d0_cut) {

          hist_map_trk_cuts["h_ladkin_time_bar"][i_track_cut][int(kin_plane_1[i_goodhit])][int(kin_paddle_1[i_goodhit])]
              ->Fill(kin_hittime_1[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_tof_bar"][i_track_cut][int(kin_plane_1[i_goodhit])][int(kin_paddle_1[i_goodhit])]
            ->Fill(kin_hit_tof_1[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_edep_bar"][i_track_cut][int(kin_plane_1[i_goodhit])][int(kin_paddle_1[i_goodhit])]
              ->Fill(kin_hitedep_1[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_theta_bar"][i_track_cut][int(kin_plane_1[i_goodhit])]
                           [int(kin_paddle_1[i_goodhit])]
                               ->Fill(kin_hittheta_1[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_phi_bar"][i_track_cut][int(kin_plane_1[i_goodhit])][int(kin_paddle_1[i_goodhit])]
              ->Fill(kin_hitphi_1[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_dTrkVert_bar"][i_track_cut][int(kin_plane_1[i_goodhit])]
                           [int(kin_paddle_1[i_goodhit])]
                               ->Fill(kin_dTrkVert_1[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_dTrkHoriz_bar"][i_track_cut][int(kin_plane_1[i_goodhit])]
                           [int(kin_paddle_1[i_goodhit])]
                               ->Fill(kin_dTrkHoriz_1[i_goodhit]);
          hist_map_trk_cuts["h_edep_vs_time_bar"][i_track_cut][int(kin_plane_1[i_goodhit])]
                           [int(kin_paddle_1[i_goodhit])]
                               ->Fill(kin_hittime_1[i_goodhit], kin_hitedep_1[i_goodhit]);
          hist_map_trk_cuts["h_edep_vs_tof_bar"][i_track_cut][int(kin_plane_1[i_goodhit])]
                           [int(kin_paddle_1[i_goodhit])]
                               ->Fill(kin_hit_tof_1[i_goodhit], kin_hitedep_1[i_goodhit]);
          hist_map_trk_cuts["h_ladkin_d0_bar"][i_track_cut][int(kin_plane_1[i_goodhit])][int(kin_paddle_1[i_goodhit])]
              ->Fill(trk_d0[int(kin_trackID_1[i_goodhit])]);
          hist_map_trk_cuts["h_projz_bar"][i_track_cut][int(kin_plane_1[i_goodhit])][int(kin_paddle_1[i_goodhit])]
              ->Fill(trk_projz[int(kin_trackID_1[i_goodhit])]);
          hist_map_trk_cuts["h_XY_GEM0_bar"][i_track_cut][int(kin_plane_1[i_goodhit])][int(kin_paddle_1[i_goodhit])]
              ->Fill(trk_xloc[0][int(kin_trackID_1[i_goodhit])], trk_yloc[0][int(kin_trackID_1[i_goodhit])]);
          hist_map_trk_cuts["h_XY_GEM1_bar"][i_track_cut][int(kin_plane_1[i_goodhit])][int(kin_paddle_1[i_goodhit])]
              ->Fill(trk_xloc[1][int(kin_trackID_1[i_goodhit])], trk_yloc[1][int(kin_trackID_1[i_goodhit])]);
        }
      }
    }

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

  // Open multiple ROOT files
  std::vector<TString> fileNames = {
      // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22591_0_2_-1.root"};

  // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22572_0_6_2000000.root"};
  // "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/first_seg_runs/LAD_COIN_22615_6_6_-1.root"  
  // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22382_0_21_-1.root",
  // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22382_0_21_-1_1.root"};

  // "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22562_0_0_-1.root",
  // "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22563_0_0_-1.root",
  // "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22564_0_0_-1.root", 
  // "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22565_0_0_-1.root",
  // "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22566_0_0_-1.root",
  // "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22567_0_0_-1.root", 
  // "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22568_0_0_-1.root", 

  "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22611_0_6_-1.root", 


};

  TString outputFileName = Form("files/hodo_timing_plots/hodo_timing_plots_22611_%c.root", spec_prefix);

  // Create a TChain to combine the trees from multiple files
  TChain *chain = new TChain("T");
  for (const auto &fileName : fileNames) {
    chain->Add(fileName);
  }

  if (chain->GetEntries() == 0) {
    std::cerr << "Error: Cannot open the ROOT files or no entries in the TChain!" << std::endl;
    delete chain;
    return;
  }

  // Get the TTree
  TTree *T = chain;
  if (!T) {
    std::cerr << "Error: Cannot find the TTree named 'T'!" << std::endl;
    delete chain;
    return;
  }

  // Number of entries in the TTree
  int nEntries = T->GetEntries();

  // Number of threads to use
  int numThreads = std::thread::hardware_concurrency();
  // numThreads     = 1;
  int chunkSize = nEntries / numThreads;

  // Adjust the number of threads if the chunk size is too small
  if (chunkSize < MINT_EVTS_PER_THREAD) {
    numThreads = std::max(1, nEntries / MINT_EVTS_PER_THREAD);
    chunkSize  = nEntries / numThreads;
  }

  std::vector<map<string, TH1 *[N_PLANES][N_SIDES][N_PADDLES]>> hist_map_vec(numThreads);
  std::vector<map<string, TH1 *[nTrackCuts][N_PLANES][N_PADDLES]>> hist_map_trk_cuts_vec(numThreads);

  for (int thread = 0; thread < numThreads; ++thread) {
    for (int plane = 0; plane < N_PLANES; ++plane) {
      for (int bar = 0; bar < N_PADDLES; ++bar) {
        for (int side = 0; side < N_SIDES; ++side) {
          hist_map_vec[thread]["h_tdc_bar"][plane][side][bar] =
              new TH1F(Form("h_tdc_%s_plane_%s_bar_%d_thread_%d", side_names[side].c_str(), plane_names[plane].c_str(),
                            bar, thread),
                       Form("TDC Raw Time %s Plane %s Bar %d Thread %d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       raw_tdc_time_params.NBINS, raw_tdc_time_params.MIN, raw_tdc_time_params.MAX);
          hist_map_vec[thread]["h_tdc_bar"][plane][side][bar]->GetXaxis()->SetTitle("TDC Time (ns)");
          hist_map_vec[thread]["h_tdc_bar"][plane][side][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_adc_bar"][plane][side][bar] =
              new TH1F(Form("h_adc_%s_plane_%s_bar_%d_thread_%d", side_names[side].c_str(), plane_names[plane].c_str(),
                            bar, thread),
                       Form("ADC Raw Time %s Plane %s Bar %d Thread %d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       raw_adc_time_params.NBINS, raw_adc_time_params.MIN, raw_adc_time_params.MAX);
          hist_map_vec[thread]["h_adc_bar"][plane][side][bar]->GetXaxis()->SetTitle("ADC Time (ns)");
          hist_map_vec[thread]["h_adc_bar"][plane][side][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_adc_good_bar"][plane][side][bar] =
              new TH1F(Form("h_adc_good_%s_plane_%s_bar_%d_thread_%d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       Form("ADC Time Good %s Plane %s Bar %d Thread %d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       time_params.NBINS, time_params.MIN, time_params.MAX);
          hist_map_vec[thread]["h_adc_good_bar"][plane][side][bar]->GetXaxis()->SetTitle("ADC Time (ns)");
          hist_map_vec[thread]["h_adc_good_bar"][plane][side][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_adc_amp_bar"][plane][side][bar] =
              new TH1F(Form("h_adc_amp_%s_plane_%s_bar_%d_thread_%d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       Form("ADC Amp %s Plane %s Bar %d Thread %d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       amp_params.NBINS, amp_params.MIN, amp_params.MAX);
          hist_map_vec[thread]["h_adc_amp_bar"][plane][side][bar]->GetXaxis()->SetTitle("Amplitude (mV)");
          hist_map_vec[thread]["h_adc_amp_bar"][plane][side][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_adc_amp_good_bar"][plane][side][bar] =
              new TH1F(Form("h_adc_amp_good_%s_plane_%s_bar_%d_thread_%d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       Form("ADC Amp Good %s Plane %s Bar %d Thread %d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       amp_params.NBINS, amp_params.MIN, amp_params.MAX);
          hist_map_vec[thread]["h_adc_amp_good_bar"][plane][side][bar]->GetXaxis()->SetTitle("Amplitude (mV)");
          hist_map_vec[thread]["h_adc_amp_good_bar"][plane][side][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_adc_int_bar"][plane][side][bar] =
              new TH1F(Form("h_adc_int_%s_plane_%s_bar_%d_thread_%d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       Form("ADC Int %s Plane %s Bar %d Thread %d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       int_params.NBINS, int_params.MIN, int_params.MAX);
          hist_map_vec[thread]["h_adc_int_bar"][plane][side][bar]->GetXaxis()->SetTitle("Integral (pC)");
          hist_map_vec[thread]["h_adc_int_bar"][plane][side][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_adc_int_good_bar"][plane][side][bar] =
              new TH1F(Form("h_adc_int_good_%s_plane_%s_bar_%d_thread_%d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       Form("ADC Int Good %s Plane %s Bar %d Thread %d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       int_params.NBINS, int_params.MIN, int_params.MAX);
          hist_map_vec[thread]["h_adc_int_good_bar"][plane][side][bar]->GetXaxis()->SetTitle("Integral (pC)");
          hist_map_vec[thread]["h_adc_int_good_bar"][plane][side][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_tdc_UnCorr_bar"][plane][side][bar] =
              new TH1F(Form("h_tdc_UnCorr_%s_plane_%s_bar_%d_thread_%d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       Form("TDC UnCorr %s Plane %s Bar %d Thread %d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       time_params.NBINS, time_params.MIN, time_params.MAX);
          hist_map_vec[thread]["h_tdc_UnCorr_bar"][plane][side][bar]->GetXaxis()->SetTitle("TDC Time (ns)");
          hist_map_vec[thread]["h_tdc_UnCorr_bar"][plane][side][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_tdc_TWCorr_bar"][plane][side][bar] =
              new TH1F(Form("h_tdc_TWCorr_%s_plane_%s_bar_%d_thread_%d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       Form("TDC TWCorr %s Plane %s Bar %d Thread %d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       time_params.NBINS, time_params.MIN, time_params.MAX);
          hist_map_vec[thread]["h_tdc_TWCorr_bar"][plane][side][bar]->GetXaxis()->SetTitle("TDC Time (ns)");
          hist_map_vec[thread]["h_tdc_TWCorr_bar"][plane][side][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_vec[thread]["h_tdc_Corr_bar"][plane][side][bar] =
              new TH1F(Form("h_tdc_Corr_%s_plane_%s_bar_%d_thread_%d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       Form("TDC Corr %s Plane %s Bar %d Thread %d", side_names[side].c_str(),
                            plane_names[plane].c_str(), bar, thread),
                       time_params.NBINS, time_params.MIN, time_params.MAX);
          hist_map_vec[thread]["h_tdc_Corr_bar"][plane][side][bar]->GetXaxis()->SetTitle("TDC Time (ns)");
          hist_map_vec[thread]["h_tdc_Corr_bar"][plane][side][bar]->GetYaxis()->SetTitle("Counts");

          // Initializing side 0 and 1 for following histograms to prevent segfaults. Only side 0 is used (data values
          // are avg of top and btm)
          if (side == 1) {
            hist_map_vec[thread]["h_tdc_avg_bar"][plane][side][bar] =
                new TH1F(Form("h_tdc_avg_dummy_plane_%d_bar_%d_thread_%d", plane, bar, thread),
                         Form("Dummy TDC Avg Plane %d Bar %d Thread %d", plane, bar, thread), 1, 0, 1);
            hist_map_vec[thread]["h_edep_avg_bar"][plane][side][bar] =
                new TH1F(Form("h_edep_avg_dummy_plane_%d_bar_%d_thread_%d", plane, bar, thread),
                         Form("Dummy EDEP Avg Plane %d Bar %d Thread %d", plane, bar, thread), 1, 0, 1);
            hist_map_vec[thread]["h_tdc_avg_vs_edep_bar"][plane][side][bar] =
                new TH2F(Form("h_tdc_avg_vs_edep_dummy_plane_%d_bar_%d_thread_%d", plane, bar, thread),
                         Form("Dummy TDC Avg vs EDEP Plane %d Bar %d Thread %d", plane, bar, thread), 1, 0, 1, 1, 0, 1);
          } else {
            hist_map_vec[thread]["h_tdc_avg_bar"][plane][side][bar] =
                new TH1F(Form("h_tdc_avg_plane_%s_bar_%d_thread_%d", plane_names[plane].c_str(), bar, thread),
                         Form("TDC Avg Plane %s Bar %d Thread %d", plane_names[plane].c_str(), bar, thread),
                         time_params.NBINS, time_params.MIN, time_params.MAX);
            hist_map_vec[thread]["h_tdc_avg_bar"][plane][side][bar]->GetXaxis()->SetTitle("TDC Time (ns)");
            hist_map_vec[thread]["h_tdc_avg_bar"][plane][side][bar]->GetYaxis()->SetTitle("Counts");

            hist_map_vec[thread]["h_edep_avg_bar"][plane][side][bar] =
                new TH1F(Form("h_edep_avg_plane_%s_bar_%d_thread_%d", plane_names[plane].c_str(), bar, thread),
                         Form("EDEP Avg Plane %s Bar %d Thread %d", plane_names[plane].c_str(), bar, thread),
                         edep_params.NBINS, edep_params.MIN, edep_params.MAX);
            hist_map_vec[thread]["h_edep_avg_bar"][plane][side][bar]->GetXaxis()->SetTitle("Energy Deposition (a.u.)");
            hist_map_vec[thread]["h_edep_avg_bar"][plane][side][bar]->GetYaxis()->SetTitle("Counts");

            hist_map_vec[thread]["h_tdc_avg_vs_edep_bar"][plane][side][bar] =
                new TH2F(Form("h_tdc_avg_vs_edep_plane_%s_bar_%d_thread_%d", plane_names[plane].c_str(), bar, thread),
                         Form("TDC Avg vs EDEP Plane %s Bar %d Thread %d", plane_names[plane].c_str(), bar, thread),
                         time_params.NBINS, time_params.MIN, time_params.MAX, edep_params.NBINS, edep_params.MIN,
                         edep_params.MAX);
            hist_map_vec[thread]["h_tdc_avg_vs_edep_bar"][plane][side][bar]->GetXaxis()->SetTitle("TDC Time (ns)");
            hist_map_vec[thread]["h_tdc_avg_vs_edep_bar"][plane][side][bar]->GetYaxis()->SetTitle(
                "Energy Deposition (a.u.)");
          }
        }
        for (int i_trk_cut = 0; i_trk_cut < nTrackCuts; i_trk_cut++) {
          hist_map_trk_cuts_vec[thread]["h_ladkin_time_bar"][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_time_plane_%s_bar_%d_%s_thread_%d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN Time Plane %s Bar %d %s Thread %d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       time_params.NBINS, time_params.MIN, time_params.MAX);
          hist_map_trk_cuts_vec[thread]["h_ladkin_time_bar"][i_trk_cut][plane][bar]->GetXaxis()->SetTitle("Time (ns)");
          hist_map_trk_cuts_vec[thread]["h_ladkin_time_bar"][i_trk_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_trk_cuts_vec[thread]["h_ladkin_tof_bar"][i_trk_cut][plane][bar] =
            new TH1F(Form("h_ladkin_tof_plane_%s_bar_%d_%s_thread_%d", plane_names[plane].c_str(), bar,
                    track_cuts[i_trk_cut].cut_name.c_str(), thread),
                 Form("LADKIN TOF Plane %s Bar %d %s Thread %d", plane_names[plane].c_str(), bar,
                    track_cuts[i_trk_cut].cut_name.c_str(), thread),
                 tof_params.NBINS, tof_params.MIN, tof_params.MAX);
          hist_map_trk_cuts_vec[thread]["h_ladkin_tof_bar"][i_trk_cut][plane][bar]->GetXaxis()->SetTitle("TOF (ns)");
          hist_map_trk_cuts_vec[thread]["h_ladkin_tof_bar"][i_trk_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_trk_cuts_vec[thread]["h_ladkin_edep_bar"][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_edep_plane_%s_bar_%d_%s_thread_%d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN EDEP Plane %s Bar %d %s Thread %d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       edep_params.NBINS, edep_params.MIN, edep_params.MAX);
          hist_map_trk_cuts_vec[thread]["h_ladkin_edep_bar"][i_trk_cut][plane][bar]->GetXaxis()->SetTitle(
              "Energy Deposition (a.u.)");
          hist_map_trk_cuts_vec[thread]["h_ladkin_edep_bar"][i_trk_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_trk_cuts_vec[thread]["h_ladkin_theta_bar"][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_theta_plane_%s_bar_%d_%s_thread_%d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN Theta Plane %s Bar %d %s Thread %d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       theta_params.NBINS, theta_params.MIN, theta_params.MAX);
          hist_map_trk_cuts_vec[thread]["h_ladkin_theta_bar"][i_trk_cut][plane][bar]->GetXaxis()->SetTitle(
              "Theta (degrees)");
          hist_map_trk_cuts_vec[thread]["h_ladkin_theta_bar"][i_trk_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_trk_cuts_vec[thread]["h_ladkin_phi_bar"][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_phi_plane_%s_bar_%d_%s_thread_%d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN Phi Plane %s Bar %d %s Thread %d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       phi_params.NBINS, phi_params.MIN, phi_params.MAX);
          hist_map_trk_cuts_vec[thread]["h_ladkin_phi_bar"][i_trk_cut][plane][bar]->GetXaxis()->SetTitle(
              "Phi (degrees)");
          hist_map_trk_cuts_vec[thread]["h_ladkin_phi_bar"][i_trk_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_trk_cuts_vec[thread]["h_ladkin_dTrkVert_bar"][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_dTrkVert_plane_%s_bar_%d_%s_thread_%d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN Delta Pos Long Plane %s Bar %d %s Thread %d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       dTrk_params.NBINS, dTrk_params.MIN, dTrk_params.MAX);
          hist_map_trk_cuts_vec[thread]["h_ladkin_dTrkVert_bar"][i_trk_cut][plane][bar]->GetXaxis()->SetTitle(
              "Delta Pos Long (cm)");
          hist_map_trk_cuts_vec[thread]["h_ladkin_dTrkVert_bar"][i_trk_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_trk_cuts_vec[thread]["h_ladkin_dTrkHoriz_bar"][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_dTrkHoriz_plane_%s_bar_%d_%s_thread_%d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN Delta Pos Trans Plane %s Bar %d %s Thread %d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       dTrk_params.NBINS, dTrk_params.MIN, dTrk_params.MAX);
          hist_map_trk_cuts_vec[thread]["h_ladkin_dTrkHoriz_bar"][i_trk_cut][plane][bar]->GetXaxis()->SetTitle(
              "Delta Pos Trans (cm)");
          hist_map_trk_cuts_vec[thread]["h_ladkin_dTrkHoriz_bar"][i_trk_cut][plane][bar]->GetYaxis()->SetTitle(
              "Counts");

          hist_map_trk_cuts_vec[thread]["h_edep_vs_time_bar"][i_trk_cut][plane][bar] = new TH2F(
              Form("h_edep_vs_time_plane_%s_bar_%d_%s_thread_%d", plane_names[plane].c_str(), bar,
                   track_cuts[i_trk_cut].cut_name.c_str(), thread),
              Form("EDEP vs Time Plane %s Bar %d %s Thread %d", plane_names[plane].c_str(), bar,
                   track_cuts[i_trk_cut].cut_name.c_str(), thread),
              time_params.NBINS, time_params.MIN, time_params.MAX, edep_params.NBINS, edep_params.MIN, edep_params.MAX);
          hist_map_trk_cuts_vec[thread]["h_edep_vs_time_bar"][i_trk_cut][plane][bar]->GetXaxis()->SetTitle("Time (ns)");
          hist_map_trk_cuts_vec[thread]["h_edep_vs_time_bar"][i_trk_cut][plane][bar]->GetYaxis()->SetTitle(
              "Energy Deposition Int (a.u.)");

          hist_map_trk_cuts_vec[thread]["h_edep_vs_tof_bar"][i_trk_cut][plane][bar] = new TH2F(
              Form("h_edep_vs_tof_plane_%s_bar_%d_%s_thread_%d", plane_names[plane].c_str(), bar,
                   track_cuts[i_trk_cut].cut_name.c_str(), thread),
              Form("EDEP vs TOF Plane %s Bar %d %s Thread %d", plane_names[plane].c_str(), bar,
                   track_cuts[i_trk_cut].cut_name.c_str(), thread),
              tof_params.NBINS, tof_params.MIN, tof_params.MAX, edep_params.NBINS, edep_params.MIN, edep_params.MAX);
          hist_map_trk_cuts_vec[thread]["h_edep_vs_tof_bar"][i_trk_cut][plane][bar]->GetXaxis()->SetTitle("TOF (ns)");
          hist_map_trk_cuts_vec[thread]["h_edep_vs_tof_bar"][i_trk_cut][plane][bar]->GetYaxis()->SetTitle(
              "Energy Deposition Int (a.u.)");

          hist_map_trk_cuts_vec[thread]["h_ladkin_d0_bar"][i_trk_cut][plane][bar] =
              new TH1F(Form("h_ladkin_d0_plane_%s_bar_%d_%s_thread_%d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN D0 Plane %s Bar %d %s Thread %d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       d0_params.NBINS, d0_params.MIN, d0_params.MAX);
          hist_map_trk_cuts_vec[thread]["h_ladkin_d0_bar"][i_trk_cut][plane][bar]->GetXaxis()->SetTitle("D0 (cm)");
          hist_map_trk_cuts_vec[thread]["h_ladkin_d0_bar"][i_trk_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_trk_cuts_vec[thread]["h_projz_bar"][i_trk_cut][plane][bar] =
              new TH1F(Form("h_projz_bar_plane_%s_bar_%d_%s_thread_%d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN ProjZ Bar Plane %s Bar %d %s Thread %d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       projz_params.NBINS, projz_params.MIN, projz_params.MAX);
          hist_map_trk_cuts_vec[thread]["h_projz_bar"][i_trk_cut][plane][bar]->GetXaxis()->SetTitle("ProjZ (cm)");
          hist_map_trk_cuts_vec[thread]["h_projz_bar"][i_trk_cut][plane][bar]->GetYaxis()->SetTitle("Counts");

          hist_map_trk_cuts_vec[thread]["h_XY_GEM0_bar"][i_trk_cut][plane][bar] =
              new TH2F(Form("h_XY_GEM0_bar_plane_%s_bar_%d_%s_thread_%d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN XY GEM0 Bar Plane %s Bar %d %s Thread %d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       gem_X_params.NBINS, gem_X_params.MIN, gem_X_params.MAX, gem_Y_params.NBINS, gem_Y_params.MIN,
                       gem_Y_params.MAX);
          hist_map_trk_cuts_vec[thread]["h_XY_GEM0_bar"][i_trk_cut][plane][bar]->GetXaxis()->SetTitle("X (cm)");
          hist_map_trk_cuts_vec[thread]["h_XY_GEM0_bar"][i_trk_cut][plane][bar]->GetYaxis()->SetTitle("Y (cm)");

          hist_map_trk_cuts_vec[thread]["h_XY_GEM1_bar"][i_trk_cut][plane][bar] =
              new TH2F(Form("h_XY_GEM1_bar_plane_%s_bar_%d_%s_thread_%d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       Form("LADKIN XY GEM1 Bar Plane %s Bar %d %s Thread %d", plane_names[plane].c_str(), bar,
                            track_cuts[i_trk_cut].cut_name.c_str(), thread),
                       gem_X_params.NBINS, gem_X_params.MIN, gem_X_params.MAX, gem_Y_params.NBINS, gem_Y_params.MIN,
                       gem_Y_params.MAX);
          hist_map_trk_cuts_vec[thread]["h_XY_GEM1_bar"][i_trk_cut][plane][bar]->GetXaxis()->SetTitle("X (cm)");
          hist_map_trk_cuts_vec[thread]["h_XY_GEM1_bar"][i_trk_cut][plane][bar]->GetYaxis()->SetTitle("Y (cm)");
        }
      }
    }
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
    threads.emplace_back(process_chunk, i_thread, start, end, ref(fileNames), ref(hist_map_vec[i_thread]),
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
              (hist_map_vec[0][hist_pair.first][plane][side][bar])
                  ->SetTitle(TString(hist_map_vec[0][hist_pair.first][plane][side][bar]->GetTitle())
                                 .ReplaceAll(Form(" Thread %d", i_thread), ""));
              (hist_map_vec[0][hist_pair.first][plane][side][bar])
                  ->SetName(TString(hist_map_vec[0][hist_pair.first][plane][side][bar]->GetName())
                                .ReplaceAll(Form("_thread_%d", i_thread), ""));
            } else {

              (hist_map_vec[0][hist_pair.first][plane][side][bar])->Add(hist_pair.second[plane][side][bar]);
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
              (hist_map_trk_cuts_vec[0][hist_pair.first][i_trk_cut][plane][bar])
                  ->SetTitle(TString(hist_map_trk_cuts_vec[0][hist_pair.first][i_trk_cut][plane][bar]->GetTitle())
                                 .ReplaceAll(Form(" Thread %d", i_thread), ""));
              (hist_map_trk_cuts_vec[0][hist_pair.first][i_trk_cut][plane][bar])
                  ->SetName(TString(hist_map_trk_cuts_vec[0][hist_pair.first][i_trk_cut][plane][bar]->GetName())
                                .ReplaceAll(Form("_thread_%d", i_thread), ""));
            } else {
              (hist_map_trk_cuts_vec[0][hist_pair.first][i_trk_cut][plane][bar])
                  ->Add(hist_pair.second[i_trk_cut][plane][bar]);
            }
          }
        }
      }
    }
  }

  cout << "Writing histograms to file..." << endl;
  // Write histograms to the output file
  outputFile->cd();
  write_to_canvas(hist_map_vec[0]["h_tdc_bar"], outputFile, "Raw/TDC_Time", "TDC");
  write_to_canvas(hist_map_vec[0]["h_adc_bar"], outputFile, "Raw/ADC_Time", "ADC");
  write_to_canvas(hist_map_vec[0]["h_adc_amp_bar"], outputFile, "Raw/Amp", "ADC_Amp");
  write_to_canvas(hist_map_vec[0]["h_adc_int_bar"], outputFile, "Raw/Int", "ADC_Int");
  write_to_canvas(hist_map_vec[0]["h_adc_good_bar"], outputFile, "Good/ADC_Time", "ADC_Good");
  write_to_canvas(hist_map_vec[0]["h_adc_amp_good_bar"], outputFile, "Good/Amp", "ADC_Amp_Good");
  write_to_canvas(hist_map_vec[0]["h_adc_int_good_bar"], outputFile, "Good/Int", "ADC_Int_Good");
  write_to_canvas(hist_map_vec[0]["h_tdc_UnCorr_bar"], outputFile, "Good/TDC_Time/UnCorr", "TDC_UnCorr");
  write_to_canvas(hist_map_vec[0]["h_tdc_TWCorr_bar"], outputFile, "Good/TDC_Time/TWCorr", "TDC_TWCorr");
  write_to_canvas(hist_map_vec[0]["h_tdc_Corr_bar"], outputFile, "Good/TDC_Time/Corr", "TDC_Corr");
  write_to_canvas(hist_map_vec[0]["h_tdc_avg_bar"], outputFile, "Good/TDC_Time/Average", "TDC_Avg", true);
  write_to_canvas(hist_map_vec[0]["h_edep_avg_bar"], outputFile, "Good/edep/", "EDEP_Avg", true);
  write_to_canvas(hist_map_vec[0]["h_tdc_avg_vs_edep_bar"], outputFile, "Good/TDC_vs_Edep", "TDC_vs_Edep", true);

  for (int i_trk_cut = 0; i_trk_cut < nTrackCuts; ++i_trk_cut) {
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_time_bar"][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_Time/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_Time_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_tof_bar"][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_TOF/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_TOF_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_edep_bar"][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_Edep/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_Edep_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_theta_bar"][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_Theta/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_Theta_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_phi_bar"][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_Phi/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_Phi_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_dTrkVert_bar"][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_dTrkVert/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_dTrkVert_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_dTrkHoriz_bar"][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_dTrkHoriz/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_dTrkHoriz_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_edep_vs_time_bar"][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_Edep_vs_Time/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_Edep_vs_Time_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_edep_vs_tof_bar"][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_Edep_vs_TOF/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_Edep_vs_TOF_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_ladkin_d0_bar"][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_D0/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_D0_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_projz_bar"][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_ProjZ/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_ProjZ_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_XY_GEM0_bar"][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_XY_GEM0/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_XY_GEM0_%s", track_cuts[i_trk_cut].cut_name.c_str()));
    write_to_canvas_plane(hist_map_trk_cuts_vec[0]["h_XY_GEM1_bar"][i_trk_cut], outputFile,
                          Form("KIN/LADKIN_XY_GEM1/%s", track_cuts[i_trk_cut].cut_name.c_str()),
                          Form("LADKIN_XY_GEM1_%s", track_cuts[i_trk_cut].cut_name.c_str()));
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

  std::_Exit(0);
  // return;
  
  // gROOT->Reset();

  // return;
  // Clean up
  // Start threads to delete histograms in parallel
  std::vector<std::thread> cleanup_threads;
  for (int thread = 0; thread < numThreads; ++thread) {
    cleanup_threads.emplace_back([&, thread]() {
      size_t total_hist = 0;
      size_t cleaned_hist = 0;
      // Count total histograms for progress bar (only for thread 0)
      if (thread == 0) {
        for (auto &hist_pair : hist_map_vec[thread]) {
          total_hist += N_PLANES * N_SIDES * N_PADDLES;
        }
        for (auto &hist_pair : hist_map_trk_cuts_vec[thread]) {
          total_hist += nTrackCuts * N_PLANES * N_PADDLES;
        }
      }
      for (auto &hist_pair : hist_map_vec[thread]) {
        for (int plane = 0; plane < N_PLANES; ++plane) {
          for (int side = 0; side < N_SIDES; ++side) {
            for (int bar = 0; bar < N_PADDLES; ++bar) {
              delete hist_pair.second[plane][side][bar];
              if (thread == 0) {
                ++cleaned_hist;
                if (cleaned_hist % 100 == 0 || cleaned_hist == total_hist) {
                  int percent = int(100.0 * cleaned_hist / total_hist);
                  std::cout << "\rCleanup progress: " << percent << "% completed." << std::flush;
                }
              }
            }
          }
        }
      }
      for (auto &hist_pair : hist_map_trk_cuts_vec[thread]) {
        for (int i_trk_cut = 0; i_trk_cut < nTrackCuts; ++i_trk_cut) {
          for (int plane = 0; plane < N_PLANES; ++plane) {
            for (int bar = 0; bar < N_PADDLES; ++bar) {
              delete hist_pair.second[i_trk_cut][plane][bar];
              if (thread == 0) {
                ++cleaned_hist;
                if (cleaned_hist % 100 == 0 || cleaned_hist == total_hist) {
                  int percent = int(100.0 * cleaned_hist / total_hist);
                  std::cout << "\rCleanup progress: " << percent << "% completed." << std::flush;
                }
              }
            }
          }
        }
      }
      if (thread == 0) {
        std::cout << "\rCleanup progress: 100% completed." << std::endl;
      }
    });
  }
  for (auto &t : cleanup_threads) {
    t.join();
  }

  cout << "Done!" << endl;
  // Close the output file

  return;
}
