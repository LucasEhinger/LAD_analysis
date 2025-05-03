#include <TBranch.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <algorithm>
#include <string>

const static int NDATAMAX       = 1000;
const int NUM_BARS              = 11;
const int NUM_BARS_REF          = 1;
const int NUM_PLANES            = 6;
const std::string plane_names[] = {"000", "001", "100", "101", "200", "REFBAR"};

const int MAX_N_TRACKS                = 10;
const int MAX_N_GOODHITS              = 40;
const static int MAX_TOTAL_HITS       = 3;
const static int MAX_N_HITS_PER_PLANE = 2;

const int NBINS    = 200;
const double N_MIN = -500;
const double N_MAX = 500;

const static int NBINS_T_DIFF  = 100;
const static double T_DIFF_MIN = -50;
const static double T_DIFF_MAX = 50;
void get_good_evts(int settingNumber) {
  // int settingNumber = 9; // Replace this with the desired setting number
  const char *inputFileName = Form(
      "/u/home/ehingerl/hallc/software/lad_replay/ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_296_-1_setting%d.root",
      settingNumber);
  const char *outputFileName = Form(
      "/u/home/ehingerl/hallc/software/lad_replay/ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_296_-1_setting%d_filtered.root",
      settingNumber);

  // Open the input ROOT file
  TFile *inputFile = TFile::Open(inputFileName, "READ");
  if (!inputFile || inputFile->IsZombie()) {
    printf("Error opening input file: %s\n", inputFileName);
    return;
  }

  // Get the tree "T" from the input file
  TTree *inputTree = (TTree *)inputFile->Get("T");
  if (!inputTree) {
    printf("Error: Tree 'T' not found in input file\n");
    inputFile->Close();
    return;
  }

  // Create the output ROOT file
  TFile *outputFile = TFile::Open(outputFileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    printf("Error creating output file: %s\n", outputFileName);
    inputFile->Close();
    return;
  }

  // Clone the input tree structure to the output tree
  TTree *outputTree = inputTree->CloneTree(0);

  // Define the criteria for selecting events
  Double_t nTracks;
  Double_t L_gem_trk_x1_local[NDATAMAX];

  int Ndata_L_gem_hit_nclustu_layer;
  Double_t L_gem_hit_nclustu_layer[NDATAMAX];

  int Ndata_L_gem_hit_nclustv_layer;
  Double_t L_gem_hit_nclustv_layer[NDATAMAX];

  Double_t L_ladhod_goodhit_n;
  Double_t L_ladhod_goodhit_delta_pos_long[NDATAMAX];
  Double_t L_ladhod_goodhit_delta_pos_trans[NDATAMAX];
  Double_t L_ladhod_goodhit_plane[NDATAMAX];
  Double_t L_ladhod_goodhit_paddle[NDATAMAX];

  Double_t nHits[NUM_PLANES]                       = {0};
  Double_t Btm_ADC_Int[NUM_PLANES][NUM_BARS]       = {0};
  Double_t Top_ADC_Int[NUM_PLANES][NUM_BARS]       = {0};
  Double_t Btm_ADC_Amp[NUM_PLANES][NUM_BARS]       = {0};
  Double_t Top_ADC_Amp[NUM_PLANES][NUM_BARS]       = {0};
  Double_t Btm_TDC_Time[NUM_PLANES][NUM_BARS]      = {0};
  Double_t Top_TDC_Time[NUM_PLANES][NUM_BARS]      = {0};
  Double_t Btm_ADC_Ped[NUM_PLANES][NUM_BARS]       = {0};
  Double_t Top_ADC_Ped[NUM_PLANES][NUM_BARS]       = {0};
  Double_t Good_Hit_Time_Avg[NUM_PLANES][NUM_BARS] = {0};

  inputTree->SetBranchAddress("L.gem.trk.ntracks", &nTracks);
  inputTree->SetBranchAddress("L.gem.trk.x1_local", L_gem_trk_x1_local);

  inputTree->SetBranchAddress("Ndata.L.gem.hit.nclustu_layer", &Ndata_L_gem_hit_nclustu_layer);
  inputTree->SetBranchAddress("L.gem.hit.nclustu_layer", L_gem_hit_nclustu_layer);

  inputTree->SetBranchAddress("Ndata.L.gem.hit.nclustv_layer", &Ndata_L_gem_hit_nclustv_layer);
  inputTree->SetBranchAddress("L.gem.hit.nclustv_layer", L_gem_hit_nclustv_layer);

  inputTree->SetBranchAddress("L.ladhod.goodhit_n", &L_ladhod_goodhit_n);
  inputTree->SetBranchAddress("L.ladhod.goodhit_delta_pos_long", L_ladhod_goodhit_delta_pos_long);
  inputTree->SetBranchAddress("L.ladhod.goodhit_delta_pos_trans", L_ladhod_goodhit_delta_pos_trans);
  inputTree->SetBranchAddress("L.ladhod.goodhit_plane", L_ladhod_goodhit_plane);
  inputTree->SetBranchAddress("L.ladhod.goodhit_paddle", L_ladhod_goodhit_paddle);

  for (int plane = 0; plane < NUM_PLANES; ++plane) {
    inputTree->SetBranchAddress(Form("L.ladhod.%s.nhits", plane_names[plane].c_str()), &nHits[plane]);
    inputTree->SetBranchAddress(Form("L.ladhod.%s.GoodBtmAdcPulseInt", plane_names[plane].c_str()), Btm_ADC_Int[plane]);
    inputTree->SetBranchAddress(Form("L.ladhod.%s.GoodTopAdcPulseInt", plane_names[plane].c_str()), Top_ADC_Int[plane]);
    inputTree->SetBranchAddress(Form("L.ladhod.%s.GoodBtmAdcPulseAmp", plane_names[plane].c_str()), Btm_ADC_Amp[plane]);
    inputTree->SetBranchAddress(Form("L.ladhod.%s.GoodTopAdcPulseAmp", plane_names[plane].c_str()), Top_ADC_Amp[plane]);
    inputTree->SetBranchAddress(Form("L.ladhod.%s.GoodBtmTdcTimeWalkCorr", plane_names[plane].c_str()),
                                Btm_TDC_Time[plane]);
    inputTree->SetBranchAddress(Form("L.ladhod.%s.GoodTopTdcTimeWalkCorr", plane_names[plane].c_str()),
                                Top_TDC_Time[plane]);
    inputTree->SetBranchAddress(Form("L.ladhod.%s.GoodBtmAdcPulsePed", plane_names[plane].c_str()), Btm_ADC_Ped[plane]);
    inputTree->SetBranchAddress(Form("L.ladhod.%s.GoodTopAdcPulsePed", plane_names[plane].c_str()), Top_ADC_Ped[plane]);
    inputTree->SetBranchAddress(Form("L.ladhod.%s.GoodHitTimeAvg", plane_names[plane].c_str()),
                                Good_Hit_Time_Avg[plane]);
  }

  Long64_t nEntries = inputTree->GetEntries();

  TH1D *hist_delta_pos_long  = new TH1D("hist_delta_pos_long", "Delta Pos Long;Value;Entries", NBINS, N_MIN, N_MAX);
  TH1D *hist_delta_pos_trans = new TH1D("hist_delta_pos_trans", "Delta Pos Trans;Value;Entries", NBINS, N_MIN, N_MAX);
  TH1D *hist_min_delta_pos_long =
      new TH1D("hist_min_delta_pos_long", "Min Delta Pos Long;Value;Entries", NBINS, N_MIN, N_MAX);
  TH1D *hist_min_delta_pos_trans =
      new TH1D("hist_min_delta_pos_trans", "Min Delta Pos Trans;Value;Entries", NBINS, N_MIN, N_MAX);
  TH1D *hist_delta_pos_long_vs_plane[NUM_PLANES];
  TH1D *hist_delta_pos_trans_vs_plane[NUM_PLANES];
  TH1D *hist_delta_pos_trans_vs_plane_unique[NUM_PLANES];

  for (int plane = 0; plane < NUM_PLANES; ++plane) {
    hist_delta_pos_long_vs_plane[plane] =
        new TH1D(Form("hist_delta_pos_long_vs_plane_%d", plane),
                 Form("Delta Pos Long for Plane %d;Delta Pos Long;Entries", plane), NBINS, N_MIN, N_MAX);
    hist_delta_pos_trans_vs_plane[plane] =
        new TH1D(Form("hist_delta_pos_trans_vs_plane_%d", plane),
                 Form("Delta Pos Trans for Plane %d;Delta Pos Trans;Entries", plane), NBINS, N_MIN, N_MAX);
    hist_delta_pos_trans_vs_plane_unique[plane] =
        new TH1D(Form("hist_delta_pos_trans_vs_plane_unique_%d", plane),
                 Form("Delta Pos Trans for Unique Plane %d;Delta Pos Trans;Entries", plane), NBINS, N_MIN, N_MAX);
  }

  TH1D *hist_time_diff_refbar_200 =
      new TH1D("hist_time_diff_refbar_200", "Time Diff Refbar f- Hodo 200;Time Difference;Entries", NBINS_T_DIFF,
               T_DIFF_MIN, T_DIFF_MAX);
  TH1D *hist_time_diff_refbar_100 =
      new TH1D("hist_time_diff_refbar_100", "Time Diff Refbar - Hodo 100;Time Difference;Entries", NBINS_T_DIFF,
               T_DIFF_MIN, T_DIFF_MAX);
  TH1D *hist_time_diff_refbar_101 =
      new TH1D("hist_time_diff_refbar_101", "Time Diff Refbar - Hodo 101;Time Difference;Entries", NBINS_T_DIFF,
               T_DIFF_MIN, T_DIFF_MAX);
  TH1D *hist_delta_pos_long_vs_bar[NUM_PLANES][NUM_BARS];
  TH1D *hist_delta_pos_trans_vs_bar[NUM_PLANES][NUM_BARS];
  TH1D *hist_delta_pos_trans_vs_bar_unique[NUM_PLANES][NUM_BARS];

  for (int plane = 0; plane < NUM_PLANES; ++plane) {
    for (int bar = 0; bar < NUM_BARS; ++bar) {
      hist_delta_pos_long_vs_bar[plane][bar] =
          new TH1D(Form("hist_delta_pos_long_vs_bar_plane_%d_bar_%d", plane, bar),
                   Form("Delta Pos Long for Plane %d Bar %d;Delta Pos Long;Entries", plane, bar), NBINS, N_MIN, N_MAX);
      hist_delta_pos_trans_vs_bar[plane][bar] = new TH1D(
          Form("hist_delta_pos_trans_vs_bar_plane_%d_bar_%d", plane, bar),
          Form("Delta Pos Trans for Plane %d Bar %d;Delta Pos Trans;Entries", plane, bar), NBINS, N_MIN, N_MAX);
      hist_delta_pos_trans_vs_bar_unique[plane][bar] = new TH1D(
          Form("hist_delta_pos_trans_vs_bar_unique_plane_%d_bar_%d", plane, bar),
          Form("Delta Pos Trans for Unique Plane %d Bar %d;Delta Pos Trans;Entries", plane, bar), NBINS, N_MIN, N_MAX);
    }
  }

  // Loop over all events in the input tree
  for (Long64_t i = 0; i < nEntries; ++i) {
    inputTree->GetEntry(i);

    // Calculate what percentage of events meet certain qualities
    int refbar_hits    = nHits[5];
    int plane_100_hits = nHits[2];
    int plane_101_hits = nHits[3];
    int plane_200_hits = nHits[4];

    static int total_events                      = 0;
    static int refbar_events                     = 0;
    static int refbar_100_101_events             = 0;
    static int refbar_100_101_less_n_hits_events = 0;
    static int refbar_200_events                 = 0;
    static int refbar_200_less_n_hits_events     = 0;

    total_events++;

    if (refbar_hits > 0) {
      refbar_events++;
      if (plane_100_hits > 0 && plane_101_hits > 0) {
        refbar_100_101_events++;
        if (plane_100_hits < MAX_N_HITS_PER_PLANE && plane_101_hits < MAX_N_HITS_PER_PLANE) {
          refbar_100_101_less_n_hits_events++;
        }
      }
      if (plane_200_hits > 0) {
        refbar_200_events++;
        if (plane_200_hits < MAX_N_HITS_PER_PLANE) {
          refbar_200_less_n_hits_events++;
        }
      }
    }

    if (i == nEntries - 1) {
      printf("Percentage of events that hit the refbar: %.2f%%\n", (double)refbar_events / total_events * 100);
      printf("Percentage of events that hit the refbar and plane 100 and plane 101: %.2f%%\n",
             (double)refbar_100_101_events / total_events * 100);
      printf("Percentage of events that hit the refbar and plane 100 and plane 101, but each plane has less than %d "
             "hits: %.2f%%\n",
             MAX_N_HITS_PER_PLANE, (double)refbar_100_101_less_n_hits_events / total_events * 100);
      printf("Percentage of events that hit the refbar and plane 200: %.2f%%\n",
             (double)refbar_200_events / total_events * 100);
      printf("Percentage of events that hit the refbar and plane 200, but each plane has less than %d hits: %.2f%%\n",
             MAX_N_HITS_PER_PLANE, (double)refbar_200_less_n_hits_events / total_events * 100);
    }
    // End calculating percentage of events that meet certain qualities

    // Fill refbar histos
    for (int j = 0; j < NUM_BARS; ++j) {
      for (int k = 0; k < NUM_BARS_REF; ++k) {
        if ((Good_Hit_Time_Avg[4][j] - Good_Hit_Time_Avg[5][k]) != 0) {
          hist_time_diff_refbar_200->Fill(Good_Hit_Time_Avg[4][j] - Good_Hit_Time_Avg[5][k]);
        }
        if ((Good_Hit_Time_Avg[2][j] - Good_Hit_Time_Avg[5][k]) != 0) {
          hist_time_diff_refbar_100->Fill(Good_Hit_Time_Avg[2][j] - Good_Hit_Time_Avg[5][k]);
        }

        if ((Good_Hit_Time_Avg[3][j] - Good_Hit_Time_Avg[5][k]) != 0) {
          hist_time_diff_refbar_101->Fill(Good_Hit_Time_Avg[3][j] - Good_Hit_Time_Avg[5][k]);
        }
      }
    }

    // Fill histograms for all "good" entries
    bool goodEvent = true;
    if (nTracks > MAX_N_TRACKS || L_ladhod_goodhit_n > MAX_N_GOODHITS ||
        (refbar_hits + plane_100_hits + plane_101_hits + plane_200_hits) > MAX_TOTAL_HITS) {
      goodEvent = false;
    }
    if (nTracks == 0) {
      goodEvent = false;
    }

    if (goodEvent) {
      // Fill histograms for all entries
      for (int j = 0; j < L_ladhod_goodhit_n; ++j) {
        hist_delta_pos_long->Fill(L_ladhod_goodhit_delta_pos_long[j]);
        hist_delta_pos_trans->Fill(L_ladhod_goodhit_delta_pos_trans[j]);
        hist_delta_pos_long_vs_plane[int(L_ladhod_goodhit_plane[j])]->Fill(L_ladhod_goodhit_delta_pos_long[j]);
        hist_delta_pos_trans_vs_plane[int(L_ladhod_goodhit_plane[j])]->Fill(L_ladhod_goodhit_delta_pos_trans[j]);
        hist_delta_pos_long_vs_bar[int(L_ladhod_goodhit_plane[j])][int(L_ladhod_goodhit_paddle[j])]->Fill(
            L_ladhod_goodhit_delta_pos_long[j]);
        hist_delta_pos_trans_vs_bar[int(L_ladhod_goodhit_plane[j])][int(L_ladhod_goodhit_paddle[j])]->Fill(
            L_ladhod_goodhit_delta_pos_trans[j]);
      }

      // Find and fill histograms for minimum entries per paddle
      if (L_ladhod_goodhit_n > 0) {
        std::map<int, Double_t> minDeltaPosTransPerPaddle;

        for (int j = 0; j < L_ladhod_goodhit_n; ++j) {
          int plane     = int(L_ladhod_goodhit_plane[j]);
          int paddle    = int(L_ladhod_goodhit_paddle[j]);
          int uniqueKey = plane * 100 + paddle; // Combine plane and paddle into a unique key

          // Update the minimum delta_pos_trans for this paddle
          if (minDeltaPosTransPerPaddle.find(uniqueKey) == minDeltaPosTransPerPaddle.end() ||
              std::abs(L_ladhod_goodhit_delta_pos_trans[j]) < std::abs(minDeltaPosTransPerPaddle[uniqueKey])) {
            minDeltaPosTransPerPaddle[uniqueKey] = L_ladhod_goodhit_delta_pos_trans[j];
          }

          // // Fill unique hits into hist_delta_pos_trans_vs_bar_unique
          // if (plane < NUM_PLANES && paddle < NUM_BARS) {
          //   hist_delta_pos_trans_vs_bar_unique[plane][paddle]->Fill(L_ladhod_goodhit_delta_pos_trans[j]);
          // }
        }

        // Fill histograms with the minimum delta_pos_trans for each paddle
        for (const auto &entry : minDeltaPosTransPerPaddle) {
          hist_min_delta_pos_trans->Fill(entry.second);
          hist_delta_pos_trans_vs_plane_unique[int(entry.first / 100)]->Fill(entry.second);
          hist_delta_pos_trans_vs_bar_unique[entry.first / 100][entry.first % 100]->Fill(entry.second);
        }
      }
      // Map for unique hit multiplicity
      std::map<int, int> uniqueHitMultiplicity;
      for (int j = 0; j < L_ladhod_goodhit_n; ++j) {
        int plane     = int(L_ladhod_goodhit_plane[j]);
        int paddle    = int(L_ladhod_goodhit_paddle[j]);
        int uniqueKey = plane * 100 + paddle; // Combine plane and paddle into a unique key
        uniqueHitMultiplicity[uniqueKey]++;
      }

      // // Print or process the unique hit multiplicity if needed
      // for (const auto &entry : uniqueHitMultiplicity) {
      //   printf("Plane %d has %d unique hits\n", entry.first, entry.second);
      // }
      outputTree->Fill();
    }
  }

  // Write the output tree and histograms to the output file
  outputTree->Write();
  hist_delta_pos_long->Write();
  hist_delta_pos_trans->Write();
  hist_min_delta_pos_long->Write();
  hist_min_delta_pos_trans->Write();
  hist_time_diff_refbar_200->Write();
  hist_time_diff_refbar_100->Write();
  hist_time_diff_refbar_101->Write();
  outputFile->mkdir("per_plane");
  outputFile->cd("per_plane");
  outputFile->mkdir("per_plane/long");
  outputFile->mkdir("per_plane/trans");
  outputFile->mkdir("per_plane/trans_unique");

  for (int plane = 0; plane < NUM_PLANES; ++plane) {
    if (hist_delta_pos_long_vs_plane[plane]->GetEntries() > 0) {
      outputFile->cd("per_plane/long");
      hist_delta_pos_long_vs_plane[plane]->Write();
    }
    if (hist_delta_pos_trans_vs_plane[plane]->GetEntries() > 0) {
      outputFile->cd("per_plane/trans");
      hist_delta_pos_trans_vs_plane[plane]->Write();
    }
    if (hist_delta_pos_trans_vs_plane_unique[plane]->GetEntries() > 0) {
      outputFile->cd("per_plane/trans_unique");
      hist_delta_pos_trans_vs_plane_unique[plane]->Write();
    }
  }
  outputFile->mkdir("per_paddle");
  outputFile->cd("per_paddle");
  outputFile->mkdir("per_paddle/long");
  outputFile->mkdir("per_paddle/trans");
  outputFile->mkdir("per_paddle/trans_unique");

  for (int plane = 0; plane < NUM_PLANES; ++plane) {
    for (int bar = 0; bar < NUM_BARS; ++bar) {
      if (hist_delta_pos_long_vs_bar[plane][bar]->GetEntries() > 0) {
        outputFile->cd("per_paddle/long");
        hist_delta_pos_long_vs_bar[plane][bar]->Write();
      }
      if (hist_delta_pos_trans_vs_bar[plane][bar]->GetEntries() > 0) {
        outputFile->cd("per_paddle/trans");
        hist_delta_pos_trans_vs_bar[plane][bar]->Write();
      }
      if (hist_delta_pos_trans_vs_bar_unique[plane][bar]->GetEntries() > 0) {
        outputFile->cd("per_paddle/trans_unique");
        hist_delta_pos_trans_vs_bar_unique[plane][bar]->Write();
      }
    }
  }

  // Close the files
  inputFile->Close();
  outputFile->Close();

  printf("Filtered events and histograms written to output file: %s\n", outputFileName);
}