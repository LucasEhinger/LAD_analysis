#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TMinuit.h>
#include <TTree.h>
#include <TVector3.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace std;

const double min_projy = -5.0;
const double max_projy = 0.0;
const double max_d0    = 30.0;
const double max_dHoriz = 80.0; // Maximum horizontal distance from LAD
const double max_dVert = 80.0; // Maximum vertical distance from LAD

// const double min_projz = -13.0;
// const double max_projz = 13.0;
const double projz_sigma       = 5.0;
const int nFixedz              = 3;
const double target_z[nFixedz] = {-10.0, 0.0, 10.0}; // Fixed z positions for the planes

const double TDC2NS = 0.09766; // ns per TDC channel

struct hist_range {
  double min;
  double max;
  int nbins;

  hist_range(double min_val, double max_val, int bins) : min(min_val), max(max_val), nbins(bins) {}
};

const hist_range TDC_TIME_RAW(0, 4000, 200);
const hist_range TARGED_POS(-20, 20, 100);
const hist_range LAD_POS(-100, 100, 100);

int run_number = 300000;

const int MAX_DATA                = 10000;
const int maxTracks               = 30;
const int nPlanes                 = 5;
const string plane_names[nPlanes] = {"000", "001", "100", "101", "200"};

const double DCA_XZ_MAX = 100.0; // Maximum DCA in XZ plane
const double DCA_YZ_MAX = 200.0; // Maximum DCA in YZ plane

const double plane_theta[nPlanes] = {150.0, 150.0, 127.0, 127.0, 104.0}; // Angle in degrees
const double plane_r[nPlanes]     = {615.0, 655.6, 523.0, 563.6, 615.0}; // Radius of the second point

const int NSTRIP_U = 128*12; // Number of strips in U plane
const int NSTRIP_V = 128*12*2; // Number of strips in V plane

const bool use_projz = true; // Fix z position to GEM projz

const double dx_min = -30.0;
const double dx_max = 30.0;
const int dx_NBINS  = 60;
const double dy_min = -30.0;
const double dy_max = 30.0;
const int dy_NBINS  = 60;



void gem_good_track_properties() {

  TString fileName = Form("/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22812_0_0_-1.root");

                     //  "LAD_COIN_22282_-1_inverted.root";
                     //  "LAD_COIN_22282_-1_500trks_good_timing.root";
                     //  "LAD_COIN_22383_0_0_500002.root";

  TString outputFileName = Form("files/good_track_histos/gem_good_track_histos_%d_-1_P.root", 22812);
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
  //Tracks
  Double_t trk_d0[MAX_DATA], trk_d0_good[MAX_DATA];
  Double_t trk_projz[MAX_DATA], trk_projy[MAX_DATA];
  Double_t trk_t[MAX_DATA], trk_dt[MAX_DATA];
  Double_t trk_x[2][MAX_DATA], trk_y[2][MAX_DATA], trk_z[2][MAX_DATA];
  Double_t trk_x_local[2][MAX_DATA], trk_y_local[2][MAX_DATA];
  Double_t trk_spID_0[2][MAX_DATA]; // Space Point ID for each track
  Double_t trk_spID_1[2][MAX_DATA]; // Space Point ID for each track
  // Space Points
  Double_t sp_time[MAX_DATA], sp_layer[MAX_DATA], sp_dt[MAX_DATA], sp_asym[MAX_DATA];
  Double_t sp_clusID1[MAX_DATA], sp_clusID2[MAX_DATA], sp_adc[MAX_DATA];
  // Clusters
  Double_t clust_time[MAX_DATA], clust_nstrip[MAX_DATA], clust_module[MAX_DATA];
  Double_t clust_maxstrip[MAX_DATA], clust_maxsamp[MAX_DATA], clust_maxadc[MAX_DATA], clust_adc[MAX_DATA];
  Double_t clust_index[MAX_DATA], clust_axis[MAX_DATA]; // Axis of the cluster
  Int_t nClusters;

  // TDC
  Double_t tdc_time_btm[nPlanes][MAX_DATA], tdc_time_top[nPlanes][MAX_DATA];
  Double_t tdc_counter_btm[nPlanes][MAX_DATA], tdc_counter_top[nPlanes][MAX_DATA];
  Int_t nTracks, nTdcTopHits[nPlanes], nTdcBtmHits[nPlanes];
  Double_t vertex_x, vertex_y, vertex_z;

  // Good Hits
  Double_t good_hit_time[MAX_DATA], good_hit_plane[MAX_DATA], good_hit_paddle[MAX_DATA];
  Double_t good_hit_track_ID[MAX_DATA], good_hit_dTrkHoriz[MAX_DATA], good_hit_dTrkVert[MAX_DATA];
  Int_t nGoodHits = 0;
  char spect_prefix = 'H'; // Default to 'H', can be changed to 'P' if needed

  T->SetBranchAddress(Form("Ndata.%c.gem.trk.d0", spect_prefix), &nTracks);
  T->SetBranchAddress(Form("%c.gem.trk.d0", spect_prefix), &trk_d0);
  T->SetBranchAddress(Form("%c.gem.trk.d0_good", spect_prefix), &trk_d0_good);
  T->SetBranchAddress(Form("%c.gem.trk.projz", spect_prefix), &trk_projz);
  T->SetBranchAddress(Form("%c.gem.trk.projy", spect_prefix), &trk_projy);
  T->SetBranchAddress(Form("%c.gem.trk.x1", spect_prefix), &trk_x[0]);
  T->SetBranchAddress(Form("%c.gem.trk.y1", spect_prefix), &trk_y[0]);
  T->SetBranchAddress(Form("%c.gem.trk.z1", spect_prefix), &trk_z[0]);
  T->SetBranchAddress(Form("%c.gem.trk.x2", spect_prefix), &trk_x[1]);
  T->SetBranchAddress(Form("%c.gem.trk.y2", spect_prefix), &trk_y[1]);
  T->SetBranchAddress(Form("%c.gem.trk.z2", spect_prefix), &trk_z[1]);
  T->SetBranchAddress(Form("%c.gem.trk.x1_local", spect_prefix), &trk_x_local[0]);
  T->SetBranchAddress(Form("%c.gem.trk.y1_local", spect_prefix), &trk_y_local[0]);
  T->SetBranchAddress(Form("%c.gem.trk.x2_local", spect_prefix), &trk_x_local[1]);
  T->SetBranchAddress(Form("%c.gem.trk.y2_local", spect_prefix), &trk_y_local[1]);
  T->SetBranchAddress(Form("%c.gem.trk.t", spect_prefix), &trk_t);
  T->SetBranchAddress(Form("%c.gem.trk.dt", spect_prefix), &trk_dt);
  T->SetBranchAddress(Form("%c.gem.trk.spID_0u", spect_prefix), &trk_spID_0[0]);
  T->SetBranchAddress(Form("%c.gem.trk.spID_0v", spect_prefix), &trk_spID_0[1]);
  T->SetBranchAddress(Form("%c.gem.trk.spID_1u", spect_prefix), &trk_spID_1[0]);
  T->SetBranchAddress(Form("%c.gem.trk.spID_1v", spect_prefix), &trk_spID_1[1]);

  T->SetBranchAddress(Form("%c.gem.sp.time", spect_prefix), &sp_time);
  T->SetBranchAddress(Form("%c.gem.sp.layer", spect_prefix), &sp_layer);
  T->SetBranchAddress(Form("%c.gem.sp.dt", spect_prefix), &sp_dt);
  T->SetBranchAddress(Form("%c.gem.sp.asym", spect_prefix), &sp_asym);
  T->SetBranchAddress(Form("%c.gem.sp.clusID1", spect_prefix), &sp_clusID1);
  T->SetBranchAddress(Form("%c.gem.sp.clusID2", spect_prefix), &sp_clusID2);
  T->SetBranchAddress(Form("%c.gem.sp.adc", spect_prefix), &sp_adc);

  T->SetBranchAddress(Form("Ndata.%c.gem.clust.time", spect_prefix), &nClusters);
  T->SetBranchAddress(Form("%c.gem.clust.time", spect_prefix), &clust_time);
  T->SetBranchAddress(Form("%c.gem.clust.nstrip", spect_prefix), &clust_nstrip);
  T->SetBranchAddress(Form("%c.gem.clust.module", spect_prefix), &clust_module);
  T->SetBranchAddress(Form("%c.gem.clust.maxstrip", spect_prefix), &clust_maxstrip);
  T->SetBranchAddress(Form("%c.gem.clust.maxsamp", spect_prefix), &clust_maxsamp);
  T->SetBranchAddress(Form("%c.gem.clust.maxadc", spect_prefix), &clust_maxadc);
  T->SetBranchAddress(Form("%c.gem.clust.adc", spect_prefix), &clust_adc);
  T->SetBranchAddress(Form("%c.gem.clust.index", spect_prefix), &clust_index);
  T->SetBranchAddress(Form("%c.gem.clust.axis", spect_prefix), &clust_axis);


  T->SetBranchAddress(Form("%c.react.x", spect_prefix), &vertex_x);
  T->SetBranchAddress(Form("%c.react.y", spect_prefix), &vertex_y);
  T->SetBranchAddress(Form("%c.react.z", spect_prefix), &vertex_z);

  T->SetBranchAddress(Form("Ndata.%c.ladkin.goodhit_dTrkHoriz_0", spect_prefix), &nGoodHits);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_hittime_0", spect_prefix), &good_hit_time);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_plane_0", spect_prefix), &good_hit_plane);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_paddle_0", spect_prefix), &good_hit_paddle);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_trackid_0", spect_prefix), &good_hit_track_ID);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_dTrkHoriz_0", spect_prefix), &good_hit_dTrkHoriz);
  T->SetBranchAddress(Form("%c.ladkin.goodhit_dTrkVert_0", spect_prefix), &good_hit_dTrkVert);


  for (int i = 0; i < nPlanes; ++i) {
    T->SetBranchAddress(Form("%c.ladhod.%s.BtmTdcTimeRaw", spect_prefix, plane_names[i].c_str()), &tdc_time_btm[i]);
    T->SetBranchAddress(Form("%c.ladhod.%s.BtmTdcCounter", spect_prefix, plane_names[i].c_str()), &tdc_counter_btm[i]);
    T->SetBranchAddress(Form("%c.ladhod.%s.TopTdcTimeRaw", spect_prefix, plane_names[i].c_str()), &tdc_time_top[i]);
    T->SetBranchAddress(Form("%c.ladhod.%s.TopTdcCounter", spect_prefix, plane_names[i].c_str()), &tdc_counter_top[i]);
    T->SetBranchAddress(Form("Ndata.%c.ladhod.%s.BtmTdcTime", spect_prefix, plane_names[i].c_str()), &nTdcBtmHits[i]);
    T->SetBranchAddress(Form("Ndata.%c.ladhod.%s.TopTdcTime", spect_prefix, plane_names[i].c_str()), &nTdcTopHits[i]);
  }

  //////////////////////////////////////////////////////////////////
  // Create histograms
  TH1F *h_d0 = new TH1F("h_d0", "DCA in XZ plane; DCA (cm); Counts", 100, -max_d0, max_d0);
  TH1F *h_projy = new TH1F("h_projy", "DCA in YZ plane; DCA (cm); Counts", 100, min_projy, max_projy);
  TH1F *h_projz = new TH1F("h_projz", "Projection in Z; Projection (cm); Counts", 100, -20.0, 20.0);

  // Cluster-level histograms
  TH1F *h_clust_time_u[2];
  TH1F *h_clust_time_v[2];
  TH1F *h_clust_nstrip_u[2];
  TH1F *h_clust_nstrip_v[2];
  TH1F *h_clust_maxstrip_u[2];
  TH1F *h_clust_maxstrip_v[2];
  TH1F *h_clust_maxsamp_u[2];
  TH1F *h_clust_maxsamp_v[2];
  TH1F *h_clust_maxadc_u[2];
  TH1F *h_clust_maxadc_v[2];
  TH1F *h_clust_adc_u[2];
  TH1F *h_clust_adc_v[2];

  for (int i = 0; i < 2; ++i) {
    h_clust_time_u[i]     = new TH1F(Form("h_clust_time_u_%d", i), Form("Cluster Time U GEM %d; Time (ns); Counts", i), 100, 0, 150);
    h_clust_time_v[i]     = new TH1F(Form("h_clust_time_v_%d", i), Form("Cluster Time V GEM %d; Time (ns); Counts", i), 100, 0, 150);
    h_clust_nstrip_u[i]   = new TH1F(Form("h_clust_nstrip_u_%d", i), Form("Cluster NStrip U GEM %d; NStrip; Counts", i), 100, 0, NSTRIP_U);
    h_clust_nstrip_v[i]   = new TH1F(Form("h_clust_nstrip_v_%d", i), Form("Cluster NStrip V GEM %d; NStrip; Counts", i), 100, 0, NSTRIP_V);
    h_clust_maxstrip_u[i] = new TH1F(Form("h_clust_maxstrip_u_%d", i), Form("Cluster MaxStrip U GEM %d; MaxStrip; Counts", i), 100, 0, NSTRIP_U);
    h_clust_maxstrip_v[i] = new TH1F(Form("h_clust_maxstrip_v_%d", i), Form("Cluster MaxStrip V GEM %d; MaxStrip; Counts", i), 100, 0, NSTRIP_V);
    h_clust_maxsamp_u[i]  = new TH1F(Form("h_clust_maxsamp_u_%d", i), Form("Cluster MaxSamp U GEM %d; MaxSamp; Counts", i), 6, -0.5, 5.5);
    h_clust_maxsamp_v[i]  = new TH1F(Form("h_clust_maxsamp_v_%d", i), Form("Cluster MaxSamp V GEM %d; MaxSamp; Counts", i), 6, -0.5, 5.5);
    h_clust_maxadc_u[i]   = new TH1F(Form("h_clust_maxadc_u_%d", i), Form("Cluster MaxADC U GEM %d; MaxADC; Counts", i), 100, 0, 4096);
    h_clust_maxadc_v[i]   = new TH1F(Form("h_clust_maxadc_v_%d", i), Form("Cluster MaxADC V GEM %d; MaxADC; Counts", i), 100, 0, 4096);
    h_clust_adc_u[i]      = new TH1F(Form("h_clust_adc_u_%d", i), Form("Cluster ADC U GEM %d; ADC; Counts", i), 100, 0, 4096);
    h_clust_adc_v[i]      = new TH1F(Form("h_clust_adc_v_%d", i), Form("Cluster ADC V GEM %d; ADC; Counts", i), 100, 0, 4096);
  }

  ///////////////////////////////////////////////////////////////////

  // Loop through the tree entries
  Long64_t nEntries = T->GetEntries();
  // nEntries          = 10000; // For testing purposes, limit to 1000 entries
  for (Long64_t i = 0; i < nEntries; ++i) {
    T->GetEntry(i);
    // Loop through the good hits
    for (int j = 0; j < nGoodHits; ++j) {
      //skip good hit if it doesn't point back to target
      if( abs(good_hit_dTrkHoriz[j]) > max_dHoriz || 
          abs(good_hit_dTrkVert[j]) > max_dVert ) {
        continue;
      }
      int track_id = int(good_hit_track_ID[j]);
      if(track_id<0 || trk_d0[track_id] > max_d0) {
        continue; // Skip if track ID is invalid or d0 is out of range
      }

      // Fill the histograms for DCA and projection
      h_d0->Fill(trk_d0[track_id]);
      h_projy->Fill(trk_projy[track_id]);
      h_projz->Fill(trk_projz[track_id]);

      int cluster_id_0u = int(trk_spID_0[0][track_id]);
      int cluster_id_0v = int(trk_spID_0[1][track_id]);
      int cluster_id_1u = int(trk_spID_1[0][track_id]);
      int cluster_id_1v = int(trk_spID_1[1][track_id]);

      // Loop through all clusters
      for (int k = 0; k < nClusters; ++k) {
        // Check if the cluster belongs to the current track
        // Only fill if axis is 0 (U) or 1 (V) as appropriate
        if (clust_module[k] == 0 && clust_axis[k] == 0 && cluster_id_0u == clust_index[k]) {
          h_clust_time_u[0]->Fill(clust_time[k]);
          h_clust_nstrip_u[0]->Fill(clust_nstrip[k]);
          h_clust_maxstrip_u[0]->Fill(clust_maxstrip[k]);
          h_clust_maxsamp_u[0]->Fill(clust_maxsamp[k]);
          h_clust_maxadc_u[0]->Fill(clust_maxadc[k]);
          h_clust_adc_u[0]->Fill(clust_adc[k]);
        } else if (clust_module[k] == 0 && clust_axis[k] == 1 && cluster_id_0v == clust_index[k]) {
          h_clust_time_v[0]->Fill(clust_time[k]);
          h_clust_nstrip_v[0]->Fill(clust_nstrip[k]);
          h_clust_maxstrip_v[0]->Fill(clust_maxstrip[k]);
          h_clust_maxsamp_v[0]->Fill(clust_maxsamp[k]);
          h_clust_maxadc_v[0]->Fill(clust_maxadc[k]);
          h_clust_adc_v[0]->Fill(clust_adc[k]);
        } else if (clust_module[k] == 1 && clust_axis[k] == 0 && cluster_id_1u == clust_index[k]) {
          h_clust_time_u[1]->Fill(clust_time[k]);
          h_clust_nstrip_u[1]->Fill(clust_nstrip[k]);
          h_clust_maxstrip_u[1]->Fill(clust_maxstrip[k]);
          h_clust_maxsamp_u[1]->Fill(clust_maxsamp[k]);
          h_clust_maxadc_u[1]->Fill(clust_maxadc[k]);
          h_clust_adc_u[1]->Fill(clust_adc[k]);
        } else if (clust_module[k] == 1 && clust_axis[k] == 1 && cluster_id_1v == clust_index[k]) {
          h_clust_time_v[1]->Fill(clust_time[k]);
          h_clust_nstrip_v[1]->Fill(clust_nstrip[k]);
          h_clust_maxstrip_v[1]->Fill(clust_maxstrip[k]);
          h_clust_maxsamp_v[1]->Fill(clust_maxsamp[k]);
          h_clust_maxadc_v[1]->Fill(clust_maxadc[k]);
          h_clust_adc_v[1]->Fill(clust_adc[k]);
        }

      }
    }
    
    
    // Print the status as a percentage
    if (i % (nEntries / 100) == 0) {
      std::cout << "\rProcessing: " << int(i * 100.0 / nEntries) << "% completed." << std::flush;
    }

  } // End Event Loop
  std::cout << "\nProcessing completed." << std::endl;

  // Create the output file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
    return;
  }
  // Write the histograms and canvases to the output file
  outputFile->cd();
  h_d0->Write();
  h_projy->Write();
  h_projz->Write();
  for (int i = 0; i < 2; ++i) {
    h_clust_time_u[i]->Write();
    h_clust_time_v[i]->Write();
    h_clust_nstrip_u[i]->Write();
    h_clust_nstrip_v[i]->Write();
    h_clust_maxstrip_u[i]->Write();
    h_clust_maxstrip_v[i]->Write();
    h_clust_maxsamp_u[i]->Write();
    h_clust_maxsamp_v[i]->Write();
    h_clust_maxadc_u[i]->Write();
    h_clust_maxadc_v[i]->Write();
    h_clust_adc_u[i]->Write();
    h_clust_adc_v[i]->Write();
  }

  // Create directories in the output file
  
  // Create histograms for all planes added together
  // Close the files
  file->Close();
  outputFile->Close();
  std::cout << "Output file created: " << outputFileName << std::endl;
  return;
}