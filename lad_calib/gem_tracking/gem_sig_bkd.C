#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TTree.h>
#include <iostream>
#include <vector>

const int MAX_DATA     = 10000;
const int d0_NBINS     = 100;
const double d0_MIN    = 0.0;
const double d0_MAX    = 100.0;
const int projz_NBINS  = 100;
const double projz_MIN = -20.0;
const double projz_MAX = 20.0;
const int dTrk_NBINS   = 100;
const double dTrk_MIN  = -500.0;
const double dTrk_MAX  = 500.0;
const int time_NBINS   = 100;
const double time_MIN  = 1700.0;
const double time_MAX  = 2000.0;
const double theta_MIN = 0.0;
const double theta_MAX = 180.0;
const int theta_NBINS  = 180;

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

const double delta_pos_trans_sig = 80.0;
const double delta_pos_long_sig  = 80.0;
const double d0_cut_sig          = 30.0;
const double delta_pos_trans_bkd = 400.0;
const double delta_pos_long_bkd  = 400.0;
const double d0_cut_bkd          = 25.0;

const int maxTracks = 30;

void run_gem_sig_bkd(int run_number, int num_replayed) {
  // Set batch mode to suppress graphical output
  gROOT->SetBatch(kTRUE);
  // Open the ROOT file
  TString fileName = Form("/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_%d_%d.root",
                          run_number, num_replayed);
  TString outputFileName = Form("files/sig_bkd/gem_sig_bkd_plots_%d_%d.root", run_number, num_replayed);
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
  Double_t trk_projz[MAX_DATA];
  Double_t kin_trackID_0[MAX_DATA], kin_trackID_1[MAX_DATA];
  Double_t kin_plane_0[MAX_DATA], kin_plane_1[MAX_DATA];
  Double_t kin_paddle_0[MAX_DATA], kin_paddle_1[MAX_DATA];
  Double_t kin_hittime_0[MAX_DATA], kin_hittime_1[MAX_DATA];
  Double_t kin_hittheta_0[MAX_DATA], kin_hittheta_1[MAX_DATA];
  Double_t kin_hitphi_0[MAX_DATA], kin_hitphi_1[MAX_DATA];
  Double_t kin_hitedep_0[MAX_DATA], kin_hitedep_1[MAX_DATA];
  Double_t kin_dTrkhoriz_0[MAX_DATA], kin_dTrkhoriz_1[MAX_DATA];
  Double_t kin_dTrkVert_0[MAX_DATA], kin_dTrkVert_1[MAX_DATA];
  Int_t nTracks, nGoodHits;

  T->SetBranchAddress("Ndata.H.gem.trk.d0", &nTracks);
  T->SetBranchAddress("Ndata.H.ladkin.goodhit_trackid_0", &nGoodHits);
  T->SetBranchAddress("H.gem.trk.d0", &trk_d0);
  T->SetBranchAddress("H.gem.trk.d0_good", &trk_d0_good);
  T->SetBranchAddress("H.gem.trk.projz", &trk_projz);

  std::string branchPrefix = "H.ladhod"; // Change this to "H.ladhod" if needed

  T->SetBranchAddress((branchPrefix + ".goodhit_trackid_0").c_str(), &kin_trackID_0);
  T->SetBranchAddress((branchPrefix + ".goodhit_trackid_1").c_str(), &kin_trackID_1);
  T->SetBranchAddress((branchPrefix + ".goodhit_plane_0").c_str(), &kin_plane_0);
  T->SetBranchAddress((branchPrefix + ".goodhit_plane_1").c_str(), &kin_plane_1);
  T->SetBranchAddress((branchPrefix + ".goodhit_paddle_0").c_str(), &kin_paddle_0);
  T->SetBranchAddress((branchPrefix + ".goodhit_paddle_1").c_str(), &kin_paddle_1);
  T->SetBranchAddress((branchPrefix + ".goodhit_hittime_0").c_str(), &kin_hittime_0);
  T->SetBranchAddress((branchPrefix + ".goodhit_hittime_1").c_str(), &kin_hittime_1);
  T->SetBranchAddress((branchPrefix + ".goodhit_hittheta_0").c_str(), &kin_hittheta_0);
  T->SetBranchAddress((branchPrefix + ".goodhit_hittheta_1").c_str(), &kin_hittheta_1);
  T->SetBranchAddress((branchPrefix + ".goodhit_hitphi_0").c_str(), &kin_hitphi_0);
  T->SetBranchAddress((branchPrefix + ".goodhit_hitphi_1").c_str(), &kin_hitphi_1);
  T->SetBranchAddress((branchPrefix + ".goodhit_hitedep_0").c_str(), &kin_hitedep_0);
  T->SetBranchAddress((branchPrefix + ".goodhit_hitedep_1").c_str(), &kin_hitedep_1);
  T->SetBranchAddress((branchPrefix + ".goodhit_deltapostrans_0").c_str(), &kin_dTrkhoriz_0);
  T->SetBranchAddress((branchPrefix + ".goodhit_deltapostrans_1").c_str(), &kin_dTrkhoriz_1);
  T->SetBranchAddress((branchPrefix + ".goodhit_deltaposlong_0").c_str(), &kin_dTrkVert_0);
  T->SetBranchAddress((branchPrefix + ".goodhit_deltaposlong_1").c_str(), &kin_dTrkVert_1);

  // Create TH2D histograms for delta_pos_long vs delta_pos_trans
  TH2D *h_dTrkVert_vs_dTrkHoriz_0 =
      new TH2D("h_dTrkVert_vs_dTrkHoriz_0", "dTrkVert vs dTrkHoriz (0);dTrkHoriz (cm);dTrkVert (cm)", dTrk_NBINS,
               dTrk_MIN, dTrk_MAX, dTrk_NBINS, dTrk_MIN, dTrk_MAX);

  TH2D *h_dTrkVert_vs_dTrkHoriz_1 =
      new TH2D("h_dTrkVert_vs_dTrkHoriz_1", "dTrkVert vs dTrkHoriz (1);dTrkHoriz (cm);dTrkVert (cm)", dTrk_NBINS,
               dTrk_MIN, dTrk_MAX, dTrk_NBINS, dTrk_MIN, dTrk_MAX);

  // Create TH2D histograms for delta_pos_long vs delta_pos_trans with d0 cut
  TH2D *h_dTrkVert_vs_dTrkHoriz_d0Cut_0 =
      new TH2D("h_dTrkVert_vs_dTrkHoriz_d0Cut_0", "dTrkVert vs dTrkHoriz with d0 Cut (0);dTrkHoriz (cm);dTrkVert (cm)",
               dTrk_NBINS, dTrk_MIN, dTrk_MAX, dTrk_NBINS, dTrk_MIN, dTrk_MAX);

  TH2D *h_dTrkVert_vs_dTrkHoriz_d0Cut_1 =
      new TH2D("h_dTrkVert_vs_dTrkHoriz_d0Cut_1", "dTrkVert vs dTrkHoriz with d0 Cut (1);dTrkHoriz (cm);dTrkVert (cm)",
               dTrk_NBINS, dTrk_MIN, dTrk_MAX, dTrk_NBINS, dTrk_MIN, dTrk_MAX);

  // Create TH2D histograms for delta_pos_long vs delta_pos_trans with track cut
  TH2D *h_dTrkVert_vs_dTrkHoriz_trkCut_0 =
      new TH2D("h_dTrkVert_vs_dTrkHoriz_trkCut_0", "dTrkVert vs dTrkHoriz with track Cut (0);dTrkHoriz (cm);dTrkVert (cm)", dTrk_NBINS,
               dTrk_MIN, dTrk_MAX, dTrk_NBINS, dTrk_MIN, dTrk_MAX);
  TH2D *h_dTrkVert_vs_dTrkHoriz_trkCut_1 =
      new TH2D("h_dTrkVert_vs_dTrkHoriz_trkCut_1", "dTrkVert vs dTrkHoriz with track Cut (1);dTrkHoriz (cm);dTrkVert (cm)", dTrk_NBINS,
               dTrk_MIN, dTrk_MAX, dTrk_NBINS, dTrk_MIN, dTrk_MAX);

  // Create TH2D histograms for delta_pos_long vs delta_pos_trans with D0 cut
  TH2D *h_dTrkVert_vs_dTrkHoriz_trkCut_D0Cut_0 = new TH2D(
      "h_dTrkVert_vs_dTrkHoriz_trkCut_D0Cut_0", "dTrkVert vs dTrkHoriz with track & d0 Cut (0);dTrkHoriz (cm);dTrkVert (cm)",
      dTrk_NBINS, dTrk_MIN, dTrk_MAX, dTrk_NBINS, dTrk_MIN, dTrk_MAX);
  TH2D *h_dTrkVert_vs_dTrkHoriz_trkCut_D0Cut_1 = new TH2D(
      "h_dTrkVert_vs_dTrkHoriz_trkCut_D0Cut_1", "dTrkVert vs dTrkHoriz with track & d0 Cut (1);dTrkHoriz (cm);dTrkVert (cm)",
      dTrk_NBINS, dTrk_MIN, dTrk_MAX, dTrk_NBINS, dTrk_MIN, dTrk_MAX);

  int n_signal             = 0;
  int n_signal_trk_cut     = 0;
  int n_background         = 0;
  int n_background_trk_cut = 0;
  // Loop over entries and fill arrays
  Long64_t nEntries = T->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    T->GetEntry(i);

    for (int k = 0; k < nGoodHits; ++k) {
      h_dTrkVert_vs_dTrkHoriz_0->Fill(kin_dTrkhoriz_0[k], kin_dTrkVert_0[k]);
      h_dTrkVert_vs_dTrkHoriz_1->Fill(kin_dTrkhoriz_1[k], kin_dTrkVert_1[k]);
      if (nTracks < maxTracks) {
        h_dTrkVert_vs_dTrkHoriz_trkCut_0->Fill(kin_dTrkhoriz_0[k], kin_dTrkVert_0[k]);
        h_dTrkVert_vs_dTrkHoriz_trkCut_1->Fill(kin_dTrkhoriz_1[k], kin_dTrkVert_1[k]);
      }
      if (abs(kin_dTrkVert_0[k]) < delta_pos_long_sig && abs(kin_dTrkhoriz_0[k]) < delta_pos_trans_sig &&
          trk_d0[int(kin_trackID_0[k])] < d0_cut_sig) {
        n_signal++;
        h_dTrkVert_vs_dTrkHoriz_d0Cut_0->Fill(kin_dTrkhoriz_0[k], kin_dTrkVert_0[k]);
        if (nTracks < maxTracks) {
          n_signal_trk_cut++;
          h_dTrkVert_vs_dTrkHoriz_trkCut_D0Cut_0->Fill(kin_dTrkhoriz_0[k], kin_dTrkVert_0[k]);
        }
      }
      if (abs(kin_dTrkVert_0[k]) < delta_pos_long_bkd && abs(kin_dTrkhoriz_0[k]) < delta_pos_trans_bkd &&
          trk_d0[int(kin_trackID_0[k])] < d0_cut_bkd) {
        n_background++;
        h_dTrkVert_vs_dTrkHoriz_d0Cut_0->Fill(kin_dTrkhoriz_0[k], kin_dTrkVert_0[k]);
        if (nTracks < maxTracks) {
          n_background_trk_cut++;
          h_dTrkVert_vs_dTrkHoriz_trkCut_D0Cut_0->Fill(kin_dTrkhoriz_0[k], kin_dTrkVert_0[k]);
        }
      }
      // if (abs(kin_dTrkVert_1[k]) < delta_pos_long_sig && abs(kin_dTrkhoriz_1[k]) < delta_pos_trans_sig) {
      //   n_signal++;
      //   if (nTracks < maxTracks) {
      //     n_signal_trk_cut++;
      //   }
      // }
      // if (abs(kin_dTrkVert_1[k]) < delta_pos_long_bkd && abs(kin_dTrkhoriz_1[k]) < delta_pos_trans_bkd) {
      //   n_background++;
      //   if (nTracks < maxTracks) {
      //     n_background_trk_cut++;
      //   }
      // }
    }
    // Print the status as a percentage
    int progress = static_cast<int>((100.0 * i) / nEntries);
    std::cout << "\rProcessing: " << progress << "% completed" << std::flush;
  }

  double bdk_norm =
      n_background / (delta_pos_trans_bkd * delta_pos_long_bkd) * (delta_pos_trans_sig * delta_pos_long_sig);
  double sig_norm = n_signal - bdk_norm;

  double bdk_norm_trk_cut =
      n_background_trk_cut / (delta_pos_trans_bkd * delta_pos_long_bkd) * (delta_pos_trans_sig * delta_pos_long_sig);
  double sig_norm_trk_cut = n_signal_trk_cut - bdk_norm_trk_cut;

  std::cout << "\nProcessing completed." << std::endl;
  // End Event Loop

  // Create output file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
    return;
  }

  // Write histograms to the output file
  outputFile->cd();
  h_dTrkVert_vs_dTrkHoriz_0->Write();
  // h_dTrkVert_vs_dTrkHoriz_1->Write();
  h_dTrkVert_vs_dTrkHoriz_trkCut_0->Write();
  // h_dTrkVert_vs_dTrkHoriz_trkCut_1->Write();
  h_dTrkVert_vs_dTrkHoriz_d0Cut_0->Write();
  // h_dTrkVert_vs_dTrkHoriz_d0Cut_1->Write();
  h_dTrkVert_vs_dTrkHoriz_trkCut_D0Cut_0->Write();
  // h_dTrkVert_vs_dTrkHoriz_trkCut_D0Cut_1->Write();

  // Draw red and green boxes on the histograms
  TCanvas *c_boxes = new TCanvas("c_boxes", "Boxes on Histograms", 800, 600);
  c_boxes->Divide(2, 2);

  c_boxes->cd(1);
  h_dTrkVert_vs_dTrkHoriz_0->Draw("COLZ");
  TBox *redBox_0 = new TBox(-delta_pos_trans_sig, -delta_pos_long_sig, delta_pos_trans_sig, delta_pos_long_sig);
  redBox_0->SetLineColor(kRed);
  redBox_0->SetLineWidth(2);
  redBox_0->SetFillStyle(0);
  redBox_0->Draw("same");

  TBox *greenBox_0 = new TBox(-delta_pos_trans_bkd, -delta_pos_long_bkd, delta_pos_trans_bkd, delta_pos_long_bkd);
  greenBox_0->SetLineColor(kGreen);
  greenBox_0->SetLineWidth(2);
  greenBox_0->SetFillStyle(0);
  greenBox_0->Draw("same");

  // c_boxes->cd(2);
  // h_deltaposlong_vs_deltapostrans_1->Draw("COLZ");
  // TBox *redBox_1 = new TBox(-delta_pos_trans_sig, -delta_pos_long_sig, delta_pos_trans_sig, delta_pos_long_sig);
  // redBox_1->SetLineColor(kRed);
  // redBox_1->SetLineWidth(2);
  // redBox_1->SetFillStyle(0);
  // redBox_1->Draw("same");

  // TBox *greenBox_1 = new TBox(-delta_pos_trans_bkd, -delta_pos_long_bkd, delta_pos_trans_bkd, delta_pos_long_bkd);
  // greenBox_1->SetLineColor(kGreen);
  // greenBox_1->SetLineWidth(2);
  // greenBox_1->SetFillStyle(0);
  // greenBox_1->Draw("same");

  c_boxes->cd(2);
  h_dTrkVert_vs_dTrkHoriz_trkCut_0->Draw("COLZ");
  redBox_0->Draw("same");
  greenBox_0->Draw("same");

  // c_boxes->cd(4);
  // h_deltaposlong_vs_deltapostrans_trackcut_1->Draw("COLZ");
  // redBox_1->Draw("same");
  // greenBox_1->Draw("same");

  // Draw the boxes on the second histogram
  c_boxes->cd(3);
  h_dTrkVert_vs_dTrkHoriz_d0Cut_0->Draw("COLZ");
  redBox_0->Draw("same");
  greenBox_0->Draw("same");
  c_boxes->cd(4);
  h_dTrkVert_vs_dTrkHoriz_trkCut_D0Cut_0->Draw("COLZ");
  redBox_0->Draw("same");
  greenBox_0->Draw("same");

  c_boxes->Write();

  // Write the number of signal and background events to the output file
  // Create histograms for signal, background, and signal-to-background ratio
  // Create a canvas to display the text
  TCanvas *c_signal_background = new TCanvas("c_signal_background", "Signal and Background", 800, 600);

  // Draw text for signal count, background count, and signal-to-background ratio
  c_signal_background->cd();
  TLatex latex;
  latex.SetTextSize(0.04);
  latex.DrawLatex(0.1, 0.8, Form("Bkd Sub Signal Count (No Track Cut): %0.1f", sig_norm));
  latex.DrawLatex(0.1, 0.75, Form("Normalized Background Count (No Track Cut): %0.1f", bdk_norm));
  if (n_background > 0) {
    latex.DrawLatex(0.1, 0.7, Form("Signal to Background Ratio (No Track Cut): %.2f", sig_norm / bdk_norm));
  } else {
    latex.DrawLatex(0.1, 0.7, "Signal to Background Ratio (No Track Cut): Undefined (Background = 0)");
  }

  latex.DrawLatex(0.1, 0.6, Form("Bkd Sub Signal Count (With Track Cut): %0.1f", sig_norm_trk_cut));
  latex.DrawLatex(0.1, 0.55, Form("Normalized Background Count (With Track Cut): %0.1f", bdk_norm_trk_cut));
  if (n_background_trk_cut > 0) {
    latex.DrawLatex(0.1, 0.5,
                    Form("Signal to Background Ratio (With Track Cut): %.2f", sig_norm_trk_cut / bdk_norm_trk_cut));
  } else {
    latex.DrawLatex(0.1, 0.5, "Signal to Background Ratio (With Track Cut): Undefined (Background = 0)");
  }

  // Write the canvas to the output file
  c_signal_background->Write();

  delete c_signal_background;
  delete c_boxes;
  // Clean up
  outputFile->Close();
  delete outputFile;
  file->Close();
  delete file;
}

void gem_sig_bkd() {
  int run_number = 22282;
  for (int num_replayed = 400000; num_replayed <= 400010; ++num_replayed) {
    run_gem_sig_bkd(run_number, num_replayed);
  }
}