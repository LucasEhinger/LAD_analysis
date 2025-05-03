#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <iostream>

void good_track_info() {
  gROOT->SetBatch(kTRUE);
  // Open the ROOT file
  const char *filename = "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/CALIB/LAD_COIN_22281_-1.root"; // Replace with your file name
  TFile *file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open the file!" << std::endl;
    return;
  }

  // Get the tree "T"
  TTree *tree = (TTree*)file->Get("T");
  if (!tree) {
    std::cerr << "Error: Cannot find the tree 'T' in the file!" << std::endl;
    file->Close();
    return;
  }

  // Define variables to hold branch data
  const int maxSize = 1000; // Adjust size as needed
  double goodhit_deltapostrans_0[maxSize];
  double goodhit_deltapostrans_1[maxSize];
  double goodhit_plane_0[maxSize];
  double goodhit_plane_1[maxSize];
  double goodhit_paddle_0[maxSize];
  double goodhit_paddle_1[maxSize];
  double goodhit_track_id_0[maxSize];
  double goodhit_track_id_1[maxSize];
  int Ndata_0;
  int Ndata_1;

  // Set branch addresses
  tree->SetBranchAddress("H.ladkin.goodhit_deltapostrans_0", goodhit_deltapostrans_0);
  tree->SetBranchAddress("H.ladkin.goodhit_deltapostrans_1", goodhit_deltapostrans_1);
  tree->SetBranchAddress("H.ladkin.goodhit_plane_0", goodhit_plane_0);
  tree->SetBranchAddress("H.ladkin.goodhit_plane_1", goodhit_plane_1);
  tree->SetBranchAddress("H.ladkin.goodhit_paddle_0", goodhit_paddle_0);
  tree->SetBranchAddress("H.ladkin.goodhit_paddle_1", goodhit_paddle_1);
  tree->SetBranchAddress("H.ladkin.goodhit_trackid_0", goodhit_track_id_0);
  tree->SetBranchAddress("H.ladkin.goodhit_trackid_1", goodhit_track_id_1);
  tree->SetBranchAddress("Ndata.H.ladkin.goodhit_deltapostrans_0", &Ndata_0);
  tree->SetBranchAddress("Ndata.H.ladkin.goodhit_deltapostrans_1", &Ndata_1);

  // Create histograms for paddle numbers for each plane
  TH1D *hist_plane_0[3];
  TH1D *hist_plane_1[2];

  for (int i = 0; i < 3; ++i) {
    hist_plane_0[i] = new TH1D(Form("hist_plane_0_%d", i), Form("Paddle Numbers for Plane 0 - Histogram %d;Paddle Number;Counts", i), 11, -0.5, 10.5);
  }

  for (int i = 0; i < 2; ++i) {
    hist_plane_1[i] = new TH1D(Form("hist_plane_1_%d", i), Form("Paddle Numbers for Plane 1 - Histogram %d;Paddle Number;Counts", i), 11, -0.5, 10.5);
  }

  // Define the range for deltapostrans
  const double deltapostrans_min = -100.0; // Adjust as needed
  const double deltapostrans_max = 100.0;  // Adjust as needed

  // Loop through all entries
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    // Loop through data for each event
    for (int j = 0; j < Ndata_0; ++j) {
      // Apply the cut on deltapostrans for plane 0
      if (goodhit_deltapostrans_0[j] >= deltapostrans_min && goodhit_deltapostrans_0[j] <= deltapostrans_max) {
        hist_plane_0[(int)goodhit_plane_0[j]/2]->Fill(goodhit_paddle_0[j]);
      }
    }
    for (int j = 0; j < Ndata_1; ++j) {
      // Apply the cut on deltapostrans for plane 1
      if (goodhit_deltapostrans_1[j] >= deltapostrans_min && goodhit_deltapostrans_1[j] <= deltapostrans_max) {
        hist_plane_1[(int)goodhit_plane_1[j]/2]->Fill(goodhit_paddle_1[j]);
      }
    }
    // Print status percentage
    if (i % (nEntries / 100) == 0 && nEntries > 100) {
      std::cout << "\rProcessing: " << (i * 100 / nEntries) << "% completed." << std::flush;
    }
  }

  // Draw the histograms
  TCanvas *canvas = new TCanvas("canvas", "Paddle Number Histograms", 800, 600);
  canvas->Divide(3, 2); // Adjusted to 3 columns and 2 rows to fit all histograms

  for (int i = 0; i < 3; ++i) {
    canvas->cd(i + 1);
    hist_plane_0[i]->Draw();
  }

  for (int i = 0; i < 2; ++i) {
    canvas->cd(i + 4);
    hist_plane_1[i]->Draw();
  }

  // Save the canvas to a PDF file
  canvas->SaveAs("good_track_info.pdf");

  // Clean up
  file->Close();
  delete file;
}