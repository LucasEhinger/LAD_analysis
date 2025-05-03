#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TTree.h>

const static int MAXDATA = 1000;
const static int nTracksCut = 30; // Set the cut for the number of tracks
const static int Z_NBINS = 40;
const static double Z_MIN = -20.0;
const static double Z_MAX = 20.0;

void gem_track_vertex() {
  gROOT->SetBatch(kTRUE);
  const char *filename =
      "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/CALIB/LAD_COIN_22281_6002.root";
  // Open the ROOT file
  TFile *file = TFile::Open(filename);
  if (!file || file->IsZombie()) {
    printf("Error: Cannot open file %s\n", filename);
    return;
  }

  // Access the tree
  TTree *tree = (TTree *)file->Get("T");
  if (!tree) {
    printf("Error: Cannot find tree T in file %s\n", filename);
    file->Close();
    return;
  }

  // Variables to hold data
  Double_t react_z;
  Double_t ntracks;
  Double_t projz[MAXDATA] = {0};
  int ndata_projz         = 0;

  // Set branch addresses
  tree->SetBranchAddress("P.react.z", &react_z);
  tree->SetBranchAddress("Ndata.P.gem.trk.projz", &ndata_projz);
  tree->SetBranchAddress("P.gem.trk.projz", projz);
  tree->SetBranchAddress("P.gem.trk.ntracks", &ntracks);

  // Histograms for P.react.z and P.gem.trk.projz
  TH1F *h_react_z = new TH1F("h_react_z", "P.react.z;Z Position;Counts", 100, -20, 20);
  TH1F *h_projz   = new TH1F("h_projz", "P.gem.trk.projz;Projected Z;Counts", Z_NBINS, Z_MIN, Z_MAX);
  // Histograms for P.gem.trk.projz with react.z cuts
  TH1F *h_projz_cut1 =
      new TH1F("h_projz_cut1", "P.gem.trk.projz (-15 < react.z < -5);Projected Z;Counts", Z_NBINS, Z_MIN, Z_MAX);
  TH1F *h_projz_cut2 = new TH1F("h_projz_cut2", "P.gem.trk.projz (-5 < react.z < 5);Projected Z;Counts", Z_NBINS, Z_MIN, Z_MAX);
  TH1F *h_projz_cut3 = new TH1F("h_projz_cut3", "P.gem.trk.projz (5 < react.z < 15);Projected Z;Counts", Z_NBINS, Z_MIN, Z_MAX);

  // Loop over entries in the tree
  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
    // Load the entry
    tree->GetEntry(i);

    if(ntracks > nTracksCut)
      continue;
    // Fill histograms for all entries
    h_react_z->Fill(react_z);

    // Loop over tracks
    for (Int_t j = 0; j < ndata_projz; j++) {
      Float_t projz_val = projz[j];
      h_projz->Fill(projz_val);

      // Apply cuts on react.z and fill corresponding histograms
      if (react_z > -15 && react_z < -5) {
        h_projz_cut1->Fill(projz_val);
      } else if (react_z > -5 && react_z < 5) {
        h_projz_cut2->Fill(projz_val);
      } else if (react_z > 5 && react_z < 15) {
        h_projz_cut3->Fill(projz_val);
      }
    }
  }

  // Create canvas for histograms
  TCanvas *c1 = new TCanvas("c1", "Histograms", 800, 600);
  c1->Divide(2, 1);

  // Draw P.react.z and P.gem.trk.projz histograms
  c1->cd(1);
  h_react_z->Draw();

  c1->cd(2);
  h_projz->Draw();

  // Create canvas for overlaid plots
  TCanvas *c2 = new TCanvas("c2", "Overlaid Plots", 800, 600);
  h_projz_cut1->SetLineColor(kRed);
  h_projz_cut2->SetLineColor(kBlue);
  h_projz_cut3->SetLineColor(kGreen);

  h_projz_cut1->Draw();
  h_projz_cut2->Draw("SAME");
  h_projz_cut3->Draw("SAME");

  // Add legend
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(h_projz_cut1, "-15 < react.z < -5", "l");
  legend->AddEntry(h_projz_cut2, "-5 < react.z < 5", "l");
  legend->AddEntry(h_projz_cut3, "5 < react.z < 15", "l");
  legend->Draw();

  // Save all canvases into a single PDF
  TCanvas *c_combined = new TCanvas("c_combined", "Combined PDF", 800, 600);
  c_combined->Print("output.pdf["); // Open the PDF file
  c1->Print("output.pdf");
  c2->Print("output.pdf");
  c_combined->Print("output.pdf]"); // Close the PDF file

  // Clean up
  file->Close();
  delete file;
}
