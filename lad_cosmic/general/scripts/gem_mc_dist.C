#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TString.h>
#include <TTree.h>
#include <string>

void gem_mc_dist() {
  std::string filename = "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/ROOTfiles/gem_mc.root";
  // Open the ROOT file
  TFile *file = TFile::Open(filename.c_str());
  if (!file || file->IsZombie()) {
    printf("Error opening file %s\n", filename.c_str());
    return;
  }

  // Get the tree T
  TTree *tree = (TTree *)file->Get("T");
  if (!tree) {
    printf("Error getting tree T from file %s\n", filename.c_str());
    file->Close();
    return;
  }

  const double gem_width  = 0.4*3072/1000*4;
  const double gem_height = 0.4*1536/1000*4; // 0.4mm spacing * n strips in m
  // Set branch addresses
  const int maxNdata = 1000; // Adjust this value as needed
  float maxpos[maxNdata];
  double axis[maxNdata];
  double moduleArr[maxNdata];
  int Ndata;

  tree->SetBranchAddress("L.gem.clust.maxpos", maxpos);
  tree->SetBranchAddress("L.gem.clust.axis", axis);
  tree->SetBranchAddress("L.gem.clust.module", moduleArr);
  tree->SetBranchAddress("Ndata.L.gem.clust.module", &Ndata);

  // Open the output file
  TFile outFile("../histos/gem_maxpos_plots.root", "RECREATE");

  // Loop over modules
  for (int module = 0; module < 2; module++) { // Adjust the range as needed
    // Create a histogram for the current module
    TH2D *hist = new TH2D(Form("hist_module_%d", module), Form("GEM Maxpos Plot Module %d;X;Y", module), 100, -gem_width , gem_width, 100, -gem_height,gem_height);

    // Loop over entries in the tree
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
      tree->GetEntry(i);
      for (int j = 0; j < Ndata; j++) {
        if (moduleArr[j] == module) {
          hist->Fill(maxpos[j][0], maxpos[j][1]);
        }
      }
    }

    // Write the histogram to the output file
    hist->Write();

    // Clean up
    delete hist;
  }

  // Close the output file
  outFile.Close();

  // Clean up
  file->Close();
}