#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <tuple>
#include <cmath>
#include <TMath.h>
#include <string>
#include <TH2D.h>
#include <TH1D.h>

using namespace std;

void check_gem_positions(string energy= "400") {
  // Open the ROOT file
  string filepath     = "/volatile/hallc/c-lad/ehingerl/G4_LAD/carlos_proton/";
  string infile_name  = filepath + "raw/" + "ScanLAD_proton_" + energy + "MeV_10k_20240205.root";
  string outfile_name = "../hists/GEM_2D_position_plots.root";

  TFile *file = TFile::Open(infile_name.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file" << std::endl;
    return;
  }

  // Get the tree from the file
  TTree *tree;
  file->GetObject("gemana", tree);
  if (!tree) {
    std::cerr << "Error getting tree" << std::endl;
    file->Close();
    return;
  }

  // Define vectors to hold the branch data
  vector<double> *vXglo = nullptr;
  vector<double> *vYglo = nullptr;
  vector<double> *vZglo = nullptr;

  // Set branch addresses
  tree->SetBranchAddress("vXglo", &vXglo);
  tree->SetBranchAddress("vYglo", &vYglo);
  tree->SetBranchAddress("vZglo", &vZglo);

  // Create histograms
  TH2D *hXZ = new TH2D("hXZ", "X vs Z", 100, 500, 1500, 100, -1000, 0);
  TH1D *hY = new TH1D("hY", "Y", 100, -70, 70);
  TH2D *h_rXZ = new TH2D("h_rXZ", "X vs rXZ", 100, 0, 1500, 100, 0, 2000);


  // Loop over the entries in the tree
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    // Loop over the vectors and fill the histograms
    for (size_t j = 0; j < vXglo->size(); ++j) {
      hXZ->Fill(vXglo->at(j), vZglo->at(j));
      hY->Fill(vYglo->at(j));
      h_rXZ->Fill(vXglo->at(j), sqrt(vXglo->at(j)*vXglo->at(j) + vZglo->at(j)*vZglo->at(j)));
    }
  }



  // Save histograms to a new ROOT file
  TFile *outfile = new TFile(outfile_name.c_str(), "RECREATE");
  hXZ->Write();
  hY->Write();
  h_rXZ->Write();
  outfile->Close();

  // Clean up
  file->Close();
}