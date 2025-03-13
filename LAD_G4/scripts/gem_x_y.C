
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2.h>
#include <TTree.h>
#include <iostream>
#include <vector>

using namespace std;

const int nXStrips = 128 * 24 * 2;
const int nYStrips = 128 * 12 * 2;

void gem_x_y() {
  // string input_file =
  // "/lustre24/expphy/volatile/hallc/c-lad/ehingerl/G4_LAD/carlos_proton/dig/ScanLAD_proton_400MeV_10k_20240205_dig.root";
  // string input_file  = "/u/home/ehingerl/hallc/software/libLADdig/test_scripts/lad_hodo_gem_sim.root";
  string input_file  = "/volatile/hallc/c-lad/ehingerl/G4_LAD/carlos_proton/dig/ScanLAD_proton_400MeV_10k_20240205_dig_test.root";
  string output_file = "/work/hallc/c-lad/ehingerl/analysis/LAD_G4/hists/GEM_2D_plots_fake_g4_400MeV_test.root";

  TFile *file = TFile::Open(input_file.c_str());
  if (!file || file->IsZombie()) {
    cerr << "Error opening file: " << input_file << endl;
    return;
  }

  TTree *tree = (TTree *)file->Get("T");
  if (!tree) {
    cerr << "Error: Tree 'T' not found in file: " << input_file << endl;
    file->Close();
    return;
  }

  vector<int> *module = nullptr;
  vector<int> *strip  = nullptr;

  tree->SetBranchAddress("LAD.GEM.dighit.module", &module);
  tree->SetBranchAddress("LAD.GEM.dighit.strip", &strip);

  TH2F *hist2D_01 = new TH2F("hist2D_01", "2D Histogram for Modules 0 and 1;X Strip;Y Strip", nXStrips, 0, nXStrips,
                             nYStrips, 0, nYStrips);
  TH2F *hist2D_23 = new TH2F("hist2D_23", "2D Histogram for Modules 2 and 3;X Strip;Y Strip", nXStrips, 0, nXStrips,
                             nYStrips, 0, nYStrips);

  for (Long64_t i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    for (size_t j = 0; j < module->size(); j++) {
      if (module->at(j) == 0) {
        for (size_t k = 0; k < module->size(); k++) {
          if (module->at(k) == 1) {
            hist2D_01->Fill(strip->at(j), strip->at(k));
          }
        }
      } else if (module->at(j) == 2) {
        for (size_t k = 0; k < module->size(); k++) {
          if (module->at(k) == 3) {
            hist2D_23->Fill(strip->at(j), strip->at(k));
          }
        }
      }
    }
  }

  TFile *output = new TFile(output_file.c_str(), "RECREATE");
  hist2D_01->Write();
  hist2D_23->Write();
  output->Close();

  file->Close();
}