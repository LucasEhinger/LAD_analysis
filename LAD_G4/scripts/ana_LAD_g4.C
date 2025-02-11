#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void ana_LAD_g4(string energy= "400") {
  // Open the ROOT file
  string filepath     = "/volatile/hallc/c-lad/ehingerl/G4_LAD/carlos_proton/";
  string infile_name  = filepath + "raw/" + "ScanLAD_proton_" + energy + "MeV_10k_20240205.root";
  string outfile_name = filepath + "ana/" + "ScanLAD_proton_" + energy + "MeV_10k_20240205_ana.root";

  TFile *file = TFile::Open(infile_name.c_str(), "READ");
  if (!file || file->IsZombie()) {
    cerr << "Error opening file" << endl;
    return;
  }

  // Get the hodoposition tree
  TTree *hodoposition = (TTree *)file->Get("hodoposition");
  // Get the hodoenergy tree
  TTree *hodoenergy = (TTree *)file->Get("hodoenergy");
  if (!hodoposition || !hodoenergy) {
    cerr << "Error getting trees" << endl;
    file->Close();
    return;
  }

  // Set branch addresses for hodoposition
  vector<double> *vXbar = nullptr;
  vector<double> *vYbar = nullptr;
  vector<double> *vZbar = nullptr;
  vector<double> *vTbar = nullptr;
  vector<int> *vPaddle  = nullptr;
  vector<int> *vTrackID = nullptr;

  hodoposition->SetBranchAddress("vXbar", &vXbar);
  hodoposition->SetBranchAddress("vYbar", &vYbar);
  hodoposition->SetBranchAddress("vZbar", &vZbar);
  hodoposition->SetBranchAddress("vTbar", &vTbar);
  hodoposition->SetBranchAddress("vPaddle", &vPaddle);
  hodoposition->SetBranchAddress("vTrackID", &vTrackID);

  // Set branch addresses for hodoenergy
  int EventID;
  vector<int> *vPadCopy   = nullptr;
  vector<double> *vEneDep = nullptr;
  vector<int> *vPDG       = nullptr;
  vector<int> *vLevel     = nullptr;

  hodoenergy->SetBranchAddress("EventID", &EventID);
  hodoenergy->SetBranchAddress("vPadCopy", &vPadCopy);
  hodoenergy->SetBranchAddress("vEneDep", &vEneDep);
  hodoenergy->SetBranchAddress("vPDG", &vPDG);
  hodoenergy->SetBranchAddress("vLevel", &vLevel);

  TFile *outfile = TFile::Open(outfile_name.c_str(), "RECREATE");
  if (!outfile || outfile->IsZombie()) {
    std::cerr << "Error creating output file" << std::endl;
    return;
  }
  // Create new tree with only the necessary branches
  TTree *outTree                = new TTree("T", "T");
  vector<double> *out_vXbar     = new vector<double>();
  vector<double> *out_vYbar     = new vector<double>();
  vector<double> *out_vZbar     = new vector<double>();
  vector<double> *out_vTbar     = new vector<double>();
  vector<int> *out_vPaddle      = new vector<int>();
  vector<double> *out_vEDep     = new vector<double>();
  vector<double> *out_vEDep_all = new vector<double>();

  outTree->Branch("vXbar", &out_vXbar);
  outTree->Branch("vYbar", &out_vYbar);
  outTree->Branch("vZbar", &out_vZbar);
  outTree->Branch("vTbar", &out_vTbar);
  outTree->Branch("vPaddle", &out_vPaddle);
  outTree->Branch("vEDep", &out_vEDep);
  outTree->Branch("vEDep_all", &out_vEDep_all);

  // Loop over entries and print the vectors
  Long64_t nEntries = hodoposition->GetEntries();
  // nEntries          = 15; // For testing purposes
  for (Long64_t i = 0; i < nEntries; ++i) {
    out_vXbar->clear();
    out_vYbar->clear();
    out_vZbar->clear();
    out_vTbar->clear();
    out_vPaddle->clear();
    out_vEDep->clear();
    out_vEDep_all->clear();
    hodoposition->GetEntry(i);
    hodoenergy->GetEntry(i);

    // Loop over the hits
    for (int j = 0; j < vXbar->size(); ++j) {
      if (vTrackID->at(j) != 1) {
        continue;
      }
      out_vXbar->push_back(vXbar->at(j));
      out_vYbar->push_back(vYbar->at(j));
      out_vZbar->push_back(vZbar->at(j));
      out_vTbar->push_back(vTbar->at(j));
      out_vPaddle->push_back(vPaddle->at(j));
    }
    for (int j = 0; j < vPadCopy->size(); ++j) {
      int indx = -1;
      for (int k = 0; k < out_vPaddle->size(); ++k) {
        if (out_vPaddle->at(k) == vPadCopy->at(j)) {
          indx = k;
          break;
        }
      }
      if (indx == -1) {
        continue;
      }

      if (vPDG->at(j) == 2212 && vLevel->at(j) == 1) {
        if (out_vEDep->size() <= indx) {
          out_vEDep->resize(indx + 1, 0.0);
        }
        (*out_vEDep)[indx] += vEneDep->at(j);
      }

      if (out_vEDep_all->size() <= indx) {
        out_vEDep_all->resize(indx + 1, 0.0);
      }
      (*out_vEDep_all)[indx] += vEneDep->at(j);
    }

    outTree->Fill();
  }

  // Write the new tree to a new file

  outfile->cd();
  outTree->Write();
  // Close the files
  outfile->Close();
  file->Close();
}