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
  // Get the gemana tree
  TTree *gemana = (TTree *)file->Get("gemana");
  if (!hodoposition || !hodoenergy || !gemana) {
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

// Set branch addresses for gemana
  vector<double> *vXloc     = nullptr;
  vector<double> *vYloc     = nullptr;
  vector<double> *vZloc     = nullptr;
  vector<double> *vXglo     = nullptr;
  vector<double> *vYglo     = nullptr;
  vector<double> *vZglo     = nullptr;
  vector<double> *vTchamber = nullptr;
  vector<int> *vChamber     = nullptr;
  vector<int> *vgPDG        = nullptr;
  vector<int> *vgLevel      = nullptr;

  gemana->SetBranchAddress("vXloc", &vXloc);
  gemana->SetBranchAddress("vYloc", &vYloc);
  gemana->SetBranchAddress("vZloc", &vZloc);
  gemana->SetBranchAddress("vXglo", &vXglo);
  gemana->SetBranchAddress("vYglo", &vYglo);
  gemana->SetBranchAddress("vZglo", &vZglo);
  gemana->SetBranchAddress("vTchamber", &vTchamber);
  gemana->SetBranchAddress("vChamber", &vChamber);
  gemana->SetBranchAddress("vgPDG", &vgPDG);
  gemana->SetBranchAddress("vgLevel", &vgLevel);

  // Create a new ROOT file

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

  // Create branches for GEM hits
  int LAD_GEM_hit_nhits;
  vector<double> *LAD_GEM_hit_xin = new vector<double>();
  vector<double> *LAD_GEM_hit_yin = new vector<double>();
  vector<double> *LAD_GEM_hit_zin = new vector<double>();
  vector<double> *LAD_GEM_hit_xout = new vector<double>();
  vector<double> *LAD_GEM_hit_yout = new vector<double>();
  vector<double> *LAD_GEM_hit_zout = new vector<double>();
  vector<double> *LAD_GEM_hit_t = new vector<double>();
  vector<double> *LAD_GEM_hit_t_out = new vector<double>();
  vector<double> *LAD_GEM_hit_edep = new vector<double>();
  vector<int> *LAD_GEM_hit_plane = new vector<int>();

  outTree->Branch("LAD_GEM_hit_nhits", &LAD_GEM_hit_nhits);
  outTree->Branch("LAD_GEM_hit_xin", &LAD_GEM_hit_xin);
  outTree->Branch("LAD_GEM_hit_yin", &LAD_GEM_hit_yin);
  outTree->Branch("LAD_GEM_hit_zin", &LAD_GEM_hit_zin);
  outTree->Branch("LAD_GEM_hit_xout", &LAD_GEM_hit_xout);
  outTree->Branch("LAD_GEM_hit_yout", &LAD_GEM_hit_yout);
  outTree->Branch("LAD_GEM_hit_zout", &LAD_GEM_hit_zout);
  outTree->Branch("LAD_GEM_hit_t", &LAD_GEM_hit_t);
  outTree->Branch("LAD_GEM_hit_t_out", &LAD_GEM_hit_t_out);
  outTree->Branch("LAD_GEM_hit_edep", &LAD_GEM_hit_edep);
  outTree->Branch("LAD_GEM_hit_plane", &LAD_GEM_hit_plane);


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

    // Loop over the gem hits
    for (int j = 0; j < vXloc->size(); ++j) {

      if (vgLevel->at(j) != 1 || vgPDG->at(j) != 2212) { // 2212 is the PDG ID for protons
        continue;
      }
      if(j != 0 && vChamber->at(j) == vChamber->at(j-1) && vgPDG->at(j) == vgPDG->at(j-1) && vgLevel->at(j) == vgLevel->at(j-1)) {
        continue;
      }
      LAD_GEM_hit_xin->push_back(vXloc->at(j) / 1000);
      LAD_GEM_hit_yin->push_back(vYloc->at(j) / 1000);
      LAD_GEM_hit_zin->push_back(vZloc->at(j) / 1000);
      LAD_GEM_hit_xout->push_back(vXloc->at(j) / 1000);
      LAD_GEM_hit_yout->push_back(vYloc->at(j) / 1000);
      LAD_GEM_hit_zout->push_back((vZloc->at(j) + 5.0) / 1000);
      LAD_GEM_hit_t->push_back(vTchamber->at(j));
      LAD_GEM_hit_t_out->push_back(vTchamber->at(j));
      LAD_GEM_hit_edep->push_back(0);
      LAD_GEM_hit_plane->push_back(vChamber->at(j));
    }
    LAD_GEM_hit_nhits = LAD_GEM_hit_xin->size();

    outTree->Fill();
  }

  // Write the new tree to a new file

  outfile->cd();
  outTree->Write();
  // Close the files
  outfile->Close();
  file->Close();
}