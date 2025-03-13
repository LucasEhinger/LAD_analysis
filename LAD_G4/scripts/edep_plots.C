#include <TFile.h>
#include <TH2D.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

const int t_bins = 100;
const int e_bins = 100;
const double t_min = 0;
const double t_max = 20;
const double e_min = 0;
const double e_max = 30;


void edep_plots() {
  double momenta[]    = {300, 400, 500, 600, 700};
  std::string fileloc = "/volatile/hallc/c-lad/ehingerl/G4_LAD/carlos_proton/ana/";
  std::vector<std::string> filenames;
  for (const auto &energy : momenta) {
    filenames.push_back(fileloc + "ScanLAD_proton_" + std::to_string(static_cast<int>(energy)) +
                        "MeV_10k_20240205_ana.root");
  }

  int colors[]    = {kRed, kBlue, kGreen, kMagenta, kCyan};

  double min_time = std::numeric_limits<double>::max();
  double max_time = std::numeric_limits<double>::lowest();
  double min_edep = std::numeric_limits<double>::max();
  double max_edep = std::numeric_limits<double>::lowest();

  std::vector<TH2D *> histograms;
  //Add histograms here. Adding inside loop caused unknown seg fault.
  for (size_t idx = 0; idx < filenames.size(); ++idx) {
    TH2D *hist_edep_vs_time =
      new TH2D(Form("hist_%d", static_cast<int>(momenta[idx])),
           "Energy Deposition vs. Time of Flight;Time of Flight (ns);Energy Deposition / 5 cm", t_bins, t_min, t_max, e_bins, e_min, e_max);
    histograms.push_back(hist_edep_vs_time);
  }


  for (size_t idx = 0; idx < filenames.size(); ++idx) {
    const auto &filename = filenames[idx];
    TFile file(filename.c_str(), "READ");
    if (file.IsZombie()) {
      std::cerr << "Error: Could not open file " << filename << std::endl;
      continue;
    }

    TTree *tree = nullptr;
    file.GetObject("T", tree);
    if (!tree) {
      std::cerr << "Error: Could not find tree in file " << filename << std::endl;
      file.Close();
      continue;
    }

    std::vector<double> *edep    = nullptr;
    std::vector<double> *time    = nullptr;
    std::vector<double> *vPaddle = nullptr;
    std::vector<double> *vXbar   = nullptr;
    std::vector<double> *vYbar   = nullptr;
    std::vector<double> *vZbar   = nullptr;
    tree->SetBranchAddress("vEDep", &edep);
    tree->SetBranchAddress("vTbar", &time);
    tree->SetBranchAddress("vPaddle", &vPaddle);
    tree->SetBranchAddress("vXbar", &vXbar);
    tree->SetBranchAddress("vYbar", &vYbar);
    tree->SetBranchAddress("vZbar", &vZbar);

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
      tree->GetEntry(i);
      std::map<int, double> total_edep;
      std::map<int, double> event_time;
      std::map<int, double> event_dist;
      for (size_t j = 0; j < edep->size(); ++j) {
        int paddle = vPaddle->at(j);
        if (total_edep.find(paddle) != total_edep.end()) {
          total_edep[paddle] += edep->at(j);
        } else {
          total_edep[paddle] = edep->at(j);
          event_time[paddle] = time->at(j);
          event_dist[paddle] = sqrt(pow(vXbar->at(j), 2) + pow(vYbar->at(j), 2) + pow(vZbar->at(j), 2)) * pow(10, -3);
        }
      }

      for (const auto &paddle : total_edep) {
        double e_pcm  = paddle.second / 5;
        double tof_pm = event_time[paddle.first] / event_dist[paddle.first];
        histograms[idx]->Fill(tof_pm, e_pcm);
        if (tof_pm < min_time)
          min_time = tof_pm;
        if (tof_pm > max_time)
          max_time = tof_pm;
        if (e_pcm < min_edep)
          min_edep = e_pcm;
        if (e_pcm > max_edep)
          max_edep = e_pcm;
      }
    }
    file.Close();
  }

  // Create a new ROOT file to save the histograms
  TFile *outputFile = new TFile("../hists/edep_vs_time_of_flight.root", "RECREATE");

  // Create a histogram to store the sum of all histograms
  TH2D *hist_sum = new TH2D("hist_sum", "Sum of Energy Deposition vs. Time of Flight;Time of Flight (ns);Energy Deposition / 5 cm", t_bins, t_min, t_max, e_bins, e_min, e_max);

  for (size_t idx = 0; idx < histograms.size(); ++idx) {
    histograms[idx]->Write();
    cout << "Address of histogram " << idx << ": " << histograms[idx] << endl;
    hist_sum->Add(histograms[idx]);
  }
  // Create a copy of hist_sum with bins less than min_bin set to zero
  TH2D *hist_sum_min_bin = (TH2D *)hist_sum->Clone("hist_sum_min_bin");
  int min_bin = 50; // Example value for min_bin

  for (int x = 1; x <= hist_sum_min_bin->GetNbinsX(); ++x) {
    for (int y = 1; y <= hist_sum_min_bin->GetNbinsY(); ++y) {
      if (hist_sum_min_bin->GetBinContent(x, y) < min_bin) {
        hist_sum_min_bin->SetBinContent(x, y, 0);
      }
    }
  }

  hist_sum_min_bin->Write();
  hist_sum->Write();
  outputFile->Close();

  // Clean up
  for (auto hist : histograms) {
    delete hist;
  }
  delete outputFile;
}
