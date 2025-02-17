#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TTree.h>
#include <TText.h>
#include <TLegend.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

const int paddle_num = 10000;

void edep_plots() {
  double energies[] = {300, 400, 500, 600, 700};
  std::string fileloc = "/volatile/hallc/c-lad/ehingerl/G4_LAD/carlos_proton/ana/";
  std::vector<std::string> filenames;
  for (const auto &energy : energies) {
    filenames.push_back(fileloc + "ScanLAD_proton_" + std::to_string(static_cast<int>(energy)) +
                        "MeV_10k_20240205_ana.root");
  }

  TCanvas *c2 = new TCanvas("c2", "Energy Deposition vs. Time of Flight", 800, 600);
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  int colors[] = {kRed, kBlue, kGreen, kMagenta, kCyan};

  double min_time = std::numeric_limits<double>::max();
  double max_time = std::numeric_limits<double>::lowest();
  double min_edep = std::numeric_limits<double>::max();
  double max_edep = std::numeric_limits<double>::lowest();

  std::vector<TGraph*> graphs;

  for (size_t idx = 0; idx < filenames.size(); ++idx) {
    const auto &filename = filenames[idx];
    TFile file(filename.c_str(), "READ");
    if (file.IsZombie()) {
      std::cerr << "Error: Could not open file " << filename << std::endl;
      continue;
    }

    TTree *tree;
    file.GetObject("T", tree); // Replace "tree_name" with the actual name of the tree
    if (!tree) {
      std::cerr << "Error: Could not find tree in file " << filename << std::endl;
      file.Close();
      continue;
    }

    std::vector<double> *edep = nullptr;
    std::vector<double> *time = nullptr;
    std::vector<double> *vPaddle = nullptr;
    tree->SetBranchAddress("vEDep", &edep);
    tree->SetBranchAddress("vTbar", &time);
    tree->SetBranchAddress("vPaddle", &vPaddle);

    TGraph *graph_edep_vs_time = new TGraph();
    Long64_t nentries = tree->GetEntries();
    int point_idx = 0;
    for (Long64_t i = 0; i < nentries; ++i) {
      tree->GetEntry(i);
      double total_edep = 0;
      double event_time = 0;
      bool valid_event = false;
      for (size_t j = 0; j < edep->size(); ++j) {
        if (vPaddle->at(j) == paddle_num) {
          total_edep += edep->at(j);
          event_time = time->at(j);
          valid_event = true;
        }
      }
      if (valid_event) {
        double e = total_edep / 5;
        graph_edep_vs_time->SetPoint(point_idx++, event_time, e);
        if (event_time < min_time) min_time = event_time;
        if (event_time > max_time) max_time = event_time;
        if (e < min_edep) min_edep = e;
        if (e > max_edep) max_edep = e;
      }
    }

    graph_edep_vs_time->SetMarkerStyle(20 + idx);
    graph_edep_vs_time->SetMarkerColor(colors[idx]);
    graph_edep_vs_time->SetLineColor(colors[idx]);
    graph_edep_vs_time->SetTitle("Energy Deposition vs. Time of Flight;Time of Flight (ns);Energy Deposition / 5 cm");
    graphs.push_back(graph_edep_vs_time);
    legend->AddEntry(graph_edep_vs_time, (std::to_string(static_cast<int>(energies[idx])) + " MeV").c_str(), "p");

    file.Close();
  }

  for (size_t idx = 0; idx < graphs.size(); ++idx) {
    if (idx == 0) {
      graphs[idx]->Draw("AP");
      graphs[idx]->GetXaxis()->SetLimits(min_time, max_time);
      graphs[idx]->GetYaxis()->SetRangeUser(min_edep, max_edep);
    } else {
      graphs[idx]->Draw("P SAME");
    }
  }

  legend->Draw();
  c2->SaveAs("../hists/edep_vs_time_of_flight.pdf");
  delete c2;
}
