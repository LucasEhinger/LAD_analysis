#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TTree.h>
#include <TText.h>
#include <iostream>
#include <string>
#include <vector>

void edep_plots() {
  double energies[] = {300, 400, 500, 600, 700};

  std::vector<double> avg_edep;
  std::vector<double> avg_time;
  std::string fileloc = "/volatile/hallc/c-lad/ehingerl/G4_LAD/carlos_proton/ana/";
  std::vector<std::string> filenames;
  for (const auto &energy : energies) {
    filenames.push_back(fileloc + "ScanLAD_proton_" + std::to_string(static_cast<int>(energy)) +
                        "MeV_10k_20240205_ana.root");
  }
  for (const auto &filename : filenames) {
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
    tree->SetBranchAddress("vEDep", &edep);
    tree->SetBranchAddress("vTbar", &time);

    double sum_edep   = 0;
    double sum_time   = 0;
    int nhits         = 0;
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
      tree->GetEntry(i);
      for (size_t j = 0; j < edep->size(); ++j) {
        sum_edep += edep->at(j);
        sum_time += time->at(j);
        nhits++;
      }
    }

    avg_edep.push_back(sum_edep / nhits);
    avg_time.push_back(sum_time / nhits);
    file.Close();
  }

  // TCanvas *c1 = new TCanvas("c1", "Average Energy Deposition and Time", 800, 600);
  // TGraph *graph_edep = new TGraph(sizeof(energies) / sizeof(energies[0]), energies, avg_edep.data());
  // graph_edep->SetTitle("Average Energy Deposition vs. Energy;Energy (MeV);Average Energy Deposition (MeV)");
  // graph_edep->SetMarkerStyle(20);
  // graph_edep->SetMarkerColor(kRed);
  // graph_edep->Draw("APL");

  // TGraph *graph_time = new TGraph(sizeof(energies) / sizeof(energies[0]), energies, avg_time.data());
  // graph_time->SetTitle("Average Time vs. Energy;Energy (MeV);Average Time (ns)");
  // graph_time->SetMarkerStyle(21);
  // graph_time->SetMarkerColor(kBlue);
  // graph_time->Draw("PL SAME");

  // c1->BuildLegend();
  // c1->SaveAs("../hists/avg_edep_time_vs_energy.pdf");
  // delete c1;
  TCanvas *c2                = new TCanvas("c2", "Energy Deposition vs. Time of Flight", 800, 600);
  TGraph *graph_edep_vs_time = new TGraph(avg_time.size());
  for (size_t i = 0; i < avg_time.size(); ++i) {
    graph_edep_vs_time->SetPoint(i, avg_time[i], avg_edep[i] / 5);
  }
  graph_edep_vs_time->SetTitle("Energy Deposition vs. Time of Flight;Time of Flight (ns);Energy Deposition / 5 cm");
  graph_edep_vs_time->SetMarkerStyle(22);
  graph_edep_vs_time->SetMarkerColor(kRed);
  graph_edep_vs_time->SetLineColor(kRed);
  graph_edep_vs_time->SetLineStyle(2); // Dashed line
  graph_edep_vs_time->Draw("AP");

  for (size_t i = 0; i < avg_time.size(); ++i) {
    TText *label = new TText(avg_time[i] + 1, avg_edep[i] / 5, std::to_string(static_cast<int>(energies[i])).c_str());
    label->SetTextSize(0.02);
    label->SetTextAlign(12); // Align text to the left
    label->Draw();
  }

  c2->SaveAs("../hists/edep_vs_time_of_flight.pdf");
  delete c2;
}