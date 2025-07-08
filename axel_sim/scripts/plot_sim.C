#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <vector>

const int N_PLANES                  = 5;
const std::string PLANE_NAMES[N_PLANES] = {"000", "001", "100", "101", "200"}; // Names of the planes
const int N_MOMENTA                 = 9;
const double MOMENTA[N_MOMENTA + 1] = {0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}; // Example momenta values
const int N_PADDLES                 = 11;
const double FILE_CHARGE            = 60 * 60 * 0.3 / 1000; // Charge in mC
const double CHARGE_NORM            = 10;                   // Target charge normalization in mC
const char *data_filename = "/home/ehingerl/hallc/analysis/lad_rates/files/lad_edep_plots_LD2_setting3_P.root";

void plot_sim() {
  gROOT->SetBatch(kTRUE);
  const char *filename = "../source_files/lad_sim_shms_13.5.root";
  TFile *file          = TFile::Open(filename, "READ");
  TTree *T             = (TTree *)file->Get("T");

  // Example: Assume tree has a branch "value" (double) and indices plane, mom, paddle (int)
  double weight, theta_p_recon, phi_p_recon, mom_p_recon;
  int lad_plane, lad_bar;
  T->SetBranchAddress("weight", &weight);
  T->SetBranchAddress("theta_p_recon", &theta_p_recon);
  T->SetBranchAddress("phi_p_recon", &phi_p_recon);
  T->SetBranchAddress("mom_p_recon", &mom_p_recon);
  T->SetBranchAddress("lad_plane", &lad_plane);
  T->SetBranchAddress("lad_bar", &lad_bar);

  // Create histograms
  TH1D *hists[N_PLANES][N_MOMENTA];
  for (int i = 0; i < N_PLANES; ++i) {
    for (int j = 0; j < N_MOMENTA; ++j) {
      std::string hname = Form("h_p%d_%.1f_%.1f", i, MOMENTA[j], MOMENTA[j + 1]);
      hists[i][j]       = new TH1D(hname.c_str(), hname.c_str(), N_PADDLES, -0.5, N_PADDLES - 0.5);
      hists[i][j]->GetXaxis()->SetTitle("Bar number");
      hists[i][j]->GetYaxis()->SetTitle("Counts");
      std::string title = Form("Plane %d, %.0f - %.0f MeV/c", i, MOMENTA[j] * 1000, MOMENTA[j + 1] * 1000);
      hists[i][j]->SetTitle(title.c_str());
    }
  }

  Long64_t nentries = T->GetEntries();
  for (Long64_t entry = 0; entry < nentries; ++entry) {
    T->GetEntry(entry);

    int plane  = lad_plane;
    int paddle = N_PADDLES - lad_bar - 1;
    // Find momentum bin
    int mom_bin = -1;
    for (int j = 0; j < N_MOMENTA; ++j) {
      if (mom_p_recon >= MOMENTA[j] && mom_p_recon < MOMENTA[j + 1]) {
        mom_bin = j;
        break;
      }
    }
    if (mom_bin == -1)
      continue;

    int true_plane1 = 4 - plane * 2;
    int true_plane2 = 4 - plane * 2 + 1;

    hists[true_plane1][mom_bin]->Fill(paddle, weight);
    if (true_plane2 < N_PLANES) {
      hists[true_plane2][mom_bin]->Fill(paddle, weight);
    }
  }
  // Create summed histograms over all momenta for each plane
  TH1D *hists_sum[N_PLANES];
  TH1D *hists_stoped[N_PLANES];
  TH1D *hists_punch_through[N_PLANES];
  for (int i = 0; i < N_PLANES; ++i) {
    std::string hname_sum = Form("h_p%d_sum", i);
    hists_sum[i]          = (TH1D *)hists[i][0]->Clone(hname_sum.c_str());
    hists_sum[i]->SetTitle(hname_sum.c_str());
    hists_sum[i]->Reset();
    std::string hname_stoped = Form("h_p%d_stoped", i);
    hists_stoped[i]          = (TH1D *)hists[i][0]->Clone(hname_stoped.c_str());
    hists_stoped[i]->SetTitle(hname_stoped.c_str());
    hists_stoped[i]->Reset();
    std::string hname_punch_through = Form("h_p%d_punch_through", i);
    hists_punch_through[i]          = (TH1D *)hists[i][0]->Clone(hname_punch_through.c_str());
    hists_punch_through[i]->SetTitle(hname_punch_through.c_str());
    hists_punch_through[i]->Reset();
    for (int j = 0; j < N_MOMENTA; ++j) {
      // Scale all histograms by the charge normalization
      hists[i][j]->Scale(CHARGE_NORM / FILE_CHARGE);
      hists_sum[i]->Add(hists[i][j]);
      if (MOMENTA[j] < 0.4 && MOMENTA[j] >= 0.3) {
        // For momenta below 0.4 GeV/c, consider it as stopped
        hists_stoped[i]->Add(hists[i][j]);
      } else if (MOMENTA[j] >= 0.4) {
        // For momenta between 0.5 and 1.0 GeV/c, consider it as punch-through
        hists_punch_through[i]->Add(hists[i][j]);
      }
    }
  }

  TFile *outfile = new TFile("../files/histograms.root", "RECREATE");
  for (int i = 0; i < N_PLANES; ++i) {
    std::string folder_name = Form("sim_only/plane_%d", i);
    outfile->mkdir(folder_name.c_str());
    outfile->cd(folder_name.c_str());
    for (int j = 0; j < N_MOMENTA; ++j)
      hists[i][j]->Write();
    hists_sum[i]->Write();
    hists_stoped[i]->Write();
    hists_punch_through[i]->Write();
    outfile->cd(); // Go back to the root directory for the next plane
  }

  // Grab plots from data, and make comparison plots
  TFile *datafile = TFile::Open(data_filename, "READ");
  TCanvas *c      = (TCanvas *)datafile->Get("KIN/LAD_RATES/c_PROTON_RATES_all_planes_rates");
  if (!c) {
    std::cerr << "Could not find canvas KIN/LAD_RATES/c_PROTON_RATES_all_planes_rates in data file." << std::endl;
  } else {
    std::vector<TGraph *> data_graphs;
    for (int i = 0; i < N_PLANES; ++i) {
      // TObject *prim = c->GetListOfPrimitives()->At(i);
      // if (prim)
      //   std::cout << "Primitive " << i << " is of type: " << prim->ClassName() << std::endl;
      TPad *pad = dynamic_cast<TPad *>(c->GetListOfPrimitives()->At(i));
      if (pad) {
        // Assume the first TGraph in the pad is the data graph
        TGraph *g = nullptr;
        TIter next(pad->GetListOfPrimitives());
        TObject *obj = nullptr;
        while ((obj = next())) {
          g = dynamic_cast<TGraph *>(obj);
          if (g)
            break;
        }
        if (g) {
          data_graphs.push_back(g);
        } else {
          std::cerr << "Could not find TGraph in pad for plane " << i << std::endl;
          data_graphs.push_back(nullptr);
        }
      } else {
        std::cerr << "Could not find TPad for plane " << i << std::endl;
        data_graphs.push_back(nullptr);
      }
      // Create a canvas for comparison
      std::string comp_folder = "comp";
      // Check if the folder already exists
      if (!outfile->GetDirectory(comp_folder.c_str())) {
        outfile->mkdir(comp_folder.c_str());
      }
      outfile->cd(comp_folder.c_str());

      TCanvas *comp_c = new TCanvas(Form("c_comp_plane_%d", i), Form("Comparison Plane %d", i), 800, 600);

      // Determine y-axis range to fit both simulation and data
      double sim_max  = hists_punch_through[i]->GetMaximum();
      double data_max = 0;
      if (data_graphs[i]) {
        double *y = data_graphs[i]->GetY();
        int n     = data_graphs[i]->GetN();
        for (int k = 0; k < n; ++k) {
          if (y[k] > data_max)
            data_max = y[k];
        }
      }
      double y_max = std::max(sim_max, data_max) * 1.2; // Add 20% headroom
      hists_punch_through[i]->SetMaximum(y_max);
      hists_punch_through[i]->SetMinimum(0);

      // Draw simulation histogram (summed over momenta)
      hists_punch_through[i]->SetLineColor(kRed);
      hists_punch_through[i]->SetMarkerColor(kRed);
      hists_punch_through[i]->SetMarkerStyle(20);
      hists_punch_through[i]->SetStats(0);
      hists_punch_through[i]->Draw("E");

      // Draw data graph
      if (data_graphs[i]) {
        data_graphs[i]->SetLineColor(kBlue + 2);
        data_graphs[i]->SetMarkerColor(kBlue + 2);
        data_graphs[i]->SetMarkerStyle(21);
        data_graphs[i]->Draw("P SAME");
      }

      auto legend = new TLegend(0.7, 0.75, 0.9, 0.9);
      legend->AddEntry(hists_punch_through[i], "Sim", "lep");
      if (data_graphs[i])
        legend->AddEntry(data_graphs[i], "Data", "lep");
      legend->Draw();
      comp_c->Write();
      outfile->cd();
    }
  }

  outfile->Close();

  file->Close();
}