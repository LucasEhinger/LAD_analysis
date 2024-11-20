#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TChain.h>
#include <ROOT/RDataFrame.hxx>

using namespace std;
using namespace ROOT;

double get_Momentum(int n_hits, RVec<Int_t> pdgid, RVec<Int_t> detid, RVec<Int_t> trackid, RVec<double> p,
                    double p_init, int gem_ID) {
  for (int i = 0; i < n_hits; i++) {
    if (detid[i] == gem_ID && pdgid[i] == 2212 && trackid[i] == 1) {
      return p[i];
    }
  }
  return 0;
}

double get_Edep(int n_hits, RVec<Int_t> pdgid, RVec<Int_t> detid, RVec<Int_t> trackid, RVec<double> edep, int gem_ID) {
  for (int i = 0; i < n_hits; i++) {
    if (detid[i] == gem_ID && pdgid[i] == 2212 && trackid[i] == 1) {
      return edep[i];
    }
  }
  return 0;
}

double get_delta_p(int n_hits, RVec<Int_t> pdgid, RVec<Int_t> detid, RVec<Int_t> trackid, RVec<double> p, double p_init,
                   int gem_ID) {
  for (int i = 0; i < n_hits; i++) {
    if (detid[i] == gem_ID && pdgid[i] == 2212 && trackid[i] == 1) {
      return p_init - p[i];
    }
  }
  return 0;
}

double get_delta_E(int n_hits, RVec<Int_t> pdgid, RVec<Int_t> detid, RVec<Int_t> trackid, RVec<double> p, double p_init,
                   int gem_ID) {
  for (int i = 0; i < n_hits; i++) {
    if (detid[i] == gem_ID && pdgid[i] == 2212 && trackid[i] == 1) {
      return sqrt(pow(0.938, 2) + pow(p_init, 2)) - sqrt(pow(0.938, 2) + pow(p[i], 2));
    }
  }
  return 0;
}

double get_Beta(int n_hits, RVec<Int_t> pdgid, RVec<Int_t> detid, RVec<Int_t> trackid, RVec<double> p,
                RVec<double> x_vec, RVec<double> y_vec, RVec<double> z_vec, RVec<double> t_vec, int gem_ID) {
  double x, y, z, t;
  for (int i = 0; i < n_hits; i++) {
    if (detid[i] == gem_ID && pdgid[i] == 2212 && trackid[i] == 1) {
      x = x_vec[i];
      y = y_vec[i];
      z = z_vec[i];
      t = t_vec[i];
      break;
    }
  }
  double r    = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  double c    = 3 * pow(10, 8);
  double beta = r / (c * t) * pow(10, 9);
  return beta;
}

double get_Beta_2det(int n_hits, RVec<Int_t> pdgid, RVec<Int_t> detid, RVec<Int_t> trackid, RVec<double> p,
                     RVec<double> x_vec, RVec<double> y_vec, RVec<double> z_vec, RVec<double> t_vec, int det_ID1,
                     int det_ID2) {
  double x, y, z, t;
  for (int i = 0; i < n_hits; i++) {
    if (detid[i] == det_ID1 && pdgid[i] == 2212 && trackid[i] == 1) {
      x = x_vec[i];
      y = y_vec[i];
      z = z_vec[i];
      t = t_vec[i];
      break;
    }
  }

  for (int i = 0; i < n_hits; i++) {
    if (detid[i] == det_ID2 && pdgid[i] == 2212 && trackid[i] == 1) {
      x -= x_vec[i];
      y -= y_vec[i];
      z -= z_vec[i];
      t -= t_vec[i];
      break;
    }
  }

  double r    = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  double c    = 3 * pow(10, 8);
  double beta = r / (c * t) * pow(10, 9);
  return abs(beta);
}

double get_tof_Momentum(double beta) {
  double m     = 0.938;
  double gamma = 1 / sqrt(1 - pow(beta, 2));
  double p     = m * gamma * beta;
  return p;
}

void proton_momentum(string inFileName, string outHistName, string treeName) {

  int n_bins                    = 1000 + 1;
  int n_events                  = pow(10, 6);
  string norm_factor_GEM        = to_string(((double)n_bins) / n_events);
  string norm_factor_GEM_Theta1 = to_string(((double)n_bins) / n_events * 7 / 2);
  string norm_factor_GEM_Theta2 = to_string(((double)n_bins) / n_events * 7 / 3);
  TFile *histFile               = new TFile(outHistName.c_str(), "RECREATE");
  histFile->cd();

  TChain chain(treeName.c_str());
  chain.Add(inFileName.c_str());

  RDataFrame rdf_raw(chain);
  auto rdf_def = rdf_raw.Define("weight", norm_factor_GEM)
                     .Define("weight_Theta1", norm_factor_GEM_Theta1)
                     .Define("weight_Theta2", norm_factor_GEM_Theta2)
                     .Define("weight_binned_momentum", "weight*5")
                     .Define("GEM1_p", "get_Momentum(hit.n,hit.pid,hit.det,hit.trid,hit.p,ev.p,101)")
                     .Define("GEM2_p", "get_Momentum(hit.n,hit.pid,hit.det,hit.trid,hit.p,ev.p,102)")
                     .Define("LAD_p", "get_Momentum(hit.n,hit.pid,hit.det,hit.trid,hit.p,ev.p,501)")
                     .Define("GEM1_edep", "get_Edep(hit.n,hit.pid,hit.det,hit.trid,hit.edep,101)")
                     .Define("GEM2_edep", "get_Edep(hit.n,hit.pid,hit.det,hit.trid,hit.edep,102)")
                     .Define("theta", "ev.theta*180/3.14159")
                     .Define("delta_p_GEM1", "get_delta_p(hit.n,hit.pid,hit.det,hit.trid,hit.p,ev.p,101)")
                     .Define("delta_p_GEM2", "get_delta_p(hit.n,hit.pid,hit.det,hit.trid,hit.p,ev.p,102)")
                     .Define("delta_p_LAD", "get_delta_p(hit.n,hit.pid,hit.det,hit.trid,hit.p,ev.p,501)")
                     .Define("delta_E_GEM1", "get_delta_E(hit.n,hit.pid,hit.det,hit.trid,hit.p,ev.p,101)")
                     .Define("delta_E_GEM2", "get_delta_E(hit.n,hit.pid,hit.det,hit.trid,hit.p,ev.p,102)")
                     .Define("delta_E_LAD", "get_delta_E(hit.n,hit.pid,hit.det,hit.trid,hit.p,ev.p,501)");

#pragma region detector_momentum_plots
  TDirectory *Detector_momentum_plots = histFile->mkdir("Detector momentum_plots");
  Detector_momentum_plots->cd();

  auto h_GEM1_p =
      *rdf_def.Filter("GEM1_p > 0")
           .Histo1D({"GEM1 Momentum", "GEM1 Momentum; p [GeV/c]; Survival Fraction", n_bins, 0, 1}, "GEM1_p", "weight");
  h_GEM1_p.Write();

  auto h_GEM2_p =
      *rdf_def.Filter("GEM2_p > 0")
           .Histo1D({"GEM2 Momentum", "GEM2 Momentum; p [GeV/c]; Survival Fraction", n_bins, 0, 1}, "GEM2_p", "weight");
  h_GEM2_p.Write();

  auto h_LAD_p =
      *rdf_def.Filter("LAD_p > 0")
           .Histo1D({"LAD Momentum", "LAD Momentum; p [GeV/c]; Survival Fraction", n_bins, 0, 1}, "LAD_p", "weight");
  h_LAD_p.Write();

  auto h_GEM1_p_theta_95_110 =
      *rdf_def.Filter("GEM1_p > 0 && theta > 95 && theta < 110")
           .Histo1D({"GEM1 Momentum #Theta #in [95, 110]",
                     "GEM1 Momentum #Theta #in [95, 110]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "GEM1_p", "weight_Theta1");
  h_GEM1_p_theta_95_110.Write();

  auto h_GEM1_p_theta_110_130 =
      *rdf_def.Filter("GEM1_p > 0 && theta > 110 && theta < 130")
           .Histo1D({"GEM1 Momentum #Theta #in [110, 130]",
                     "GEM1 Momentum #Theta #in [110, 130]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "GEM1_p", "weight_Theta1");
  h_GEM1_p_theta_110_130.Write();

  auto h_GEM1_p_theta_130_160 =
      *rdf_def.Filter("GEM1_p > 0 && theta > 130 && theta < 160")
           .Histo1D({"GEM1 Momentum #Theta #in [130, 160]",
                     "GEM1 Momentum #Theta #in [130, 160]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "GEM1_p", "weight_Theta2");
  h_GEM1_p_theta_130_160.Write();

  auto h_GEM2_p_theta_95_110 =
      *rdf_def.Filter("GEM2_p > 0 && theta > 95 && theta < 110")
           .Histo1D({"GEM2 Momentum #Theta #in [95, 110]",
                     "GEM2 Momentum #Theta #in [95, 110]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "GEM2_p", "weight_Theta1");
  h_GEM2_p_theta_95_110.Write();

  auto h_GEM2_p_theta_110_130 =
      *rdf_def.Filter("GEM2_p > 0 && theta > 110 && theta < 130")
           .Histo1D({"GEM2 Momentum #Theta #in [110, 130]",
                     "GEM2 Momentum #Theta #in [110, 130]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "GEM2_p", "weight_Theta1");
  h_GEM2_p_theta_110_130.Write();

  auto h_GEM2_p_theta_130_160 =
      *rdf_def.Filter("GEM2_p > 0 && theta > 130 && theta < 160")
           .Histo1D({"GEM2 Momentum #Theta #in [130, 160]",
                     "GEM2 Momentum #Theta #in [130, 160]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "GEM2_p", "weight_Theta2");
  h_GEM2_p_theta_130_160.Write();

  auto h_LAD_p_theta_95_110 =
      *rdf_def.Filter("LAD_p > 0 && theta > 95 && theta < 110")
           .Histo1D({"LAD Momentum #Theta #in [95, 110]",
                     "LAD Momentum #Theta #in [95, 110]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "LAD_p", "weight_Theta1");
  h_LAD_p_theta_95_110.Write();

  auto h_LAD_p_theta_110_130 =
      *rdf_def.Filter("LAD_p > 0 && theta > 110 && theta < 130")
           .Histo1D({"LAD Momentum #Theta #in [110, 130]",
                     "LAD Momentum #Theta #in [110, 130]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "LAD_p", "weight_Theta1");
  h_LAD_p_theta_110_130.Write();

  auto h_LAD_p_theta_130_160 =
      *rdf_def.Filter("LAD_p > 0 && theta > 130 && theta < 160")
           .Histo1D({"LAD Momentum #Theta #in [130, 160]",
                     "LAD Momentum #Theta #in [130, 160]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "LAD_p", "weight_Theta2");
  h_LAD_p_theta_130_160.Write();
#pragma endregion detector_momentum_plots

#pragma region Original_momentum_plots
  TDirectory *Original_momentum_plots = histFile->mkdir("Original_momentum_plots");
  Original_momentum_plots->cd();

  auto h_GEM1_p_original =
      *rdf_def.Filter("GEM1_p > 0")
           .Histo1D({"GEM1 Original Momentum", "GEM1 Original Momentum; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "ev.p", "weight");
  h_GEM1_p_original.Write();

  auto h_GEM2_p_original =
      *rdf_def.Filter("GEM2_p > 0")
           .Histo1D({"GEM2 Original Momentum", "GEM2 Original Momentum; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "ev.p", "weight");
  h_GEM2_p_original.Write();

  auto h_LAD_p_original =
      *rdf_def.Filter("LAD_p > 0")
           .Histo1D({"LAD Original Momentum", "LAD Original Momentum; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "ev.p", "weight");
  h_LAD_p_original.Write();

  auto h_GEM1_p_original_theta_95_110 =
      *rdf_def.Filter("GEM1_p > 0 && theta > 95 && theta < 110")
           .Histo1D({"GEM1 Original Momentum #Theta #in [95, 110]",
                     "GEM1 Original Momentum #Theta #in [95, 110]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "ev.p", "weight_Theta1");
  h_GEM1_p_original_theta_95_110.Write();

  auto h_GEM1_p_original_theta_110_130 =
      *rdf_def.Filter("GEM1_p > 0 && theta > 110 && theta < 130")
           .Histo1D({"GEM1 Original Momentum #Theta #in [110, 130]",
                     "GEM1 Original Momentum #Theta #in [110, 130]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "ev.p", "weight_Theta1");
  h_GEM1_p_original_theta_110_130.Write();

  auto h_GEM1_p_original_theta_130_160 =
      *rdf_def.Filter("GEM1_p > 0 && theta > 130 && theta < 160")
           .Histo1D({"GEM1 Original Momentum #Theta #in [130, 160]",
                     "GEM1 Original Momentum #Theta #in [130, 160]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "ev.p", "weight_Theta2");
  h_GEM1_p_original_theta_130_160.Write();

  auto h_GEM2_p_original_theta_95_110 =
      *rdf_def.Filter("GEM2_p > 0 && theta > 95 && theta < 110")
           .Histo1D({"GEM2 Original Momentum #Theta #in [95, 110]",
                     "GEM2 Original Momentum #Theta #in [95, 110]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "ev.p", "weight_Theta1");
  h_GEM2_p_original_theta_95_110.Write();

  auto h_GEM2_p_original_theta_110_130 =
      *rdf_def.Filter("GEM2_p > 0 && theta > 110 && theta < 130")
           .Histo1D({"GEM2 Original Momentum #Theta #in [110, 130]",
                     "GEM2 Original Momentum #Theta #in [110, 130]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "ev.p", "weight_Theta1");
  h_GEM2_p_original_theta_110_130.Write();

  auto h_GEM2_p_original_theta_130_160 =
      *rdf_def.Filter("GEM2_p > 0 && theta > 130 && theta < 160")
           .Histo1D({"GEM2 Original Momentum #Theta #in [130, 160]",
                     "GEM2 Original Momentum #Theta #in [130, 160]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "ev.p", "weight_Theta2");
  h_GEM2_p_original_theta_130_160.Write();

  auto h_LAD_p_original_theta_95_110 =
      *rdf_def.Filter("LAD_p > 0 && theta > 95 && theta < 110")
           .Histo1D({"LAD Original Momentum #Theta #in [95, 110]",
                     "LAD Original Momentum #Theta #in [95, 110]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "ev.p", "weight_Theta1");
  h_LAD_p_original_theta_95_110.Write();

  auto h_LAD_p_original_theta_110_130 =
      *rdf_def.Filter("LAD_p > 0 && theta > 110 && theta < 130")
           .Histo1D({"LAD Original Momentum #Theta #in [110, 130]",
                     "LAD Original Momentum #Theta #in [110, 130]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "ev.p", "weight_Theta1");
  h_LAD_p_original_theta_110_130.Write();

  auto h_LAD_p_original_theta_130_160 =
      *rdf_def.Filter("LAD_p > 0 && theta > 130 && theta < 160")
           .Histo1D({"LAD Original Momentum #Theta #in [130, 160]",
                     "LAD Original Momentum #Theta #in [130, 160]; p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                    "ev.p", "weight_Theta2");
  h_LAD_p_original_theta_130_160.Write();

#pragma endregion Original_momentum_plots

#pragma region Theta_plots
  TDirectory *Theta_plots = histFile->mkdir("Theta_plots");
  Theta_plots->cd();

  auto h_GEM1_theta =
      *rdf_def.Filter("GEM1_p > 0")
           .Histo1D({"GEM1 Theta", "GEM1 Theta; theta [Deg]; Survival Fraction", n_bins, 95, 160}, "theta", "weight");
  h_GEM1_theta.Write();

  auto h_GEM2_theta =
      *rdf_def.Filter("GEM2_p > 0")
           .Histo1D({"GEM2 Theta", "GEM2 Theta; theta [Deg]; Survival Fraction", n_bins, 95, 160}, "theta", "weight");
  h_GEM2_theta.Write();

  auto h_LAD_theta =
      *rdf_def.Filter("(LAD_p > 0)")
           .Histo1D({"LAD Theta", "LAD Theta; theta [Deg]; Survival Fraction", n_bins, 95, 160}, "theta", "weight");
  h_LAD_theta.Write();

  for (double p = 0; p <= 1; p += 0.05) {
    string p_str     = to_string(p);
    p_str            = p_str.substr(0, p_str.find('.') + 3); // Limit to 2 decimal places
    string p_str_end = to_string(p + 0.05);
    p_str_end        = p_str_end.substr(0, p_str_end.find('.') + 3);

    string name         = "GEM1 Theta p #in [" + p_str + ", " + p_str_end + "]";
    string title        = "GEM1 Theta p #in [" + p_str + ", " + p_str_end + "]; theta [Deg]; Survival Fraction";
    auto h_GEM1_theta_p = *rdf_def.Filter("GEM1_p > 0 && ev.p > " + p_str + " && ev.p < " + p_str_end)
                               .Define("weight_theta", to_string(((double)n_bins) / n_events / 0.05))
                               .Histo1D({name.c_str(), title.c_str(), n_bins, 95, 160}, "theta", "weight_theta");
    h_GEM1_theta_p.Write();

    name                = "GEM2 Theta p #in [" + p_str + ", " + p_str_end + "]";
    title               = "GEM2 Theta p #in [" + p_str + ", " + p_str_end + "]; theta [Deg]; Survival Fraction";
    auto h_GEM2_theta_p = *rdf_def.Filter("GEM2_p > 0 && ev.p > " + p_str + " && ev.p < " + p_str_end)
                               .Define("weight_theta", to_string(((double)n_bins) / n_events / 0.05))
                               .Histo1D({name.c_str(), title.c_str(), n_bins, 95, 160}, "theta", "weight_theta");
    h_GEM2_theta_p.Write();

    name               = "LAD Theta p #in [" + p_str + ", " + p_str_end + "]";
    title              = "LAD Theta p #in [" + p_str + ", " + p_str_end + "]; theta [Deg]; Survival Fraction";
    auto h_LAD_theta_p = *rdf_def.Filter("LAD_p > 0 && ev.p > " + p_str + " && ev.p < " + p_str_end)
                              .Define("weight_theta", to_string(((double)n_bins) / n_events / 0.05))
                              .Histo1D({name.c_str(), title.c_str(), n_bins, 95, 160}, "theta", "weight_theta");
    h_LAD_theta_p.Write();
  }

#pragma endregion Theta_plots

#pragma region EP_Loss_plots

  TDirectory *EP_Loss_plots = histFile->mkdir("EP_Loss_plots");
  EP_Loss_plots->cd();

  auto h_GEM1_delta_p =
      *rdf_def.Filter("GEM1_p > 0")
           .Histo1D({"GEM1 Delta p", "GEM1 Delta p; #Delta p [GeV/c]; Survival Fraction", n_bins, 0, 1}, "delta_p_GEM1",
                    "weight");
  h_GEM1_delta_p.Write();

  auto h_GEM2_delta_p =
      *rdf_def.Filter("GEM2_p > 0")
           .Histo1D({"GEM2 Delta p", "GEM2 Delta p; #Delta p [GeV/c]; Survival Fraction", n_bins, 0, 1}, "delta_p_GEM2",
                    "weight");
  h_GEM2_delta_p.Write();

  auto h_LAD_delta_p = *rdf_def.Filter("LAD_p > 0")
                            .Histo1D({"LAD Delta p", "LAD Delta p; #Delta p [GeV/c]; Survival Fraction", n_bins, 0, 1},
                                     "delta_p_LAD", "weight");
  h_LAD_delta_p.Write();

  auto h_GEM1_delta_E =
      *rdf_def.Filter("GEM1_p > 0")
           .Histo1D({"GEM1 Delta E", "GEM1 Delta E; #Delta E [GeV]; Survival Fraction", n_bins, 0, 0.1}, "delta_E_GEM1",
                    "weight");
  h_GEM1_delta_E.Write();

  auto h_GEM2_delta_E =
      *rdf_def.Filter("GEM2_p > 0")
           .Histo1D({"GEM2 Delta E", "GEM2 Delta E; #Delta E [GeV]; Survival Fraction", n_bins, 0, 0.1}, "delta_E_GEM2",
                    "weight");
  h_GEM2_delta_E.Write();

  auto h_LAD_delta_E = *rdf_def.Filter("LAD_p > 0")
                            .Histo1D({"LAD Delta E", "LAD Delta E; #Delta E [GeV]; Survival Fraction", n_bins, 0, 0.1},
                                     "delta_E_LAD", "weight");
  h_LAD_delta_E.Write();

  auto h_GEM1_p_delta_p =
      *rdf_def.Filter("GEM1_p > 0")
           .Histo2D({"GEM1 Momentum vs Delta p", "GEM1 Momentum vs Delta p; p [GeV/c]; #Delta p [GeV/c]", n_bins, 0, 1,
                     n_bins, 0, 1},
                    "ev.p", "delta_p_GEM1", "weight");
  h_GEM1_p_delta_p.Write();

  auto h_GEM2_p_delta_p =
      *rdf_def.Filter("GEM2_p > 0")
           .Histo2D({"GEM2 Momentum vs Delta p", "GEM2 Momentum vs Delta p; p [GeV/c]; #Delta p [GeV/c]", n_bins, 0, 1,
                     n_bins, 0, 1},
                    "ev.p", "delta_p_GEM2", "weight");
  h_GEM2_p_delta_p.Write();

  auto h_LAD_p_delta_p =
      *rdf_def.Filter("LAD_p > 0")
           .Histo2D({"LAD Momentum vs Delta p", "LAD Momentum vs Delta p; p [GeV/c]; #Delta p [GeV/c]", n_bins, 0, 1,
                     n_bins, 0, 1},
                    "ev.p", "delta_p_LAD", "weight");
  h_LAD_p_delta_p.Write();

  auto h_GEM1_p_delta_E =
      *rdf_def.Filter("GEM1_p > 0")
           .Histo2D({"GEM1 Momentum vs Delta E", "GEM1 Momentum vs Delta E; p [GeV/c]; #Delta E [GeV]", n_bins, 0, 1,
                     n_bins, 0, 0.1},
                    "ev.p", "delta_E_GEM1", "weight");
  h_GEM1_p_delta_E.Write();

  auto h_GEM2_p_delta_E =
      *rdf_def.Filter("GEM2_p > 0")
           .Histo2D({"GEM2 Momentum vs Delta E", "GEM2 Momentum vs Delta E; p [GeV/c]; #Delta E [GeV]", n_bins, 0, 1,
                     n_bins, 0, 0.1},
                    "ev.p", "delta_E_GEM2", "weight");
  h_GEM2_p_delta_E.Write();

  auto h_LAD_p_delta_E = *rdf_def.Filter("LAD_p > 0")
                              .Histo2D({"LAD Momentum vs Delta E", "LAD Momentum vs Delta E; p [GeV/c]; #Delta E [GeV]",
                                        n_bins, 0, 1, n_bins, 0, 0.1},
                                       "ev.p", "delta_E_LAD", "weight");
  h_LAD_p_delta_E.Write();

#pragma endregion EP_Loss_plots

#pragma region tof_plots

  TDirectory *tof_plots = histFile->mkdir("tof_plots");
  tof_plots->cd();

  auto rdf_tof =
      rdf_def.Define("GEM1_Beta", "get_Beta(hit.n,hit.pid,hit.det,hit.trid,hit.p,hit.x,hit.y,hit.z,hit.t,101)")
          .Define("GEM1_tof_p", "get_tof_Momentum(GEM1_Beta)")
          .Define("GEM2_Beta", "get_Beta(hit.n,hit.pid,hit.det,hit.trid,hit.p,hit.x,hit.y,hit.z,hit.t,102)")
          .Define("GEM2_tof_p", "get_tof_Momentum(GEM2_Beta)")
          .Define("LAD_Beta", "get_Beta(hit.n,hit.pid,hit.det,hit.trid,hit.p,hit.x,hit.y,hit.z,hit.t,501)")
          .Define("LAD_tof_p", "get_tof_Momentum(LAD_Beta)")
          .Define("GEM1_GEM2_Beta",
                  "get_Beta_2det(hit.n,hit.pid,hit.det,hit.trid,hit.p,hit.x,hit.y,hit.z,hit.t,101,102)")
          .Define("GEM1_GEM2_tof_p", "get_tof_Momentum(GEM1_GEM2_Beta)")
          .Define("GEM1_LAD_Beta",
                  "get_Beta_2det(hit.n,hit.pid,hit.det,hit.trid,hit.p,hit.x,hit.y,hit.z,hit.t,101,501)")
          .Define("GEM1_LAD_tof_p", "get_tof_Momentum(GEM1_LAD_Beta)")
          .Define("GEM2_LAD_Beta",
                  "get_Beta_2det(hit.n,hit.pid,hit.det,hit.trid,hit.p,hit.x,hit.y,hit.z,hit.t,102,501)")
          .Define("GEM2_LAD_tof_p", "get_tof_Momentum(GEM2_LAD_Beta)");

  auto h_GEM1_tof_p = *rdf_tof.Filter("LAD_p > 0")
                           .Define("delta_GEM1_tof_p", "ev.p-GEM1_tof_p")
                           .Histo2D({"GEM1 TOF P vs P", "GEM1 TOF P vs P; Initial p [GeV/c]; Initial - tof p [GeV/c]",
                                     n_bins, 0, 1, n_bins, -0.05, 0.3},
                                    "ev.p", "delta_GEM1_tof_p", "weight");
  h_GEM1_tof_p.Write();

  auto h_GEM2_tof_p = *rdf_tof.Filter("LAD_p > 0")
                           .Define("delta_GEM2_tof_p", "ev.p-GEM2_tof_p")
                           .Histo2D({"GEM2 TOF P vs P", "GEM2 TOF P vs P; Initial p [GeV/c]; Initial - tof p [GeV/c]",
                                     n_bins, 0, 1, n_bins, -0.05, 0.3},
                                    "ev.p", "delta_GEM2_tof_p", "weight");
  h_GEM2_tof_p.Write();

  auto h_LAD_tof_p = *rdf_tof.Filter("LAD_p > 0")
                          .Define("delta_LAD_tof_p", "ev.p-LAD_tof_p")
                          .Histo2D({"LAD TOF P vs P", "LAD TOF P vs P; Initial p [GeV/c]; Initial - tof p [GeV/c]",
                                    n_bins, 0, 1, n_bins, -0.05, 0.3},
                                   "ev.p", "delta_LAD_tof_p", "weight");
  h_LAD_tof_p.Write();

  auto h_GEM1_GEM2_tof_p =
      *rdf_tof.Filter("LAD_p > 0")
           .Define("delta_GEM1_GEM2_tof_p", "ev.p-GEM1_GEM2_tof_p")
           .Histo2D({"GEM1 GEM2 TOF P vs P", "GEM1 GEM2 TOF P vs P; Initial p [GeV/c]; Initial - tof p [GeV/c]", n_bins,
                     0, 1, n_bins, -0.05, 0.3},
                    "ev.p", "delta_GEM1_GEM2_tof_p", "weight");
  h_GEM1_GEM2_tof_p.Write();

  auto h_GEM1_LAD_tof_p =
      *rdf_tof.Filter("LAD_p > 0")
           .Define("delta_GEM1_LAD_tof_p", "ev.p-GEM1_LAD_tof_p")
           .Histo2D({"GEM1 LAD TOF P vs P", "GEM1 LAD TOF P vs P; Initial p [GeV/c]; Initial - tof p [GeV/c]", n_bins,
                     0, 1, n_bins, -0.05, 0.3},
                    "ev.p", "delta_GEM1_LAD_tof_p", "weight");
  h_GEM1_LAD_tof_p.Write();

  auto h_GEM2_LAD_tof_p =
      *rdf_tof.Filter("LAD_p > 0")
           .Define("delta_GEM2_LAD_tof_p", "ev.p-GEM2_LAD_tof_p")
           .Histo2D({"GEM2 LAD TOF P vs P", "GEM2 LAD TOF P vs P; Initial p [GeV/c]; Initial - tof p [GeV/c]", n_bins,
                     0, 1, n_bins, -0.05, 0.3},
                    "ev.p", "delta_GEM2_LAD_tof_p", "weight");
  h_GEM2_LAD_tof_p.Write();

#pragma endregion tof_plots

  return;
}