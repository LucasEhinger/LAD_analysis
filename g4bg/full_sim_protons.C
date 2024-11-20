// #include "analysis.h"
#include <iostream>
#include <TFile.h>
#include <TChain.h>
#include <ROOT/RDataFrame.hxx>

using namespace std;
using namespace ROOT;

double get_Momentum(int n_hits, RVec<Int_t> pdgid, RVec<Int_t> detid, RVec<Int_t> trackid, RVec<double> p, double p_init, int gem_ID)
{
  for (int i = 0; i < n_hits; i++)
  {
    if (detid[i] == gem_ID && pdgid[i] == 2212)
    {
      return p[i];
    }
  }
  return 0;
}

void full_sim_protons(string inFileName, string outHistName, string treeName)
{
  int n_bins = 20;

  TFile *histFile = new TFile(outHistName.c_str(), "RECREATE");
  histFile->cd();

  TChain chain(treeName.c_str());
  chain.Add(inFileName.c_str());
  RDataFrame rdf_raw(chain);
  auto rdf_def = rdf_raw.Define("GEM1_p", "get_Momentum(hit.n,hit.pid,hit.det,hit.trid,hit.p,ev.p,101)")
                     .Define("GEM2_p", "get_Momentum(hit.n,hit.pid,hit.det,hit.trid,hit.p,ev.p,102)")
                     .Define("LAD_p", "get_Momentum(hit.n,hit.pid,hit.det,hit.trid,hit.p,ev.p,501)")
                     .Define("theta", "ev.theta*180/3.14159");

  auto h_GEM1_p = *rdf_def.Filter("GEM1_p > 0").Histo1D({"GEM1 Momentum", "GEM1 Momentum; p [GeV/c]; Counts", n_bins, 0, 1}, "GEM1_p");
  h_GEM1_p.Write();

  auto h_GEM2_p = *rdf_def.Filter("GEM2_p > 0").Histo1D({"GEM2 Momentum", "GEM2 Momentum; p [GeV/c]; Counts", n_bins, 0, 1}, "GEM2_p");
  h_GEM2_p.Write();

  auto h_LAD_p = *rdf_def.Filter("LAD_p > 0").Histo1D({"LAD Momentum", "LAD Momentum; p [GeV/c]; Counts", n_bins, 0, 1}, "LAD_p");
  h_LAD_p.Write();

  auto h_GEM1_p_original = *rdf_def.Filter("GEM1_p > 0").Histo1D({"GEM1 Original Momentum", "GEM1 Original Momentum; p [GeV/c]; Counts", n_bins, 0, 1}, "ev.p");
  h_GEM1_p_original.Write();

  auto h_GEM2_p_original = *rdf_def.Filter("GEM2_p > 0").Histo1D({"GEM2 Original Momentum", "GEM2 Original Momentum; p [GeV/c]; Counts", n_bins, 0, 1}, "ev.p");
  h_GEM2_p_original.Write();

  auto h_LAD_p_original = *rdf_def.Filter("LAD_p > 0").Histo1D({"LAD Original Momentum", "LAD Original Momentum; p [GeV/c]; Counts", n_bins, 0, 1}, "ev.p");
  h_LAD_p_original.Write();

  histFile->Close();
}