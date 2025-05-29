#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TMinuit.h>
#include <TTree.h>
#include <TVector3.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace std;

const double min_projy = -5.0;
const double max_projy = 0.0;
const double max_d0    = 30.0;
// const double min_projz = -13.0;
// const double max_projz = 13.0;
const double projz_sigma       = 5.0;
const int nFixedz              = 3;
const double target_z[nFixedz] = {-10.0, 0.0, 10.0}; // Fixed z positions for the planes

const double TDC2NS = 0.09766; // ns per TDC channel

struct hist_range {
  double min;
  double max;
  int nbins;

  hist_range(double min_val, double max_val, int bins) : min(min_val), max(max_val), nbins(bins) {}
};

const hist_range TDC_TIME_RAW(0, 4000, 200);
const hist_range TARGED_POS(-20, 20, 100);
const hist_range LAD_POS(-100, 100, 100);

int run_number = 300000;

const int MAX_DATA                = 10000;
const int maxTracks               = 30;
const int nPlanes                 = 5;
const string plane_names[nPlanes] = {"000", "001", "100", "101", "200"};

const double DCA_XZ_MAX = 100.0; // Maximum DCA in XZ plane
const double DCA_YZ_MAX = 200.0; // Maximum DCA in YZ plane

const double plane_theta[nPlanes] = {150.0, 150.0, 127.0, 127.0, 104.0}; // Angle in degrees
const double plane_r[nPlanes]     = {615.0, 655.6, 523.0, 563.6, 615.0}; // Radius of the second point

const bool use_projz = true; // Fix z position to GEM projz

const double dx_min = -30.0;
const double dx_max = 30.0;
const int dx_NBINS  = 60;
const double dy_min = -30.0;
const double dy_max = 30.0;
const int dy_NBINS  = 60;

TVector3 LinePlaneIntersection(const TVector3 &p1, const TVector3 &p2, const TVector3 &p3, const TVector3 &l1,
                               const TVector3 &l2) {
  // Define the plane normal
  TVector3 planeNormal = (p2 - p1).Cross(p3 - p1);
  planeNormal          = planeNormal.Unit();

  // Line direction
  TVector3 lineDir = l2 - l1;

  // Check if the line is parallel to the plane
  double denom = planeNormal.Dot(lineDir);
  if (fabs(denom) < 1e-6) {
    std::cerr << "The line is parallel to the plane and does not intersect." << std::endl;
    return TVector3(0, 0, 0); // Return a zero vector if the line is parallel to the plane
  }

  // Calculate the intersection point
  double t              = planeNormal.Dot(p1 - l1) / denom;
  TVector3 intersection = l1 + t * lineDir;

  return intersection;
}

TVector3 GetHodoHitPosition(const int paddle, const int plane) {

  double dTrans = (5 - paddle) * 22.0;
  // Define the plane using three points
  TVector3 p_hit(plane_r[plane] * cos((plane_theta[plane] - 90) * TMath::DegToRad()), 0,
                 -plane_r[plane] * sin((plane_theta[plane] - 90) * TMath::DegToRad()));
  p_hit = p_hit + TVector3(-dTrans * cos((180 - plane_theta[plane]) * TMath::DegToRad()), 0,
                           -dTrans * sin((180 - plane_theta[plane]) * TMath::DegToRad()));

  return p_hit;
}

TVector3 GetGEMHitPosition(const double x_loc, const double y_loc, const double dx, const double dy, const double r,
                           const double gem_theta) {
  // Define the plane using three points
  TVector3 p_hit(r * cos((gem_theta - 90) * TMath::DegToRad()), 0, -r * sin((gem_theta - 90) * TMath::DegToRad()));
  p_hit = p_hit + TVector3(-(x_loc + dx) * cos((180 - gem_theta) * TMath::DegToRad()), dy,
                           -(x_loc + dx) * sin((180 - gem_theta) * TMath::DegToRad()));

  return p_hit;
}

double GetProjZ(const TVector3 &gem_hit0, const TVector3 &gem_hit1, const int x_targ) {

  // Calculate the direction vector of the line
  TVector3 direction = gem_hit1 - gem_hit0;

  // Check if the direction vector is valid
  if (fabs(direction.X()) < 1e-6) {
    std::cerr << "Warning: The line is parallel to the x-axis and cannot be projected. Setting proj_z to 100."
              << std::endl;
    return 100.0;
  }

  // Calculate the parameter t for the line equation
  double t = (x_targ - gem_hit0.X()) / direction.X();

  // Calculate the z-coordinate at the projected point
  double proj_z = gem_hit0.Z() + t * direction.Z();

  return proj_z;
}

double get_dca_xz_hodo(int hodoPlane, int hodoPaddle, TVector3 p_gem_hit0, TVector3 p_gem_hit1) {
  // Get the hodo hit position
  TVector3 p_hit = GetHodoHitPosition(hodoPaddle, hodoPlane);
  // Get the GEM hit position

  TVector3 lineDir   = p_gem_hit1 - p_gem_hit0;
  TVector3 linePoint = p_gem_hit0;

  // Vector from the line point to p_hit
  TVector3 vecToPoint = p_hit - linePoint;

  // Calculate the distance of closest approach in the x-z plane
  // Project vecToPoint onto the line direction
  double t = vecToPoint.Dot(lineDir) / lineDir.Mag2();

  // Closest point on the line to p_hit
  TVector3 closestPoint = linePoint + t * lineDir;

  // Calculate the distance of closest approach in the x-z plane
  double dca_xz = sqrt(pow(closestPoint.X() - p_hit.X(), 2) + pow(closestPoint.Z() - p_hit.Z(), 2));

  return dca_xz;
}

double get_dca_y_hodo(int hodoPlane, int hodoPaddle, TVector3 p_gem_hit0, TVector3 p_gem_hit1) {
  // Get the hodo hit position
  TVector3 p_hit = GetHodoHitPosition(hodoPaddle, hodoPlane);
  // Get the GEM hit position

  TVector3 lineDir   = p_gem_hit1 - p_gem_hit0;
  TVector3 linePoint = p_gem_hit0;

  // Vector from the line point to p_hit
  TVector3 vecToPoint = p_hit - linePoint;

  // Calculate the distance of closest approach in the y direction
  double dca_y = vecToPoint.Y();

  return dca_y;
}

void gem_lad_timing_full(int run_num) {

  TString fileName = Form("/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/"
                          "LAD_COIN_%d_0_0_-1.root", run_num);

                     //  "LAD_COIN_22282_-1_inverted.root";
                     //  "LAD_COIN_22282_-1_500trks_good_timing.root";
                     //  "LAD_COIN_22383_0_0_500002.root";

  TString outputFileName = Form("files/gem_window/gem_window_%d_-1_P.root", run_num);
  // Open the ROOT file
  TFile *file = TFile::Open(fileName);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open the ROOT file!" << std::endl;
    return;
  }

  // Get the TTree
  TTree *T = dynamic_cast<TTree *>(file->Get("T"));
  if (!T) {
    std::cerr << "Error: Cannot find the TTree named 'T'!" << std::endl;
    file->Close();
    return;
  }
  // Define arrays to hold the data
  Double_t trk_d0[MAX_DATA], trk_d0_good[MAX_DATA];
  Double_t trk_projz[MAX_DATA], trk_projy[MAX_DATA];
  Double_t trk_t[MAX_DATA], trk_dt[MAX_DATA];
  Double_t trk_x[2][MAX_DATA], trk_y[2][MAX_DATA], trk_z[2][MAX_DATA];
  Double_t trk_x_local[2][MAX_DATA], trk_y_local[2][MAX_DATA];
  Double_t tdc_time_btm[nPlanes][MAX_DATA], tdc_time_top[nPlanes][MAX_DATA];
  Double_t tdc_counter_btm[nPlanes][MAX_DATA], tdc_counter_top[nPlanes][MAX_DATA];
  Int_t nTracks, nTdcTopHits[nPlanes], nTdcBtmHits[nPlanes];
  Double_t vertex_x, vertex_y, vertex_z;
  char spect_prefix = 'P'; // Default to 'H', can be changed to 'P' if needed

  T->SetBranchAddress(Form("Ndata.%c.gem.trk.d0", spect_prefix), &nTracks);
  T->SetBranchAddress(Form("%c.gem.trk.d0", spect_prefix), &trk_d0);
  T->SetBranchAddress(Form("%c.gem.trk.d0_good", spect_prefix), &trk_d0_good);
  T->SetBranchAddress(Form("%c.gem.trk.projz", spect_prefix), &trk_projz);
  T->SetBranchAddress(Form("%c.gem.trk.projy", spect_prefix), &trk_projy);
  T->SetBranchAddress(Form("%c.gem.trk.x1", spect_prefix), &trk_x[0]);
  T->SetBranchAddress(Form("%c.gem.trk.y1", spect_prefix), &trk_y[0]);
  T->SetBranchAddress(Form("%c.gem.trk.z1", spect_prefix), &trk_z[0]);
  T->SetBranchAddress(Form("%c.gem.trk.x2", spect_prefix), &trk_x[1]);
  T->SetBranchAddress(Form("%c.gem.trk.y2", spect_prefix), &trk_y[1]);
  T->SetBranchAddress(Form("%c.gem.trk.z2", spect_prefix), &trk_z[1]);
  T->SetBranchAddress(Form("%c.gem.trk.x1_local", spect_prefix), &trk_x_local[0]);
  T->SetBranchAddress(Form("%c.gem.trk.y1_local", spect_prefix), &trk_y_local[0]);
  T->SetBranchAddress(Form("%c.gem.trk.x2_local", spect_prefix), &trk_x_local[1]);
  T->SetBranchAddress(Form("%c.gem.trk.y2_local", spect_prefix), &trk_y_local[1]);
  T->SetBranchAddress(Form("%c.gem.trk.t", spect_prefix), &trk_t);
  T->SetBranchAddress(Form("%c.gem.trk.dt", spect_prefix), &trk_dt);
  T->SetBranchAddress(Form("%c.react.x", spect_prefix), &vertex_x);
  T->SetBranchAddress(Form("%c.react.y", spect_prefix), &vertex_y);
  T->SetBranchAddress(Form("%c.react.z", spect_prefix), &vertex_z);

  for (int i = 0; i < nPlanes; ++i) {
    T->SetBranchAddress(Form("%c.ladhod.%s.BtmTdcTimeRaw", spect_prefix, plane_names[i].c_str()), &tdc_time_btm[i]);
    T->SetBranchAddress(Form("%c.ladhod.%s.BtmTdcCounter", spect_prefix, plane_names[i].c_str()), &tdc_counter_btm[i]);
    T->SetBranchAddress(Form("%c.ladhod.%s.TopTdcTimeRaw", spect_prefix, plane_names[i].c_str()), &tdc_time_top[i]);
    T->SetBranchAddress(Form("%c.ladhod.%s.TopTdcCounter", spect_prefix, plane_names[i].c_str()), &tdc_counter_top[i]);
    T->SetBranchAddress(Form("Ndata.%c.ladhod.%s.BtmTdcTime", spect_prefix, plane_names[i].c_str()), &nTdcBtmHits[i]);
    T->SetBranchAddress(Form("Ndata.%c.ladhod.%s.TopTdcTime", spect_prefix, plane_names[i].c_str()), &nTdcTopHits[i]);
  }

  //////////////////////////////////////////////////////////////////
  // Create histograms for dx and dy
  TH1D *h_time_wTrack_top[nPlanes];
  TH1D *h_time_wTrack_btm[nPlanes];
  TH1D *h_all_times_top[nPlanes];
  TH1D *h_all_times_btm[nPlanes];
  TH1D *h_lad_d0_xz_top[nPlanes];
  TH1D *h_lad_d0_xz_btm[nPlanes];
  TH1D *h_lad_d0_yz_top[nPlanes];
  TH1D *h_lad_d0_yz_btm[nPlanes];
  TH1D *h_projz_top[nPlanes];
  TH1D *h_projz_btm[nPlanes];
  TH1D *h_projy_top[nPlanes];
  TH1D *h_projy_btm[nPlanes];

  for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
    h_time_wTrack_top[i_plane] = new TH1D(Form("h_time_top_%d", i_plane), Form("Top TDC Time Plane %d", i_plane),
                                          TDC_TIME_RAW.nbins, TDC_TIME_RAW.min, TDC_TIME_RAW.max);
    h_time_wTrack_btm[i_plane] = new TH1D(Form("h_time_btm_%d", i_plane), Form("Bottom TDC Time Plane %d", i_plane),
                                          TDC_TIME_RAW.nbins, TDC_TIME_RAW.min, TDC_TIME_RAW.max);
    h_all_times_top[i_plane] = new TH1D(Form("h_all_times_top_%d", i_plane), Form("Top All TDC Time Plane %d", i_plane),
                                        TDC_TIME_RAW.nbins, TDC_TIME_RAW.min, TDC_TIME_RAW.max);
    h_all_times_btm[i_plane] =
        new TH1D(Form("h_all_times_btm_%d", i_plane), Form("Bottom All TDC Time Plane %d", i_plane), TDC_TIME_RAW.nbins,
                 TDC_TIME_RAW.min, TDC_TIME_RAW.max);
    h_lad_d0_xz_top[i_plane] = new TH1D(Form("h_lad_d0_xz_top_%d", i_plane), Form("Top DCA XZ Plane %d", i_plane),
                                        LAD_POS.nbins, LAD_POS.min, LAD_POS.max);
    h_lad_d0_xz_btm[i_plane] = new TH1D(Form("h_lad_d0_xz_btm_%d", i_plane), Form("Bottom DCA XZ Plane %d", i_plane),
                                        LAD_POS.nbins, LAD_POS.min, LAD_POS.max);
    h_lad_d0_yz_top[i_plane] = new TH1D(Form("h_lad_d0_yz_top_%d", i_plane), Form("Top DCA YZ Plane %d", i_plane),
                                        LAD_POS.nbins, LAD_POS.min, LAD_POS.max);
    h_lad_d0_yz_btm[i_plane] = new TH1D(Form("h_lad_d0_yz_btm_%d", i_plane), Form("Bottom DCA YZ Plane %d", i_plane),
                                        LAD_POS.nbins, LAD_POS.min, LAD_POS.max);
    h_projz_top[i_plane]     = new TH1D(Form("h_projz_top_%d", i_plane), Form("Top ProjZ Plane %d", i_plane),
                                        TARGED_POS.nbins, TARGED_POS.min, TARGED_POS.max);
    h_projz_btm[i_plane]     = new TH1D(Form("h_projz_btm_%d", i_plane), Form("Bottom ProjZ Plane %d", i_plane),
                                        TARGED_POS.nbins, TARGED_POS.min, TARGED_POS.max);
    h_projy_top[i_plane]     = new TH1D(Form("h_projy_top_%d", i_plane), Form("Top ProjY Plane %d", i_plane),
                                        TARGED_POS.nbins, TARGED_POS.min, TARGED_POS.max);
    h_projy_btm[i_plane]     = new TH1D(Form("h_projy_btm_%d", i_plane), Form("Bottom ProjY Plane %d", i_plane),
                                        TARGED_POS.nbins, TARGED_POS.min, TARGED_POS.max);
  }
  ///////////////////////////////////////////////////////////////////

  // Loop through the tree entries
  Long64_t nEntries = T->GetEntries();
  // nEntries          = 10000; // For testing purposes, limit to 1000 entries
  for (Long64_t i = 0; i < nEntries; ++i) {
    T->GetEntry(i);
    // Skip processing if there are too many tracks
    for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
      for (int i_top = 0; i_top < nTdcTopHits[i_plane]; ++i_top) {
        h_all_times_top[i_plane]->Fill(tdc_time_top[i_plane][i_top] * TDC2NS);
      }
      for (int i_btm = 0; i_btm < nTdcBtmHits[i_plane]; ++i_btm) {
        h_all_times_btm[i_plane]->Fill(tdc_time_btm[i_plane][i_btm] * TDC2NS);
      }
    }

    for (int i_trk = 0; i_trk < nTracks; ++i_trk) {
      // Skip if the track if it doesn't point back to the target
      bool is_on_target = false;
      for (int i_foil = 0; i_foil < nFixedz; ++i_foil) {
        if (trk_projz[i_trk] > target_z[i_foil] - projz_sigma && trk_projz[i_trk] < target_z[i_foil] + projz_sigma) {
          is_on_target = true;
          break;
        }
      }
      if (!is_on_target) {
        continue;
      }
      if (trk_projy[i_trk] < min_projy || trk_projy[i_trk] > max_projy) {
        continue;
      }
      if (trk_d0[i_trk] > max_d0) {
        continue;
      }
      if (abs(trk_dt[i_trk]) > 15.0) {
        continue;
      }

      // Loop through all TDC hits
      for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {

        for (int i_top = 0; i_top < nTdcTopHits[i_plane]; ++i_top) {
          TVector3 p_gemhit_0(trk_x[0][i_trk], trk_y[0][i_trk], trk_z[0][i_trk]);
          TVector3 p_gemhit_1(trk_x[1][i_trk], trk_y[1][i_trk], trk_z[1][i_trk]);
          int paddle = tdc_counter_top[i_plane][i_top]-1;
          double dca_xz_top = get_dca_xz_hodo(i_plane, paddle, p_gemhit_0, p_gemhit_1);
          double dca_y_top  = get_dca_y_hodo(i_plane, paddle, p_gemhit_0, p_gemhit_1);

          if (dca_xz_top < DCA_XZ_MAX) {
            h_lad_d0_xz_top[i_plane]->Fill(dca_xz_top);
            h_lad_d0_yz_top[i_plane]->Fill(dca_y_top);
            h_projz_top[i_plane]->Fill(trk_projz[i_trk]);
            h_projy_top[i_plane]->Fill(trk_projy[i_trk]);
            h_time_wTrack_top[i_plane]->Fill(tdc_time_top[i_plane][i_top] * TDC2NS);
          }
        }
        for (int i_btm = 0; i_btm < nTdcBtmHits[i_plane]; ++i_btm) {
          TVector3 p_gemhit_0(trk_x[0][i_trk], trk_y[0][i_trk], trk_z[0][i_trk]);
          TVector3 p_gemhit_1(trk_x[1][i_trk], trk_y[1][i_trk], trk_z[1][i_trk]);
          int paddle = tdc_counter_btm[i_plane][i_btm]-1;
          double dca_xz_btm = get_dca_xz_hodo(i_plane, paddle, p_gemhit_0, p_gemhit_1);
          double dca_y_btm  = get_dca_y_hodo(i_plane, paddle, p_gemhit_0, p_gemhit_1);

          if (dca_xz_btm < DCA_XZ_MAX) {
            h_lad_d0_xz_btm[i_plane]->Fill(dca_xz_btm);
            h_lad_d0_yz_btm[i_plane]->Fill(dca_y_btm);
            h_projz_btm[i_plane]->Fill(trk_projz[i_trk]);
            h_projy_btm[i_plane]->Fill(trk_projy[i_trk]);
            h_time_wTrack_btm[i_plane]->Fill(tdc_time_btm[i_plane][i_btm] * TDC2NS);
          }
        }
      }
    }
    // Print the status as a percentage
    if (i % (nEntries / 100) == 0) {
      std::cout << "\rProcessing: " << int(i * 100.0 / nEntries) << "% completed." << std::flush;
    }

  } // End Event Loop
  std::cout << "\nProcessing completed." << std::endl;

  // Create the output file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Error: Cannot create the output ROOT file!" << std::endl;
    return;
  }
  // Write the histograms and canvases to the output file
  outputFile->cd();

  // Create directories in the output file
  TDirectory *dir_time_wTrack = outputFile->mkdir("time_wTrack");
  TDirectory *dir_all_times   = outputFile->mkdir("all_times");
  TDirectory *dir_time_ratio  = outputFile->mkdir("time_ratio");
  TDirectory *dir_xz          = outputFile->mkdir("xz");
  TDirectory *dir_y           = outputFile->mkdir("y");
  TDirectory *dir_projy       = outputFile->mkdir("projy");
  TDirectory *dir_projz       = outputFile->mkdir("projz");

  // Create histograms for all planes added together
  TH1D *h_time_wTrack_top_all_planes = (TH1D *)h_time_wTrack_top[0]->Clone("h_all_time_top");
  h_time_wTrack_top_all_planes->Reset();
  h_time_wTrack_top_all_planes->SetTitle("Top TDC Time With Track (All Planes)");

  TH1D *h_time_wTrack_btm_all_planes = (TH1D *)h_time_wTrack_btm[0]->Clone("h_all_time_btm");
  h_time_wTrack_btm_all_planes->Reset();
  h_time_wTrack_btm_all_planes->SetTitle("Bottom TDC Time With Track (All Planes)");

  TH1D *h_all_times_top_all_planes = (TH1D *)h_all_times_top[0]->Clone("h_all_times_top");
  h_all_times_top_all_planes->Reset();
  h_all_times_top_all_planes->SetTitle("Top TDC Time (All Planes)");

  TH1D *h_all_times_btm_all_planes = (TH1D *)h_all_times_btm[0]->Clone("h_all_times_btm");
  h_all_times_btm_all_planes->Reset();
  h_all_times_btm_all_planes->SetTitle("Bottom TDC Time (All Planes)");

  TH1D *h_all_xz_top = (TH1D *)h_lad_d0_xz_top[0]->Clone("h_all_xz_top");
  h_all_xz_top->Reset();
  h_all_xz_top->SetTitle("All Top Planes Combined (XZ)");

  TH1D *h_all_xz_btm = (TH1D *)h_lad_d0_xz_btm[0]->Clone("h_all_xz_btm");
  h_all_xz_btm->Reset();
  h_all_xz_btm->SetTitle("All Bottom Planes Combined (XZ)");

  TH1D *h_all_y_top = (TH1D *)h_lad_d0_yz_top[0]->Clone("h_all_y_top");
  h_all_y_top->Reset();
  h_all_y_top->SetTitle("All Top Planes Combined (Y)");

  TH1D *h_all_y_btm = (TH1D *)h_lad_d0_yz_btm[0]->Clone("h_all_y_btm");
  h_all_y_btm->Reset();
  h_all_y_btm->SetTitle("All Bottom Planes Combined (Y)");

  TH1D *h_all_projy_top = (TH1D *)h_projy_top[0]->Clone("h_all_projy_top");
  h_all_projy_top->Reset();
  h_all_projy_top->SetTitle("All Top Planes Combined (ProjY)");

  TH1D *h_all_projy_btm = (TH1D *)h_projy_btm[0]->Clone("h_all_projy_btm");
  h_all_projy_btm->Reset();
  h_all_projy_btm->SetTitle("All Bottom Planes Combined (ProjY)");

  TH1D *h_all_projz_top = (TH1D *)h_projz_top[0]->Clone("h_all_projz_top");
  h_all_projz_top->Reset();
  h_all_projz_top->SetTitle("All Top Planes Combined (ProjZ)");

  TH1D *h_all_projz_btm = (TH1D *)h_projz_btm[0]->Clone("h_all_projz_btm");
  h_all_projz_btm->Reset();
  h_all_projz_btm->SetTitle("All Bottom Planes Combined (ProjZ)");

  // Combine histograms for all planes
  for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
    h_time_wTrack_top_all_planes->Add(h_time_wTrack_top[i_plane]);
    h_time_wTrack_btm_all_planes->Add(h_time_wTrack_btm[i_plane]);
    h_all_times_top_all_planes->Add(h_all_times_top[i_plane]);
    h_all_times_btm_all_planes->Add(h_all_times_btm[i_plane]);
    h_all_xz_top->Add(h_lad_d0_xz_top[i_plane]);
    h_all_xz_btm->Add(h_lad_d0_xz_btm[i_plane]);
    h_all_y_top->Add(h_lad_d0_yz_top[i_plane]);
    h_all_y_btm->Add(h_lad_d0_yz_btm[i_plane]);
    h_all_projy_top->Add(h_projy_top[i_plane]);
    h_all_projy_btm->Add(h_projy_btm[i_plane]);
    h_all_projz_top->Add(h_projz_top[i_plane]);
    h_all_projz_btm->Add(h_projz_btm[i_plane]);
  }

  // Write histograms to their respective directories
  dir_time_wTrack->cd();
  for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
    h_time_wTrack_top[i_plane]->Write();
    h_time_wTrack_btm[i_plane]->Write();
  }
  h_time_wTrack_top_all_planes->Write();
  h_time_wTrack_btm_all_planes->Write();

  dir_all_times->cd();
  for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
    h_all_times_top[i_plane]->Write();
    h_all_times_btm[i_plane]->Write();
  }
  h_all_times_top_all_planes->Write();
  h_all_times_btm_all_planes->Write();

  dir_time_ratio->cd();
  for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
    TH1D *h_time_ratio_top = (TH1D *)h_time_wTrack_top[i_plane]->Clone(Form("h_time_ratio_top_%d", i_plane));
    h_time_ratio_top->Divide(h_all_times_top[i_plane]);
    h_time_ratio_top->Write();

    TH1D *h_time_ratio_btm = (TH1D *)h_time_wTrack_btm[i_plane]->Clone(Form("h_time_ratio_btm_%d", i_plane));
    h_time_ratio_btm->Divide(h_all_times_btm[i_plane]);
    h_time_ratio_btm->Write();
  }
  TH1D *h_time_ratio_top_all_planes = (TH1D *)h_time_wTrack_top_all_planes->Clone("h_time_ratio_top_all_planes");
  h_time_ratio_top_all_planes->Divide(h_all_times_top_all_planes);
  h_time_ratio_top_all_planes->SetTitle("Bottom TDC Time With Track / All Hits (All Planes)");
  h_time_ratio_top_all_planes->Write();
  TH1D *h_time_ratio_btm_all_planes = (TH1D *)h_time_wTrack_btm_all_planes->Clone("h_time_ratio_btm_all_planes");
  h_time_ratio_btm_all_planes->Divide(h_all_times_btm_all_planes);
  h_time_ratio_btm_all_planes->SetTitle("Bottom TDC Time With Track / All Hits (All Planes)");
  h_time_ratio_btm_all_planes->Write();

  dir_xz->cd();
  for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
    h_lad_d0_xz_top[i_plane]->Write();
    h_lad_d0_xz_btm[i_plane]->Write();
  }
  h_all_xz_top->Write();
  h_all_xz_btm->Write();

  dir_y->cd();
  for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
    h_lad_d0_yz_top[i_plane]->Write();
    h_lad_d0_yz_btm[i_plane]->Write();
  }
  h_all_y_top->Write();
  h_all_y_btm->Write();

  dir_projy->cd();
  for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
    h_projy_top[i_plane]->Write();
    h_projy_btm[i_plane]->Write();
  }
  h_all_projy_top->Write();
  h_all_projy_btm->Write();

  dir_projz->cd();
  for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
    h_projz_top[i_plane]->Write();
    h_projz_btm[i_plane]->Write();
  }
  h_all_projz_top->Write();
  h_all_projz_btm->Write();
  // Close the files
  file->Close();
  outputFile->Close();
  std::cout << "Output file created: " << outputFileName << std::endl;
  return;
}