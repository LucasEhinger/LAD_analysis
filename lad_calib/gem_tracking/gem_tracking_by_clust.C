#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMinuit.h>
#include <TROOT.h>
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
const hist_range GEM_ADC(0, 10000, 100);

const hist_range lad_time(1600, 2100, 100);
const hist_range time_peak(1770, 1790, 100);
const hist_range time_before(1700, 1750, 100);
const hist_range time_after(1850, 1950, 100);

int run_number = 300000;

const int MAX_DATA                = 10000;
const int maxTracks               = 30;
const int nPlanes                 = 5;
const string plane_names[nPlanes] = {"000", "001", "100", "101", "200"};
const int nPaddles                = 11;

const double DCA_XZ_MAX = 100.0; // Maximum DCA in XZ plane
const double DCA_YZ_MAX = 200.0; // Maximum DCA in YZ plane

const double plane_theta[nPlanes] = {150.0, 150.0, 127.0, 127.0, 104.0}; // Angle in degrees
const double plane_r[nPlanes]     = {615.0, 655.6, 523.0, 563.6, 615.0}; // Radius of the second point

const bool use_projz = true; // Fix z position to GEM projz

const double dx_min = -70.0;
const double dx_max = 70.0;
const int dx_NBINS  = 40;
const double dy_min = -40.0;
const double dy_max = 40.0;
const int dy_NBINS  = 40;

const double MAX_DX_CLUST = 5.0; // Maximum dx for cluster hits
const double MAX_DY_CLUST = 5.0; // Maximum dy for cluster hits

double gem_theta = 127.0;            // Angle in degrees
double gem_phi   = 0.0;              // Angle in degrees
double gem_r[2]  = {77.571, 95.571}; // Radius of the GEM's
double gem_dx[2] = {0.0, 0.0};       // GEM dx offsets
double gem_dy[2] = {0.0, 0.0};       // GEM dy offsets

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

void gem_tracking_by_clust() {
  // Set ROOT to batch mode to avoid opening canvases
  gROOT->SetBatch(kTRUE);
  std::vector<TString> fileNames = {
      Form("/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22572_0_6_2000000.root")
      // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22565_0_0_-1.root"
      // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_296_0_0_-1.root"
      // "LAD_COIN_22282_-1_inverted.root",
      // "LAD_COIN_22282_-1_500trks_good_timing.root",
      // "LAD_COIN_22383_0_0_500002.root"
      // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22565_0_0_-1.root",
      // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22565_0_0_-1.root",
      // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22565_0_0_-1.root",
      // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22565_0_0_-1.root"

  };
  // Form("/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22572_0_6_2000000.root");
  // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_22565_0_0_-1.root";
  // "/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/LAD_COIN_296_0_0_-1.root";

  //  "LAD_COIN_22282_-1_inverted.root";
  //  "LAD_COIN_22282_-1_500trks_good_timing.root";
  //  "LAD_COIN_22383_0_0_500002.root";

  TFile *file            = new TFile(fileNames[0]);
  TString outputFileName = Form("files/tracking_gem/tracking_gem_22572_-1_P.root");
  // Create a TChain instead of a TTree
  // TChain *T = new TChain("T");
  TTree *T = (TTree *)file->Get("T");
  // Add files to the TChain
  // for (const auto &fileName : fileNames) {
  //   T->Add(fileName);
  // }
  if (T->GetEntries() == 0) {
    std::cerr << "Error: Cannot find any entries in the TChain!" << std::endl;
    return;
  }
  // Define arrays to hold the data
  Double_t clust_pos[MAX_DATA], clust_layer[MAX_DATA], clust_axis[MAX_DATA], clust_adc[MAX_DATA]; // 0:y. 1:x
  Int_t nData_clust;
  Double_t sp_X[MAX_DATA], sp_Y[MAX_DATA], sp_Z[MAX_DATA], sp_adc[MAX_DATA];
  Double_t sp_layer[MAX_DATA];
  Int_t nData_sp;

  Double_t trk_d0[MAX_DATA], trk_d0_good[MAX_DATA];
  Double_t trk_projz[MAX_DATA], trk_projy[MAX_DATA];
  Double_t trk_t[MAX_DATA], trk_dt[MAX_DATA];
  Double_t trk_x[2][MAX_DATA], trk_y[2][MAX_DATA], trk_z[2][MAX_DATA];
  Double_t trk_x_local[2][MAX_DATA], trk_y_local[2][MAX_DATA];
  Double_t tdc_time_btm[nPlanes][MAX_DATA], tdc_time_top[nPlanes][MAX_DATA];
  Double_t tdc_counter_btm[nPlanes][MAX_DATA], tdc_counter_top[nPlanes][MAX_DATA];
  Int_t nTracks, nTdcTopHits[nPlanes], nTdcBtmHits[nPlanes];
  Double_t vertex_x, vertex_y, vertex_z;

  Double_t time_avg[nPlanes][MAX_DATA], time_avg_paddle[nPlanes][MAX_DATA], time_avg_ypos[nPlanes][MAX_DATA];
  Int_t nData_hodo[nPlanes];
  char spect_prefix = 'P'; // Default to 'H', can be changed to 'P' if needed

  T->SetBranchAddress(Form("Ndata.%c.gem.clust.pos", spect_prefix), &nData_clust);
  T->SetBranchAddress(Form("%c.gem.clust.pos", spect_prefix), &clust_pos);
  T->SetBranchAddress(Form("%c.gem.clust.layer", spect_prefix), &clust_layer);
  T->SetBranchAddress(Form("%c.gem.clust.axis", spect_prefix), &clust_axis);
  T->SetBranchAddress(Form("%c.gem.clust.adc", spect_prefix), &clust_adc);
  T->SetBranchAddress(Form("Ndata.%c.gem.sp.posX", spect_prefix), &nData_sp);
  T->SetBranchAddress(Form("%c.gem.sp.posX", spect_prefix), &sp_X);
  T->SetBranchAddress(Form("%c.gem.sp.posY", spect_prefix), &sp_Y);
  T->SetBranchAddress(Form("%c.gem.sp.posZ", spect_prefix), &sp_Z);
  T->SetBranchAddress(Form("%c.gem.sp.layer", spect_prefix), &sp_layer);
  T->SetBranchAddress(Form("%c.gem.sp.adc", spect_prefix), &sp_adc);
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
    T->SetBranchAddress(Form("%c.ladhod.%s.HodoHitTime", spect_prefix, plane_names[i].c_str()), &time_avg[i]);
    T->SetBranchAddress(Form("%c.ladhod.%s.HodoHitPaddleNum", spect_prefix, plane_names[i].c_str()),
                        &time_avg_paddle[i]);
    T->SetBranchAddress(Form("Ndata.%c.ladhod.%s.HodoHitTime", spect_prefix, plane_names[i].c_str()), &nData_hodo[i]);
    T->SetBranchAddress(Form("%c.ladhod.%s.HodoHitPos", spect_prefix, plane_names[i].c_str()), &time_avg_ypos[i]);
  }
  //////////////////////////////////////////////////////////////////
  TH1F *h_gem_cust_dx[2][nPlanes], *h_gem_cust_dy[2][nPlanes];
  TH1F *h_gem_sp_dx[2][nPlanes], *h_gem_sp_dy[2][nPlanes];

  TH1F *h_gem_cust_dx_min[2][nPlanes], *h_gem_cust_dy_min[2][nPlanes];
  TH1F *h_gem_sp_dx_min[2][nPlanes], *h_gem_sp_dy_min[2][nPlanes];
  TH1F *h_gem_cust_dx_min_before[2][nPlanes], *h_gem_cust_dy_min_before[2][nPlanes];
  TH1F *h_gem_sp_dx_min_before[2][nPlanes], *h_gem_sp_dy_min_before[2][nPlanes];
  TH1F *h_gem_cust_dx_min_peak[2][nPlanes], *h_gem_cust_dy_min_peak[2][nPlanes];
  TH1F *h_gem_sp_dx_min_peak[2][nPlanes], *h_gem_sp_dy_min_peak[2][nPlanes];
  TH1F *h_gem_cust_dx_min_after[2][nPlanes], *h_gem_cust_dy_min_after[2][nPlanes];
  TH1F *h_gem_sp_dx_min_after[2][nPlanes], *h_gem_sp_dy_min_after[2][nPlanes];

  TH2D *h_gem_clust_dx_x[2][nPlanes], *h_gem_clust_dy_y[2][nPlanes];
  TH2D *h_gem_sp_dx_x[2][nPlanes], *h_gem_sp_dy_y[2][nPlanes];
  TH2D *h_gem_clust_dx_adc[2][nPlanes], *h_gem_clust_dy_adc[2][nPlanes];
  TH2D *h_gem_sp_dx_adc[2][nPlanes], *h_gem_sp_dy_adc[2][nPlanes];
  TH2D *h_gem_clust_dxy[2][nPlanes], *h_gem_sp_dxy[2][nPlanes];

  TH1F *h_gem_cust_dx_punchthrough[2][nPlanes], *h_gem_cust_dy_punchthrough[2][nPlanes];
  TH1F *h_gem_sp_dx_punchthrough[2][nPlanes], *h_gem_sp_dy_punchthrough[2][nPlanes];

  TH2D *h_gem_clust_dx_x_punchthrough[2][nPlanes], *h_gem_clust_dy_y_punchthrough[2][nPlanes];
  TH2D *h_gem_sp_dx_x_punchthrough[2][nPlanes], *h_gem_sp_dy_y_punchthrough[2][nPlanes];
  TH2D *h_gem_clust_dx_adc_punchthrough[2][nPlanes], *h_gem_clust_dy_adc_punchthrough[2][nPlanes];
  TH2D *h_gem_sp_dx_adc_punchthrough[2][nPlanes], *h_gem_sp_dy_adc_punchthrough[2][nPlanes];
  TH2D *h_gem_clust_dxy_punchthrough[2][nPlanes], *h_gem_sp_dxy_punchthrough[2][nPlanes];

  TH1F *h_gem_cust_dx_punchthrough_min[2][nPlanes], *h_gem_cust_dy_punchthrough_min[2][nPlanes];
  TH1F *h_gem_sp_dx_punchthrough_min[2][nPlanes], *h_gem_sp_dy_punchthrough_min[2][nPlanes];
  TH1F *h_gem_cust_dx_punchthrough_min_before[2][nPlanes], *h_gem_cust_dy_punchthrough_min_before[2][nPlanes];
  TH1F *h_gem_sp_dx_punchthrough_min_before[2][nPlanes], *h_gem_sp_dy_punchthrough_min_before[2][nPlanes];
  TH1F *h_gem_cust_dx_punchthrough_min_peak[2][nPlanes], *h_gem_cust_dy_punchthrough_min_peak[2][nPlanes];
  TH1F *h_gem_sp_dx_punchthrough_min_peak[2][nPlanes], *h_gem_sp_dy_punchthrough_min_peak[2][nPlanes];
  TH1F *h_gem_cust_dx_punchthrough_min_after[2][nPlanes], *h_gem_cust_dy_punchthrough_min_after[2][nPlanes];
  TH1F *h_gem_sp_dx_punchthrough_min_after[2][nPlanes], *h_gem_sp_dy_punchthrough_min_after[2][nPlanes];

  // Hodo Time Histograms
  TH1F *h_hodo_time[nPlanes];
  TH1F *h_hodo_time_GEM0_x[nPlanes];
  TH1F *h_hodo_time_GEM0_y[nPlanes];
  TH1F *h_hodo_time_GEM1_x[nPlanes];
  TH1F *h_hodo_time_GEM1_y[nPlanes];
  TH1F *h_hodo_time_GEM0_x_punchthrough[nPlanes];
  TH1F *h_hodo_time_GEM0_y_punchthrough[nPlanes];
  TH1F *h_hodo_time_GEM1_x_punchthrough[nPlanes];
  TH1F *h_hodo_time_GEM1_y_punchthrough[nPlanes];

  TH1F *h_hodo_time_GEM0_xy[nPlanes];
  TH1F *h_hodo_time_GEM1_xy[nPlanes];
  TH1F *h_hodo_time_GEM_all_x[nPlanes];
  TH1F *h_hodo_time_GEM_all_y[nPlanes];
  TH1F *h_hodo_time_GEM_all_xy[nPlanes];
  TH1F *h_hodo_time_GEM0_xy_punchthrough[nPlanes];
  TH1F *h_hodo_time_GEM1_xy_punchthrough[nPlanes];
  TH1F *h_hodo_time_GEM_all_x_punchthrough[nPlanes];
  TH1F *h_hodo_time_GEM_all_y_punchthrough[nPlanes];
  TH1F *h_hodo_time_GEM_all_xy_punchthrough[nPlanes];

  for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
    for (int i_gem = 0; i_gem < 2; ++i_gem) {
      h_gem_cust_dx[i_gem][i_plane] = new TH1F(
          Form("h_gem_cust_dx_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_cust_dx_%d_%s;Cluster x - Golden Track x (cm);Counts", i_gem, plane_names[i_plane].c_str()),
          dx_NBINS * 2, dx_min * 2, dx_max * 2);
      h_gem_cust_dy[i_gem][i_plane] = new TH1F(
          Form("h_gem_cust_dy_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_cust_dy_%d_%s;Cluster y - Golden Track y (cm);Counts", i_gem, plane_names[i_plane].c_str()),
          dy_NBINS * 2, dy_min * 2, dy_max * 2);
      h_gem_clust_dx_x[i_gem][i_plane] = new TH2D(
          Form("h_gem_clust_dx_x_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_clust_dx_x_%d_%s;Cluster x;Cluster x - Golden Track x (cm)", i_gem, plane_names[i_plane].c_str()),
          dx_NBINS, dx_min, dx_max, dx_NBINS * 2, dx_min * 2, dx_max * 2);
      h_gem_clust_dy_y[i_gem][i_plane] = new TH2D(
          Form("h_gem_clust_dy_y_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_clust_dy_y_%d_%s;Cluster y;Cluster y - Golden Track y (cm)", i_gem, plane_names[i_plane].c_str()),
          dy_NBINS, dy_min, dy_max, dy_NBINS * 2, dy_min * 2, dy_max * 2);

      h_gem_sp_dx[i_gem][i_plane] = new TH1F(
          Form("h_gem_sp_dx_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_sp_dx_%d_%s;Cluster x - Golden Track x (cm);Counts", i_gem, plane_names[i_plane].c_str()),
          dx_NBINS * 2, dx_min * 2, dx_max * 2);
      h_gem_sp_dy[i_gem][i_plane] = new TH1F(
          Form("h_gem_sp_dy_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_sp_dy_%d_%s;Cluster y - Golden Track y (cm);Counts", i_gem, plane_names[i_plane].c_str()),
          dy_NBINS * 2, dy_min * 2, dy_max * 2);
      h_gem_sp_dx_x[i_gem][i_plane] = new TH2D(
          Form("h_gem_sp_dx_x_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_sp_dx_x_%d_%s;Cluster x;Cluster x - Golden Track x (cm)", i_gem, plane_names[i_plane].c_str()),
          dx_NBINS, dx_min, dx_max, dx_NBINS * 2, dx_min * 2, dx_max * 2);
      h_gem_sp_dy_y[i_gem][i_plane] = new TH2D(
          Form("h_gem_sp_dy_y_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_sp_dy_y_%d_%s;Cluster y;Cluster y - Golden Track y (cm)", i_gem, plane_names[i_plane].c_str()),
          dy_NBINS, dy_min, dy_max, dy_NBINS * 2, dy_min * 2, dy_max * 2);

      h_gem_cust_dx_punchthrough[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dx_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dx_punchthrough_%d_%s;Cluster x - Golden Track x (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS * 2, dx_min * 2, dx_max * 2);
      h_gem_cust_dy_punchthrough[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dy_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dy_punchthrough_%d_%s;Cluster y - Golden Track y (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS * 2, dy_min * 2, dy_max * 2);
      h_gem_clust_dx_x_punchthrough[i_gem][i_plane] =
          new TH2D(Form("h_gem_clust_dx_x_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_clust_dx_x_punchthrough_%d_%s;Cluster x;Cluster x - Golden Track x (cm)", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max, dx_NBINS * 2, dx_min * 2, dx_max * 2);
      h_gem_clust_dy_y_punchthrough[i_gem][i_plane] =
          new TH2D(Form("h_gem_clust_dy_y_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_clust_dy_y_punchthrough_%d_%s;Cluster y;Cluster y - Golden Track y (cm)", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max, dy_NBINS * 2, dy_min * 2, dy_max * 2);

      h_gem_sp_dx_punchthrough[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dx_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dx_punchthrough_%d_%s;Cluster x - Golden Track x (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS * 2, dx_min * 2, dx_max * 2);
      h_gem_sp_dy_punchthrough[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dy_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dy_punchthrough_%d_%s;Cluster y - Golden Track y (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS * 2, dy_min * 2, dy_max * 2);
      h_gem_sp_dx_x_punchthrough[i_gem][i_plane] =
          new TH2D(Form("h_gem_sp_dx_x_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dx_x_punchthrough_%d_%s;Cluster x;Cluster x - Golden Track x (cm)", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max, dx_NBINS * 2, dx_min * 2, dx_max * 2);
      h_gem_sp_dy_y_punchthrough[i_gem][i_plane] =
          new TH2D(Form("h_gem_sp_dy_y_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dy_y_punchthrough_%d_%s;Cluster y;Cluster y - Golden Track y (cm)", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max, dy_NBINS * 2, dy_min * 2, dy_max * 2);

      h_gem_clust_dxy[i_gem][i_plane] =
          new TH2D(Form("h_gem_clust_dxy_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_clust_dxy_%d_%s;Cluster x - Golden Track x (cm);Cluster y - Golden Track y (cm)", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);
      h_gem_sp_dxy[i_gem][i_plane] =
          new TH2D(Form("h_gem_sp_dxy_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dxy_%d_%s;Cluster x - Golden Track x (cm);Cluster y - Golden Track y (cm)", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);
      h_gem_clust_dxy_punchthrough[i_gem][i_plane] = new TH2D(
          Form("h_gem_clust_dxy_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_clust_dxy_punchthrough_%d_%s;Cluster x - Golden Track x (cm);Cluster y - Golden Track y (cm)",
               i_gem, plane_names[i_plane].c_str()),
          dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);
      h_gem_sp_dxy_punchthrough[i_gem][i_plane] = new TH2D(
          Form("h_gem_sp_dxy_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_sp_dxy_punchthrough_%d_%s;Cluster x - Golden Track x (cm);Cluster y - Golden Track y (cm)", i_gem,
               plane_names[i_plane].c_str()),
          dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);

      h_gem_cust_dx_min[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dx_min_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dx_min_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_cust_dy_min[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dy_min_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dy_min_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      h_gem_sp_dx_min[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dx_min_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dx_min_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_sp_dy_min[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dy_min_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dy_min_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);

      h_gem_cust_dx_min_before[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dx_min_before_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dx_min_before_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_cust_dy_min_before[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dy_min_before_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dy_min_before_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      h_gem_sp_dx_min_before[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dx_min_before_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dx_min_before_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_sp_dy_min_before[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dy_min_before_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dy_min_before_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      h_gem_cust_dx_min_peak[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dx_min_peak_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dx_min_peak_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_cust_dy_min_peak[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dy_min_peak_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dy_min_peak_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      h_gem_sp_dx_min_peak[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dx_min_peak_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dx_min_peak_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_sp_dy_min_peak[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dy_min_peak_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dy_min_peak_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      h_gem_cust_dx_min_after[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dx_min_after_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dx_min_after_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_cust_dy_min_after[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dy_min_after_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dy_min_after_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      h_gem_sp_dx_min_after[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dx_min_after_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dx_min_after_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_sp_dy_min_after[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dy_min_after_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dy_min_after_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);

      h_gem_cust_dx_punchthrough_min[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dx_punchthrough_min_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dx_punchthrough_min_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_cust_dy_punchthrough_min[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dy_punchthrough_min_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dy_punchthrough_min_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      h_gem_sp_dx_punchthrough_min[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dx_punchthrough_min_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dx_punchthrough_min_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_sp_dy_punchthrough_min[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dy_punchthrough_min_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dy_punchthrough_min_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      h_gem_cust_dx_punchthrough_min_before[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dx_punchthrough_min_before_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dx_punchthrough_min_before_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts",
                        i_gem, plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_cust_dy_punchthrough_min_before[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dy_punchthrough_min_before_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dy_punchthrough_min_before_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts",
                        i_gem, plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      h_gem_sp_dx_punchthrough_min_before[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dx_punchthrough_min_before_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dx_punchthrough_min_before_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_sp_dy_punchthrough_min_before[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dy_punchthrough_min_before_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dy_punchthrough_min_before_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      h_gem_cust_dx_punchthrough_min_peak[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dx_punchthrough_min_peak_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dx_punchthrough_min_peak_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_cust_dy_punchthrough_min_peak[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dy_punchthrough_min_peak_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dy_punchthrough_min_peak_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      h_gem_sp_dx_punchthrough_min_peak[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dx_punchthrough_min_peak_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dx_punchthrough_min_peak_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_sp_dy_punchthrough_min_peak[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dy_punchthrough_min_peak_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dy_punchthrough_min_peak_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      h_gem_cust_dx_punchthrough_min_after[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dx_punchthrough_min_after_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dx_punchthrough_min_after_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts",
                        i_gem, plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_cust_dy_punchthrough_min_after[i_gem][i_plane] =
          new TH1F(Form("h_gem_cust_dy_punchthrough_min_after_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_cust_dy_punchthrough_min_after_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts",
                        i_gem, plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      h_gem_sp_dx_punchthrough_min_after[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dx_punchthrough_min_after_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dx_punchthrough_min_after_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max);
      h_gem_sp_dy_punchthrough_min_after[i_gem][i_plane] =
          new TH1F(Form("h_gem_sp_dy_punchthrough_min_after_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dy_punchthrough_min_after_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max);
      // ADC histograms
      h_gem_clust_dx_adc[i_gem][i_plane] = new TH2D(
          Form("h_gem_clust_dx_adc_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_clust_dx_adc_%d_%s;Cluster x - Golden Track x (cm);ADC", i_gem, plane_names[i_plane].c_str()),
          dx_NBINS, dx_min, dx_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
      h_gem_clust_dy_adc[i_gem][i_plane] = new TH2D(
          Form("h_gem_clust_dy_adc_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_clust_dy_adc_%d_%s;Cluster y - Golden Track y (cm);ADC", i_gem, plane_names[i_plane].c_str()),
          dy_NBINS, dy_min, dy_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
      h_gem_sp_dx_adc[i_gem][i_plane] = new TH2D(
          Form("h_gem_sp_dx_adc_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_sp_dx_adc_%d_%s;Cluster x - Golden Track x (cm);ADC", i_gem, plane_names[i_plane].c_str()),
          dx_NBINS, dx_min, dx_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
      h_gem_sp_dy_adc[i_gem][i_plane] = new TH2D(
          Form("h_gem_sp_dy_adc_%d_%s", i_gem, plane_names[i_plane].c_str()),
          Form("h_gem_sp_dy_adc_%d_%s;Cluster y - Golden Track y (cm);ADC", i_gem, plane_names[i_plane].c_str()),
          dy_NBINS, dy_min, dy_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
      h_gem_clust_dx_adc_punchthrough[i_gem][i_plane] =
          new TH2D(Form("h_gem_clust_dx_adc_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_clust_dx_adc_punchthrough_%d_%s;Cluster x - Golden Track x (cm);ADC", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
      h_gem_clust_dy_adc_punchthrough[i_gem][i_plane] =
          new TH2D(Form("h_gem_clust_dy_adc_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_clust_dy_adc_punchthrough_%d_%s;Cluster y - Golden Track y (cm);ADC", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
      h_gem_sp_dx_adc_punchthrough[i_gem][i_plane] =
          new TH2D(Form("h_gem_sp_dx_adc_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dx_adc_punchthrough_%d_%s;Cluster x - Golden Track x (cm);ADC", i_gem,
                        plane_names[i_plane].c_str()),
                   dx_NBINS, dx_min, dx_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
      h_gem_sp_dy_adc_punchthrough[i_gem][i_plane] =
          new TH2D(Form("h_gem_sp_dy_adc_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()),
                   Form("h_gem_sp_dy_adc_punchthrough_%d_%s;Cluster y - Golden Track y (cm);ADC", i_gem,
                        plane_names[i_plane].c_str()),
                   dy_NBINS, dy_min, dy_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
    }
    // Hodo time histograms

    h_hodo_time[i_plane] = new TH1F(Form("h_hodo_time_%s", plane_names[i_plane].c_str()),
                                    Form("h_hodo_time_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                                    lad_time.nbins, lad_time.min, lad_time.max);
    h_hodo_time_GEM0_x[i_plane] =
        new TH1F(Form("h_hodo_time_GEM0_x_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM0_x_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()), lad_time.nbins,
                 lad_time.min, lad_time.max);
    h_hodo_time_GEM0_y[i_plane] =
        new TH1F(Form("h_hodo_time_GEM0_y_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM0_y_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()), lad_time.nbins,
                 lad_time.min, lad_time.max);
    h_hodo_time_GEM1_x[i_plane] =
        new TH1F(Form("h_hodo_time_GEM1_x_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM1_x_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()), lad_time.nbins,
                 lad_time.min, lad_time.max);
    h_hodo_time_GEM1_y[i_plane] =
        new TH1F(Form("h_hodo_time_GEM1_y_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM1_y_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()), lad_time.nbins,
                 lad_time.min, lad_time.max);
    h_hodo_time_GEM0_x_punchthrough[i_plane] =
        new TH1F(Form("h_hodo_time_GEM0_x_punchthrough_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM0_x_punchthrough_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                 lad_time.nbins, lad_time.min, lad_time.max);
    h_hodo_time_GEM0_y_punchthrough[i_plane] =
        new TH1F(Form("h_hodo_time_GEM0_y_punchthrough_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM0_y_punchthrough_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                 lad_time.nbins, lad_time.min, lad_time.max);
    h_hodo_time_GEM1_x_punchthrough[i_plane] =
        new TH1F(Form("h_hodo_time_GEM1_x_punchthrough_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM1_x_punchthrough_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                 lad_time.nbins, lad_time.min, lad_time.max);
    h_hodo_time_GEM1_y_punchthrough[i_plane] =
        new TH1F(Form("h_hodo_time_GEM1_y_punchthrough_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM1_y_punchthrough_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                 lad_time.nbins, lad_time.min, lad_time.max);

    h_hodo_time_GEM0_xy[i_plane] =
        new TH1F(Form("h_hodo_time_GEM0_xy_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM0_xy_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                 lad_time.nbins, lad_time.min, lad_time.max);
    h_hodo_time_GEM1_xy[i_plane] =
        new TH1F(Form("h_hodo_time_GEM1_xy_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM1_xy_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                 lad_time.nbins, lad_time.min, lad_time.max);
    h_hodo_time_GEM_all_x[i_plane] =
        new TH1F(Form("h_hodo_time_GEM_all_x_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM_all_x_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                 lad_time.nbins, lad_time.min, lad_time.max);
    h_hodo_time_GEM_all_y[i_plane] =
        new TH1F(Form("h_hodo_time_GEM_all_y_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM_all_y_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                 lad_time.nbins, lad_time.min, lad_time.max);
    h_hodo_time_GEM_all_xy[i_plane] =
        new TH1F(Form("h_hodo_time_GEM_all_xy_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM_all_xy_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                 lad_time.nbins, lad_time.min, lad_time.max);

    h_hodo_time_GEM0_xy_punchthrough[i_plane] =
        new TH1F(Form("h_hodo_time_GEM0_xy_punchthrough_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM0_xy_punchthrough_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                 lad_time.nbins, lad_time.min, lad_time.max);
    h_hodo_time_GEM1_xy_punchthrough[i_plane] =
        new TH1F(Form("h_hodo_time_GEM1_xy_punchthrough_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM1_xy_punchthrough_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                 lad_time.nbins, lad_time.min, lad_time.max);
    h_hodo_time_GEM_all_x_punchthrough[i_plane] =
        new TH1F(Form("h_hodo_time_GEM_all_x_punchthrough_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM_all_x_punchthrough_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                 lad_time.nbins, lad_time.min, lad_time.max);
    h_hodo_time_GEM_all_y_punchthrough[i_plane] =
        new TH1F(Form("h_hodo_time_GEM_all_y_punchthrough_%s", plane_names[i_plane].c_str()),
                 Form("h_hodo_time_GEM_all_y_punchthrough_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
                 lad_time.nbins, lad_time.min, lad_time.max);
    h_hodo_time_GEM_all_xy_punchthrough[i_plane] = new TH1F(
        Form("h_hodo_time_GEM_all_xy_punchthrough_%s", plane_names[i_plane].c_str()),
        Form("h_hodo_time_GEM_all_xy_punchthrough_%s;Hodoscope time (ns);Counts", plane_names[i_plane].c_str()),
        lad_time.nbins, lad_time.min, lad_time.max);
  }
  ///////////////////////////////////////////////////////////////////

  TVector3 p1[2], p2[2], p3[2];
  for (int i = 0; i < 2; ++i) {
    double theta = gem_theta;
    p1[i]        = TVector3(gem_r[i] * cos((theta - 90) * TMath::DegToRad()), 0,
                            -gem_r[i] * sin((theta - 90) * TMath::DegToRad()));
    p2[i]        = p1[i] + TVector3(-gem_r[i] * cos((180 - theta) * TMath::DegToRad()), 0,
                                    -gem_r[i] * sin((180 - theta) * TMath::DegToRad()));
    p3[i]        = p1[i] + TVector3(0, 10.0, 0); // Arbitrary y value of 10.0
  }

  Bool_t has_hodo_hit_avg[nPlanes][nPaddles];
  // Loop through the tree entries
  Long64_t nEntries = T->GetEntries();
  // nEntries          = 100000; // For testing purposes, limit to 1000 entries
  for (Long64_t i = 0; i < nEntries; ++i) {
    T->GetEntry(i);
    if (abs(vertex_x) > 50 || abs(vertex_y) > 10 || abs(vertex_z) > 10) {
      continue;
    }
    /////////////////////////////////////////////////////////
    // Check for matching hits in the hodoscope
    for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
      for (int i_paddle = 0; i_paddle < nPaddles; ++i_paddle) {
        has_hodo_hit_avg[i_plane][i_paddle] = false;
      }
    }

    for (int plane = 0; plane < nPlanes; ++plane) {
      for (int i_hit = 0; i_hit < nData_hodo[plane]; ++i_hit) {
        int bar = int(time_avg_paddle[plane][i_hit]) - 1;
        if (bar < 0 || bar >= nPaddles)
          continue; // Skip invalid bars
        has_hodo_hit_avg[plane][bar] = true;
      }
    }
    ////////////////////////////////////////////////////////
    for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
      TVector3 vertex(vertex_x, vertex_y, vertex_z);
      int matching_plane;
      switch (i_plane) {
      case 0:
        matching_plane = 1;
        break;
      case 1:
        matching_plane = 0;
        break;
      case 2:
        matching_plane = 3;
        break;
      case 3:
        matching_plane = 2;
        break;
      case 4:
        matching_plane = 4;
        break;
      }
      for (int i_hit = 0; i_hit < nData_hodo[i_plane]; ++i_hit) {
        int paddle        = static_cast<int>(time_avg_paddle[i_plane][i_hit]) - 1;
        TVector3 hodo_pos = GetHodoHitPosition(paddle, i_plane);
        hodo_pos.SetY(time_avg_ypos[i_plane][i_hit]);

        bool is_punchthrough = has_hodo_hit_avg[matching_plane][paddle];
        int time_id;
        if (time_avg[i_plane][i_hit] >= time_peak.min && time_avg[i_plane][i_hit] < time_peak.max) {
          time_id = 1; // peak
        } else if (time_avg[i_plane][i_hit] >= time_before.min && time_avg[i_plane][i_hit] < time_before.max) {
          time_id = 0; // before
        } else if (time_avg[i_plane][i_hit] >= time_after.min && time_avg[i_plane][i_hit] < time_after.max) {
          time_id = 2; // after
        } else {
          time_id = -1; // out of range
        }
        // Calculate the intersection point

        // The line is from the vertex to the hodo hit position
        TVector3 intersection[2] = {LinePlaneIntersection(p1[0], p2[0], p3[0], vertex, hodo_pos),
                                    LinePlaneIntersection(p1[1], p2[1], p3[1], vertex, hodo_pos)};

        // intersection now contains the intersection point of the track (vertex->hodo) with the hodoscope plane

        double clust_min_dx[2]              = {9999, 9999};
        double clust_min_dy[2]              = {9999, 9999};
        double clust_min_dx_punchthrough[2] = {9999, 9999};
        double clust_min_dy_punchthrough[2] = {9999, 9999};

        for (int i_clust = 0; i_clust < nData_clust; ++i_clust) {
          double clust_position = -clust_pos[i_clust];
          int clust_layer_id    = clust_layer[i_clust];
          double clust_axis_id  = clust_axis[i_clust];
          double clust_adc_val  = clust_adc[i_clust];
          TVector3 gem_hit;
          if (clust_axis_id == 1) { // x axis
            gem_hit    = GetGEMHitPosition(clust_position, 0, gem_dx[clust_layer_id], gem_dy[clust_layer_id],
                                           gem_r[clust_layer_id], gem_theta);
            double dxz = sqrt(pow(gem_hit.X() - intersection[clust_layer_id].X(), 2) +
                              pow(gem_hit.Z() - intersection[clust_layer_id].Z(), 2));
            if (gem_hit.X() > intersection[clust_layer_id].X()) {
              dxz = -dxz;
            }
            h_gem_cust_dx[clust_layer_id][i_plane]->Fill(dxz);
            h_gem_clust_dx_x[clust_layer_id][i_plane]->Fill(clust_position, dxz);
            h_gem_clust_dx_adc[clust_layer_id][i_plane]->Fill(dxz, clust_adc_val);
            if (abs(clust_min_dx[clust_layer_id]) > abs(dxz)) {
              clust_min_dx[clust_layer_id] = dxz;
            }
            if (is_punchthrough) {
              h_gem_cust_dx_punchthrough[clust_layer_id][i_plane]->Fill(dxz);
              h_gem_clust_dx_x_punchthrough[clust_layer_id][i_plane]->Fill(clust_position, dxz);
              h_gem_clust_dx_adc_punchthrough[clust_layer_id][i_plane]->Fill(dxz, clust_adc_val);
              if (abs(clust_min_dx_punchthrough[clust_layer_id]) > abs(dxz)) {
                clust_min_dx_punchthrough[clust_layer_id] = dxz;
              }
            }
          } else if (clust_axis_id == 0) { // y axis
            gem_hit   = GetGEMHitPosition(0, clust_position, gem_dx[clust_layer_id], gem_dy[clust_layer_id],
                                          gem_r[clust_layer_id], gem_theta);
            double dy = gem_hit.Y() - intersection[clust_layer_id].Y();
            h_gem_cust_dy[clust_layer_id][i_plane]->Fill(dy);
            h_gem_clust_dy_y[clust_layer_id][i_plane]->Fill(clust_position, dy);
            h_gem_clust_dy_adc[clust_layer_id][i_plane]->Fill(clust_adc_val, dy);
            if (abs(clust_min_dy[clust_layer_id]) > abs(dy)) {
              clust_min_dy[clust_layer_id] = dy;
            }
            if (is_punchthrough) {
              h_gem_cust_dy_punchthrough[clust_layer_id][i_plane]->Fill(dy);
              h_gem_clust_dy_y_punchthrough[clust_layer_id][i_plane]->Fill(clust_position, dy);
              h_gem_clust_dy_adc_punchthrough[clust_layer_id][i_plane]->Fill(clust_adc_val, dy);
              if (abs(clust_min_dy_punchthrough[clust_layer_id]) > abs(dy)) {
                clust_min_dy_punchthrough[clust_layer_id] = dy;
              }
            }
          } else {
            std::cerr << "Error: Invalid cluster axis ID!" << std::endl;
            continue;
          }
        }
        // Fill the minimum dx and dy histograms
        h_gem_cust_dx_min[0][i_plane]->Fill(clust_min_dx[0]);
        h_gem_cust_dy_min[0][i_plane]->Fill(clust_min_dy[0]);
        h_gem_cust_dx_min[1][i_plane]->Fill(clust_min_dx[1]);
        h_gem_cust_dy_min[1][i_plane]->Fill(clust_min_dy[1]);
        if (time_id == 0) { // before
          h_gem_cust_dx_min_before[0][i_plane]->Fill(clust_min_dx[0]);
          h_gem_cust_dy_min_before[0][i_plane]->Fill(clust_min_dy[0]);
          h_gem_cust_dx_min_before[1][i_plane]->Fill(clust_min_dx[1]);
          h_gem_cust_dy_min_before[1][i_plane]->Fill(clust_min_dy[1]);
        } else if (time_id == 1) { // peak
          h_gem_cust_dx_min_peak[0][i_plane]->Fill(clust_min_dx[0]);
          h_gem_cust_dy_min_peak[0][i_plane]->Fill(clust_min_dy[0]);
          h_gem_cust_dx_min_peak[1][i_plane]->Fill(clust_min_dx[1]);
          h_gem_cust_dy_min_peak[1][i_plane]->Fill(clust_min_dy[1]);
        } else if (time_id == 2) { // after
          h_gem_cust_dx_min_after[0][i_plane]->Fill(clust_min_dx[0]);
          h_gem_cust_dy_min_after[0][i_plane]->Fill(clust_min_dy[0]);
          h_gem_cust_dx_min_after[1][i_plane]->Fill(clust_min_dx[1]);
          h_gem_cust_dy_min_after[1][i_plane]->Fill(clust_min_dy[1]);
        }

        // Fill the punchthrough minimum dx and dy histograms
        h_gem_cust_dx_punchthrough_min[0][i_plane]->Fill(clust_min_dx_punchthrough[0]);
        h_gem_cust_dy_punchthrough_min[0][i_plane]->Fill(clust_min_dy_punchthrough[0]);
        h_gem_cust_dx_punchthrough_min[1][i_plane]->Fill(clust_min_dx_punchthrough[1]);
        h_gem_cust_dy_punchthrough_min[1][i_plane]->Fill(clust_min_dy_punchthrough[1]);
        if (time_id == 0) { // before
          h_gem_cust_dx_punchthrough_min_before[0][i_plane]->Fill(clust_min_dx_punchthrough[0]);
          h_gem_cust_dy_punchthrough_min_before[0][i_plane]->Fill(clust_min_dy_punchthrough[0]);
          h_gem_cust_dx_punchthrough_min_before[1][i_plane]->Fill(clust_min_dx_punchthrough[1]);
          h_gem_cust_dy_punchthrough_min_before[1][i_plane]->Fill(clust_min_dy_punchthrough[1]);
        } else if (time_id == 1) { // peak
          h_gem_cust_dx_punchthrough_min_peak[0][i_plane]->Fill(clust_min_dx_punchthrough[0]);
          h_gem_cust_dy_punchthrough_min_peak[0][i_plane]->Fill(clust_min_dy_punchthrough[0]);
          h_gem_cust_dx_punchthrough_min_peak[1][i_plane]->Fill(clust_min_dx_punchthrough[1]);
          h_gem_cust_dy_punchthrough_min_peak[1][i_plane]->Fill(clust_min_dy_punchthrough[1]);
        } else if (time_id == 2) { // after
          h_gem_cust_dx_punchthrough_min_after[0][i_plane]->Fill(clust_min_dx_punchthrough[0]);
          h_gem_cust_dy_punchthrough_min_after[0][i_plane]->Fill(clust_min_dy_punchthrough[0]);
          h_gem_cust_dx_punchthrough_min_after[1][i_plane]->Fill(clust_min_dx_punchthrough[1]);
          h_gem_cust_dy_punchthrough_min_after[1][i_plane]->Fill(clust_min_dy_punchthrough[1]);
        }

        // Fill time histograms
        h_hodo_time[i_plane]->Fill(time_avg[i_plane][i_hit]);
        if (fabs(clust_min_dx[0]) < MAX_DX_CLUST && fabs(clust_min_dy[0]) < MAX_DY_CLUST) {
          h_hodo_time_GEM0_xy[i_plane]->Fill(time_avg[i_plane][i_hit]);
          if (is_punchthrough)
            h_hodo_time_GEM0_xy_punchthrough[i_plane]->Fill(time_avg[i_plane][i_hit]);
        }
        if (fabs(clust_min_dx[1]) < MAX_DX_CLUST && fabs(clust_min_dy[1]) < MAX_DY_CLUST) {
          h_hodo_time_GEM1_xy[i_plane]->Fill(time_avg[i_plane][i_hit]);
          if (is_punchthrough)
            h_hodo_time_GEM1_xy_punchthrough[i_plane]->Fill(time_avg[i_plane][i_hit]);
        }
        if (fabs(clust_min_dx[0]) < MAX_DX_CLUST) {
          h_hodo_time_GEM0_x[i_plane]->Fill(time_avg[i_plane][i_hit]);
          if (is_punchthrough)
            h_hodo_time_GEM0_x_punchthrough[i_plane]->Fill(time_avg[i_plane][i_hit]);
        }
        if (fabs(clust_min_dy[0]) < MAX_DY_CLUST) {
          h_hodo_time_GEM0_y[i_plane]->Fill(time_avg[i_plane][i_hit]);
          if (is_punchthrough)
            h_hodo_time_GEM0_y_punchthrough[i_plane]->Fill(time_avg[i_plane][i_hit]);
        }
        if (fabs(clust_min_dx[1]) < MAX_DX_CLUST) {
          h_hodo_time_GEM1_x[i_plane]->Fill(time_avg[i_plane][i_hit]);
          if (is_punchthrough)
            h_hodo_time_GEM1_x_punchthrough[i_plane]->Fill(time_avg[i_plane][i_hit]);
        }
        if (fabs(clust_min_dy[1]) < MAX_DY_CLUST) {
          h_hodo_time_GEM1_y[i_plane]->Fill(time_avg[i_plane][i_hit]);
          if (is_punchthrough)
            h_hodo_time_GEM1_y_punchthrough[i_plane]->Fill(time_avg[i_plane][i_hit]);
        }
        if (fabs(clust_min_dx[0]) < MAX_DX_CLUST && fabs(clust_min_dx[1]) < MAX_DX_CLUST) {
          h_hodo_time_GEM_all_x[i_plane]->Fill(time_avg[i_plane][i_hit]);
          if (is_punchthrough)
            h_hodo_time_GEM_all_x_punchthrough[i_plane]->Fill(time_avg[i_plane][i_hit]);
        }
        if (fabs(clust_min_dy[0]) < MAX_DY_CLUST && fabs(clust_min_dy[1]) < MAX_DY_CLUST) {
          h_hodo_time_GEM_all_y[i_plane]->Fill(time_avg[i_plane][i_hit]);
          if (is_punchthrough)
            h_hodo_time_GEM_all_y_punchthrough[i_plane]->Fill(time_avg[i_plane][i_hit]);
        }
        if (fabs(clust_min_dx[0]) < MAX_DX_CLUST && fabs(clust_min_dy[0]) < MAX_DY_CLUST &&
            fabs(clust_min_dx[1]) < MAX_DX_CLUST && fabs(clust_min_dy[1]) < MAX_DY_CLUST) {
          h_hodo_time_GEM_all_xy[i_plane]->Fill(time_avg[i_plane][i_hit]);
          if (is_punchthrough)
            h_hodo_time_GEM_all_xy_punchthrough[i_plane]->Fill(time_avg[i_plane][i_hit]);
        }

        // Loop through the space points
        for (int i_sp = 0; i_sp < nData_sp; ++i_sp) {
          TVector3 gem_hit;
          int sp_layer_id   = static_cast<int>(sp_layer[i_sp]);
          double sp_adc_val = sp_adc[i_sp];
          gem_hit           = GetGEMHitPosition(sp_X[i_sp], sp_X[i_sp], gem_dx[sp_layer_id], gem_dy[sp_layer_id],
                                                gem_r[sp_layer_id], gem_theta);
          double dxz        = sqrt(pow(gem_hit.X() - intersection[sp_layer_id].X(), 2) +
                                   pow(gem_hit.Z() - intersection[sp_layer_id].Z(), 2));
          if (gem_hit.X() > intersection[sp_layer_id].X()) {
            dxz = -dxz;
          }
          double dy = gem_hit.Y() - intersection[sp_layer_id].Y();

          h_gem_sp_dx[sp_layer_id][i_plane]->Fill(dxz);
          h_gem_sp_dy[sp_layer_id][i_plane]->Fill(dy);
          h_gem_sp_dx_x[sp_layer_id][i_plane]->Fill(sp_X[i_sp], dxz);
          h_gem_sp_dy_y[sp_layer_id][i_plane]->Fill(sp_Y[i_sp], dy);
          h_gem_sp_dx_adc[sp_layer_id][i_plane]->Fill(dxz, sp_adc_val);
          h_gem_sp_dy_adc[sp_layer_id][i_plane]->Fill(dy, sp_adc_val);
          h_gem_sp_dxy[sp_layer_id][i_plane]->Fill(dxz, dy);
          if (is_punchthrough) {
            h_gem_sp_dx_punchthrough[sp_layer_id][i_plane]->Fill(dxz);
            h_gem_sp_dy_punchthrough[sp_layer_id][i_plane]->Fill(dy);
            h_gem_sp_dx_x_punchthrough[sp_layer_id][i_plane]->Fill(sp_X[i_sp], dxz);
            h_gem_sp_dy_y_punchthrough[sp_layer_id][i_plane]->Fill(sp_Y[i_sp], dy);
            h_gem_sp_dx_adc_punchthrough[sp_layer_id][i_plane]->Fill(dxz, sp_adc_val);
            h_gem_sp_dy_adc_punchthrough[sp_layer_id][i_plane]->Fill(dy, sp_adc_val);
            h_gem_sp_dxy_punchthrough[sp_layer_id][i_plane]->Fill(dxz, dy);
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
  outputFile->mkdir("clust/x/dx");
  outputFile->mkdir("clust/x/2D");
  outputFile->mkdir("clust/x/2D_ADC");
  outputFile->mkdir("clust/y/dy");
  outputFile->mkdir("clust/y/2D");
  outputFile->mkdir("clust/y/2D_ADC");
  outputFile->mkdir("clust/y/dy_punchthrough");
  outputFile->mkdir("clust/y/2D_punchthrough");
  outputFile->mkdir("clust/y/2D_ADC_punchthrough");
  outputFile->mkdir("clust/x/dx_punchthrough");
  outputFile->mkdir("clust/x/2D_punchthrough");
  outputFile->mkdir("clust/x/2D_ADC_punchthrough");
  outputFile->mkdir("clust/dxdy/all");
  outputFile->mkdir("clust/dxdy/punchthrough");
  outputFile->mkdir("clust/x/dx/min/all");
  outputFile->mkdir("clust/x/dx/min/before");
  outputFile->mkdir("clust/x/dx/min/peak");
  outputFile->mkdir("clust/x/dx/min/after");
  outputFile->mkdir("clust/x/dx/min/comp");
  outputFile->mkdir("clust/y/dy/min/all");
  outputFile->mkdir("clust/y/dy/min/before");
  outputFile->mkdir("clust/y/dy/min/peak");
  outputFile->mkdir("clust/y/dy/min/after");
  outputFile->mkdir("clust/y/dy/min/comp");
  outputFile->mkdir("clust/x/dx/min/before_punchthrough");
  outputFile->mkdir("clust/y/dy/min/before_punchthrough");
  outputFile->mkdir("clust/x/dx/min/peak_punchthrough");
  outputFile->mkdir("clust/y/dy/min/peak_punchthrough");
  outputFile->mkdir("clust/x/dx/min/after_punchthrough");
  outputFile->mkdir("clust/y/dy/min/after_punchthrough");
  outputFile->mkdir("clust/x/dx/min/comp_punchthrough");
  outputFile->mkdir("clust/y/dy/min/comp_punchthrough");

  for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
    outputFile->mkdir(Form("time/%s", plane_names[i_plane].c_str()));
  }

  outputFile->mkdir("sp/x/dx");
  outputFile->mkdir("sp/x/2D");
  outputFile->mkdir("sp/x/2D_ADC");
  outputFile->mkdir("sp/y/dy");
  outputFile->mkdir("sp/y/2D");
  outputFile->mkdir("sp/y/2D_ADC");
  outputFile->mkdir("sp/x/dx_punchthrough");
  outputFile->mkdir("sp/x/2D_punchthrough");
  outputFile->mkdir("sp/x/2D_ADC_punchthrough");
  outputFile->mkdir("sp/y/dy_punchthrough");
  outputFile->mkdir("sp/y/2D_punchthrough");
  outputFile->mkdir("sp/y/2D_ADC_punchthrough");
  outputFile->mkdir("sp/dxdy/all");
  outputFile->mkdir("sp/dxdy/punchthrough");
  // Write histograms to the output file
  for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
    for (int i_gem = 0; i_gem < 2; ++i_gem) {
      outputFile->cd("clust/x/dx");
      h_gem_cust_dx[i_gem][i_plane]->Write();
      outputFile->cd("clust/x/2D");
      h_gem_clust_dx_x[i_gem][i_plane]->Write();
      outputFile->cd("clust/x/2D_ADC");
      h_gem_clust_dx_adc[i_gem][i_plane]->SetOption("COLZ");
      h_gem_clust_dx_adc[i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy");
      h_gem_cust_dy[i_gem][i_plane]->Write();
      outputFile->cd("clust/y/2D");
      h_gem_clust_dy_y[i_gem][i_plane]->Write();
      outputFile->cd("clust/y/2D_ADC");
      h_gem_clust_dy_adc[i_gem][i_plane]->SetOption("COLZ");
      h_gem_clust_dy_adc[i_gem][i_plane]->Write();
      outputFile->cd("clust/x/dx_punchthrough");
      h_gem_cust_dx_punchthrough[i_gem][i_plane]->Write();
      outputFile->cd("clust/x/2D_punchthrough");
      h_gem_clust_dx_x_punchthrough[i_gem][i_plane]->Write();
      outputFile->cd("clust/x/2D_ADC_punchthrough");
      h_gem_clust_dx_adc_punchthrough[i_gem][i_plane]->SetOption("COLZ");
      h_gem_clust_dx_adc_punchthrough[i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy_punchthrough");
      h_gem_cust_dy_punchthrough[i_gem][i_plane]->Write();
      outputFile->cd("clust/y/2D_punchthrough");
      h_gem_clust_dy_y_punchthrough[i_gem][i_plane]->Write();
      outputFile->cd("clust/y/2D_ADC_punchthrough");
      h_gem_clust_dy_adc_punchthrough[i_gem][i_plane]->SetOption("COLZ");
      h_gem_clust_dy_adc_punchthrough[i_gem][i_plane]->Write();
      outputFile->cd("clust/dxdy/all");
      h_gem_clust_dxy[i_gem][i_plane]->Write();
      outputFile->cd("clust/dxdy/punchthrough");
      h_gem_clust_dxy_punchthrough[i_gem][i_plane]->Write();

      outputFile->cd("clust/x/dx/min/all");
      h_gem_cust_dx_min[i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/all");
      h_gem_cust_dy_min[i_gem][i_plane]->Write();
      outputFile->cd("clust/x/dx/min/before");
      h_gem_cust_dx_min_before[i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/before");
      h_gem_cust_dy_min_before[i_gem][i_plane]->Write();
      outputFile->cd("clust/x/dx/min/peak");
      h_gem_cust_dx_min_peak[i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/peak");
      h_gem_cust_dy_min_peak[i_gem][i_plane]->Write();
      outputFile->cd("clust/x/dx/min/after");
      h_gem_cust_dx_min_after[i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/after");
      h_gem_cust_dy_min_after[i_gem][i_plane]->Write();

      outputFile->cd("clust/x/dx/min/before_punchthrough");
      h_gem_cust_dx_punchthrough_min_before[i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/before_punchthrough");
      h_gem_cust_dy_punchthrough_min_before[i_gem][i_plane]->Write();
      outputFile->cd("clust/x/dx/min/peak_punchthrough");
      h_gem_cust_dx_punchthrough_min_peak[i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/peak_punchthrough");
      h_gem_cust_dy_punchthrough_min_peak[i_gem][i_plane]->Write();
      outputFile->cd("clust/x/dx/min/after_punchthrough");
      h_gem_cust_dx_punchthrough_min_after[i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/after_punchthrough");
      h_gem_cust_dy_punchthrough_min_after[i_gem][i_plane]->Write();
      // outputFile->cd("clust/x/dx/min/comp_punchthrough");
      // h_gem_cust_dx_punchthrough_min[i_gem][i_plane]->Write();
      // outputFile->cd("clust/y/dy/min/comp_punchthrough");
      // h_gem_cust_dy_punchthrough_min[i_gem][i_plane]->Write();
      // Comparison Histograms
      TCanvas *c1;
      TLegend *legend;
      outputFile->cd("clust/x/dx/min/comp");
      c1 = new TCanvas(Form("c1_x_%d_%s", i_gem, plane_names[i_plane].c_str()), "Comparison Canvas", 800, 600);
      c1->cd();
      h_gem_cust_dx_min_before[i_gem][i_plane]->SetLineColor(kBlue);
      h_gem_cust_dx_min_peak[i_gem][i_plane]->SetLineColor(kRed);
      h_gem_cust_dx_min_after[i_gem][i_plane]->SetLineColor(kBlack);

      h_gem_cust_dx_min_before[i_gem][i_plane]->Scale(1.0 / h_gem_cust_dx_min_before[i_gem][i_plane]->Integral());
      h_gem_cust_dx_min_peak[i_gem][i_plane]->Scale(1.0 / h_gem_cust_dx_min_peak[i_gem][i_plane]->Integral());
      h_gem_cust_dx_min_after[i_gem][i_plane]->Scale(1.0 / h_gem_cust_dx_min_after[i_gem][i_plane]->Integral());
      h_gem_cust_dx_min_before[i_gem][i_plane]->Draw("hist");
      h_gem_cust_dx_min_peak[i_gem][i_plane]->Draw("hist same");
      h_gem_cust_dx_min_after[i_gem][i_plane]->Draw("hist same");
      legend = new TLegend(0.65, 0.7, 0.88, 0.88);
      legend->AddEntry(h_gem_cust_dx_min_before[i_gem][i_plane], "Before", "l");
      legend->AddEntry(h_gem_cust_dx_min_peak[i_gem][i_plane], "Peak", "l");
      legend->AddEntry(h_gem_cust_dx_min_after[i_gem][i_plane], "After", "l");
      legend->Draw();
      c1->Update();
      c1->Write();
      outputFile->cd("clust/y/dy/min/comp");
      c1 = new TCanvas(Form("c1_y_%d_%s", i_gem, plane_names[i_plane].c_str()), "Comparison Canvas", 800, 600);
      c1->cd();
      h_gem_cust_dy_min_before[i_gem][i_plane]->SetLineColor(kBlue);
      h_gem_cust_dy_min_peak[i_gem][i_plane]->SetLineColor(kRed);
      h_gem_cust_dy_min_after[i_gem][i_plane]->SetLineColor(kBlack);
      h_gem_cust_dy_min_before[i_gem][i_plane]->Scale(1.0 / h_gem_cust_dy_min_before[i_gem][i_plane]->Integral());
      h_gem_cust_dy_min_peak[i_gem][i_plane]->Scale(1.0 / h_gem_cust_dy_min_peak[i_gem][i_plane]->Integral());
      h_gem_cust_dy_min_after[i_gem][i_plane]->Scale(1.0 / h_gem_cust_dy_min_after[i_gem][i_plane]->Integral());
      h_gem_cust_dy_min_before[i_gem][i_plane]->Draw("hist");
      h_gem_cust_dy_min_peak[i_gem][i_plane]->Draw("hist same");
      h_gem_cust_dy_min_after[i_gem][i_plane]->Draw("hist same");
      legend = new TLegend(0.65, 0.7, 0.88, 0.88);
      legend->AddEntry(h_gem_cust_dy_min_before[i_gem][i_plane], "Before", "l");
      legend->AddEntry(h_gem_cust_dy_min_peak[i_gem][i_plane], "Peak", "l");
      legend->AddEntry(h_gem_cust_dy_min_after[i_gem][i_plane], "After", "l");
      legend->Draw();
      c1->Update();
      c1->Write();

      outputFile->cd("clust/x/dx/min/comp_punchthrough");
      c1 = new TCanvas(Form("c1_x_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()), "Comparison Canvas", 800,
                       600);
      c1->cd();
      h_gem_cust_dx_punchthrough_min_before[i_gem][i_plane]->SetLineColor(kBlue);
      h_gem_cust_dx_punchthrough_min_peak[i_gem][i_plane]->SetLineColor(kRed);
      h_gem_cust_dx_punchthrough_min_after[i_gem][i_plane]->SetLineColor(kBlack);
      h_gem_cust_dx_punchthrough_min_before[i_gem][i_plane]->Scale(
          1.0 / h_gem_cust_dx_punchthrough_min_before[i_gem][i_plane]->Integral());
      h_gem_cust_dx_punchthrough_min_peak[i_gem][i_plane]->Scale(
          1.0 / h_gem_cust_dx_punchthrough_min_peak[i_gem][i_plane]->Integral());
      h_gem_cust_dx_punchthrough_min_after[i_gem][i_plane]->Scale(
          1.0 / h_gem_cust_dx_punchthrough_min_after[i_gem][i_plane]->Integral());
      h_gem_cust_dx_punchthrough_min_before[i_gem][i_plane]->Draw("hist");
      h_gem_cust_dx_punchthrough_min_peak[i_gem][i_plane]->Draw("hist same");
      h_gem_cust_dx_punchthrough_min_after[i_gem][i_plane]->Draw("hist same");
      legend = new TLegend(0.65, 0.7, 0.88, 0.88);
      legend->AddEntry(h_gem_cust_dx_punchthrough_min_before[i_gem][i_plane], "Before", "l");
      legend->AddEntry(h_gem_cust_dx_punchthrough_min_peak[i_gem][i_plane], "Peak", "l");
      legend->AddEntry(h_gem_cust_dx_punchthrough_min_after[i_gem][i_plane], "After", "l");
      legend->Draw();
      c1->Update();
      c1->Write();
      outputFile->cd("clust/y/dy/min/comp_punchthrough");
      c1 = new TCanvas(Form("c1_y_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()), "Comparison Canvas", 800,
                       600);
      c1->cd();
      h_gem_cust_dy_punchthrough_min_before[i_gem][i_plane]->SetLineColor(kBlue);
      h_gem_cust_dy_punchthrough_min_peak[i_gem][i_plane]->SetLineColor(kRed);
      h_gem_cust_dy_punchthrough_min_after[i_gem][i_plane]->SetLineColor(kBlack);
      h_gem_cust_dy_punchthrough_min_before[i_gem][i_plane]->Scale(
          1.0 / h_gem_cust_dy_punchthrough_min_before[i_gem][i_plane]->Integral());
      h_gem_cust_dy_punchthrough_min_peak[i_gem][i_plane]->Scale(
          1.0 / h_gem_cust_dy_punchthrough_min_peak[i_gem][i_plane]->Integral());
      h_gem_cust_dy_punchthrough_min_after[i_gem][i_plane]->Scale(
          1.0 / h_gem_cust_dy_punchthrough_min_after[i_gem][i_plane]->Integral());
      h_gem_cust_dy_punchthrough_min_before[i_gem][i_plane]->Draw("hist");
      h_gem_cust_dy_punchthrough_min_peak[i_gem][i_plane]->Draw("hist same");
      h_gem_cust_dy_punchthrough_min_after[i_gem][i_plane]->Draw("hist same");
      legend = new TLegend(0.65, 0.7, 0.88, 0.88);
      legend->AddEntry(h_gem_cust_dy_punchthrough_min_before[i_gem][i_plane], "Before", "l");
      legend->AddEntry(h_gem_cust_dy_punchthrough_min_peak[i_gem][i_plane], "Peak", "l");
      legend->AddEntry(h_gem_cust_dy_punchthrough_min_after[i_gem][i_plane], "After", "l");
      legend->Draw();
      c1->Update();
      c1->Write();

      // Space point histograms
      outputFile->cd("sp/x/dx");
      h_gem_sp_dx[i_gem][i_plane]->Write();
      outputFile->cd("sp/x/2D");
      h_gem_sp_dx_x[i_gem][i_plane]->Write();
      outputFile->cd("sp/x/2D_ADC");
      h_gem_sp_dx_adc[i_gem][i_plane]->SetOption("COLZ");
      h_gem_sp_dx_adc[i_gem][i_plane]->Write();
      outputFile->cd("sp/y/dy");
      h_gem_sp_dy[i_gem][i_plane]->Write();
      outputFile->cd("sp/y/2D");
      h_gem_sp_dy_y[i_gem][i_plane]->Write();
      outputFile->cd("sp/y/2D_ADC");
      h_gem_sp_dy_adc[i_gem][i_plane]->SetOption("COLZ");
      h_gem_sp_dy_adc[i_gem][i_plane]->Write();
      outputFile->cd("sp/x/dx_punchthrough");
      h_gem_sp_dx_punchthrough[i_gem][i_plane]->Write();
      outputFile->cd("sp/x/2D_punchthrough");
      h_gem_sp_dx_x_punchthrough[i_gem][i_plane]->Write();
      outputFile->cd("sp/x/2D_ADC_punchthrough");
      h_gem_sp_dx_adc_punchthrough[i_gem][i_plane]->SetOption("COLZ");
      h_gem_sp_dx_adc_punchthrough[i_gem][i_plane]->Write();
      outputFile->cd("sp/y/dy_punchthrough");
      h_gem_sp_dy_punchthrough[i_gem][i_plane]->Write();
      outputFile->cd("sp/y/2D_punchthrough");
      h_gem_sp_dy_y_punchthrough[i_gem][i_plane]->Write();
      outputFile->cd("sp/y/2D_ADC_punchthrough");
      h_gem_sp_dy_adc_punchthrough[i_gem][i_plane]->SetOption("COLZ");
      h_gem_sp_dy_adc_punchthrough[i_gem][i_plane]->Write();
      outputFile->cd("sp/dxdy/all");
      h_gem_sp_dxy[i_gem][i_plane]->Write();
      outputFile->cd("sp/dxdy/punchthrough");
      h_gem_sp_dxy_punchthrough[i_gem][i_plane]->Write();
    }
    // Write time histograms
    outputFile->cd(Form("time/%s", plane_names[i_plane].c_str()));
    h_hodo_time[i_plane]->Write();
    h_hodo_time_GEM0_x[i_plane]->Write();
    h_hodo_time_GEM0_y[i_plane]->Write();
    h_hodo_time_GEM1_x[i_plane]->Write();
    h_hodo_time_GEM1_y[i_plane]->Write();
    h_hodo_time_GEM0_x_punchthrough[i_plane]->Write();
    h_hodo_time_GEM0_y_punchthrough[i_plane]->Write();
    h_hodo_time_GEM1_x_punchthrough[i_plane]->Write();
    h_hodo_time_GEM1_y_punchthrough[i_plane]->Write();
    h_hodo_time_GEM0_xy[i_plane]->Write();
    h_hodo_time_GEM1_xy[i_plane]->Write();
    h_hodo_time_GEM_all_x[i_plane]->Write();
    h_hodo_time_GEM_all_y[i_plane]->Write();
    h_hodo_time_GEM_all_xy[i_plane]->Write();
    h_hodo_time_GEM0_xy_punchthrough[i_plane]->Write();
    h_hodo_time_GEM1_xy_punchthrough[i_plane]->Write();
    h_hodo_time_GEM_all_x_punchthrough[i_plane]->Write();
    h_hodo_time_GEM_all_y_punchthrough[i_plane]->Write();
    h_hodo_time_GEM_all_xy_punchthrough[i_plane]->Write();
  }
  // Create directories for histograms

  // Combine histograms for all planes
  // for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {

  // }

  // Close the files
  file->Close();
  outputFile->Close();
  std::cout << "Output file created: " << outputFileName << std::endl;
  return;
}