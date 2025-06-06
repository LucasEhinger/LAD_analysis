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

const double TDC2NS     = 0.09766; // ns per TDC channel
const char spect_prefix = 'H';     // Spectrometer prefix, 'H' or 'P'

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
const hist_range LAD_ADC(0, 300, 100);

const hist_range lad_time(1500, 2100, 100);
const hist_range time_diff(-2, 8, 100);
const hist_range time_peak(1770, 1790, 100);
const hist_range time_before(1700, 1750, 100);
const hist_range time_after(1850, 1950, 100);
const hist_range time_window_sig(1625, 1800, 100);
const hist_range time_window_bkd(1800, 1975, 100);

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
const int MINT_EVTS_PER_THREAD    = 10000;

const bool use_projz = true; // Fix z position to GEM projz

const double dx_min = -70.0;
const double dx_max = 70.0;
const int dx_NBINS  = 40;
const double dy_min = -40.0;
const double dy_max = 40.0;
const int dy_NBINS  = 40;

const int n_dxy_cuts             = 5;                // Number of dxy cuts
const double dx_cuts[n_dxy_cuts] = {1, 2, 4, 6, 10}; // dxy cuts in mm
const double dy_cuts[n_dxy_cuts] = {1, 2, 4, 6, 10}; // dxy cuts in mm

double gem_theta = 127.0;            // Angle in degrees
double gem_phi   = 0.0;              // Angle in degrees
double gem_r[2]  = {77.571, 95.571}; // Radius of the GEM's
double gem_dx[2] = {0.0, 0.0};       // GEM dx offsets
double gem_dy[2] = {0.0, 0.0};       // GEM dy offsets

template <typename T>
void drawEfficiencyHistograms(T *h_all_sig, T *h_all_bkd, T *h_hits_sig, T *h_hits_bkd, const double bkdSub_scale,
                              const char *title, const char *name, TFile *outfile, const char *dir_nosub, const char *dir_bkdSub) {

  // Draw histogram without background subtraction
  outfile->cd(dir_nosub);
  T* h_eff_sig = (T *)h_hits_sig->Clone(Form("h_efficiency_%s", name));
  h_eff_sig->SetTitle(Form("%s Efficiency (no bkd sub)", title));
  h_eff_sig->Divide(h_all_sig);
  h_eff_sig->SetMinimum(0);
  h_eff_sig->SetMaximum(2);
  h_eff_sig->Write();


  // Draw histogram with background subtraction
  outfile->cd(dir_bkdSub);
  h_all_sig->Add(h_all_bkd, -bkdSub_scale);
  h_hits_sig->Add(h_hits_bkd, -bkdSub_scale);

  T *h_eff_sig_bkdSub = (T *)h_all_sig->Clone(Form("h_efficiency_ratio_%s", name));
  h_eff_sig_bkdSub->SetTitle(Form("%s Efficiency (bkd sub)", title));
  h_eff_sig_bkdSub->Divide(h_hits_sig);
  h_eff_sig_bkdSub->SetMinimum(0);
  h_eff_sig_bkdSub->SetMaximum(2);
  h_eff_sig_bkdSub->Write();
}

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

void process_chunk(int thread_id, int start, int end, std::vector<TString> &fileNames,
                   map<string, TH1 *[2][nPlanes]> &hist_map,
                   map<string, TH1 *[n_dxy_cuts][nPlanes]> &hist_map_dxdy_cuts) {
  TChain *T = new TChain("T");
  for (const auto &fileName : fileNames) {
    T->Add(fileName);
  }
  //////////////////////////////////////////////////////////
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

  Double_t time_avg[nPlanes][MAX_DATA], time_avg_paddle[nPlanes][MAX_DATA], time_avg_ypos[nPlanes][MAX_DATA],
      adc_amp_avg[nPlanes][MAX_DATA];
  Int_t nData_hodo[nPlanes];

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
    T->SetBranchAddress(Form("%c.ladhod.%s.HodoHitEdepAmp", spect_prefix, plane_names[i].c_str()), &adc_amp_avg[i]);
  }
  //////////////////////////////////////////////////////////
  // Loop over the entries in the chunk
  Double_t hodo_hit_time_punchthrough[nPlanes][nPaddles], hodo_hit_adc_punchthrough[nPlanes][nPaddles];
  TVector3 p1[2], p2[2], p3[2];
  for (int i = 0; i < 2; ++i) {
    double theta = gem_theta;
    p1[i]        = TVector3(gem_r[i] * cos((theta - 90) * TMath::DegToRad()), 0,
                            -gem_r[i] * sin((theta - 90) * TMath::DegToRad()));
    p2[i]        = p1[i] + TVector3(-gem_r[i] * cos((180 - theta) * TMath::DegToRad()), 0,
                                    -gem_r[i] * sin((180 - theta) * TMath::DegToRad()));
    p3[i]        = p1[i] + TVector3(0, 10.0, 0); // Arbitrary y value of 10.0
  }

  for (int i = start; i < end; ++i) {
    T->GetEntry(i);
    if (abs(vertex_x) > 50 || abs(vertex_y) > 10 || abs(vertex_z) > 10) {
      continue;
    }
    /////////////////////////////////////////////////////////
    // Check for matching hits in the hodoscope
    for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
      for (int i_paddle = 0; i_paddle < nPaddles; ++i_paddle) {
        hodo_hit_time_punchthrough[i_plane][i_paddle] = -999;
        hodo_hit_adc_punchthrough[i_plane][i_paddle]  = -999;
      }
    }

    for (int plane = 0; plane < nPlanes; ++plane) {
      for (int i_hit = 0; i_hit < nData_hodo[plane]; ++i_hit) {
        int bar = int(time_avg_paddle[plane][i_hit]) - 1;
        if (bar < 0 || bar >= nPaddles)
          continue; // Skip invalid bars
        hodo_hit_time_punchthrough[plane][bar] = time_avg[plane][i_hit];
        hodo_hit_adc_punchthrough[plane][bar]  = adc_amp_avg[plane][i_hit];
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

        bool is_punchthrough = (hodo_hit_time_punchthrough[matching_plane][paddle] != -999);
        double diff_time =
            hodo_hit_time_punchthrough[matching_plane][paddle] - hodo_hit_time_punchthrough[i_plane][paddle];
        if (matching_plane < i_plane) {
          diff_time = -diff_time; // Adjust the time difference based on the plane order
        }

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

        TVector3 gem_hit_center_0   = GetGEMHitPosition(0, 0, 0, 0, gem_r[0], gem_theta);
        TVector3 gem_hit_center_1   = GetGEMHitPosition(0, 0, 0, 0, gem_r[1], gem_theta);
        Double_t golden_track_dy[2] = {gem_hit_center_0.Y() - intersection[0].Y(),
                                       gem_hit_center_1.Y() - intersection[1].Y()};
        // Calculate the difference in the x-z plane between the GEM center and the intersection
        Double_t golden_track_dxz[2] = {(gem_hit_center_0.X() > intersection[0].X() ? -1 : 1) *
                                            sqrt(pow(gem_hit_center_0.X() - intersection[0].X(), 2) +
                                                 pow(gem_hit_center_0.Z() - intersection[0].Z(), 2)),
                                        (gem_hit_center_1.X() > intersection[1].X() ? -1 : 1) *
                                            sqrt(pow(gem_hit_center_1.X() - intersection[1].X(), 2) +
                                                 pow(gem_hit_center_1.Z() - intersection[1].Z(), 2))};

        // intersection now contains the intersection point of the track (vertex->hodo) with the hodoscope plane

        double clust_min_dx[2]              = {9999, 9999};
        double clust_min_dy[2]              = {9999, 9999};
        double clust_min_x[2]               = {9999, 9999};
        double clust_min_y[2]               = {9999, 9999};
        double clust_min_dx_punchthrough[2] = {9999, 9999};
        double clust_min_dy_punchthrough[2] = {9999, 9999};
        double clust_min_x_punchthrough[2]  = {9999, 9999};
        double clust_min_y_punchthrough[2]  = {9999, 9999};

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
            hist_map["h_gem_cust_dx"][clust_layer_id][i_plane]->Fill(dxz);
            hist_map["h_gem_clust_dx_x"][clust_layer_id][i_plane]->Fill(clust_position, dxz);
            hist_map["h_gem_clust_dx_adc"][clust_layer_id][i_plane]->Fill(dxz, clust_adc_val);
            if (abs(clust_min_dx[clust_layer_id]) > abs(dxz)) {
              clust_min_dx[clust_layer_id] = dxz;
              clust_min_x[clust_layer_id]  = clust_position;
            }
            if (is_punchthrough) {
              hist_map["h_gem_cust_dx_punchthrough"][clust_layer_id][i_plane]->Fill(dxz);
              hist_map["h_gem_clust_dx_x_punchthrough"][clust_layer_id][i_plane]->Fill(clust_position, dxz);
              hist_map["h_gem_clust_dx_adc_punchthrough"][clust_layer_id][i_plane]->Fill(dxz, clust_adc_val);
              if (abs(clust_min_dx_punchthrough[clust_layer_id]) > abs(dxz)) {
                clust_min_dx_punchthrough[clust_layer_id] = dxz;
                clust_min_x_punchthrough[clust_layer_id]  = clust_position;
              }
            }
          } else if (clust_axis_id == 0) { // y axis
            gem_hit   = GetGEMHitPosition(0, clust_position, gem_dx[clust_layer_id], gem_dy[clust_layer_id],
                                          gem_r[clust_layer_id], gem_theta);
            double dy = gem_hit.Y() - intersection[clust_layer_id].Y();
            dy        = (-clust_position) -
                 intersection[clust_layer_id].Y(); // Adjust dy to be relative to the cluster position
            // LHE. FIXME. Something seems up with dy calculation.
            hist_map["h_gem_cust_dy"][clust_layer_id][i_plane]->Fill(dy);
            hist_map["h_gem_clust_dy_y"][clust_layer_id][i_plane]->Fill(clust_position, dy);
            hist_map["h_gem_clust_dy_adc"][clust_layer_id][i_plane]->Fill(dy, clust_adc_val);
            if (abs(clust_min_dy[clust_layer_id]) > abs(dy)) {
              clust_min_dy[clust_layer_id] = dy;
              clust_min_y[clust_layer_id]  = clust_position;
            }
            if (is_punchthrough) {
              hist_map["h_gem_cust_dy_punchthrough"][clust_layer_id][i_plane]->Fill(dy);
              hist_map["h_gem_clust_dy_y_punchthrough"][clust_layer_id][i_plane]->Fill(clust_position, dy);
              hist_map["h_gem_clust_dy_adc_punchthrough"][clust_layer_id][i_plane]->Fill(dy, clust_adc_val);
              if (abs(clust_min_dy_punchthrough[clust_layer_id]) > abs(dy)) {
                clust_min_dy_punchthrough[clust_layer_id] = dy;
                clust_min_y_punchthrough[clust_layer_id]  = clust_position;
              }
            }
          } else {
            std::cerr << "Error: Invalid cluster axis ID!" << std::endl;
            continue;
          }
        }
        // Fill the minimum dx and dy histograms
        hist_map["h_gem_cust_dx_min"][0][i_plane]->Fill(clust_min_dx[0]);
        hist_map["h_gem_cust_dy_min"][0][i_plane]->Fill(clust_min_dy[0]);
        hist_map["h_gem_cust_dx_min"][1][i_plane]->Fill(clust_min_dx[1]);
        hist_map["h_gem_cust_dy_min"][1][i_plane]->Fill(clust_min_dy[1]);
        if (time_id == 0) { // before
          hist_map["h_gem_cust_dx_min_before"][0][i_plane]->Fill(clust_min_dx[0]);
          hist_map["h_gem_cust_dy_min_before"][0][i_plane]->Fill(clust_min_dy[0]);
          hist_map["h_gem_cust_dx_min_before"][1][i_plane]->Fill(clust_min_dx[1]);
          hist_map["h_gem_cust_dy_min_before"][1][i_plane]->Fill(clust_min_dy[1]);
        } else if (time_id == 1) { // peak
          hist_map["h_gem_cust_dx_min_peak"][0][i_plane]->Fill(clust_min_dx[0]);
          hist_map["h_gem_cust_dy_min_peak"][0][i_plane]->Fill(clust_min_dy[0]);
          hist_map["h_gem_cust_dx_min_peak"][1][i_plane]->Fill(clust_min_dx[1]);
          hist_map["h_gem_cust_dy_min_peak"][1][i_plane]->Fill(clust_min_dy[1]);
        } else if (time_id == 2) { // after
          hist_map["h_gem_cust_dx_min_after"][0][i_plane]->Fill(clust_min_dx[0]);
          hist_map["h_gem_cust_dy_min_after"][0][i_plane]->Fill(clust_min_dy[0]);
          hist_map["h_gem_cust_dx_min_after"][1][i_plane]->Fill(clust_min_dx[1]);
          hist_map["h_gem_cust_dy_min_after"][1][i_plane]->Fill(clust_min_dy[1]);
        }

        // Fill the punchthrough minimum dx and dy histograms
        hist_map["h_gem_cust_dx_punchthrough_min"][0][i_plane]->Fill(clust_min_dx_punchthrough[0]);
        hist_map["h_gem_cust_dy_punchthrough_min"][0][i_plane]->Fill(clust_min_dy_punchthrough[0]);
        hist_map["h_gem_cust_dx_punchthrough_min"][1][i_plane]->Fill(clust_min_dx_punchthrough[1]);
        hist_map["h_gem_cust_dy_punchthrough_min"][1][i_plane]->Fill(clust_min_dy_punchthrough[1]);
        if (time_id == 0) { // before
          hist_map["h_gem_cust_dx_punchthrough_min_before"][0][i_plane]->Fill(clust_min_dx_punchthrough[0]);
          hist_map["h_gem_cust_dy_punchthrough_min_before"][0][i_plane]->Fill(clust_min_dy_punchthrough[0]);
          hist_map["h_gem_cust_dx_punchthrough_min_before"][1][i_plane]->Fill(clust_min_dx_punchthrough[1]);
          hist_map["h_gem_cust_dy_punchthrough_min_before"][1][i_plane]->Fill(clust_min_dy_punchthrough[1]);
        } else if (time_id == 1) { // peak
          hist_map["h_gem_cust_dx_punchthrough_min_peak"][0][i_plane]->Fill(clust_min_dx_punchthrough[0]);
          hist_map["h_gem_cust_dy_punchthrough_min_peak"][0][i_plane]->Fill(clust_min_dy_punchthrough[0]);
          hist_map["h_gem_cust_dx_punchthrough_min_peak"][1][i_plane]->Fill(clust_min_dx_punchthrough[1]);
          hist_map["h_gem_cust_dy_punchthrough_min_peak"][1][i_plane]->Fill(clust_min_dy_punchthrough[1]);
        } else if (time_id == 2) { // after
          hist_map["h_gem_cust_dx_punchthrough_min_after"][0][i_plane]->Fill(clust_min_dx_punchthrough[0]);
          hist_map["h_gem_cust_dy_punchthrough_min_after"][0][i_plane]->Fill(clust_min_dy_punchthrough[0]);
          hist_map["h_gem_cust_dx_punchthrough_min_after"][1][i_plane]->Fill(clust_min_dx_punchthrough[1]);
          hist_map["h_gem_cust_dy_punchthrough_min_after"][1][i_plane]->Fill(clust_min_dy_punchthrough[1]);
        }

        // Fill time histograms
        for (int i_cut = 0; i_cut < n_dxy_cuts; ++i_cut) {
          double max_dx = dx_cuts[i_cut];
          double max_dy = dy_cuts[i_cut];
          hist_map_dxdy_cuts["h_hodo_time"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
          hist_map_dxdy_cuts["h_hodo_adc_amp"][i_cut][i_plane]->Fill(adc_amp_avg[i_plane][i_hit]);

          if (is_punchthrough) {
            hist_map_dxdy_cuts["h_hodo_time_punchthrough"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
            hist_map_dxdy_cuts["h_hodo_time_diff_punchthrough"][i_cut][i_plane]->Fill(diff_time);
            hist_map_dxdy_cuts["h_hodo_adc_amp_punchthrough"][i_cut][i_plane]->Fill(adc_amp_avg[i_plane][i_hit]);
          }
          if (fabs(clust_min_dx[0]) < max_dx && fabs(clust_min_dy[0]) < max_dy) {
            hist_map_dxdy_cuts["h_hodo_time_GEM0_xy"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
            if (is_punchthrough) {
              hist_map_dxdy_cuts["h_hodo_time_GEM0_xy_punchthrough"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
              hist_map_dxdy_cuts["h_hodo_time_diff_GEM0_xy_punchthrough"][i_cut][i_plane]->Fill(diff_time);
              hist_map_dxdy_cuts["h_hodo_adc_amp_GEM0_xy_punchthrough"][i_cut][i_plane]->Fill(
                  adc_amp_avg[i_plane][i_hit]);
            }
          }
          if (fabs(clust_min_dx[1]) < max_dx && fabs(clust_min_dy[1]) < max_dy) {
            hist_map_dxdy_cuts["h_hodo_time_GEM1_xy"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
            if (is_punchthrough) {
              hist_map_dxdy_cuts["h_hodo_time_GEM1_xy_punchthrough"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
              hist_map_dxdy_cuts["h_hodo_time_diff_GEM1_xy_punchthrough"][i_cut][i_plane]->Fill(diff_time);
              hist_map_dxdy_cuts["h_hodo_adc_amp_GEM1_xy_punchthrough"][i_cut][i_plane]->Fill(
                  adc_amp_avg[i_plane][i_hit]);
            }
          }
          if (fabs(clust_min_dx[0]) < max_dx) {
            hist_map_dxdy_cuts["h_hodo_time_GEM0_x"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
            if (is_punchthrough) {
              hist_map_dxdy_cuts["h_hodo_time_GEM0_x_punchthrough"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
              hist_map_dxdy_cuts["h_hodo_time_diff_GEM0_x_punchthrough"][i_cut][i_plane]->Fill(diff_time);
              hist_map_dxdy_cuts["h_hodo_adc_amp_GEM0_x_punchthrough"][i_cut][i_plane]->Fill(
                  adc_amp_avg[i_plane][i_hit]);
            }
          }
          if (fabs(clust_min_dy[0]) < max_dy) {
            hist_map_dxdy_cuts["h_hodo_time_GEM0_y"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
            if (is_punchthrough) {
              hist_map_dxdy_cuts["h_hodo_time_GEM0_y_punchthrough"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
              hist_map_dxdy_cuts["h_hodo_time_diff_GEM0_y_punchthrough"][i_cut][i_plane]->Fill(diff_time);
              hist_map_dxdy_cuts["h_hodo_adc_amp_GEM0_y_punchthrough"][i_cut][i_plane]->Fill(
                  adc_amp_avg[i_plane][i_hit]);
            }
          }
          if (fabs(clust_min_dx[1]) < max_dx) {
            hist_map_dxdy_cuts["h_hodo_time_GEM1_x"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
            if (is_punchthrough) {
              hist_map_dxdy_cuts["h_hodo_time_GEM1_x_punchthrough"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
              hist_map_dxdy_cuts["h_hodo_time_diff_GEM1_x_punchthrough"][i_cut][i_plane]->Fill(diff_time);
              hist_map_dxdy_cuts["h_hodo_adc_amp_GEM1_x_punchthrough"][i_cut][i_plane]->Fill(
                  adc_amp_avg[i_plane][i_hit]);
            }
          }
          if (fabs(clust_min_dy[1]) < max_dy) {
            hist_map_dxdy_cuts["h_hodo_time_GEM1_y"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
            if (is_punchthrough) {
              hist_map_dxdy_cuts["h_hodo_time_GEM1_y_punchthrough"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
              hist_map_dxdy_cuts["h_hodo_time_diff_GEM1_y_punchthrough"][i_cut][i_plane]->Fill(diff_time);
              hist_map_dxdy_cuts["h_hodo_adc_amp_GEM1_y_punchthrough"][i_cut][i_plane]->Fill(
                  adc_amp_avg[i_plane][i_hit]);
            }
          }
          if (fabs(clust_min_dx[0]) < max_dx && fabs(clust_min_dx[1]) < max_dx) {
            hist_map_dxdy_cuts["h_hodo_time_GEM_all_x"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
            if (is_punchthrough) {
              hist_map_dxdy_cuts["h_hodo_time_GEM_all_x_punchthrough"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
              hist_map_dxdy_cuts["h_hodo_time_diff_GEM_all_x_punchthrough"][i_cut][i_plane]->Fill(diff_time);
              hist_map_dxdy_cuts["h_hodo_adc_amp_GEM_all_x_punchthrough"][i_cut][i_plane]->Fill(
                  adc_amp_avg[i_plane][i_hit]);
            }
          }
          if (fabs(clust_min_dy[0]) < max_dy && fabs(clust_min_dy[1]) < max_dy) {
            hist_map_dxdy_cuts["h_hodo_time_GEM_all_y"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
            if (is_punchthrough) {
              hist_map_dxdy_cuts["h_hodo_time_GEM_all_y_punchthrough"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
              hist_map_dxdy_cuts["h_hodo_time_diff_GEM_all_y_punchthrough"][i_cut][i_plane]->Fill(diff_time);
              hist_map_dxdy_cuts["h_hodo_adc_amp_GEM_all_y_punchthrough"][i_cut][i_plane]->Fill(
                  adc_amp_avg[i_plane][i_hit]);
            }
          }
          if (fabs(clust_min_dx[0]) < max_dx && fabs(clust_min_dy[0]) < max_dy && fabs(clust_min_dx[1]) < max_dx &&
              fabs(clust_min_dy[1]) < max_dy) {
            hist_map_dxdy_cuts["h_hodo_time_GEM_all_xy"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
            if (is_punchthrough) {
              hist_map_dxdy_cuts["h_hodo_time_GEM_all_xy_punchthrough"][i_cut][i_plane]->Fill(time_avg[i_plane][i_hit]);
              hist_map_dxdy_cuts["h_hodo_time_diff_GEM_all_xy_punchthrough"][i_cut][i_plane]->Fill(diff_time);
              hist_map_dxdy_cuts["h_hodo_adc_amp_GEM_all_xy_punchthrough"][i_cut][i_plane]->Fill(
                  adc_amp_avg[i_plane][i_hit]);
            }
          }

          // GEM Efficiency histograms
          if (fabs(clust_min_dx[0]) < max_dx) {
            if (time_avg[i_plane][i_hit] >= time_window_sig.min && time_avg[i_plane][i_hit] < time_window_sig.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_all_GEM1_x"][i_cut][i_plane]->Fill(golden_track_dxz[1]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_all_GEM1_x_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[1]);
              }
            }
            if (time_avg[i_plane][i_hit] >= time_window_bkd.min && time_avg[i_plane][i_hit] < time_window_bkd.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_all_GEM1_x_bkd"][i_cut][i_plane]->Fill(golden_track_dxz[1]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_all_GEM1_x_bkd_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[1]);
              }
            }
          }

          if (fabs(clust_min_dy[0]) < max_dy) {
            if (time_avg[i_plane][i_hit] >= time_window_sig.min && time_avg[i_plane][i_hit] < time_window_sig.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_all_GEM1_y"][i_cut][i_plane]->Fill(golden_track_dy[1]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_all_GEM1_y_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dy[1]);
              }
            }
            if (time_avg[i_plane][i_hit] >= time_window_bkd.min && time_avg[i_plane][i_hit] < time_window_bkd.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_all_GEM1_y_bkd"][i_cut][i_plane]->Fill(golden_track_dy[1]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_all_GEM1_y_bkd_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dy[1]);
              }
            }
          }

          if (fabs(clust_min_dx[1]) < max_dx) {
            if (time_avg[i_plane][i_hit] >= time_window_sig.min && time_avg[i_plane][i_hit] < time_window_sig.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_all_GEM0_x"][i_cut][i_plane]->Fill(golden_track_dxz[0]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_all_GEM0_x_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[0]);
              }
            }
            if (time_avg[i_plane][i_hit] >= time_window_bkd.min && time_avg[i_plane][i_hit] < time_window_bkd.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_all_GEM0_x_bkd"][i_cut][i_plane]->Fill(golden_track_dxz[0]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_all_GEM0_x_bkd_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[0]);
              }
            }
          }

          if (fabs(clust_min_dy[1]) < max_dy) {
            if (time_avg[i_plane][i_hit] >= time_window_sig.min && time_avg[i_plane][i_hit] < time_window_sig.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_all_GEM0_y"][i_cut][i_plane]->Fill(golden_track_dy[0]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_all_GEM0_y_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dy[0]);
              }
            }
            if (time_avg[i_plane][i_hit] >= time_window_bkd.min && time_avg[i_plane][i_hit] < time_window_bkd.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_all_GEM0_y_bkd"][i_cut][i_plane]->Fill(golden_track_dy[0]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_all_GEM0_y_bkd_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dy[0]);
              }
            }
          }

          if (fabs(clust_min_dx[0]) < max_dx && fabs(clust_min_dx[1]) < max_dx) {
            if (time_avg[i_plane][i_hit] >= time_window_sig.min && time_avg[i_plane][i_hit] < time_window_sig.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM0_x"][i_cut][i_plane]->Fill(golden_track_dxz[0]);
              hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM1_x"][i_cut][i_plane]->Fill(golden_track_dxz[1]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM0_x_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[0]);
                hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM1_x_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[1]);
              }
            }
            if (time_avg[i_plane][i_hit] >= time_window_bkd.min && time_avg[i_plane][i_hit] < time_window_bkd.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM0_x_bkd"][i_cut][i_plane]->Fill(golden_track_dxz[0]);
              hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM1_x_bkd"][i_cut][i_plane]->Fill(golden_track_dxz[1]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM0_x_bkd_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[0]);
                hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM1_x_bkd_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[1]);
              }
            }
          }
          if (fabs(clust_min_dy[0]) < max_dy && fabs(clust_min_dy[1]) < max_dy) {
            if (time_avg[i_plane][i_hit] >= time_window_sig.min && time_avg[i_plane][i_hit] < time_window_sig.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM0_y"][i_cut][i_plane]->Fill(golden_track_dy[0]);
              hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM1_y"][i_cut][i_plane]->Fill(golden_track_dy[1]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM0_y_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dy[0]);
                hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM1_y_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dy[1]);
              }
            }
            if (time_avg[i_plane][i_hit] >= time_window_bkd.min && time_avg[i_plane][i_hit] < time_window_bkd.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM0_y_bkd"][i_cut][i_plane]->Fill(golden_track_dy[0]);
              hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM1_y_bkd"][i_cut][i_plane]->Fill(golden_track_dy[1]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM0_y_bkd_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dy[0]);
                hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM1_y_bkd_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dy[1]);
              }
            }
          }

          if (fabs(clust_min_dx[0]) < max_dx && fabs(clust_min_dy[0]) < max_dy) {
            if (time_avg[i_plane][i_hit] >= time_window_sig.min && time_avg[i_plane][i_hit] < time_window_sig.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_all_GEM1_2D"][i_cut][i_plane]->Fill(golden_track_dxz[0],
                                                                                       golden_track_dy[0]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_all_GEM1_2D_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[0], golden_track_dy[0]);
              }
            }
            if (time_avg[i_plane][i_hit] >= time_window_bkd.min && time_avg[i_plane][i_hit] < time_window_bkd.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_all_GEM1_2D_bkd"][i_cut][i_plane]->Fill(golden_track_dxz[0],
                                                                                           golden_track_dy[0]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_all_GEM1_2D_bkd_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[0], golden_track_dy[0]);
              }
            }
          }
          if (fabs(clust_min_dx[1]) < max_dx && fabs(clust_min_dy[1]) < max_dy) {
            if (time_avg[i_plane][i_hit] >= time_window_sig.min && time_avg[i_plane][i_hit] < time_window_sig.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_all_GEM0_2D"][i_cut][i_plane]->Fill(golden_track_dxz[1],
                                                                                       golden_track_dy[1]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_all_GEM0_2D_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[1], golden_track_dy[1]);
              }
            }
            if (time_avg[i_plane][i_hit] >= time_window_bkd.min && time_avg[i_plane][i_hit] < time_window_bkd.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_all_GEM0_2D_bkd"][i_cut][i_plane]->Fill(golden_track_dxz[1],
                                                                                           golden_track_dy[1]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_all_GEM0_2D_bkd_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[1], golden_track_dy[1]);
              }
            }
          }

          if (fabs(clust_min_dx[0]) < max_dx && fabs(clust_min_dx[1]) < max_dx && fabs(clust_min_dy[0]) < max_dy &&
              fabs(clust_min_dy[1]) < max_dy) {
            if (time_avg[i_plane][i_hit] >= time_window_sig.min && time_avg[i_plane][i_hit] < time_window_sig.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM0_2D"][i_cut][i_plane]->Fill(golden_track_dxz[0],
                                                                                            golden_track_dy[0]);
              hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM1_2D"][i_cut][i_plane]->Fill(golden_track_dxz[1],
                                                                                            golden_track_dy[1]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM0_2D_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[0], golden_track_dy[0]);
                hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM1_2D_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[1], golden_track_dy[1]);
              }
            }
            if (time_avg[i_plane][i_hit] >= time_window_bkd.min && time_avg[i_plane][i_hit] < time_window_bkd.max) {
              hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM0_2D_bkd"][i_cut][i_plane]->Fill(golden_track_dxz[0],
                                                                                                golden_track_dy[0]);
              hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM1_2D_bkd"][i_cut][i_plane]->Fill(golden_track_dxz[1],
                                                                                                golden_track_dy[1]);
              if (is_punchthrough) {
                hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM0_2D_bkd_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[0], golden_track_dy[0]);
                hist_map_dxdy_cuts["h_gem_efficiency_good_hit_GEM1_2D_bkd_punchthrough"][i_cut][i_plane]->Fill(
                    golden_track_dxz[1], golden_track_dy[1]);
              }
            }
          }
        } // End Loop over cuts

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

          hist_map["h_gem_sp_dx"][sp_layer_id][i_plane]->Fill(dxz);
          hist_map["h_gem_sp_dy"][sp_layer_id][i_plane]->Fill(dy);
          hist_map["h_gem_sp_dx_x"][sp_layer_id][i_plane]->Fill(sp_X[i_sp], dxz);
          hist_map["h_gem_sp_dy_y"][sp_layer_id][i_plane]->Fill(sp_Y[i_sp], dy);
          hist_map["h_gem_sp_dx_adc"][sp_layer_id][i_plane]->Fill(dxz, sp_adc_val);
          hist_map["h_gem_sp_dy_adc"][sp_layer_id][i_plane]->Fill(dy, sp_adc_val);
          hist_map["h_gem_sp_dxy"][sp_layer_id][i_plane]->Fill(dxz, dy);
          if (is_punchthrough) {
            hist_map["h_gem_sp_dx_punchthrough"][sp_layer_id][i_plane]->Fill(dxz);
            hist_map["h_gem_sp_dy_punchthrough"][sp_layer_id][i_plane]->Fill(dy);
            hist_map["h_gem_sp_dx_x_punchthrough"][sp_layer_id][i_plane]->Fill(sp_X[i_sp], dxz);
            hist_map["h_gem_sp_dy_y_punchthrough"][sp_layer_id][i_plane]->Fill(sp_Y[i_sp], dy);
            hist_map["h_gem_sp_dx_adc_punchthrough"][sp_layer_id][i_plane]->Fill(dxz, sp_adc_val);
            hist_map["h_gem_sp_dy_adc_punchthrough"][sp_layer_id][i_plane]->Fill(dy, sp_adc_val);
            hist_map["h_gem_sp_dxy_punchthrough"][sp_layer_id][i_plane]->Fill(dxz, dy);
          }
        }
      }
    }

    if (i % ((end - start) / 100) == 0 && thread_id == 0) {
      std::cout << "\rProcessing: " << int((i - start) * 100.0 / (end - start)) << "% completed." << std::flush;
    }

  } // End Event Loop
}

void gem_tracking_by_clust() {
  // Set ROOT to batch mode to avoid opening canvases
  gROOT->SetBatch(kTRUE);
  ROOT::EnableThreadSafety();
  std::vector<TString> fileNames = {
      Form("/volatile/hallc/c-lad/ehingerl/lad_replay/ROOTfiles/LAD_COIN/PRODUCTION/bad_timing/"
           "LAD_COIN_22615_0_6_-1.root")
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
  TString outputFileName = Form("files/tracking_gem/tracking_gem_22615_newtiming_-1_%c.root", spect_prefix);
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

  int nEntries = T->GetEntries();
  // Number of threads to use
  int numThreads = std::thread::hardware_concurrency();
  // numThreads     = 1;
  int chunkSize = nEntries / numThreads;

  // Adjust the number of threads if the chunk size is too small
  if (chunkSize < MINT_EVTS_PER_THREAD) {
    numThreads = std::max(1, nEntries / MINT_EVTS_PER_THREAD);
    chunkSize  = nEntries / numThreads;
  }
  // numThreads = 1;
  // chunkSize  = nEntries / numThreads;

  vector<map<string, TH1 *[2][nPlanes]>> hist_map_vec(numThreads);
  vector<map<string, TH1 *[n_dxy_cuts][nPlanes]>> hist_map_dxdy_cuts_vec(numThreads);

  //////////////////////////////////////////////////////////////////
  for (int thread = 0; thread < numThreads; ++thread) {
    for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
      for (int i_gem = 0; i_gem < 2; ++i_gem) {
        hist_map_vec[thread]["h_gem_cust_dx"][i_gem][i_plane] = new TH1F(
            Form("h_gem_cust_dx_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_cust_dx_%d_%s;Cluster x - Golden Track x (cm);Counts", i_gem, plane_names[i_plane].c_str()),
            dx_NBINS * 2, dx_min * 2, dx_max * 2);
        hist_map_vec[thread]["h_gem_cust_dy"][i_gem][i_plane] = new TH1F(
            Form("h_gem_cust_dy_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_cust_dy_%d_%s;Cluster y - Golden Track y (cm);Counts", i_gem, plane_names[i_plane].c_str()),
            dy_NBINS * 2, dy_min * 2, dy_max * 2);

        hist_map_vec[thread]["h_gem_clust_dx_x"][i_gem][i_plane] =
            new TH2D(Form("h_gem_clust_dx_x_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_clust_dx_x_%d_%s;Cluster x;Cluster x - Golden Track x (cm)", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max, dx_NBINS * 2, dx_min * 2, dx_max * 2);

        hist_map_vec[thread]["h_gem_clust_dy_y"][i_gem][i_plane] =
            new TH2D(Form("h_gem_clust_dy_y_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_clust_dy_y_%d_%s;Cluster y;Cluster y - Golden Track y (cm)", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max, dy_NBINS * 2, dy_min * 2, dy_max * 2);

        hist_map_vec[thread]["h_gem_sp_dx"][i_gem][i_plane] = new TH1F(
            Form("h_gem_sp_dx_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_sp_dx_%d_%s;Cluster x - Golden Track x (cm);Counts", i_gem, plane_names[i_plane].c_str()),
            dx_NBINS * 2, dx_min * 2, dx_max * 2);
        hist_map_vec[thread]["h_gem_sp_dy"][i_gem][i_plane] = new TH1F(
            Form("h_gem_sp_dy_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_sp_dy_%d_%s;Cluster y - Golden Track y (cm);Counts", i_gem, plane_names[i_plane].c_str()),
            dy_NBINS * 2, dy_min * 2, dy_max * 2);
        hist_map_vec[thread]["h_gem_sp_dx_x"][i_gem][i_plane] = new TH2D(
            Form("h_gem_sp_dx_x_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_sp_dx_x_%d_%s;Cluster x;Cluster x - Golden Track x (cm)", i_gem, plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max, dx_NBINS * 2, dx_min * 2, dx_max * 2);
        hist_map_vec[thread]["h_gem_sp_dy_y"][i_gem][i_plane] = new TH2D(
            Form("h_gem_sp_dy_y_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_sp_dy_y_%d_%s;Cluster y;Cluster y - Golden Track y (cm)", i_gem, plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max, dy_NBINS * 2, dy_min * 2, dy_max * 2);

        hist_map_vec[thread]["h_gem_cust_dx_punchthrough"][i_gem][i_plane] =
            new TH1F(Form("h_gem_cust_dx_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_cust_dx_punchthrough_%d_%s;Cluster x - Golden Track x (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS * 2, dx_min * 2, dx_max * 2);
        hist_map_vec[thread]["h_gem_cust_dy_punchthrough"][i_gem][i_plane] =
            new TH1F(Form("h_gem_cust_dy_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_cust_dy_punchthrough_%d_%s;Cluster y - Golden Track y (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS * 2, dy_min * 2, dy_max * 2);
        hist_map_vec[thread]["h_gem_clust_dx_x_punchthrough"][i_gem][i_plane] =
            new TH2D(Form("h_gem_clust_dx_x_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_clust_dx_x_punchthrough_%d_%s;Cluster x;Cluster x - Golden Track x (cm)", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max, dx_NBINS * 2, dx_min * 2, dx_max * 2);
        hist_map_vec[thread]["h_gem_clust_dy_y_punchthrough"][i_gem][i_plane] =
            new TH2D(Form("h_gem_clust_dy_y_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_clust_dy_y_punchthrough_%d_%s;Cluster y;Cluster y - Golden Track y (cm)", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max, dy_NBINS * 2, dy_min * 2, dy_max * 2);

        hist_map_vec[thread]["h_gem_sp_dx_punchthrough"][i_gem][i_plane] =
            new TH1F(Form("h_gem_sp_dx_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dx_punchthrough_%d_%s;Cluster x - Golden Track x (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS * 2, dx_min * 2, dx_max * 2);
        hist_map_vec[thread]["h_gem_sp_dy_punchthrough"][i_gem][i_plane] =
            new TH1F(Form("h_gem_sp_dy_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dy_punchthrough_%d_%s;Cluster y - Golden Track y (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS * 2, dy_min * 2, dy_max * 2);
        hist_map_vec[thread]["h_gem_sp_dx_x_punchthrough"][i_gem][i_plane] =
            new TH2D(Form("h_gem_sp_dx_x_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dx_x_punchthrough_%d_%s;Cluster x;Cluster x - Golden Track x (cm)", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max, dx_NBINS * 2, dx_min * 2, dx_max * 2);
        hist_map_vec[thread]["h_gem_sp_dy_y_punchthrough"][i_gem][i_plane] =
            new TH2D(Form("h_gem_sp_dy_y_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dy_y_punchthrough_%d_%s;Cluster y;Cluster y - Golden Track y (cm)", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max, dy_NBINS * 2, dy_min * 2, dy_max * 2);

        hist_map_vec[thread]["h_gem_clust_dxy"][i_gem][i_plane] =
            new TH2D(Form("h_gem_clust_dxy_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_clust_dxy_%d_%s;Cluster x - Golden Track x (cm);Cluster y - Golden Track y (cm)",
                          i_gem, plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_sp_dxy"][i_gem][i_plane] =
            new TH2D(Form("h_gem_sp_dxy_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dxy_%d_%s;Cluster x - Golden Track x (cm);Cluster y - Golden Track y (cm)", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_clust_dxy_punchthrough"][i_gem][i_plane] = new TH2D(
            Form("h_gem_clust_dxy_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_clust_dxy_punchthrough_%d_%s;Cluster x - Golden Track x (cm);Cluster y - Golden Track y (cm)",
                 i_gem, plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_sp_dxy_punchthrough"][i_gem][i_plane] = new TH2D(
            Form("h_gem_sp_dxy_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_sp_dxy_punchthrough_%d_%s;Cluster x - Golden Track x (cm);Cluster y - Golden Track y (cm)",
                 i_gem, plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);

        hist_map_vec[thread]["h_gem_cust_dx_min"][i_gem][i_plane] =
            new TH1F(Form("h_gem_cust_dx_min_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_cust_dx_min_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_cust_dy_min"][i_gem][i_plane] =
            new TH1F(Form("h_gem_cust_dy_min_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_cust_dy_min_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_sp_dx_min"][i_gem][i_plane] =
            new TH1F(Form("h_gem_sp_dx_min_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dx_min_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_sp_dy_min"][i_gem][i_plane] =
            new TH1F(Form("h_gem_sp_dy_min_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dy_min_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);

        hist_map_vec[thread]["h_gem_cust_dx_min_before"][i_gem][i_plane] =
            new TH1F(Form("h_gem_cust_dx_min_before_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_cust_dx_min_before_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_cust_dy_min_before"][i_gem][i_plane] =
            new TH1F(Form("h_gem_cust_dy_min_before_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_cust_dy_min_before_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_sp_dx_min_before"][i_gem][i_plane] =
            new TH1F(Form("h_gem_sp_dx_min_before_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dx_min_before_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_sp_dy_min_before"][i_gem][i_plane] =
            new TH1F(Form("h_gem_sp_dy_min_before_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dy_min_before_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_cust_dx_min_peak"][i_gem][i_plane] =
            new TH1F(Form("h_gem_cust_dx_min_peak_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_cust_dx_min_peak_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_cust_dy_min_peak"][i_gem][i_plane] =
            new TH1F(Form("h_gem_cust_dy_min_peak_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_cust_dy_min_peak_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_sp_dx_min_peak"][i_gem][i_plane] =
            new TH1F(Form("h_gem_sp_dx_min_peak_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dx_min_peak_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_sp_dy_min_peak"][i_gem][i_plane] =
            new TH1F(Form("h_gem_sp_dy_min_peak_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dy_min_peak_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_cust_dx_min_after"][i_gem][i_plane] =
            new TH1F(Form("h_gem_cust_dx_min_after_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_cust_dx_min_after_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_cust_dy_min_after"][i_gem][i_plane] =
            new TH1F(Form("h_gem_cust_dy_min_after_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_cust_dy_min_after_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_sp_dx_min_after"][i_gem][i_plane] =
            new TH1F(Form("h_gem_sp_dx_min_after_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dx_min_after_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_sp_dy_min_after"][i_gem][i_plane] =
            new TH1F(Form("h_gem_sp_dy_min_after_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dy_min_after_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);

        hist_map_vec[thread]["h_gem_cust_dx_punchthrough_min"][i_gem][i_plane] = new TH1F(
            Form("h_gem_cust_dx_punchthrough_min_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_cust_dx_punchthrough_min_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_cust_dy_punchthrough_min"][i_gem][i_plane] = new TH1F(
            Form("h_gem_cust_dy_punchthrough_min_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_cust_dy_punchthrough_min_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_sp_dx_punchthrough_min"][i_gem][i_plane] =
            new TH1F(Form("h_gem_sp_dx_punchthrough_min_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dx_punchthrough_min_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_sp_dy_punchthrough_min"][i_gem][i_plane] =
            new TH1F(Form("h_gem_sp_dy_punchthrough_min_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dy_punchthrough_min_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_cust_dx_punchthrough_min_before"][i_gem][i_plane] = new TH1F(
            Form("h_gem_cust_dx_punchthrough_min_before_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_cust_dx_punchthrough_min_before_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_cust_dy_punchthrough_min_before"][i_gem][i_plane] = new TH1F(
            Form("h_gem_cust_dy_punchthrough_min_before_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_cust_dy_punchthrough_min_before_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_sp_dx_punchthrough_min_before"][i_gem][i_plane] = new TH1F(
            Form("h_gem_sp_dx_punchthrough_min_before_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_sp_dx_punchthrough_min_before_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_sp_dy_punchthrough_min_before"][i_gem][i_plane] = new TH1F(
            Form("h_gem_sp_dy_punchthrough_min_before_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_sp_dy_punchthrough_min_before_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_cust_dx_punchthrough_min_peak"][i_gem][i_plane] = new TH1F(
            Form("h_gem_cust_dx_punchthrough_min_peak_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_cust_dx_punchthrough_min_peak_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_cust_dy_punchthrough_min_peak"][i_gem][i_plane] = new TH1F(
            Form("h_gem_cust_dy_punchthrough_min_peak_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_cust_dy_punchthrough_min_peak_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_sp_dx_punchthrough_min_peak"][i_gem][i_plane] = new TH1F(
            Form("h_gem_sp_dx_punchthrough_min_peak_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_sp_dx_punchthrough_min_peak_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_sp_dy_punchthrough_min_peak"][i_gem][i_plane] = new TH1F(
            Form("h_gem_sp_dy_punchthrough_min_peak_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_sp_dy_punchthrough_min_peak_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_cust_dx_punchthrough_min_after"][i_gem][i_plane] = new TH1F(
            Form("h_gem_cust_dx_punchthrough_min_after_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_cust_dx_punchthrough_min_after_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_cust_dy_punchthrough_min_after"][i_gem][i_plane] = new TH1F(
            Form("h_gem_cust_dy_punchthrough_min_after_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_cust_dy_punchthrough_min_after_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);
        hist_map_vec[thread]["h_gem_sp_dx_punchthrough_min_after"][i_gem][i_plane] = new TH1F(
            Form("h_gem_sp_dx_punchthrough_min_after_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_sp_dx_punchthrough_min_after_%d_%s;Min |Cluster x - Golden Track x| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_vec[thread]["h_gem_sp_dy_punchthrough_min_after"][i_gem][i_plane] = new TH1F(
            Form("h_gem_sp_dy_punchthrough_min_after_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_sp_dy_punchthrough_min_after_%d_%s;Min |Cluster y - Golden Track y| (cm);Counts", i_gem,
                 plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);
        // ADC histograms
        hist_map_vec[thread]["h_gem_clust_dx_adc"][i_gem][i_plane] = new TH2D(
            Form("h_gem_clust_dx_adc_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_clust_dx_adc_%d_%s;Cluster x - Golden Track x (cm);ADC", i_gem, plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
        hist_map_vec[thread]["h_gem_clust_dy_adc"][i_gem][i_plane] = new TH2D(
            Form("h_gem_clust_dy_adc_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_clust_dy_adc_%d_%s;Cluster y - Golden Track y (cm);ADC", i_gem, plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
        hist_map_vec[thread]["h_gem_sp_dx_adc"][i_gem][i_plane] = new TH2D(
            Form("h_gem_sp_dx_adc_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_sp_dx_adc_%d_%s;Cluster x - Golden Track x (cm);ADC", i_gem, plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
        hist_map_vec[thread]["h_gem_sp_dy_adc"][i_gem][i_plane] = new TH2D(
            Form("h_gem_sp_dy_adc_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_sp_dy_adc_%d_%s;Cluster y - Golden Track y (cm);ADC", i_gem, plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
        hist_map_vec[thread]["h_gem_clust_dx_adc_punchthrough"][i_gem][i_plane] = new TH2D(
            Form("h_gem_clust_dx_adc_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_clust_dx_adc_punchthrough_%d_%s;Cluster x - Golden Track x (cm);ADC", i_gem,
                 plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
        hist_map_vec[thread]["h_gem_clust_dy_adc_punchthrough"][i_gem][i_plane] = new TH2D(
            Form("h_gem_clust_dy_adc_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
            Form("h_gem_clust_dy_adc_punchthrough_%d_%s;Cluster y - Golden Track y (cm);ADC", i_gem,
                 plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
        hist_map_vec[thread]["h_gem_sp_dx_adc_punchthrough"][i_gem][i_plane] =
            new TH2D(Form("h_gem_sp_dx_adc_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dx_adc_punchthrough_%d_%s;Cluster x - Golden Track x (cm);ADC", i_gem,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
        hist_map_vec[thread]["h_gem_sp_dy_adc_punchthrough"][i_gem][i_plane] =
            new TH2D(Form("h_gem_sp_dy_adc_punchthrough_%d_%s_thread_%d", i_gem, plane_names[i_plane].c_str(), thread),
                     Form("h_gem_sp_dy_adc_punchthrough_%d_%s;Cluster y - Golden Track y (cm);ADC", i_gem,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max, GEM_ADC.nbins, GEM_ADC.min, GEM_ADC.max);
      }
      // Hodo time histograms

      for (int i_cut = 0; i_cut < n_dxy_cuts; ++i_cut) {
        hist_map_dxdy_cuts_vec[thread]["h_hodo_time"][i_cut][i_plane] =
            new TH1F(Form("h_hodo_time_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_time_cut%d_%s;Hodoscope time (ns);Counts", i_cut, plane_names[i_plane].c_str()),
                     lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_punchthrough"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_punchthrough_cut%d_%s;Hodoscope time (ns);Counts", i_cut, plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_diff_punchthrough"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_diff_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_diff_punchthrough_cut%d_%s;Hodoscope time Difference (ns);Counts", i_cut,
                 plane_names[i_plane].c_str()),
            time_diff.nbins, time_diff.min, time_diff.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM0_x"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM0_x_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM0_x_cut%d_%s;Hodoscope time (ns);Counts", i_cut, plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM0_y"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM0_y_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM0_y_cut%d_%s;Hodoscope time (ns);Counts", i_cut, plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM1_x"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM1_x_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM1_x_cut%d_%s;Hodoscope time (ns);Counts", i_cut, plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM1_y"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM1_y_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM1_y_cut%d_%s;Hodoscope time (ns);Counts", i_cut, plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM0_x_punchthrough"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM0_x_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM0_x_punchthrough_cut%d_%s;Hodoscope time (ns);Counts", i_cut,
                 plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_diff_GEM0_x_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_hodo_time_diff_GEM0_x_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_time_diff_GEM0_x_punchthrough_cut%d_%s;Hodoscope time Difference (ns);Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     time_diff.nbins, time_diff.min, time_diff.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM0_y_punchthrough"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM0_y_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM0_y_punchthrough_cut%d_%s;Hodoscope time (ns);Counts", i_cut,
                 plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_diff_GEM0_y_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_hodo_time_diff_GEM0_y_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_time_diff_GEM0_y_punchthrough_cut%d_%s;Hodoscope time Difference (ns);Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     time_diff.nbins, time_diff.min, time_diff.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM1_x_punchthrough"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM1_x_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM1_x_punchthrough_cut%d_%s;Hodoscope time (ns);Counts", i_cut,
                 plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_diff_GEM1_x_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_hodo_time_diff_GEM1_x_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_time_diff_GEM1_x_punchthrough_cut%d_%s;Hodoscope time Difference (ns);Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     time_diff.nbins, time_diff.min, time_diff.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM1_y_punchthrough"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM1_y_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM1_y_punchthrough_cut%d_%s;Hodoscope time (ns);Counts", i_cut,
                 plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_diff_GEM1_y_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_hodo_time_diff_GEM1_y_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_time_diff_GEM1_y_punchthrough_cut%d_%s;Hodoscope time Difference (ns);Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     time_diff.nbins, time_diff.min, time_diff.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM0_xy"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM0_xy_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM0_xy_cut%d_%s;Hodoscope time (ns);Counts", i_cut, plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM1_xy"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM1_xy_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM1_xy_cut%d_%s;Hodoscope time (ns);Counts", i_cut, plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM_all_x"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM_all_x_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM_all_x_cut%d_%s;Hodoscope time (ns);Counts", i_cut, plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM_all_y"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM_all_y_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM_all_y_cut%d_%s;Hodoscope time (ns);Counts", i_cut, plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM_all_xy"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM_all_xy_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM_all_xy_cut%d_%s;Hodoscope time (ns);Counts", i_cut, plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM0_xy_punchthrough"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM0_xy_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM0_xy_punchthrough_cut%d_%s;Hodoscope time (ns);Counts", i_cut,
                 plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_diff_GEM0_xy_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_hodo_time_diff_GEM0_xy_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_time_diff_GEM0_xy_punchthrough_cut%d_%s;Hodoscope time Difference (ns);Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     time_diff.nbins, time_diff.min, time_diff.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM1_xy_punchthrough"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM1_xy_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM1_xy_punchthrough_cut%d_%s;Hodoscope time (ns);Counts", i_cut,
                 plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_diff_GEM1_xy_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_hodo_time_diff_GEM1_xy_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_time_diff_GEM1_xy_punchthrough_cut%d_%s;Hodoscope time Difference (ns);Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     time_diff.nbins, time_diff.min, time_diff.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM_all_x_punchthrough"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM_all_x_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM_all_x_punchthrough_cut%d_%s;Hodoscope time (ns);Counts", i_cut,
                 plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_diff_GEM_all_x_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_hodo_time_diff_GEM_all_x_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_time_diff_GEM_all_x_punchthrough_cut%d_%s;Hodoscope time Difference (ns);Counts",
                          i_cut, plane_names[i_plane].c_str()),
                     time_diff.nbins, time_diff.min, time_diff.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM_all_y_punchthrough"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM_all_y_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM_all_y_punchthrough_cut%d_%s;Hodoscope time (ns);Counts", i_cut,
                 plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_diff_GEM_all_y_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_hodo_time_diff_GEM_all_y_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_time_diff_GEM_all_y_punchthrough_cut%d_%s;Hodoscope time Difference (ns);Counts",
                          i_cut, plane_names[i_plane].c_str()),
                     time_diff.nbins, time_diff.min, time_diff.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_GEM_all_xy_punchthrough"][i_cut][i_plane] = new TH1F(
            Form("h_hodo_time_GEM_all_xy_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_time_GEM_all_xy_punchthrough_cut%d_%s;Hodoscope time (ns);Counts", i_cut,
                 plane_names[i_plane].c_str()),
            lad_time.nbins, lad_time.min, lad_time.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_time_diff_GEM_all_xy_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_hodo_time_diff_GEM_all_xy_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_time_diff_GEM_all_xy_punchthrough_cut%d_%s;Hodoscope time Difference (ns);Counts",
                          i_cut, plane_names[i_plane].c_str()),
                     time_diff.nbins, time_diff.min, time_diff.max);

        // ADC amplitude histograms
        hist_map_dxdy_cuts_vec[thread]["h_hodo_adc_amp"][i_cut][i_plane] =
            new TH1D(Form("h_hodo_adc_amp_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_adc_amp_cut%d_%s;ADC Amplitude;Counts", i_cut, plane_names[i_plane].c_str()),
                     LAD_ADC.nbins, LAD_ADC.min, LAD_ADC.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_adc_amp_punchthrough"][i_cut][i_plane] = new TH1D(
            Form("h_hodo_adc_amp_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_adc_amp_punchthrough_cut%d_%s;ADC Amplitude;Counts", i_cut, plane_names[i_plane].c_str()),
            LAD_ADC.nbins, LAD_ADC.min, LAD_ADC.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_adc_amp_GEM0_x_punchthrough"][i_cut][i_plane] = new TH1D(
            Form("h_hodo_adc_amp_GEM0_x_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_adc_amp_GEM0_x_punchthrough_cut%d_%s;ADC Amplitude;Counts", i_cut,
                 plane_names[i_plane].c_str()),
            LAD_ADC.nbins, LAD_ADC.min, LAD_ADC.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_adc_amp_GEM0_y_punchthrough"][i_cut][i_plane] = new TH1D(
            Form("h_hodo_adc_amp_GEM0_y_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_adc_amp_GEM0_y_punchthrough_cut%d_%s;ADC Amplitude;Counts", i_cut,
                 plane_names[i_plane].c_str()),
            LAD_ADC.nbins, LAD_ADC.min, LAD_ADC.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_adc_amp_GEM1_x_punchthrough"][i_cut][i_plane] = new TH1D(
            Form("h_hodo_adc_amp_GEM1_x_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_adc_amp_GEM1_x_punchthrough_cut%d_%s;ADC Amplitude;Counts", i_cut,
                 plane_names[i_plane].c_str()),
            LAD_ADC.nbins, LAD_ADC.min, LAD_ADC.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_adc_amp_GEM1_y_punchthrough"][i_cut][i_plane] = new TH1D(
            Form("h_hodo_adc_amp_GEM1_y_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_adc_amp_GEM1_y_punchthrough_cut%d_%s;ADC Amplitude;Counts", i_cut,
                 plane_names[i_plane].c_str()),
            LAD_ADC.nbins, LAD_ADC.min, LAD_ADC.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_adc_amp_GEM0_xy_punchthrough"][i_cut][i_plane] = new TH1D(
            Form("h_hodo_adc_amp_GEM0_xy_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_adc_amp_GEM0_xy_punchthrough_cut%d_%s;ADC Amplitude;Counts", i_cut,
                 plane_names[i_plane].c_str()),
            LAD_ADC.nbins, LAD_ADC.min, LAD_ADC.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_adc_amp_GEM1_xy_punchthrough"][i_cut][i_plane] = new TH1D(
            Form("h_hodo_adc_amp_GEM1_xy_punchthrough_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_hodo_adc_amp_GEM1_xy_punchthrough_cut%d_%s;ADC Amplitude;Counts", i_cut,
                 plane_names[i_plane].c_str()),
            LAD_ADC.nbins, LAD_ADC.min, LAD_ADC.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_adc_amp_GEM_all_x_punchthrough"][i_cut][i_plane] =
            new TH1D(Form("h_hodo_adc_amp_GEM_all_x_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_adc_amp_GEM_all_x_punchthrough_cut%d_%s;ADC Amplitude;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     LAD_ADC.nbins, LAD_ADC.min, LAD_ADC.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_adc_amp_GEM_all_y_punchthrough"][i_cut][i_plane] =
            new TH1D(Form("h_hodo_adc_amp_GEM_all_y_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_adc_amp_GEM_all_y_punchthrough_cut%d_%s;ADC Amplitude;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     LAD_ADC.nbins, LAD_ADC.min, LAD_ADC.max);

        hist_map_dxdy_cuts_vec[thread]["h_hodo_adc_amp_GEM_all_xy_punchthrough"][i_cut][i_plane] =
            new TH1D(Form("h_hodo_adc_amp_GEM_all_xy_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_hodo_adc_amp_GEM_all_xy_punchthrough_cut%d_%s;ADC Amplitude;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     LAD_ADC.nbins, LAD_ADC.min, LAD_ADC.max);

        // GEM Efficiency histograms (1D for x and y, for both GEM0 and GEM1)
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM0_x"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_all_GEM0_x_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_gem_efficiency_all_GEM0_x_cut%d_%s;Cluster x;Counts", i_cut, plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM0_x"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_good_hit_GEM0_x_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_gem_efficiency_good_hit_GEM0_x_cut%d_%s;Cluster x;Counts", i_cut, plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM0_y"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_all_GEM0_y_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_gem_efficiency_all_GEM0_y_cut%d_%s;Cluster y;Counts", i_cut, plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM0_y"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_good_hit_GEM0_y_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_gem_efficiency_good_hit_GEM0_y_cut%d_%s;Cluster y;Counts", i_cut, plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);

        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM1_x"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_all_GEM1_x_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_gem_efficiency_all_GEM1_x_cut%d_%s;Cluster x;Counts", i_cut, plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM1_x"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_good_hit_GEM1_x_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_gem_efficiency_good_hit_GEM1_x_cut%d_%s;Cluster x;Counts", i_cut, plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM1_y"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_all_GEM1_y_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_gem_efficiency_all_GEM1_y_cut%d_%s;Cluster y;Counts", i_cut, plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM1_y"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_good_hit_GEM1_y_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_gem_efficiency_good_hit_GEM1_y_cut%d_%s;Cluster y;Counts", i_cut, plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM0_2D"][i_cut][i_plane] = new TH2F(
            Form("h_gem_efficiency_all_GEM0_2D_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("GEM0 Efficiency All (x vs y);Golden Track x (cm);Golden Track y (cm)"), dx_NBINS, dx_min, dx_max,
            dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM0_2D"][i_cut][i_plane] = new TH2F(
            Form("h_gem_efficiency_good_hit_GEM0_2D_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("GEM0 Efficiency Good Hit (x vs y);Golden Track x (cm);Golden Track y (cm)"), dx_NBINS, dx_min, dx_max,
            dy_NBINS, dy_min, dy_max);

        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM1_2D"][i_cut][i_plane] = new TH2F(
            Form("h_gem_efficiency_all_GEM1_2D_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("GEM1 Efficiency All (x vs y);Golden Track x (cm);Golden Track y (cm)"), dx_NBINS, dx_min, dx_max,
            dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM1_2D"][i_cut][i_plane] = new TH2F(
            Form("h_gem_efficiency_good_hit_GEM1_2D_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("GEM1 Efficiency Good Hit (x vs y);Golden Track x (cm);Golden Track y (cm)"), dx_NBINS, dx_min, dx_max,
            dy_NBINS, dy_min, dy_max);

        // _bkd versions
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM0_x_bkd"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_all_GEM0_x_bkd_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_gem_efficiency_all_GEM0_x_bkd_cut%d_%s;Cluster x;Counts", i_cut, plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM0_x_bkd"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_good_hit_GEM0_x_bkd_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(),
                 thread),
            Form("h_gem_efficiency_good_hit_GEM0_x_bkd_cut%d_%s;Cluster x;Counts", i_cut, plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM0_y_bkd"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_all_GEM0_y_bkd_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_gem_efficiency_all_GEM0_y_bkd_cut%d_%s;Cluster y;Counts", i_cut, plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM0_y_bkd"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_good_hit_GEM0_y_bkd_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(),
                 thread),
            Form("h_gem_efficiency_good_hit_GEM0_y_bkd_cut%d_%s;Cluster y;Counts", i_cut, plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);

        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM1_x_bkd"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_all_GEM1_x_bkd_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_gem_efficiency_all_GEM1_x_bkd_cut%d_%s;Cluster x;Counts", i_cut, plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM1_x_bkd"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_good_hit_GEM1_x_bkd_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(),
                 thread),
            Form("h_gem_efficiency_good_hit_GEM1_x_bkd_cut%d_%s;Cluster x;Counts", i_cut, plane_names[i_plane].c_str()),
            dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM1_y_bkd"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_all_GEM1_y_bkd_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("h_gem_efficiency_all_GEM1_y_bkd_cut%d_%s;Cluster y;Counts", i_cut, plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM1_y_bkd"][i_cut][i_plane] = new TH1F(
            Form("h_gem_efficiency_good_hit_GEM1_y_bkd_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(),
                 thread),
            Form("h_gem_efficiency_good_hit_GEM1_y_bkd_cut%d_%s;Cluster y;Counts", i_cut, plane_names[i_plane].c_str()),
            dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM0_2D_bkd"][i_cut][i_plane] = new TH2F(
            Form("h_gem_efficiency_all_GEM0_2D_bkd_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("GEM0 Efficiency All BKD (x vs y);Golden Track x (cm);Golden Track y (cm)"), dx_NBINS, dx_min, dx_max,
            dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM0_2D_bkd"][i_cut][i_plane] =
            new TH2F(Form("h_gem_efficiency_good_hit_GEM0_2D_bkd_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("GEM0 Efficiency Good Hit BKD (x vs y);Golden Track x (cm);Golden Track y (cm)"), dx_NBINS,
                     dx_min, dx_max, dy_NBINS, dy_min, dy_max);

        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM1_2D_bkd"][i_cut][i_plane] = new TH2F(
            Form("h_gem_efficiency_all_GEM1_2D_bkd_cut%d_%s_thread_%d", i_cut, plane_names[i_plane].c_str(), thread),
            Form("GEM1 Efficiency All BKD (x vs y);Golden Track x (cm);Golden Track y (cm)"), dx_NBINS, dx_min, dx_max,
            dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM1_2D_bkd"][i_cut][i_plane] =
            new TH2F(Form("h_gem_efficiency_good_hit_GEM1_2D_bkd_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("GEM1 Efficiency Good Hit BKD (x vs y);Golden Track x (cm);Golden Track y (cm)"), dx_NBINS,
                     dx_min, dx_max, dy_NBINS, dy_min, dy_max);

        // Punchthrough versions
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM0_x_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_all_GEM0_x_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_all_GEM0_x_punchthrough_cut%d_%s;Cluster x;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM0_x_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_good_hit_GEM0_x_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_good_hit_GEM0_x_punchthrough_cut%d_%s;Cluster x;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM0_y_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_all_GEM0_y_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_all_GEM0_y_punchthrough_cut%d_%s;Cluster y;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM0_y_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_good_hit_GEM0_y_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_good_hit_GEM0_y_punchthrough_cut%d_%s;Cluster y;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);

        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM1_x_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_all_GEM1_x_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_all_GEM1_x_punchthrough_cut%d_%s;Cluster x;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM1_x_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_good_hit_GEM1_x_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_good_hit_GEM1_x_punchthrough_cut%d_%s;Cluster x;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM1_y_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_all_GEM1_y_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_all_GEM1_y_punchthrough_cut%d_%s;Cluster y;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM1_y_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_good_hit_GEM1_y_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_good_hit_GEM1_y_punchthrough_cut%d_%s;Cluster y;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM0_2D_punchthrough"][i_cut][i_plane] =
            new TH2F(Form("h_gem_efficiency_all_GEM0_2D_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("GEM0 Efficiency All Punchthrough (x vs y);Golden Track x (cm);Golden Track y (cm)"),
                     dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM0_2D_punchthrough"][i_cut][i_plane] =
            new TH2F(Form("h_gem_efficiency_good_hit_GEM0_2D_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("GEM0 Efficiency Good Hit Punchthrough (x vs y);Golden Track x (cm);Golden Track y (cm)"),
                     dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM1_2D_punchthrough"][i_cut][i_plane] =
            new TH2F(Form("h_gem_efficiency_all_GEM1_2D_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("GEM1 Efficiency All Punchthrough (x vs y);Golden Track x (cm);Golden Track y (cm)"),
                     dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM1_2D_punchthrough"][i_cut][i_plane] =
            new TH2F(Form("h_gem_efficiency_good_hit_GEM1_2D_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("GEM1 Efficiency Good Hit Punchthrough (x vs y);Golden Track x (cm);Golden Track y (cm)"),
                     dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);

        // _bkd punchthrough versions
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM0_x_bkd_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_all_GEM0_x_bkd_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_all_GEM0_x_bkd_punchthrough_cut%d_%s;Cluster x;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM0_x_bkd_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_good_hit_GEM0_x_bkd_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_good_hit_GEM0_x_bkd_punchthrough_cut%d_%s;Cluster x;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM0_y_bkd_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_all_GEM0_y_bkd_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_all_GEM0_y_bkd_punchthrough_cut%d_%s;Cluster y;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM0_y_bkd_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_good_hit_GEM0_y_bkd_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_good_hit_GEM0_y_bkd_punchthrough_cut%d_%s;Cluster y;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);

        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM1_x_bkd_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_all_GEM1_x_bkd_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_all_GEM1_x_bkd_punchthrough_cut%d_%s;Cluster x;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM1_x_bkd_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_good_hit_GEM1_x_bkd_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_good_hit_GEM1_x_bkd_punchthrough_cut%d_%s;Cluster x;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dx_NBINS, dx_min, dx_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM1_y_bkd_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_all_GEM1_y_bkd_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_all_GEM1_y_bkd_punchthrough_cut%d_%s;Cluster y;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM1_y_bkd_punchthrough"][i_cut][i_plane] =
            new TH1F(Form("h_gem_efficiency_good_hit_GEM1_y_bkd_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("h_gem_efficiency_good_hit_GEM1_y_bkd_punchthrough_cut%d_%s;Cluster y;Counts", i_cut,
                          plane_names[i_plane].c_str()),
                     dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM0_2D_bkd_punchthrough"][i_cut][i_plane] =
            new TH2F(Form("h_gem_efficiency_all_GEM0_2D_bkd_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("GEM0 Efficiency All BKD Punchthrough (x vs y);Golden Track x (cm);Golden Track y (cm)"),
                     dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM0_2D_bkd_punchthrough"][i_cut][i_plane] =
            new TH2F(Form("h_gem_efficiency_good_hit_GEM0_2D_bkd_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("GEM0 Efficiency Good Hit BKD Punchthrough (x vs y);Golden Track x (cm);Golden Track y (cm)"),
                     dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_all_GEM1_2D_bkd_punchthrough"][i_cut][i_plane] =
            new TH2F(Form("h_gem_efficiency_all_GEM1_2D_bkd_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("GEM1 Efficiency All BKD Punchthrough (x vs y);Golden Track x (cm);Golden Track y (cm)"),
                     dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);
        hist_map_dxdy_cuts_vec[thread]["h_gem_efficiency_good_hit_GEM1_2D_bkd_punchthrough"][i_cut][i_plane] =
            new TH2F(Form("h_gem_efficiency_good_hit_GEM1_2D_bkd_punchthrough_cut%d_%s_thread_%d", i_cut,
                          plane_names[i_plane].c_str(), thread),
                     Form("GEM1 Efficiency Good Hit BKD Punchthrough (x vs y);Golden Track x (cm);Golden Track y (cm)"),
                     dx_NBINS, dx_min, dx_max, dy_NBINS, dy_min, dy_max);
      }
    }
  }
  ///////////////////////////////////////////////////////////////////

  // start threads
  cout << "Starting " << numThreads << " threads..." << endl;
  std::vector<std::thread> threads;
  for (int i_thread = 0; i_thread < numThreads; ++i_thread) {
    int start = i_thread * chunkSize;
    int end   = (i_thread == numThreads - 1) ? nEntries : start + chunkSize;
    threads.emplace_back(process_chunk, i_thread, start, end, ref(fileNames), ref(hist_map_vec[i_thread]),
                         ref(hist_map_dxdy_cuts_vec[i_thread]));
  }
  // Wait for all threads to finish
  for (auto &thread : threads) {
    thread.join();
  }
  std::cout << "\rProcessing: 100% completed. \nMerging histograms" << std::endl;

  // Merge histograms from all threads by adding to the first thread
  for (int i_thread = 0; i_thread < numThreads; ++i_thread) {
    for (const auto &pair : hist_map_vec[i_thread]) {
      const std::string &histName = pair.first;
      for (int i_gem = 0; i_gem < 2; ++i_gem) {
        for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
          if (hist_map_vec[0].count(histName) && hist_map_vec[0][histName][i_gem][i_plane]) {
            if (i_thread == 0) {
              // Set the title and name for the histogram in the zero thread
              // Remove "thread_%d" from the title and name
              TString title = hist_map_vec[i_thread][histName][i_gem][i_plane]->GetTitle();
              TString name  = hist_map_vec[i_thread][histName][i_gem][i_plane]->GetName();
              title.ReplaceAll(Form("_thread_%d", i_thread), "");
              name.ReplaceAll(Form("_thread_%d", i_thread), "");
              hist_map_vec[0][histName][i_gem][i_plane]->SetTitle(title);
              hist_map_vec[0][histName][i_gem][i_plane]->SetName(name);
            } else {
              // Add the histogram to the zero thread
              hist_map_vec[0][histName][i_gem][i_plane]->Add(hist_map_vec[i_thread][histName][i_gem][i_plane]);
            }
          }
        }
      }
    }
    for (const auto &pair : hist_map_dxdy_cuts_vec[i_thread]) {
      const std::string &histName = pair.first;
      for (int i_cut = 0; i_cut < n_dxy_cuts; ++i_cut) {
        for (int i_plane = 0; i_plane < nPlanes; ++i_plane) {
          if (hist_map_dxdy_cuts_vec[0].count(histName) && hist_map_dxdy_cuts_vec[0][histName][i_cut][i_plane]) {
            if (i_thread == 0) {
              // Set the title and name for the histogram in the zero thread
              // Remove "thread_%d" from the title and name
              TString title = hist_map_dxdy_cuts_vec[i_thread][histName][i_cut][i_plane]->GetTitle();
              TString name  = hist_map_dxdy_cuts_vec[i_thread][histName][i_cut][i_plane]->GetName();
              title.ReplaceAll(Form("_thread_%d", i_thread), "");
              name.ReplaceAll(Form("_thread_%d", i_thread), "");
              hist_map_dxdy_cuts_vec[0][histName][i_cut][i_plane]->SetTitle(title);
              hist_map_dxdy_cuts_vec[0][histName][i_cut][i_plane]->SetName(name);
            } else {
              // Add the histogram to the zero thread
              hist_map_dxdy_cuts_vec[0][histName][i_cut][i_plane]->Add(
                  hist_map_dxdy_cuts_vec[i_thread][histName][i_cut][i_plane]);
            }
          }
        }
      }
    }
  }

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
    for (int i_cut = 0; i_cut < n_dxy_cuts; ++i_cut) {
      outputFile->mkdir(Form("time/%s/cut%d/Times", plane_names[i_plane].c_str(), i_cut));
      outputFile->mkdir(Form("time/%s/cut%d/Diff_Times", plane_names[i_plane].c_str(), i_cut));
      outputFile->mkdir(Form("time/%s/cut%d/ADC_Amp", plane_names[i_plane].c_str(), i_cut));
      outputFile->mkdir(Form("GEM_Efficiency_BkdSub/%s/cut%d", plane_names[i_plane].c_str(), i_cut));
      outputFile->mkdir(Form("GEM_Efficiency/%s/cut%d", plane_names[i_plane].c_str(), i_cut));
    }
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
      hist_map_vec[0]["h_gem_cust_dx"][i_gem][i_plane]->Write();
      outputFile->cd("clust/x/2D");
      hist_map_vec[0]["h_gem_clust_dx_x"][i_gem][i_plane]->Write();
      outputFile->cd("clust/x/2D_ADC");
      hist_map_vec[0]["h_gem_clust_dx_adc"][i_gem][i_plane]->SetOption("COLZ");
      hist_map_vec[0]["h_gem_clust_dx_adc"][i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy");
      hist_map_vec[0]["h_gem_cust_dy"][i_gem][i_plane]->Write();
      outputFile->cd("clust/y/2D");
      hist_map_vec[0]["h_gem_clust_dy_y"][i_gem][i_plane]->Write();
      outputFile->cd("clust/y/2D_ADC");
      hist_map_vec[0]["h_gem_clust_dy_adc"][i_gem][i_plane]->SetOption("COLZ");
      hist_map_vec[0]["h_gem_clust_dy_adc"][i_gem][i_plane]->Write();
      outputFile->cd("clust/x/dx_punchthrough");
      hist_map_vec[0]["h_gem_cust_dx_punchthrough"][i_gem][i_plane]->Write();
      outputFile->cd("clust/x/2D_punchthrough");
      hist_map_vec[0]["h_gem_clust_dx_x_punchthrough"][i_gem][i_plane]->Write();
      outputFile->cd("clust/x/2D_ADC_punchthrough");
      hist_map_vec[0]["h_gem_clust_dx_adc_punchthrough"][i_gem][i_plane]->SetOption("COLZ");
      hist_map_vec[0]["h_gem_clust_dx_adc_punchthrough"][i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy_punchthrough");
      hist_map_vec[0]["h_gem_cust_dy_punchthrough"][i_gem][i_plane]->Write();
      outputFile->cd("clust/y/2D_punchthrough");
      hist_map_vec[0]["h_gem_clust_dy_y_punchthrough"][i_gem][i_plane]->Write();
      outputFile->cd("clust/y/2D_ADC_punchthrough");
      hist_map_vec[0]["h_gem_clust_dy_adc_punchthrough"][i_gem][i_plane]->SetOption("COLZ");
      hist_map_vec[0]["h_gem_clust_dy_adc_punchthrough"][i_gem][i_plane]->Write();
      outputFile->cd("clust/dxdy/all");
      hist_map_vec[0]["h_gem_clust_dxy"][i_gem][i_plane]->Write();
      outputFile->cd("clust/dxdy/punchthrough");
      hist_map_vec[0]["h_gem_clust_dxy_punchthrough"][i_gem][i_plane]->Write();

      outputFile->cd("clust/x/dx/min/all");
      hist_map_vec[0]["h_gem_cust_dx_min"][i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/all");
      hist_map_vec[0]["h_gem_cust_dy_min"][i_gem][i_plane]->Write();
      outputFile->cd("clust/x/dx/min/before");
      hist_map_vec[0]["h_gem_cust_dx_min_before"][i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/before");
      hist_map_vec[0]["h_gem_cust_dy_min_before"][i_gem][i_plane]->Write();
      outputFile->cd("clust/x/dx/min/peak");
      hist_map_vec[0]["h_gem_cust_dx_min_peak"][i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/peak");
      hist_map_vec[0]["h_gem_cust_dy_min_peak"][i_gem][i_plane]->Write();
      outputFile->cd("clust/x/dx/min/after");
      hist_map_vec[0]["h_gem_cust_dx_min_after"][i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/after");
      hist_map_vec[0]["h_gem_cust_dy_min_after"][i_gem][i_plane]->Write();

      outputFile->cd("clust/x/dx/min/before_punchthrough");
      hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_before"][i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/before_punchthrough");
      hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_before"][i_gem][i_plane]->Write();
      outputFile->cd("clust/x/dx/min/peak_punchthrough");
      hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_peak"][i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/peak_punchthrough");
      hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_peak"][i_gem][i_plane]->Write();
      outputFile->cd("clust/x/dx/min/after_punchthrough");
      hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_after"][i_gem][i_plane]->Write();
      outputFile->cd("clust/y/dy/min/after_punchthrough");
      hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_after"][i_gem][i_plane]->Write();
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
      hist_map_vec[0]["h_gem_cust_dx_min_before"][i_gem][i_plane]->SetLineColor(kBlue);
      hist_map_vec[0]["h_gem_cust_dx_min_peak"][i_gem][i_plane]->SetLineColor(kRed);
      hist_map_vec[0]["h_gem_cust_dx_min_after"][i_gem][i_plane]->SetLineColor(kBlack);

      hist_map_vec[0]["h_gem_cust_dx_min_before"][i_gem][i_plane]->Scale(
          1.0 / hist_map_vec[0]["h_gem_cust_dx_min_before"][i_gem][i_plane]->Integral());
      hist_map_vec[0]["h_gem_cust_dx_min_peak"][i_gem][i_plane]->Scale(
          1.0 / hist_map_vec[0]["h_gem_cust_dx_min_peak"][i_gem][i_plane]->Integral());
      hist_map_vec[0]["h_gem_cust_dx_min_after"][i_gem][i_plane]->Scale(
          1.0 / hist_map_vec[0]["h_gem_cust_dx_min_after"][i_gem][i_plane]->Integral());
      hist_map_vec[0]["h_gem_cust_dx_min_before"][i_gem][i_plane]->Draw("hist");
      hist_map_vec[0]["h_gem_cust_dx_min_peak"][i_gem][i_plane]->Draw("hist same");
      hist_map_vec[0]["h_gem_cust_dx_min_after"][i_gem][i_plane]->Draw("hist same");
      legend = new TLegend(0.65, 0.7, 0.88, 0.88);
      legend->AddEntry(hist_map_vec[0]["h_gem_cust_dx_min_before"][i_gem][i_plane], "Before", "l");
      legend->AddEntry(hist_map_vec[0]["h_gem_cust_dx_min_peak"][i_gem][i_plane], "Peak", "l");
      legend->AddEntry(hist_map_vec[0]["h_gem_cust_dx_min_after"][i_gem][i_plane], "After", "l");
      legend->Draw();
      c1->Update();
      c1->Write();
      outputFile->cd("clust/y/dy/min/comp");
      c1 = new TCanvas(Form("c1_y_%d_%s", i_gem, plane_names[i_plane].c_str()), "Comparison Canvas", 800, 600);
      c1->cd();
      hist_map_vec[0]["h_gem_cust_dy_min_before"][i_gem][i_plane]->SetLineColor(kBlue);
      hist_map_vec[0]["h_gem_cust_dy_min_peak"][i_gem][i_plane]->SetLineColor(kRed);
      hist_map_vec[0]["h_gem_cust_dy_min_after"][i_gem][i_plane]->SetLineColor(kBlack);
      hist_map_vec[0]["h_gem_cust_dy_min_before"][i_gem][i_plane]->Scale(
          1.0 / hist_map_vec[0]["h_gem_cust_dy_min_before"][i_gem][i_plane]->Integral());
      hist_map_vec[0]["h_gem_cust_dy_min_peak"][i_gem][i_plane]->Scale(
          1.0 / hist_map_vec[0]["h_gem_cust_dy_min_peak"][i_gem][i_plane]->Integral());
      hist_map_vec[0]["h_gem_cust_dy_min_after"][i_gem][i_plane]->Scale(
          1.0 / hist_map_vec[0]["h_gem_cust_dy_min_after"][i_gem][i_plane]->Integral());
      hist_map_vec[0]["h_gem_cust_dy_min_before"][i_gem][i_plane]->Draw("hist");
      hist_map_vec[0]["h_gem_cust_dy_min_peak"][i_gem][i_plane]->Draw("hist same");
      hist_map_vec[0]["h_gem_cust_dy_min_after"][i_gem][i_plane]->Draw("hist same");
      legend = new TLegend(0.65, 0.7, 0.88, 0.88);
      legend->AddEntry(hist_map_vec[0]["h_gem_cust_dy_min_before"][i_gem][i_plane], "Before", "l");
      legend->AddEntry(hist_map_vec[0]["h_gem_cust_dy_min_peak"][i_gem][i_plane], "Peak", "l");
      legend->AddEntry(hist_map_vec[0]["h_gem_cust_dy_min_after"][i_gem][i_plane], "After", "l");
      legend->Draw();
      c1->Update();
      c1->Write();

      outputFile->cd("clust/x/dx/min/comp_punchthrough");
      c1 = new TCanvas(Form("c1_x_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()), "Comparison Canvas", 800,
                       600);
      c1->cd();
      hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_before"][i_gem][i_plane]->SetLineColor(kBlue);
      hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_peak"][i_gem][i_plane]->SetLineColor(kRed);
      hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_after"][i_gem][i_plane]->SetLineColor(kBlack);
      hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_before"][i_gem][i_plane]->Scale(
          1.0 / hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_before"][i_gem][i_plane]->Integral());
      hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_peak"][i_gem][i_plane]->Scale(
          1.0 / hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_peak"][i_gem][i_plane]->Integral());
      hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_after"][i_gem][i_plane]->Scale(
          1.0 / hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_after"][i_gem][i_plane]->Integral());
      hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_before"][i_gem][i_plane]->Draw("hist");
      hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_peak"][i_gem][i_plane]->Draw("hist same");
      hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_after"][i_gem][i_plane]->Draw("hist same");
      legend = new TLegend(0.65, 0.7, 0.88, 0.88);
      legend->AddEntry(hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_before"][i_gem][i_plane], "Before", "l");
      legend->AddEntry(hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_peak"][i_gem][i_plane], "Peak", "l");
      legend->AddEntry(hist_map_vec[0]["h_gem_cust_dx_punchthrough_min_after"][i_gem][i_plane], "After", "l");
      legend->Draw();
      c1->Update();
      c1->Write();
      outputFile->cd("clust/y/dy/min/comp_punchthrough");
      c1 = new TCanvas(Form("c1_y_punchthrough_%d_%s", i_gem, plane_names[i_plane].c_str()), "Comparison Canvas", 800,
                       600);
      c1->cd();
      hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_before"][i_gem][i_plane]->SetLineColor(kBlue);
      hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_peak"][i_gem][i_plane]->SetLineColor(kRed);
      hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_after"][i_gem][i_plane]->SetLineColor(kBlack);
      hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_before"][i_gem][i_plane]->Scale(
          1.0 / hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_before"][i_gem][i_plane]->Integral());
      hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_peak"][i_gem][i_plane]->Scale(
          1.0 / hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_peak"][i_gem][i_plane]->Integral());
      hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_after"][i_gem][i_plane]->Scale(
          1.0 / hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_after"][i_gem][i_plane]->Integral());
      hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_before"][i_gem][i_plane]->Draw("hist");
      hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_peak"][i_gem][i_plane]->Draw("hist same");
      hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_after"][i_gem][i_plane]->Draw("hist same");
      legend = new TLegend(0.65, 0.7, 0.88, 0.88);
      legend->AddEntry(hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_before"][i_gem][i_plane], "Before", "l");
      legend->AddEntry(hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_peak"][i_gem][i_plane], "Peak", "l");
      legend->AddEntry(hist_map_vec[0]["h_gem_cust_dy_punchthrough_min_after"][i_gem][i_plane], "After", "l");
      legend->Draw();
      c1->Update();
      c1->Write();

      // Space point histograms
      outputFile->cd("sp/x/dx");
      hist_map_vec[0]["h_gem_sp_dx"][i_gem][i_plane]->Write();
      outputFile->cd("sp/x/2D");
      hist_map_vec[0]["h_gem_sp_dx_x"][i_gem][i_plane]->Write();
      outputFile->cd("sp/x/2D_ADC");
      hist_map_vec[0]["h_gem_sp_dx_adc"][i_gem][i_plane]->SetOption("COLZ");
      hist_map_vec[0]["h_gem_sp_dx_adc"][i_gem][i_plane]->Write();
      outputFile->cd("sp/y/dy");
      hist_map_vec[0]["h_gem_sp_dy"][i_gem][i_plane]->Write();
      outputFile->cd("sp/y/2D");
      hist_map_vec[0]["h_gem_sp_dy_y"][i_gem][i_plane]->Write();
      outputFile->cd("sp/y/2D_ADC");
      hist_map_vec[0]["h_gem_sp_dy_adc"][i_gem][i_plane]->SetOption("COLZ");
      hist_map_vec[0]["h_gem_sp_dy_adc"][i_gem][i_plane]->Write();
      outputFile->cd("sp/x/dx_punchthrough");
      hist_map_vec[0]["h_gem_sp_dx_punchthrough"][i_gem][i_plane]->Write();
      outputFile->cd("sp/x/2D_punchthrough");
      hist_map_vec[0]["h_gem_sp_dx_x_punchthrough"][i_gem][i_plane]->Write();
      outputFile->cd("sp/x/2D_ADC_punchthrough");
      hist_map_vec[0]["h_gem_sp_dx_adc_punchthrough"][i_gem][i_plane]->SetOption("COLZ");
      hist_map_vec[0]["h_gem_sp_dx_adc_punchthrough"][i_gem][i_plane]->Write();
      outputFile->cd("sp/y/dy_punchthrough");
      hist_map_vec[0]["h_gem_sp_dy_punchthrough"][i_gem][i_plane]->Write();
      outputFile->cd("sp/y/2D_punchthrough");
      hist_map_vec[0]["h_gem_sp_dy_y_punchthrough"][i_gem][i_plane]->Write();
      outputFile->cd("sp/y/2D_ADC_punchthrough");
      hist_map_vec[0]["h_gem_sp_dy_adc_punchthrough"][i_gem][i_plane]->SetOption("COLZ");
      hist_map_vec[0]["h_gem_sp_dy_adc_punchthrough"][i_gem][i_plane]->Write();
      outputFile->cd("sp/dxdy/all");
      hist_map_vec[0]["h_gem_sp_dxy"][i_gem][i_plane]->Write();
      outputFile->cd("sp/dxdy/punchthrough");
      hist_map_vec[0]["h_gem_sp_dxy_punchthrough"][i_gem][i_plane]->Write();
    }
    // Write time histograms
    for (int i_cut = 0; i_cut < n_dxy_cuts; ++i_cut) {
      outputFile->cd(Form("time/%s/cut%d/Times", plane_names[i_plane].c_str(), i_cut));
      hist_map_dxdy_cuts_vec[0]["h_hodo_time"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM0_x"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM0_y"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM1_x"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM1_y"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM0_x_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM0_y_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM1_x_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM1_y_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM0_xy"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM1_xy"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM_all_x"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM_all_y"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM_all_xy"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM0_xy_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM1_xy_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM_all_x_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM_all_y_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_GEM_all_xy_punchthrough"][i_cut][i_plane]->Write();
      outputFile->cd(Form("time/%s/cut%d/Diff_Times", plane_names[i_plane].c_str(), i_cut));
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_diff_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_diff_GEM0_x_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_diff_GEM0_y_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_diff_GEM1_x_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_diff_GEM1_y_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_diff_GEM0_xy_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_diff_GEM1_xy_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_diff_GEM_all_x_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_diff_GEM_all_y_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_time_diff_GEM_all_xy_punchthrough"][i_cut][i_plane]->Write();
      outputFile->cd(Form("time/%s/cut%d/ADC_Amp", plane_names[i_plane].c_str(), i_cut));
      hist_map_dxdy_cuts_vec[0]["h_hodo_adc_amp"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_adc_amp_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_adc_amp_GEM0_x_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_adc_amp_GEM0_y_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_adc_amp_GEM1_x_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_adc_amp_GEM1_y_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_adc_amp_GEM0_xy_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_adc_amp_GEM1_xy_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_adc_amp_GEM_all_x_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_adc_amp_GEM_all_y_punchthrough"][i_cut][i_plane]->Write();
      hist_map_dxdy_cuts_vec[0]["h_hodo_adc_amp_GEM_all_xy_punchthrough"][i_cut][i_plane]->Write();
      outputFile->cd(Form("GEM_Efficiency_BkdSub/%s/cut%d", plane_names[i_plane].c_str(), i_cut));
      string bkdSubDir = Form("GEM_Efficiency_BkdSub/%s/cut%d", plane_names[i_plane].c_str(), i_cut);
      string noSubDir = Form("GEM_Efficiency/%s/cut%d", plane_names[i_plane].c_str(), i_cut);
      // Draw efficiency histograms: good_hit / all for GEM0_x, GEM0_y, GEM1_x, GEM1_y
      // Compute signal - background for each histogram, then divide by all
      double bkdSub_scale = (time_window_sig.max - time_window_sig.min) / (time_window_bkd.max - time_window_bkd.min);

      // GEM0_x
      drawEfficiencyHistograms(
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_x"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_x"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_x_bkd"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_x_bkd"][i_cut][i_plane], bkdSub_scale,
          Form("GEM0_x Efficiency cut%d_%s", i_cut, plane_names[i_plane].c_str()),
          "GEM0_x", outputFile,noSubDir.c_str(),bkdSubDir.c_str());

      // GEM0_y
      drawEfficiencyHistograms(
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_y"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_y"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_y_bkd"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_y_bkd"][i_cut][i_plane], bkdSub_scale,
          Form("GEM0_y Efficiency cut%d_%s", i_cut, plane_names[i_plane].c_str()),
          "GEM0_y", outputFile,noSubDir.c_str(),bkdSubDir.c_str());

      // GEM1_x
      drawEfficiencyHistograms(
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_x"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_x"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_x_bkd"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_x_bkd"][i_cut][i_plane], bkdSub_scale,
          Form("GEM1_x Efficiency cut%d_%s", i_cut, plane_names[i_plane].c_str()),
          "GEM1_x", outputFile,noSubDir.c_str(),bkdSubDir.c_str());

      // GEM1_y
      drawEfficiencyHistograms(
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_y"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_y"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_y_bkd"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_y_bkd"][i_cut][i_plane], bkdSub_scale,
          Form("GEM1_y Efficiency cut%d_%s", i_cut, plane_names[i_plane].c_str()),
          "GEM1_y", outputFile,noSubDir.c_str(),bkdSubDir.c_str());

      // GEM0_2D
      drawEfficiencyHistograms(
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_2D"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_2D"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_2D_bkd"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_2D_bkd"][i_cut][i_plane], bkdSub_scale,
          Form("GEM0_2D Efficiency cut%d_%s", i_cut, plane_names[i_plane].c_str()),
          "GEM0_2D", outputFile,noSubDir.c_str(),bkdSubDir.c_str());

      // GEM1_2D
      drawEfficiencyHistograms(
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_2D"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_2D"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_2D_bkd"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_2D_bkd"][i_cut][i_plane], bkdSub_scale,
          Form("GEM1_2D Efficiency cut%d_%s", i_cut, plane_names[i_plane].c_str()),
          "GEM1_2D", outputFile,noSubDir.c_str(),bkdSubDir.c_str());

      // Punchthrough versions
      // GEM0_x
      drawEfficiencyHistograms(
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_x_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_x_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_x_bkd_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_x_bkd_punchthrough"][i_cut][i_plane], bkdSub_scale,
          Form("GEM0_x Efficiency Punchthrough cut%d_%s", i_cut, plane_names[i_plane].c_str()),
          "GEM0_x_punchthrough", outputFile,noSubDir.c_str(),bkdSubDir.c_str());

      // GEM0_y
      drawEfficiencyHistograms(
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_y_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_y_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_y_bkd_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_y_bkd_punchthrough"][i_cut][i_plane], bkdSub_scale,
          Form("GEM0_y Efficiency Punchthrough cut%d_%s", i_cut, plane_names[i_plane].c_str()),
          "GEM0_y_punchthrough", outputFile,noSubDir.c_str(),bkdSubDir.c_str());

      // GEM1_x
      drawEfficiencyHistograms(
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_x_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_x_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_x_bkd_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_x_bkd_punchthrough"][i_cut][i_plane], bkdSub_scale,
          Form("GEM1_x Efficiency Punchthrough cut%d_%s", i_cut, plane_names[i_plane].c_str()),
          "GEM1_x_punchthrough", outputFile,noSubDir.c_str(),bkdSubDir.c_str());

      // GEM1_y
      drawEfficiencyHistograms(
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_y_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_y_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_y_bkd_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_y_bkd_punchthrough"][i_cut][i_plane], bkdSub_scale,
          Form("GEM1_y Efficiency Punchthrough cut%d_%s", i_cut, plane_names[i_plane].c_str()),
          "GEM1_y_punchthrough", outputFile,noSubDir.c_str(),bkdSubDir.c_str());

      // GEM0 2D punchthrough
      drawEfficiencyHistograms(
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_2D_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_2D_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_2D_bkd_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_2D_bkd_punchthrough"][i_cut][i_plane], bkdSub_scale,
          Form("GEM0 2D Efficiency Punchthrough cut%d_%s", i_cut, plane_names[i_plane].c_str()),
          "GEM0_2D_punchthrough", outputFile,noSubDir.c_str(),bkdSubDir.c_str());

      // GEM1 2D punchthrough
      drawEfficiencyHistograms(
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_2D_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_2D_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_2D_bkd_punchthrough"][i_cut][i_plane],
          hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_2D_bkd_punchthrough"][i_cut][i_plane], bkdSub_scale,
          Form("GEM1 2D Efficiency Punchthrough cut%d_%s", i_cut, plane_names[i_plane].c_str()),
          "GEM1_2D_punchthrough", outputFile,noSubDir.c_str(),bkdSubDir.c_str());

      // Write the efficiency histograms
      // Write GEM Efficiency histograms for signal region only (no background subtraction)
      // Efficiency = good_hit / all
      // outputFile->cd(Form("GEM_Efficiency/%s/cut%d", plane_names[i_plane].c_str(), i_cut));

      // TH1F *h_eff_GEM0_x = (TH1F *)hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_x"][i_cut][i_plane]->Clone(
      //     Form("h_gem_efficiency_ratio_GEM0_x_cut%d_%s", i_cut, plane_names[i_plane].c_str()));
      // h_eff_GEM0_x->SetTitle("GEM0_x Efficiency (good_hit/all, signal region only)");
      // h_eff_GEM0_x->Divide(hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_x"][i_cut][i_plane]);
      // h_eff_GEM0_x->Write();

      // TH1F *h_eff_GEM0_y = (TH1F *)hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_y"][i_cut][i_plane]->Clone(
      //     Form("h_gem_efficiency_ratio_GEM0_y_cut%d_%s", i_cut, plane_names[i_plane].c_str()));
      // h_eff_GEM0_y->SetTitle("GEM0_y Efficiency (good_hit/all, signal region only)");
      // h_eff_GEM0_y->Divide(hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_y"][i_cut][i_plane]);
      // h_eff_GEM0_y->Write();

      // TH1F *h_eff_GEM1_x = (TH1F *)hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_x"][i_cut][i_plane]->Clone(
      //     Form("h_gem_efficiency_ratio_GEM1_x_cut%d_%s", i_cut, plane_names[i_plane].c_str()));
      // h_eff_GEM1_x->SetTitle("GEM1_x Efficiency (good_hit/all, signal region only)");
      // h_eff_GEM1_x->Divide(hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_x"][i_cut][i_plane]);
      // h_eff_GEM1_x->Write();

      // TH1F *h_eff_GEM1_y = (TH1F *)hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_y"][i_cut][i_plane]->Clone(
      //     Form("h_gem_efficiency_ratio_GEM1_y_cut%d_%s", i_cut, plane_names[i_plane].c_str()));
      // h_eff_GEM1_y->SetTitle("GEM1_y Efficiency (good_hit/all, signal region only)");
      // h_eff_GEM1_y->Divide(hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_y"][i_cut][i_plane]);
      // h_eff_GEM1_y->Write();

      // TH2F *h_eff_GEM0_2D =
      //     (TH2F *)hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM0_2D"][i_cut][i_plane]->Clone(
      //         Form("h_gem_efficiency_ratio_GEM0_2D_cut%d_%s", i_cut, plane_names[i_plane].c_str()));
      // h_eff_GEM0_2D->SetTitle("GEM0_2D Efficiency (good_hit/all, signal region only)");
      // h_eff_GEM0_2D->Divide(hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM0_2D"][i_cut][i_plane]);
      // h_eff_GEM0_2D->Write();

      // TH2F *h_eff_GEM1_2D =
      //     (TH2F *)hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_good_hit_GEM1_2D"][i_cut][i_plane]->Clone(
      //         Form("h_gem_efficiency_ratio_GEM1_2D_cut%d_%s", i_cut, plane_names[i_plane].c_str()));
      // h_eff_GEM1_2D->SetTitle("GEM1_2D Efficiency (good_hit/all, signal region only)");
      // h_eff_GEM1_2D->Divide(hist_map_dxdy_cuts_vec[0]["h_gem_efficiency_all_GEM1_2D"][i_cut][i_plane]);
      // h_eff_GEM1_2D->Write();
    }
  }

  // Close the files
  file->Close();
  outputFile->Close();
  std::cout << "Output file created: " << outputFileName << std::endl;
  std::_Exit(0);
  return;
}