#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

void change_gem_param(double m0_x, double m0_y, double m0_z, double m1_x, double m1_y, double m1_z) {
  std::ifstream infile("gem_param");
  if (!infile.is_open()) {
    std::cerr << "Error: Could not open file 'gem_param' for reading." << std::endl;
    return;
  }

  std::ofstream outfile("gem_param_modified");
  if (!outfile.is_open()) {
    std::cerr << "Error: Could not open file 'gem_param_modified' for writing." << std::endl;
    return;
  }

  std::string line;
  while (std::getline(infile, line)) {
    if (line.find("lgem_m0_position") != std::string::npos) {
      outfile << "lgem_m0_position = 0, " << m0_x << ", " << m0_y << ", " << m0_z << "; Position from survey" << std::endl;
    } else if (line.find("lgem_m1_position") != std::string::npos) {
      outfile << "lgem_m1_position = 1, " << m1_x << ", " << m1_y << ", " << m1_z << "; Position from survey" << std::endl;
    } else {
      outfile << line << std::endl;
    }
  }

  infile.close();
  outfile.close();

  // Optionally, replace the original file with the modified one
  if (std::rename("gem_param_modified", "gem_param") != 0) {
    std::cerr << "Error: Could not replace the original file with the modified one." << std::endl;
  }
}