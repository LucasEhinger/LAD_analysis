#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

const int voltage_word_index = 7;

void modifyVoltage(const std::string &inputFilename, const std::string &outputFilename, int changeAmount) {
  std::ifstream inputFile(inputFilename);
  if (!inputFile.is_open()) {
    std::cerr << "Error opening file: " << inputFilename << std::endl;
    return;
  }

  std::vector<std::string> lines;
  std::string line;

  while (std::getline(inputFile, line)) {
    if (line.empty()) {
      lines.push_back(line);
      continue;
    }

    if (line[0] == '#') {
      lines.push_back(line);
      continue;
    }

    std::istringstream iss(line);
    std::string word;
    std::vector<std::string> words;
    std::vector<std::string> separators;

    while (iss >> word) {
      words.push_back(word);
      char separator;
      if (iss.get(separator)) {
        separators.push_back(std::string(1, separator));
      } else {
        separators.push_back("");
      }
    }

    if (!words.empty() &&
        (words[0].find("HVG") != std::string::npos || words[0].find("Ref-Bar") != std::string::npos) &&
        words.size() > voltage_word_index) {
      try {
        int voltage = std::stoi(words[voltage_word_index]);
        voltage += changeAmount;
        words[voltage_word_index] = std::to_string(voltage);
      } catch (const std::invalid_argument &e) {
        std::cerr << "Invalid voltage value on line: " << line << std::endl;
      }
    }

    std::ostringstream oss;
    for (size_t i = 0; i < words.size(); ++i) {
      oss << words[i];
      if (i < separators.size()) {
        oss << separators[i];
      }
    }

    lines.push_back(oss.str());
  }

  inputFile.close();

  std::ofstream outputFile(outputFilename);
  if (!outputFile.is_open()) {
    std::cerr << "Error opening file for writing: " << outputFilename << std::endl;
    return;
  }

  for (const auto &modifiedLine : lines) {
    outputFile << modifiedLine << std::endl;
  }

  outputFile.close();
}

int voltage_edit() {
  std::vector<int> changeAmounts;
  for (int i = -450; i <= 150; i += 30) {
    changeAmounts.push_back(i);
  }
  for (int changeAmount : changeAmounts) {
    std::string inputFilename  = "files/HV_Files/HV-backup_2025-03-26_orig.sav"; // Replace with your input file name
    std::string outputFilename = "files/HV_Files/HV-backup_2025-03-26_change" + std::to_string(changeAmount) + ".sav";

    modifyVoltage(inputFilename, outputFilename, changeAmount);
  }

  return 0;
}