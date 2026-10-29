#include "TmpDir.h"

#include <filesystem>
#include <iostream>
#include <random>
#include <sstream>

void TmpDir::createTempFile(std::string& tcf, std::string& tcd) {
  namespace fs = std::filesystem;
  try {
    // a new directory with a random name in the system's temporary directory
    std::random_device seed;
    std::mt19937_64 generator(seed());
    fs::path dir;
    do {
      std::ostringstream name;
      name << "percolator-" << std::hex << generator();
      dir = fs::temp_directory_path() / name.str();
    } while (!fs::create_directory(dir));
    tcd = dir.string();
    tcf = (dir / "converters-tmp.tcb").string();
  } catch (fs::filesystem_error& e) {
    std::cerr << e.what() << std::endl;
  }
}
