#include "mattak/Monitoring.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

// Include mattak classes
#include "mattak/Waveforms.h"
#include "mattak/Header.h"
#include "mattak/DAQStatus.h"
#include "mattak/Constants.h"
// #include "mattak/RunInfo.h"
// #include "mattak/Pedestals.h"
// #include "mattak/Converter.h" 

#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TBranch.h>
#include <TLeaf.h>


ClassImp(mattak::Monitoring); 
// non-member function to check if a file exists
// bool fileExists(const std::string& filename) {
//     std::ifstream file(filename);
//     bool fstatus = file.good();
//     if (fstatus) {
//         std::cout << "File " << filename << " exists." << std::endl;
//     } else {
//         std::cout << "File " << filename << " does not exist." << std::endl;
//     }
//     file.close();
//     return fstatus;
// }
mattak::Monitoring::Monitoring(uint32_t run, uint16_t station)
{
    this->run_number = run;
    this->station_number = station;
};



