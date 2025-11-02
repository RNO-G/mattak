#include "mattak/Monitoring.h"
#include <iostream>
// #include <format>
#include <string>
#include <sstream>
#include <cstdlib>

// Include mattak classes
#include "mattak/Waveforms.h"
#include "mattak/Converter.h"
#include "mattak/Header.h"
#include "mattak/DAQStatus.h"
#include "mattak/Pedestals.h"
#include "mattak/RunInfo.h"

#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TBranch.h>
#include <TLeaf.h>


ClassImp(mattak::Monitoring); 

mattak::Monitoring::Monitoring(uint32_t run, uint16_t station, std::string RNOG_DATA_DIR)
{
    this->run_number = run;
    this->station_number = station;
    
    std::stringstream  headers_ss; 
    headers_ss << RNOG_DATA_DIR << "/station" << station << "/run" << run << "/headers.root";
    std::string headers_filename = headers_ss.str();
    std::cout << "Header filename: " << headers_filename << std::endl;
    TFile* headers_file = TFile::Open(headers_filename.c_str(), "READ");
    if (!headers_file || headers_file->IsZombie()) 
    {
        std::cerr << "Error opening file: " << headers_filename << std::endl;
    }
    else 
    {
        std::cout << "Successfully opened file: " << headers_filename << std::endl;
        headers_file->Close();
    }
    std::cout<<"\n"<<std::endl;

    std::stringstream  daqstatus_ss; 
    daqstatus_ss << RNOG_DATA_DIR << "/station" << station << "/run" << run << "/daqstatus.root";
    std::string daqstatus_filename = daqstatus_ss.str();
    std::cout << "DAQStatus filename: " << daqstatus_filename << std::endl;
    TFile* daqstatus_file = TFile::Open(daqstatus_filename.c_str(), "READ");
    if (!daqstatus_file || daqstatus_file->IsZombie()) 
    {
        std::cerr << "Error opening file: " << daqstatus_file << std::  endl;
    }    
    else {
        std::cout << "Successfully opened file: " << daqstatus_file << std::endl;
        daqstatus_file->Close();
    }
    std::cout<<"\n"<<std::endl;

    std::stringstream combined_ss;;
    combined_ss << RNOG_DATA_DIR << "/station" << station << "/run" << run << "/combined.root";
    std::string combined_filename = combined_ss.str();
    std::cout << "Combined filename: " << combined_filename << std::endl;
    TFile* combined_file = TFile::Open(combined_filename.c_str(), "READ");
    if (!combined_file || combined_file->IsZombie()) 
    {
        std::cerr << "Error opening file: " << combined_filename << std::  endl;
    }    
    else {
        std::cout << "Successfully opened file: " << combined_filename << std::endl;
        combined_file->Close();
    }
    std::cout<<"\n"<<std::endl;

    std::stringstream pedestal_ss;;
    pedestal_ss << RNOG_DATA_DIR << "/station" << station << "/run" << run << "/pedestal.root";
    std::string pedestal_filename = pedestal_ss.str();
    std::cout << "Pedestal filename: " << pedestal_filename << std::endl;
    TFile* pedestal_file = TFile::Open(pedestal_filename.c_str(), "READ");
    if (!pedestal_file || pedestal_file->IsZombie()) 
    {
        std::cerr << "Error opening file: " << pedestal_filename << std::  endl;
    }    
    else {
        std::cout << "Successfully opened file: " << pedestal_filename << std::endl;
        pedestal_file->Close();
    }
    std::cout<<"\n"<<std::endl;

    // header_filename = std::format("{}/run_{}/station_{}/Header.root", RNOG_DATA_DIR, run, station);
    // DAQStatus_filename = std::format("{}/run_{}/station_{}/DAQStatus.root", RNOG_DATA_DIR, run, station);
    
    // TFile* header_file = TFile::Open(header_filename.c_str(), "READ");
    // if (!file || file->IsZombie()) 
    // {
    //     std::cerr << "Error opening file: " << header_filename << std::endl;
    // }
    // else 
    // {
    //     std::cout << "Successfully opened file: " << header_filename << std::endl;
    //     header_file->Close();
    // }
    // TFile* DAQStatus_file = TFile::Open(DAQStatus_filename.c_str(), "READ");
    // if (!file || file->IsZombie()) 
    // {
    //     std::cerr << "Error opening file: " << DAQStatus_filename << std::endl;
    // }    
    // else {
    //     std::cout << "Successfully opened file: " << DAQStatus_filename << std::endl;
    //     DAQStatus_file->Close();
    // }

  // TTree* tree = dynamic_cast<TTree*>(file->Get("monitoring_tree"));
  // if (!tree) {
  //   std::cerr << "Error: Tree 'monitoring_tree' not found in file: " << filename << std::endl;
  //   file->Close();
  //   return;
  // }

  // MonitoringData data;
  // tree->SetBranchAddress("run_number", &data.run_number);
  // tree->SetBranchAddress("radiant_voltage", &data.radiant_voltage);
  // tree->SetBranchAddress("battery_voltage", &data.battery_voltage);
  // tree->SetBranchAddress("wind_speed", &data.wind_speed);

  // Long64_t nentries = tree->GetEntries();
  // for (Long64_t i = 0; i < nentries; ++i) {
  //   tree->GetEntry(i);
  //   // Process the data as needed
  //   std::cout << "Run: " << data.run_number 
  //             << ", Radiant Voltage: " << data.radiant_voltage 
  //             << ", Battery Voltage: " << data.battery_voltage 
  //             << ", Wind Speed: " << data.wind_speed << std::endl;
  // }

  // file->Close();
};

// int main(){
    // const char* variable = "RNOG_DATA_DIR";
    // const char* data_dir = std::getenv(variable);
    // monitoring_obj = mattak.Monitoring(23, 999, *data_dir)
    // std::cout << "Run number: " << monitoring_obj.run_number << std::endl;
    // std::cout << "Station number: " << monitoring_obj.station_number << std::endl;
    // return 0;
// }

