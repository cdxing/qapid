#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "ConfigReader.h"

// I have to initialize the maps like this because c++98 sucks (or maybe just the CINT does)
void ConfigReader::initialize()
{
  intValCuts["fixed_target"] = -999;
  intValCuts["epd_max_weight"] = -999;
  intValCuts["nHits"] = -999;
  intValCuts["nHits_dEdx"] = -999;
  intValCuts["min_tracks"] = -999;
  intValCuts["shift_terms"] = -999;
  intValCuts["epdA_inner_row"] = -999;
  intValCuts["epdA_outer_row"] = -999;
  intValCuts["epdB_inner_row"] = -999;
  intValCuts["epdB_outer_row"] = -999;

  dblValCuts["sqrt_s_NN"] = -999.0;
  dblValCuts["order_n"] = -999.0; 
  dblValCuts["order_m"] = -999.0; 
  dblValCuts["epd_threshold"] = -999.0; 
  dblValCuts["nHits_ratio"] = -999.0; 
  dblValCuts["dca"] = -999.0; 
  dblValCuts["tpc_A_low_eta"] = -999.0; 
  dblValCuts["tpc_A_high_eta"] = -999.0; 
  dblValCuts["tpc_B_low_eta"] = -999.0; 
  dblValCuts["tpc_B_high_eta"] = -999.0; 
  dblValCuts["r_vtx"] = -999.0; 
  dblValCuts["z_vtx_low"] = -999.0; 
  dblValCuts["z_vtx_high"] = -999.0; 
  dblValCuts["y_mid"] = -999.0; 
  dblValCuts["nSig_pi_low"] = -999.0; 
  dblValCuts["nSig_pi_high"] = -999.0; 
  dblValCuts["nSig_ka_low"] = -999.0; 
  dblValCuts["nSig_ka_high"] = -999.0; 
  dblValCuts["nSig_pr_low"] = -999.0; 
  dblValCuts["nSig_pr_high"] = -999.0; 
  dblValCuts["z_de_low"] = -999.0; 
  dblValCuts["z_de_high"] = -999.0; 
  dblValCuts["z_tr_low"] = -999.0; 
  dblValCuts["z_tr_high"] = -999.0; 
  dblValCuts["m2_pi_low"] = -999.0; 
  dblValCuts["m2_pi_high"] = -999.0; 
  dblValCuts["m2_ka_low"] = -999.0; 
  dblValCuts["m2_ka_high"] = -999.0; 
  dblValCuts["m2_de_low"] = -999.0; 
  dblValCuts["m2_de_high"] = -999.0; 
  dblValCuts["m2_tr_low"] = -999.0; 
  dblValCuts["m2_tr_high"] = -999.0; 
  dblValCuts["yCM_pid_pi_low"] = -999.0;
  dblValCuts["yCM_pid_pi_high"] = -999.0;
  dblValCuts["yCM_flow_pi_low"] = -999.0;
  dblValCuts["yCM_flow_pi_high"] = -999.0;
  dblValCuts["yCM_ext_flow_pi_low"] = -999.0;
  dblValCuts["yCM_ext_flow_pi_high"] = -999.0;
  dblValCuts["yCM_pid_de_low"] = -999.0;
  dblValCuts["yCM_pid_de_high"] = -999.0;
  dblValCuts["yCM_flow_de_low"] = -999.0;
  dblValCuts["yCM_flow_de_high"] = -999.0;
  dblValCuts["yCM_ext_flow_de_low"] = -999.0;
  dblValCuts["yCM_ext_flow_de_high"] = -999.0;
  dblValCuts["yCM_pid_tr_low"] = -999.0;
  dblValCuts["yCM_pid_tr_high"] = -999.0;
  dblValCuts["yCM_flow_tr_low"] = -999.0;
  dblValCuts["yCM_flow_tr_high"] = -999.0;
  dblValCuts["yCM_ext_flow_tr_low"] = -999.0;
  dblValCuts["yCM_ext_flow_tr_high"] = -999.0;
  dblValCuts["yCM_pid_ka_low"] = -999.0;
  dblValCuts["yCM_pid_ka_high"] = -999.0;
  dblValCuts["yCM_flow_ka_low"] = -999.0;
  dblValCuts["yCM_flow_ka_high"] = -999.0;
  dblValCuts["yCM_ext_flow_ka_low"] = -999.0;
  dblValCuts["yCM_ext_flow_ka_high"] = -999.0;
  dblValCuts["yCM_pid_pr_low"] = -999.0;
  dblValCuts["yCM_pid_pr_high"] = -999.0;
  dblValCuts["yCM_flow_pr_low"] = -999.0;
  dblValCuts["yCM_flow_pr_high"] = -999.0;
  dblValCuts["yCM_dep_flow_pr_low"] = -999.0;
  dblValCuts["yCM_dep_flow_pr_high"] = -999.0;
  dblValCuts["yCM_ext_flow_pr_low"] = -999.0;
  dblValCuts["yCM_ext_flow_pr_high"] = -999.0;
  dblValCuts["yCM_sym_flow_pr_low"] = -999.0;
  dblValCuts["yCM_sym_flow_pr_high"] = -999.0;
  dblValCuts["yCM_for_flow_pr_low"] = -999.0;
  dblValCuts["yCM_for_flow_pr_high"] = -999.0;
  dblValCuts["pt_pid_pi_low"] = -999.0; 
  dblValCuts["pt_pid_pi_high"] = -999.0; 
  dblValCuts["pt_pid_ka_low"] = -999.0; 
  dblValCuts["pt_pid_ka_high"] = -999.0; 
  dblValCuts["pt_pid_pr_low"] = -999.0; 
  dblValCuts["pt_pid_pr_high"] = -999.0; 
  dblValCuts["pt_flow_pr_low"] = -999.0; 
  dblValCuts["pt_flow_pr_high"] = -999.0; 
  dblValCuts["pt_ydep_flow_pr_low"] = -999.0;
  dblValCuts["pt_ydep_flow_pr_high"] = -999.0;
  dblValCuts["pt_ext_flow_pr_low"] = -999.0;
  dblValCuts["pt_ext_flow_pr_high"] = -999.0;
  dblValCuts["pt_sym_flow_pr_low"] = -999.0;
  dblValCuts["pt_sym_flow_pr_high"] = -999.0;
  dblValCuts["pt_for_flow_pr_low"] = -999.0;
  dblValCuts["pt_for_flow_pr_high"] = -999.0;
  dblValCuts["pt_pid_de_low"] = -999.0; 
  dblValCuts["pt_pid_de_high"] = -999.0; 
  dblValCuts["pt_pid_tr_low"] = -999.0; 
  dblValCuts["pt_pid_tr_high"] = -999.0; 
}

void ConfigReader::setAllCuts()
{
  fixed_target = intValCuts["fixed_target"];
  epd_max_weight = intValCuts["epd_max_weight"];
  nHits = intValCuts["nHits"];
  nHits_dEdx = intValCuts["nHits_dEdx"];
  min_tracks = intValCuts["min_tracks"];
  shift_terms = intValCuts["shift_terms"];
  epdA_inner_row = intValCuts["epdA_inner_row"];
  epdA_outer_row = intValCuts["epdA_outer_row"];
  epdB_inner_row = intValCuts["epdB_inner_row"];
  epdB_outer_row = intValCuts["epdB_outer_row"];

  sqrt_s_NN = dblValCuts["sqrt_s_NN"];
  order_n = dblValCuts["order_n"]; 
  order_m = dblValCuts["order_m"]; 
  epd_threshold = dblValCuts["epd_threshold"]; 
  nHits_ratio = dblValCuts["nHits_ratio"]; 
  dca = dblValCuts["dca"]; 
  tpc_A_low_eta = dblValCuts["tpc_A_low_eta"]; 
  tpc_A_high_eta = dblValCuts["tpc_A_high_eta"]; 
  tpc_B_low_eta = dblValCuts["tpc_B_low_eta"]; 
  tpc_B_high_eta = dblValCuts["tpc_B_high_eta"]; 
  r_vtx = dblValCuts["r_vtx"]; 
  z_vtx_low = dblValCuts["z_vtx_low"]; 
  z_vtx_high = dblValCuts["z_vtx_high"]; 
  y_mid = dblValCuts["y_mid"]; 
  nSig_pi_low = dblValCuts["nSig_pi_low"]; 
  nSig_pi_high = dblValCuts["nSig_pi_high"]; 
  nSig_ka_low = dblValCuts["nSig_ka_low"]; 
  nSig_ka_high = dblValCuts["nSig_ka_high"]; 
  nSig_pr_low = dblValCuts["nSig_pr_low"]; 
  nSig_pr_high = dblValCuts["nSig_pr_high"]; 
  z_de_low = dblValCuts["z_de_low"]; 
  z_de_high = dblValCuts["z_de_high"]; 
  z_tr_low = dblValCuts["z_tr_low"]; 
  z_tr_high = dblValCuts["z_tr_high"]; 
  m2_pi_low = dblValCuts["m2_pi_low"]; 
  m2_pi_high = dblValCuts["m2_pi_high"]; 
  m2_ka_low = dblValCuts["m2_ka_low"]; 
  m2_ka_high = dblValCuts["m2_ka_high"]; 
  m2_de_low = dblValCuts["m2_de_low"]; 
  m2_de_high = dblValCuts["m2_de_high"]; 
  m2_tr_low = dblValCuts["m2_tr_low"]; 
  m2_tr_high = dblValCuts["m2_tr_high"]; 
  yCM_pid_pi_low = dblValCuts["yCM_pid_pi_low"];
  yCM_pid_pi_high = dblValCuts["yCM_pid_pi_high"];
  yCM_flow_pi_low = dblValCuts["yCM_flow_pi_low"];
  yCM_flow_pi_high = dblValCuts["yCM_flow_pi_high"];
  yCM_ext_flow_pi_low = dblValCuts["yCM_ext_flow_pi_low"];
  yCM_ext_flow_pi_high = dblValCuts["yCM_ext_flow_pi_high"];
  yCM_pid_de_low = dblValCuts["yCM_pid_de_low"];
  yCM_pid_de_high = dblValCuts["yCM_pid_de_high"];
  yCM_flow_de_low = dblValCuts["yCM_flow_de_low"];
  yCM_flow_de_high = dblValCuts["yCM_flow_de_high"];
  yCM_ext_flow_de_low = dblValCuts["yCM_ext_flow_de_low"];
  yCM_ext_flow_de_high = dblValCuts["yCM_ext_flow_de_high"];
  yCM_pid_tr_low = dblValCuts["yCM_pid_tr_low"];
  yCM_pid_tr_high = dblValCuts["yCM_pid_tr_high"];
  yCM_flow_tr_low = dblValCuts["yCM_flow_tr_low"];
  yCM_flow_tr_high = dblValCuts["yCM_flow_tr_high"];
  yCM_ext_flow_tr_low = dblValCuts["yCM_ext_flow_tr_low"];
  yCM_ext_flow_tr_high = dblValCuts["yCM_ext_flow_tr_high"];
  yCM_pid_ka_low = dblValCuts["yCM_pid_ka_low"];
  yCM_pid_ka_high = dblValCuts["yCM_pid_ka_high"];
  yCM_flow_ka_low = dblValCuts["yCM_flow_ka_low"];
  yCM_flow_ka_high = dblValCuts["yCM_flow_ka_high"];
  yCM_ext_flow_ka_low = dblValCuts["yCM_ext_flow_ka_low"];
  yCM_ext_flow_ka_high = dblValCuts["yCM_ext_flow_ka_high"];
  yCM_pid_pr_low = dblValCuts["yCM_pid_pr_low"];
  yCM_pid_pr_high = dblValCuts["yCM_pid_pr_high"];
  yCM_flow_pr_low = dblValCuts["yCM_flow_pr_low"];
  yCM_flow_pr_high = dblValCuts["yCM_flow_pr_high"];
  yCM_dep_flow_pr_low = dblValCuts["yCM_dep_flow_pr_low"];
  yCM_dep_flow_pr_high = dblValCuts["yCM_dep_flow_pr_high"];
  yCM_ext_flow_pr_low = dblValCuts["yCM_ext_flow_pr_low"];
  yCM_ext_flow_pr_high = dblValCuts["yCM_ext_flow_pr_high"];
  yCM_sym_flow_pr_low = dblValCuts["yCM_sym_flow_pr_low"];
  yCM_sym_flow_pr_high = dblValCuts["yCM_sym_flow_pr_high"];
  yCM_for_flow_pr_low = dblValCuts["yCM_for_flow_pr_low"];
  yCM_for_flow_pr_high = dblValCuts["yCM_for_flow_pr_high"];
  pt_pid_pi_low = dblValCuts["pt_pid_pi_low"]; 
  pt_pid_pi_high = dblValCuts["pt_pid_pi_high"]; 
  pt_pid_ka_low = dblValCuts["pt_pid_ka_low"]; 
  pt_pid_ka_high = dblValCuts["pt_pid_ka_high"]; 
  pt_pid_pr_low = dblValCuts["pt_pid_pr_low"]; 
  pt_pid_pr_high = dblValCuts["pt_pid_pr_high"]; 
  pt_flow_pr_low = dblValCuts["pt_flow_pr_low"]; 
  pt_flow_pr_high = dblValCuts["pt_flow_pr_high"]; 
  pt_ydep_flow_pr_low = dblValCuts["pt_ydep_flow_pr_low"];
  pt_ydep_flow_pr_high = dblValCuts["pt_ydep_flow_pr_high"];
  pt_ext_flow_pr_low = dblValCuts["pt_ext_flow_pr_low"];
  pt_ext_flow_pr_high = dblValCuts["pt_ext_flow_pr_high"];
  pt_sym_flow_pr_low = dblValCuts["pt_sym_flow_pr_low"];
  pt_sym_flow_pr_high = dblValCuts["pt_sym_flow_pr_high"];
  pt_for_flow_pr_low = dblValCuts["pt_for_flow_pr_low"];
  pt_for_flow_pr_high = dblValCuts["pt_for_flow_pr_high"];
  pt_pid_de_low = dblValCuts["pt_pid_de_low"]; 
  pt_pid_de_high = dblValCuts["pt_pid_de_high"]; 
  pt_pid_tr_low = dblValCuts["pt_pid_tr_low"]; 
  pt_pid_tr_high = dblValCuts["pt_pid_tr_high"]; 
}

ConfigReader::ConfigReader() 
{ 
  initialize();
  errorFlag = false; 
}

ConfigReader::~ConfigReader()
{}

bool ConfigReader::errorFound()
{ return errorFlag; }

void ConfigReader::notifyError()
{
  std::cout << std::endl;
  std::cout << "There was an error in reading the config file." << std::endl;

  if (lastKey != "")
    { std::cout << "The last key and value read were: " << lastKey << ", " << lastValue << std::endl; }
  else
    { std::cout << "There were no keys or values read successfully." << std::endl; }
}

Bool_t ConfigReader::triggersMatch(UInt_t readTrigger)
{
  Bool_t triggerMatchFound = false;
  for (unsigned int i = 0; i < triggers.size(); i++)
    { 
      if (readTrigger == triggers[i]) 
	{ 
	  triggerMatchFound = true; 
	  break;
	} 
    }
  return triggerMatchFound;
}

void ConfigReader::read(std::string fileName)
{
  std::ifstream inputStream(fileName.c_str());

  std::string line;
  std::getline(inputStream, line);  // Get the first line

  // Loop over lines of input in current file
  while (inputStream.good())
    {
      // Skip the text lines and the empty lines
      if (line[0] == '#' || line.empty())
	{
	  std::getline(inputStream, line);
	  continue;
	}

      // Split string by delimeter '='
      std::string delimeter = "=";
      size_t pos = line.find(delimeter);

      if (pos == std::string::npos) 
	{ 
	  std::cout << "Missing \'=\' in a line of the config file." << std::endl; 
	  errorFlag = true;
	  break;
	}

      std::string key   = line.substr(0, pos);
      std::string value = line.substr(pos+1, std::string::npos);
      
      lastKey   = key;
      lastValue = value;

      // Now check if "value" is actually a list of values, if so, split it and try to save the triggers.
      size_t commaPos = value.find(",");

      if (commaPos != std::string::npos && key.compare("triggers") == 0)
	{
	  std::stringstream valueStream(value);

	  while(valueStream.good())
	    {
	      std::string subString;
	      std::getline(valueStream, subString, ',');
	      try
		{ 
		  UInt_t triggerValue = (UInt_t)std::atoi(subString.c_str()); 
		  triggers.push_back(triggerValue);
		}
	      catch (...)
		{
		  std::cout << "Error parsing this value: " << value << std::endl;
		  errorFlag = true;
		  break;
		}
	    }
	}
      else if (key.compare("triggers") == 0) // No commas, so only one trigger
	{
	  try
	    { 
	      UInt_t triggerValue = (UInt_t)std::atoi(value.c_str()); 
	      triggers.push_back(triggerValue); 
	    }
	  catch (...)
	    {
	      std::cout << "Error in this value for triggers: " << value << std::endl;
	      errorFlag = true;
	    }
	}
      if (errorFlag) break;

      if (key.compare("triggers") != 0)
	{
	  try
	    { intValCuts.at(key) = std::atoi(value.c_str()); }
	  catch (...)//(const std::out_of_range& oorInt) 
	    {
	      try
		{ dblValCuts.at(key) = std::atof(value.c_str()); }
	      catch (...)//(const std::out_of_range& oorDbl) 
		{
		  std::cout << "Unknown config key: " << key << std::endl;
		  errorFlag = true;
		  break;	      
		}
	    }
	}

      std::getline(inputStream, line); // Get the next line
    }// End while(inputStream.good())

  if (errorFlag) { notifyError(); }
  setAllCuts();
}// End function read()
