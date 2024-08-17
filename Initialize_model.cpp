#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <limits> 
#include <unordered_map>
#include "header.h"
#include "ECAgent.h"


using namespace std;


void Initialize_model(	std::unordered_map<std::string, double> &p,
						std::unordered_map<std::string, double> &v,
						std::vector<ECAgent> &Agentvector,
						double k)
{
   
//*********************************************************************************************************
//                                      Initialize ECs
//*********************************************************************************************************   
     
	// parameters
	p["k1"] = 0.1;
  	p["k_1"] = 0.001;
  	p["k2"] = 0.001;
  	p["k_2"] = 0.1;
  	p["k3"] = 0.005;
  	p["k4"] = 0.1;
  	p["k5"] = 0.1; 
  	p["k6"] = 0.005; 
  	p["beta"] = 0.001;
  	p["gamma"] = 0.005;
  	p["phi"] = 0.005;
	p["theta1"] = 0.2;
	p["theta2"] = 0.15;
  	p["phi_f"] = 0.001;
  	p["n"] = 2;
  
	// parameters stiffness-dependent functions
	p["k"] = k; 
	p["c1"] = 0.12; 
	p["c2"] = 52; 
	p["c"] = 1.26;
	p["n2"] = 0.92;
  	
	// variables
    v["dll4"] = 0;
	v["notch"] = 1;
	v["nicd"] = 0;
	v["notchbound"] = 0;
    v["vegfr2"] = 1;
	v["vegfr2bound"] = 0;
	v["hey"] = 0;
	v["filo"] = 0;
	v["dll4_saved"] = 0;
	v["notch_saved"] = 1;
	v["notchbound_saved"] = 0;
	v["filo_saved"] = 0;
	v["vegf_sensed"] = 0;
	v["vegf_sensed_tot"] = 0;

//*********************************************************************************************************
//                                     Define set-up
//									-------------------
//									| 0,0 0,1 ... 0,n |
//									| 1,0 	  ...	  |
//									|				  |
//									|				  |
//									|				  |
//									| n,0     ... n,n |
//									------------------
//								Matlab plot starts from 1,1
//*********************************************************************************************************   

	// Place EC Agents
	
	int row = 68;
	
	ECAgent Agent0 = ECAgent(row, 45, 1, {2,10}, &v, &p); 
	ECAgent Agent1 = ECAgent(row, 46, 2, {1,3}, &v, &p); 
    ECAgent Agent2 = ECAgent(row, 47, 3, {2,4}, &v, &p); 
    ECAgent Agent3 = ECAgent(row, 48, 4, {3,5}, &v, &p); 
    ECAgent Agent4 = ECAgent(row, 49, 5, {4,6}, &v, &p);
    ECAgent Agent5 = ECAgent(row, 50, 6, {5,7}, &v, &p); 
    ECAgent Agent6 = ECAgent(row, 51, 7, {6,8}, &v, &p); 
    ECAgent Agent7 = ECAgent(row, 52, 8, {7,9}, &v, &p); 
    ECAgent Agent8 = ECAgent(row, 53, 9, {8,10}, &v, &p); 
    ECAgent Agent9 = ECAgent(row, 54, 10, {9,11}, &v, &p); 
    ECAgent Agent10 = ECAgent(row, 55, 11, {10,1}, &v, &p); 
    
    Agentvector = {Agent0, Agent1, Agent2, Agent3, Agent4, Agent5, Agent6, Agent7, Agent8, Agent9, Agent10}; 

	
}

