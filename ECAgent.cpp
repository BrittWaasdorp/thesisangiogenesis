#include <iostream>
#include <cmath>
#include "ECAgent.h"

ECAgent::ECAgent(int x_val, int y_val, int id_val, std::vector<int> neighbours,
            std::unordered_map<std::string, double> *var_vals,
            std::unordered_map<std::string, double> *par_vals){
            
 	x       = x_val;
  	y       = y_val;
  	id      = id_val;
  	nb		= neighbours;
  	v       = std::unordered_map<std::string, double>(*var_vals);
  	p       = std::unordered_map<std::string, double>(*par_vals);
  	
};

void ECAgent::update(double dt, double sensed_dll4, double sensed_notch, double sensed_notchbound, VEGFgradient Vfield){
	
  	// Update VEGF!
  	double V_available = Vfield.value(x,y); // is actually (row,col) 
  	double a = 0.06;
  	double V_sensed = V_available*( 2/(1+exp(-a*v["filo"])) - 1); 

  	// Save concentrations before updating
  	v["dll4_saved"] = v["dll4"];
  	v["notch_saved"] = v["notch"];
  	v["notchbound_saved"] = v["notchbound"];
  	v["filo_saved"] = v["filo"]; // needed for stopping criterion
  	
  	// Determine rate eqs
  	double dll4_basal = p["beta"] * dt;
  	double notch_production = p["gamma"] * dt;
  	double vegfr2_production = p["gamma"] * dt;
  	double hey_basal = p["beta"] * dt;
  	double filo_production = p["gamma"] * dt;
  	
  	double dll4_degradation = (p["phi"] * v["dll4"]) * dt; 
  	double notch_degradation = (p["phi"] * v["notch"]) * dt;
  	double nicd_degradation = (p["phi"] * v["nicd"]) * dt;
  	double notchbound_degradation = (p["phi"] * v["notchbound"]) * dt;
  	double vegfr2_degradation = (p["phi"] * v["vegfr2"]) * dt;
  	double vegfr2bound_degradation = (p["phi"] * v["vegfr2bound"]) * dt;
  	double hey_degradation = (p["phi"] * v["hey"]) * dt;
  	double filo_degradation = (p["phi_f"] * v["filo"]) * dt;
  
  	double dll4_binding = (p["k2"] * ( p["c"] / (p["c"]+pow(p["k"]-1,p["n2"]) )) * v["dll4"] * sensed_notch ) * dt; // includes stiffness through suppressive Hill function with manual n
  	double notch_binding = (p["k2"] * ( p["c"] / (p["c"]+pow(p["k"]-1,p["n2"]) )) * sensed_dll4 * v["notch"] ) * dt; // includes stiffness through suppressive Hill function with manual n
  	double vegfr2_binding = (p["k1"] * V_sensed * v["vegfr2"]) * dt;
   
	double dll4_dissociation = (p["k_2"] * sensed_notchbound) * dt; 
  	double notch_dissociation = (p["k_2"] * v["notchbound"]) * dt;
  	double vegfr2_dissociation = (p["k_1"] * v["vegfr2bound"]) * dt;
  	
  	double nicd_catalysis = (p["k4"] * v["notchbound"]) * dt;
  	double filo_feedback = (p["k5"] * pow(v["vegfr2bound"], p["n"])) * dt;
  	double vegfr2_inhibition = (p["k3"] * v["vegfr2"] * pow(v["hey"], p["n"])) * dt;
  	double dll4_stiff_inhibition =  ( p["c1"] * v["dll4"] * ( (p["k"]-1) / (p["c2"]+p["k"]-1) ) ) * dt; // includes stiffness through Hill function
  	
  	double hey_expression = (p["theta2"] * pow(v["nicd"], p["n"])/(1 + pow(v["nicd"], p["n"]))) * dt;
  	double dll4_expression = (p["theta1"] * pow(v["vegfr2bound"], p["n"])/(1 + pow(v["vegfr2bound"], p["n"]))) * dt;
  
  	// Update internal variables
  	v["dll4"] = 		std::max(0.0, v["dll4"]			+ dll4_basal 		- dll4_degradation 			- dll4_binding 		+ dll4_dissociation		- dll4_stiff_inhibition		+ dll4_expression);
  	v["notch"] = 		std::max(0.0, v["notch"]		+ notch_production 	- notch_degradation 		- notch_binding 	+ notch_dissociation);
  	v["nicd"] = 		std::max(0.0, v["nicd"]								- nicd_degradation														+ nicd_catalysis);
  	v["notchbound"] =	std::max(0.0, v["notchbound"] 						- notchbound_degradation 	+ notch_binding 	- notch_dissociation	- nicd_catalysis); 
  	v["vegfr2"] = 		std::max(0.0, v["vegfr2"]		+ vegfr2_production - vegfr2_degradation 		- vegfr2_binding 	+ vegfr2_dissociation 	- vegfr2_inhibition);
  	v["vegfr2bound"] =	std::max(0.0, v["vegfr2bound"] 						- vegfr2bound_degradation 	+ vegfr2_binding 	- vegfr2_dissociation);
  	v["hey"] = 			std::max(0.0, v["hey"]			+ hey_basal			- hey_degradation																				+ hey_expression);
  	v["filo"] =			std::max(0.0, v["filo"]			+ filo_production	- filo_degradation														+ filo_feedback); 	
};

void ECAgent::printvars(){
        std::cout << "{dll4: " << v["dll4"] << "}, " << "{notch: " << v["notch"] << "}, " << "{nicd: " << v["nicd"] << "}, " << "{notchbound: " << v["notchbound"] << "}, " << "{vegfr2: " << v["vegfr2"] << "}, " << "{vegfr2bound: " << v["vegfr2bound"] << "}, " << "{hey: " << v["hey"] << "}, " << "{filo: " << v["filo"] << "}\n";
		std::cout << "{dll4_saved: " << v["dll4_saved"] << "}, " << "{notch_saved: " << v["notch_saved"] << "}, " << "{notchbound_saved: " << v["notchbound_saved"] << "}\n";
};

void ECAgent::printpars(){
	for (auto it = p.cbegin(); it != p.cend(); ++it) {
        std::cout << "{" << (*it).first << ": " << (*it).second << "}\n";
    }
};





