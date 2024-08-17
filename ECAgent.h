#ifndef ECAGENT_H
#define ECAGENT_H

#include <algorithm>    // using, at least, std::max_element
#include <random>
#include <ctime>
#include <vector>
#include <string>
#include <unordered_map>
#include "sys/time.h"
#include "VEGFgradient.h"


class ECAgent{
public:
	int x, y, id;
	std::vector<int> nb;
	std::unordered_map<std::string, double> v; // variables
	std::unordered_map<std::string, double> p; // parameters
	
    // ___________________________________________________________________________
    // constructors:
    ECAgent(int x_val, int y_val, int id_val, std::vector<int> neighbours,
            std::unordered_map<std::string, double> *var_vals,
            std::unordered_map<std::string, double> *par_vals);

    // ___________________________________________________________________________
    // member functions:
  	~ECAgent() {}; 

    void update(double dt, double sensed_dll4, double sensed_notch, double sensed_notchbound, VEGFgradient Vfield);
  	
	void printvars();
	
	void printpars();
};

#endif


