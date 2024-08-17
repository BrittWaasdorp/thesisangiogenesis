#ifndef VEGFGRADIENT_H
#define VEGFGRADIENT_H

#include <cstdlib>
#include <iostream>
#include <vector>

class VEGFgradient{
public:
	double matrix[101][101], prev_matrix[101][101];
	int X, Y, xp, yp;
	double h, D, d, s;
	
    // ___________________________________________________________________________
    // constructors:
    VEGFgradient();

    // ___________________________________________________________________________
    // member functions:
  	~VEGFgradient() {};
	
	void initialize();
	
	void update(double dt);
	
	double value(int x, int y);
};

#endif
