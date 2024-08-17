#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

#include "VEGFgradient.h"

VEGFgradient::VEGFgradient()
{   
	// parameters
	h = 10; 					// <-- 10 microns between gridpoints
	D = 10;		 				// <-- 10 microns^2 per second, Prokopiou2016
	d = 1.6*pow(10,-4);			// <-- 1.6*10^(-4) per second, Prokopiou2016
	s = 0.5; 					// <-- 0.5 / 1.5
	
	// dimensions
	X = 101; 
	Y = 101; 
	
	// point source coordinates
	xp = 50;
	yp = 50;
    
};

void VEGFgradient::initialize()
{
	for(int i=0; i<X; i++) 
	{
		for(int j=0; j<Y; j++) 
		{
			matrix[i][j] = 0;
			prev_matrix[i][j] = 0;
//			// nonzero initial conditions
//			if(i == xp && j == yp){ 
//				prev_matrix[i][j] = ...;
//			}
		}	
	}
};


void VEGFgradient::update(double dt){
	
	// Update with finite difference scheme
	for(int i=0; i<X; i++)
	{
		for(int j=0; j<Y; j++)
		{
			if(i == 0)
			{
				if(j == 0) 			// (i,j) = (0,0)
				{
					matrix[i][j] = 	prev_matrix[i][j] + dt*( 1/(h*h) * D*(													- 4*prev_matrix[i][j] 	+ prev_matrix[i+1][j] 	+ prev_matrix[i][j+1]) - d*prev_matrix[i][j]);
				}else if(j == Y-1) 	// (i,j) = (0,Y-1)
				{
					matrix[i][j] = 	prev_matrix[i][j] + dt*( 1/(h*h) * D*(							  prev_matrix[i][j-1] 	- 4*prev_matrix[i][j] 	+ prev_matrix[i+1][j]						 ) - d*prev_matrix[i][j]);
				}else 				// (i,j) = (0,j)
				{
					matrix[i][j] = 	prev_matrix[i][j] + dt*( 1/(h*h) * D*(							  prev_matrix[i][j-1] 	- 4*prev_matrix[i][j] 	+ prev_matrix[i+1][j] 	+ prev_matrix[i][j+1]) - d*prev_matrix[i][j]);
				}
			}else if(i == X-1)
			{
				if(j == 0) 			// (i,j) = (X-1,0)
				{
					matrix[i][j] = 	prev_matrix[i][j] + dt*( 1/(h*h) * D*(	prev_matrix[i-1][j] 							- 4*prev_matrix[i][j] 							+ prev_matrix[i][j+1]) - d*prev_matrix[i][j]);
				}else if(j == Y-1) 	// (i,j) = (X-1,Y-1)
				{
					matrix[i][j] = 	prev_matrix[i][j] + dt*( 1/(h*h) * D*(	prev_matrix[i-1][j] 	+ prev_matrix[i][j-1] 	- 4*prev_matrix[i][j]												 ) - d*prev_matrix[i][j]);
				}else 				// (i,j) = (X-1,j)
				{
					matrix[i][j] = 	prev_matrix[i][j] + dt*( 1/(h*h) * D*(	prev_matrix[i-1][j] 	+ prev_matrix[i][j-1] 	- 4*prev_matrix[i][j] 							+ prev_matrix[i][j+1]) - d*prev_matrix[i][j]);
				}
			}else					
			{
				if(j == 0)			// (i,j) = (i,0)
				{
					matrix[i][j] = 	prev_matrix[i][j] + dt*( 1/(h*h) * D*(	prev_matrix[i-1][j] 							- 4*prev_matrix[i][j] 	+ prev_matrix[i+1][j]	+ prev_matrix[i][j+1]) - d*prev_matrix[i][j]);
				}else if(j == Y-1) 	// (i,j) = (i,Y-1)
				{
					matrix[i][j] = 	prev_matrix[i][j] + dt*( 1/(h*h) * D*(	prev_matrix[i-1][j] 	+ prev_matrix[i][j-1] 	- 4*prev_matrix[i][j]	+ prev_matrix[i+1][j]						 ) - d*prev_matrix[i][j]);
				}else if(i == xp && j == yp) // (i,j) = (xp,yp) point source
				{
					matrix[i][j] = 	prev_matrix[i][j] + dt*( 1/(h*h) * D*(	prev_matrix[i-1][j] 	+ prev_matrix[i][j-1] 	- 4*prev_matrix[i][j] 	+ prev_matrix[i+1][j] 	+ prev_matrix[i][j+1]) - d*prev_matrix[i][j] + s);
				}else 				// (i,j) = (i,j)
				{
					matrix[i][j] = 	prev_matrix[i][j] + dt*( 1/(h*h) * D*(	prev_matrix[i-1][j] 	+ prev_matrix[i][j-1] 	- 4*prev_matrix[i][j] 	+ prev_matrix[i+1][j] 	+ prev_matrix[i][j+1]) - d*prev_matrix[i][j]);
				}
			}
		}	
	}
	
	for(int i=0; i<X; i++)
	{
		for(int j=0; j<Y; j++)
		{
			prev_matrix[i][j] = matrix[i][j];
		}	
	} 
	
};


double VEGFgradient::value(int x, int y){	// actually (row,col)
	return matrix[x][y];
}




