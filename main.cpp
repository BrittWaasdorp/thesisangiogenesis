#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>		// std::cout
#include <iomanip>      // std::setprecision
#include <fstream>
#include <sstream>
#include <time.h>
#include <unordered_map>
#include <vector>
#include <chrono>
#include "header.h" 
#include "ECAgent.h"
#include "VEGFgradient.h"



using namespace std;

std::unordered_map<std::string, double> p;
std::unordered_map<std::string, double> v;
std::vector<ECAgent> Agentvector;
double dt = 1; // timestep of 1 second per iteration 
double TimeMatrix[40*60*60][99]; // Saves filopodia value per agent over time (for stopping criterium) // Adjust size to at least ['hourvector size'*60*60, 'Agentvector size'] 


int main()
{
    srand((unsigned)time(NULL));
    
    
	// Write down patterning per protein over stiffness
    ofstream kdll4file; 		// to save Dll4 over stiffness
	kdll4file.open("kdll4_matrix.txt");
	
	ofstream kfilofile; 		// to save filopodia over stiffness
	kfilofile.open("kfilo_matrix.txt");
	
	ofstream knotchfile; 		// to save Notch over stiffness
	knotchfile.open("knotch_matrix.txt");
	
	ofstream knotchboundfile; 	// to save notchbound over stiffness
	knotchboundfile.open("knotchbound_matrix.txt");
	
	ofstream knicdfile; 		// to save nicd over stiffness
	knicdfile.open("knicd_matrix.txt");
	
	ofstream kvegfr2file; 		// to save vegfr2 over stiffness
	kvegfr2file.open("kvegfr2_matrix.txt");
	
	ofstream kvegfr2boundfile; 	// to save vegfr2bound over stiffness
	kvegfr2boundfile.open("kvegfr2bound_matrix.txt");
	
	ofstream kheyfile; 			// to save hey over stiffness
	kheyfile.open("khey_matrix.txt");
	
	ofstream ittimefile; 		// to save real times
	ittimefile.open("ittime_matrix.txt");
	
	ofstream totalfile; 		// to save total protein values of all stiffness
	totalfile.open("total_matrix.txt");
	
	
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now(); // Tracks computational time
	
	
	VEGFgradient Vfield = VEGFgradient();
	
	std::vector<double> hourvector = {40}; // simulation time (real time in hours)
	for (int it_hour = 0; it_hour != hourvector.size(); it_hour++)
	{	
	
		// Stiffness 1 till 25 kPa, steps of 1 kPa
		std::vector<double> kvector = {1,2};
//		for (double k = 2; k != 26; k++) // stiffness 1 till 25 kPa
//		{
//			kvector.insert(kvector.end(), k);
//		}
//		// Stiffness 1 till 7 kPa, steps of 0.33 kPa
//		std::vector<double> kvector = {1, 1.33, 1.67, 2, 2.33, 2.67, 3, 3.33, 3.67, 4, 4.33, 4.67, 5, 5.33, 5.67, 6, 6.33, 6.67, 7};

		for (int it_k = 0; it_k != kvector.size(); it_k++)
		{
			double ittime;
		    int iteration;
	
			int itend = hourvector[it_hour]*60*60/dt; // Time in seconds is iterations*dt
		
		    for (iteration=0;iteration<itend;iteration++) 
		    {
			 	
			    cout<<"Iteration number: "<<iteration<<endl;
			
			    /**************************************************************************************
					Initialize vessel sprout and VEGF field
				*****************************************************************************************/ 
				int stopcode = 0; // gives info on which stopping criterium was met
				
				if (iteration==0)
			    {
				    Initialize_model(p, v, Agentvector, kvector[it_k]);
					Vfield.initialize();
			    }        

				/**************************************************************************************
		        	Update vessel sprout and VEGF field
		        *****************************************************************************************/ 
				if ((iteration)!=0) 
				{	
					// Update VEGF gradient field
					Vfield.update(dt); 
					
					// Update vessel sprout
					for(int it = 0; it != Agentvector.size(); it++)
					{
						TimeMatrix[iteration][it] = Agentvector[it].v["filo"]; // needed for stopping criterium
						
						// step 1: retrieve dll4, notch and notchbound from neighbouring cells 
						double sensed_dll4 = 0;
						double sensed_notch = 0;
						double sensed_notchbound = 0;
						for(int n = 0; n != Agentvector[it].nb.size(); n++) // iterate over neighbours of the cell
						{
							int neighbour = Agentvector[it].nb[n];
							if ( neighbour < Agentvector[it].id ) // if neighbour is updated already, take saved values from previous iteration
							{
								for (int it2 = 0; it2 != Agentvector.size(); it2++) // find iteration number of neighbour
								{
									if (Agentvector[it2].id == neighbour)
									{
										sensed_dll4 += Agentvector[it2].v["dll4_saved"];
										sensed_notch += Agentvector[it2].v["notch_saved"];
										sensed_notchbound += Agentvector[it2].v["notchbound_saved"];
									}
								}
							}
							else // if neighbour is not updated already
							{
								for (int it2 = 0; it2 != Agentvector.size(); it2++) // find iteration number of neighbour
								{
									if (Agentvector[it2].id == neighbour)
									{
										sensed_dll4 += Agentvector[it2].v["dll4"];
										sensed_notch += Agentvector[it2].v["notch"];
										sensed_notchbound += Agentvector[it2].v["notchbound"];
									}
								}
							}
						}
						
						// step 2: update cell
						Agentvector[it].update(dt, sensed_dll4, sensed_notch, sensed_notchbound, Vfield); 
						
					}
					
					/**************************************************************************************
			        	Stopping criterion
			    	*****************************************************************************************/
		    	
					ittime = iteration*dt;
					stopcode = 1;
					
					// Stopping criterion: check for all cells if filopodia changed in the last 2 hours
					
					// First 2 hours
					if (iteration > 1 && iteration*dt <= 2*60*60)
					{
						int check = 0;
						int it = 0; // counts the cells that satisfy the stopping criterium
						while (check != 1)
						{
							if ( abs(Agentvector[it].v["filo"]-Agentvector[it].v["filo_saved"]) < 10e-7 )
							{
								it += 1;
							}else // stopping criterion is not met
							{
								check = 1;
							}
							if ( it==Agentvector.size() ) // stopping criterion is met for all cells
							{
								stopcode = 2;
								iteration = itend-1;
								check = 1;
							}
						}	
					}
					// After 2 hours
					if (iteration*dt > 2*60*60)
					{
						int check = 0;
						int it = 0; // counts the cells that satisfy the stopping criterium
						while (check != 1)
						{
							int index = iteration-2*60*60/dt;
							if ( abs(Agentvector[it].v["filo"]-TimeMatrix[index][it]) < 10e-4 ) 
							{
								it += 1;
							}else // stopping criterion is not met
							{
								check = 1;
							}
							if ( it==Agentvector.size() ) // stopping criterion is met for all cells
							{
								stopcode = 3;
								ittime = iteration*dt - 2*60*60;
								iteration = itend-1;
								check = 1;
							}
						}	
					}


					/**************************************************************************************
			            Last iteration
			    	*****************************************************************************************/
					if (iteration==itend-1)
					{
						
						ittimefile << ittime << " " << stopcode << "\n"; // saves real pattern time and stopcode for stopping criterium

						double nicd_tot = 0; 
						double dll4_tot = 0;
						double notch_tot = 0; 
						double vegfr2_tot = 0; 
						double vegfr2b_tot = 0; 
						double dnb_tot = 0;
						double filo_tot = 0; 
						double he_tot = 0; 
						for(int it = 0; it != Agentvector.size(); it++)
						{
							// Adds up protein concentrations for all cells
							nicd_tot += Agentvector[it].v["nicd"];
							dll4_tot += Agentvector[it].v["dll4"];
							notch_tot += Agentvector[it].v["notch"];
							vegfr2_tot += Agentvector[it].v["vegfr2"];
							he_tot += Agentvector[it].v["hey"];
							vegfr2b_tot += Agentvector[it].v["vegfr2bound"]; 
							dnb_tot += Agentvector[it].v["notchbound"];
							filo_tot += Agentvector[it].v["filo"]; 
							
							// Saves protein concentrations in files
							kfilofile << Agentvector[it].v["filo"] << " ";
							kdll4file << Agentvector[it].v["dll4"] << " ";
							knotchfile << Agentvector[it].v["notch"] << " ";
							knotchboundfile << Agentvector[it].v["notchbound"] << " ";
							knicdfile << Agentvector[it].v["nicd"] << " ";
							kvegfr2file << Agentvector[it].v["vegfr2"] << " ";
							kvegfr2boundfile << Agentvector[it].v["vegfr2bound"] << " ";
							kheyfile << Agentvector[it].v["hey"] << " ";
						}
						// Writes down total protein concentrations
						totalfile << dll4_tot << " " << notch_tot << " " << nicd_tot << " " << dnb_tot << " " << vegfr2_tot << " " << vegfr2b_tot << " " << he_tot << " " << filo_tot << "\n";
					}
					
				}	
		    }
			
		    kdll4file << "\n";
		    kfilofile << "\n";
		    knotchfile << "\n";
		    knotchboundfile << "\n";
		    knicdfile << "\n";
		    kvegfr2file << "\n";
		    kvegfr2boundfile << "\n";
		    kheyfile << "\n";
		}
	}

	kdll4file.close();
	kfilofile.close();
	knotchfile.close();
	knotchboundfile.close();
	knicdfile.close();
	kvegfr2file.close();
	kvegfr2boundfile.close();
	kheyfile.close();
	
	ittimefile.close();
	totalfile.close();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Computation time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " [s]" << std::endl;
    
	return 0;
}


