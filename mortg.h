// class mortg.h
// Contains code to compute mortgage payments and balances                                            
// Copyright A. Michael Sharifi, 2016

#ifndef MORTG_H
#define MORTG_H

class mortg {
	int i_m;   // index for mortgage states
	int i_r2;  // index for mortgage rates
	int t_yr;    // index for time periods
	int N_term = 30; // mortgage term 
	
	double loan_init = 1.0; // initialize mortgage bal at 100k 


public:
	vector<vector<double>> bal;  //mortgage balances indexed by rate then year 
	vector<vector<double>> pmt; // mortgage payments indexed by state then year
	vector<double> pmt0; // original mortgage payment

	mortg(); // constructor
	double fpmt(double loan_bal_in, double rm_in, int t_left_in); // mortgage payment formula
};

#endif
