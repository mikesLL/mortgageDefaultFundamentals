// class mortg.h
// Contains code to compute mortgage payments and balances                                            
// Copyright A. Michael Sharifi, 2016

#ifndef MORTG_H
#define MORTG_H

class mortg {
	int i_m;   // index for mortgage states
	int i_r2;  // index for mortgage rates
	int t_yr;    // index for time periods
	int N_term;        // mortgage term 
	double loan_init;    // initialize mortgage bal
	double ph0;          // initial home price
	
	vector<double> rm_grid;  // mortgage rate grid

public:
	vector<vector<double>> bal;  //mortgage balances indexed by rate then year 
	vector<vector<double>> pmt; // mortgage payments indexed by state then year
	vector<double> pmt0; // original mortgage payment

	// mortgage state maps
	vector<int> m2rcurr_map;
	vector<int> m2rpmt_map;
	vector<int> m2rlb_map;
	//vector<vector<vector<int>>> r2m_map;
	vector<vector<int>> r2m_map;

	//vector<int> m2mrefi_map; // state map: m to m refi

	mortg( void *snodes_in, double ph0_in, double ltv0_in ); // constructor
	double fpmt(double loan_bal_in, double rm_in, int t_left_in); // mortgage payment formula
};

#endif
