// class mortg.h
// Contains code to compute mortgage payments and balances                                            
// Copyright A. Michael Sharifi, 2016

#ifndef MORTG_H
#define MORTG_H

class mortg {
	
	double mapr;          // mortgage apr
	int N_term;           // mortgage term 
	double ph0;          // initial home price


public:
	double mpmt;                        // mortgage payment
	vector<double> loan_bal;            // mortgage loan balance (by year)

	mortg( void *snodes_in, double ph0_in, double ltv0_in, double mapr_in, int N_term_in); // constructor

	double fpmt(double loan_bal_in, double rm_in, int t_left_in);                         // mortgage payment formula
};

#endif


//int i_m;   // index for mortgage states
//int i_r2;  // index for mortgage rates
//int t_yr;    // index for time periods

//vector<vector<double>> bal;  //mortgage balances indexed by rate then year 
//vector<vector<double>> pmt; // mortgage payments indexed by state then year
//vector<double> pmt0; // original mortgage payment

// mortgage state maps
//vector<int> m2rcurr_map;
//vector<int> m2rpmt_map;
//vector<int> m2rlb_map;

//vector<vector<vector<int>>> r2m_map;
//vector<vector<int>> r2m_map;

//vector<int> m2mrefi_map; // state map: m to m refi