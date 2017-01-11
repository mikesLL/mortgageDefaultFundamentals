// Copyright A. Michael Sharifi, 2016
#include "headers.h"

mortg::mortg(){

	// store zeros associated with each rate and then each time period
	vector<vector<double>> zeros_RMN_T(rm_n, vector<double>(N_term, 0.0));
	vector<vector<double>> zeros_MN_T(m_n, vector<double>(N_term, 0.0));

	vector<double> zeros_RMN(rm_n, 0.0);

	pmt0 = zeros_RMN;  // initialize original mortgage pmt vector
	bal = zeros_RMN_T; // initialize mortgage balance vector
	pmt = zeros_MN_T;  // initialize mortgage pmt vector

	// initialize maps
	// states: current rate, payment rate, loan balance rate
	/*
	i2s_map = zeros_NPH_NRENT_NYI;
	s2i_ph = zeros_NS;
	s2i_rent = zeros_NS;
	s2i_yi = zeros_NS;

	i_s = 0;

	for (i_ph = 0; i_ph < n_ph; i_ph++) {
		for (i_rent = 0; i_rent < n_rent; i_rent++) {
			for (i_yi = 0; i_yi < n_yi; i_yi++) {
				i2s_map[i_ph][i_rent][i_yi] = i_s;   // initialize i2s map
				s2i_ph[i_s] = i_ph;                 // initialize s2_i maps
				s2i_rent[i_s] = i_rent;
				s2i_yi[i_s] = i_yi;
				i_s++;
			}
		}
	}*/
	
	// compute initial mortgage payment
	for (i_r2 = 0; i_r2 < rm_n; i_r2++) {
		pmt0[i_r2] = fpmt(loan_init, rm_store[i_r2], N_term); 
	}

	// compute mortgage balances (assuming paid off at prevailing rate)
	for (i_r2 = 0; i_r2 < rm_n; i_r2++) { 
		t_yr = 0;
		bal[i_r2][t_yr] = loan_init;                                      // mortgage balance at origination

		for (t_yr = 1; t_yr < N_term; t_yr++) {
			bal[i_r2][t_yr] = bal[i_r2][t_yr - 1] - pmt0[i_r2];        // mortgage balance in later periods
		}
	}

	// compute mortgage payments in all states
	for (t_yr = 0; t_yr < N_term; t_yr++) {
		for (i_m = 0; i_m < m_n; i_m++ ) {
			// given i_m and year, know: 
			// 1) current rate
			// 2) payment rate
			// 3) loan balance (assuming paid off at original rate)
			

		}
	}
}


// mortgage payment formula
double mortg::fpmt(double loan_bal_in, double rm_in, int t_left_in) {
	double loan_bal0 = loan_bal_in;                           // loan balance
	double rm0 = rm_in;                                       // mortgage rate
	int t_left = t_left_in;                                   // years left on mortgage
	double t_leftd = double(t_left);                          

	double mnum = rm0*pow( (1.0 + rm0), t_leftd);
	double mden = pow( (1.0 + rm0), t_leftd) - 1.0;
	double mfac = mnum / mden;                                // mortgage factor

	double pmt_ret = loan_bal0 * mfac;                        // yearly mortgage payment

	return pmt_ret;
}

// const int m_n = rm_n*rm_n*rm_n; 
// Ex:
// i_m = 0 (curr rate 1, pmt rate 1, loan bal 1) loan bal 1: loan balance if ammort at rate 1
// i_m = 1 (curr rate 1, pmt rate 1, loan bal 2) loan bal 2: loan balance if ammort at rate 2 (higher rate)

// i_m = 2 (curr rate 1, pmt rate 2, loan bal 1) increase the payment rate
// i_m = 3 (curr rate 1, pmt rate 2, loan bal 2) 

// i_m = 4 (curr rate 2, pmt rate 1, loan bal 1) increase the current rate
// i_m = 5 (curr rate 2, pmt rate 1, loan bal 2) 

// i_m = 6 (curr rate 2, pmt rate 2, loan bal 1) 
// i_m = 7 (curr rate 2, pmt rate 2, loan bal 2) 

// Refinance example:
// Suppose HH begins the problem with pmt rate 2. Then loan begins with bal 2.
// Then consider the next period.
// If curr rate = 2, no benefit to refinancing. 
// If curr rate = 1, value of refinancing is value of pmt rate 1, loan bal 2

// Increase the loan balance but change the origination rate to accomodate refinancing at a lower rate
// orig rate is important because it determines payment


//const int m_n = 8; // mortgage states
// mortgage parameters: rate, adjustability, duration (term)
// TODO: write a function: input i_m, output i_m_refi 
// MOD HERE ********************************