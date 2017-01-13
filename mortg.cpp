// Copyright A. Michael Sharifi, 2016
#include "headers.h"

mortg::mortg(){

	cout << "large object created" << endl; 
	cin.get();

	// store zeros associated with each rate and then each time period
	vector<double> zeros_RMN(rm_n, 0.0);
	vector<int> zeros_int_MN(m_n, 0);
	vector<vector<double>> zeros_RMN_T(rm_n, vector<double>(N_term, 0.0));
	vector<vector<double>> zeros_MN_T(m_n, vector<double>(N_term, 0.0));
	vector<vector<vector<int>>> zeros_int_RMN_RMN_RMN(rm_n, vector<vector<int>>(rm_n, vector<int>(rm_n, 0)));

	pmt0 = zeros_RMN;  // initialize original mortgage pmt vector
	bal = zeros_RMN_T; // initialize mortgage balance vector
	pmt = zeros_MN_T;  // initialize mortgage pmt vector
	
	m2rcurr_map = zeros_int_MN;
	m2rpmt_map = zeros_int_MN;
	m2rlb_map = zeros_int_MN;

	r2m_map = zeros_int_RMN_RMN_RMN;
	m2mrefi_map = zeros_int_MN;

	//vector<vector<vector<int>>> r2m_map(rm_n, vector<vector<int>> (rm_n, vector<int> (rm_n, 0)));
	// state maps
	//vector<int> m2rcurr_map(m_n, 0);
	//vector<int> m2rpmt_map(m_n, 0);
	//vector<int> m2rlb_map(m_n, 0);
	//vector<vector<vector<int>>> r2m_map(rm_n, vector<vector<int>> (rm_n, vector<int> (rm_n, 0)));

	//vector<int> m2mrefi_map(m_n, 0); // state map: m to m refi

	// initialize maps
	// states: current rate, payment rate, loan balance rate
	int i_m_refi;
	i_m = 0;

	int i_rcurr; // current rate
	int i_rpmt;  // payment rate
	int i_rlb;   // loan balance assuming payment at this rate

	//r2m_map[1][1][1] = 0;
	
	for (i_rcurr = 0; i_rcurr < rm_n; i_rcurr++) {
		for (i_rpmt = 0; i_rpmt < rm_n; i_rpmt++) {
			for (i_rlb = 0; i_rlb < rm_n; i_rlb++) {
				m2rcurr_map[i_m] = i_rcurr; // pass in i_m, retrieve i_rcurr state
				m2rpmt_map[i_m] = i_rpmt; // pass in i_m, retrieve i_pmt state
				m2rlb_map[i_m] = i_rlb; // pass in i_m, retrieve i_rlb state

				r2m_map[i_rcurr][i_rpmt][i_rlb] = i_m;

				i_m++;
			}
		}
	}


	int foo = r2m_map[0][0][0];
	// initialize refinance map
	i_m = 0;
	// int i_m_refi; // i_m associated with refinancing
	// note that due to indexing r2m_map must be pre-computed 
	for (i_rcurr = 0; i_rcurr < rm_n; i_rcurr++) {
		for (i_rpmt = 0; i_rpmt < rm_n; i_rpmt++) {
			for (i_rlb = 0; i_rlb < rm_n; i_rlb++) {
				// refinance: retrieve i_m_refi where i_lb, i_rcurr fixed but pmt rate = curr rate
				m2mrefi_map[i_m] = r2m_map[i_rcurr][i_rcurr][i_rlb]; 
				i_m++;
			}
		}
	}

	// Q: will the value fn at the lower rate have already been computed?
	// Q: how to compute the difference in loan balances at time of refinance?
	// bal[i_r_refi] = ...
	// bal[i_r] = ...
	// difference modeled as a change in liquid assets

	// Value of refinancing: eval value fn at i_m_refi
	
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
			i_rcurr = m2rcurr_map[i_m];
			i_rpmt = m2rpmt_map[i_m];
			i_rlb = m2rlb_map[i_m];		

			// compute payment given loan balance, payment rate, and years left
			// TODO: check index / boundary conditions
			pmt[i_m][t_yr] = fpmt(bal[i_rlb][t_yr], rm_store[i_rpmt], N_term - t_yr);  
		}
	}
}

/*
Now, given mortgage payment in each state can compute the value function
Write a method to determine the value of refinancing
*/

// int mortg::

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