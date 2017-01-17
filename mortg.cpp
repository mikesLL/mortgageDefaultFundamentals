// Copyright A. Michael Sharifi, 2016
#include "headers.h"

mortg::mortg( void *snodes_in, double ph0_in ){

	snodes *snodes1 = (snodes *)snodes_in;

	int i_rcurr, i_rpmt, i_rlb; // current, payment, loan balance rates

	// TODO: get rid of rm_store
	double sr_mort_prem = 0.0322;                             // TODO: move this somewhere else
	rm_grid = (*snodes1).rm_gridt[0];

	for (i_rpmt = 0; i_rpmt < rm_n; i_rpmt++) {
		rm_grid[i_rpmt] = rm_grid[i_rpmt] + sr_mort_prem;
	}

	ph0 = ph0_in;                  // Load in home price
	loan_init = ltv_init*ph0;      // Compute initial loan balance
	N_term = 30;                   // mortgage term

	vector<double> zeros_RMN(rm_n, 0.0);
	vector<int> zeros_int_MN(m_n, 0);
	vector<vector<double>> zeros_RMN_T(rm_n, vector<double>(N_term, 0.0));
	vector<vector<double>> zeros_MN_T(m_n, vector<double>(N_term, 0.0));

	vector<vector<int>> zeros_int_NRM_NRM (n_rm, vector<int>(n_rm, 0));
	vector<vector<vector<int>>> zeros_int_RMN_RMN_RMN(rm_n, vector<vector<int>>(rm_n, vector<int>(rm_n, 0)));

	pmt0 = zeros_RMN;                      // initialize original mortgage pmt vector
	pmt = zeros_MN_T;                      // initialize mortgage pmt vector
	bal = zeros_RMN_T;                     // initialize mortgage balance vector
	
	m2rpmt_map = zeros_int_MN;           // mortgage index to index under repayment
	m2rlb_map = zeros_int_MN;            // mortgage index to remaining loan balance

	r2m_map = zeros_int_NRM_NRM;
	
	// initialize mortgage index maps
	i_m = 0;
	for (i_rpmt = 0; i_rpmt < n_rm; i_rpmt++) {
		for (i_rlb = 0; i_rlb < n_rm; i_rlb++) {
			m2rpmt_map[i_m] = i_rpmt;                   // pass in i_m, retrieve i_pmt state
			m2rlb_map[i_m] = i_rlb;                     // pass in i_m, retrieve i_rlb state
			r2m_map[i_rpmt][i_rlb] = i_m;
			i_m++;
		}
	}
		
	// compute initial mortgage payment
	for (i_r2 = 0; i_r2 < rm_n; i_r2++) {
		//pmt0[i_r2] = fpmt(loan_init, rm_store[i_r2], N_term);        // fmpt: computes mortgage payments
		pmt0[i_r2] = fpmt(loan_init, rm_grid[i_r2], N_term);        // fmpt: computes mortgage payments
	}

	// compute mortgage balances (assuming paid off at prevailing rate)
	for (i_r2 = 0; i_r2 < rm_n; i_r2++) { 
		t_yr = 0;
		bal[i_r2][t_yr] = loan_init;                                      // mortgage balance at origination

		for (t_yr = 1; t_yr < N_term; t_yr++) {
			bal[i_r2][t_yr] = (1.0 + rm_grid[i_r2] ) *bal[i_r2][t_yr - 1] - pmt0[i_r2];        // mortgage balance in later periods
		}
	}

	// compute mortgage payments in all states
	for (t_yr = 0; t_yr < N_term; t_yr++) {
		for (i_m = 0; i_m < m_n; i_m++ ) {
			
			i_rpmt = m2rpmt_map[i_m];  // compute payment given loan balance, payment rate, and years left
			i_rlb = m2rlb_map[i_m];		

			//pmt[i_m][t_yr] = fpmt(bal[i_rlb][t_yr], rm_store[i_rpmt], N_term - t_yr);  
			pmt[i_m][t_yr] = fpmt(bal[i_rlb][t_yr], rm_grid[i_rpmt], N_term - t_yr);  
		}
	}

	// NOTE:
	// i_m contains: i_rcurr current rate, i_r
	// given i_m and year, know: 
	// 1) current rate
	// 2) payment rate
	// 3) loan balance (assuming paid off at original rate)
	//i_rcurr = m2rcurr_map[i_m];  // error is here!
	// TODO: check index / boundary conditions
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