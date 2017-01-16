// Copyright A. Michael Sharifi, 2016
#include "headers.h"

snodes::snodes(int age0_in, int T_max_in, int city_id_in) {

	city_id = city_id_in;
	age0 = age0_in;
	T_max = T_max_in;

	csfLevSn = csfLev_pidxw[param_id] * 1.0 / csfmarg_store[city_id];

	// initialize price, rent, income grids
	// use max time horizon, n_ph, n_rent, n_yi gridpoints
	rm_gridt = vector<vector<double>>(T_max + 1, vector<double>(n_rm, 0.0));
	p_gridt = vector<vector<double>>(T_max + 1, vector<double>(n_ph, 0.0));
	rent_gridt = vector<vector<double>>(T_max + 1, vector<double>(n_rent, 0.0));
	yi_gridt = vector<vector<double>>(T_max + 1, vector<double>(n_yi, 0.0));
	yi_gridt_btax = vector<vector<double>>(T_max + 1, vector<double>(n_yi, 0.0));
	

	// initialize gammat state-state transition matrix
	vector<vector<vector<double>>> zeros_T_NS_NS(T_max + 1, vector<vector<double>>(n_s, vector<double>(n_s, 0.0)));
	gammat = zeros_T_NS_NS;

	// initialize i2s and s2i mappings
	vector<vector<vector<int>>> zeros_NPH_NRENT_NYI(n_ph, vector<vector<int>>(n_rent, vector<int>(n_yi, 0)));

	vector<vector<vector<vector<int>>>> zeros_NPH_NRENT_NYI_NRM(n_ph, vector<vector<vector<int>>>(n_rent, vector<vector<int>>(n_yi, vector<int>(n_rm, 0))));
	
	vector<int> zeros_NS(n_s, 0);

	i2s_map = zeros_NPH_NRENT_NYI_NRM;
	//i2s_map = zeros_NPH_NRENT_NYI;
	s2i_ph = zeros_NS;
	s2i_rent = zeros_NS;
	s2i_yi = zeros_NS;
	s2i_rm = zeros_NS;

	i_s = 0;

	/*
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

	for (i_ph = 0; i_ph < n_ph; i_ph++) {
		for (i_rent = 0; i_rent < n_rent; i_rent++) {
			for (i_yi = 0; i_yi < n_yi; i_yi++) {
				for (i_rm = 0; i_rm < n_rm; i_rm++) {
					i2s_map[i_ph][i_rent][i_yi][i_rm] = i_s;      // initialize i2s map
					s2i_ph[i_s] = i_ph;                           // initialize s2_i maps
					s2i_rent[i_s] = i_rent;     
					s2i_yi[i_s] = i_yi;
					s2i_rm[i_s] = i_rm; 
					i_s++;
				}
			}
		}
	}

	// compute i_s_mid;
	i_ph_mid = (int) floor( 0.5*n_ph );
	i_rent_mid = (int) floor( 0.5*n_rent );
	i_yi_mid = (int) floor( 0.5*n_yi );
	i_rm_mid = (int)floor(0.5*n_rm);
	
	i_s_mid = i2s_map[i_ph_mid][i_rent_mid][i_yi_mid][i_rm_mid];

    for( i_ph = 0; i_ph < n_ph; i_ph++ ){
		s_ph_midry[i_ph] = i2s_map[i_ph][i_rent_mid][i_yi_mid][i_rm_mid];
    }

	int i_t;

	for (i_t = 0; i_t < t_n; i_t++) {
		hu_ten[i_t] = hu_ten_store[city_id][i_t];
	}

	ten_w[0] = 0.0;

	for (i_t = 1; i_t < t_n; i_t++) {
		ten_w[i_t] = hu_ten[i_t] / hu_med[city_id];
	}
	rent_adj = hu_ten[0] / hu_med[city_id];

	 
	// work here
	vector<vector<double>> zeros_WN_RX(w_n, vector<double>(retxn, 0.0));
	vector<vector<vector<vector<double>>>> zeros_NS_NS_WN_RX(n_s, vector<vector<vector<double>>>(n_s, zeros_WN_RX));
	vector<vector<vector<vector<vector<double>>>>> zeros_T_NS_NS_WN_RX(T_max, zeros_NS_NS_WN_RX);

	w_t2_state = zeros_T_NS_NS_WN_RX; 
	// w_t2_state
	// given s1, s2, plots where wealth will be in next period (based on equity returns)
	// Given: (t_hor, i_s1, i_s2, i_w1), 
	// Last two arguments store w_low, w_high at beginning of next period

	// work here 
	int own_state0 = 1;  // assume household is homeowner by default
	vector<vector<int>> ones_int_NS_WN(n_s, vector<int>(w_n, own_state0));
	vector<vector<vector<int>>> ones_int_T_NS_WN(T_max, ones_int_NS_WN);

	own_state = ones_int_T_NS_WN;
	// own_state
	// Given: (t_hor, i_s1) and assume HH begins period as owner
	// i_w1 value = {0 (HH defaults), 1 (HH owns) }

	// in order to go further, will have to assume the wealth distribution
}


// This function is used to model wealth when household refinances or defaults
// i_w1_new: new wealth index once costs of refinace/default are accounted for
void snodes::w_state_swap( int i_w1_swap_in) {

	int i_w1_swap = i_w1_swap_in;
	int i_s2, i_x2;
	
	for (i_s2 = 0; i_s2 < n_s; i_s2++) {
		for (i_x2 = 0; i_x2 < retxn; i_x2++) {
			w_t2_state[t_hor][i_s1][i_s2][i_w1][i_x2] = w_t2_state[t_hor][i_s1][i_s2][i_w1_swap][i_x2];
		}
	}

	 
}

// Initialize state distribution, wealth distribution
void snodes::init_dist() {
	vector<double> zeros_NS(n_s, 0.0);
	vector<double> zeros_WN(w_n, 0.0); 
	sdist = zeros_NS;                                     // Initial state prob distribution
	wdist = zeros_WN;                                     // Initial Wealth distribution

	// Initial conditions
	sdist[floor(n_s / 2.0)] = 1.0;                                  // First Year: impose in median state
	double w_med = 0.4;                                             // impose median wealth = 40k
	int i_w_med = int(floor((w_med - w_min) / (w_max - w_min)));    // convert to wealth index
	wdist[i_w_med] = 1.0;                                           // impose mass of wealth is in median of w_n
}

// wtrans: given a vector representing the wealth distribution, calculate wealth distribution in next period
// Inputs
// 1) initial state distribution (if no input, impose in median state)
// 2) initial wealth distribution (if no input, all HH's at median wealth) 
// 3) state transition matrix
// Output
// 1) next period wealth distribution (used to calculate probability of default)
// 2) hazard rate (conditional upon owning a home this year, state dist, and wealth dist,
// probability of defaulting in next period)
// NOTES: 
// wdist2: prob(in state) * prob trans to next state * wealth in curr state * wealth in next state
// hazard rate = probability of being in current state * wealth dist in current state * pdef | state, wealth
void snodes::wtrans() {
	//vector<double> sdist;
	// Main inputs
	//vector<double> sdist(n_s, 0.0);                    // Initial state prob distribution
	//vector<double> wdist(w_n, 0.0);                    // Initial Wealth distribution
	vector<vector<double>> gammap = gammat[t_hor];     // State transition matrix

	// Main outputs
	vector<double> wdist2(w_n, 0.0);            // Next period wealth distribution
	vector<double> def_state_dist(n_s, 0.0);    // proportion of HH's who default in each state
	double hazard_rate = 0.0;                   // hazard rate

	// Indices
	int i_s1p, i_s2p, i_w1p; // current state, next state, current wealth
	int i_w2l, i_w2h;        // wealth in low state, high state (next period)
	double w2_l, w2_h;       // wealth in low state, high state (next period)
	
	// COMPUTE: wdist2: wealth distribution in the next period
	for (i_s1p = 0; i_s1p < n_s; i_s1p++) {                         // cycle through possible states
		for (i_s2p = 0; i_s2p < n_s; i_s2p++) {
			for (i_w1p = 0; i_w1p < w_n; i_w1p++) {                
				w2_l = w_t2_state[t_hor][i_s1p][i_s2p][i_w1p][0];    // load in wealth in next period given state
				w2_h = w_t2_state[t_hor][i_s1p][i_s2p][i_w1p][1];    // (low and high wealth realizations)
	
				i_w2l = round( (w2_l - w_min) / (w_max - w_min) );   // round wealth realization to closest wealth index
				i_w2l = min( max(i_w2l, 0) , w_n - 1);               // TODO: double check round fn

				i_w2h = round( (w2_h - w_min) / (w_max - w_min) );
				i_w2h = min( max(i_w2h, 0), w_n - 1);

				wdist2[i_w2l] = wdist2[i_w2l] + 0.5*sdist[i_s1p]*gammap[i_s1p][i_s2p]*wdist[i_w1p];  // add outcome probability mass
				wdist2[i_w2h] = wdist2[i_w2h] + 0.5*sdist[i_s1p]*gammap[i_s1p][i_s2p]*wdist[i_w1p];  // 0.5: equity return
				
			}
		}
	}

	// COMPUTE: default hazard rate
	for (i_s1p = 0; i_s1p < n_s; i_s1p++) {                       // cycle through current state
		for (i_w1p = 0; i_w1p < w_n; i_w1p++) {                   // cycle through wealth
			hazard_rate = hazard_rate + sdist[i_s1p] * wdist[i_w1p] * (1.0 - double(own_state[t_hor][i_s1p][i_w1p]));
		}
	}
}

void snodes::adj_tax() {

	// variable of interest: 
	// double yi_gridt[T_max + 1][n_yi]; 

	int i_t2 = 0;
	int i_y2 = 0;

	double tax_brack[] = {0.0, 0.4385, 1.0595, 1.61, 2.88 };
	double tax_rate[] = { 0.15, 0.28, 0.31, 0.36, 0.396 };

	double y_btax2 = 0.0;
	double y_tax_bill = 0.0;
	double y_atax2 = 0.0;
	double y_diff = 0.0;

	for (i_t2 = 0; i_t2 < (T_max + 1); i_t2++) {
		for (i_y2 = 0; i_y2 < n_yi; i_y2++){ 

			yi_gridt_btax[i_t2][i_y2] = yi_gridt[i_t2][i_y2];
			
			y_btax2 = yi_gridt[i_t2][i_y2];           // load pre-tax income
			y_tax_bill = 0.0;

			if (y_btax2 >= tax_brack[0]) {
				y_diff = min(y_btax2 - tax_brack[0], tax_brack[1] - tax_brack[0] );
				y_tax_bill = y_tax_bill +  tax_rate[0] * y_diff;
			}

			if (y_btax2 >= tax_brack[1]) {
				y_diff = min( y_btax2 - tax_brack[1], tax_brack[2] - tax_brack[1] );
				y_tax_bill = y_tax_bill + tax_rate[1] * y_diff;
			}

			if (y_btax2 >= tax_brack[2]) {
				y_diff = min(y_btax2 - tax_brack[2], tax_brack[3] - tax_brack[2]);
				y_tax_bill = y_tax_bill + tax_rate[2] * y_diff;
			}

			if (y_btax2 >= tax_brack[3]) {
				y_diff = min( y_btax2 - tax_brack[3], tax_brack[4] - tax_brack[3]);
				y_tax_bill = y_tax_bill + tax_rate[3] * y_diff;
			}

			if (y_btax2 >= tax_brack[4]) {
				y_diff = y_btax2 - tax_brack[4];
				y_tax_bill = y_tax_bill + tax_rate[4] * y_diff;
			}

			yi_gridt[i_t2][i_y2] = y_btax2 - y_tax_bill;

		}
	}

}
