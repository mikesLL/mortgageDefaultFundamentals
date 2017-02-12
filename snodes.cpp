// Copyright A. Michael Sharifi, 2016
#include "headers.h"

snodes::snodes(int age0_in, int T_max_in, int city_id_in) {

	city_id = city_id_in;
	age0 = age0_in;
	T_max = T_max_in;

	//csfLevSn = csfLev_pidxw[param_id] * 1.0 / csfmarg_store[city_id];

	// initialize home price, price level, urate, and fedfunds grids
	p_gridt = vector<vector<double>>(T_max + 1, vector<double>(n_ph, 0.0));
	rent_gridt = vector<vector<double>>(T_max + 1, vector<double>(n_rent, 0.0));
	yi_gridt = vector<vector<double>>(T_max + 1, vector<double>(n_yi, 0.0));

	plevel_gridt = vector<vector<double>>(T_max + 1, vector<double>(n_plevel, 0.0));
	urate_gridt = vector<vector<double>>(T_max + 1, vector<double>(n_urate, 0.0));
	fedfunds_gridt = vector<vector<double>>(T_max + 1, vector<double>(n_fedfunds, 0.0));

	// fedfunds store: pass in t_hor, i_s
	// retrieve fedfunds rate
	fedfunds_store = vector<vector<double>>(T_max + 1, vector<double>(n_s, 0.0));            
	
	// initialize gammat state-state transition matrix
	vector<vector<vector<double>>> zeros_T_NS_NS(T_max + 1, vector<vector<double>>(n_s, vector<double>(n_s, 0.0)));
	gammat = zeros_T_NS_NS;

	// initialize i2s and s2i mappings
	vector<int> zeros_NS(n_s, 0);

	// states: home price, price level, urate, fed funds
	// [n_ph][n_plevel][n_urate][n_fedfunds]
	// i2s_map2 will capture information on all states
	int i_plevel, i_urate, i_fedfunds;
	
	i_s = 0;
	s2i_ph.resize(n_s);
	s2i_plevel.resize(n_s);
	s2i_urate.resize(n_s);
	s2i_fedfunds.resize(n_s);
	//s2i_yi.resize(n_s);

	i2s_map.resize(n_ph);

	for (i_ph = 0; i_ph < n_ph; i_ph++) {

		i2s_map[i_ph].resize(n_plevel);
		for (i_plevel = 0; i_plevel < n_plevel; i_plevel++) {

			i2s_map[i_ph][i_plevel].resize(n_urate);
			for (i_urate = 0; i_urate < n_urate; i_urate++) {

				i2s_map[i_ph][i_plevel][i_urate].resize(n_fedfunds);
				for (i_fedfunds = 0; i_fedfunds < n_fedfunds; i_fedfunds++){

					i2s_map[i_ph][i_plevel][i_urate][i_fedfunds] = i_s;          // initialize i2s map

					s2i_ph[i_s] = i_ph;                                           // initialize s2_i maps
					s2i_plevel[i_s] = i_plevel;                            // maps state to current price level
					s2i_urate[i_s] = i_urate;                              // maps state to current unemployment rate
					s2i_fedfunds[i_s] = i_fedfunds;                        // maps state to current fedfunds
					i_s++;
				}
			}
		}
	}


	// compute i_s_mid;
	int i_plevel_mid, i_urate_mid, i_fedfunds_mid;

	i_ph_mid = (int) floor( 0.5*n_ph );
	i_plevel_mid = (int)floor(0.5*n_plevel);
	i_urate_mid = (int)floor(0.5*n_urate);
	i_fedfunds_mid = (int)floor(0.5*n_fedfunds);

	i_s_mid = i2s_map[i_ph_mid][i_plevel_mid][i_urate_mid][i_fedfunds_mid];

    for( i_ph = 0; i_ph < n_ph; i_ph++ ){
		s_ph_midry[i_ph] = i2s_map[i_ph][i_plevel_mid][i_urate_mid][i_fedfunds_mid];
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
	int i1, i2, i3, i4, i5;
	w_t2_state_low.resize(T_max);
	w_t2_state_high.resize(T_max);

	for (i1 = 0; i1 < T_max; i1++) {
		w_t2_state_low[i1].resize(n_s);
		w_t2_state_high[i1].resize(n_s);

		for (i2 = 0; i2 < n_s; i2++) {
			w_t2_state_low[i1][i2].resize(n_s);
			w_t2_state_high[i1][i2].resize(n_s);

			for (i3 = 0; i3 < n_s; i3++) {
				w_t2_state_low[i1][i2][i3].resize(w_n);
				w_t2_state_high[i1][i2][i3].resize(w_n);

				w_t2_state_low[i1][i2][i3].assign(w_n, 0.0);
				w_t2_state_high[i1][i2][i3].assign(w_n, 0.0);
			}
		}
	}
	
	// work here 
	int own_state0 = 1, def_state0 = 0;  // assume household is homeowner by default
	vector<vector<int>> ones_int_NS_WN(n_s, vector<int>(w_n, own_state0));
	vector<vector<int>> zeros_int_NS_WN(n_s, vector<int>(w_n, def_state0));

	vector<vector<vector<int>>> ones_int_T_NS_WN(T_max, ones_int_NS_WN);
	vector<vector<vector<int>>> zeros_int_T_NS_WN(T_max, zeros_int_NS_WN);

	own_state = ones_int_T_NS_WN;
	def_state = zeros_int_T_NS_WN;

	double def_ltv_state0 = 0.0;    // default ltv state
	def_ltv_state.resize(T_max);
	for (i1 = 0; i1 < T_max; i1++) {
		def_ltv_state[i1] = vector<vector<double>>(n_s, vector<double>(w_n, 0.0)); 
	}



}

// own_state
// Given: (t_hor, i_s1) and assume HH begins period as owner
// i_w1 value = {0 (HH defaults), 1 (HH owns) }

// in order to go further, will have to assume the wealth distribution

// This function is used to model wealth when household refinances or defaults
// i_w1_new: new wealth index once costs of refinace/default are accounted for
void snodes::w_state_swap( int i_w1_swap_in) {

	int i_w1_swap = i_w1_swap_in;
	int i_s2;
	
	for (i_s2 = 0; i_s2 < n_s; i_s2++) {
		w_t2_state_low[t_hor][i_s1][i_s2][i_w1] = w_t2_state_low[t_hor][i_s1][i_s2][i_w1_swap];
		w_t2_state_high[t_hor][i_s1][i_s2][i_w1] = w_t2_state_high[t_hor][i_s1][i_s2][i_w1_swap];
	}
}

//for (i_x2 = 0; i_x2 < retxn; i_x2++) {
//	w_t2_state[t_hor][i_s1][i_s2][i_w1][i_x2] = w_t2_state[t_hor][i_s1][i_s2][i_w1_swap][i_x2];
//}

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
