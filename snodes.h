// class snodes.h
// contains nodes used to discretize home prices, rents, and income 
// as well as transition matrix                                                
// Copyright A. Michael Sharifi, 2016

#ifndef SNODES_H
#define SNODES_H

class snodes {
	int i_s, i_rent, i_ph, i_yi;
	int i_rm;
	
public:

	// compute i_s_mid;
	int i_ph_mid, i_rent_mid, i_yi_mid, i_rm_mid;
	double loan_bal1; 

	int city_id;
	double hu_ten[t_n]; //= hu_ten_store[N_cities][t_n];
	double ten_w[t_n]; 
	double rent_adj;
	double csfLevSn;

	int age0;
	int T_max;
	int i_s_mid;
	int t_hor;
	int s_ph_midry[n_ph];                               // states where home prices index from low to high, but rent and yi are always median


	// State Grids
	vector<vector<double>> p_gridt, rent_gridt, yi_gridt, yi_gridt_btax;
	vector<vector<double>> rm_gridt;

	vector<vector<double>> plevel_gridt;     // price level grid
	vector<vector<double>> urate_gridt;       // unemployment grid
	vector<vector<double>> fedfunds_gridt;     // fedfunds
	
	vector<vector<double>> fedfunds_store;   // fed funds value

	vector<vector<vector<double>>> gammat;              // state transition matrix for each time period

	//vector<vector<vector<int>>> i2s_map;                // maps individual dimensions to state
	vector<vector<vector<vector<int>>>> i2s_map;                // maps individual dimensions to state

	vector<vector<vector<vector<int>>>> i2s_map2;              // maps individual dimension to state (version2)

	//vector<int> s2i_ph;                                 // maps state to current home price
	vector<int> s2i_plevel;                             // maps state to current price level
	vector<int> s2i_urate;                              // maps state to current unemployment rate
	vector<int> s2i_fedfunds;                           // maps state to current fedfunds
	
	vector<int> s2i_ph;                                 // maps state to current price
	vector<int> s2i_rent;                               // maps state to current rent
	vector<int> s2i_yi;                                 // maps state to current income
	vector<int> s2i_rm;

	snodes(int age0_in, int T_max_in, int city_id_in );
	
	void adj_tax();

	// adding default variables here
	int i_s1; // current state
	int i_w1; // current wealth state

	// w_t2_state: wealth outcome matrix
	// (t_hor, i_s1, i_s2, w_i1, w2 {low, high} )
    // given above values, w_t2_state returns wealth in next period (either low or high)
	//vector<vector<vector<vector<vector<double>>>>> w_t2_state;

	vector<vector<vector<vector<double>>>> w_t2_state_low;
	vector<vector<vector<vector<double>>>> w_t2_state_high;

	
	vector<vector<vector<int>>> own_state;  // for each (t_hor, i_s1, w_i1 ) track whether HH defaults
	vector<vector<vector<int>>> def_state;  // default outcome matrix

	void w_state_swap(int i_w1_new_in);

	//void wtrans();  // compute wealth transition and hazard rat
	//void init_dist();

	vector<double> sdist, wdist;

};

#endif
