/*
ufnEV2.h
Class ufnEV2 contains data related to value function evaluation and stores data in a
way which speeds up computing home price and value function expecations
Copyright A. Michael Sharifi, 2016
*/

using namespace std;


#ifndef UFNEV2_H
#define UFNEV2_H

class vfn;

class ufnEV2 {
	eval_res res2;
	int w_i_low;
	double w_i_d;
	double alpha1;
	double lambda1;
	double v_tlower;
	double v_tupper;
	double w_low, w_high;
	double csfLevSn;
	double loan_bal1;

	int N_s2p;           // number of states with positive probability
	int i_s2p;            // state positive probability index
	int i_s2p_vec[n_s];       // store positive probability indices

	vector<vector<double>> vw3_grid_ti2;

	vector<vector<double>> vw3_grid_move;

	vector<vector<double>> vw3_d_grid_ti2;
	vector<vector<double>> vw3_dd_grid_ti2;

	int i_s2, i_ph2, i_x2;                         // state index, price index, equity return index
	double rb_eff = rb;                            // effective return on bonds	(default = rb)       
	double Evw_2 = 0.0;                            // Value Function Expectation
	double w2, w2_move, vw2;

	double uc; // = ufn(x[0], hu, (*vf2).pref);  // composite utility
	double rb_eff_agg;                       // aggregated gross return on bonds

	double b_unsec; // = 0.0;
	double b_sec; // = 0.0;
	double rb_unsec; //  = rb + credit_prem;
	double rb_sec;

	int i_urate1, i_fedfunds1, i_plevel1, i_yi1;

public:
	int i_s1;             // current state
	int t_hor;

	vfn *vf2;              // initialize next period structure 
	snodes *snodes1;       // state space

	int t_i2;              // tenure in next period
	double hu;             // housing utility
	double ph_2e, ph1;     // home price expectation, current home price


	double mort_add_equity; // MOD HERE: increasing wealth by home equity
	// equity increases by this amount

	eval_res res1;
	eval_res res1_move;

	double csf_net2[n_s];   // home price futures payout

	void enter_data(void *snodes_in, void *vf2_in); // enter data into utilty fn 
	double eval(vector<double> x_in);               // evaluate utility fn

	void store_wlh2(vector<double> x_in);               // evaluate utility fn

	inline eval_res eval_v( int i_s_in, double w_in);
	inline eval_res eval_v_move(int i_s_in, double w_in);
};

#endif

