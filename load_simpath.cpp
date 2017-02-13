/*
load_simpath.cpp
given intial parameters, this program simulates a home price path
and returns interpolating prices as well as a transition matrix

Copyright A. Michael Sharifi, 2016

*ps0: pointer to price structure
*gs0: pointer to transition matrix stucture
rent_in: initial rent
ph0_in: initial home price
ret0_in: most recently observed hpa
csf_1_yr_in: currently observed futures price
t_id: initial year
city_init_in: initial city

impose forecasted 1st year return equals that implied by csf
ph0 : current price
csf1_yr: forecasted 1yr price
ret_tn = (csf_1yr - ph0) / ph0 + sigma_ret*randn_draw;

In actual regression:
1. Estimate cointegrating relation ship with rent to price ratio
2. regress nominal returns on lagged nominal returns
3. Done; parameters can be entered here

Numeraire: non-durable consumption
Rent grows in real terms
Income grows in real terms as the agent ages and city-wide income increases
*/

#include "headers.h"

void load_simpath(void *snodes_in, int grent_id_in, double grent_in, double rent_in, double ph0_in, double ret0_in, double csf_1yr_in, int t_id, string city_init_in, int city_id, int age_begin_in) {

	snodes *snodes1 = (snodes*)snodes_in;

	int s1, s2, s_test;                             // state in current period, state in next period

	// MOD HERE: limiting number of simulations for now 
	int N_print =  1000;                          // number of observations to print to file
	int N_sim = 1000000;                                            // number of simulations

	int i_ph, i_rent, i_yi, i_rm, i_s;                                      // state and individual dimension indices

	cout << "load ppath: ph0_in:  " << ph0_in << endl;
	
	city_id = 1;

	double inf_mult = 104.5; // inflation multiplier

	//determnistic (real) rent growth
	double g_rent = grent_in; // 0.0; //0.02;

	ph0_in = ph0_in / inf_mult * 100.0;
	rent_in = rent_in / inf_mult * 100.0;
	csf_1yr_in = csf_1yr_in / inf_mult * 100.0; 

	double y_city_mult[] = { 1.23, 1.54, 1.02, 1.34, 1.02, 1.00,  0.95, 1.39 };

	double rho_rm = 0.4;    // short-rate: autocorrelation
	double theta_rm = 1.0;  // short-rate: mean-reversion
	double rm_mu = 0.04; // long-run short rate mean

	double alpha_hat = 0.02;
	double rhof_hat = 0.46;
	double theta_hat = 0.34;
	double sigma_ret = 0.08;
	double gamma0_hat = -2.27;
	double gamma1_hat = 0.59;

	int grent_id = grent_id_in;   // set = 0 for low rent growth, set = 1 fo high rent growth

	// set rent growth parameters here
	if (grent_id) {
		alpha_hat = 0.0052;
		rhof_hat = 0.5780;
		theta_hat = 0.3526;
		sigma_ret = 0.0347;
		gamma0_hat = -2.27;
		gamma1_hat = 0.59;
		g_rent = 0.0029;
	}
	else {
		alpha_hat = 0.0053;
		rhof_hat = 0.6885;
		theta_hat = 0.2990;
		sigma_ret = 0.0390;
		gamma0_hat = -2.27;
		gamma1_hat = 0.59;
		g_rent = -0.0106;
	}

	/* Parameters from MATLAB (real terms)
	rho0_est_high: 0.5780
	rho0_est_low : 0.6885
	alpha0_est_high : 0.0052
	alpha0_est_low : 0.0053
	theta0_est_high : 0.3526
	theta0_est_low : 0.2990
	sigma1_est_high : 0.0347
	sigma1_est_low : 0.0390
	rp_avg_high : 0.0574
	rp_avg_low : 0.0551
	ret_avg_high : 0.0062
	ret_avg_low : 0.0047

	g_rent_avg_low : -0.0106
	g_rent_avg_high : 0.0029
	*/

	 
	double p_min = 0.01;                                                         // lower bound on home prices	

	double rm_str0 = 0.0124;                                                      // initial short rate
	//double mort_rate = rm_str0 + sr_mort_prem;                                 // initial mortgage rate
	
	// stdev centers (tauchen-style discretization) 
    double rm_nd_std[] = { -2.0, -1.0, 0.0, 1.0, 2.0 };
	//double ph_nd_std[] = { -2.0, -1.0, 0.0, 1.0, 2.0 };     
	double ph_nd_std[] = { -2.0, -1.5, -1.0, -0.5, 0.0, .5, 1.0, 1.5, 2.0 };     
	double rent_nd_std[] = { 0.0 };                                            // can also set { -1.0, 0.0, 1.0 };
	
	double plevel_nd_std[] = { -1.0, 0.0, 1.0 };            // store price-level std nodes
	double urate_nd_std[] = { -1.0, 0.0, 1.0 };             // store unemployment rate std nodes
	double fedfunds_nd_std[] = { 0.0 };                     // fed funds mean income std nodes
	double yi_nd_std[] = { 0.0, 0.0, 0.0 };                 // store mean income 
	
	string city_init = city_init_in;                                          // initial conditions for simulation
	    
	double rent0 = rent_in;

	cout << "load ppath; rent = " << rent0 << endl;
	
	double ph0 = ph0_in;
	double ret0 = ret0_in;
	
	int T_max = (*snodes1).T_max;

	cout << "begin load_simpath" << endl ;
    double duration;
	clock_t start;                          // instantiate timer
	start = clock();
	
	random_device rd;                       // instantiate random generator
	mt19937 gen(rd());
	normal_distribution<> dist(0.0, 1.0);   // standard normal distribution

	int i = 0, j = 0;
	cout << " N_sim = " << N_sim << endl;
	
	int T_sim = T_max + 1;                                                              // T_max is a measure of the investor's horizon
	
	cout << "T_sim = " << T_sim << endl;
	
	double ret_tn, ret_lag, ecm;

	// store observations
	vector<vector<double>> plevel_str(T_sim + 1, vector<double>(N_sim, 0.0));               // prive level
	vector<vector<double>> urate_str(T_sim + 1, vector<double>(N_sim, 0.0));               // unemployment rate
	vector<vector<double>> fedfunds_str(T_sim + 1, vector<double>(N_sim, 0.0));             // fed-funds rate	
	vector<vector<double>> ph_str(T_sim + 1, vector<double>(N_sim, 0.0) );              // store city-wide home prices
	vector<vector<double>> rent_str(T_sim + 1, vector<double>(N_sim, 0.0));             // store city-wide rent

	// NODES 
	vector<vector<double>> plevel_str_nds(T_sim + 1, vector<double>(n_plevel, 0.0));           // price-level nodes
	vector<vector<double>> urate_str_nds(T_sim + 1, vector<double>(n_urate, 0.0));           // unemployment rate nodes
	vector<vector<double>> fedfunds_str_nds(T_sim + 1, vector<double>(n_fedfunds, 0.0));           // fed funds nodes
	vector<vector<double>> ph_str_nds(T_sim + 1, vector<double> (n_ph, 0.0) );         // home price nds
	vector<vector<double>> rent_str_nds(T_sim + 1, vector<double> (n_rent, 0.0) );     // rent nds

	double z1, z2, z3;
	double mult_2005 =   .56 / .51; // 1.5;                              // multiplier: convert $1992 to $2005
	double mult_unit = 0.01;                             // multiplier: convert $1000s from parameters to $100000s in paper


	int t = 0, n = 0, n1 = 0, n2 = 0;                                          // time and simulation indices
	double res_sum, res_sq_sum, res_mean, res_std, row_sum;                    // sums for computing means and standard deviations

 
	double eps_h; // housing shock

	// Macro initial conditions
	double pinf0 = 0.3;
	double urate0 = 0.05;
	double fedfunds0 = 0.01; 
	double y_inc0 = 0.8; // TODO: let be a fn of MTI
	double g_y = 0.01; // Real income growth

	double var_a[] = { 0.0025, 0.0395, 0.0474 };                       // VAR: constants

	double var_b[3][3] = { {-0.2505,    0.2879,    0.3663},            // VAR: coefficients
	                       {0.0528,    0.4644, -0.1822},
						   {-0.5300, -0.4020, 0.7285} };
	
	//double var_b[3][3] = { { 1.0, 0.0, 0.0 },
	//						{ 1.0, 1.0, 1.0 },
	//						{ 0.0, 0.0, 1.0 } };

	double plevel0 = 1.0;
	double cholQ[2][2] = { { 0.0099, 0.0}, { -0.0020,  0.0104 } };

	double v0[] = { pinf0, urate0, fedfunds0 };             // inflation rate unemployment rate
	double v1[] = { 0.0, 0.0, 0.0 };
	double v2[] = { 0.0, 0.0, 0.0 };

	double vz1, vz2, vz3;
	double vu1, vu2, vu3;

	double pinf, pinf_lag;
	
	// loop through simulations
	cout << "load_simpath.cpp: Begin Simulations" << endl;
	for (n = 0; n < N_sim; n++){                  
		if ( n % 10000 == 0) {
			cout << "n = " << n << endl;
		}
		
		// Initial year
		t = 0;            

		pinf_lag = pinf0;
		
		plevel_str[t][n] = 1.0;             // original price level
		urate_str[t][n] = urate0;           // unemployment rate
		fedfunds_str[t][n] = fedfunds0;     // fed-funds rate

		ph_str[t][n] = log(ph0);
		rent_str[t][n] = rent0;

		// simulate later time periods
		for (t = 1; t < (T_sim + 1); t++){
			
			z1 = dist(gen);                    // draw shocks
			z2 = dist(gen);
			z3 = dist(gen);

			vz1 = dist(gen);                   // VAR shocks (independent)
			vz2 = dist(gen);
			vz3 = dist(gen);

			vu1 = cholQ[0][0] * vz1;          //  VAR shocks (reduced form)
			vu2 = cholQ[1][0] * vz1 + cholQ[1][1] * vz2;
 
			eps_h = sigma_ret*z1;                // home price innovation
			
			// compute inflation, unemp, fedfunds using the VAR + shocks
			v0[0] = pinf_lag; // pinf_str[t - 1][n];
			v0[1] = urate_str[t - 1][n];
			v0[2] = fedfunds_str[t - 1][n];
			matrix_mult(var_a, var_b, v0, v1, 1);
			pinf_lag = v1[0] + vu1;
			plevel_str[t][n] = plevel_str[t - 1][n] * (1.0 + pinf_lag);
			urate_str[t][n] = v1[1] + vu2;
			fedfunds_str[t][n] = v1[2];
			
			//matrix_mult(var_a, var_b, v0, v2, 10);          // compute the mortgage rate (if were to add FRM)
			//rm_str[t][n] = max(v2[2] + 0.03, 0.0);                  // Mortgage premium!
		
			// compute home price
			ret_lag = ph_str[t - 1][n] - ph_str[max(t - 2, 0)][n];             // ph_str is in logs 
			
			ecm = log(rent_str[t - 1][n]) - gamma0_hat - gamma1_hat*(ph_str[t - 1][n]);          // cointegrate rents, prices

			// cointegrate interest rates, rents, and prices
			ret_tn = alpha_hat + rhof_hat*ret_lag + theta_hat*ecm + eps_h;         // return series
			
			rent_str[t][n] = exp(g_rent + pinf_lag)*rent_str[t-1][n];               // update rent, price
			ph_str[t][n] = ret_tn + ph_str[t - 1][n] + pinf_lag;
			
		}
	}

	cout << "load_simpath.cpp: Finished Running Simulations" << endl;

	cout << "load_simpath.cpp: Compute Price-Level Nodes " << endl;
	for (t = 0; t < T_sim; t++) {
		res_sum = accumulate(plevel_str[t].begin(), plevel_str[t].end(), 0.0);
		res_mean = res_sum / (double)plevel_str[t].size();
		res_sq_sum = inner_product(plevel_str[t].begin(), plevel_str[t].end(), plevel_str[t].begin(), 0.0);
		res_std = sqrt(res_sq_sum / (double) plevel_str[t].size() - res_mean * res_mean);

		if (t == 0) {
			res_std = 0.01;  // initial period: manually set stdev 
		}

		for (n = 0; n < n_plevel; n++) {
			plevel_str_nds[t][n] = res_mean + plevel_nd_std[n] * res_std;   // compute mortgage rate nodes 
			plevel_str_nds[t][n] = max(plevel_str_nds[t][n], 0.0);          // set non-negative short-rate
			(*snodes1).plevel_gridt[t][n] = plevel_str_nds[t][n];

			cout << (*snodes1).plevel_gridt[t][n] << "...";              // print out interest rate grid entry
		}
		cout << endl;
	}

	cout << "load_simpath.cpp: Compute Unemployment Rates Nodes " << endl;
	for (t = 0; t < T_sim; t++) {
		res_sum = accumulate(urate_str[t].begin(), urate_str[t].end(), 0.0);
		res_mean = res_sum / (double) urate_str[t].size();
		res_sq_sum = inner_product(urate_str[t].begin(), urate_str[t].end(), urate_str[t].begin(), 0.0);
		res_std = sqrt(res_sq_sum / (double)urate_str[t].size() - res_mean * res_mean);

		if (t == 0) {
			res_std = 0.01;  // initial period: manually set stdev 
		}

		for (n = 0; n < n_urate; n++) {
			urate_str_nds[t][n] = res_mean + urate_nd_std[n] * res_std;   // compute mortgage rate nodes 
			urate_str_nds[t][n] = max(urate_str_nds[t][n], 0.0);          // set non-negative short-rate
			(*snodes1).urate_gridt[t][n] = urate_str_nds[t][n];

			cout << (*snodes1).urate_gridt[t][n] << "...";              // print out interest rate grid entry
		}
		cout << endl;
	}

	cout << "load_simpath.cpp: Compute Fedfunds Nodes " << endl;
	for (t = 0; t < T_sim; t++) {
		res_sum = accumulate(fedfunds_str[t].begin(), fedfunds_str[t].end(), 0.0);
		res_mean = res_sum / (double)fedfunds_str[t].size();
		res_sq_sum = inner_product(fedfunds_str[t].begin(), fedfunds_str[t].end(), fedfunds_str[t].begin(), 0.0);
		res_std = sqrt(res_sq_sum / (double)fedfunds_str[t].size() - res_mean * res_mean);

		if (t >= 0) {
			res_std = 0.01;  // initial period: manually set stdev 
		}

		for (n = 0; n < n_fedfunds; n++) {
			fedfunds_str_nds[t][n] = res_mean + fedfunds_nd_std[n] * res_std;   // compute mortgage rate nodes 
			fedfunds_str_nds[t][n] = max(fedfunds_str_nds[t][n], 0.0);          // set non-negative short-rate
			(*snodes1).fedfunds_gridt[t][n] = fedfunds_str_nds[t][n];

			cout << (*snodes1).fedfunds_gridt[t][n] << "...";              // print out interest rate grid entry
		}
		cout << endl;
	}

	
	cout << "load_simpath.cpp: Compute Home Price Nodes " << endl;
	for (t = 0; t < T_sim; t++) {
		res_sum = accumulate(ph_str[t].begin(), ph_str[t].end(), 0.0);
		res_mean = res_sum / (double)ph_str[t].size();
		res_sq_sum = inner_product(ph_str[t].begin(), ph_str[t].end(), ph_str[t].begin(), 0.0);
		res_std = sqrt(res_sq_sum / (double) ph_str[t].size() - res_mean * res_mean);

		if (t == 0) {
			res_std = 0.01;  // initial period: manually set stdev 
		}

		for (n = 0; n < n_ph; n++) {
			ph_str_nds[t][n] = res_mean + ph_nd_std[n] * res_std;  // compute home price nodes (log scale)
			(*snodes1).p_gridt[t][n] = exp(ph_str_nds[t][n]);    // load home price nodes into p_gridt
			cout << (*snodes1).p_gridt[t][n] << "...";      // print out price grid entry
		}
		cout << endl;    
	}

	
	cout << "load_simpath.cpp: Compute Rent Nodes " << endl;
	for (t = 0; t < T_sim; t++) {
		res_sum = accumulate(rent_str[t].begin(), rent_str[t].end(), 0.0);
		res_mean = res_sum / (double) rent_str[t].size();

		res_sq_sum = inner_product(rent_str[t].begin(), rent_str[t].end(), rent_str[t].begin(), 0.0);
		res_std = sqrt(res_sq_sum / (double) rent_str[t].size() - res_mean * res_mean);

		// check null stdev case
		if ( (res_std != res_std ) || (res_std <= 0.01) ) {
			res_std = 0.01;  // initial period: manually set stdev
		}

		for (n = 0; n < n_rent; n++) {
			rent_str_nds[t][n] = res_mean + rent_nd_std[n] * res_std;
			(*snodes1).rent_gridt[t][n] = rent_str_nds[t][n];    // load rent nodes into rent_gridt
		}
	}


	cout << "load_simpath.cpp: Compute Y income (Median City-wide income) Nodes " << endl;
	for (t = 0; t < T_sim; t++) {
		for (n = 0; n < n_yi; n++) {
			(*snodes1).yi_gridt[t][n] = y_inc0 * pow( 1.0 + g_y, t );    // load rent nodes into rent_gridt
		}
	}
	
	
	cout << "load_simpath.cpp: Begin Computing State-Space Transition Matrix" << endl;

	vector<vector<int>> zeros_NS_NS(n_s, vector<int>(n_s, 0)); // default matrix for state to state transition
	vector<vector<int>> gamma_store = zeros_NS_NS;             // gamma matrix tracks state to state transition

	t = 0;
	double ph_step;
	double rent_step;
	double yi_step;
	double rm_step;

	double plevel_step;
	double urate_step;
	double fedfunds_step;

	double ph_step1;
	double rent_step1;
	double yi_step1;
	double rm_step1;

	double plevel_step1;
	double urate_step1;
	double fedfunds_step1;

	int i_plevel, i_urate, i_fedfunds;

	vector<double> fedfunds_store_sum(n_s, 0.0);
	vector<int> fedfunds_store_count(n_s, 0);

	cout << "load_simpath.cpp: Cycle Through Observations" << endl;
	for (t = 0; t < T_sim; t++) {
		fedfunds_store_sum = vector<double>(n_s, 0.0);
		fedfunds_store_count = vector<int>(n_s, 0);

		gamma_store = zeros_NS_NS;                                  // initialize current-period transition matrix

		ph_step = ph_str_nds[t][1] - ph_str_nds[t][0];                         // step-sizes
		plevel_step = plevel_str_nds[t][1] - plevel_str_nds[t][0];             // step-sizes
		urate_step = urate_str_nds[t][1] - urate_str_nds[t][0];                // step-sizes
		fedfunds_step = 0.001; // fedfunds_str_nds[t][1] - fedfunds_str_nds[t][0];    

		ph_step1 = ph_str_nds[t + 1][1] - ph_str_nds[t + 1][0];
		plevel_step1 = plevel_str_nds[t+1][1] - plevel_str_nds[t+1][0];          // step-sizes
		urate_step1 = urate_str_nds[t+1][1] - urate_str_nds[t+1][0];             // step-sizes
		fedfunds_step1 = 0.001; // fedfunds_str_nds[t + 1][1] - fedfunds_str_nds[t + 1][0];    

		for (n = 0; n < N_sim; n++) {

			// current period: map observation to closest node
			i_ph = (int) round( (ph_str[t][n] - ph_str_nds[t][0]) / ph_step);
			i_plevel = (int) round( (plevel_str[t][n] - plevel_str_nds[t][0]) / plevel_step);
			i_urate = (int) round( (urate_str[t][n] - urate_str_nds[t][0]) / urate_step);
			i_fedfunds = (int) round( (fedfunds_str[t][n] - fedfunds_str_nds[t][0]) / fedfunds_step);
			
			i_ph = min(max(i_ph, 0), n_ph - 1);              // bound extremes
			i_plevel = min(max(i_plevel, 0), n_plevel - 1);              // bound extremes
			i_urate = min(max(i_urate, 0), n_urate - 1);              // bound extremes
			i_fedfunds = min(max(i_fedfunds, 0), n_fedfunds - 1);              // bound extremes
			 
			s1 = (*snodes1).i2s_map[i_ph][i_plevel][i_urate][i_fedfunds];             // map node to state

			// fedfunds: store observation
			fedfunds_store_sum[s1] += fedfunds_str[t][n];
			fedfunds_store_count[s1] += 1;
			
			// next period: map observation to closest node
			i_ph = (int)round((ph_str[t + 1][n] - ph_str_nds[t + 1][0]) / ph_step1);
			i_plevel = (int)round((plevel_str[t+1][n] - plevel_str_nds[t+1][0]) / plevel_step1);
			i_urate = (int)round((urate_str[t+1][n] - urate_str_nds[t+1][0]) / urate_step1);
			i_fedfunds = (int)round((fedfunds_str[t+1][n] - fedfunds_str_nds[t+1][0]) / fedfunds_step1);

			i_ph = min(max(i_ph, 0), n_ph - 1);              // bound extremes
			i_plevel = min(max(i_plevel, 0), n_plevel - 1);              // bound extremes
			i_urate = min(max(i_urate, 0), n_urate - 1);              // bound extremes
			i_fedfunds = min(max(i_fedfunds, 0), n_fedfunds - 1);              // bound extremes

			s2 = (*snodes1).i2s_map[i_ph][i_plevel][i_urate][i_fedfunds];          // map node to state
			
			gamma_store[s1][s2] += 1;                                // store  in transition matrix
		}

		cout << "load_simpath.cpp: t = " << t  << "  load transition matrix into data structure "<< endl;

		// load transition matrix into data structure
		for (i = 0; i < n_s; i++) {
			row_sum = (double) accumulate(gamma_store[i].begin(), gamma_store[i].end(), 0);

			for (j = 0; j < n_s; j++) {
				if (row_sum >= 1.0) {
					(*snodes1).gammat[t][i][j] = (double)gamma_store[i][j] / row_sum;   // load transition matrix into structure
				} else {
					//(*snodes1).gammat[t][i][j] = 0.0;
					(*snodes1).gammat[t][i][j] = 1.0 / double(n_s); // default setting if no obs in space
				}
			}
		}

		// load fedfunds into data structure
		for (i = 0; i < n_s; i++) {
			fedfunds_store_count[i] = max(fedfunds_store_count[i], 1);
			(*snodes1).fedfunds_store[t][i] = max( fedfunds_store_sum[i] / (double)fedfunds_store_count[i], 0.0 );
		}

	}

	cout << "load_simpath.cpp: Finished Computing Transition Matrix" << endl;

	cout << "load_simpath.cpp: gamma T max:" << endl;
	
	string fn_beg = "price_results/" + city_init + "_age" + to_string((*snodes1).age0) + "_yr" + to_string(t_id);

	// PRINT: HOME PRICE, PLEVEL, URATE< FEDFUNDS PATHS
	// print sample HOME PRICE paths
	ofstream price_obs_file;                                        
	price_obs_file.open(fn_beg + "_price_obs_file.csv", ios::out | ios::trunc);

	for (n = 0; n < N_print; n++) {
		for (t = 0; t < T_sim; t++) {
			price_obs_file << exp(ph_str[t][n]) << ",";
		}
		price_obs_file << endl;
	}
	price_obs_file.close();

	// print sample PLEVEL paths
	ofstream plevel_obs_file;                                     
	plevel_obs_file.open( fn_beg +  "_plevel_obs_file.csv", ios::out | ios::trunc);
	for (n = 0; n < N_print; n++) {
		for (t = 0; t < T_sim; t++) {
			plevel_obs_file << plevel_str[t][n] << ",";
		}
		plevel_obs_file << endl;
	}
	plevel_obs_file.close();
	
	// print sample URATE paths
	ofstream urate_obs_file;
	urate_obs_file.open(fn_beg + "_urate_obs_file.csv", ios::out | ios::trunc);
	for (n = 0; n < N_print; n++) {
		for (t = 0; t < T_sim; t++) {
			urate_obs_file << urate_str[t][n] << ",";
		}
		urate_obs_file << endl;
	}
	urate_obs_file.close();

	// print sample fedfunds paths
	ofstream fedfunds_obs_file;
	fedfunds_obs_file.open(fn_beg + "_fedfunds_obs_file.csv", ios::out | ios::trunc);
	for (n = 0; n < N_print; n++) {
		for (t = 0; t < T_sim; t++) {
			fedfunds_obs_file << fedfunds_str[t][n] << ",";
		}
		fedfunds_obs_file << endl;
	}
	fedfunds_obs_file.close();


	// print out HOME PRICE grid
	ofstream pstruct_file;
	pstruct_file.open(fn_beg + "_pstruct_file.csv", ios::out | ios::trunc);
	for (t = 0; t < T_sim; t++) {
		for (n = 0; n < n_ph; n++) {
			pstruct_file << (*snodes1).p_gridt[t][n] << ",";         // alt: ph_str_nds[t][n];
		}
		pstruct_file << endl;
	}

	// print out PLEVEL grid
	ofstream plevel_struct_file;
	plevel_struct_file.open(fn_beg + "_plevel_struct_file.csv", ios::out | ios::trunc);
	for (t = 0; t < T_sim; t++) {
		for (n = 0; n < n_plevel; n++) {
			plevel_struct_file << (*snodes1).plevel_gridt[t][n] << ",";         // alt: ph_str_nds[t][n];
		}
		plevel_struct_file << endl;
	}

	// print out URATE grid
	ofstream urate_struct_file;
	urate_struct_file.open(fn_beg + "_urate_struct_file.csv", ios::out | ios::trunc);
	for (t = 0; t < T_sim; t++) {
		for (n = 0; n < n_urate; n++) {
			urate_struct_file << (*snodes1).urate_gridt[t][n] << ",";         // alt: ph_str_nds[t][n];
		}
		urate_struct_file << endl;
	}

	// print out FEDFUNDS grid
	ofstream fedfunds_struct_file;
	fedfunds_struct_file.open(fn_beg + "_fedfunds_struct_file.csv", ios::out | ios::trunc);
	for (t = 0; t < T_sim; t++) {
		for (n = 0; n < n_fedfunds; n++) {
			fedfunds_struct_file << (*snodes1).fedfunds_gridt[t][n] << ",";         // alt: ph_str_nds[t][n];
		}
		fedfunds_struct_file << endl;
	}


	// print out rent grid
	ofstream rstruct_file;
	rstruct_file.open(fn_beg + "_rstruct_file.csv", ios::out | ios::trunc);
	for (t = 0; t < T_sim; t++) {
		for (n = 0; n < n_rent; n++) {
			rstruct_file << (*snodes1).rent_gridt[t][n] << ",";         // alt: ph_str_nds[t][n];
		}
		rstruct_file << endl;
	}

	// print sample rent paths
	ofstream rent_obs_file;
	rent_obs_file.open(fn_beg + "_rent_obs_file.csv", ios::out | ios::trunc);
	for (n = 0; n < N_print; n++) {
		for (t = 0; t < T_sim; t++) {
			rent_obs_file << rent_str[t][n] << ",";
		}
		rent_obs_file << endl;
	}
	rent_obs_file.close();

	// print out transition matrix
	ofstream gstruct_file;                                                           
	gstruct_file.open(fn_beg + "_gstruct_file.csv", ios::out | ios::trunc);

	for (t = 0; t < T_sim; t++) {
		for (i = 0; i < n_s; i++) {
			for (j = 0; j < n_s; j++) {
				gstruct_file << (*snodes1).gammat[t][i][j] << ",";
			}
			gstruct_file << endl;
		}
		gstruct_file << endl;
	}

	duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	cout<<"Finished load_rand. Time elapsed: "<< duration <<'\n';


}


// print out income grid
/*
ofstream ystruct_file;
ystruct_file.open( fn_beg + "_ystruct_file.csv", ios::out | ios::trunc);
for (t = 0; t < T_sim; t++) {
for (n = 0; n < n_yi; n++) {
ystruct_file << (*snodes1).yi_gridt[t][n] << ",";         // alt: ph_str_nds[t][n];
}
ystruct_file << endl;
}*/

/*
// print sample income paths
ofstream yi_obs_file;
yi_obs_file.open( fn_beg +  "_yi_obs_file.csv", ios::out | ios::trunc);
for (n = 0; n < N_print; n++) {
for (t = 0; t < T_sim; t++) {
yi_obs_file << exp(yi_str[t][n]) << ",";
}
yi_obs_file << endl;
}
yi_obs_file.close();
*/


/*
cout << "load_simpath.cpp: Compute Rate Nodes " << endl;
for (t = 0; t < T_sim; t++) {
res_sum = accumulate(rm_str[t].begin(), rm_str[t].end(), 0.0);
res_mean = res_sum / (double)rm_str[t].size();
res_sq_sum = inner_product(rm_str[t].begin(), rm_str[t].end(), rm_str[t].begin(), 0.0);
res_std = sqrt(res_sq_sum / (double)rm_str[t].size() - res_mean * res_mean);

if (t == 0) {
res_std = 0.01;  // initial period: manually set stdev
}

for (n = 0; n < n_rm; n++) {
rm_str_nds[t][n] = res_mean + rm_nd_std[n] * res_std;   // compute mortgage rate nodes
rm_str_nds[t][n] = max(rm_str_nds[t][n], 0.0);          // set non-negative short-rate
(*snodes1).rm_gridt[t][n] = rm_str_nds[t][n];

cout << (*snodes1).rm_gridt[t][n] << "...";              // print out interest rate grid entry
}
cout << endl;
}
*/

/*
cout << "load_simpath.cpp: Compute Income Nodes " << endl;
for (t = 0; t < T_sim; t++) {
res_sum = accumulate(yi_str[t].begin(), yi_str[t].end(), 0.0);
res_mean = res_sum / (double)yi_str[t].size();

res_sq_sum = inner_product(yi_str[t].begin(), yi_str[t].end(), yi_str[t].begin(), 0.0);
res_std = sqrt(res_sq_sum / (double)yi_str[t].size() - res_mean * res_mean);

if (t == 0) {
res_std = 0.01;
}
for (n = 0; n < n_yi; n++) {
yi_str_nds[t][n] = res_mean + yi_nd_std[n] * res_std; // compute income nodes (log scale)
//yi_str_nds[t][n] = res_mean + yi_nd_std[n] * res_std; // compute income nodes (log scale)
(*snodes1).yi_gridt[t][n] = exp(yi_str_nds[t][n]);    // load individual-income nodes into income_gridt
cout << t << "..." << n << "..." << (*snodes1).yi_gridt[t][n] << "...";

//if (t <= 5) {
//	(*snodes1).yi_gridt[t][n] = 0.0;
//}
}
cout << endl;
}
*/