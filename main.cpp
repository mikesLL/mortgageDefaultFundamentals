/*
This program numerically solves a household's (HH's) lifecycle portfolio optimization 
problem when the HH is able to take S&P Case-Shiller Index Futures positions. The HH's value fn 
and status in t = 11 (11 year horizon) is modeled as vf_11. Given vf_11, program solves for 
vf_10,..., vf_1 by value function iteration using gen_v1.

An integral part of the portfolio problem is estimating the future 
of home prices. I simulate the home price path in load_rand.

pstruct contains possible home prices in given number of time periods
gstruct contains a transition matrix for each time period.
sd_read_in.csv, sf_read_in.csv, and lax_read_in.csv include housing related data for 
that particular year (median price, median rent, lagged home price appreciation, 
Case-Shiller Inex Futures Price, etc...)
(not included on GitHub due to proprietary data)
Currently works for on data from 2007-2014
Cities: San Diego, San Francisco, Los Angeles, Boston, Chicago, 
Denver, Miami, New York 

Copyright A. Michael Sharifi, 2016
*/

#include "headers.h"


int main(){
	string city_init, city_filename;                                   // city name
	hdata city_data;                                     // housing data structure stores previous rents and lagged returns
	
	const string city_init_vec[] = {"sd", "sf", "lax", "bos", "chi", "den", "mia", "nym"  };
	const string city_filename_vec[] = { "sd_read_in.csv", "sf_read_in.csv", "lax_read_in.csv",
		                                 "bos_read_in.csv", "chi_read_in.csv", "den_read_in.csv",
		                                 "mia_read_in.csv", "nym_read_in.csv" };

	int age_begin_store2[5][3] = {  {0, 0, 0}, { 30, 60, 45 }, { 30, 60, 45 }, { 60, 45, 30 }, { 60, 45, 30 } };

	age_begin_store2[0][0] = age_begin_store[0];
	age_begin_store2[0][1] = age_begin_store[1];
	age_begin_store2[0][2] = age_begin_store[2];

	int city_id;
	for (city_id = city_begin; city_id <= city_end; city_id++) {

		city_filename = city_filename_vec[city_id];
		city_init = city_init_vec[city_id];
		cout << "city_init = " << city_init << endl;
		load_csv(&city_data, city_filename);

		int i_age, t, age0;
	
		double grent_store[] = { 0.02, 0.0 };  // rent growth 
		double ltv0_store[] = {0.8, 0.95};
		double rp0_store[] = { 0.06, 0.04 };   // rp = 0.045;
			
		int id = 0, id_grent, id_ltv0, id_rp0;
		
		t = 11;

		double gr = 0.025, gr1 = 0.0;
		double ltv = 0.8, ltv1 = 0.9;
		double rp = 0.06, rp1 = 0.045; 

		// order: rent growth, LTV0, rp0 mult
		/*
		double param_store[8][3] = { { 0.0, 0.8, 1.0 },     // baseline case,
		                             { 1.0, 0.8, 1.0 },     //
							         { 0.0, 0.8, 1.5 },     // low R/P
							         { 1.0, 0.8, 1.5 },     // 
							         { 0.0, 0.90, 1.0 },     // high ltv
							         { 1.0, 0.90, 1.0 },     // 
							         { 0.0, 0.90, 1.5 },     // low R/P and High LTV
							         { 1.0, 0.90, 1.5 } };
		*/

		// as you cycle through parameterizations, the r/p will decrease
		double param_store[8][3] = { { 1.0, 0.8, 0.6 },     // baseline case,
		                             { 1.0, 0.8, 0.8 },     //
		                             { 1.0, 0.8, 1.0 },     // low R/P
		                             { 1.0, 0.8, 1.2 },     // 
		                             { 1.0, 0.8, 1.4 },     // high ltv
		                             { 1.0, 0.8, 1.6 },     // 
		                             { 1.0, 0.8, 1.8 },     // low R/P and High LTV
		                             { 1.0, 0.8, 2.0 } };

		/*
		Initial rent-price ratio: 
		rp + 0.02, rp - 0.02? Nice
		*/
		#pragma omp parallel for
		
		for (id = 0; id <= 3; id++) {                      // id <= 7

			int grent_id = round(param_store[id][0]);     // set = 0 for low rent growth, = 1 for high rent growth
			double grent = 0.0;                           //param_store[id][0];      //grent_store[id_grent];
			double ltv0 = param_store[id][1];             //ltv0_store[id_ltv0];

			double phr_in, gamma0_est, gamma1_est, ph0;

			double ph0_low = 1.3188, ph0_high = 1.5695, ph0_def;
			double rp0;

			if (false){
				if (grent_id) {
					ph0_def = ph0_high;
					ph0 = ph0_high * param_store[id][2];
					gamma0_est = -2.8092;
					gamma1_est = 0.7155;
				}
				else {
					ph0_def = ph0_low;
					ph0 = ph0_low * param_store[id][2];
					gamma0_est = -2.8164;
					gamma1_est = 0.5796;
				}
				phr_in = exp(gamma0_est)*pow(ph0_def, gamma1_est); // I guess do this for now using average home price
		    }
			
			if (true) {
				grent_id = 1; 
				ph0_def = ph0_high;
				phr_in = 0.057*ph0_def;                 // holds rent fixed 
				ph0 = ph0_high * param_store[id][2];    // ph0 can be over or undervalued
				gamma0_est = -2.8092;
				gamma1_est = 0.7155;
			}

			/*
			rent_avg_high : 0.0773
			rent_avg_low : 0.0680
			price_avg_high : 1.5695
			price_avg_low : 1.3188
			rp_avg_high : 0.0552
			rp_avg_low : 0.0564
			*/
			
			age0 = 30;
			int T_max = age_max - age0;            // optimization problem horizon
			int t_hor = T_max;                     // household's planning horizon time index
			
			double duration;			
			//double ph0 = 1.0 / rp0 * phr_in;       // rent to price parameter imposes current price
			
			clock_t start = clock();                // discretized states including home prices, rents, incomes
			snodes snodes1(age0, T_max, city_id);
			
			cout << "main.cpp: snodes1 finished" << endl;

			def_stats def_stats1(&snodes1); 
				
			cout << "main.cpp : finished running def_stats" << endl;

			cout << "snodes1.rent_adj = " << snodes1.rent_adj << endl; 

			cout << "housing tenure: " << endl;
			cout << snodes1.hu_ten[0] << "..." << snodes1.hu_ten[1] << "..." << snodes1.hu_ten[2] << "..." << endl;

			cout << "housing wealth weight: " << endl;
			cout << snodes1.ten_w[0] << "..." << snodes1.ten_w[1] << "..." << snodes1.ten_w[2] << "..." << endl;

			vfn vf_F;
			vfn vf_P;

			double mapr = 0.06;
			int N_term = 30;
		
			cout << "Compute Mortgage Path" << endl;
			mortg mortg1(&snodes1, ph0, ltv0, mapr, N_term); 

			cout << "Simulate/Discretize Macro and Housing State-space" << endl;
		
			load_simpath(&snodes1, grent_id, grent, phr_in, ph0, city_data.ret_lag[t],
					city_data.csf_1yr[t], t, city_init, city_id, age0);

			cout << "main.cpp: begin enter data" << endl;
			vf_F.enter_data(&snodes1, phr_in, t, t_hor, city_data.csf_1yr[t], pref, T_max);

			cout << "main.cpp: set terminal" << endl;
			vf_F.set_terminal(&mortg1, phr_in );

			cout << "main.cpp: store_data" << endl;
			store_data(&snodes1, &vf_F, city_init, t, t_hor, T_max);

			for (t_hor = (T_max - 1); t_hor >= 0; t_hor--) {
					snodes1.t_hor = t_hor;
					
					cout << "main.cpp: vf_P: enter data" << endl;
					vf_P.enter_data(&snodes1, phr_in, t, t_hor, city_data.csf_1yr[t], pref, T_max);

					cout << " main.cpp: gen_VP" << endl;
					gen_VP( &snodes1, &mortg1, &vf_P, &vf_F);

					cout << "main.cpp: begin store_data" << endl;
					store_data(&snodes1, &vf_P, city_init, t, t_hor, T_max);

					vf_F = vf_P;
				}
				def_stats1.wtrans_iterate(T_max);
				def_stats1.print_def_stats(T_max, id);   // print out transition matrices

				t_hor = 0;

				cout << "\n Value Function Iteration Completed \n";
				duration = (clock() - start) / (double)CLOCKS_PER_SEC;
				cout << "total run time: " << duration << '\n';
		}
		
	}

	cin.get();
	return 0;
}


