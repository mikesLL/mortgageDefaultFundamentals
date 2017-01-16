#include "headers.h"

// Initialize state distribution, wealth distribution
def_stats::def_stats( void *snodes_in ) {

	snodes1 = (snodes *)snodes_in;
	//t_hor = (*snodes1).t_hor;

	vector<double> zeros_NS(n_s, 0.0);
	vector<double> zeros_WN(w_n, 0.0);
	sdist = zeros_NS;                                     // Initial state prob distribution
	wdist = zeros_WN;                                     // Initial Wealth distribution


	int t_max = 40; // TODO: set = something from calibration
	vector<double> zeros_tmax(t_max, 0.0);
	hazard_store = zeros_tmax; 

	vector<vector<double>> zeros_tmax_NS( t_max, vector<double> ( n_s, 0.0) );
	vector<vector<double>> zeros_tmax_WN(t_max, vector<double>( w_n, 0.0) );

	sdist_store = zeros_tmax_NS;
	wdist_store = zeros_tmax_WN;

	// Initial conditions
	int t_hor = 0;
	sdist[floor(n_s / 2.0)] = 1.0;                                  // First Year: impose in median state
	double w_med = 0.4;                                             // impose median wealth = 40k
	int i_w_med = int(floor((w_med - w_min) / (w_max - w_min)));    // convert to wealth index
	wdist[i_w_med] = 1.0;                                           // impose mass of wealth is in median of w_n

	sdist_store[t_hor] = sdist;
	wdist_store[t_hor] = wdist;
}


void def_stats::wtrans_iterate(int t_hor_in) {
	int t_hor = t_hor_in;
	int i_t_hor;
	double hazard_rate;                   // hazard rate

	vector<vector<double>> gammap;                      // State transition matrix															
	vector<double> zeros_WN(w_n, 0.0);
	vector<double> zeros_NS(n_s, 0.0); 

	// Indices
	int i_s1p, i_s2p, i_w1p; // current state, next state, current wealth
	int i_w2l, i_w2h;        // wealth in low state, high state (next period)
	double w2_l, w2_h;       // wealth in low state, high state (next period)

	for (i_t_hor = 0; i_t_hor < t_hor; i_t_hor++) {

		hazard_rate = 0.0;
		wdist = wdist_store[i_t_hor];     // set wealth distribution
		sdist = sdist_store[i_t_hor];     // set state distribution
		wdist2 = zeros_WN;
		sdist2 = zeros_NS;

		gammap = (*snodes1).gammat[i_t_hor];     // State transition matrix

		// COMPUTE: wdist2: wealth distribution in the next period
		for (i_s1p = 0; i_s1p < n_s; i_s1p++) {                         // cycle through possible states
			for (i_s2p = 0; i_s2p < n_s; i_s2p++) {
				for (i_w1p = 0; i_w1p < w_n; i_w1p++) {
					w2_l = (*snodes1).w_t2_state[t_hor][i_s1p][i_s2p][i_w1p][0];    // load in wealth in next period given state
					w2_h = (*snodes1).w_t2_state[t_hor][i_s1p][i_s2p][i_w1p][1];    // (low and high wealth realizations)

					i_w2l = round((w2_l - w_min) / (w_max - w_min));   // round wealth realization to closest wealth index
					i_w2l = min(max(i_w2l, 0), w_n - 1);               // TODO: double check round fn

					i_w2h = round((w2_h - w_min) / (w_max - w_min));
					i_w2h = min(max(i_w2h, 0), w_n - 1);

					wdist2[i_w2l] = wdist2[i_w2l] + 0.5*sdist[i_s1p] * gammap[i_s1p][i_s2p] * wdist[i_w1p];  // add outcome probability mass
					wdist2[i_w2h] = wdist2[i_w2h] + 0.5*sdist[i_s1p] * gammap[i_s1p][i_s2p] * wdist[i_w1p];  // 0.5: equity return
				}

				sdist2[i_s2p] = sdist2[i_s2p] + sdist[i_s1p] * gammap[i_s1p][i_s2p];
			}

			hazard_rate = hazard_rate + sdist[i_s1p] * wdist[i_w1p] * (1.0 - double((*snodes1).own_state[t_hor][i_s1p][i_w1p]));
		}

		wdist_store[i_t_hor + 1] = wdist2;
		sdist_store[i_t_hor + 1] = sdist2;
		hazard_store[i_t_hor] = hazard_rate;
	}
}

// COMPUTE: default hazard rate
//for (i_s1p = 0; i_s1p < n_s; i_s1p++) {                       // cycle through current state
//	for (i_w1p = 0; i_w1p < w_n; i_w1p++) {                   // cycle through wealth
//		hazard_rate = hazard_rate + sdist[i_s1p] * wdist[i_w1p] * (1.0 - double((*snodes1).own_state[t_hor][i_s1p][i_w1p]));
//	}
//}

//wdist_store[t_hor + 1] = wdist2;

void def_stats::print_def_stats( int t_hor_in) {
	int t_hor = t_hor_in;
	int i_t_hor;

	// PRINT: hazard store
	ofstream v1_file;                                           // open output file stream 
	string file_name1 = "def_results/hazard_rate.csv";
	v1_file.open(file_name1, ios::out | ios::trunc);              // outstream, truncate

	// print headers
	v1_file << "t_hor,hazardRate" << endl;
	for (i_t_hor = 0; i_t_hor < t_hor; i_t_hor++) {
		v1_file << i_t_hor << "," << hazard_store[i_t_hor] << "," << endl;
	}
	v1_file << endl;
	v1_file.close();


	// PRINT: sdist store
	ofstream v2_file;                                           // open output file stream 
	string file_name2 = "def_results/sdist_store.csv";
	v2_file.open(file_name2, ios::out | ios::trunc);              // outstream, truncate
	int i_sf;
																 // print headers
	v2_file << "t_hor,sdistStore" << endl;
	for (i_t_hor = 0; i_t_hor < t_hor; i_t_hor++) {
		v2_file << i_t_hor << ",";
		
		for(i_sf = 0; i_sf < n_s; i_sf++) {
			v2_file << sdist_store[i_t_hor][i_sf] << "," ;
		}
		v2_file << "," << endl; 
	}
	v2_file << endl;
	v2_file.close();

	// PRINT: wdist store
	ofstream v3_file;                                           // open output file stream 
	string file_name3 = "def_results/wdist_store.csv";
	v3_file.open(file_name3, ios::out | ios::trunc);              // outstream, truncate
	int i_wf;
	// print headers
	v3_file << "t_hor,wdistStore" << endl;
	for (i_t_hor = 0; i_t_hor < t_hor; i_t_hor++) {
		v3_file << i_t_hor << ",";

		for (i_wf = 0; i_wf < w_n; i_wf++) {
			v3_file << wdist_store[i_t_hor][i_wf] << ",";
		}
		v3_file << "," << endl;
	}
	v3_file << endl;
	v3_file.close();

	// RUN CODE TO PRINT OUT 
	// hazard_store, sdist_store, wdist_store
};





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
/*
void def_stats::wtrans( int t_hor_in ) {
t_hor = t_hor_in;

wdist = wdist_store[t_hor];     // set wealth distribution
sdist = sdist_store[t_hor];     // set state distribution

vector<vector<double>> gammap = (*snodes1).gammat[t_hor];     // State transition matrix

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
w2_l = (*snodes1).w_t2_state[t_hor][i_s1p][i_s2p][i_w1p][0];    // load in wealth in next period given state
w2_h = (*snodes1).w_t2_state[t_hor][i_s1p][i_s2p][i_w1p][1];    // (low and high wealth realizations)

i_w2l = round((w2_l - w_min) / (w_max - w_min));   // round wealth realization to closest wealth index
i_w2l = min(max(i_w2l, 0), w_n - 1);               // TODO: double check round fn

i_w2h = round((w2_h - w_min) / (w_max - w_min));
i_w2h = min(max(i_w2h, 0), w_n - 1);

wdist2[i_w2l] = wdist2[i_w2l] + 0.5*sdist[i_s1p] * gammap[i_s1p][i_s2p] * wdist[i_w1p];  // add outcome probability mass
wdist2[i_w2h] = wdist2[i_w2h] + 0.5*sdist[i_s1p] * gammap[i_s1p][i_s2p] * wdist[i_w1p];  // 0.5: equity return

}
}
}

// COMPUTE: default hazard rate
for (i_s1p = 0; i_s1p < n_s; i_s1p++) {                       // cycle through current state
for (i_w1p = 0; i_w1p < w_n; i_w1p++) {                   // cycle through wealth
hazard_rate = hazard_rate + sdist[i_s1p] * wdist[i_w1p] * (1.0 - double((*snodes1).own_state[t_hor][i_s1p][i_w1p]));
}
}

wdist_store[t_hor + 1] = wdist2;
}
*/