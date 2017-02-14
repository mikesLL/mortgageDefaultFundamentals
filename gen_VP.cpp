/*
Given the current value fn VFN_3d_2, this program solve for the preceeding value fn VFN_3d_1

Start with the renter problem: t_i= 0 and loop through home price states ph_i and wealth states w_i
Given each state ph_i and w_i, cycle through next period tenure t_i2 = 0,1,2 
Solve optimization problem for that inidividual state
A solution includes policies x[5] = {C, B, X, CSFp, CSFn}, t_i2 in {0,1,2}
and a value of the current state V1()

C: Consumption 
B: Bond (Mortgage) holdings
X: Equities
CSFp: Case-Shiller Index Futures, Long Position, Margin
CSFn: Case-Shiller Index Futures, Short Position, Margin
t_i, t_i2 in {0: rent, 1: starter home, 2: full-size home}

Copyright A. Michael Sharifi, 2016

TODO: add an algorithm to use the concavity of the
value function to generate quick guesses at the begining of the state
also exploits fact that policies / strategies should be smooth
Would be kind of like the MCMC analog to VFI
*/

#include "headers.h"

void gen_VP(void *snodes_in, void *mortg_in, void *VFN_3d_1, void *VFN_3d_2 ){
	
double duration;
clock_t start;

snodes *snodes1 = (snodes *)snodes_in;
mortg *mortg1 = (mortg *)mortg_in;

vfn *rr2 = (vfn *) VFN_3d_2;          // address to initialized V2
vfn *rr1 = (vfn *) VFN_3d_1;          // address to initialized V1

gen_res res1;
eval_res res0;
eval_res res_t_lag;
eval_res res_t_0;
	
double b_min;
	
double w_refi;
eval_res res_refi;
int t_refi = 1;

vector<double> theta1;
	
int y_i; 
int t_i, t_i2 , ph_i, w_i, t_adj;
int i_rent, i_ph, i_yi;

int i_s;
int t_hor = (*snodes1).t_hor;

if (t_hor == 0) {
	cout << "first horizon" << endl; 
}

double coh, w_adj, v_adj, v_i_floor, v0_opt, v0, beg_equity, mpmt, b_min2;
double v1;
double v_lag_w; // value function guess from previous wealth level
double v_lag_t; // value function guess from previous tenure computation

vector<double> x(5,0.0);                     // policy / strategy
vector<double> x_lag_w(5, 0.0);
vector<double> x0(5, 0.0);
vector<double> x1(5, 0.0);

vector<double> x_guess(5, 0.0);

vector<double> x_ti1(5, 0.0);       // policy: owns smallest home
double v_ti1 = -1.0e20;              // value: owns smallest home

double h_mult = 0.0;
double bond_est = 0.0;
	
vector<vector<double>> x_lag_wt(t_n, vector<double> (5, 0.0) );

for (i_ph = 0; i_ph < n_ph; i_ph++) {
	cout << (*snodes1).p_gridt[t_hor][i_ph] << "..." << endl;
}

t_i = 0; // consider t_i = 0 first; case: begin with renter 
int res1_flag = 1;
int t_i2_lag_w = 0;

// MODS HERE: add cycle for i_m
int i_m = 0; // mortgage state
int i_rcurr, i_rpmt, i_rlb;  // mortgage state: current rate, payment rate, loan balance

// assume the household cannot sell the home this year and then buy back in the future
// permanent renter case: set t_i2 = 0!
t_i2 = 0;
	
int i_urate, i_fedfunds, i_plevel;

double ltv;


// COMPUTE Renter problem
for (i_s = 0; i_s < n_s; i_s++) {


	// /i_yi = (*snodes1).s2i_yi[i_s];              // load in states
	// pass in URATE, PLEVEL, and FEDFUNDS
	(*rr2).urate1 = (*snodes1).urate_gridt[t_hor][(*snodes1).s2i_urate[i_s]];                   // pass in initial unemployment rate
	(*rr2).yinc1 = (*snodes1).yi_gridt[t_hor][0];                                               // pass in income
	(*rr2).fed_funds1 = (*snodes1).fedfunds_store[t_hor][i_s];                                  // pass in fed funds
	(*rr2).plevel1 = (*snodes1).plevel_gridt[t_hor][(*snodes1).s2i_plevel[i_s]];

	i_yi = 0; //i_yi = (*snodes1).s2i_yi[i_s];              // load in states
	i_rent = 0;                                   // (*snodes1).s2i_rent[i_s];
	i_ph = (*snodes1).s2i_ph[i_s];                // load in states

	start = clock();
	cout << "i_m = " << i_m << "  i_s = " << i_s << "  t_i = " << t_i << "  i_ph = " << i_ph << "  begin renter problem" << endl;
	
	for (w_i = 0; w_i < w_n; w_i++) {
		i_m = 0;                                                       // No Mortgage
		t_i = 0;                                                       // Renter

		if (w_i % 100 == 0) {
			cout << "w_i = " << w_i << endl;
		}
			
		b_min = b_min_unsec;
		beg_equity = -1.0e6;
		mpmt = -1.0e6;             // TODO: check on this

		// load previous w_i policy as a benchmark
		(*rr1).get_pol(t_i, i_m, i_s, w_i - 1, x_lag_w);
		t_i2_lag_w = t_i2; // tenure choice: previous wealth  
		v_lag_w = (*rr1).vw3_grid[t_i][i_m][i_s][max(w_i - 1, 0)];

		// load t_i2-restricted guess as an initial starting point
		x_guess = x_lag_wt[t_i2];

		// compute cash on hand
		coh = (*rr1).w_grid[w_i] - (*snodes1).rent_gridt[t_hor][i_rent] * (*snodes1).rent_adj +
			(*snodes1).yi_gridt[t_hor][(*snodes1).s2i_yi[i_s]];
	    //(*snodes1).y_gridt[t_hor][i_yi = (*snodes1).s2i_yi[i_s];]              // load in states

		(*rr2).i_s1 = i_s;          // pass in current state to next-period value function
		(*rr2).t_i1 = t_i;
		(*rr2).i_m1 = i_m; 
		(*rr2).w_i1 = w_i;
		(*rr2).t_i2 = t_i2;         //  t_i2 is a choice variable, so adding it to fn pointer 
		(*rr2).def_flag = 0;

		(*snodes1).i_s1 = i_s;      // load states into snodes
		(*snodes1).i_w1 = w_i;
		(*snodes1).loan_bal1 = 0.0; 
		
		res1 = gen_VPw(snodes1, rr1, rr2, coh, x_guess, b_min, beg_equity, mpmt);                               // pass problem into gen_V1_w

		v1 = res1.v_opt;                 // guess for the current vfn given current t_i2
		x1 = res1.x_opt;                 // guess for current policy given current t_i2

		(*rr1).v_move[w_i] = v1;         // set moving value fn
	
		x_lag_wt[t_i2] = x1;             // store x_opt as a guess for the next wealth under same tenure

		for (i_m = 0; i_m < m_n; i_m++) {                          // set opt result for all mort states
			(*rr1).set_pol_ten_v(t_i, i_m, i_s, w_i, x1, t_i2, v1);      

			if (res1.valid_flag == 0) {
				(*rr1).set_pol_ten_v(t_i, i_m, i_s, w_i, x1, 0, v1);
			}
		}
	}

	for (i_m = 0; i_m < m_n; i_m++) {
		(*rr1).interp_vw3(t_i, i_m, i_s);  // clean and interpolate grid	
	}
	duration = (clock() - start) / (double)CLOCKS_PER_SEC;
	cout << "time elapsed: " << duration << '\n';
}


// MODS HERE
int i_m_refi;  // mortgage state associated with refinancing
double mortg_pmt3;  // mortgage payment
double loan_diff;  // loan balance difference
double loan_bal; // current loan balance


cout << "gen_VP.cpp: begin homeowner problem" << endl;

double v_refi;
int t_def, i_m_def;
double w_def;
eval_res res_def;
// TODO: think about ordering computation cases

int t_sell, i_m_sell;
double w_sell;
eval_res res_sell;
int t_i2_nonvalid = 0;

for (t_i = 1; t_i < t_n; t_i++) {                        // Cycle through homeowner problem

	for (i_s = 0; i_s < n_s; i_s++) {                  // Cycle through macro states

		// pass in URATE, FEDFUNDS, YINC
		(*rr2).urate1 = (*snodes1).urate_gridt[t_hor][(*snodes1).s2i_urate[i_s]];            // pass in initial unemployment rate
		(*rr2).fed_funds1 = (*snodes1).fedfunds_store[t_hor][i_s];                            // pass in fed funds
		(*rr2).yinc1 = (*snodes1).yi_gridt[t_hor][0];                                         // pass in income
		(*rr2).plevel1 = (*snodes1).plevel_gridt[t_hor][(*snodes1).s2i_plevel[i_s]];          // pass in plevel

		t_i2 = t_i;                                                                // impose t_i2 = t_i (homeowner)

		i_ph = (*snodes1).s2i_ph[i_s];
		
		start = clock();
		cout << "i_s = " << i_s << "   t_i = " << t_i << "i_yi = " << i_yi
			<< " i_ph = " << i_ph << " begin homeowner problem " << endl;

		for (i_m = 0; i_m < m_n; i_m++) {                        // Cycle through mortgage states (m_n = 2)        

			if (i_m > 0) {         
				mpmt = (*mortg1).mpmt;                           // Compute mortgage payment, loan balance
				loan_bal = (*mortg1).loan_bal[t_hor];

				if ((*snodes1).s2i_yi[i_s] >= 1) {
					mpmt = (*mortg1).mpmt - (*mortg1).mapr*loan_bal; 
				}
				
			} else {
				mpmt = 0.0;
				loan_bal = 0.0;
			}

			for (w_i = 0; w_i < w_n; w_i++) {

				b_min = -max_ltv*((*snodes1).p_gridt[t_hor][i_ph] - loan_bal) + 0.0*b_min_unsec;  // get rid of b_min_unsec for now

				coh = (*rr1).w_grid[w_i] - mpmt + (*snodes1).yi_gridt[t_hor][(*snodes1).s2i_yi[i_s]];
					// COH = liquid assets + income - mortgage payment)

				// load previous w_i policy as an initial guess
				(*rr1).get_pol(t_i, i_m, i_s, w_i - 1, x_lag_w);                     // get x pol sol from previous w_i and assign to x
				t_i2_lag_w = (*rr1).xt_grid[t_i][i_m][i_s][max(w_i - 1, 0)];
				v_lag_w = (*rr1).vw3_grid[t_i][i_m][i_s][max(w_i - 1, 0)];
				(*rr2).i_m1 = i_m;

				beg_equity = -1.0e6;

				(*rr2).i_s1 = i_s;
				(*rr2).t_i1 = t_i;
				(*rr2).w_i1 = w_i;
				(*rr2).t_i2 = t_i2;

				(*snodes1).i_s1 = i_s;      // load states into snodes
				(*snodes1).i_w1 = w_i;
				(*snodes1).loan_bal1 = loan_bal;

				// solve optimization problem
				res1 = gen_VPw(snodes1, rr1, rr2, coh, x_lag_w, b_min, beg_equity, mpmt);
				v1 = res1.v_opt;                                                 // guess for current solution: vfn
				x1 = res1.x_opt;                                                 // guess for current policy

				
				(*rr1).set_pol_ten_v(t_i, i_m, i_s, w_i, x1, t_i2, v1);          // store opt result

				if (res1.valid_flag <= 0){
					(*rr1).set_pol_ten_v(t_i, i_m, i_s, w_i, x1, t_i2_nonvalid, v1);          // check valid flag
				}
					
				// CASE: HH SELLS
				t_sell = 0;
				i_m_sell = 0;

				w_sell = (*rr1).w_grid[w_i] + (1.0 - phi)*(*snodes1).p_gridt[t_hor][i_ph] - loan_bal;        // calc wealth if HH sells
				res_sell = (*rr1).eval_v(t_sell, i_m_sell, i_s, w_sell);                                    // eval value fn if HH sells
			
				// CASE: value of selling > value of owning
				if ( ( (w_sell >= 0.0) && (res_sell.v_i_floor > v1) && (res_sell.w_i_floor >= 0) ) || (res1.valid_flag <= 0) ) {

					(*rr1).get_pol(t_sell, i_m_sell, i_s, res_sell.w_i_floor, x);                // submit x as reference and load in x pol from t1 = 0
					(*rr1).set_pol_ten_v(t_i, i_m, i_s, w_i, x, t_sell, res_sell.v_i_floor);     // first arguments are current state variables, x containts updated policy

					v1 = res_sell.v_i_floor;                                                     // update v1 with value of default
					(*snodes1).own_state[t_hor][i_s][w_i] = 0;                                   // update ownership state
					(*snodes1).w_state_swap(res_sell.w_i_floor);                                  // update wealth transition state
				}

				if(i_m >= 1){
					// CASE: HH REFINANCES
					t_refi = 1;                                                          // Refi: still owns
					i_m_refi = 0;                                                        // Refi: mortgage index state
					w_refi = (*rr1).w_grid[w_i] - loan_bal;                              // If HH refinances, liquid wealth declines by loan_bal
					res_refi = (*rr1).eval_v(t_refi, i_m_refi, i_s, w_refi);                  // Evaluate value fn if household refinances

																							  // check if HH refinances (compare value fn's and check validity)
					if ( (w_refi >= 0.0) && (res_refi.v_i_floor > v1) && (res_refi.w_i_floor >= 0) ) {
						(*rr1).get_pol(t_refi, i_m_refi, i_s, res_refi.w_i_floor, x);              // load in policy associated with refinance
						(*rr1).set_pol_ten_v(t_i, i_m, i_s, w_i, x, t_refi, res_refi.v_i_floor);   // set policy associated with refinance
						v1 = res_refi.v_i_floor;                                                   // upadate value fn guess
						(*snodes1).w_state_swap(res_refi.w_i_floor);                               // update wealth transition state
					}

					// COMPUTE CASE: HH DEFAULTS
					t_def = 0;
					i_m_def = 0;
					w_def = (*rr1).w_grid[w_i];                                // calc wealth if HH defaults
					res_def = (*rr1).eval_v(t_def, i_m_def, i_s, w_def);               // eval value fn if HH defaults

					// CASE: value of default > value of owning
					if ( ( (w_def >= 0.0) && (res_def.v_i_floor > v1) && (res_def.w_i_floor >= 0) ) ) {
						(*rr1).get_pol(t_def, i_m_def, i_s, res_def.w_i_floor, x);                         // submit x as reference and load in x pol from t1 = 0
						(*rr1).set_pol_ten_v(t_i, i_m, i_s, w_i, x, t_def, res_def.v_i_floor);     // first arguments are current state variables, x containts updated policy

						v1 = res_def.v_i_floor;                        // update v1 with value of default
						(*snodes1).own_state[t_hor][i_s][w_i] = 0;     // update ownership state
						(*snodes1).def_state[t_hor][i_s][w_i] = 1;     // update default state
						(*snodes1).w_state_swap(res_def.w_i_floor);    // update wealth transition state

						ltv = loan_bal / (*snodes1).p_gridt[t_hor][(*snodes1).s2i_ph[i_s]];  // compute ltv
						(*snodes1).def_ltv_state[t_hor][i_s][w_i] = ltv;                     // store

						if (ltv <= 0.9) {
							cout << "gen_VP.cpp: ltv error" << endl; 
						}

					}
				}

			}

			(*rr1).interp_vw3(t_i, i_m, i_s);  // clean grid

		}

		duration = (clock() - start) / (double)CLOCKS_PER_SEC;
		cout << "time elapsed: " << duration << "  lcount = " << (*rr2).lcount << endl;
	}
}


// clean grid again
for (t_i = 0; t_i < t_n; t_i++) {  
	for (i_m = 0; i_m < m_n; i_m++) {
		for (i_s = 0; i_s < n_s; i_s++) {
			(*rr1).interp_vw3(t_i, i_m, i_s);
		}
	}
}

cout << " gen" << to_string( (*rr1).t_num) << "completed" << endl;
}


//b_min = -max_ltv*(*snodes1).ten_w[t_i2] * (*snodes1).p_gridt[t_hor][i_ph] + b_min_unsec;
//b_min2 = -max_lti * (*snodes1).yi_gridt_btax[t_hor][i_yi] / (rb + mort_spread - 1.0);


//************************************ 
// OLD VERSION
//************************************
/*
cout << "gen_VP.cpp: begin homeowner problem" << endl;
for (t_i = 1; t_i < t_n; t_i++) {                        // Cycle through homeowner problem
for (i_s = 0; i_s < n_s; i_s++) {

t_i2 = t_i;                                          // impose t_i2 = t_i
i_yi = (*snodes1).s2i_yi[i_s];
i_rent = (*snodes1).s2i_rent[i_s];
i_ph = (*snodes1).s2i_ph[i_s];
i_rcurr = (*snodes1).s2i_rm[i_s];                    // short rate is a state var

start = clock();
cout << "i_s = " << i_s << "   t_i = " << t_i << "i_yi = " << i_yi
<< " i_ph = " << i_ph << " begin homeowner problem " << endl;

for (i_m = 0; i_m < m_n; i_m++) {                       // Cycle through mortgage states

// Load mortgage state

i_rpmt = (*mortg1).m2rpmt_map[i_m];        // current payment rate on mortgage
i_rlb = (*mortg1).m2rlb_map[i_m];          // loan balance (given ammort rate)

// TODO: verify this works
//mortg_pmt3 = (*mortg1).pmt[i_m][t_hor]; // compute mortgage payment in each state
mpmt = (*mortg1).pmt[i_m][t_hor];                   // mortgage payment

for (w_i = 0; w_i < w_n; w_i++) {

loan_bal = (*mortg1).bal[i_rlb][t_hor];                // retrieve current loan balance
b_min = -max_ltv*((*snodes1).p_gridt[t_hor][i_ph] - loan_bal) + 0.0*b_min_unsec;  // get rid of b_min_unsec for now

coh = (*rr1).w_grid[w_i] + y_atax*(*snodes1).yi_gridt[t_hor][i_yi] - mpmt;  // COH = liquid assets + income - mortgage payment)

// load previous w_i policy as an initial guess
(*rr1).get_pol(t_i, i_m, i_s, w_i - 1, x_lag_w);                     // get x pol sol from previous w_i and assign to x
t_i2_lag_w = (*rr1).xt_grid[t_i][i_m][i_s][max(w_i - 1, 0)];
v_lag_w = (*rr1).vw3_grid[t_i][i_m][i_s][max(w_i - 1, 0)];
(*rr2).m_i1 = i_m;

//mpmt = 0.0; //-1.0e6;                // FOR NOW: set = 0.0; think about getting rid of
beg_equity = -1.0e6;

(*rr2).i_s1 = i_s;
(*rr2).t_i1 = t_i;
(*rr2).w_i1 = w_i;
(*rr2).t_i2 = t_i2;

(*snodes1).i_s1 = i_s;      // load states into snodes
(*snodes1).i_w1 = w_i;
(*snodes1).loan_bal1 = loan_bal;

// solve optimization problem
res1 = gen_VPw(snodes1, rr1, rr2, coh, x_lag_w, b_min, beg_equity, mpmt);
v1 = res1.v_opt;                                                 // guess for current solution: vfn
x1 = res1.x_opt;                                                 // guess for current policy

(*rr1).set_pol_ten_v(t_i, i_m, i_s, w_i, x1, t_i2, v1);          // store opt result

// COMPUTE CASE: HH REFINANCES
if (i_rcurr < i_rpmt) {
t_refi = 1;                                                      // tenure index in event of refi

i_m_refi = (*mortg1).r2m_map[i_rcurr][i_rlb];                // refi: i_rpmt = i_rcurr, i_rlb is the same

// loan balance: ammort at repayment rate - ammort at current market rate
loan_diff = (*mortg1).bal[i_rpmt][t_hor] - (*mortg1).bal[i_rcurr][t_hor];

w_refi = (*rr1).w_grid[w_i] - loan_diff;                                        // compute liquid wealth if HH refinances
res_refi = (*rr1).eval_v(1, i_m_refi, i_s, w_refi);                             // evaluate value fn if household refinances

if ( (w_refi >= 0.0) && (res_refi.v_i_floor > v1) && (res_refi.w_i_floor >= 0)) {

(*rr1).get_pol(t_refi, i_m_refi, i_s, res_refi.w_i_floor, x);              // load in policy associated with refinance
(*rr1).set_pol_ten_v(t_i, i_m, i_s, w_i, x, t_refi, res_refi.v_i_floor);   // set policy associated with refinance
v1 = res_refi.v_i_floor;                                                   // upadate value fn guess

(*snodes1).w_state_swap(res_refi.w_i_floor);                               // update wealth transition state
}
}

// COMPUTE CASE: HH SELLS
w_adj = (*rr1).w_grid[w_i] + (*snodes1).p_gridt[t_hor][i_ph] - loan_bal;           // calc wealth if HH sells
res_t_0 = (*rr1).eval_v(0, i_m, i_s, w_adj);               // eval value fn if HH sells

// CASE: value of selling > value of owning
if ((w_adj >= 0.0) && (res_t_0.v_i_floor > v1) && (res_t_0.w_i_floor >= 0)) {
(*rr1).get_pol(0, i_m, i_s, res_t_0.w_i_floor, x);                         // submit x as reference and load in x pol from t1 = 0
t_adj = (*rr1).xt_grid[0][i_m][i_s][res_t_0.w_i_floor];                    // get t2 pol from t1 = 0; simulated sale
(*rr1).set_pol_ten_v(t_i, i_m, i_s, w_i, x, t_adj, res_t_0.v_i_floor);     // first arguments are current state variables, x containts updated policy

v1 = res_t_0.v_i_floor;                        // update v1 with value of default
(*snodes1).own_state[t_hor][i_s][w_i] = 0;     // update ownership state
(*snodes1).w_state_swap(res_t_0.w_i_floor);    // update wealth transition state
}

// COMPUTE CASE: HH DEFAULTS
w_adj = (*rr1).w_grid[w_i];                                // calc wealth if HH defaults
res_t_0 = (*rr1).eval_v(0, i_m, i_s, w_adj);               // eval value fn if HH defaults

// CASE: value of default > value of owning
if ((w_adj >= 0.0) && (res_t_0.v_i_floor > v1) && (res_t_0.w_i_floor >= 0)) {
(*rr1).get_pol(0, i_m, i_s, res_t_0.w_i_floor, x);                         // submit x as reference and load in x pol from t1 = 0
t_adj = (*rr1).xt_grid[0][i_m][i_s][res_t_0.w_i_floor];                    // get t2 pol from t1 = 0; simulated sale
(*rr1).set_pol_ten_v(t_i, i_m, i_s, w_i, x, t_adj, res_t_0.v_i_floor);     // first arguments are current state variables, x containts updated policy

v1 = res_t_0.v_i_floor;                        // update v1 with value of default
(*snodes1).own_state[t_hor][i_s][w_i] = 0;     // update ownership state
(*snodes1).def_state[t_hor][i_s][w_i] = 1;     // update default state
(*snodes1).w_state_swap(res_t_0.w_i_floor);    // update wealth transition state
}
}

(*rr1).interp_vw3(t_i, i_m, i_s);  // clean grid

}

duration = (clock() - start) / (double)CLOCKS_PER_SEC;
cout << "time elapsed: " << duration << "  lcount = " << (*rr2).lcount << endl;
}
}
*/