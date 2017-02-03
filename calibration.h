// calibration.h
// Includes parameters used throughout the program
// Copyright A. Michael Sharifi, 2016
// 
// cities: San Diego, San Francisco, Los Angeles, Boston, Chicago, Denver,  Miami, New York
// city_id: 0,1,2,3,4,5,6,7,8
// *.read_in files begin at year 5 (2007) and end at year 11 (2013)


// MOD HERE *******************************
const double loan_amt = 1.0;
const double apr_frm = 0.06;
const double mort_term = 30.0;
const int mort_term_int = 3; //30;

// state variable: mortgage rates
const double unemp_mult = 0.5; // Proportion of income received if unemployment shock

const double mti0 = 0.3;        // Mortgage-to-income at origination

const int m_n = 2; // 2 mortgage states: i_m = 0 for no mortgage, i_m = 1 for a mortgage

const int n_plevel = 3;
const int n_urate = 3;
const int n_fedfunds = 1;

const int city_begin = 2;
const int city_end = 2;

const int t_begin = 11;    // MOD HERE: begin in final time period, compute this only
const int t_end = 11;                   

// MOD HERE: Fix param_id = 0 for manual selection
const int param_id = 0;  // set = 0 to define parameters here manually; set = 1, 2, 3, 4 for presets and load in main

const int age_begin_store[] = { 30, 30, 30 }; // manual age settings here
const int n_age_store[] = { 1, 1, 1, 2, 2 };
const int n_age = n_age_store[param_id]; 

const double csfLevStore[] = {1.0/0.055, 1.0/0.055, 0.0, 1.0/0.055, 0.0}; // manually set futures leverage
const double csfLev = csfLevStore[param_id];
const int w_n = 100; // 100; // Grid points in wealth; set = 200 for fast computation, = 2000 for precision

const int age_max = 31; //35; //65;                  // age at which household retires / annuitizes wealth  

const int csfLevi = int(floor(csfLev));   // Floor for identification

// MOD HERE: set = 0 for renting, = 1 for owning
const int t_n = 2; // 5;                        // possible tenure states
const int pref = 0;                       // set pref = 0 for Cobb-Douglas, = 1 for CES
const int N_control = 6;
const int N_cities = 8;                    // number of cities

const int n_ph = 9; //9      // possible home price states
const int n_rent = 1; // 3;  possible rent states
const int n_yi = 1; // 3; //3; // 3;  // labor income states

const int n_s = n_ph * n_plevel * n_urate * n_fedfunds;   // impose fedfunds determined by plevel, urate

// Labor income related parameters
const double maint_mult = 0.98;
const double y_tax = 0.0;                 // 0.3;  // taxation is handled in snodes.cpp
const double y_atax = 1.0 - y_tax;
const double y_replace = 0.9388;            // From Cocco, Gomes, Maenhout (2005)

const double w_max = 40.0; //16.05;         // maximum wealth (on grid) (100's thousands)             
const double w_min = -2.0; // 0.0; // 0.05;          // minimum wealth (on grid) (100's thousands) 

const int w_i_zero = (int)ceil(-w_min * double(w_n) / (w_max - w_min));


const double rho = 5.0;
//const double rho = 1.8;                  // Power: curvature parameter; governs risk-aversion  also = 1.0, 2.0, 4.0
const int rhoi = int(floor(rho));        // Floor for identification

const double beta = .97;   // alt: =.95               // Beta: time preferences
const double phi = 0.08;                  // moving / transaction costs in event of home sale; also =.15
const double phi_sell = 0.08; //0.10;
const double phi_buy = 0.00;  //0.02;

//double var_a[] = { 0.0025, 0.0395, 0.0474 };                       // VAR: constants

// Housing-service related parameters
// median square footage by city
//const double hu_med[N_cities] = { 1.6, 1.7, 1.5, 1.9, 1.8, 1.6, 1.5, 1.83 };   // san fran updated

const double hu_med[N_cities] = { 1.6, 1.7, 1.5, 1.9, 1.8, 1.6, 1.5, 1.83 };   // san fran updated


// MOD HERE: assume representative house
const double hu_ten_store[N_cities][t_n] = {
	{ 1.0, 1.0 },
	{ 1.0, 1.0 },
	{ 1.0, 1.0 },
	{ 1.0, 1.0 },
	{ 1.0, 1.0 },
	{ 1.0, 1.0 },
	{ 1.0, 1.0 },
	{ 1.0, 1.0 }
};

const double hu_ten_def =  .5;  // square footage in default case

// If Cobb-Douglas Preferences:
const double alpha_cd = 0.7;                         // Prefence weight for C: Non-durable consumption
const double calpha_cd = 0.3;                        // =1.0 - alpha_sd; Preference weight for H: Housing services

const double rb = 1.0304;                             // Gross return on bonds / mortgage rate; sometimes = 1.04

// IF CES Preferences:
const double alpha_ces = -6.485; // -6.7; // .75;                // Low substitutability between C and H
const double gammad = 1.01;
const int gammai = int(floor(gammad));

//Equity approximation (Gaussian-Hermite quadrature, 2 node approximation
const double x_mu = 0.04 + (rb - 1.0); // Cocco, Gomes Maenhout (2005)
const double x_std = 0.157; 
const int retxn = 2;
const double retxv[] = { x_mu - x_std, x_mu + x_std }; //{ -0.10, 0.25 };
const double retxp[] = { 0.5, 0.5 };

const double p_move = 0.0;  // also: set = 0.0        // Probability receiving an exogenous moving shock

const double b_min_const = -20.0;                    
const double b_motive = 1.0;                         // Strength of bequest motive

const double c_fs = .06;                             // Minimum baseline consumption (Gov Assistance)
const double coh_fs = .05;                           // Cash on hand (Gov Asssistance: Non-durable Consumption + Housing)

// down-payment criteria
const double delta = .20;                            // Minimum down payment
const double min_dpmt =  .20;    //0.2;                     // minimum down payment
const double max_ltv = 0.80; // .95;  0.8;                        // max loan to value

// mortgage risk criteria
const double max_lti = 0.3;
const double mort_spread = 0.00;                       // mortgage spread above risk-free rate
const double pmi_dpmt = 0.00;                          // if down payment below this amount, add to mortgage spread
const double pmi_prem = 0.0; //0.01;                         // pmi premium
const double credit_prem = .15;                       // unsecured credit apr
const double b_min_unsec = -2.0; //-1.0;               // unsecured borrowing limit

const double num_small = -1.0e20;
