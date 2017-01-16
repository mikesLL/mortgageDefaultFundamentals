// class def_stats.h
// Contains code to compute default statistics                                     
// Copyright A. Michael Sharifi, 2016

#ifndef DEF_STATS_H
#define DEF_STATS_H

class def_stats {
	snodes *snodes1;       // pointer to snodes 
	//int t_hor;             // current planning period


public:

	def_stats(void *snodes_in);

	//void wtrans( int t_hor_in );  // compute wealth transition and hazard rate
	void wtrans_iterate(int t_hor_in);  // compute wealth transition and hazard rate

	vector<double> sdist, wdist, wdist2, sdist2; // state, wealth distributions

	vector<double> hazard_store;                       // hazard rate in each year
	vector<vector<double>> sdist_store, wdist_store;  // sdist, wdist in each year

	void print_def_stats( int t_hor_in );


};

#endif

