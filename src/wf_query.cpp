/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file wf_query.cpp
 *  \brief Application of the Window Farming Pattern
 *  
 *  The application consists in a pipeline of two stages where the first one
 *  is the Window Farming pattern while the second is a Client node that
 *  collects the results received from the Window Farming. The Window Farming
 *  can execute several types of continuous queries depending whether some
 *  macros are enabled or not. To date, two queries are supported:
 *  1- Skyline Query: macro SKYLINE;
 *  2- Top-δ Dominant Skyline Query: macro TOPD.
 *  
 *  \note The application receives stream elements from the Generator through
 *  a TCP/IP socket.
 *  
 */

/* ***************************************************************************
 *  
 *  This program is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License version 3 as
 *  published by the Free Software Foundation.
 *  
 *  This program is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 *  License for more details.
 *  
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software Foundation,
 *  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *  
 ****************************************************************************
 */

/*  Author: Gabriele Mencagli <mencagli@di.unipi.it>
 *  July 2016
 */

// include
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <limits>
#include <iostream>
#include <ff/node.hpp>
#include <ff/mapper.hpp>
#include <ff/pipeline.hpp>
#include <wf.hpp>
#include <pane.hpp>
#include <general.hpp>
#include <tuple_t.hpp>
#include <queries.hpp>
#include <socket_utils.hpp>

using namespace std;
using namespace ff;

// set the types of the Pane and Window according to the query executed
#if defined(SKYLINE)
typedef SkyPane Pane_t;
typedef SkyWin Win_t;
#elif defined(TOPD)
typedef TopdPane Pane_t;
typedef TopdWin Win_t;
#endif

/*! 
 *  \class Client
 *  
 *  \brief Client class of the Window Farming
 *  
 *  Class that implements the generic client that receives completed
 *  windows from the Window Farming pattern.
 *  
 *  This class is defined in \ref Pane_Farming/include/wf_query.cpp
 */
class Client: public ff_node_t<Pane_t> {
private:
	volatile ticks last_print, start_ticks;
	unsigned long rcv_panes = 0; // number of received windows over the entire execution
	unsigned long rcv_instances = 0; // number of received window instances over the entire execution
	vector<int> cnt_ins_vect; // counters of instances to receive for each window
	vector<int> sizes_pane_vect; // number of result tuples for each window
	vector<double> starting_pane_vect; // starting_pane_vect[i] is the starting time in usec of the i-th window
	vector<double> ending_pane_vect; // ending_pane_vect[i] is the ending time in usec of the i-th window
	vector<double> closing_pane_vect; // closing_pane_vect[i] is the closing time in usec of the i-th window
	double pane_rate = 0; // number of window received per second over the entire execution
	double avg_ins_per_pane = 0; // average number of instances per window over the entire execution
	double tot_avg_pane_size = 0; // average number of result tuples per window over the entire execution
	double tot_Sn_pane_size = 0; // n*\sigma^{2}_{n} of result tuples per window over the entire execution
	double tot_avg_closing_delay = 0; // average closing delay (usec) of windows over the entire execution
	double tot_avg_pane_latency = 0; // average window latency (usec) over the entire execution
	// statistics of a sampling interval
	unsigned long panes_sample = 0; // number of received windows per sample
	unsigned long ins_sample = 0; // number of received window instances per sample
	double avg_ins_size = 0; // average number of result tuples per window instance measured during a sampling interval
	double Sn_ins = 0; // n*\sigma^{2}_{n} of result tuples per window instance measured during a sampling interval
	double avg_pane_size = 0; // average number of result tuples per windows (sum of the sizes of the instances of the window) measured during a sampling interval
	double avg_closing_delay = 0; // average closing delay (usec) of windows measured during a sampling interval
	double avg_pane_latency = 0; // average pane latency (usec) of windows measured during a sampling interval

public:
	// constructor
	Client(): cnt_ins_vect(RESERVE_V_SPACE, -1), sizes_pane_vect(RESERVE_V_SPACE, 0), starting_pane_vect(RESERVE_V_SPACE, numeric_limits<size_t>::max()), ending_pane_vect(RESERVE_V_SPACE, 0),
				  closing_pane_vect(RESERVE_V_SPACE, 0) {}

	// destructor
	~Client() {}

	// svc_init method
	int svc_init() {
		// entering the sampling barrier
		pthread_barrier_wait(&sampleBarr);
		start_ticks = getticks();
		last_print = start_ticks;
		return 0;
	}

	// svc_end method
	void svc_end() {
		// print statistics of the last results (if there are any)
		volatile double elapsed = FROM_TICKS_TO_USECS(getticks()-last_print);
		if(panes_sample > 0) {
			double time = ((double)(getticks()-start_ticks)/((long long)((unsigned long)(FREQ)*1000000))); // in seconds
			printf("Time: %6.3fs, recvd wins: %ld, results (win instance) [avg: %6.2f std: %6.2f], results (win) [avg: %6.2f], closing delay (win) [avg: %6.2f us], latency (win) [avg: %6.2f us]\n", time, panes_sample, avg_ins_size, sqrt(Sn_ins/ins_sample), avg_pane_size, avg_closing_delay, avg_pane_latency);
		}
		// compute the rate of windows computed per second
		pane_rate = ((double) rcv_panes) / (FROM_TICKS_TO_USECS(getticks() - start_ticks) / 1000000);
		// entering the printing barrier
		pthread_barrier_wait(&printBarr);
		// entering the precedence barrier 1
		pthread_barrier_wait(&precBarr1);
		MY_PRINT("*******************************Final Statistics*******************************\n" <<\
			      "Received: " << pane_rate << " wins/sec, total " << rcv_panes << " windows\n" <<\
			      "Average closing delay per window: " << tot_avg_closing_delay << " usec\n" <<\
			      "Average latency per window: " << tot_avg_pane_latency/1000 << " msec\n" <<\
			      "Results per window: avg " << tot_avg_pane_size << " std " << sqrt(tot_Sn_pane_size/rcv_panes) << "\n"\
			      "******************************************************************************\n");
	}

	// svc method
	Pane_t *svc(Pane_t *r) {
		// increase the number of window instances
		ins_sample++;
		rcv_instances++;
		// increase the vectors if needed
		size_t pane_id = r->getId();
		if(pane_id >= cnt_ins_vect.size()) {
			cnt_ins_vect.resize(pane_id+1, -1);
			sizes_pane_vect.resize(pane_id+1, 0);
			starting_pane_vect.resize(pane_id+1, numeric_limits<size_t>::max());
			ending_pane_vect.resize(pane_id+1, 0);
			closing_pane_vect.resize(pane_id+1, 0);
		}
		// if it is the first instance of pane_id we set the number of instances to receive
		if(cnt_ins_vect[pane_id] == -1) {
			cnt_ins_vect[pane_id] = r->getNoInstances();
			avg_ins_per_pane += (1/((double) rcv_panes+1)) * (r->getNoInstances() - avg_ins_per_pane);
		}
		cnt_ins_vect[pane_id]--;
		// update the size of pane_id
		sizes_pane_vect[pane_id] += r->getSize();
		// update the starting time of the window
		if(starting_pane_vect[pane_id] > FROM_TICKS_TO_USECS(r->getStartingTime()))
			starting_pane_vect[pane_id] = FROM_TICKS_TO_USECS(r->getStartingTime());
		// update the ending time of the window
		if(ending_pane_vect[pane_id] < FROM_TICKS_TO_USECS(r->getEndingTime()))
			ending_pane_vect[pane_id] = FROM_TICKS_TO_USECS(r->getEndingTime());
		// update the closing time of the window
		if(closing_pane_vect[pane_id] < FROM_TICKS_TO_USECS(r->getClosingTime()))
			closing_pane_vect[pane_id] = FROM_TICKS_TO_USECS(r->getClosingTime());
		// statistics of the window instance size
		double old_avg = avg_ins_size;
		avg_ins_size += (1/((double) ins_sample)) * (r->getSize() - avg_ins_size);
		Sn_ins += (r->getSize() - old_avg) * (r->getSize() - avg_ins_size);
		// a window is complete when we have received all its instances
		if(cnt_ins_vect[pane_id] == 0) {
			rcv_panes++;
			panes_sample++;
			// compute the average latency of the windows
			double latency = ending_pane_vect[pane_id] - starting_pane_vect[pane_id];
			avg_pane_latency += (1/((double) panes_sample)) * (latency - avg_pane_latency);		
			// compute the total average latency of windows
			tot_avg_pane_latency += (1/((double) rcv_panes)) * (latency - tot_avg_pane_latency);		
			// compute the average delay of the windows
			double closing_delay = closing_pane_vect[pane_id] - ending_pane_vect[pane_id];
			avg_closing_delay += (1/((double) panes_sample)) * (closing_delay - avg_closing_delay);
			// compute the total average closing delay of windows
			tot_avg_closing_delay += (1/((double) rcv_panes)) * (closing_delay - tot_avg_closing_delay);
			// compute the average window size
			avg_pane_size += (1/((double) panes_sample)) * (sizes_pane_vect[pane_id] - avg_pane_size);
			// compute the total average window size
			double old_tot_avg = tot_avg_pane_size;
			tot_avg_pane_size += (1/((double) rcv_panes)) * (sizes_pane_vect[pane_id] - tot_avg_pane_size);
			tot_Sn_pane_size += (sizes_pane_vect[pane_id] - old_tot_avg) * (sizes_pane_vect[pane_id] - tot_avg_pane_size);
		}
		volatile double elapsed = FROM_TICKS_TO_USECS(getticks()-last_print);
		if(elapsed > PRINT_RATE_USEC) {
			double time = ((double)(getticks()-start_ticks)/((long long)((unsigned long)(FREQ)*1000000))); // in seconds
			printf("Time: %6.3fs, recvd wins: %ld, results (win instance) [avg: %6.2f std: %6.2f], results (win) [avg: %6.2f], closing delay (win) [avg: %6.2f us], latency (win) [avg: %6.2f us]\n", time, panes_sample, avg_ins_size, sqrt(Sn_ins/ins_sample), avg_pane_size, avg_closing_delay, avg_pane_latency);
			panes_sample = avg_ins_size = Sn_ins = avg_closing_delay = avg_pane_latency = ins_sample = avg_pane_size = 0;
			last_print = getticks();
		}
		delete r;
		return ff_node_t<Pane_t>::GO_ON;
	}
};

// main
int main(int argc, char *argv[]) {
	int option = 0;
	size_t port = 10000;
	double drop_prob = 0;
	unsigned long win_len = 0; // window length (in ms)
	unsigned long win_slide = 0; // window slide (in ms)
	size_t nw = 1;
	// arguments from command line
	if(argc != 11) {
		cout << argv[0] << " -p <port> -d <dropping probability> -w <win length ms> -s <win slide ms> -n <workers>" << endl;
		exit(0);
	}
	while ((option = getopt(argc, argv, "p:d:w:s:n:")) != -1) {
		switch (option) {
			case 'p': port = atoi(optarg);
				break;
			case 'd': drop_prob = stod(optarg);
				break;
			case 'w': win_len = atoi(optarg);
				break;
			case 's': win_slide = atoi(optarg);
				break;
			case 'n': nw = atoi(optarg);
				break;
			default: {
				cout << argv[0] << " -p <port> -d <dropping probability> -w <win length ms> -s <win slide ms> -n <workers>" << endl;
				exit(0);
			}
        }
    }
    cout << "The frequency of the CPU is " << FREQ << " Mhz" << endl;
    // set the thread mapping (only for RePhrase)
#if defined(REPHRASE)
    const char worker_mapping[]="0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 1, \
                                 9, 17, 25, 33, 41, 49, 57, 65, 73, 81, 89, 97, 105, 113, 121, 129, 137, 145, 153, 2, 10, \
                                 18, 26, 34, 42, 50, 58, 66, 74, 82, 90, 98, 106, 114, 122, ,130, 138, 146, 154, 3, 11, \
                                 19, 27, 35, 43, 51, 59, 67, 75, 83, 91, 99, 107, 115, 123, 131, 139, 147, 155, 4, 12, 20, \
                                 28, 36, 44, 52, 60, 68, 76, 84, 92, 100, 108, 116, 124, 132, 140, 158, 156, 5, 13, 21, 29, \
                                 37, 45, 53, 61, 69, 77, 85, 93, 101, 109, 117, 125, 133, 141, 149, 157, 6, 14, 22, 30, 38, \
                                 46, 54, 62, 70, 78, 86, 94, 102, 110, 118, 126, 134, 142, 150, 158, 7, 15, 23, 31, 39, 47, \
                                 55, 63, 71, 79, 87, 95, 103, 111, 119, 127, 135, 143, 151, 159";
    threadMapper::instance()->setMappingList(worker_mapping);
#endif
    // if LOG mode is enabled we initialize the log files
#if defined(LOG)
    LOG_WF_INIT();
#endif
    // create a pipeline with unbounded queues
    ff::ff_pipeline pipe;
	// initialize the sampling barrier between WF_Emitter and Client
	pthread_barrier_init(&sampleBarr, NULL, 2);
	// initialize the printing barrier between WF_Emitter and Client
	pthread_barrier_init(&printBarr, NULL, 2);
	// initialize the precedence barrier 1 between WF_Emitter and Client
	pthread_barrier_init(&precBarr1, NULL, 2);
    WF<Pane_t> wf(win_len, win_slide, nw, drop_prob, port);
    pipe.add_stage(&wf); // add the WF pattern to the pipeline
    Client client;
    pipe.add_stage(&client); // add the Client node to the pipeline
    pipe.setFixedSize(false); // the pipeline uses unbounded queues!
#if defined(SKYLINE)
    cout << "Starting pipeline (Skyline Query)..." << endl;
#elif defined(TOPD)
    cout << "Starting pipeline (Top-δ Dominant Skyline Query)..." << endl;
#endif
    if(pipe.run_and_wait_end() < 0) {
    	cerr << "Error executing the pipeline" << endl;
    	return -1;
    }
    else {
    	cout << "...done" << endl;
    }
    // if LOG mode is enabled we close the log files
#if defined(LOG)
    LOG_WF_CLOSE();
#endif
    //pipe.ffStats(cout);
    return 0;
}
