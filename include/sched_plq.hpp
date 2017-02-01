/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file sched_plq.hpp
 *  \brief Scheduling Strategies used by the PLQ Stage
 *  
 *  Different implementations of the scheduling strategies adopted by the
 *  Emitter functionality of the PLQ stage of the Pane Farming pattern.
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
 *  April 2016
 */

#ifndef SCHED_H
#define SCHED_H

// include
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <ff/lb.hpp>
#include <ff/node.hpp>
#include <ff/utils.hpp>
#include <pid.hpp>
#include <tuple_t.hpp>
#include <general.hpp>

using namespace ff;
using namespace std;

/* This file provides the following classes:
 *  -PB_RR_Sched: pane-based round robin strategy;
 *  -TB_RR_Sched: tuple-based round robin scheduling strategy;
 *  -Adaptive_Sched: adaptive scheduling with splitting strategy of large panes;
 *  -Adaptive_Sched_PID: adaptive scheduling with splitting strategy of large panes and controlled by a PID.
 */

// function to compute the actual utilization factor of the PLQ stage
inline double computePLQUtilization(ff_loadbalancer * const lb, vector<size_t> &tuples, vector<ticks> &time_ws, vector<size_t> &tasks_ws, size_t n_punc) {
	size_t nw = tuples.size(); // parallelism degree
	// statistics per sampling interval (SAMPLING_SCHED_PLQ_USEC)
	double lambdas[nw]; // lambdas[i] is the amount of tuples/punctuations to the i-th worker
	double mus[nw]; // mus[i] is the ideal number of tuples/punctuations that the i-th worker can process
	double time_idle[nw]; // time_idle[i] is the length (usec) of the idle fraction of the sampling interval of the i-th worker
	double lambda = 0; // total amount of tuples/punctuations managed by the PLQ stage
	double avg_ts = 0; // estimated average processing time per tuple/punctuation (usec)
	double rho = 0; // utilization factor of the PLQ stage
	double tot_time_working = 0; // sum of the time intervals where the workers are effectively processing tuples/punctuations
	double tot_tasks = 0; // sum of the number of tuples/punctuations processed by the workers
	// compute lambdas and lambda
	for(size_t i=0; i<nw; i++) {
		lambdas[i] = tuples[i] + n_punc;
		lambda += tuples[i] + n_punc;
	}
	// get the pointers to the PLQ workers
	const svector<ff_node*> &workers = lb->getWorkers();
	for(size_t i=0; i<nw; i++) {
		// total time (in ticks) passed in the svc by the i-th worker from the beginning of the execution
		ticks active_time_ticks = workers[i]->getsvcticks();
		// time passed in processing tuples/punctuations by the i-th worker during the last sampling interval
		double time_working = FROM_TICKS_TO_USECS(active_time_ticks - time_ws[i]);
		tot_time_working += time_working;
		// idle time of the i-th worker
		time_idle[i] = std::max(0.0, SAMPLING_SCHED_PLQ_USEC - time_working);
		// total number of tuples/punctuations processed by the i-th worker from the beginning of the execution
		int tasks = workers[i]->getnumtask();
		// number of tuples/punctuations processed by the i-th worker during the last sampling interval
		mus[i] = tasks - tasks_ws[i];
		tot_tasks += mus[i];
		// save the last statistics
		tasks_ws[i] = tasks;
		time_ws[i] = active_time_ticks;
	}
	// compute the estimated average processing time per tuple
	avg_ts = tot_time_working / tot_tasks;
	// compute the ideal service rate of the PLQ stage
	for(size_t i=0; i<nw; i++) {
		mus[i] += time_idle[i] / avg_ts;
		rho += (lambdas[i] * lambdas[i]) / (lambda * mus[i]);
	}
	if(std::isnan(rho) || std::isinf(rho)) return 0;
	else return rho;
}

/*! 
 *  \class PB_RR_Sched
 *  
 *  \brief Pane-based Round-Robin Scheduling Strategy
 *  
 *  Strategy that assigns each new pane to a worker in a round-robin
 *  fashion.
 *  
 *  This class is defined in \ref Pane_Farming/include/sched_plq.hpp
 */
class PB_RR_Sched {
private:
	size_t nw; // number of workers
	ff_loadbalancer * const lb; // load balancer object
	vector<size_t> tuples_to_w; // tuples_to_w[i] is the number of tuples scheduled to the i-th worker
	vector<int> panes_ws; // panes_ws[i] is the id of the worker assigned to the i-th pane (-1 not assigned yet)
	vector<int> pane_sizes; // pane_sizes[i] is the size (no. of tuples) of the pane with identifier i
	size_t last_w=0; // id of the last worker assigned to a new pane
	vector<ticks> time_ws; // vector of ticks for each worker (getsvcticks statistics)
	vector<size_t> tasks_ws; // vector of tasks for each worker (getnumtask statistics)
	volatile ticks last_sample_ticks; // time (in ticks) of the last sample triggering
	size_t next_pane_id = 0; // identifier of the next pane to complete
	double avg_rho = 0; // average utilization factor of the PLQ stage
	size_t n_samples = 0; // number of sampling intervals
	size_t n_punc = 0; // number of punctuations generated by the PLQ Emitter in the current sampling interval
	bool iop; // true if the PLQ runs in IOP mode (otherwise we use the OOP model)

public:
	// constructor
	PB_RR_Sched(size_t _nw, ff_loadbalancer * const _lb, ticks _starting_time, bool _iop=false): nw(_nw), lb(_lb), last_sample_ticks(_starting_time), tuples_to_w(nw, 0), time_ws(nw, 0), tasks_ws(nw, 0), iop(_iop) {}

	// destructor
	~PB_RR_Sched() {}

	// method to schedule a new tuple to a worker
	inline void sched_to_worker(tuple_t *t) {
		// check if a new sampling interval is triggered
		double elapsed = FROM_TICKS_TO_USECS(getticks() - last_sample_ticks);
		if(elapsed >= SAMPLING_SCHED_PLQ_USEC) {
			last_sample_ticks = getticks();
			// get the utilization factor
			double rho = computePLQUtilization(lb, tuples_to_w, time_ws, tasks_ws, n_punc);
			// update the average utilization factor of the PLQ stage
			if(rho < 50 && rho > 0) {
				n_samples++; // increase the number of completed samples
				avg_rho += (1/(double) n_samples) * (rho - avg_rho);
			}
			// reset the counters of scheduled tuples
			fill(tuples_to_w.begin(), tuples_to_w.end(), 0);
			// reset the number of punctuations
			n_punc = 0;
		}
		// get the id of the pane
		size_t pane_id = t->pane_id;
		if(pane_id >= panes_ws.size()) {
			panes_ws.resize(pane_id+1, -1);
			pane_sizes.resize(pane_id+1, 0);
		}
		// increase the counter of the number of tuples per pane
		pane_sizes[pane_id]++;
		// if the pane has not been assigned to a worker yet
		if(panes_ws[pane_id] == -1) {
			// select the worker in a round-robin fashion
			size_t w = last_w;
			last_w = (last_w + 1) % nw;
			panes_ws[pane_id] = w;
			tuples_to_w[w]++;
			while(!lb->ff_send_out_to(t, w));
		}
		// the pane already has an assigned worker
		else {
			size_t w = panes_ws[pane_id];
			tuples_to_w[w]++;
			while(!lb->ff_send_out_to(t, w));
		}
	}

	// method to generate a punctuation
	inline void generatePunctuation(double ts, size_t p_id) {
		if(p_id < next_pane_id) return;
		// check for lulls
		for(int i=next_pane_id; i<p_id; i++) {
			if(pane_sizes[i] == 0) {
				// emit a special tuple indicating a lull (empty pane)
				tuple_t *t = new tuple_t();
				t->isLull = true;
				t->pane_id = i;
				t->no_instances = 1;
				while(!lb->ff_send_out_to(t, 0)); // always to worker 0 (little additional overhead for it)
			}
		}
		next_pane_id = p_id;
		// generate punctuation
		if(!iop) {
			n_punc++;
			// create the punctuation
			tuple_t *p = new tuple_t();
			p->set_P(ts, p_id, nw);
			// broadcast the punctuation to the workers
			lb->broadcast_task(p);
		}
	}

	// method to get the average utilization of the PLQ stage
	double getAvgUtilization() const {
		return avg_rho;
	}

	// method to get the average splitting factor
	double getSplitFactor() const {
		return 1;
	}

	// method to compute the load balancing ratio of the PLQ
	double getLBRatio() {
		double min_time = numeric_limits<double>::max();
		double max_time = 0;
		double total_time[nw];
		// get the pointers to the PLQ workers
		const svector<ff_node*> &workers = lb->getWorkers();
		// get the working time for each worker
		for(size_t i=0; i<nw; i++) {
			// get the total execution time of the i-th worker
			total_time[i] = FROM_TICKS_TO_USECS(workers[i]->getsvcticks());
			if(total_time[i] < min_time) min_time = total_time[i];
			if(total_time[i] > max_time) max_time = total_time[i];
		}
		return (1-(min_time/max_time)) * 100;
	}

	// method to obtain the mean and std of the pane size
	pair<double,double> getPaneSizeStats() {
		double avg = 0;
		double Sn = 0;
		int count = 0;
		for(auto const&s: pane_sizes) {
			count++;
			double old_avg = avg;
			avg += (1/((double) count)) * (s - avg);
			Sn += (s - old_avg) * (s - avg);
		}
		return pair<double,double>(avg, sqrt(Sn/count));
	}
};

/*! 
 *  \class TB_RR_Sched
 *  
 *  \brief Tuple-based Round-Robin Scheduling Strategy
 *  
 *  Strategy that assigns each new tuple to a worker selected in a
 *  round-robin fashion.
 *  
 *  This class is defined in \ref Pane_Farming/include/sched_plq.hpp
 */
class TB_RR_Sched {
private:
	size_t nw; // number of workers
	ff_loadbalancer * const lb; // load balancer object
	vector<size_t> tuples_to_w; // tuples_to_w[i] is the number of tuples scheduled to the i-th worker
	vector<size_t> p_instances; // p_instances[i] is the number of instances of the i-th pane
	vector<int> pane_sizes; // pane_sizes[i] is the size (no. of tuples) of the pane with identifier i
	vector<vector<size_t>> p_workers; // p_worker[i] is a vector containing the ids of the workers having an instance of the i-th pane
	size_t last_w=0; // id of the last worker assigned to a new tuple
	vector<ticks> time_ws; // vector of ticks for each worker (getsvcticks statistics)
	vector<size_t> tasks_ws; // vector of tasks for each worker (getnumtask statistics)
	volatile ticks last_sample_ticks; // time (in ticks) of the last sample triggering
	size_t next_pane_id = 0; // identifier of the next pane to complete
	double avg_rho = 0; // average utilization factor of the PLQ stage
	size_t n_samples = 0; // number of sampling intervals
	size_t n_punc = 0; // number of punctuations generated by the PLQ Emitter in the current sampling interval
	bool iop; // true if the PLQ runs in IOP mode (otherwise we use the OOP model)

public:
	// constructor
	TB_RR_Sched(size_t _nw, ff_loadbalancer * const _lb, ticks _starting_time, bool _iop=false): nw(_nw), lb(_lb), last_sample_ticks(_starting_time), tuples_to_w(nw, 0), time_ws(nw, 0), tasks_ws(nw, 0), iop(_iop) {}

	// destructor
	~TB_RR_Sched() {}

	// method to schedule a new tuple to a worker
	inline void sched_to_worker(tuple_t *t) {
		// check if a new sampling interval is triggered
		double elapsed = FROM_TICKS_TO_USECS(getticks() - last_sample_ticks);
		if(elapsed >= SAMPLING_SCHED_PLQ_USEC) {
			last_sample_ticks = getticks();
			// get the utilization factor
			double rho = computePLQUtilization(lb, tuples_to_w, time_ws, tasks_ws, n_punc);
			// update the average utilization factor of the PLQ stage
			if(rho < 50 && rho > 0) {
				n_samples++; // increase the number of completed samples
				avg_rho += (1/(double) n_samples) * (rho - avg_rho);
			}
			// reset the counters of scheduled tuples
			fill(tuples_to_w.begin(), tuples_to_w.end(), 0);
			// reset the number of punctuations
			n_punc = 0;
		}
		// get the id of the pane
		size_t pane_id = t->pane_id;
		if(pane_id >= p_instances.size()) {
			p_instances.resize(pane_id+1, 0);
			p_workers.resize(pane_id+1);
			pane_sizes.resize(pane_id+1, 0);
		}
		// increase the counter of the number of tuples per pane
		pane_sizes[pane_id]++;
		// get the id of the next worker in a round-robin fashion
		size_t w = last_w;
		last_w = (last_w + 1) % nw;
		// if the pane has less than nw instances and the selected worker does not still have an instance of pane_id
		if((p_instances[pane_id] < nw) && (find(p_workers[pane_id].begin(), p_workers[pane_id].end(), w) == p_workers[pane_id].end())) {
			// multicast the presence of a new instance to the interested workers
			p_instances[pane_id]++;
			multicast_instance(p_workers[pane_id], pane_id, p_instances[pane_id]);
			t->no_instances = p_instances[pane_id];
			// add the worker to the list of workers having an instance of pane_id
			p_workers[pane_id].push_back(w);
		}
		// send the tuple to the selected worker
		tuples_to_w[w]++;
		while(!lb->ff_send_out_to(t, w));
	}

	// method to generate a punctuation
	inline void generatePunctuation(double ts, size_t p_id) {
		if(p_id < next_pane_id) return;
		// check for lulls
		for(int i=next_pane_id; i<p_id; i++) {
			if(pane_sizes[i] == 0) {
				// emit a special tuple indicating a lull (empty pane)
				tuple_t *t = new tuple_t();
				t->isLull = true;
				t->pane_id = i;
				t->no_instances = 1;
				while(!lb->ff_send_out_to(t, 0)); // always to worker 0 (little additional overhead for it)
			}
		}
		next_pane_id = p_id;
		// generate punctuation
		if(!iop) {
			n_punc++;
			// create the punctuation
			tuple_t *p = new tuple_t();
			p->set_P(ts, p_id, nw);
			// broadcast the punctuation to the workers
			lb->broadcast_task(p);
		}
	}

	// method to multicast a new_pane_instance message to a subset of the workers
	inline void multicast_instance(vector<size_t> &ws, size_t pane_id, size_t count_instances) {
		// without destination workers we simply return without effects
		if(ws.size() == 0) return;
		// create the message
		tuple_t *t = new tuple_t();
		t->set_NewInstance(pane_id, count_instances, ws.size());
		std::vector<size_t> retry;
		for(auto const &w_id: ws) {
			if(!lb->ff_send_out_to(t, w_id)) retry.push_back(w_id);
		}
		while(retry.size()) {
			if(lb->ff_send_out_to(t, retry.back())) retry.pop_back();
			else ticks_wait(TICKS2WAIT); // wait a short time interval
   		}
	}

	// method to get the average utilization of the PLQ stage
	double getAvgUtilization() const {
		return avg_rho;
	}

	// method to get the average splitting factor
	double getSplitFactor() const {
		double avg = 0;
		for(auto const &c: p_instances) avg += c;
		avg = avg / p_instances.size();
		return avg;
	}

	// method to compute the load balancing ratio of the PLQ
	double getLBRatio() {
		double min_time = numeric_limits<double>::max();
		double max_time = 0;
		double total_time[nw];
		// get the pointers to the PLQ workers
		const svector<ff_node*> &workers = lb->getWorkers();
		// get the working time for each worker
		for(size_t i=0; i<nw; i++) {
			// get the total execution time of the i-th worker
			total_time[i] = FROM_TICKS_TO_USECS(workers[i]->getsvcticks());
			if(total_time[i] < min_time) min_time = total_time[i];
			if(total_time[i] > max_time) max_time = total_time[i];
		}
		return (1-(min_time/max_time)) * 100;
	}
};

/*! 
 *  \class Adaptive_Sched
 *  
 *  \brief Adaptive Scheduling Strategy
 *  
 *  First implementation of the adaptive scheduling strategy used by the PLQ
 *  stage of the Pane Farming pattern. The idea is to balance the workload
 *  by assigning to the workers roughly the same number of input tuples.
 *  To this end, if the actual pane instance size is greater than a historical
 *  measure, we split a pane into another instance assigned to a distinct worker.
 *  The goal is to split only the panes large enough.
 *  
 *  This class is defined in \ref Pane_Farming/include/sched_plq.hpp
 */
class Adaptive_Sched {
private:
	size_t nw; // number of workers
	ff_loadbalancer * const lb; // load balancer object
	vector<size_t> tuples_to_w; // number of tuples distributed per worker
	// Pane Descriptor: struct with useful information of a pane
	struct pane_des_t {
		int a_worker; // id of the worker assigned to the active instance of the pane
		size_t n_instances; // number of existing instances of the pane
		vector<size_t> w_ids; // ids of the workers that have an instance of the pane
		size_t a_tuples; // number of tuples of the active instance of the pane
		vector<size_t> i_sizes; // i_sizes[i] is the size of the instance assigned to worker w_ids[i]
		// default constructor
		pane_des_t(size_t _nw): i_sizes(_nw, 0) { a_worker = -1; n_instances = 0; a_tuples = 0; }
	};
	vector<pane_des_t> pv; // vector of pane descriptors
	vector<int> pane_sizes; // pane_sizes[i] is the size of the pane with identifier i
	size_t next_pane_id = 0; // identifier of the next pane to complete
	double avg_size = 0; // average size of the pane instances
	double Sn_size = 0; // n*\sigma_{n}^{2} of the pane instances size
	double std_size = 0; // standard deviation of the pane instances size
	double tot_instances = 0; // number of pane instances produced
	size_t last_w = 0; // id of the last worker (used by get_llw_id with queue lengths)
	vector<ticks> time_ws; // vector of ticks for each worker (getsvcticks statistics)
	vector<size_t> tasks_ws; // vector of tasks for each worker (getnumtask statistics)
	volatile ticks last_sample_ticks; // time (in ticks) of the last sample triggering
	double avg_rho = 0; // average utilization factor of the PLQ stage
	size_t n_samples = 0; // number of sampling intervals
	size_t n_punc = 0; // number of punctuations generated by the PLQ Emitter in the current sampling interval
	bool warmup = true; // true if we are in the initial warmup period in order to have good estimates of avg_size and std_size
	bool iop; // true if the PLQ runs in IOP mode (otherwise we use the OOP model)

public:
	// constructor
	Adaptive_Sched(size_t _nw, ff_loadbalancer * const _lb, ticks _starting_time, bool _iop=false): nw(_nw), lb(_lb), last_sample_ticks(_starting_time), tuples_to_w(nw, 0), time_ws(nw, 0), tasks_ws(nw, 0), iop(_iop) {}

	// destructor
	~Adaptive_Sched() {}

	// method to schedule a new tuple to a worker
	inline void sched_to_worker(tuple_t *t) {
		// check if a new sampling interval is triggered
		double elapsed = FROM_TICKS_TO_USECS(getticks() - last_sample_ticks);
		if(elapsed >= SAMPLING_SCHED_PLQ_USEC) {
			last_sample_ticks = getticks();
			// get the utilization factor
			double rho = computePLQUtilization(lb, tuples_to_w, time_ws, tasks_ws, n_punc);
			// update the average utilization factor of the PLQ stage
			if(rho < 50 && rho > 0) {
				n_samples++; // increase the number of completed samples
				avg_rho += (1/(double) n_samples) * (rho - avg_rho);
			}
			// reset the counters of scheduled tuples
			fill(tuples_to_w.begin(), tuples_to_w.end(), 0);
			// reset the number of punctuations
			n_punc = 0;
		}
		// create the pane descriptor if it does not exist
		size_t pane_id = t->pane_id;
		if(pane_id >= pv.size()) {
			pv.resize(pane_id+1, pane_des_t(nw));
			pane_sizes.resize(pane_id+1, 0);
		}
		// increase the counter of the number of tuples per pane
		pane_sizes[pane_id]++;
		// Case 1: the tuple is the first of a new pane
		if(pv[pane_id].a_tuples == 0) {
			// we schedule to the least loaded worker
			int w = get_llw_id();
			tuples_to_w[w]++;
			pv[pane_id].a_tuples = 1;
			pv[pane_id].a_worker = w;
			pv[pane_id].n_instances = 1;
			pv[pane_id].w_ids.push_back(w);
			(pv[pane_id].i_sizes)[w] = 1;
			// send the tuple to the selected worker
			while(!lb->ff_send_out_to(t, w));
		}
		// Case 2: the tuple belongs to an existing pane
		else {
			// if the pane istance is not so large or we are in the warmup period
			if((pv[pane_id].a_tuples <= (avg_size + std_size)) || warmup) {
				// schedule the tuple as before
				int w = pv[pane_id].a_worker;
				tuples_to_w[w]++;
				pv[pane_id].a_tuples++;
				(pv[pane_id].i_sizes)[w]++;
				while(!lb->ff_send_out_to(t, w));
			}
			// otherwise we create another instance of the pane
			else {
				// we find the least loaded worker
				int w = get_llw_id();
				// if this worker already has an instance, it becomes the active one again
				vector<size_t> &w_ids = pv[pane_id].w_ids;
				if(find(w_ids.begin(), w_ids.end(), w) != w_ids.end()) {
					tuples_to_w[w]++;
					pv[pane_id].a_tuples = (pv[pane_id].i_sizes)[w] + 1;
					pv[pane_id].a_worker = w;
					(pv[pane_id].i_sizes)[w]++;
					while(!lb->ff_send_out_to(t, w));
				}
				// else we create a new instance of the pane that becomes the active one
				else {
					tuples_to_w[w]++;
					pv[pane_id].a_tuples = 1;
					pv[pane_id].a_worker = w;
					pv[pane_id].n_instances++;
					(pv[pane_id].i_sizes)[w] = 1;
					// multicast the presence of a new instance to the interested workers
					multicast_instance(pv[pane_id].w_ids, pane_id, pv[pane_id].n_instances);
					pv[pane_id].w_ids.push_back(w);
					t->no_instances = pv[pane_id].n_instances;
					// send the tuple to the new selected worker
					while(!lb->ff_send_out_to(t, w));
				}
			}
		}
	}

	// method to generate a punctuation
	inline void generatePunctuation(double ts, size_t p_id) {
		if(p_id < next_pane_id) return;
		else {
			if(p_id >= pv.size()) {
				pv.resize(p_id+1, pane_des_t(nw));
				pane_sizes.resize(p_id+1, 0);
			}
		}
		// check for lulls
		for(int i=next_pane_id; i<p_id; i++) {
			if(pane_sizes[i] == 0) {
				// emit a special tuple indicating a lull (empty pane)
				tuple_t *t = new tuple_t();
				t->isLull = true;
				t->pane_id = i;
				t->no_instances = 1;
				while(!lb->ff_send_out_to(t, 0)); // always to worker 0 (little additional overhead for it)
			}
		}
		// update the statistics of the completed panes (avg size and std)
		for(int i=next_pane_id; i<p_id; i++) {
			for(auto const& s: pv[i].i_sizes) {
				if(s > 0) {
					tot_instances++;
					double old_avg_size = avg_size;
					avg_size += (((double) 1)/tot_instances) * (s - old_avg_size);
					Sn_size += (s - old_avg_size) * (s - avg_size);
				}
			}
		}
		std_size= sqrt(((double) Sn_size)/tot_instances);
		// we conclude the warmup period after the completion of the first five pane (heuristic)
		if(p_id > 5) warmup = false;
		next_pane_id = p_id;
		// generate punctuation
		if(!iop) {
			n_punc++;
			// create the punctuation
			tuple_t *p = new tuple_t();
			p->set_P(ts, p_id, nw);
			// broadcast the punctuation to the workers
			lb->broadcast_task(p);
		}
	}

	// method to obtain the id of the least loaded worker (llw)
	inline int get_llw_id() {
		// the llw is the worker with the shortest input queue
		size_t min = numeric_limits<size_t>::max();
		int llw_id = -1;
		const svector<ff_node*> &workers = lb->getWorkers();
		size_t w = last_w;
		for(size_t i=0; i<workers.size(); i++) {
			FFBUFFER *q = workers[w]->get_in_buffer();
			size_t len = q->length();
			if(len < min) {
				min = len;
				llw_id = w;
			}
			w = (w+1) % workers.size();
		}
		last_w = (llw_id + 1) % workers.size();
		return llw_id;
	}

	// method to multicast a new_pane_instance message to a subset of the workers
	inline void multicast_instance(vector<size_t> &ws, size_t pane_id, size_t count_instances) {
		// without destination workers we simply return without effects
		if(ws.size() == 0) return;
		// create the message
		tuple_t *t = new tuple_t();
		t->set_NewInstance(pane_id, count_instances, ws.size());
		std::vector<size_t> retry;
		for(auto const &w_id: ws) {
			if(!lb->ff_send_out_to(t, w_id)) retry.push_back(w_id);
		}
		while(retry.size()) {
			if(lb->ff_send_out_to(t, retry.back())) retry.pop_back();
			else ticks_wait(TICKS2WAIT); // wait a short time interval
   		}
	}

	// method to get the average utilization of the PLQ stage
	double getAvgUtilization() const {
		return avg_rho;
	}

	// method to get the average splitting factor
	double getSplitFactor() const {
		double avg = 0;
		for(auto const &d: pv) avg += d.n_instances;
		avg = avg / pv.size();
		return avg;
	}

	// method to compute the load balancing ratio of the PLQ
	double getLBRatio() {
		double min_time = numeric_limits<double>::max();
		double max_time = 0;
		double total_time[nw];
		// get the pointers to the PLQ workers
		const svector<ff_node*> &workers = lb->getWorkers();
		// get the working time for each worker
		for(size_t i=0; i<nw; i++) {
			// get the total execution time of the i-th worker
			total_time[i] = FROM_TICKS_TO_USECS(workers[i]->getsvcticks());
			if(total_time[i] < min_time) min_time = total_time[i];
			if(total_time[i] > max_time) max_time = total_time[i];
		}
		return (1-(min_time/max_time)) * 100;
	}
};

/*! 
 *  \class Adaptive_Sched_PID
 *  
 *  \brief Adaptive Scheduling Strategy with PID-based Control
 *  
 *  Second implementation of the adaptive scheduling strategy used by the PLQ
 *  stage of the Pane Farming pattern. This implementation is an extension of
 *  Adaptive_Sched in which we dynamically adapt the threshold value used to
 *  determine whether a pane instance must be divided or not. The threshold
 *  is modified whenever the  utilization is not closer to a desired value.
 *  This dynamic adaptation is achieved using a PID controller.
 *  
 *  This class is defined in \ref Pane_Farming/include/sched_plq.hpp
 */
class Adaptive_Sched_PID {
private:
	size_t nw; // number of workers
	ff_loadbalancer * const lb; // load balancer object
	vector<size_t> tuples_to_w; // number of tuples distributed per worker
	// Pane Descriptor: struct with useful information of a pane
	struct pane_des_t {
		int a_worker; // id of the worker assigned to the active instance of the pane
		size_t n_instances; // number of existing instances of the pane
		vector<size_t> w_ids; // ids of the workers that have an instance of the pane
		size_t a_tuples; // number of tuples of the active instance of the pane
		vector<size_t> i_sizes; // i_sizes[i] is the size of the instance assigned to worker w_ids[i]
		// default constructor
		pane_des_t(size_t _nw): i_sizes(_nw, 0) { a_worker = -1; n_instances = 0; a_tuples = 0; }
	};
	vector<pane_des_t> pv; // vector of pane descriptors
	vector<int> pane_sizes; // pane_sizes[i] is the size of the pane with identifier i
	size_t next_pane_id = 0; // identifier of the next pane to complete
	double avg_size = 0; // average size of the pane instances
	double Sn_size = 0; // n*\sigma_{n}^{2} of the pane instances size
	double std_size = 0; // standard deviation of the pane instances size
	double tot_instances = 0; // number of pane instances produced
	size_t last_w = 0; // id of the last worker (used by get_llw_id with queue lengths)
	vector<ticks> time_ws; // vector of ticks for each worker (getsvcticks statistics)
	vector<size_t> tasks_ws; // vector of tasks for each worker (getnumtask statistics)
	PID *pid; // PID controller used to dynamically adapt the threshold value
	double alpha = 1; // corrective factor of the threshold used to split the pane instances
	volatile ticks last_sample_ticks; // time (in ticks) of the last sample triggering
	double avg_rho = 0; // average utilization factor of the PLQ stage
	size_t n_samples = 0; // number of sampling intervals of PID
	size_t n_punc = 0; // number of punctuations generated by the PLQ Emitter in the current sampling interval
	bool warmup = true; // true if we are in the initial warmup period in order to have good estimates of avg_size and std_size
	bool iop; // true if the PLQ runs in IOP mode (otherwise we use the OOP model)

public:
	// constructor
	Adaptive_Sched_PID(size_t _nw, ff_loadbalancer * const _lb, ticks _starting_time, double _desired_rho, bool _iop=false): nw(_nw), lb(_lb), last_sample_ticks(_starting_time), tuples_to_w(nw, 0), time_ws(nw, 0), tasks_ws(nw, 0), iop(_iop) {
		// create the controller
		double Kp = 1.8; double Ki = 0; double Kd = 2; // <----------(ad hoc values - to be checked)
		pid = new PID(Kp, Ki, Kd, _desired_rho, -0.2, 0.2, true, true /* isReversed=true */);
	}

	// destructor
	~Adaptive_Sched_PID() {
		delete pid;
	}

	// method to schedule a new tuple to a worker
	inline void sched_to_worker(tuple_t *t) {
		// check if a new sampling interval is triggered
		double elapsed = FROM_TICKS_TO_USECS(getticks() - last_sample_ticks);
		if(elapsed >= SAMPLING_SCHED_PLQ_USEC) {
			last_sample_ticks = getticks();
			// get the utilization factor
			double rho = computePLQUtilization(lb, tuples_to_w, time_ws, tasks_ws, n_punc);
			// update the average utilization factor of the PLQ stage
			if(rho < 50 && rho > 0) {
				n_samples++; // increase the number of completed samples
				avg_rho += (1/(double) n_samples) * (rho - avg_rho);
			}
 			// compute the PID output to adapt the threshold
			pid->setInput(rho);
			if(!pid->process()) cout << "PID is disabled" << endl;
			double output = pid->getOutput();
			alpha += output;
			if(alpha < 0.01) alpha = 0.01;
			//cout << "Utilization: " << rho << " alpha: " << alpha << " threshold: " << (1/alpha) * (avg_size + std_size) << endl;
			fill(tuples_to_w.begin(), tuples_to_w.end(), 0); // reset the counters of scheduled tuples
			// reset the number of punctuations
			n_punc = 0;
		}
		// create the pane descriptor if it does not exist
		size_t pane_id = t->pane_id;
		if(pane_id >= pv.size()) {
			pv.resize(pane_id+1, pane_des_t(nw));
			pane_sizes.resize(pane_id+1, 0);
		}
		// increase the counter of the number of tuples per pane
		pane_sizes[pane_id]++;
		// Case 1: the tuple is the first of a new pane
		if(pv[pane_id].a_tuples == 0) {
			// we schedule to the least loaded worker
			int w = get_llw_id();
			tuples_to_w[w]++;
			pv[pane_id].a_tuples = 1;
			pv[pane_id].a_worker = w;
			pv[pane_id].n_instances = 1;
			pv[pane_id].w_ids.push_back(w);
			(pv[pane_id].i_sizes)[w] = 1;
			// send the tuple to the selected worker
			while(!lb->ff_send_out_to(t, w));
		}
		// Case 2: the tuple belongs to an existing pane
		else {
			// compute the actual threshold
			double threshold = (1/alpha) * (avg_size + std_size);
			// if the active pane instance is not large enough or we are in the warmup period
			if((pv[pane_id].a_tuples <= threshold) || warmup) {
				// schedule the tuple as before
				int w = pv[pane_id].a_worker;
				tuples_to_w[w]++;
				pv[pane_id].a_tuples++;
				(pv[pane_id].i_sizes)[w]++;
				while(!lb->ff_send_out_to(t, w));
			}
			// otherwise we create another instance of the pane
			else {
				// we find the least loaded worker
				int w = get_llw_id();
				// if this worker already has an instance, it becomes the active one again
				vector<size_t> &w_ids = pv[pane_id].w_ids;
				if(find(w_ids.begin(), w_ids.end(), w) != w_ids.end()) {
					tuples_to_w[w]++;
					pv[pane_id].a_tuples = (pv[pane_id].i_sizes)[w] + 1;
					pv[pane_id].a_worker = w;
					(pv[pane_id].i_sizes)[w]++;
					while(!lb->ff_send_out_to(t, w));
				}
				// else we create a new instance of the pane that becomes the active one
				else {
					tuples_to_w[w]++;
					pv[pane_id].a_tuples = 1;
					pv[pane_id].a_worker = w;
					pv[pane_id].n_instances++;
					(pv[pane_id].i_sizes)[w] = 1;
					// multicast the presence of a new instance to the interested workers
					multicast_instance(pv[pane_id].w_ids, pane_id, pv[pane_id].n_instances);
					pv[pane_id].w_ids.push_back(w);
					t->no_instances = pv[pane_id].n_instances;
					// send the tuple to the new selected worker
					while(!lb->ff_send_out_to(t, w));
				}
			}
		}
	}

	// method to generate a punctuation
	inline void generatePunctuation(double ts, size_t p_id) {
		if(p_id < next_pane_id) return;
		else {
			if(p_id >= pv.size()) {
				pv.resize(p_id+1, pane_des_t(nw));
				pane_sizes.resize(p_id+1, 0);
			}
		}
		// check for lulls
		for(int i=next_pane_id; i<p_id; i++) {
			if(pane_sizes[i] == 0) {
				// emit a special tuple indicating a lull (empty pane)
				tuple_t *t = new tuple_t();
				t->isLull = true;
				t->pane_id = i;
				t->no_instances = 1;
				while(!lb->ff_send_out_to(t, 0)); // always to worker 0
			}
		}
		// update the statistics of the completed pane (avg size and std)
		for(int i=next_pane_id; i<p_id; i++) {
			for(auto const& s: pv[i].i_sizes) {
				if(s > 0) {
					tot_instances++;
					double old_avg_size = avg_size;
					avg_size += (((double) 1)/tot_instances) * (s - old_avg_size);
					Sn_size += (s - old_avg_size) * (s - avg_size);
				}
			}
		}
		std_size= sqrt(((double) Sn_size)/tot_instances);
		next_pane_id = p_id;
		// we conclude the warmup period after the completion of the first five pane (heuristic)
		if(p_id > 5) warmup = false;
		// generate punctuation
		if(!iop) {
			n_punc++;
			// create the punctuation
			tuple_t *p = new tuple_t();
			p->set_P(ts, p_id, nw);
			// broadcast the punctuation to the workers
			lb->broadcast_task(p);
		}
	}

	// method to obtain the id of the least loaded worker (llw)
	inline int get_llw_id() {
		// the llw is the worker with the shortest input queue
		size_t min = numeric_limits<size_t>::max();
		int llw_id = -1;
		const svector<ff_node*> &workers = lb->getWorkers();
		size_t w = last_w;
		for(size_t i=0; i<workers.size(); i++) {
			FFBUFFER *q = workers[w]->get_in_buffer();
			size_t len = q->length();
			if(len < min) {
				min = len;
				llw_id = w;
			}
			w = (w+1) % workers.size();
		}
		last_w = (llw_id + 1) % workers.size();
		return llw_id;
	}

	// method to multicast a new_pane_instance message to a subset of the workers
	inline void multicast_instance(vector<size_t> &ws, size_t pane_id, size_t count_instances) {
		// without destination workers we simply return without effects
		if(ws.size() == 0) return;
		// create the message
		tuple_t *t = new tuple_t();
		t->set_NewInstance(pane_id, count_instances, ws.size());
		std::vector<size_t> retry;
		for(auto const &w_id: ws) {
			if(!lb->ff_send_out_to(t, w_id)) retry.push_back(w_id);
		}
		while(retry.size()) {
			if(lb->ff_send_out_to(t, retry.back())) retry.pop_back();
			else ticks_wait(TICKS2WAIT); // wait a short time interval
   		}
	}

	// method to get the average utilization of the PLQ stage
	double getAvgUtilization() const {
		return avg_rho;
	}

	// method to get the average splitting factor
	double getSplitFactor() const {
		double avg = 0;
		for(auto const &d: pv) avg += d.n_instances;
		avg = avg / pv.size();
		return avg;
	}

	// method to compute the load balancing ratio of the PLQ
	double getLBRatio() {
		double min_time = numeric_limits<double>::max();
		double max_time = 0;
		double total_time[nw];
		// get the pointers to the PLQ workers
		const svector<ff_node*> &workers = lb->getWorkers();
		// get the working time for each worker
		for(size_t i=0; i<nw; i++) {
			// get the total execution time of the i-th worker
			total_time[i] = FROM_TICKS_TO_USECS(workers[i]->getsvcticks());
			if(total_time[i] < min_time) min_time = total_time[i];
			if(total_time[i] > max_time) max_time = total_time[i];
		}
		return (1-(min_time/max_time)) * 100;
	}
};

#endif
