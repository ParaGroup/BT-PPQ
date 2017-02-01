/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file wlq.hpp
 *  \brief Window-level Query Stage
 *  
 *  Implementation of the Window-level Query (WLQ) stage, which is the
 *  second stage of the Pane Farming pattern. WLQ is a farm in which the
 *  emitter is responsible to schedule the received pane results form the first stage
 *  to the workers that merge them into the corresponding window results. The
 *  computations on different windows are independent and can be executed by
 *  distinct workers in parallel.
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
 *  May 2016
 */

#ifndef WLQ_H
#define WLQ_H

// include
#include <deque>
#include <ff/lb.hpp>
#include <ff/node.hpp>
#include <ff/farm.hpp>
#include <ff/mapping_utils.hpp>
#include <general.hpp>
#include <queries.hpp>

using namespace ff;
using namespace std;

/* This file provides the following classes:
 *  -WLQ_Worker: class that implements the worker of the WLQ stage;
 *  -WLQ_Emitter: class that implements the emitter of the WLQ stage (round-robin strategy);
 *  -WLQ_Emitter_OD: class that implements the emitter of the WLQ stage (on-demand strategy);
 *  -WLQ_Collector: class that implements the collector of the WLQ stage.
 */

// struct of a window task scheduled by the WLQ Emitter to the WLQ Workers
template <typename Pane_t, typename Win_t>
struct Win_Task_t {
	Pane_t *p; // pane result to be added to the window result
	Win_t *w; // window result to be updated
	size_t worker_id; // id (starting from 0) of the WLQ worker receiving this task
	/* If this flag is true, the window result must be transmitted
	   to the collector (finalized) after the completion of this task. */
	bool finalize;

	// constructor
	Win_Task_t(Pane_t *_p, Win_t *_w, size_t _w_id, bool _finalize=false): p(_p), w(_w), worker_id(_w_id), finalize(_finalize) {}
};

/*! 
 *  \class WLQ_Worker
 *  
 *  \brief Worker of the Window-level Query Stage
 *  
 *  Worker functionality of the WLQ stage of the Pane Farming. It is
 *  responsible to receive windows tasks, each one consisting in a pair
 *  of a window and a pane result: the pane result will be inserted into
 *  the window by updating the corresponding window result.
 *  
 *  This class is defined in \ref Pane_Farming/include/wlq.hpp
 */
template<typename Pane_t, typename Win_t>
#if !defined(WLQ_OD)
class WLQ_Worker: public ff_node_t<Win_Task_t<Pane_t, Win_t>, Win_t> {
#else
class WLQ_Worker: public ff_monode_t<Win_Task_t<Pane_t, Win_t>, Win_t> {
#endif
private:
	typedef Win_Task_t<Pane_t, Win_t> w_task_t;
	size_t rcv_wtasks = 0; // number of processed window tasks by the worker
	size_t rcv_win_finalized = 0; // number of finalized windows
	double avg_wtask_time = 0; // average time to process a window task by the worker (usec)
	double avg_win_panes = 0; // average number of pane instances per finalized window
	double avg_win_results = 0; // average number of result tuples per finalized window
	volatile ticks time_elapsed_svc = 0; // global ticks passed in the svc by the worker

public:
	// constructor
	WLQ_Worker() {}

	// destructor
	~WLQ_Worker() {}

	// svc_init method
	int svc_init() {
		return 0;
	}

	// return the total amount of ticks spent in the svc from the beginning of the execution
	ticks getsvcticks() { return time_elapsed_svc; }

	// return the total number of svc calls from the beginning of the execution
	size_t getnumtask() { return rcv_wtasks; }

	// svc_end method
	void svc_end() {
		// if LOG mode is enabled we write on the log file
#if defined(LOG)
		ostringstream stream;
 		stream << "****************WLQ Worker " << this->get_my_id() << "****************\n";
 		stream << "Received window tasks: " << rcv_wtasks << "\n";
 		stream << "Average number of results per window: " << avg_win_results << "\n";
 		stream << "Average number of pane instances per window: " << avg_win_panes << "\n";
 		stream << "Average processing time per window task (usec): " << avg_wtask_time << "\n";
    	LOG_WLQ_WRITE(stream);
#endif
	}

	// svc method
	Win_t *svc(w_task_t *w_task) {
		// get the initial time of the svc
		volatile ticks start_ticks = getticks();
		Win_t *w = w_task->w; // window to update
		Pane_t *p = w_task->p; // pane instance to add
		// update the window with the new pane instance
		w->addComputePane(*p);
		// deallocate the pane if it is not needed by any other window
		if(p->decreaseRefCounter() == 1) delete p;
		// get the final time of the svc
		volatile ticks end_ticks = getticks();
		double elapsed = FROM_TICKS_TO_USECS(end_ticks - start_ticks); // processing time of the window task
		// update statistics
		rcv_wtasks++;
		avg_wtask_time += (1/((double) rcv_wtasks)) * (elapsed - avg_wtask_time);
		time_elapsed_svc += (end_ticks - start_ticks);
		// if the task must be finalized, we emit it to the collector
		if(w_task->finalize) {
			rcv_win_finalized++;
			avg_win_panes += (1/((double) rcv_win_finalized)) * (w->getNoInstances() - avg_win_panes);
			avg_win_results += (1/((double) rcv_win_finalized)) * (w->getSize() - avg_win_results);
			w_task->p = nullptr; // the pane pointer is no longer meaningful
			Win_t *result_win = w_task->w;
#if !defined(WLQ_OD)
			delete w_task;
			return result_win;
#else
			ff_monode_t<w_task_t, Win_t>::ff_send_out_to((Win_t *) w_task, 0); // to the emitter to notify the availability to receive new tasks
			ff_monode_t<w_task_t, Win_t>::ff_send_out_to(result_win, 1); // to the collector
			return ff_monode_t<w_task_t, Win_t>::GO_ON;
#endif
		}
		else {
#if !defined(WLQ_OD)
			delete w_task;
			return ff_node_t<w_task_t, Win_t>::GO_ON;
#else
			ff_monode_t<w_task_t, Win_t>::ff_send_out_to((Win_t *) w_task, 0); // to the emitter to notify the availability to receive new tasks
			return ff_monode_t<w_task_t, Win_t>::GO_ON;
#endif
		}
	}
};

/*! 
 *  \class WLQ_Emitter
 *  
 *  \brief Emitter of the Window-level Query Stage
 *  
 *  Emitter functionality of the WLQ stage of the Pane Farming. It is
 *  responsible to receive the pane results from the PLQ stage and to
 *  distribute window tasks to the workers. This version uses a standard
 *  round-robin assignment, in which all the tasks of the same window are
 *  always scheduled to the same worker.
 *  
 *  This class is defined in \ref Pane_Farming/include/wlq.hpp
 */
template<typename Pane_t, typename Win_t>
class WLQ_Emitter: public ff_node_t<Pane_t, Win_Task_t<Pane_t, Win_t>> {
private:
	typedef Win_Task_t<Pane_t, Win_t> w_task_t;
	// struct of a window descriptor used for scheduling
	struct Win_Descr_t {
		Win_t *w; // pointer to the window
		size_t first_pane_id; // id of the first pane of the window (starting from 0)
		size_t wp; // number of panes of the window
		/* full[i] counts the number of panes instances that we still need
		   to receive to complete the i-th pane of the window. -1 if we have
		   not received any pane instances of the i-th pane yet. */
		int *full;

		// constructor
		Win_Descr_t(size_t _id, size_t _wp, size_t _first) {
			w = new Win_t(_id);
			wp = _wp;
			first_pane_id = _first;
			full = new int[wp];
			fill_n(full, wp, -1);
		}

		// destructor
		~Win_Descr_t() {
			// remember: not delete the window!
			delete full;
		}

		// method to check if the window must be finalized
		bool isComplete() {
			bool complete = true;
			for(int i=0; i<wp; i++) complete = complete && (full[i] == 0);
			return complete;
		}

		// method to update the counter of the received pane instances of the window
		void updateInstanceCount(Pane_t *p) {
			// find the local id of the pane
			size_t lp_id = p->getId() - first_pane_id;
			// if it is the first instance received for that pane
			if(full[lp_id] == -1) full[lp_id] = p->getNoInstances() - 1;
			// otherwise we decrease the counter
			else full[lp_id]--;
		}
	};
	size_t wp; // number of panes in the window
	size_t sp; // number of panes in the window slide
	size_t mw; // number of WLQ workers
	size_t first_win_id = 0; // id of the window at position 0 of the deque descrSet
	deque<Win_Descr_t *> descrSet; // deque of open window descriptors
	ff_loadbalancer * const lb; // load balancer object for the scheduling to workers
	volatile ticks last_sample_ticks; // time (in ticks) of the last sample triggering
	double avg_rho = 0; // average utilization factor of the WLQ stage
	size_t n_samples = 0; // number of sampling intervals
	vector<size_t> tasks_to_w; // tasks_to_w[i] is the number of window tasks scheduled to the i-th worker during the last sampling interval
	vector<ticks> time_ws; // vector of ticks for each worker (getsvcticks statistics)
	vector<size_t> tasks_ws; // vector of window tasks for each worker (getnumtask statistics)

	// method to compute the load balancing ratio of the WLQ
	double getLBRatio() {
		double min_time = numeric_limits<double>::max();
		double max_time = 0;
		double total_time[mw];
		// get the pointers to the WLQ workers
		const svector<ff_node*> &workers = lb->getWorkers();
		// get the working time for each worker
		for(size_t i=0; i<mw; i++) {
			// cast to the WLQ_Worker type
			WLQ_Worker<Pane_t, Win_t> *worker_casted = static_cast<WLQ_Worker<Pane_t, Win_t> *>(workers[i]);
			// get the total execution time of the i-th worker
			total_time[i] = FROM_TICKS_TO_USECS(worker_casted->getsvcticks());
			if(total_time[i] < min_time) min_time = total_time[i];
			if(total_time[i] > max_time) max_time = total_time[i];
		}
		return (1-(min_time/max_time)) * 100;
	}

	// function to compute the actual utilization factor of the WLQ stage
	inline double computeWLQUtilization(ff_loadbalancer * const lb) {
		// statistics per sampling interval (SAMPLING_SCHED_WLQ_USEC)
		double lambdas[mw]; // lambdas[i] is the amount of window tasks to the i-th worker
		double mus[mw]; // mus[i] is the ideal number of window tasks that the i-th worker can process
		double time_idle[mw]; // time_idle[i] is the length (usec) of the idle fraction of the sampling interval of the i-th worker
		double lambda = 0; // total amount of window tasks generated by the WLQ Emitter
		double avg_ts = 0; // estimated average processing time per window task (usec)
		double rho = 0; // utilization factor of the WLQ stage
		double tot_time_working = 0; // sum of the time intervals where the workers are effectively processing window tasks
		double tot_tasks = 0; // sum of the number of window tasks processed by the workers
		// compute lambdas and lambda
		for(size_t i=0; i<mw; i++) {
			lambdas[i] = tasks_to_w[i];
			lambda += tasks_to_w[i];
		}
		// get the pointers to the WLQ workers
		const svector<ff_node*> &workers = lb->getWorkers();
		for(size_t i=0; i<mw; i++) {
			WLQ_Worker<Pane_t, Win_t> *worker_casted = static_cast<WLQ_Worker<Pane_t, Win_t> *>(workers[i]);
			// total time (in ticks) passed in the svc by the i-th worker from the beginning of the execution
			ticks active_time_ticks = worker_casted->getsvcticks();
			// time passed in processing window tasks by the i-th worker during the last sampling interval
			double time_working = FROM_TICKS_TO_USECS(active_time_ticks - time_ws[i]);
			tot_time_working += time_working;
			// idle time of the i-th worker
			time_idle[i] = std::max(0.0, SAMPLING_SCHED_WLQ_USEC - time_working);
			//cout << "Worker " << i << " Total sample: " << SAMPLING_SCHED_WLQ_USEC/1000 << " ms, Idle time: " << time_idle[i]/1000 << " ms, Work time: " << time_working/1000 << " ms" << endl;
			// total number of window tasks processed by the i-th worker from the beginning of the execution
			int tasks = worker_casted->getnumtask();
			// number of window tasks processed by the i-th worker during the last sampling interval
			mus[i] = tasks - tasks_ws[i];
			tot_tasks += mus[i];
			// save the last statistics
			tasks_ws[i] = tasks;
			time_ws[i] = active_time_ticks;
		}
		// compute the estimated average processing time per window task
		avg_ts = tot_time_working / tot_tasks;
		// compute the ideal service rate of the WLQ stage
		for(size_t i=0; i<mw; i++) {
			mus[i] += (time_idle[i] / avg_ts);
			rho += (lambdas[i] * lambdas[i]) / (lambda * mus[i]);
		}
		if(std::isnan(rho) || std::isinf(rho)) return 0;
		else return rho;
	}

public:
	// constructor
	WLQ_Emitter(size_t _wp, size_t _sp, size_t _mw, ff_loadbalancer * const _lb): wp(_wp), sp(_sp), mw(_mw), lb(_lb), tasks_to_w(mw, 0), time_ws(mw, 0), tasks_ws(mw, 0) {}

	// destructor
	~WLQ_Emitter() {}

	// svc_init method
	int svc_init() {
		// entering the sampling barrier
		pthread_barrier_wait(&sampleBarr);
		// initialize the time of the last sample triggering
		last_sample_ticks = getticks();
		return 0;
	}

	// svc_end method
	void svc_end() {
		// entering the printing barrier
		pthread_barrier_wait(&printBarr);
		// entering the precedence barrier 1
		pthread_barrier_wait(&precBarr1);
		// check whether there are ready windows still present in deque
		for(auto const& wd: descrSet) {
			if(wd != nullptr) {
				if(wd->isComplete()) MY_PRINT("WLQ_Emitter has some pending descriptors of ready windows in deque!\n");
				break;
			}
		}
		MY_PRINT("***********************WLQ Emitter Statistics***********************\n" <<\
			      "Average utilization (WLQ stage): " << avg_rho << "\n" <<\
			      "Load balancing ratio (WLQ stage): " << this->getLBRatio() << "%%\n" <<\
			      "********************************************************************\n");
		// entering the printing barrier 2
		pthread_barrier_wait(&precBarr2);
	}

	// svc method
	w_task_t *svc(Pane_t *p) {
		// check if a new sampling interval is triggered
		double elapsed = FROM_TICKS_TO_USECS(getticks() - last_sample_ticks);
		if(elapsed >= SAMPLING_SCHED_WLQ_USEC) {
			last_sample_ticks = getticks();
			// get the utilization factor
			double rho = computeWLQUtilization(lb);
			// update the average utilization factor of the WLQ stage
			if(rho < 50 && rho > 0) {
				n_samples++; // increase the number of completed samples
				avg_rho += (1/(double) n_samples) * (rho - avg_rho);
			}
			// reset the counters of scheduled tuples
			fill(tasks_to_w.begin(), tasks_to_w.end(), 0);
		}
		// get the id of the pane (starting from 0)
		size_t pane_id = p->getId();
		// determine the id of the first and last window that the pane belongs
		size_t first_win;
		if(pane_id < wp) first_win = 0;
		else first_win = ceil(((double) ((pane_id + 1) - wp)) / sp);
		size_t last_win = ceil(((double) (pane_id + 1)) / sp) - 1;
		// if needed, we have to allocate new window descriptors and store their pointers in the deque
		if(last_win - first_win_id >= descrSet.size()) {
			size_t old_size = descrSet.size();
			size_t new_size = last_win - first_win_id + 1;
			descrSet.resize(new_size);
			for(int i=old_size; i<new_size; i++) {
				// compute the id of the window to insert in the new descriptor
				size_t w_id = i + first_win_id;
				// compute the id of the first pane of the window
				size_t first_pane_id = sp * w_id;
				descrSet[i] = new Win_Descr_t(w_id, wp, first_pane_id);
			}
		}
		// for all the windows to which the pane p belongs
		for(size_t i=first_win - first_win_id; i<=last_win - first_win_id; i++) {
			/* We do the following actions:
			   1- create the corresponding task;
			   2- schedule the task to a worker by setting the finalized flag properly.
			*/
			if(descrSet[i] == nullptr) abort();
			// update the number of instances of the window
			descrSet[i]->updateInstanceCount(p);
			// create the window task to be scheduled
			bool toFinalize = descrSet[i]->isComplete();
			size_t toWorker = (descrSet[i]->w)->getId() % mw; // round-robin assignment
			tasks_to_w[toWorker]++;
			w_task_t *w_task = new w_task_t(p, descrSet[i]->w, toWorker, toFinalize);
			while(!lb->ff_send_out_to(w_task, toWorker));
			// if the window must be finalized we can destroy the descriptor
			if(toFinalize) {
				delete descrSet[i];
				descrSet[i] = nullptr; // nullptr corresponds to a finalized window
			}
		}
		/* Check the presence of a contiguous set of window descriptors,
		   starting from the beginning of the deque, to be deleted. */
		auto it = descrSet.begin();
		for(; it<descrSet.end(); it++) {
			if((*it) == nullptr) first_win_id++;
			else break;
		}
		descrSet.erase(descrSet.begin(), it);
		return ff_node_t<Pane_t, w_task_t>::GO_ON;
	}
};

/*! 
 *  \class WLQ_Emitter_OD
 *  
 *  \brief Emitter of the Window-level Query Stage
 *  
 *  Emitter functionality of the WLQ stage of the Pane Farming. It is
 *  responsible to receive pane results from the PLQ stage and to distribute
 *  window tasks to the worker. This version uses an on-demand scheduling
 *  of window tasks to the available workers by exploiting proper feedback
 *  queues from the WLQ workers to the WLQ emitter.
 *  
 *  This class is defined in \ref Pane_Farming/include/wlq.hpp
 */
template<typename Pane_t, typename Win_t>
class WLQ_Emitter_OD: public ff_node_t<Pane_t, Win_Task_t<Pane_t, Win_t>> {
private:
	typedef Win_Task_t<Pane_t, Win_t> w_task_t;
	// struct of a window descriptor used for scheduling
	struct Win_Descr_t {
		Win_t *w; // pointer to the window
		size_t first_pane_id; // id of the first pane of the window (starting from 0)
		size_t wp; // number of panes of the window
		/* full[i] counts the number of panes instances that we still need
		   to receive to complete the i-th pane of the window. -1 if we have
		   not received any pane instances of the i-th pane yet. */
		int *full;
		vector<Pane_t *> pendingPanes; // queue (actually a vector) of pointers to pending pane instances of this window
		bool isBusy = false; // if true a window task of this window is currently in execution in one of the WLQ workers

		// constructor
		Win_Descr_t(size_t _id, size_t _wp, size_t _first) {
			w = new Win_t(_id);
			wp = _wp;
			first_pane_id = _first;
			full = new int[wp];
			fill_n(full, wp, -1);
		}

		// destructor
		~Win_Descr_t() {
			// remember: not delete the window!
			delete full;
		}

		// method to check if the window must be finalized
		bool isComplete() {
			bool complete = true;
			for(int i=0; i<wp; i++) complete = complete && (full[i] == 0);
			return complete;
		}

		// method to update the counter of the received pane instances of the window
		void updateInstanceCount(Pane_t *p) {
			// find the local id of the pane
			size_t l_id = p->getId() - first_pane_id;
			// if it is the first pane instance received for that pane
			if(full[l_id] == -1) full[l_id] = p->getNoInstances() - 1;
			// otherwise we decrease the counter
			else full[l_id]--;
		}
	};
	size_t wp; // number of panes in each window
	size_t sp; // number of panes in the window slide
	size_t mw; // number of WLQ workers
	size_t plq_workers; // number of PLQ workers
	size_t last_w; // id of the last WLQ worker served
	size_t first_win_id = 0; // id of the first window in the deque descrSet
	deque<Win_Descr_t *> descrSet; // deque of open window descriptors
	size_t n_pendingWindows = 0; // number of windows with at least one pending pane instance
	bool *workerReady; // workerReady[i] = true means that the i-th WLQ worker is ready to receive a new window task
	size_t n_readyWorkers; // number of ready workers
	ff_loadbalancer * const lb; // load balancer object for the scheduling to workers
	size_t eos_received = 0; // number of EOS messages received
	volatile ticks last_sample_ticks; // time (in ticks) of the last sample triggering
	double avg_rho = 0; // average utilization factor of the WLQ stage
	size_t n_samples = 0; // number of sampling intervals
	unsigned long gen_tasks_sample = 0; // number of tasks generated during the last sampling interval
	vector<size_t> tasks_to_w; // tasks_to_w[i] is the number of window tasks scheduled to the i-th worker during the last sampling interval
	vector<ticks> time_ws; // vector of ticks for each worker (getsvcticks statistics)
	vector<size_t> tasks_ws; // vector of window tasks for each worker (getnumtask statistics)

	// method to select a ready worker for a new window task, -1 if no worker is ready
	int selectReadyWorker(int active_workers=-1) {
		// if active_workers is not specified we use all the WLQ workers
		if(active_workers == -1) active_workers = mw;
		for(size_t i=last_w+1; i<active_workers; i++) {
			if(workerReady[i]) {
				last_w = i;
				return i;
			}
		}
		for(size_t i=0; i<=last_w; i++) {
			if(workerReady[i]) {
				last_w = i;
				return i;
			}
		}
		return -1;
	}

	// method to send a new window task to a selected worker
	void sendTaskToWorker(Pane_t *p, size_t lw_id, size_t worker_id) {
		tasks_to_w[worker_id]++;
		// update the number of pane instances of the window
		descrSet[lw_id]->updateInstanceCount(p);
		// create the window task to be scheduled
		bool toFinalize = descrSet[lw_id]->isComplete();
		w_task_t *w_task = new w_task_t(p, descrSet[lw_id]->w, worker_id, toFinalize);
		while(!lb->ff_send_out_to(w_task, worker_id));
		// if the window must be finalized we can destroy the descriptor
		if(toFinalize) {
			delete descrSet[lw_id];
			descrSet[lw_id] = nullptr; // nullptr corresponds to a finalized window
		}
	}

	// method to compute the load balancing ratio of the WLQ
	double getLBRatio() {
		double min_time = numeric_limits<double>::max();
		double max_time = 0;
		double total_time[mw];
		// get the pointers to the WLQ workers
		const svector<ff_node*> &workers = lb->getWorkers();
		// get the working time for each worker
		for(size_t i=0; i<mw; i++) {
			// cast to the WLQ_Worker type
			WLQ_Worker<Pane_t, Win_t> *worker_casted = static_cast<WLQ_Worker<Pane_t, Win_t> *>(workers[i]);
			// get the total execution time of the i-th worker
			total_time[i] = FROM_TICKS_TO_USECS(worker_casted->getsvcticks());
			if(total_time[i] < min_time) min_time = total_time[i];
			if(total_time[i] > max_time) max_time = total_time[i];
		}
		return (1-(min_time/max_time)) * 100;
	}

	// function to compute the actual utilization factor of the WLQ stage
	inline double computeWLQUtilization(ff_loadbalancer * const lb) {
		// statistics per sampling interval (SAMPLING_SCHED_WLQ_USEC)
		double probs[mw]; // probs[i] is the probability to schedule a window task to the i-th worker
		double mus[mw]; // mus[i] is the ideal number of window tasks that the i-th worker can process
		double time_idle[mw]; // time_idle[i] is the length (usec) of the idle fraction of the sampling interval of the i-th worker
		double avg_ts = 0; // estimated average processing time per window task (usec)
		double rho = 0; // utilization factor of the WLQ stage
		double tot_time_working = 0; // sum of the time intervals where the workers are effectively processing window tasks
		double tot_tasks = 0; // sum of the number of window tasks processed by the workers
		// compute probs
		double tot_tasks_scheduled = 0;
		for(size_t i=0; i<mw; i++) tot_tasks_scheduled += tasks_to_w[i];
		for(size_t i=0; i<mw; i++) probs[i] = tasks_to_w[i] / tot_tasks_scheduled;
		// get the pointers to the WLQ workers
		const svector<ff_node*> &workers = lb->getWorkers();
		for(size_t i=0; i<mw; i++) {
			WLQ_Worker<Pane_t, Win_t> *worker_casted = static_cast<WLQ_Worker<Pane_t, Win_t> *>(workers[i]);
			// total time (in ticks) passed in the svc by the i-th worker from the beginning of the execution
			ticks active_time_ticks = worker_casted->getsvcticks();
			// time passed in processing window tasks by the i-th worker during the last sampling interval
			double time_working = FROM_TICKS_TO_USECS(active_time_ticks - time_ws[i]);
			tot_time_working += time_working;
			// idle time of the i-th worker
			time_idle[i] = std::max(0.0, SAMPLING_SCHED_WLQ_USEC - time_working);
			//cout << "Worker " << i << " Total sample: " << SAMPLING_SCHED_WLQ_USEC/1000 << " ms, Idle time: " << time_idle[i]/1000 << " ms, Work time: " << time_working/1000 << " ms" << endl;
			// total number of window tasks processed by the i-th worker from the beginning of the execution
			int tasks = worker_casted->getnumtask();
			// number of window tasks processed by the i-th worker during the last sampling interval
			mus[i] = tasks - tasks_ws[i];
			tot_tasks += mus[i];
			// save the last statistics
			tasks_ws[i] = tasks;
			time_ws[i] = active_time_ticks;
		}
		// compute the estimated average processing time per window task
		avg_ts = tot_time_working / tot_tasks;
		// compute the ideal service rate of the WLQ stage
		for(size_t i=0; i<mw; i++) {
			mus[i] += time_idle[i] / avg_ts;
			rho += probs[i] * ((gen_tasks_sample * probs[i]) / (mus[i]));
		}
		if(std::isnan(rho) || std::isinf(rho)) return 0;
		else return rho;
	}

public:
	// constructor
	WLQ_Emitter_OD(size_t _wp, size_t _sp, size_t _mw, size_t _plq_workers, ff_loadbalancer * const _lb): wp(_wp), sp(_sp), mw(_mw), last_w(_mw), plq_workers(_plq_workers), n_readyWorkers(_mw), lb(_lb), tasks_to_w(mw, 0), time_ws(mw, 0), tasks_ws(mw, 0) {
		workerReady = new bool[mw];
		fill_n(workerReady, mw, true);
	}

	// destructor
	~WLQ_Emitter_OD() {
		delete workerReady;
	}

	// svc_init method
	int svc_init() {
		// entering the sampling barrier
		pthread_barrier_wait(&sampleBarr);
		// initialize the time of the last sample triggering
		last_sample_ticks = getticks();
		return 0;
	}

	// svc_end method
	void svc_end() {
		// entering the printing barrier
		pthread_barrier_wait(&printBarr);
		// entering the precedence barrier 1
		pthread_barrier_wait(&precBarr1);
		// check whether there are ready windows still present in deque
		for(auto const& wd: descrSet) {
			if(wd != nullptr) {
				if(wd->isComplete()) MY_PRINT("WLQ_Emitter_OD has some pending descriptors of ready windows in deque!\n");
				break;
			}
		}
		MY_PRINT("***********************WLQ Emitter Statistics***********************\n" <<\
			      "Average utilization (WLQ stage): " << avg_rho << "\n" <<\
			      "Load balancing ratio (WLQ stage): " << this->getLBRatio() << "%%\n" <<\
			      "********************************************************************\n");
		// entering the precedence barrier 2
		pthread_barrier_wait(&precBarr2);
	}

	// svc method
	w_task_t *svc(Pane_t *p) {
		// check if a new sampling interval is triggered
		double elapsed = FROM_TICKS_TO_USECS(getticks() - last_sample_ticks);
		if(elapsed >= SAMPLING_SCHED_WLQ_USEC) {
			last_sample_ticks = getticks();
			// get the utilization factor
			double rho = computeWLQUtilization(lb);
			// update the average utilization factor of the WLQ stage
			if(rho < 50 && rho > 0) {
				n_samples++; // increase the number of completed samples
				avg_rho += (1/(double) n_samples) * (rho - avg_rho);
			}
			// reset the counters per sampling interval
			fill(tasks_to_w.begin(), tasks_to_w.end(), 0);
			gen_tasks_sample = 0;
		}
		// identify the source of the input message
		int msg_source = lb->get_channel_id();
        // the input message is a new pane instance from the PLQ stage
        if(msg_source == -1) {
        	// get the id of the pane (starting from 0)
			size_t pane_id = p->getId();
			// determine the id of the first and last window that the pane belongs
			size_t first_win;
			if(pane_id < wp) first_win = 0;
			else first_win = ceil(((double) ((pane_id + 1) - wp)) / sp);
			size_t last_win = ceil(((double) (pane_id + 1)) / sp) - 1;
			// if needed, we have to allocate new window descriptors and store their pointers in the deque
			if(last_win - first_win_id >= descrSet.size()) {
				size_t old_size = descrSet.size();
				size_t new_size = last_win - first_win_id + 1;
				descrSet.resize(new_size);
				for(int i=old_size; i<new_size; i++) {
					// compute the id of the window to insert in the new descriptor
					size_t w_id = i + first_win_id;
					// compute the id of the first pane of the window
					size_t first_pane_id = sp * w_id;
					descrSet[i] = new Win_Descr_t(w_id, wp, first_pane_id);
				}
			}
 			// for all the windows to which the pane p belongs
			for(size_t i=first_win - first_win_id; i<=last_win - first_win_id; i++) {
				/* We do the following actions:
				   1- if the window is not busy, we check the presence of an available worker:
				      1.1- if present, we schedule a new task to that worker and the window becomes busy;
				      1.2- otherwise, we add the pane instance to the pending queue of the window;
				   2- otherwise we add the pane instance to the pending queue of the window.
				*/
				gen_tasks_sample++;
				if(descrSet[i] == nullptr) abort();
				// if the window is not busy
				if(!descrSet[i]->isBusy) {
					// get a ready worker
					size_t worker_id = selectReadyWorker();
					// if an available worker exists
					if(worker_id != -1) {
						workerReady[worker_id] = false;
						n_readyWorkers--;
						// schedule a new task to this worker
						descrSet[i]->isBusy = true;
						sendTaskToWorker(p, i, worker_id);
					}
					// if no available worker exists
					else {
						(descrSet[i]->pendingPanes).push_back(p);
						if((descrSet[i]->pendingPanes).size() == 1) n_pendingWindows++;
					}
				}
				// the window is busy
				else {
					(descrSet[i]->pendingPanes).push_back(p);
					if((descrSet[i]->pendingPanes).size() == 1) n_pendingWindows++;
				}
			}
			/* Check the presence of a contiguous set of window descriptors,
			   starting from the beginning of the deque, to be deleted. */
			auto it = descrSet.begin();
			for(; it<descrSet.end(); it++) {
				if((*it) == nullptr) first_win_id++;
				else break;
			}
			descrSet.erase(descrSet.begin(), it);
			return ff_node_t<Pane_t, w_task_t>::GO_ON;
        }
        // the input message is a feedback from one of the WLQ workers
        else {
        	/* We do the following actions:
				1- if there are other pending pane instances of that window, we assign
				   a new task to the same worker;
				2- otherwise, we check the presence of a not busy window with some pending
				   pane instances:
				   2.1- if it exists, we schedule a new task to the worker and we mark the window as busy;
				   2.2- otherwise, the worker is marked as available to receive a new task.
			*/
        	w_task_t *fb = reinterpret_cast<w_task_t *>(p); // ugly but needed!
        	// get the local index in deque of the window of the task
        	size_t lw_id = (fb->w)->getId() - first_win_id;
        	if((!fb->finalize) && ((descrSet[lw_id]->pendingPanes).size() > 0)) {
        		// get the last pending pane instance
        		Pane_t *p = (descrSet[lw_id]->pendingPanes).back();
        		(descrSet[lw_id]->pendingPanes).pop_back();
        		if((descrSet[lw_id]->pendingPanes).size() == 0) n_pendingWindows--;
        		//send a new task to the same worker
        		sendTaskToWorker(p, lw_id, fb->worker_id);
        	}
        	// check if there is a not busy window with some pending pane instances to execute
        	else {
        		if((!fb->finalize)) descrSet[lw_id]->isBusy = false;
        		bool served = false;
        		for(auto const &descr: descrSet) {
        			if(descr != nullptr) {
        				if(!descr->isBusy && (descr->pendingPanes).size() > 0) {
        					descr->isBusy = true;
        					// get the last pending pane instance
        					Pane_t *p = (descr->pendingPanes).back();
        					(descr->pendingPanes).pop_back();
        					if((descr->pendingPanes).size() == 0) n_pendingWindows--;
        					// send a new task to the same worker
        					lw_id = (descr->w)->getId() - first_win_id;
        					sendTaskToWorker(p, lw_id, fb->worker_id);
							served = true;
							break;
        				}
        			}
        		}
        		// otherwise the worker is marked as available
        		if(!served) {
        			workerReady[fb->worker_id] = true;
        			n_readyWorkers++;
        		}
        	}
        	delete fb;
        	// check the termination
        	if((eos_received == plq_workers) && (n_pendingWindows == 0) && (n_readyWorkers == mw)) return ff_node_t<Pane_t, w_task_t>::EOS;
        	else return ff_node_t<Pane_t, w_task_t>::GO_ON;
        }
	}

	// method to properly manage the EOS
	void eosnotify(ssize_t id) {
		// we have to receive all EOS from the PLQ workers
		if(id == -1) {
            eos_received++;
            /* We propagate the EOS if we have received all the EOS from the previous stage,
               there is no pending pane instance in any window, and all the workers are ready .*/
            if((eos_received == plq_workers) && (n_pendingWindows == 0) && (n_readyWorkers == mw)) {
                lb->broadcast_task(EOS);
            }
        }
    }
};

/*! 
 *  \class WLQ_Collector
 *  
 *  \brief Collector of the Window-level Query Stage
 *  
 *  Collector functionality of the WLQ stage of the Pane Farming. It is
 *  in charge of receiving the window results from the WLQ_Worker instances
 *  and producing them outside in the increasing order of their identifier.
 *  
 *  This class is defined in \ref Pane_Farming/include/wlq.hpp
 */
template<typename Win_t>
class WLQ_Collector: public ff_node_t<Win_t> {
private:
	// identifier of the next window to transmit
	size_t next_win_id = 0; // window id starting from 0
	// deque of the windows received and buffered for ordering
	deque<Win_t *> winSet;

public:
	// constructor
	WLQ_Collector() {}

	// destructor
	~WLQ_Collector() {}

	// svc_init method
	int svc_init() {
		return 0;
	}

	// svc_end method
	void svc_end() {
		// check whether there are some pending windows in deque
		if(winSet.size() != 0) {
			MY_PRINT("WLQ_Collector has some pending windows in deque!\n");
		}
	}

	// svc method
	Win_t *svc(Win_t *w) {
#if !defined(UNORDERING_COLLECTOR)
		// if the received window is the next one to transmit
		if(w->getId() == next_win_id) {
			next_win_id++;
			ff_node_t<Win_t>::ff_send_out(w);
			if(winSet.size() > 0) {
				// check whether there are other windows to transmit
				auto it = winSet.begin() + 1;
				for(; it<winSet.end(); it++) {
					if((*it) == nullptr) break;
					else {
						next_win_id++;
						ff_node_t<Win_t>::ff_send_out(*it);
					}
				}
				// erase the pointers to the transmitted windows in deque
				winSet.erase(winSet.begin(), it);
			}
		}
		// otherwise the window must be buffered
		else {
			if(w->getId() - next_win_id >= winSet.size()) {
				size_t old_size = winSet.size();
				size_t new_size = w->getId() - next_win_id + 1;
				winSet.resize(new_size, nullptr);
			}
			winSet[w->getId() - next_win_id] = w;
		}
		return ff_node_t<Win_t>::GO_ON;
#else
		// simply forward the window to the next stage in the arrival order
		return w;
#endif
	}
};

#endif
