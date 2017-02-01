/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file kslack.hpp
 *  \brief Punctuation and Reordering Buffer Mechanisms based on K-Slack
 *  
 *  We implement two punctuation mechanisms. The first is based on the standard
 *  K-slack while the second one adopts an adaptive K-slack approach controller by a PID. The
 *  idea is to generate punctuations. A punctuation with timestamp t means that
 *  almost all the tuples with application timestamp lower than t have been received.
 *  
 *  We also implement the standard and adaptive (PID based) K-slack reordering buffers,
 *  where tuples are buffered and emitted amost in order.
 *  
 *  The two different approaches (punctuation mechanisms and reordering buffers) are
 *  adopted by the OOP and IOP versions of PLQ respectively)
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
 *  March 2016
 */

#ifndef KSLACK_H
#define KSLACK_H

// include
#include <deque>
#include <vector>
#include <algorithm>
#include <pid.hpp>
#include <general.hpp>
#include <tuple_t.hpp>

using namespace std;

/* This file provides the following classes:
 *  -KSlack_Punc: class that implements a punctuation mechanism based on standard K-slack;
 *  -KSlackPID_Punc: class that implements a punctuation mechanism based on adaptive K-slack
 *   with a PID control law to adapt the delay used by the algorithm at runtime;
 *  -KSlack_Buffer: class that implements a reordering buffer of tuples based on standard K-slack;
 *  -KSlackPID_Buffer: class that implements a reordering buffer of tuples based on adaptive K-slack
 *  with a PID control law to adapt the delay used by the algorithm at runtime.
 */

/*! 
 *  \class KSlack_Punc
 *  
 *  \brief Punctuation Mechanism based on Standard K-slack
 *  
 *  We implement a mechanism to generate punctuations according to the K-slack
 #  algorithm. K-slack is a well-known algorithm to re-order input tuples before
 *  presenting them to the computation. In this case we use this algorithm to
 #  derive a stream progress notion: a punctuation with timestamp t means that
 *  "almost" all the tuples with application timestamp < t have been arrived and
 *  already transmitted to the operator.
 *  
 *  This class is defined in \ref Pane_Farming/include/kslack.hpp
 */
class KSlack_Punc {
private:
	double K = 0; // K parameter of the slack (of the same time unit of the application timestamps)
	double tcurr = 0; // highest application timestamp of the tuples seen so far
	vector<double> ts_vect; // vector of the last application timestamps of the tuples

public:
	// constructor
	KSlack_Punc() {}

	// destructor
	~KSlack_Punc() {}

	// method to update the K-slack. It returns true in the case of a new punctuation
	inline bool checkProgress(tuple_t &t1, double &punctuation) {
		// insert the application timestamp of the tuple in the vector
		ts_vect.push_back(t1.app_ts);
		if(t1.app_ts <= tcurr) return false;
		else {
			tcurr = t1.app_ts; // update tcurr
			// find the maximum delay
			double max_d = 0;
			for (auto const &ts: ts_vect) {
    			double diff = tcurr-ts;
    			if(max_d < diff) max_d = diff;
			}
			if(max_d > K) K = max_d; // update K;
			/* All the tuples with application timestamp lower
			than tcurr-K have been likely received */
			double new_punc = tcurr - K;
			ts_vect.clear(); // empty the vector of application timestamps
			// emit the punctuation
			punctuation = new_punc;
			return true;
		}
	}
};

/*! 
 *  \class KSlackPID_Punc
 *  
 *  \brief Punctuation Mechanism based on Adaptive K-slack Controlled by a PID
 *  
 *  This class is an adaptive implementation of the punctuation mechanism based on K-slack
 *  in which the delay used is adapted according to the actual dropping probability. We
 *  recall that all the tuples with application timestamp t greater than the timestamp of
 *  the last punctuation are dropped. The idea is the following: try to adapt the K-slack
 *  delay length in order to keep the dropping probability close to a user-specified value.
 *  To achieve this goal we use a PID (Proportional-Integrative-Derivative) control law.
 *  
 *  This class is defined in \ref Pane_Farming/include/kslack.hpp
 */
class KSlackPID_Punc {
private:
	double K = 0; // K parameter of the slack (of the same time unit of the application timestamps)
	double tcurr = 0; // highest application timestamp of the tuples seen so far
	vector<double> ts_vect; // vector of the last application timestamps of the tuples
	double alpha = 0; // attenuation parameter of the K-slack delay (it is the output of the PID)
	PID *pid; // PID controller used to adapt the K-slack delay
	bool autoTuning; // flag: true if the auto-tuning mechanism is enabled, false otherwise
	PIDAutotuning *pid_autotuning; // PID auto-tuning mechanism for setting the gains of the PID
	volatile ticks last_sample_ticks; // time (in ticks) of the last sample triggering

public:
	// constructor
	KSlackPID_Punc(double _drop_target=0.1, bool _autoTuning=false): autoTuning(_autoTuning) {
		// create the controller
		double Kp = 6; double Ki = 0.6; double Kd = 0.2; // <----------(ad hoc values - to be checked)
		double coverage_setpoint = 1-_drop_target;
		pid = new PID(Kp, Ki, Kd, coverage_setpoint, -1, 1, true, true /* isReversed=true */);
		if(autoTuning) pid_autotuning = new PIDAutotuning(coverage_setpoint, true /* isReversed=true */);
		// initialize the time of the last sample triggering
		last_sample_ticks = getticks();
	}

	// destructor
	~KSlackPID_Punc() {
		delete pid;
		if(autoTuning) delete pid_autotuning;
	}

	// method to update the K-slack. It returns true in the case of a new punctuation
	inline bool checkProgress(tuple_t &t1, double &punctuation) {
		// insert the application timestamp of the tuple in the vector
		ts_vect.push_back(t1.app_ts);
		if(t1.app_ts <= tcurr) return false;
		else {
			tcurr = t1.app_ts; // update tcurr
			// find the maximum delay
			double max_d = 0;
			for (auto const &ts: ts_vect) {
    			double diff = tcurr-ts;
    			if(max_d < diff) max_d = diff;
			}
			if(max_d > K) K = max_d; // update K;
			/* All the tuples with application timestamp lower
			than tcurr-K have been likely received */
			double new_punc = tcurr - ((1-alpha) * K); // the real delay is (1-alpha)*K
			ts_vect.clear(); // empty the vector of application timestamps
			// emit the punctuation
			punctuation = new_punc;
			return true;
		}
	}

	// method to update the alpha parameter using the PID. It returns true if a new sampling interval is triggered
	inline bool adapt(double measured_drop_prob) {
		// check if a new sampling interval is triggered
		double elapsed = FROM_TICKS_TO_USECS(getticks() - last_sample_ticks);
		if(elapsed >= SAMPLING_PUNCTUATION_USEC) {
			last_sample_ticks = getticks();
			double measured_coverage = 1 - measured_drop_prob;
			double output;
			// if auto-tuning is enabled
			if(autoTuning) {
				pid_autotuning->setInput(measured_coverage);
				if(pid_autotuning->tuning()) {
					// the auto-tuning phase is complete
					autoTuning = false;
					// use the new gains obtained by auto-tuning
					cout << "Auto-tuning gains  = > Kp: " << pid_autotuning->getKp() << " Ki: " << pid_autotuning->getKi() << " Kd: " << pid_autotuning->getKd() << endl;
					pid->setGains(pid_autotuning->getKp(), pid_autotuning->getKi(), pid_autotuning->getKd());
				}
				output = pid_autotuning->getOutput();
			}
			// if auto-tuning is not enabled or it has been completed
			else {
				pid->setInput(measured_coverage);
				if(!pid->process()) cout << "PID is disabled" << endl;
				output = pid->getOutput();
			}
			// application of the output to the system
			alpha = alpha + output;
			// boundaries alpha is in [0,1]
			if(alpha > 1) alpha = 1;
			if(alpha < 0) alpha = 0;
			return true;
		}
		else return false;
	}
};

/*! 
 *  \class KSlack_Buffer
 *  
 *  \brief Reordering Buffer of Tuples based on Standard K-slack
 *  
 *  This class implements a classic reordering buffer of tuples based
 *  on the K-slack algorithm. This buffer is used to implement a IOP
 *  (In-Order-Processing) version of the Pane Farming pattern (used
 *  for comparison).
 *  
 *  This class is defined in \ref Pane_Farming/include/kslack.hpp
 */
class KSlack_Buffer {
private:
	double K = 0; // K parameter of the slack (of the same time unit of the application timestamps)
	double tcurr = 0; // highest application timestamp of the tuples seen so far
	deque<tuple_t *> bufferedTuples; // buffer of tuples waiting to be emitted
	vector<double> ts_vect; // vector of the last application timestamps of the tuples
	bool toEmit = false; // flag: it is true if some tuples must be emitted, false otherwise
	std::deque<tuple_t *>::iterator firstTuple; // iterator to the first tuple to emit
	std::deque<tuple_t *>::iterator lastTuple; // iterator to the last tuple to emit

public:
	// constructor
	KSlack_Buffer() {}

	// destructor
	~KSlack_Buffer() {}

	// method to insert a tuple into the buffer
	bool insertTuple(tuple_t *t) {
		// insert the application timestamp of the tuple in the vector
		ts_vect.push_back(t->app_ts);
		// insertion of the tuple in the buffer
		auto it = std::lower_bound(bufferedTuples.begin(), bufferedTuples.end(), t,
								   [](tuple_t *t1, tuple_t *t2)
								   { return t1->app_ts < t2->app_ts; });
		// if the tuple must be inserted at the end of the buffer
		if(it == bufferedTuples.end()) bufferedTuples.push_back(t);
		// if the tuple must be inserted in the middle of the buffer
		else bufferedTuples.insert(it, t);
		// check if we can emit some buffered tuples or not
		if(t->app_ts <= tcurr) return false;
		else {
			tcurr = t->app_ts; // update tcurr
			// find the maximum delay
			double max_d = 0;
			for (auto const &ts: ts_vect) {
    			double diff = tcurr-ts;
    			if(max_d < diff) max_d = diff;
			}
			if(max_d > K) K = max_d; // update K;
			// determine the iterator to the last tuple to emit
			tuple_t *tmp = new tuple_t();
			tmp->app_ts = tcurr - K;
			ts_vect.clear(); // empty the vector of application timestamps
			firstTuple = bufferedTuples.begin();
			lastTuple = std::lower_bound(bufferedTuples.begin(), bufferedTuples.end(), tmp,
								         [](tuple_t *t1, tuple_t *t2)
								         { return t1->app_ts <= t2->app_ts; });
			delete tmp;
			toEmit = true;
			return true;
		}
	}

	// method to get the tuples in the correct order from the buffer
	tuple_t *extractTuple() {
		if(!toEmit) return NULL;
		else if(firstTuple == lastTuple) {
			// remove the emitted tuples from the buffer
			bufferedTuples.erase(bufferedTuples.begin(), lastTuple);
			toEmit = false;
			return NULL;
		}
		else {
			tuple_t *t = *firstTuple;
			firstTuple++;
			return t;
		}
	}

	// method to prepare the flush of the buffer
	void prepareToFlush() {
		toEmit = true;
		firstTuple = bufferedTuples.begin();
		lastTuple = bufferedTuples.end();
	}

	// method to get the size of the buffer
	size_t getSize() const {
		return bufferedTuples.size();
	}
};

/*! 
 *  \class KSlackPID_Buffer
 *  
 *  \brief Reordering Buffer of Tuples based on Adaptive K-slack Controlled by a PID
 *  
 *  This class implements an adaptive reordering buffer of tuples based
 *  on the K-slack algorithm controlled by a PID. This buffer is used to
 *  implement a IOP (In-Order-Processing) version of the Pane Farming pattern
 *  (used for comparison).
 *  
 *  This class is defined in \ref Pane_Farming/include/kslack.hpp
 */
class KSlackPID_Buffer {
private:
	double K = 0; // K parameter of the slack (of the same time unit of the application timestamps)
	double tcurr = 0; // highest application timestamp of the tuples seen so far
	deque<tuple_t *> bufferedTuples; // buffer of tuples waiting to be emitted
	vector<double> ts_vect; // vector of the last application timestamps of the tuples
	bool toEmit = false; // flag: it is true if some tuples must be emitted, false otherwise
	std::deque<tuple_t *>::iterator firstTuple; // iterator to the first tuple to emit
	std::deque<tuple_t *>::iterator lastTuple; // iterator to the last tuple to emit
	double alpha = 0; // attenuation parameter of the K-slack delay (it is the output of the PID)
	PID *pid; // PID controller used to adapt the K-slack delay
	bool autoTuning; // flag: true if the auto-tuning mechanism is enabled, false otherwise
	PIDAutotuning *pid_autotuning; // PID auto-tuning mechanism for setting the gains of the PID
	volatile ticks last_sample_ticks; // time (in ticks) of the last sample triggering

public:
	// constructor
	KSlackPID_Buffer(double _drop_target=0.1, bool _autoTuning=false): autoTuning(_autoTuning) {
		// create the controller
		double Kp = 6; double Ki = 0.6; double Kd = 0.2; // <----------(ad hoc values - to be checked)
		double coverage_setpoint = 1-_drop_target;
		pid = new PID(Kp, Ki, Kd, coverage_setpoint, -1, 1, true, true /* isReversed=true */);
		if(autoTuning) pid_autotuning = new PIDAutotuning(coverage_setpoint, true /* isReversed=true */);
		// initialize the time of the last sample triggering
		last_sample_ticks = getticks();
	}

	// destructor
	~KSlackPID_Buffer() {
		delete pid;
		if(autoTuning) delete pid_autotuning;
	}

	// method to insert a tuple into the buffer
	bool insertTuple(tuple_t *t) {
		// insert the application timestamp of the tuple in the vector
		ts_vect.push_back(t->app_ts);
		// insertion of the tuple in the buffer
		auto it = std::lower_bound(bufferedTuples.begin(), bufferedTuples.end(), t,
								   [](tuple_t *t1, tuple_t *t2)
								   { return t1->app_ts < t2->app_ts; });
		// if the tuple must be inserted at the end of the buffer
		if(it == bufferedTuples.end()) bufferedTuples.push_back(t);
		// if the tuple must be inserted in the middle of the buffer
		else bufferedTuples.insert(it, t);
		// check if we can emit some buffered tuples or not
		if(t->app_ts <= tcurr) return false;
		else {
			tcurr = t->app_ts; // update tcurr
			// find the maximum delay
			double max_d = 0;
			for (auto const &ts: ts_vect) {
    			double diff = tcurr-ts;
    			if(max_d < diff) max_d = diff;
			}
			if(max_d > K) K = max_d; // update K;
			// determine the iterator to the last tuple to emit
			tuple_t *tmp = new tuple_t();
			tmp->app_ts = tcurr - ((1-alpha) * K); // the real delay is (1-alpha)*K
			ts_vect.clear(); // empty the vector of application timestamps
			firstTuple = bufferedTuples.begin();
			lastTuple = std::lower_bound(bufferedTuples.begin(), bufferedTuples.end(), tmp,
								         [](tuple_t *t1, tuple_t *t2)
								         { return t1->app_ts <= t2->app_ts; });
			delete tmp;
			toEmit = true;
			return true;
		}
	}

	// method to update the alpha parameter using the PID. It returns true if a new sampling interval is triggered
	inline bool adapt(double measured_drop_prob) {
		// check if a new sampling interval is triggered
		double elapsed = FROM_TICKS_TO_USECS(getticks() - last_sample_ticks);
		if(elapsed >= SAMPLING_PUNCTUATION_USEC) {
			last_sample_ticks = getticks();
			double measured_coverage = 1 - measured_drop_prob;
			double output;
			// if auto-tuning is enabled
			if(autoTuning) {
				pid_autotuning->setInput(measured_coverage);
				if(pid_autotuning->tuning()) {
					// the auto-tuning phase is complete
					autoTuning = false;
					// use the new gains obtained by auto-tuning
					cout << "Auto-tuning gains  = > Kp: " << pid_autotuning->getKp() << " Ki: " << pid_autotuning->getKi() << " Kd: " << pid_autotuning->getKd() << endl;
					pid->setGains(pid_autotuning->getKp(), pid_autotuning->getKi(), pid_autotuning->getKd());
				}
				output = pid_autotuning->getOutput();
			}
			// if auto-tuning is not enabled or it has been completed
			else {
				pid->setInput(measured_coverage);
				if(!pid->process()) cout << "PID is disabled" << endl;
				output = pid->getOutput();
			}
			// application of the output to the system
			alpha = alpha + output;
			// boundaries alpha is in [0,1]
			if(alpha > 1) alpha = 1;
			if(alpha < 0) alpha = 0;
			return true;
		}
		else return false;
	}

	// method to get the tuples in the correct order from the buffer
	tuple_t *extractTuple() {
		if(!toEmit) return NULL;
		else if(firstTuple == lastTuple) {
			// remove the emitted tuples from the buffer
			bufferedTuples.erase(bufferedTuples.begin(), lastTuple);
			toEmit = false;
			return NULL;
		}
		else {
			tuple_t *t = *firstTuple;
			firstTuple++;
			return t;
		}
	}

	// method to prepare the flush of the buffer
	void prepareToFlush() {
		toEmit = true;
		firstTuple = bufferedTuples.begin();
		lastTuple = bufferedTuples.end();
	}

	// method to get the size of the buffer
	size_t getSize() const {
		return bufferedTuples.size();
	}
};

#endif
