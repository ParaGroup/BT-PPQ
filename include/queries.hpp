/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file queries.hpp
 *  \brief Implementation of the Supported Continuous Preference Queries
 *  
 *  In this file we provide a set of classes for supporting various kinds of
 *  continuous preference queries. To date, the supported queries are:
 *  1- Skyline query: it returns the set of non-dominated points;
 *  2- Top-δ dominant skyline query: it returns the first δ>=1 points with the
 *     lowest dominance factor, where the dominance factor of a point p is the
 *     maximum k in [0, DIM] such that there exists a p' that k-dominates
 *     p (it is equal or better than p in at least k dimensions and better in
 *     at least one of these dimensions).
 *  
 *  Recall that not every preference query can be executed by the Pane
 *  Farming parallel pattern.
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

#ifndef QUERIES_H
#define QUERIES_H

// include
#include <vector>
#include <pane.hpp>
#include <window.hpp>
#include <general.hpp>
#include <tuple_t.hpp>

using namespace std;

/* This file provides the following classes:
 *  -TData: class that implements a special container of tuples (column-based layout);
 *  -SkyPane: class of the pane implementation for skyline queries;
 *  -SkyWin: class of the window implementation for skyline queries;
 *  -TopdPane: class of the pane implementation for top-δ dominant skyline query;
 *  -TopdWin: class of the window implementation for top-δ dominant skyline query.
 */

/*! 
 *  \class TData
 *  
 *  \brief Tuple Container
 *  
 *  Container for the efficient storage of tuples for the query processing.
 *  It uses a column-oriented memory layout.
 *  
 *  This class is defined in \ref Pane_Farming/include/queries.hpp
 */
class TData {
public:
	// index of the last valid tuple (-1 the container is empty)
	int last_idx = -1;
	// vectors for the tuples' coordinates + timestamps
	std::vector<double> v_dim[DIM+2];

    // constructor
    TData() {
        // reserve space for the vectors
        for(size_t d=0; d<DIM+2; d++) v_dim[d].reserve(RESERVE_V_SPACE);
    }

    // destructor
    ~TData() {}

    // method to append a tuple at the end of the container
    inline void append(tuple_t const &t) {
        last_idx++;
        for(size_t d=0; d<DIM; d++) v_dim[d].push_back(t.d[d]);
        v_dim[DIM].push_back(t.app_ts);
        v_dim[DIM+1].push_back(t.gen_ts);
    }

    // method to remove all the tuples from position pos (included) to the end of the container
    inline void removeFromIndex(size_t pos) {
        // if position is not valid return;
        if(pos > last_idx) return;
        // if the tuple to remove is the last one
        if(pos == last_idx) {
            for(size_t d=0; d<DIM+2; d++) v_dim[d].pop_back();
        }
        else {
            // maintain only the tuples from the one with index 0 to pos-1
            for(size_t d=0; d<DIM+2; d++) v_dim[d].resize(pos);
        }
        last_idx=pos-1;
    }

	// method to return the number of tuples in the container
	inline size_t getSize() const {
		return (size_t) (last_idx+1);
	}
};

/*! 
 *  \class SkyPane
 *  
 *  \brief Pane Implementation for Skyline Queries
 *  
 *  Pane implementation for continuous skyline queries.
 *  
 *  This class is defined in \ref Pane_Farming/include/queries.hpp
 */
class SkyPane: public Pane {
private:
	// pointer to the tuple container
	TData *s;

public:
	// constructor I
	SkyPane(size_t w_cnt=1, size_t _id=0, size_t _inst=0, size_t _len=0): Pane(w_cnt, _id, _inst, _len) {
		// create the TData container
		s = new TData();
	}

	// destructor
	~SkyPane() {
		delete s;
	}

	/* 
	 * Method to add a tuple to the pane and incrementally update the internal
	 * pane result.
	 */
	inline void addComputeTuple(tuple_t const &t) {
		// set the starting time of the pane
		if(this->getStartingTime() == 0) this->setStartingTime(t.arrival_ticks);
		else if(this->getStartingTime() > t.arrival_ticks) this->setStartingTime(t.arrival_ticks);
		int pos;
		bool isDominated = updateSkyline(t, pos);
		// pos is the index of the last point in the container not dominated by t
		if(pos < s->last_idx) {
			// remove all the points in the container dominated by t
			s->removeFromIndex(pos+1);
		}
		// if t does not dominate any point in the container and it is not dominated, we append it to the container
		if((pos == s->last_idx) && !isDominated) s->append(t);
		// set the ending time of the pane
		this->setEndingTime(getticks());
	}

	/* 
	 * Method to shift all the tuples dominated by t at the end of the TData container of the pane.
	 * The index of the last tuple not dominated by t (-1 if all the elements are dominated) is written
	 * in the input argument pos passed by reference. The method also returns whether the point t is
	 * dominated by any point in the container or not.
	 * Worst case complexity: O(n*d).
	 */
	inline bool updateSkyline(tuple_t const &t, int &pos) {
		// if Tdom[i] >= 1, the i-th point is dominated by t and can be shifted
		int Tdom[s->last_idx+1];
		// if dom[i] >= 1, the i-th point dominates t
		int domT[s->last_idx+1];
		// flag is true if t is dominated by at least one point in the container
		bool isDominated = false;
		// first dimension d=0
		size_t d = 0;
		for(int i=0; i<=s->last_idx; i++) {
			Tdom[i] = (s->v_dim[d][i] > t.d[d]) - (DIM+1)*(s->v_dim[d][i] < t.d[d]);
			domT[i] = (s->v_dim[d][i] < t.d[d]) - (DIM+1)*(s->v_dim[d][i] > t.d[d]);
		}
		// dimensions from d=1 to d=DIM-2
		for(d=1; d <= DIM-2; d++) {
			for(int i=0; i<=s->last_idx; i++) {
				Tdom[i] += (s->v_dim[d][i] > t.d[d]) - (DIM+1)*(s->v_dim[d][i] < t.d[d]);
				domT[i] += (s->v_dim[d][i] < t.d[d]) - (DIM+1)*(s->v_dim[d][i] > t.d[d]);
			}
		}
		// last dimension d=DIM-1
		for(int i=0; i<=s->last_idx; i++) {
			Tdom[i] += (s->v_dim[d][i] > t.d[d]) - (DIM+1)*(s->v_dim[d][i] < t.d[d]);
			domT[i] += domT[i] += (s->v_dim[d][i] < t.d[d]) - (DIM+1)*(s->v_dim[d][i] > t.d[d]);
			if(domT[i] >= 1) isDominated = true;
		}
		// scan all the tuples and shift the dominated ones to the end of the TData container
		int i = 0;
		pos = s->last_idx;
		while(i <= pos) {
			// if the i-th tuple is dominated by t
			if(Tdom[i] >= 1) {
				// copy the tuple in position s->last_idx in the position i
				for(size_t d=0; d<DIM+2; d++) s->v_dim[d][i] = s->v_dim[d][pos];
				Tdom[i] = Tdom[pos];
				domT[i] = domT[pos];
				pos--;
			}
			else i++;
		}
		return isDominated;
	}

	// method to return the actual number of tuples in the pane
	inline size_t getSize() const {
		return s->getSize();
	}

	// method to return the pointer to the TData of the pane
	inline TData *getData() const {
		return s;
	}

	// method to transform the skyline tuples of the pane into a vector of tuples
	vector<tuple_t> *toVector() const {
		vector<tuple_t> *result_set = new vector<tuple_t>(s->last_idx+1);
		// for all the tuples in the container
		for(int i=0; i<=s->last_idx; i++) {
			// for all the dimensions and timestamps
			for(size_t d=0; d<DIM; d++) (*result_set)[i].d[d] = s->v_dim[d][i];
			(*result_set)[i].app_ts = s->v_dim[DIM][i];
			(*result_set)[i].gen_ts = s->v_dim[DIM+1][i];
		}
		return result_set;
	}
};

/*! 
 *  \class SkyWin
 *  
 *  \brief Window Implementation for Skyline Queries
 *  
 *  Window implementation for continuous skyline queries.
 *  
 *  This class is defined in \ref Pane_Farming/include/queries.hpp
 */
class SkyWin: public Window<SkyPane> {
private:
	// pointer to the tuple container
	TData *w;

public:
	// constructor I
	SkyWin(size_t _id=0): Window<SkyPane>(_id) {
		// create the TData container
		w = new TData();
		// initialize the ready time to a large enough value
		this->setReadyTime(getticks());
	}

	// destructor
	~SkyWin() {
		delete w;
	}

	/* 
	 * Method to add a pane to the window.
	 */
	inline void addPane(SkyPane const &p) {
		// update the number of pane instances of the window
		n_pane_instances++;
		// update the ready time of the window
		if(p.getStartingTime() < this->getReadyTime()) this->setReadyTime(p.getStartingTime());
		// insert all the tuples of the pane in window
		TData *l_sky = p.getData(); // local skyline of the input pane
		for(int i=0; i<=l_sky->last_idx; i++) {
			// create the tuple
			tuple_t t;
			for(size_t d=0; d<DIM; d++) t.d[d] = l_sky->v_dim[d][i];
			t.app_ts = l_sky->v_dim[DIM][i];
			t.gen_ts = l_sky->v_dim[DIM+1][i];
			w->append(t);
		}
	}

	/* 
	 * Method to add a pane to the window and incrementally update the internal
	 * window result.
	 */
	inline void addComputePane(SkyPane const &p) {
		// update the number of pane instances of the window
		n_pane_instances++;
		// update the ready time of the window
		if(p.getStartingTime() < this->getReadyTime()) this->setReadyTime(p.getStartingTime());
		// update the window result (skyline) by trying to add all the tuples in the input pane
		TData *l_sky = p.getData(); // local skyline of the input pane
		int endpoint = this->getSize() - 1;
		for(int i=0; i<=l_sky->last_idx; i++) {
			// create the tuple
			tuple_t t;
			for(size_t d=0; d<DIM; d++) t.d[d] = l_sky->v_dim[d][i];
			t.app_ts = l_sky->v_dim[DIM][i];
			t.gen_ts = l_sky->v_dim[DIM+1][i];
			int pos;
			bool isDominated = updateSkyline(t, pos, endpoint);
			// pos is the index of the last point in the container not dominated by t
			if(pos < w->last_idx) {
				// remove all the points in the container dominated by t
				w->removeFromIndex(pos+1);
			}
			// if t does not dominate any point in the container and it is not dominated, we append it to the container
			if((pos == w->last_idx) && !isDominated) w->append(t);
		}
	}

	/* 
	 * Method to shift all the tuples dominated by t at the end of the TData container of the window.
	 * The index of the last tuple not dominated by t (-1 if all the elements are dominated) is written
	 * in the input argument pos passed by reference. The method also returns whether the point t is
	 * dominated by any point in the container or not.
	 * Worst case complexity: O(n*d).
	 */
	inline bool updateSkyline(tuple_t const &t, int &pos, int &endpoint) {
		// if Tdom[i] >= 1, the i-th point is dominated by t and can be shifted
		int Tdom[endpoint+1];
		// if dom[i] >= 1, the i-th point dominates t
		int domT[endpoint+1];
		// flag is true if t is dominated by at least one point in the container
		bool isDominated = false;
		// first dimension d=0
		size_t d = 0;
		for(int i=0; i<=endpoint; i++) {
			Tdom[i] = (w->v_dim[d][i] > t.d[d]) - (DIM+1)*(w->v_dim[d][i] < t.d[d]);
			domT[i] = (w->v_dim[d][i] < t.d[d]) - (DIM+1)*(w->v_dim[d][i] > t.d[d]);
		}
		// dimensions from d=1 to d=DIM-2
		for(d=1; d <= DIM-2; d++) {
			for(int i=0; i<=endpoint; i++) {
				Tdom[i] += (w->v_dim[d][i] > t.d[d]) - (DIM+1)*(w->v_dim[d][i] < t.d[d]);
				domT[i] += (w->v_dim[d][i] < t.d[d]) - (DIM+1)*(w->v_dim[d][i] > t.d[d]);
			}
		}
		// last dimension d=DIM-1
		for(int i=0; i<=endpoint; i++) {
			Tdom[i] += (w->v_dim[d][i] > t.d[d]) - (DIM+1)*(w->v_dim[d][i] < t.d[d]);
			domT[i] += domT[i] += (w->v_dim[d][i] < t.d[d]) - (DIM+1)*(w->v_dim[d][i] > t.d[d]);
			if(domT[i] >= 1) isDominated = true;
		}
		// scan all the tuples and shift the dominated ones to the end of the TData container
		int i = 0;
		pos = w->last_idx;
		while(i <= endpoint) {
			// if the i-th tuple is dominated by t
			if(Tdom[i] >= 1) {
				// copy the tuple in position endpoint in the position i
				for(size_t d=0; d<DIM+2; d++) w->v_dim[d][i] = w->v_dim[d][endpoint];
				Tdom[i] = Tdom[endpoint];
				domT[i] = domT[endpoint];
				// copy the tuple in position pos in position endpoint
				for(size_t d=0; d<DIM+2; d++) w->v_dim[d][endpoint] = w->v_dim[d][pos];
				endpoint--;
				pos--;
			}
			else i++;
		}
		return isDominated;
	}

	// method to return the actual number of tuples in the window
	inline size_t getSize() const {
		return w->getSize();
	}

	// method to transform the skyline tuples of the window into a vector of tuples
	vector<tuple_t> *toVector() const {
		vector<tuple_t> *result_set = new vector<tuple_t>(w->last_idx+1);
		// for all the tuples in the container
		for(int i=0; i<=w->last_idx; i++) {
			// for all the dimensions and timestamps
			for(size_t d=0; d<DIM; d++) (*result_set)[i].d[d] = w->v_dim[d][i];
			(*result_set)[i].app_ts = w->v_dim[DIM][i];
			(*result_set)[i].gen_ts = w->v_dim[DIM+1][i];
		}
		return result_set;
	}
};

/*! 
 *  \class TopdPane
 *  
 *  \brief Pane Implementation for Top-δ dominant skyline query
 *  
 *  Pane implementation for continuous top-δ dominant skyline query.
 *  
 *  This class is defined in \ref Pane_Farming/include/queries.hpp
 */
class TopdPane: public Pane {
private:
	// pointer to the tuple container
	TData *s;
	// vectors of the dominance factors of the tuples in the container
	vector<size_t> domFactors;

public:
	// constructor I
	TopdPane(size_t w_cnt=1, size_t _id=0, size_t _inst=0, size_t _len=0): Pane(w_cnt, _id, _inst, _len) {
		// create the TData container
		s = new TData();
		domFactors.reserve(RESERVE_V_SPACE);
	}

	// destructor
	~TopdPane() {
		delete s;
	}

	/* 
	 * Method to add a tuple to the pane and incrementally update the internal
	 * pane result.
	 */
	inline void addComputeTuple(tuple_t const &t) {
		// set the starting time of the pane
		if(this->getStartingTime() == 0) this->setStartingTime(t.arrival_ticks);
		else if(this->getStartingTime() > t.arrival_ticks) this->setStartingTime(t.arrival_ticks);
		int pos;
		size_t dom = updateTopd(t, pos);
		if(pos < s->last_idx) {
			s->removeFromIndex(pos+1);
			domFactors.resize(pos+1);
		}
		if(dom < DIM) {
			// t is a skyline point, we append it in the container
			s->append(t);
			domFactors.push_back(dom);
		}
		// set the ending time of the pane
		this->setEndingTime(getticks());
	}

	/* 
	 * Method to return the dominance factor of a tuple t given the tuples stored
	 * in the container of the pane. This method also updates the dominance factor
	 * of all the tuples in the container. The points with dominance factor
	 * equal to DIM (non skyline points) are shifted in the last positions.
	 * Worst case complexity: O(n*d).
	 */
	inline size_t updateTopd(tuple_t const &t, int &index) {
		size_t domFactor = 0; // dominance factor of t
		int domT[s->last_idx+1]; // domT[i] = K if the i-th tuple is equal/better than t in exactly K dimensions
		int noBetterT[s->last_idx+1];
		// at the end all the tuples with id in [0, idx] will be the ones not dominated by t
		int idx = s->last_idx;
		int Tdom[s->last_idx+1]; // dom[i] = K if t is equal/better than the i-th tuple in exactly K dimensions
		int TnoBetter[s->last_idx+1];
		// first dimension d=0
		size_t d = 0;
		for(int i=0; i<=s->last_idx; i++) {
			domT[i] = (s->v_dim[d][i] <= t.d[d]);
			noBetterT[i] = !(s->v_dim[d][i] <= t.d[d]) || (s->v_dim[d][i] == t.d[d]);
			Tdom[i] = (t.d[d] <= s->v_dim[d][i]);
			TnoBetter[i] = !(t.d[d] <= s->v_dim[d][i]) || (t.d[d] == s->v_dim[d][i]);
		}
		// dimensions from d=1 to d=DIM-2
		for(d=1; d <= DIM-2; d++) {
			for(int i=0; i<=s->last_idx; i++) {
				domT[i] += (s->v_dim[d][i] <= t.d[d]);
				noBetterT[i] *= !(s->v_dim[d][i] <= t.d[d]) || (s->v_dim[d][i] == t.d[d]);
				Tdom[i] += (t.d[d] <= s->v_dim[d][i]);
				TnoBetter[i] = !(t.d[d] <= s->v_dim[d][i]) || (t.d[d] == s->v_dim[d][i]);
			}
		}
		// last dimension d=DIM-1
		for(int i=0; i<=s->last_idx; i++) {
			domT[i] += (s->v_dim[d][i] <= t.d[d]);
			noBetterT[i] *= !(s->v_dim[d][i] <= t.d[d]) || (s->v_dim[d][i] == t.d[d]);
			domT[i] *= (1 - noBetterT[i]);
			// update the dominance factor of t
			if(domT[i] > domFactor) domFactor = domT[i];
			Tdom[i] += (t.d[d] <= s->v_dim[d][i]);
			TnoBetter[i] = !(t.d[d] <= s->v_dim[d][i]) || (t.d[d] == s->v_dim[d][i]);
			Tdom[i] *= (1 - TnoBetter[i]);
			if(Tdom[i] > domFactors[i]) domFactors[i] = Tdom[i];
		}
		// scan all the tuples and shift the dominated ones to the end of the TData container
		int i = 0;
		while(i <= idx) {
			// if the i-th tuple is dominated by t
			if(domFactors[i] == DIM) {
				// switch the tuple with the one in the last position
				for(size_t d=0; d<DIM+2; d++) s->v_dim[d][i] = s->v_dim[d][idx];
				domFactors[i] = domFactors[idx];
				idx--;
			}
			else i++;
		}
		index = idx;
		return domFactor;
	}

	// method to return the actual number of tuples in the pane
	inline size_t getSize() const {
		return s->getSize();
	}

	// method to return the pointer to the TData of the pane
	inline TData *getData() const {
		return s;
	}

	// method to return the vector of dominance factors of the tuples in the TData container
	inline vector<size_t> const &getDomVector() const {
		return domFactors;
	}

	// method to transform the top-δ tuples of the pane into a vector of tuples
	vector<tuple_t> *toVector() const {
		// determine the number of tuples per dominance factor
		size_t n_tuples[DIM+1]; // n_tuples[i] = no. of tuples with dominance factor i
		fill_n(n_tuples, DIM+1, 0);
		for(size_t i=0; i<=s->last_idx; i++) {
			n_tuples[domFactors[i]]++;
		}
		/* Determine the smallest dominance factor kdom for which
	   	   there exist at least δ=TOPD_NUMBER tuples with dominance factor
	   	   lower or equal than kdom. */
		size_t kdom = 0;
		size_t needed_tuples = TOPD_NUMBER; // we use δ=TOPD_NUMBER
		size_t counters[DIM];
		fill_n(counters, DIM, 0);
		while((kdom < DIM) && (needed_tuples > 0)) {
			if(needed_tuples <= n_tuples[kdom]) {
				counters[kdom] = needed_tuples;
				needed_tuples = 0;
			}
			else {
				needed_tuples -= n_tuples[kdom];
				counters[kdom] = n_tuples[kdom];
				kdom++;
			}
		}
		// create the final vector
		needed_tuples = 0;
		vector<tuple_t> *result = new vector<tuple_t>();
		for(int i=0; i<=s->last_idx; i++) {
			if(counters[domFactors[i]] > 0) {
				counters[domFactors[i]]--;
				needed_tuples++;
				// create the tuple
				tuple_t t;
				for(size_t d=0; d<DIM; d++) t.d[d] = s->v_dim[d][i];
				t.app_ts = s->v_dim[DIM][i];
				t.gen_ts = s->v_dim[DIM+1][i];
				result->push_back(t);
				if(needed_tuples >= TOPD_NUMBER) break;
			}
		}
		return result;
	}
};

/*! 
 *  \class TopdWin
 *  
 *  \brief Window Implementation for Top-δ Dominant Skyline Query
 *  
 *  Window implementation for continuous Top-δ dominant skyline query.
 *  
 *  This class is defined in \ref Pane_Farming/include/queries.hpp
 */
class TopdWin: public Window<TopdPane> {
private:
	// pointer to the tuple container
	TData *w;
	// vectors of the dominance factors of the tuples in the container
	vector<size_t> domFactors;

public:
	// constructor I
	TopdWin(size_t _id=0): Window<TopdPane>(_id) {
		// create the TData container
		w = new TData();
		// initialize the ready time to a large enough value
		this->setReadyTime(getticks());
	}

	// destructor
	~TopdWin() {
		delete w;
	}

	/* 
	 * Method to add a pane to the window.
	 */
	inline void addPane(TopdPane const &p) {
		// update the number of pane instances of the window
		n_pane_instances++;
		// update the ready time of the window
		if(p.getStartingTime() < this->getReadyTime()) this->setReadyTime(p.getStartingTime());
		// insert all the tuples of the pane in window
		TData *l_sky = p.getData(); // local skyline of the input pane
		vector<size_t> const &doms = p.getDomVector(); // dominance vector of the input pane
		for(int i=0; i<=l_sky->last_idx; i++) {
			// create the tuple
			tuple_t t;
			for(size_t d=0; d<DIM; d++) t.d[d] = l_sky->v_dim[d][i];
			t.app_ts = l_sky->v_dim[DIM][i];
			t.gen_ts = l_sky->v_dim[DIM+1][i];
			w->append(t);
			domFactors.push_back(doms[i]);
		}
	}

	/* 
	 * Method to add a pane to the window and incrementally update the internal
	 * window result.
	 */
	inline void addComputePane(TopdPane const &p) {
		// update the number of pane instances of the window
		n_pane_instances++;
		// update the ready time of the window
		if(p.getStartingTime() < this->getReadyTime()) this->setReadyTime(p.getStartingTime());
		// update the window result (skyline) by trying to add all the tuples in the input pane
		TData *l_sky = p.getData(); // local skyline of the input pane
		vector<size_t> const &doms = p.getDomVector(); // dominance vector of the input pane
		for(int i=0; i<=l_sky->last_idx; i++) {
			// create the tuple
			tuple_t t;
			for(size_t d=0; d<DIM; d++) t.d[d] = l_sky->v_dim[d][i];
			t.app_ts = l_sky->v_dim[DIM][i];
			t.gen_ts = l_sky->v_dim[DIM+1][i];
			int pos;
			size_t dom = updateTopd(t, pos);
			if(pos < w->last_idx) {
				w->removeFromIndex(pos+1);
				domFactors.resize(pos+1);
			}
			if(dom < DIM) {
				// t is a skyline point, we append it in the container
				w->append(t);
				domFactors.push_back(max(dom, doms[i]));
			}
		}
	}

	/* 
	 * Method to return the dominance factor of a tuple t given the tuples stored
	 * in the container of the window. This method also updates the dominance factor
	 * of all the tuples in the container. The points with dominance factor
	 * equal to DIM (non skyline points) are shifted in the last positions.
	 * Worst case complexity: O(n*d).
	 */
	inline size_t updateTopd(tuple_t const &t, int &index) {
		size_t domFactor = 0; // dominance factor of t
		int domT[w->last_idx+1]; // domT[i] = K if the i-th tuple is equal/better than t in exactly K dimensions
		int noBetterT[w->last_idx+1];
		// at the end all the tuples with id in [0, idx] will be the ones not dominated by t
		int idx = w->last_idx;
		int Tdom[w->last_idx+1]; // dom[i] = K if t is equal/better than the i-th tuple in exactly K dimensions
		int TnoBetter[w->last_idx+1];
		// first dimension d=0
		size_t d = 0;
		for(int i=0; i<=w->last_idx; i++) {
			domT[i] = (w->v_dim[d][i] <= t.d[d]);
			noBetterT[i] = !(w->v_dim[d][i] <= t.d[d]) || (w->v_dim[d][i] == t.d[d]);
			Tdom[i] = (t.d[d] <= w->v_dim[d][i]);
			TnoBetter[i] = !(t.d[d] <= w->v_dim[d][i]) || (t.d[d] == w->v_dim[d][i]);
		}
		// dimensions from d=1 to d=DIM-2
		for(d=1; d <= DIM-2; d++) {
			for(int i=0; i<=w->last_idx; i++) {
				domT[i] += (w->v_dim[d][i] <= t.d[d]);
				noBetterT[i] *= !(w->v_dim[d][i] <= t.d[d]) || (w->v_dim[d][i] == t.d[d]);
				Tdom[i] += (t.d[d] <= w->v_dim[d][i]);
				TnoBetter[i] = !(t.d[d] <= w->v_dim[d][i]) || (t.d[d] == w->v_dim[d][i]);
			}
		}
		// last dimension d=DIM-1
		for(int i=0; i<=w->last_idx; i++) {
			domT[i] += (w->v_dim[d][i] <= t.d[d]);
			noBetterT[i] *= !(w->v_dim[d][i] <= t.d[d]) || (w->v_dim[d][i] == t.d[d]);
			domT[i] *= (1 - noBetterT[i]);
			// update the dominance factor of t
			if(domT[i] > domFactor) domFactor = domT[i];
			Tdom[i] += (t.d[d] <= w->v_dim[d][i]);
			TnoBetter[i] = !(t.d[d] <= w->v_dim[d][i]) || (t.d[d] == w->v_dim[d][i]);
			Tdom[i] *= (1 - TnoBetter[i]);
			if(Tdom[i] > domFactors[i]) domFactors[i] = Tdom[i];
		}
		// scan all the tuples and shift the dominated ones to the end of the TData container
		int i = 0;
		while(i <= idx) {
			// if the i-th tuple is dominated by t
			if(domFactors[i] == DIM) {
				// switch the tuple with the one in the last position
				for(size_t d=0; d<DIM+2; d++) w->v_dim[d][i] = w->v_dim[d][idx];
				domFactors[i] = domFactors[idx];
				idx--;
			}
			else i++;
		}
		index = idx;
		return domFactor;
	}

	// method to return the actual number of tuples in the window
	inline size_t getSize() const {
		return w->getSize();
	}

	// method to transform the top-δ tuples of the window into a vector of tuples
	vector<tuple_t> *toVector() const {
		// determine the number of tuples per dominance factor
		size_t n_tuples[DIM+1]; // n_tuples[i] = no. of tuples with dominance factor i
		fill_n(n_tuples, DIM+1, 0);
		for(size_t i=0; i<=w->last_idx; i++) {
			n_tuples[domFactors[i]]++;
		}
		/* Determine the smallest dominance factor kdom for which
	   	   there exist at least δ=TOPD_NUMBER tuples with dominance factor
	   	   lower or equal than kdom. */
		size_t kdom = 0;
		size_t needed_tuples = TOPD_NUMBER; // we use δ=TOPD_NUMBER
		size_t counters[DIM];
		fill_n(counters, DIM, 0);
		while((kdom < DIM) && (needed_tuples > 0)) {
			if(needed_tuples <= n_tuples[kdom]) {
				counters[kdom] = needed_tuples;
				needed_tuples = 0;
			}
			else {
				needed_tuples -= n_tuples[kdom];
				counters[kdom] = n_tuples[kdom];
				kdom++;
			}
		}
		// create the final vector
		needed_tuples = 0;
		vector<tuple_t> *result = new vector<tuple_t>();
		for(int i=0; i<=w->last_idx; i++) {
			if(counters[domFactors[i]] > 0) {
				counters[domFactors[i]]--;
				needed_tuples++;
				// create the tuple
				tuple_t t;
				for(size_t d=0; d<DIM; d++) t.d[d] = w->v_dim[d][i];
				t.app_ts = w->v_dim[DIM][i];
				t.gen_ts = w->v_dim[DIM+1][i];
				result->push_back(t);
				if(needed_tuples >= TOPD_NUMBER) break;
			}
		}
		return result;
	}
};

#endif
