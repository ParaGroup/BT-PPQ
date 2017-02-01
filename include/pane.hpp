/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file pane.hpp
 *  \brief Pane
 *  
 *  A pane is a tumbling window with a specified ts_low and ts_high, which
 *  represent the lower and upper bounds of the application timestamps of
 *  the pane. A pane with ts_low and ts_high will contain all the tuples
 *  with application timestamps ts such that ts_low <= ts < ts_high.
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

#ifndef PANE_H
#define PANE_H

// include
#include <atomic>
#include <general.hpp>
#include <tuple_t.hpp>

using namespace std;

/* This file provides the following class:
 *  -Pane: class of a generic pane. It must be extended to implement a specific query.
 */

/*! 
 *  \class Pane
 *  
 *  \brief Pane
 *  
 *  This class must be extended to provide the specific implementation of the
 *  pane for a desired continuous query executed by the Pane Farming pattern.
 *  
 *  This class is defined in \ref Pane_Farming/include/pane.hpp
 */
class Pane {
protected:
	/* Atomic counter of the number of windows that
	   include this pane. */
	std::atomic_size_t refCounter;
	/* Boundaries (application timestamps). The pane contains
	 * all the tuples with timestamp ts in [ts_low, ts_high). */
	double ts_low;
	double ts_high;
	// identifier of the pane (starting from 0)
	size_t id;
	// counter of the existing instances of the pane
	size_t no_instances;
	/* Starting time (in ticks) of the pane: it is the time at which the
	   PLQ emitter receives the first tuple of the pane. */
	volatile ticks starting_time = 0;
	/* Ending time (in ticks) of the pane: it is the time at which a
	   PLQ worker completes the computation on the last received tuple
	   of the pane. */
	volatile ticks ending_time = 0;
	/* Closing time (in ticks) of the pane: it is the time at which the
	   pane is complete and its result can be produced to the next stage */
	volatile ticks closing_time = 0;

	// constructor
	Pane(size_t _c=0, size_t _id=0, size_t _inst=1, size_t _len=0): refCounter(_c), id(_id), no_instances(_inst) {
		// set the boundaries of the pane
		ts_low = (double) (id * _len);
		ts_high = (double) ((id+1) * _len);
	}

public:
	// method to add a tuple to the pane and update the pane result
	virtual void addComputeTuple(tuple_t const &t) = 0;
	// method to return the actual pane size
	virtual size_t getSize() const = 0;
	// method to return the identifier of the pane (starting from 0)
	inline size_t getId() const { return id; }
	// method to return the number of existing instances of the pane
	inline size_t getNoInstances() const { return no_instances; }
	// method to return the starting time of the pane (in ticks)
	inline ticks getStartingTime() const { return starting_time; }
	// method to return the ending time of the pane (in ticks)
	inline ticks getEndingTime() const { return ending_time; }
	// method to return the closing time of the pane (in ticks)
	inline ticks getClosingTime() const { return closing_time; }
	// method to decrease the reference counter of the pane and atomically return the past value
	inline size_t decreaseRefCounter() { return refCounter.fetch_sub(1); }
	// method to set the identifier of the pane (starting from 0)
	inline void setId(size_t _id) { id = _id; }
	// method to set the number of existing instances of the pane
	inline void setNoInstances(size_t ins) { no_instances = ins; }
	// method to set the starting time (in ticks) of the pane
	inline void setStartingTime(ticks _starting_time) { starting_time = _starting_time; }
	// method to set the ending time (in ticks) of the pane
	inline void setEndingTime(ticks _ending_time) { ending_time = _ending_time; }
	// method to set the closing time (in ticks) of the pane
	inline void setClosingTime(ticks _closing_time) { closing_time = _closing_time; }
};

#endif
