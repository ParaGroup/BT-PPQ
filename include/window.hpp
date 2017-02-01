/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file window.hpp
 *  \brief Window
 *  
 *  This file provides a generic window class. A window is a data structure
 *  in which we accomulate all the panes belonging to the window and we
 *  update the corresponding window result. Each window is composed by a fixed
 *  number of panes depending on the window specifications (window length and
 *  slide parameters).
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

#ifndef WIN_H
#define WIN_H

// include
#include <general.hpp>

using namespace std;

/* This file provides the following class:
 *  -Window: class of a generic window. It must be extended to implement a specific query.
 */

/*! 
 *  \class Window
 *  
 *  \brief Window
 *  
 *  This class must be extended to provide the specific implementation of the
 *  window for a desired continuous query.
 *  
 *  This class is defined in \ref Pane_Farming/include/window.hpp
 */
template <typename Pane_t>
class Window {
protected:
	size_t id; // identifier of the window (starting from 0)
	/* Ready time (in ticks) of the window: it is the smallest
	   starting time of a pane in the window. */
	volatile ticks ready_time = 0;
	size_t n_pane_instances = 0; // number of pane instances of the window

	// constructor
	Window(size_t _id): id(_id) {}

public:
	// method to add a pane to the window and update the window content
	virtual void addComputePane(const Pane_t &pane) = 0;
	// method to return the size of the window
	virtual size_t getSize() const = 0;
	// method to return the identifier of the window (starting from 0)
	inline size_t getId() const { return id; }
	// method to return the number of pane instances of the window
	inline size_t getNoInstances() const { return n_pane_instances; }
	// method to return the ready time of the window (in ticks)
	inline ticks getReadyTime() const { return ready_time; }
	// method to set the identifier of the window (starting from 0)
	inline void setId(size_t _id) { id = _id; }
	// method to set the ready time of the window (in ticks)
	inline void setReadyTime(ticks _time) { ready_time = _time; }
};

#endif
