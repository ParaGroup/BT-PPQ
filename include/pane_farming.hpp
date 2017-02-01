/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file pane_farming.hpp
 *  \brief Pane Farming
 *  
 *  Implementation of the Pane Farming parallel pattern. It is implemented as
 *  a pipeline of two stages, each one potentially parallel (farm). The first
 *  stage is called Pane-level Query (PLQ) and processes each pane independently
 *  (a pane is a tumbling window of tuples). Once processed, panes are transmitted
 *  to a second stage called Window-level Query (WLQ) that merges all the panes
 *  belonging to the same window into a final result. Panes in common to more than
 *  one window can be re-used and not computed multiple times.
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

#ifndef PF_H
#define PF_H

// include
#include <functional>
#include <ff/farm.hpp>
#include <ff/pipeline.hpp>
#include <plq.hpp>
#include <wlq.hpp>
#include <queries.hpp>

using namespace ff;
using namespace std;

/* This file provides the following class:
 *  -PF: class that implements the Pane Farming.
 */

/*! 
 *  \class PF
 *  
 *  \brief Pane Farming
 *  
 *  Pane Farming is a parallel pattern for window-based streaming computations.
 *  It consists in a pipeline of two stages (PLQ and WLQ). The first processes
 *  panes. The second works on count-based sliding windows of computed panes
 *  and merges them by producing a final result for each window.
 *  
 *  This class is defined in \ref Pane_Farming/include/pane_farming.hpp
 */
#if defined(PLQ_ONLY) // <----- version with the PLQ stage only!
template<typename Pane_t>
class PF: public ff_pipeline {
private:
#if !defined(IOP)
    typedef PLQ_Emitter plq_emitter_t;
    typedef PLQ_Worker<Pane_t> plq_worker_t;
#else
    typedef PLQ_Emitter_IOP plq_emitter_t;
    typedef PLQ_Worker_IOP<Pane_t> plq_worker_t;
#endif
	size_t nw_plq; // pardegree PLQ stage
	size_t nw_wlq; // pardgree WLQ stage
    plq_emitter_t *plq_emitter; // pointer to the emitter of the first stage (PLQ)
    vector<ff_node *> plq_w_vect; // vector of pointers to workers of the first stage (PLQ)
    ff_farm<> *plq_farm; // pointer to the first stage (PLQ)

public:
    /**
     * \brief constructor I
     * 
     * \param _wlen window length (in ms)
     * \param _wslide window slide (in ms)
     * \param _nw_plq pardegree of the PLQ stage
     * \param _nw_wlq pardegree of the WLQ stage
     * \param _drop_prob desired dropping probability (0 uses the standard slack mechanism)
     * \param _port port to receive connection from the Generator
     * 
     */
    PF(size_t _wlen, size_t _wslide, size_t _nw_plq, size_t _nw_wlq, double _drop_prob, size_t _port): ff_pipeline(), nw_plq(_nw_plq), nw_wlq(_nw_wlq) {
        // calculate the length in ms of the pane
        auto gcd = [](size_t u, size_t v) {
            while(v != 0) {
                unsigned long r = u % v;
                u = v;
                v = r;
            }
            return u;
        };
        size_t pane_len = gcd(_wlen, _wslide);
        // calculate the window parameters in terms of panes
        size_t wp = _wlen / pane_len; // window length (in no. panes)
        size_t sp = _wslide / pane_len; // window slide (in no. panes)
        cout << "Window parameters [wp: " << wp << " sp: " << sp << " pane length: " << pane_len << " ms], dimensions per tuple: " << DIM << endl;
        // create the first stage (PLQ): it is a farm without collector
        for(size_t i=0; i<nw_plq; i++) plq_w_vect.push_back(new plq_worker_t(wp, sp, pane_len));
        plq_farm = new ff_farm<>();
        plq_farm->add_workers(plq_w_vect);
#if !defined(IOP)
        plq_emitter = new plq_emitter_t(pane_len, nw_plq, _port, _drop_prob, plq_farm->getlb());
#else
        plq_emitter = new plq_emitter_t(pane_len, nw_plq, _port, _drop_prob, plq_farm->getlb());
#endif
        plq_farm->add_emitter(plq_emitter);
        plq_farm->remove_collector();
        this->add_stage(plq_farm);
        this->setFixedSize(false); // the pipeline uses unbounded queues!
    }

    // destructor
    ~PF() {
        delete plq_emitter;
        for(auto const &node: plq_w_vect) {
            plq_worker_t *w = static_cast<plq_worker_t *>(node);
            delete w;
        }
        delete plq_farm;
    }
};
#else // <----- complete version of the pattern
template<typename Pane_t, typename Win_t>
class PF: public ff_pipeline {
private:
#if !defined(IOP)
    typedef PLQ_Emitter plq_emitter_t;
    typedef PLQ_Worker<Pane_t> plq_worker_t;
#else
    typedef PLQ_Emitter_IOP plq_emitter_t;
    typedef PLQ_Worker_IOP<Pane_t> plq_worker_t;
#endif
#if !defined(WLQ_OD)
    typedef WLQ_Emitter<Pane_t, Win_t> wlq_emitter_t;
#else
    typedef WLQ_Emitter_OD<Pane_t, Win_t> wlq_emitter_t;
#endif
    typedef WLQ_Worker<Pane_t, Win_t> wlq_worker_t;
    typedef WLQ_Collector<Win_t> wlq_collector_t;
    size_t nw_plq; // pardegree PLQ stage
    size_t nw_wlq; // pardgree WLQ stage
    plq_emitter_t *plq_emitter; // pointer to the emitter of the first stage (PLQ)
    vector<ff_node *> plq_w_vect; // vector of pointers to workers of the first stage (PLQ)
    wlq_emitter_t *wlq_emitter; // pointer to the emitter of the second stage (WLQ)
    vector<ff_node *> wlq_w_vect; // vector of pointers to workers of the second stage (WLQ)
    wlq_collector_t *wlq_collector; // pointer to the collector of the second stage (WLQ)
    ff_farm<> *plq_farm; // pointer to the first stage (PLQ)
    ff_farm<> *wlq_farm; // pointer to the second stage (WLQ)

public:
    /**
     * \brief constructor I
     * 
     * \param _wlen window length (in ms)
     * \param _wslide window slide (in ms)
     * \param _nw_plq pardegree of the PLQ stage
     * \param _nw_wlq pardegree of the WLQ stage
     * \param _drop_prob desired dropping probability (0 uses the base slack mechanism)
     * \param _port port to receive connection from the Generator
     * 
     */
    PF(size_t _wlen, size_t _wslide, size_t _nw_plq, size_t _nw_wlq, double _drop_prob, size_t _port): ff_pipeline(), nw_plq(_nw_plq), nw_wlq(_nw_wlq) {
        assert(nw_plq > 0 && nw_wlq > 0);
        assert(_wlen >= _wslide);
        // calculate the length in ms of the pane
        auto gcd = [](size_t u, size_t v) {
            while(v != 0) {
                unsigned long r = u % v;
                u = v;
                v = r;
            }
            return u;
        };
        size_t pane_len = gcd(_wlen, _wslide);
        // calculate the window parameters in terms of panes
        size_t wp = _wlen / pane_len; // window length (in no. panes)
        size_t sp = _wslide / pane_len; // window slide (in no. panes)
        cout << "Window parameters [wp: " << wp << " sp: " << sp << " pane length: " << pane_len << " ms], dimensions per tuple: " << DIM << endl;
        // create the first stage (PLQ): it is a farm without collector
        for(size_t i=0; i<nw_plq; i++) plq_w_vect.push_back(new plq_worker_t(wp, sp, pane_len));
        plq_farm = new ff_farm<>();
        plq_farm->add_workers(plq_w_vect);
#if !defined(IOP)
        plq_emitter = new plq_emitter_t(pane_len, nw_plq, _port, _drop_prob, plq_farm->getlb());
#else
        plq_emitter = new plq_emitter_t(pane_len, nw_plq, _port, _drop_prob, plq_farm->getlb());
#endif
        plq_farm->add_emitter(plq_emitter);
        plq_farm->remove_collector();
        this->add_stage(plq_farm);
        // create the second stage (WLQ): it is a farm with collector
        for(size_t i=0; i<nw_wlq; i++) wlq_w_vect.push_back(new wlq_worker_t());
        wlq_farm = new ff_farm<>();
        wlq_farm->add_workers(wlq_w_vect);
#if !defined(WLQ_OD)
        wlq_emitter = new wlq_emitter_t(wp, sp, nw_wlq, wlq_farm->getlb());
#else
        wlq_emitter = new wlq_emitter_t(wp, sp, nw_wlq, nw_plq, wlq_farm->getlb());
#endif
        wlq_farm->add_emitter(wlq_emitter);
        wlq_farm->remove_collector();
        wlq_farm->setMultiInput(); // the second farm receives directly from the PLQ workers
#if defined(WLQ_OD)
        wlq_farm->wrap_around(); // feedback channels between the WLQ workers and the WLQ emitter
#endif
        wlq_collector = new wlq_collector_t();
        wlq_farm->add_collector(wlq_collector);
        this->add_stage(wlq_farm);
        this->setFixedSize(false); // the pipeline uses unbounded queues!
    }

    // destructor
    ~PF() {
        delete plq_emitter;
        for(auto const &node: plq_w_vect) {
            plq_worker_t *w = static_cast<plq_worker_t *>(node);
            delete w;
        }
        delete plq_farm;
        delete wlq_emitter;
        delete wlq_collector;
        for(auto const &node: wlq_w_vect) {
            wlq_worker_t *w = static_cast<wlq_worker_t *>(node);
            delete w;
        }
        delete wlq_farm;
    }
};
#endif

#endif
