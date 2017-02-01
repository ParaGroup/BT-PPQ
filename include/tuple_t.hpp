/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file tuple_t.hpp
 *  \brief Tuple Message
 *  
 *  This file provides the implementation of the general message exchanged
 *  between the Emitter and the Workers of the PLQ stage of Pane Farming.
 *  A message of this type can represent:
 *  1- a EOS (End-Of-Stream) marker that notifies the end of the input stream;
 *  2- a punctuation message;
 *  3- a new_pane_instance message "special message for pane-splitting strategies";
 *  4- a regular tuple with DIM double-precision numerical attributes and two
 *     timestamps (the application and generator timestamps).
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

#ifndef _TUPLE_H
#define _TUPLE_H

// include
#include <string.h>
#include <atomic>
#include <sstream>
#include <iostream>

using namespace std;

/* This file provides the following class:
 *  -tuple_t: type of the message exchanged between Emitter/Workers of the PLQ stage.
 */

/*! 
 *  \class tuple_t
 *  
 *  \brief Implementation of the Tuple Message
 *  
 *  This file provides the implementation of the general message exchanged
 *  between the Emitter and the Workers of the PLQ stage of Pane Farming.
 *  A message of this type can represent:
 *  1- a EOS (End-Of-Stream) marker that notifies the end of the input stream;
 *  2- a punctuation message;
 *  3- a new_pane_instance message "special message for pane-splitting strategies";
 *  4- a regular tuple with DIM double-precision numerical attributes and two
 *     timestamps (the application and generator timestamps).
 *  
 *  This class is defined in \ref Pane_Farming/include/tuple_t.hpp
 */
class tuple_t {
public:
    // if this flag is true the tuple is a EOS marker
    bool isEOS;
    // otherwise, if this flag is true the tuple is a punctuation
    bool isPunctuation;
    double punc_value; // punctuation value
    // otherwise, if this flag is true the tuple is a new_pane_instance message
    bool isNewInstance;
    size_t no_instances; // number of instances of the pane
    atomic_size_t refCounter; // number of workers that will receive the message
    // otherwise the following fields are meaningful (regular tuple)
    size_t pane_id; // id of the pane of the tuple (starting from 0)
	double d[DIM]; // array of DIM attributes
	double app_ts; // application timestamp (in usecs)
	double gen_ts; // generator timestamp (in usecs)
    bool isLull; // if true the tuple will be used to build an empty pane
    volatile ticks arrival_ticks; // time in ticks of the arrival of this tuple to the PLQ Emitter

    // empty constructor
    tuple_t() {
        isEOS = false;
        isPunctuation = false;
        punc_value = 0;
        isNewInstance = false;
        no_instances = 0;
        refCounter = 0;
        pane_id = 0;
        fill_n(d, DIM, 0);
        app_ts = 0;
        gen_ts = 0;
        isLull = false;
        arrival_ticks = 0;
    }

    // copy constructor
    tuple_t(tuple_t const &t) {
        isEOS = t.isEOS;
        isPunctuation = t.isPunctuation;
        punc_value = t.punc_value;
        isNewInstance = t.isNewInstance;
        no_instances = t.no_instances;
        refCounter.store((t.refCounter).load());
        pane_id = t.pane_id;
        memcpy(d, t.d, sizeof(double)*DIM);
        app_ts = t.app_ts;
        gen_ts = t.gen_ts;
        isLull = t.isLull;
        arrival_ticks = t.arrival_ticks;
    }

    // copy operator
    tuple_t &operator= (const tuple_t& other) {
        isEOS = other.isEOS;
        isPunctuation = other.isPunctuation;
        punc_value = other.punc_value;
        isNewInstance = other.isNewInstance;
        no_instances = other.no_instances;
        refCounter.store((other.refCounter).load());
        pane_id = other.pane_id;
        memcpy(d, other.d, sizeof(double)*DIM);
        app_ts = other.app_ts;
        gen_ts = other.gen_ts;
        isLull = other.isLull;
        arrival_ticks = other.arrival_ticks;
        return *this;
    }

    // get the size of the tuple's fields in bytes (only for regular tuples)
    static inline size_t getSize() {
        // only the fields d, app_ts, gen_ts
        return sizeof(double)*(DIM+2);
    }

    // serialization method (only for regular tuples)
    inline int serialize(char *buffer) {
        if(isEOS || isPunctuation || isNewInstance) return -1;
    	char *data = (char *) d;
    	int b = 0;
    	for(size_t i = 0; i<DIM; i++) {
    		for(size_t j = 0; j<sizeof(double); j++) buffer[b++] = data[(i*sizeof(double))+j];
    	}
    	data = (char *) &app_ts;
    	for(size_t j = 0; j<sizeof(double); j++) buffer[b++] = data[j];
    	data = (char *) &gen_ts;
    	for(size_t j = 0; j<sizeof(double); j++) buffer[b++] = data[j];
        return 0;
    }

    // deserialization method (only for regular tuples)
    inline int deserialize(char *buffer) {
        if(isEOS || isPunctuation || isNewInstance) return -1;
    	char *data = (char *) d;
    	int b = 0;
    	for(size_t i=0; i<DIM; i++) {
    		for(size_t j=0; j<sizeof(double); j++) data[(i*sizeof(double))+j] = buffer[b++];
    	}
    	data = (char *) &app_ts;
    	for(size_t j=0; j<sizeof(double); j++) data[j] = buffer[b++];
    	data = (char *) &gen_ts;
    	for(size_t j=0; j<sizeof(double); j++) data[j] = buffer[b++];
        return 0;
    }

    // method to set the tuple as a punctuation
    inline void set_P(double ts, size_t p_id, size_t nw) {
        isPunctuation = true;
        punc_value = ts; // punctuation value
        pane_id = p_id; // pane id corresponding to the punctuation value
        refCounter = nw; // number of receivers of the punctuation
    }

    // method to set the tuple as a EOS marker
    inline void set_EOS(size_t nw) {
        isEOS = true;
        refCounter = nw; // number of receivers of the EOS marker
    }

    // method to set the tuple as a new_pane_instance message
    inline void set_NewInstance(size_t p_id, size_t no_inst, size_t nw) {
        isNewInstance = true;
        pane_id = p_id; // pane id
        no_instances = no_inst; // new number of instances of pane id
        refCounter = nw; // number of receivers of the message
    }

    // print method (only for regular tuples)
    void print() const {
        if(isEOS || isPunctuation || isNewInstance) return;
        std::ostringstream print;
        print << "Tuple: [";
        for(size_t i=0; i<DIM-1; i++) {
            print << d[i] << ", ";
        }
        print << d[DIM-1] << "] app_ts " << app_ts << " gen_ts " << gen_ts;
        cout << print.str() << endl;
    }
};

#endif
