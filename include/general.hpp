/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file general.hpp
 *  \brief General Definitions and Macros
 *  
 *  This file contains a set of useful definitions and macros used by other
 *  files of the project.
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

#ifndef _GENERAL_H
#define _GENERAL_H

// include
#include <mutex>
#include <fstream>
#include <iostream>
#include <ff/cycle.h>

// define
#if (defined(PIANOSA) || defined(PIANOSAU))
    #define FREQ 2000
#endif

#if defined(REPARA)
    #define FREQ 2400
#endif

#if defined(REPHRASE)
    #define FREQ 512 // magic value, trust me!
#endif

#define FROM_TICKS_TO_USECS(ticks) (((double) ticks)/FREQ)
#define FROM_USECS_TO_TICKS(usecs) (ticks) (FREQ * usecs)
#define SAMPLING_PUNCTUATION_USEC 100000 // sampling interval in usecs (used by the Adaptive K-slack Punctuation Mechanism)
#define SAMPLING_SCHED_PLQ_USEC 250000 // sampling interval in usecs used by the PLQ emitter
#define SAMPLING_SCHED_WLQ_USEC 2000000 // sampling interval in usecs used by the WLQ emitter
#define MAX_TIME_BETWEEN_PUNC_USEC 10000 // maximum time in usec between two consecutive punctuations (used to limit the no. of punctuations)
#define RESERVE_V_SPACE 500 // pre-allocated size of vectors and deques throughout the project
#define PRINT_RATE_USEC 10000000 // frequency of prints adopted by the Client in usecs
#define TICKS2WAIT 10000 // number of ticks of the busy-waiting phase
#define TOPD_NUMBER 100 // number of δ>0 tuples to produce in the top-δ dominant skyline query
#define PLQ_WORKERS_LOG_FILENAME "log/workers_plq.log" // filename of the PLQ workers log file
#define WLQ_WORKERS_LOG_FILENAME "log/workers_wlq.log" // filename of the WLQ workers log file
#define WF_WORKERS_LOG_FILENAME "log/workers_wf.log" // filename of the WF workers log file
#define STORED_TUPLES_LOG_FILENAME "log/plq_stored_tuples.log" // filename of the PLQ stored tuples log file

// define to atomically write a string in the stdout
std::mutex mutex_screen; // mutex to avoid mixing couts

#define MY_PRINT(...) { \
 	std::ostringstream stream; \
 	stream << __VA_ARGS__; \
 	mutex_screen.lock(); \
	std::cout << stream.str(); \
	mutex_screen.unlock(); \
}

// set of defines in the case of LOG mode enabled
#if defined(LOG)
	std::mutex mutex_log_plq; // mutex to access the log file of the PLQ workers
	std::mutex mutex_log_wlq; // mutex to access the log file of the WLQ workers
	std::mutex mutex_log_wf; // mutex to access the log file of the WF workers
	std::ofstream logfile_plq; // log file of the PLQ workers
	std::ofstream logfile_wlq; // log file of the WLQ workers
	std::ofstream logfile_wf; // log file of the WF workers

	// main thread must call this macro (PLQ)
	#define LOG_PLQ_INIT() { \
		logfile_plq.open(PLQ_WORKERS_LOG_FILENAME); \
	}

	// main thread must call this macro (PLQ)
	#define LOG_PLQ_CLOSE() { \
		logfile_plq.close(); \
	}

	// macro to write atomically in the log file (PLQ)
	#define LOG_PLQ_WRITE(stream) { \
 		mutex_log_plq.lock(); \
		logfile_plq << stream.str(); \
		mutex_log_plq.unlock(); \
	}

	// main thread must call this macro (WLQ)
	#define LOG_WLQ_INIT() { \
		logfile_wlq.open(WLQ_WORKERS_LOG_FILENAME); \
	}

	// main thread must call this macro (WLQ)
	#define LOG_WLQ_CLOSE() { \
		logfile_wlq.close(); \
	}

	// macro to write atomically in the log file (WLQ)
	#define LOG_WLQ_WRITE(stream) { \
 		mutex_log_wlq.lock(); \
		logfile_wlq << stream.str(); \
		mutex_log_wlq.unlock(); \
	}

	// main thread must call this macro (WF)
	#define LOG_WF_INIT() { \
		logfile_wf.open(WF_WORKERS_LOG_FILENAME); \
	}

	// main thread must call this macro (WF)
	#define LOG_WF_CLOSE() { \
		logfile_wf.close(); \
	}

	// macro to write atomically in the log file (WF)
	#define LOG_WF_WRITE(stream) { \
 		mutex_log_wf.lock(); \
		logfile_wf << stream.str(); \
		mutex_log_wf.unlock(); \
	}
#endif

// sampling barrier to start the sampling synchronously
pthread_barrier_t sampleBarr;

// printing barrier to avoid interleaved final prints
pthread_barrier_t printBarr;

// precedence barriers to print the statistics in a given order
pthread_barrier_t precBarr1;
pthread_barrier_t precBarr2;

// global variable with the setpoint of the utilization factor (used if ADAPTIVE_SCHED_PID is enabled)
double desired_rho;

#endif
