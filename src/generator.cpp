/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file generator.cpp
 *  \brief Input Stream Generator
 *  
 *  Generator is a process reading tuples from a dataset with predefined
 *  timestamps, properly disordered and possibly with bursty characteristics.
 *  The process transmits the tuples to the Pane Farming pattern following
 *  the speed dictated by the generator timestamps of the tuples.
 *  
 *  \note This process communicates with the Pane Farming pattern through a
 *  TCP/IP socket.
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

// include
#include <stdio.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <random>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <ff/cycle.h>
#include <general.hpp>
#include <tuple_t.hpp>
#include <socket_utils.hpp>

using namespace std;

/*! 
 *  \class DatasetReader
 *  
 *  \brief Class to Parse the Dataset Text File
 *  
 *  This class is used by the Generator to read the dataset of tuples from
 *  a text file. Each tuple is a row of the file with the following format:
 *  
 *  d1 d2 ... dn app_timestamp gen_timestamp
 *  
 *  where di is the i-th dimension (double fp number). The timestamps are
 *  expressed in microseconds starting from 0.
 *  
 *  This class is defined in \ref Pane_Farming/src/generator.cpp
 */
class DatasetReader {
private:
	ifstream *file;
	size_t num_tuples; // total number fo tuples in the dataset
	size_t next_id; // next tuple to read
	double *tuples; // array of data (each tuple is DIM+2 contiguous doubles)

	// method to count the number of tuples in the file
	inline unsigned long countNoTuples() {
		// we have a tuple for each line in the text file
		string line;
		size_t n_lines = 0;
		while(getline(*file, line)) n_lines++;
		file->clear();
		file->seekg(0, ios::beg); // restart from the beginning
		return n_lines;
	}

	// method to read all the tuples from the file
	inline size_t readTuples() {
		size_t id_tuple = 0;
		// read the file row by row
		string line;
		while(getline(*file, line)) {
			int i = 0;
   			stringstream s(line);
   			while(!s.eof()) {
      			string tmp;
      			s >> tmp;
      			tuples[(id_tuple*(DIM+2))+i] = stod(tmp);
      			i++;
   			}
   			id_tuple++;
   		}
   		return id_tuple;
	}

public:
	// constructor
	DatasetReader(char *filename) {
		file = new ifstream(filename);
		num_tuples = countNoTuples();
		tuples = new double[num_tuples * (DIM+2)];
		readTuples();
		next_id = 0;
	}

	// destructor
	~DatasetReader() {
		file->close();
		delete file;
		delete[] tuples;
	}

	// method to get the number of tuples
	size_t getNoTuples() {
		return num_tuples;
	}

	// method to get the next tuple
	double *nextTuple() {
		if(next_id < num_tuples) {
			double *p;
			p = &(tuples[next_id * (DIM+2)]);
			next_id++;
			return p;
		}
		else return nullptr;
	}
};

/*! 
 *  \class DatasetReaderBin
 *  
 *  \brief Class to Parse the Dataset Binary File
 *  
 *  This class is used by the Generator to read the dataset of tuples from
 *  a binary file.
 *  
 *  This class is defined in \ref Pane_Farming/src/generator.cpp
 */
class DatasetReaderBin {
private:
	int fd;
	size_t num_tuples; // total number fo tuples in the dataset
	size_t next_id; // next tuple to read
	double *map; // mapped file as an array of doubles (each tuple is DIM+2 contiguous doubles)
	struct stat stbuf;
	size_t file_size; // size in bytes of the mapped file

	// method to count the number of tuples in the file
	inline size_t countNoTuples() {
    	if((fstat(fd, &stbuf) != 0) || (!S_ISREG(stbuf.st_mode))) {
    		perror("Error fstat on the dataset binary file");
    		exit(-1);
		}
		size_t file_size = stbuf.st_size;
		return file_size / (sizeof(double) * (DIM+2));
	}

	// method to read all the tuples from the file
	inline size_t readTuples() {
		// stub
	}

public:
	// constructor
	DatasetReaderBin(char *filename) {
		fd = open(filename, O_RDONLY);
		if(fd == -1) {
			perror("Error opening the dataset binary file");
    		exit(-1);
		}
		num_tuples = countNoTuples();
		next_id = 0;
		// map the file in memory
		file_size = sizeof(double) * (DIM+2) * num_tuples;
		map = (double *) mmap(0, file_size, PROT_READ, MAP_SHARED, fd, 0);
		if(map == MAP_FAILED) {
			close(fd);
			perror("Error mmapping the file");
			exit(-1);
    	}
	}

	// destructor
	~DatasetReaderBin() {
		if(munmap(map, file_size) == -1) {
			perror("Error un-mmapping the dataset binary file");
    	}
    	close(fd);
	}

	// method to get the number of tuples
	inline size_t getNoTuples() {
		return num_tuples;
	}

	// method to get the next tuple
	inline double* nextTuple() {
		if(next_id < num_tuples) {
			double *p;
			p = &(map[next_id * (DIM+2)]);
			next_id++;
			return p;
		}
		else return nullptr;
	}
};

// main
int main(int argc, char *argv[]) {
	int option = 0;
	cpu_set_t cpuset;
	size_t port = 10000;
	char *hostname = NULL;
	char *file_dataset = NULL;
	/* Max number of tuples to generate, if -1 we generate
	   all the tuples of the input dataset. */
	int max_tuples = -1;
	// arguments from command line
	if(argc == 7) {
		while ((option = getopt(argc, argv, "h:p:f:")) != -1) {
			switch (option) {
				case 'h': hostname = optarg;
					break;
				case 'p': port = atoi(optarg);
					break;
				case 'f': file_dataset = optarg;
					break;
				default: {
					cout << argv[0] << " -h <hostname/ip> -p <port> -f <dataset_file> [-n <max no. of tuples>]" << endl;
					exit(0);
				}
        	}
    	}
	}
	else if(argc == 9) {
		while ((option = getopt(argc, argv, "h:p:f:n:")) != -1) {
			switch (option) {
				case 'h': hostname = optarg;
					break;
				case 'p': port = atoi(optarg);
					break;
				case 'f': file_dataset = optarg;
					break;
				case 'n': max_tuples = atoi(optarg);
					break;
				default: {
					cout << argv[0] << " -h <hostname/ip> -p <port> -f <dataset_file> [-n <max no. of tuples>]" << endl;
					exit(0);
				}
        	}
    	}
	}
	else {
		cout << argv[0] << " -h <hostname/ip> -p <port> -f <dataset_file> [-n <max no. of tuples>]" << endl;
		exit(0);
	}
    cout << "The frequency of the CPU is " << FREQ << " Mhz" << endl;
    // set affinity of the process on a specific thread context
	CPU_ZERO(&cpuset);
	#if (defined(PIANOSA) || defined(PIANOSAU))
    	CPU_SET(31, &cpuset);
	#endif
	#if defined(REPARA)
    	CPU_SET(47, &cpuset);
	#endif
	#if defined(REPHRASE)
    	CPU_SET(159, &cpuset);
	#endif
	pthread_t tid = pthread_self();
	if(pthread_setaffinity_np(tid, sizeof(cpu_set_t), &cpuset)) {
		perror("Generator process cannot be binded on the chosen CPU\n");
	}
	else cout << "Generator process is executed on CPU " << sched_getcpu() << endl;
	// read the dataset
#if 0
	// dataset read in text mode
	DatasetReader dataset(file_dataset);
#else
	// dataset read in binary mode
	DatasetReaderBin dataset(file_dataset);
#endif
	cout << "Dataset acquired successfully" << endl;
	// connect to the pane farming pattern
	int socket = socket_connect(hostname, port);
	if(socket < 0) {
		cout << "Connection to Pane Farming: refused" << endl;
		exit(-1);
	}
	else cout << "Connection to Pane Farming: OK" << endl;
	// send all the tuples in the dataset following their generator timestamps
	double *tuple;
	ticks start_ticks = getticks();
	size_t tuples_to_receive = 0; // no. of tuples to read from the dataset
	if(max_tuples == -1) tuples_to_receive = dataset.getNoTuples();
	else tuples_to_receive = max_tuples;
	volatile ticks gen_starting_ticks = getticks();
	for(size_t i=0; i<tuples_to_receive; i++) {
		// get next tuple (already serialized)
		tuple = dataset.nextTuple();
		// wait according to the generator timestamps of the tuple
		volatile ticks end_ticks = start_ticks+FROM_USECS_TO_TICKS(tuple[DIM+1]); // DIM+1 is the last field
		volatile ticks cycles = getticks();
		while(cycles < end_ticks) cycles = getticks();
		// we can send the tuple
		if(socket_send(socket, (char *) tuple, tuple_t::getSize()) != tuple_t::getSize()) {
			cout << "Error in transmitting a tuple" << endl;
			exit(-1);
		}
	}
	double gen_time_us = FROM_TICKS_TO_USECS(getticks() - gen_starting_ticks);
	cout << "Generation done in " << gen_time_us/1000000 << " seconds" << endl;
	return 0;
}
