/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file pane_skyline_checker.cpp
 *  \brief Utility for testing the correctness of the SkyPane implementation
 *  
 *  This program tests the correctness of the skyline pane. The idea is
 *  to compare the results of the algorithm in the SkyPane class with the ones
 *  using the R-tree data structure from the SpatialIndex library.
 *  
 *  \note This process communicates with the Generator through a TCP/IP socket.
 *        This tool works correctly only with dataset files with integer
 *        dimensions for the tuples.
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

// include
#include <stdio.h>
#include <unistd.h>
#include <spatialindex/SpatialIndex.h>
#include <spatialindex/Region.h>
#include <queue>
#include <iostream>
#include <queries.hpp>
#include <socket_utils.hpp>

// define
#define FILLFACTOR 0.4
#define INDEXCAPACITY 5
#define LEAFCAPACITY 10
#define MAXCORNER 1
#define MINCORNER 2
#define DOMINATES 0
#define ISDOMINATED 1

using namespace std;
using namespace SpatialIndex;

/*! 
 *  \class SkyPoint
 *  
 *  \brief Skyline Point
 *  
 *  Class that implements the generic skyline point used by the SpatialIndex
 *  library.
 *  
 *  This class is defined in \ref Pane_Farming/include/pane_skyline_checker.cpp
 */
struct SkyPoint {
public:
	// tuple encapsulated into the SkyPoint
	tuple_t t;
	// region describing a bounding box, in our case is simply a point
	Region *pRegion;

	// constructor
	SkyPoint(tuple_t &_t): t(_t) {
		// initialization of the pRegion
		pRegion = new Region(t.d, t.d, DIM);
	}

	// destructor
	~SkyPoint(){
		delete pRegion;
	}

	// check the dominance
	inline bool dominance(const Region &other, int relation, int what) {
		return dominance((what == MINCORNER) ? other.m_pLow : other.m_pHigh, relation);
	}

	// check the dominance
	inline bool dominance(const SkyPoint &other, int relation) {
		return dominance((other.t).d, relation);
	}

	// check the dominance
	inline bool dominance(const double *pcoords, int relation) {
		int tocheck, a, b;
		bool sless = false;
		for(size_t i=0; i<DIM; i++) {
			a = (relation == ISDOMINATED) ? pcoords[i] : t.d[i];
			b = (relation == ISDOMINATED) ? t.d[i] : pcoords[i];
			if(a > b) return false;
			if(a < b) sless = true;
		}
		return sless;
	}
};

/*! 
 *  \class GenericQuery
 *  
 *  \brief Class of a Generic Query to the SkyRTree.
 *  
 *  Class that implements a generic query of the SkyRTree.
 *  
 *  This class is defined in \ref Pane_Farming/include/pane_skyline_checker.cpp
 */
class GenericQuery : public SpatialIndex::IQueryStrategy {
protected:
	bool stop = false;
	queue<id_type> ids;

	// constructor (protected)
	GenericQuery() {}

	// methods to be implemented
	virtual void handleIndex(id_type id, Region *region) = 0;
	virtual void handleLeaf(SkyPoint * candidate) = 0;

public:
	// method to get the next entry
	void getNextEntry(const IEntry &entry, id_type &nextEntry, bool &hasNext) {
		const INode *n = dynamic_cast<const INode*>(&entry);
		//IShape * ps;
		Region *region;
		if(n->isLeaf()) {
			for(uint32_t child = 0; child < n->getChildrenCount(); child++) {
				SkyPoint *candidate = (SkyPoint *) n->getChildIdentifier(child);
				handleLeaf(candidate);
			}
		}
		else {
			for(uint32_t child=0; child < n->getChildrenCount(); child++) {
				n->getChildRegion(child, &region);
				handleIndex(n->getChildIdentifier(child), region);
			}
		}
		if(!ids.empty() && !stop) {
			nextEntry = ids.front();
			ids.pop();
			hasNext = true;
		}
		else hasNext = false;
	}
};

/*! 
 *  \class DominanceQuery
 *  
 *  \brief Dominance Query
 *  
 *  Class that implements a query to find all the tuples dominated by the input one.
 *  
 *  This class is defined in \ref Pane_Farming/include/pane_skyline_checker.cpp
 */
class DominanceQuery : public GenericQuery {
private:
	SkyPoint * pt;
	int mode;
	int corner;
public:
	list<SkyPoint*> pts;

	// constructor
	DominanceQuery(SkyPoint *pt, int mode, int corner): GenericQuery(), pt(pt), mode(mode), corner(corner) {}

	// method to handle the index
	void handleIndex(id_type id, Region *region) {
		if(pt->dominance(*region, mode, corner)) ids.push(id);
	}

	// method to handle a leaf
	void handleLeaf(SkyPoint *candidate) {
		if (pt->dominance(*candidate, mode)) pts.push_back(candidate);
	}
};

/*! 
 *  \class SkylineQuery
 *  
 *  \brief SkylineQuery Query
 *  
 *  Class that implements a query to find the skyline in the SkyRTree.
 *  
 *  This class is defined in \ref Pane_Farming/include/pane_skyline_checker.cpp
 */
class SkylineQuery: public GenericQuery {
private:
	SkyPoint *pt;

public:
	bool isSkyline = true;

	// constructor
	SkylineQuery(SkyPoint *pt): GenericQuery(), pt(pt) {}

	// method to handle the index
	void handleIndex(id_type id, Region *region) {
		if(pt->dominance(*region, ISDOMINATED, MAXCORNER)) {
			isSkyline=false;
			stop=true;
		}
		else if(pt->dominance(*region, ISDOMINATED, MINCORNER)) ids.push(id);
	}

	// method to handle a leaf
	void handleLeaf(SkyPoint *candidate) {
		if(pt->dominance(*candidate, ISDOMINATED)) {
			isSkyline=false;
			stop=true;
		}
	}
};

/*! 
 *  \class SkyRTree
 *  
 *  \brief R-Tree Index for Skyline Queries
 *  
 *  Class that implements a R-Tree index structure to perform Skyline Queries.
 *  It uses the SpatialIndex library.
 *  
 *  This class is defined in \ref Pane_Farming/include/pane_skyline_checker.cpp
 */
class SkyRTree {
private:
	IStorageManager *storage;
	ISpatialIndex *rtree;
	// id_type is a typedef for int_64_t
	id_type indexIdentifier;
	size_t size = 0;
	//Region mbr;

public:
	// constructor
	SkyRTree() {
		storage = StorageManager::createNewMemoryStorageManager();
		rtree = RTree::createNewRTree(*storage, FILLFACTOR, INDEXCAPACITY, LEAFCAPACITY, DIM, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
	}

	// destructor
	~SkyRTree() {
		delete rtree;
		delete storage;
	}

	// method to add a point
	void insertPoint(SkyPoint *x) {
		size++;
		rtree->insertData(0, 0, *(x->pRegion), (id_type) x);
	}

	// method to remove a point
	void removePoint(SkyPoint *x) {
		size--;
		rtree->deleteData(*(x->pRegion), (id_type) x);
	}

	// method to return the current size of the R-tree
	size_t getSize() const {
		return size;
	}

	// method to return the list of points dominated by x
	list<SkyPoint *> getDominatedPoints(SkyPoint *x) {
		DominanceQuery qs(x, DOMINATES, MAXCORNER);
		rtree->queryStrategy(qs);
		return qs.pts;
	}

	// method to remove the points in the R-Tree dominated by x
	size_t removeDominatedPoints(SkyPoint *x) {
		size_t removed = 0;
		// get the dominated points
		list<SkyPoint *> removed_list = this->getDominatedPoints(x);
		// remove the dominated points
		for(list<SkyPoint *>::iterator it=removed_list.begin() ; it != removed_list.end(); ++it) {
			this->removePoint(*it);
			removed++;
		}
		return removed;
	}

	// method to check if a point x is in the skyline
	bool isSkyline(SkyPoint * x) {
		SkylineQuery qs(x);
		rtree->queryStrategy(qs);
		return qs.isSkyline;
	}
};

// method used to print a set of tuples in the same order
bool operator<(const tuple_t &t1, const tuple_t &t2) {
	for(size_t i=0; i<DIM; i++) {
		if(t1.d[i] < t2.d[i]) return true;
		if(t1.d[i] > t2.d[i]) return false;
	}
	if(t1.app_ts < t2.app_ts) return true;
	if(t1.app_ts > t2.app_ts) return false;
	if(t1.gen_ts < t2.gen_ts) return true;
	if(t1.gen_ts > t2.gen_ts) return false;
	return false;
}

// main
int main(int argc, char *argv[]) {
	int option = 0;
	size_t port = 10000;
	size_t start_idx = 0;
	size_t end_idx = 0;
	int socket;
	SkyPane pane;
	SkyRTree db;
	tuple_t tuple_min;
	for(size_t i=0; i<DIM; i++) tuple_min.d[i] = -10000000; // large enough negative value
	SkyPoint point_min(tuple_min);
	// arguments from command line
	if(argc != 7) {
		cout << argv[0] << " -p <port> -s <starting index> -e <ending index>" << endl;
		exit(0);
	}
	while ((option = getopt(argc, argv, "p:s:e:")) != -1) {
		switch (option) {
			case 'p': port = atoi(optarg);
				break;
			case 's': start_idx = atoi(optarg);
				break;
			case 'e': end_idx = atoi(optarg);
				break;
			default: {
				cout << argv[0] << " -p <port> -s <starting index> -e <ending index>" << endl;
				exit(0);
			}
        }
    }
	// wait for the connection from the generator
	int *s = socket_accept(1, port);
	socket = *s;
	size_t rcv_count = 0;
	tuple_t t;
	char tuple_buffer[tuple_t::getSize()];
	// receive at most num_tuples tuples from the generator
	while(socket_receive(socket, tuple_buffer, tuple_t::getSize()) > 0) {
		rcv_count++;
		if(rcv_count < start_idx) continue;
		t.deserialize(tuple_buffer);
		pane.addComputeTuple(t);
		SkyPoint *p = new SkyPoint(t);
		if(db.removeDominatedPoints(p) > 0) {
			db.insertPoint(p);
		}
		else {
			if(db.isSkyline(p)) db.insertPoint(p);
		}
		if(rcv_count > end_idx) break;
	}
	// obtain the two results
	vector<tuple_t> *result1 = pane.toVector();
	list<SkyPoint *> removed = db.getDominatedPoints(&point_min);
	// if they have different sizes return false
	bool areEqual = true;
	if(result1->size() != removed.size()) {
		cout << "False" << endl;
		cout << "(My) size: " << result1->size() << " (R-Tree) size: " << removed.size() << endl;
		/*cout << "Result 1" << endl;
		for(auto const &t: *result1) t.print();
		cout << "Result 2" << endl;
		for(auto const &p: removed) (p->t).print();*/
		exit(0);
	}
	else {
		// sort the first result
		std::sort ((*result1).begin(), (*result1).end());
		// sort the second result
		vector<tuple_t> *result2 = new vector<tuple_t>();
		list<SkyPoint *>::const_iterator iterator;
		for (iterator = removed.begin(); iterator != removed.end(); ++iterator) {
    		result2->push_back((*iterator)->t);
		}
		std::sort ((*result2).begin(), (*result2).end());
		vector<tuple_t>::iterator compareIterator1 = result1->begin();
		vector<tuple_t>::iterator compareIterator2 = result2->begin();
		// compare the two results
		for(size_t i=0; i<result1->size(); i++) {
			if(*compareIterator1 != *compareIterator2) {
				areEqual=false;
				break;
			}
			//(*compareIterator1).print();
			//(*compareIterator2).print();
			compareIterator1++;
			compareIterator2++;
		}
	}
	if(areEqual) {
		cout << "True" << endl;
		cout << "Skyline size: " << result1->size() << endl;
	}
	else {
		cout << "False" << endl;
		cout << "(My) size: " << result1->size() << " (R-Tree) size: " << removed.size() << endl;
		/* cout << "Result 1" << endl;
		for(auto const &t: *result1) t.print();
		cout << "Result 2" << endl;
		for(auto const &p: removed) (p->t).print(); */
	}	
	return 0;
}
