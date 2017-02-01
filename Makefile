# --------------------------------------------------------------------------
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License version 2 as
#  published by the Free Software Foundation.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#  
#  As a special exception, you may use this file as part of a free software
#  library without restriction.  Specifically, if other files instantiate
#  templates or use macros or inline functions from this file, or you compile
#  this file and link it with other files to produce an executable, this
#  file does not by itself cause the resulting executable to be covered by
#  the GNU General Public License.  This exception does not however
#  invalidate any other reasons why the executable file might be covered by
#  the GNU General Public License.
# ---------------------------------------------------------------------------

# Author: Gabriele Mencagli <mencagli@di.unipi.it>
# Date:   March 2016

FF_ROOT		= $(HOME)/fastflow

CXX			= g++
INCLUDES	= -I $(FF_ROOT) -I $(PWD)/include
CXXFLAGS  	= -std=c++11

LDFLAGS 	= -pthread
MACROS		= -DNDEBUG -DTRACE_FASTFLOW -DLOG -DREPHRASE -DDIM=8
OPTFLAGS	= -g -O3 -finline-functions

TARGETS		= generator plq_skyline_pb plq_skyline_tb plq_skyline_as plq_skyline_aspid pf_skyline pf_skyline_od pf_topd pf_topd_od wf_skyline wf_topd

.DEFAULT_GOAL := all
.PHONY: all clean cleanall
.SUFFIXES: .cpp

generator: src/generator.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) $(OPTFLAGS) -o bin/generator src/generator.cpp $(LDFLAGS)

plq_skyline_pb: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DSKYLINE -DPLQ_ONLY -DPB_RR_SCHED $(OPTFLAGS) -o bin/plq_skyline_pb src/pf_query.cpp $(LDFLAGS)

plq_skyline_pb_iop: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DSKYLINE -DPLQ_ONLY -DIOP -DPB_RR_SCHED $(OPTFLAGS) -o bin/plq_skyline_pb_iop src/pf_query.cpp $(LDFLAGS)

plq_skyline_tb: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DSKYLINE -DPLQ_ONLY -DTB_RR_SCHED $(OPTFLAGS) -o bin/plq_skyline_tb src/pf_query.cpp $(LDFLAGS)

plq_skyline_tb_iop: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DSKYLINE -DPLQ_ONLY -DIOP -DTB_RR_SCHED $(OPTFLAGS) -o bin/plq_skyline_tb_iop src/pf_query.cpp $(LDFLAGS)

plq_skyline_as: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DSKYLINE -DPLQ_ONLY -DADAPTIVE_SCHED $(OPTFLAGS) -o bin/plq_skyline_as src/pf_query.cpp $(LDFLAGS)

plq_skyline_as_iop: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DSKYLINE -DPLQ_ONLY -DIOP -DADAPTIVE_SCHED $(OPTFLAGS) -o bin/plq_skyline_as_iop src/pf_query.cpp $(LDFLAGS)

plq_skyline_aspid: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DSKYLINE -DPLQ_ONLY -DADAPTIVE_SCHED_PID $(OPTFLAGS) -o bin/plq_skyline_aspid src/pf_query.cpp $(LDFLAGS)

plq_skyline_aspid_iop: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DSKYLINE -DPLQ_ONLY -DIOP -DADAPTIVE_SCHED_PID $(OPTFLAGS) -o bin/plq_skyline_aspid_iop src/pf_query.cpp $(LDFLAGS)

plq_topd_pb: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DTOPD -DPLQ_ONLY -DPB_RR_SCHED $(OPTFLAGS) -o bin/plq_topd_pb src/pf_query.cpp $(LDFLAGS)

plq_topd_pb_iop: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DTOPD -DPLQ_ONLY -DIOP -DPB_RR_SCHED $(OPTFLAGS) -o bin/plq_topd_pb_iop src/pf_query.cpp $(LDFLAGS)

plq_topd_tb: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DTOPD -DPLQ_ONLY -DTB_RR_SCHED $(OPTFLAGS) -o bin/plq_topd_tb src/pf_query.cpp $(LDFLAGS)

plq_topd_tb_iop: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DTOPD -DPLQ_ONLY -DIOP -DTB_RR_SCHED $(OPTFLAGS) -o bin/plq_topd_tb_iop src/pf_query.cpp $(LDFLAGS)

plq_topd_as: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DTOPD -DPLQ_ONLY -DADAPTIVE_SCHED $(OPTFLAGS) -o bin/plq_topd_as src/pf_query.cpp $(LDFLAGS)

plq_topd_as_iop: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DTOPD -DPLQ_ONLY -DIOP -DADAPTIVE_SCHED $(OPTFLAGS) -o bin/plq_topd_as_iop src/pf_query.cpp $(LDFLAGS)

plq_topd_aspid: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DTOPD -DPLQ_ONLY -DADAPTIVE_SCHED_PID $(OPTFLAGS) -o bin/plq_topd_aspid src/pf_query.cpp $(LDFLAGS)

plq_topd_aspid_iop: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DTOPD -DPLQ_ONLY -DIOP -DADAPTIVE_SCHED_PID $(OPTFLAGS) -o bin/plq_topd_aspid_iop src/pf_query.cpp $(LDFLAGS)

pf_skyline: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DSKYLINE -DADAPTIVE_SCHED_PID $(OPTFLAGS) -o bin/pf_skyline src/pf_query.cpp $(LDFLAGS)

pf_skyline_od: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DSKYLINE -DADAPTIVE_SCHED_PID -DWLQ_OD $(OPTFLAGS) -o bin/pf_skyline_od src/pf_query.cpp $(LDFLAGS)

pf_topd: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DTOPD -DADAPTIVE_SCHED_PID $(OPTFLAGS) -o bin/pf_topd src/pf_query.cpp $(LDFLAGS)

pf_topd_od: src/pf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DTOPD -DADAPTIVE_SCHED_PID -DWLQ_OD $(OPTFLAGS) -o bin/pf_topd_od src/pf_query.cpp $(LDFLAGS)

wf_skyline: src/wf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DSKYLINE $(OPTFLAGS) -o bin/wf_skyline src/wf_query.cpp $(LDFLAGS)

wf_topd: src/wf_query.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MACROS) -DTOPD $(OPTFLAGS) -o bin/wf_topd src/wf_query.cpp $(LDFLAGS)

all: $(TARGETS)

clean:
	rm -f bin/*

cleanall:
	\rm -f bin/*.o bin/*~
