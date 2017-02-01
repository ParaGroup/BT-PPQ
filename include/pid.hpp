/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*! 
 *  \file pid.hpp
 *  \brief PID Controller
 *  
 *  This file implements a PID (Proportional-Integrative-Derivative) controller.
 *  The approach is inspired to the implementation in the Arduino PID Library as
 *  described at the following web page:
 *  
 *  http://brettbeauregard.com/blog/2011/04/improving-the-beginners-pid-introduction/
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

#ifndef PID_H
#define PID_H

// include
#include <cmath>
#include <general.hpp>

// define
#define MANUAL 0
#define AUTOMATIC 1

using namespace std;

/* This file provides the following classes:
 *  -PID: class that implements a PID controller;
 *  -PIDAutotuning: class that implements an auto-tuning mechanism for PIDs.
 */

/*! 
 *  \class PID
 *  
 *  \brief PID Controller
 *  
 *  This file implements a PID (Proportional-Integrative-Derivative) controller.
 *  The approach is inspired to the implementation in the Arduino PID Library.
 *  
 *  This class is defined in \ref Pane_Farming/include/pid.hpp
 */
class PID {
private:
	double setpoint; // desired value of the system output
	double iTerm = 0, lastInput = -1;
	double input = 0, output = 0; // current input and output of the PID
	double kp, ki, kd; // gains
	double outMin, outMax; // boundaries of the controller output
	bool inAuto; // true if the controller is enabled, false otherwise
	bool isReversed; // false=DIRECT, true=REVERSE

public:
	// constructor
	PID(double _kp, double _ki, double _kd, double _setpoint, double _outMin, double _outMax, bool _inAuto=true, bool _isReversed=false) {
		isReversed = _isReversed;
		setGains(_kp, _ki, _kd);
		setpoint = _setpoint;
		outMin = _outMin;
		outMax = _outMax;
		inAuto = _inAuto;
	}

	// destructor
	~PID() {}

	// method to initialize the PID
	void initialize() {
   		lastInput = input;
   		iTerm = output;
   		if(iTerm > outMax) iTerm = outMax;
   		else if(iTerm < outMin) iTerm = outMin;
	}

	/* 
	 * PID behavior: setInput() -> process() -> getOutput().
	 * Set the current input of the controller.
	 */
	void setInput(double _input) {
		input = _input;
	}

	/* 
	 * PID behavior: setInput() -> process() -> getOutput().
	 * Return true if the controller computes a new output from the present input, false otherwise.
	 */
	bool process() {
		if(lastInput == -1) lastInput = input; // only for the first sampling interval
		// if the PID is disabled return false
		if(!inAuto) return false;
		else {
			// the controller is enabled
      		double error = setpoint-input;
      		iTerm += (ki * error);
      		if(iTerm > outMax) iTerm = outMax;
      		else if(iTerm < outMin) iTerm = outMin;
      		double dInput = (input - lastInput);
      		// compute the PID output
      		output = kp * error + iTerm - kd * dInput;
      		if(output > outMax) output = outMax;
      		else if(output < outMin) output = outMin;
      		// save the current input as the last one
      		lastInput = input;
      		//cout << "Setpoint: " << setpoint << " sampled coverage: " << input << " error: " << error << " new output: " << output << endl;
      		return true;
		}
	}

	/* 
	 * PID behavior: setInput() -> process() -> getOutput().
	 * Get the current output from the controller.
	 */
	double getOutput() {
		return output;
	}

	// method to set the gain constants of the PID
	void setGains(double _kp, double _ki, double _kd) {
		if (_kp<0 || _ki<0 || _kd<0) return;
		if(!isReversed) {
			kp = _kp;
			ki = _ki;
			kd = _kd;
		}
		else {
      		kp = (0 - _kp);
      		ki = (0 - _ki);
      		kd = (0 - _kd);
		}
	}

	// method to set the output boundaries
	void setOutputLimits(double min, double max) {
   		if(min > max) return;
   		outMin = min;
   		outMax = max;
   		if(output > outMax) output = outMax;
   		else if(output < outMin) output = outMin;
   		if(iTerm > outMax) iTerm = outMax;
   		else if(iTerm < outMin) iTerm = outMin;
	}

	// method to set the mode of the controller
	void setMode(int mode) {
    	bool newAuto = (mode == AUTOMATIC);
    	if(newAuto == !inAuto) {
    		// we go from manual to auto
        	initialize();
    	}
    	inAuto = newAuto;
	}

	// method to set the new setpoint
	void setSetPoint(double _s) {
		setpoint = _s;
		iTerm = 0, lastInput = -1;
		input = 0, output = 0;
	}

	// method to control the PID direction
	void setReversed(bool _isReversed) {
		isReversed = _isReversed;
	}
};

/*! 
 *  \class PIDAutotuning
 *  
 *  \brief Auto-tuning Mechanism (Relay-based) for PIDs
 *  
 *  This file implements an auto-tuning mechanism (Relay Method) for PIDs. This method has been
 *  ripped from the Arduino Auto-tuning Library.
 *  
 *  This class is defined in \ref Pane_Farming/include/pid.hpp
 */
class PIDAutotuning {
private:
	bool isMax, isMin;
	double input, output; // current input and output
	double setpoint;
	double noiseBand; // noise band on the controller's input (system output)
	int controlType;
	bool running;
	unsigned long peak1, peak2;
	int sampleTime = SAMPLING_PUNCTUATION_USEC/1000; // sampling time in ms
	int nLookBack;
	int peakType;
	double lastInputs[101];
    double peaks[10];
	int peakCount;
	bool justchanged;
	bool justevaled;
	double absMax, absMin;
	double oStep;
	double outputStart;
	double Ku, Pu;
	bool isReversed;

	// method to finalize the tuning phase
	void finishUp() {
	  output = outputStart;
      // we can generate tuning parameters!
      Ku = 4*(2*oStep)/((absMax-absMin)*3.14159);
      Pu = (double)(peak1-peak2) / 1000;
	}

public:
	// construtor
    PIDAutotuning(double _setpoint, bool _isReversed=false) {
    	setpoint = _setpoint;
		controlType = 1 ; //default to PID (1) otherwise PI (0)
		noiseBand = 0;
		nLookBack = 4;
		running = false;
		oStep = 0.15;
		isReversed = _isReversed;
		//SetLookbackSec(10);
	}

	// destructor
	~PIDAutotuning() {}

	// method to stop the auto-tuning mechanism
	void cancel() {
		running = false;
	}

	/* 
	 * PIDAutotuning behavior: setInput() -> tuning() -> getOutput().
	 * Set the current input of the auto-tuning mechanism.
	 */
	void setInput(double _input) {
		input = _input;
	}

	/* 
	 * PIDAutotuning behavior: setInput() -> tuning() -> getOutput().
	 * Tuning method (return true if the tuning phase is complete, false otherwise).
	 */
	bool tuning() {
		justevaled = false;
		if(peakCount > 9 && running) {
			running = false;
			finishUp();
			return true;
		}
		volatile unsigned long now = (unsigned long) FROM_TICKS_TO_USECS(getticks())/1000; // in ms
		double refVal = input;
		justevaled = true;
		if(!running) {
			// initialize working variables the first time around
			peakType = 0;
			peakCount = 0;
			justchanged = false;
			absMax = refVal;
			absMin = refVal;
			setpoint = refVal;
			running = true;
			outputStart = output;
			output = outputStart+oStep;
		}
		else {
			if(refVal > absMax) absMax = refVal;
			if(refVal < absMin) absMin = refVal;
		}
		// oscillate the output base on the input's relation to the setpoint
		if(!isReversed) {
			if(refVal > setpoint+noiseBand) output = outputStart-oStep;
			else if(refVal < setpoint-noiseBand) output = outputStart+oStep;
		}
		else {
    		if(refVal > setpoint+noiseBand) output = outputStart+oStep;
    		else if(refVal < setpoint-noiseBand) output = outputStart-oStep;
		}
  		// bool isMax = true, isMin = true;
  		isMax = true;
  		isMin = true;
  		// id peaks
  		for(size_t i=nLookBack-1; i>=0; i--) {
    		double val = lastInputs[i];
    		if(isMax) isMax = refVal>val;
    		if(isMin) isMin = refVal<val;
    		lastInputs[i+1] = lastInputs[i];
  		}
  		lastInputs[0] = refVal;
  		/*if(nLookBack < 9) {
  			// we don't want to trust the maxes or mins until the inputs array has been filled
			return false;
		}*/
  		if(isMax) {
    		if(peakType == 0) peakType = 1;
    		if(peakType == -1) {
      			peakType = 1;
      			justchanged = true;
      			peak2 = peak1;
    		}
    		peak1 = now;
    		peaks[peakCount] = refVal;
  		}
  		else if(isMin) {
    		if(peakType == 0) peakType = -1;
    		if(peakType == 1) {
      			peakType = -1;
      			peakCount++;
      			justchanged = true;
      		}
    		if(peakCount < 10) peaks[peakCount] = refVal;
  		}
  		if(justchanged && peakCount > 2) {
  			// we've transitioned: check if we can autotune based on the last peaks
    		double avgSeparation = (abs(peaks[peakCount-1]-peaks[peakCount-2])+abs(peaks[peakCount-2]-peaks[peakCount-3]))/2;
    		if(avgSeparation < 0.05*(absMax-absMin)) {
				finishUp();
      			running = false;
	  			return true;
    		}
  		}
   		justchanged = false;
		return false;
	}

	/* 
	 * PIDAutotuning behavior: setInput() -> tuning() -> getOutput().
	 * Get the current output from the auto-tuning mechanism.
	 */
	double getOutput() {
		return output;
	}

	// method to get the output step
	double getOutputStep() {
		return oStep;
	}

	// method to get the control type
	int getControlType() {
		return controlType;
	}

	// method to get the noise band
	double getNoiseBand() {
		// this should be acurately set
		return noiseBand;
	}

	// method to get how far back are we looking to identify peaks
	int getLookbackSec() {
		return nLookBack * sampleTime / 1000;
	}

	// method to get the proportional gain
	double getKp() {
		return controlType == 1 ? 0.6 * Ku : 0.4 * Ku;
	}

	// method to get the integrative gain
	double getKi() {
		return controlType == 1? 1.2*Ku / Pu : 0.48 * Ku / Pu;  // Ki = Kc/Ti
	}

	// method to get the derivative gain
	double getKd() {
		return controlType == 1? 0.075 * Ku * Pu : 0;  // Kd = Kc * Td
	}

	// method to set the control type
	void setControlType(int _type) {
		// 0 = PI, 1 = PID
		controlType = _type;
	}

	// method to set the noise band
	void setNoiseBand(double _band) {
		// the autotune will ignore signal chatter smaller than this value
		noiseBand = _band;
	}

	// method to set the output step
	void setOutputStep(double _step) {
		// how far above and below the starting value will the output step?
		oStep = _step;
	}

	// method to set how far back are we looking to identify peaks
	void setLookbackSec(int _value) {
    	if(_value < 1) _value = 1;
		if(_value < 25) {
			nLookBack = _value * 4;
			sampleTime = 250;
		}
		else {
			nLookBack = 100;
			sampleTime = _value*10;
		}
	}
};

#endif
