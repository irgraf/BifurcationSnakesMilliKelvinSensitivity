//
//  main.cpp
//  Thermosensing_Voltage_AP_Firing_Simple_Model_Scaling
//
//  Created by Isabella Graf
//

#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <cstdio>
#include <fstream>
#include <cmath>
#include "RandomNumberGenerator.h"


using namespace std;

random_device rd;           // to seed the pseudo-random number generator
unsigned seed = rd();       // seed for pseudo-random number generator
mt19937 mt(seed);           // pseudo-random number generator
const int NumberRuns = pow(2,21);   // number of recorded interspike times; for Fig. 5A, we chose NumberRuns = pow(2,21)
const vector<double> alphaVal={0.16768, 11.48};   // vector of alpha values for which interspike times are recorded
const int LengthalphaVal=2;   // length of alphaVal
const double DeltaTime=0.001;   // time step
const double UpperThresholdVoltage =500;   // threshold when AP is fired
const double LowerThresholdVoltage=-500;   // starting value for (rescaled) voltage
const bool PrintFile=0;   // set to 1 if interspike times should be saved in a file
const bool PrintSteps=0;   // set to 1 if intermediate steps should be printed to output device
const int SimulationNumber=1;   // identifies simulation - to avoid that file is overwritten if the simulation is run again
const string FileNameHelp = "_NumberRuns_"+to_string(NumberRuns)+"_TimeStep_"+to_string(double(DeltaTime))+"_UpperThreshold_" +to_string(int(UpperThresholdVoltage))+"_LowerThreshold_"+ to_string(int(LowerThresholdVoltage))+"_NumberAlphaVals_"+ to_string(int(LengthalphaVal))+"_SimNumber_"+ to_string(int(SimulationNumber));   // define file name
const string FileName = "SpikeTimes"+FileNameHelp+".txt";   // define file name

void OneSimulationStep (double& Voltage, const double& alpha){   // update voltage during one simulation step for a simulation with scaling variable alpha
    double noise = StandardNormalDist();   // draw noise according to a Standard Normal distribution
    if constexpr(PrintSteps){
        cout << "noise is " << noise << endl;
    }
    Voltage += DeltaTime*(alpha + pow(Voltage,2))+ sqrt(DeltaTime)*noise;   // voltage incremented by determinstic and stochastic part (first and second term, respectively)
}

void PrintToFile(const string& filename, const double& alpha, const vector<double>& SpikeTimes)  {   // print a vector of interspike times to a file with file name filename (and record scaling variable alpha)
    if(!PrintFile)
        return;
    ofstream file(filename, ios_base::app);
    if(file.is_open())  {
        cout << "opened file" << endl;
        file << "alpha is:" << endl;   // print value of scaling variable alpha
        file << alpha << endl;
        file << "spike times: " << endl;   // print vector of interspike times
        for (int i=0; i<SpikeTimes.size(); ++i) {
            file << SpikeTimes[i] << " ";
        }
        file << endl;
        file << endl;
    }
    else {
        cout << "could not open file/n";
    }
    file.close();
}

void PrepareFile(const string& filename)  {   // print basic information to file
    if(!PrintFile)
        return;
    ofstream file(filename);
    if(file.is_open())  {
        cout << "opened file for file preparation" << endl;
        file << "parameters are: " << endl;   // print parameters:
        file << "Number of runs per alpha are: " << NumberRuns << endl;   // number of recorded interspike times per value of scaling variable alpha
        file << "Time step is: " << DeltaTime << endl;   // time step
        file << "Lower Threshold value (starting value) is: " << LowerThresholdVoltage << endl;   // starting value of (rescaled) voltage
        file << "Upper Threshold value is: " << UpperThresholdVoltage << endl;   // value of (rescaled) voltage at which AP is fired
        file << endl;
    }
    else {
        cout << "could not open file/n";
    }
    file.close();
}


void ResetSystem (double& Volt, double& time){   // reset of system after firing of AP
    Volt = LowerThresholdVoltage;   // reset voltage back to starting value
    time=0.;   // reset time to zero
}

double SpikeTimeSingleRun (double& Voltage, double& time, const double& alpha){   // for specific value of alpha, record time until AP is fired (after system has been reset), i.e., the interspike time
    while (Voltage < UpperThresholdVoltage) {   // as long as voltage is below the threshold value:
        OneSimulationStep(Voltage, alpha);   // update voltage
        time += DeltaTime;   // update time
    }
    return time;   // return time when threshold is crossed first
}

void MeasureSpikeTimes (double& Voltage, double& time, vector<double>& spiketimes, const double& alpha){   // record NumberRuns interspike times
    cout << "measure spike times" << endl;
    for (int i=0; i<NumberRuns; ++i) {
        double SpikeTimeRun= SpikeTimeSingleRun(Voltage, time, alpha);   // calculate interspike time
        if constexpr (PrintSteps){
            cout << "time is " << time << endl;
        }
        spiketimes.push_back(SpikeTimeRun);   // record last interspike time (append to vector of interspike times)
        ResetSystem(Voltage, time);   // reset system
    }
    PrintToFile(FileName, alpha, spiketimes);   // print scaling variable and vector of interspike times to file
}



int main() {
    double CurrentVoltage=LowerThresholdVoltage;   // tracks voltage over time
    double CurrentTime=0.;   // tracks time
    vector<double> SpikeTimes;   // vector of interspike times
    SpikeTimes.reserve(NumberRuns+1);   // reserve space for vector of interspike times
    cout << "Upper threshold value is " << UpperThresholdVoltage << endl;
    cout << "Lower threshold value (starting value) is " << LowerThresholdVoltage << endl;
    PrepareFile(FileName);   // print basic information to file before start of simulation
    for (int j=0; j<alphaVal.size(); ++j) {   // for all values of the rescaling variable, determine interspike times
        cout << "current alpha value is " << alphaVal[j] << endl;
        MeasureSpikeTimes(CurrentVoltage, CurrentTime, SpikeTimes, alphaVal[j]);   // record interspike times and print vector to file
        SpikeTimes.clear();   // clear vector of interspike times for next value of alpha
    }
    
    return 0;
}



