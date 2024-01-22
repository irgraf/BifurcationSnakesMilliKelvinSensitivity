//
//  main.cpp
//  Thermosensing_Voltage_AP_Firing
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
const int VoltageUpdatesPerChannelOpening=pow(2,8);     // number of batches
const int SimulationSteps = VoltageUpdatesPerChannelOpening*pow(2,15);     // number of simulation steps
const int SimulationStepsStartRecording = 0.1*SimulationSteps;     // number of simulation steps after which recording starts
const int NumberRecordings = (SimulationSteps-SimulationStepsStartRecording)/VoltageUpdatesPerChannelOpening;     // number of recordings: here, record roughly once per channel opening
const int SimulationStepsBetweenRecording = (SimulationSteps-SimulationStepsStartRecording)/NumberRecordings;    // number of simulation steps between two recordings: here, corresponds to number of batches
const int NumberChannels = pow(2,23);       // total number of channels
const int NumberChannelsPerBatch = NumberChannels/VoltageUpdatesPerChannelOpening;       // number of channels in each batch
const double WidthSigmoidVoltage = 30;       // width of voltage sigmoid, DeltaV, in mV
const double WidthSigmoidTemp = 1;       // width of temperature sigmoid, DeltaT, in K
const double EquilibriumVoltage = -70;       // resting potential
const double ChannelOpeningTime = 7.*pow(10,-4);       // channel opening time, tau_c, in s
const double VoltageTimeScale = 0.00004*VoltageUpdatesPerChannelOpening/ChannelOpeningTime;       // membrane timescale tau_rest / (channel opening time tau_c/number of batches): corresponds to membrane timescale in units of timescale of one simulation step (=tauc/number of batches)
const double MaxVoltageChange= 3000;       // corresponds to (ic/cmem)*tau_rest, in mV
const double UpperThresholdVoltage = EquilibriumVoltage+MaxVoltageChange*(1.+0.5*sqrt(1.-4.*WidthSigmoidVoltage/MaxVoltageChange))/2.;        // threshold at which AP is fired, in mV
const double ChangeHalfVoltageReset = 0.01;       // if feedback for control parameter is on: amount by which half voltage is increased each time an AP is fired, d+, in mV
const double TimeBetweenSpikes = 0.12*VoltageUpdatesPerChannelOpening/ChannelOpeningTime;          // if feedback for control parameter is on: target time between two consecutive spikes (1/frequency) in units of timescale of one simulation step (=tauc/number of batches)
const double ChangeHalfVoltageDrift = ChangeHalfVoltageReset/TimeBetweenSpikes;      // if feedback for control parameter is on: amount by which half voltage is decreased every simulation step, d-, in mV
const double VoltageReset = EquilibriumVoltage;       // voltage after AP firing and reset
const bool PrintFile=1;       // set to 1 if data should be saved to file
const bool SingleTimeSeriesSimulation=1;       // set to 1 if a single time series should be recorded
const bool FixedHalfVoltageSweep = 0;      // set to 1 if spike times should be recorded for different values of the voltage at half maximum (control parameter not subject to feedback)
const int NumberRuns = 1;     // how often simulation should be run for a single value of the sweep parameter (set to 1 for simulations in the paper)
double VoltageShiftTemperature =0.0;    // initial value for external shift in the voltage at half maximum; if temperature is changed: will be updated during the simulation (in units of voltage in mV)
double NewTemperature =0.001*WidthSigmoidVoltage/WidthSigmoidTemp;   // shift in the voltage at half maximum introduced in the middle of the simulation, in mV; a positive value corresponds to an increase in temperature (set to zero for Fig. 2a)
double HalfVoltageOffset=-0.0001*WidthSigmoidVoltage/WidthSigmoidTemp;     // for single time series simulation: distance of voltage at half maximum relative to the bifurcation (if VoltageShiftTemperature = 0); a positive value corresponds to an increase in the voltage at half maximum (i.e., towards the irregular regime)
const bool HalfVoltageFeebackAfterVoltageShift=1;     // for single time series simulation: set to 0 (1) if feedback after temperature shift is off (on)
const bool HalfVoltageFeebackBeforeVoltageShift=1;     // for single time series simulation: set to 0 (1) if feedback before temperature shift is off (on)
int SimStepShiftTemperature = SimulationStepsStartRecording+(SimulationSteps-SimulationStepsStartRecording)/2.;  // time (in units of simulation steps) when shift in temperature (or, equivalently, voltage at half maximum) is set to NewTemperature (here, in the middle of the recorded simulation time)
const vector<double> HalfVoltageOffsetSweep = {0.0359195, 0.0269396, 0.0179597, 0.00897987, 0., -0.00897987, -0.0179597, -0.0269396, -0.0359195};  // for FixedHalfVoltageSweep: HalfVoltageOffsetSweep corresponds to offset of voltage at half maximum with respect to the bifurcation; positive values correspond to an increase in the voltage at half maximum (i.e., towards the irregular regime)
const vector<double> SweepParameters = {146848., 88771.5, 58983.9, 42568.6, 32874.6, 26767.8, 22689.2, 19821.5, 17713.7};  // for FixedHalfVoltageSweep: corresponds to the expected average interspike time (in units of one simulation step, i.e., tau_c/number of batches); this number is not too relevant since it is used only to reserve an appropriate amount of space for the vector of interspike times. However, the vector has to have the same length as HalfVoltageOffsetSweep
const string FileNameHelp = "_"+to_string(NewTemperature)+"_SimSteps2e"+to_string(int(log2(SimulationSteps)))+"_numchannel2e" +to_string(int(log2(NumberChannels)))+"_ChangeHalfVol"+ to_string(ChangeHalfVoltageReset)+"_tauve"+to_string(int(log10(VoltageTimeScale)))+"_TimeBetwSpikese"+ to_string(int(log10(TimeBetweenSpikes))) +"_FeedbackAfter"+to_string(int(HalfVoltageFeebackAfterVoltageShift))+"_FeedbackBefore"+to_string(int(HalfVoltageFeebackBeforeVoltageShift))+ "VhalfOffset_" + to_string(double(HalfVoltageOffset))+ "VolUpdatesPerChannelOp_" + to_string(double(VoltageUpdatesPerChannelOpening))+ "tauc_" + to_string(double(ChannelOpeningTime));  // used to make file name for single time series
const string FileNameSingleTimeSeries = "timeseries" + FileNameHelp +"_1.txt";  // make file name for single time series

vector<int> ListNumberOpenChannelsBatches(VoltageUpdatesPerChannelOpening,NumberChannelsPerBatch*0.1);  // initialize vector listing the number of open channels in every batch

double OpenProb (double Voltage, double Vhalf, double Width, double Shift){ // channel opening probability as a function of the current value of the voltage, the voltage at half maximum, the width of the sigmoid, and a shift of the voltage at half maximum (e.g., due to a change in temperature)
    return 1./(1.+exp(-(Voltage-Vhalf+Shift)/Width));
}

void OneSimulationStep (double &Voltage, double &Vhalf, double & VoltageShift, vector<int> &VectorOpenChannels, int LabelBatch){     // one simulation step if there is feedback in the control parameter; arguments are the current value of the voltage, the current value of the voltage at half maximum, a potential shift in the voltage at half maximum (the temperature), the vector of the number of open channels for all batches, and the label of the batch that is updated this simulation step
    if (Voltage < UpperThresholdVoltage) { // if voltage is below the threshold
        VectorOpenChannels[LabelBatch]=BinomialDist(NumberChannelsPerBatch, OpenProb(Voltage, Vhalf, WidthSigmoidVoltage, VoltageShift)); // randomly redraw the number of open channels in batch LabelBatch, according to the current value of the opening probability
        Voltage += (EquilibriumVoltage-Voltage)/VoltageTimeScale + MaxVoltageChange/VoltageTimeScale/NumberChannels*accumulate(VectorOpenChannels.begin(), VectorOpenChannels.end(), 0.); // increase the voltage according to Eq. 2 in the main text for a timestep (tauc_c/number of batches); accumulate() yields total number of open channels in all batches
        Vhalf -= ChangeHalfVoltageDrift;  // decrese voltage at half maximum
    } else { // if voltage has exceeded the threshold for AP firing
        Voltage = VoltageReset; // reset the voltage to the membrane resting potential
        Vhalf += ChangeHalfVoltageReset;  // increase the voltage at half maximum by d+
        for (int i=0; i<VoltageUpdatesPerChannelOpening; i++){  // redraw the number of open channels per batch according to a binomial distribution with opening probability in the system reset state
            VectorOpenChannels[i]=BinomialDist(NumberChannelsPerBatch, OpenProb(Voltage, Vhalf, WidthSigmoidVoltage, VoltageShift));
        }
    }
    
}

void OneSimulationStepWithoutHalfVoltageFeedback (double &Voltage, double &Vhalf, double VoltageShift, vector<int> &VectorOpenChannels, int LabelBatch){  //same as OneSimulationStep, only that the voltage at half maximum is not updated (control parameter is not subject to feedback)
    if (Voltage < UpperThresholdVoltage) {
        VectorOpenChannels[LabelBatch]=BinomialDist(NumberChannelsPerBatch, OpenProb(Voltage, Vhalf, WidthSigmoidVoltage, VoltageShift));
        Voltage += (EquilibriumVoltage-Voltage)/VoltageTimeScale + MaxVoltageChange/VoltageTimeScale/NumberChannels*accumulate(VectorOpenChannels.begin(), VectorOpenChannels.end(), 0.);
    } else {
        Voltage = VoltageReset;
        for (int i=0; i<VoltageUpdatesPerChannelOpening; i++){
            VectorOpenChannels[i]=BinomialDist(NumberChannelsPerBatch, OpenProb(Voltage, Vhalf, WidthSigmoidVoltage, VoltageShift));
        }
    }
    
}


void PrintToFile(const string& filename, const vector<double>& VoltageSeries, const vector<double>& HalfVoltageSeries, const vector<double>& VoltageShiftTemperatureSeries, const vector<double>& TimeSeries)  { // print simulation data to file (e.g., for time series)
    if(!PrintFile)
        return;
    ofstream file(filename);
    if(file.is_open())  {
        cout << "opened file" << endl;
        file << "voltage: " << endl; // print recorded voltage over time
        for (int i=0; i<VoltageSeries.size(); ++i) {
            file << VoltageSeries[i] << " ";
        }
        file << endl;
        file << endl;
        file << "half voltage: " << endl; // print recorded voltage at half maximum over time
        for (int i=0; i<HalfVoltageSeries.size(); ++i) {
            file << HalfVoltageSeries[i] << " ";
        }
        file << endl;
        file << endl;
        file << "opening probability: " << endl; // print channel opening probability over time
        for (int i=0; i<HalfVoltageSeries.size(); ++i) {
            file << OpenProb(VoltageSeries[i], HalfVoltageSeries[i], WidthSigmoidVoltage, VoltageShiftTemperatureSeries[i]) << " ";
        }
        file << endl;
        file << endl;
        file << "times: " << endl; // print corresponding times
        for (int i=0; i<TimeSeries.size(); ++i) {
            file << TimeSeries[i] << " ";
        }
    }
    else {
        cout << "could not open file/n";
    }
    file.close();
}

void PrintToFileSweep(const string& filename, const vector<double>& SpikeTimesBefore, const vector<double>& SpikeTimesAfter)  { // print simulation data to file (e.g., for calculating spike time distribution)
    if(!PrintFile)
        return;
    ofstream file;
    file.open(filename, std::ios_base::app);
    if(file.is_open())  {
        cout << "opened file" << endl; // print vector of spike times before the temperature is varied
        for (int i=0; i<SpikeTimesBefore.size(); ++i) {
            file << fixed << int(SpikeTimesBefore[i]) << " ";
        }
        file << endl; // print vector of spike times after the temperature is varied
        for (int i=0; i<SpikeTimesAfter.size(); ++i) {
            file << fixed << int(SpikeTimesAfter[i]) << " ";
        }
        file << endl;
        file << endl;
    }
    else {
        cout << "could not open file/n";
    }
    file.close();
}


void PrintToFileBetweenSweeps(const string& filename)  { // between sweeps, add empty line to the file
    if(!PrintFile)
        return;
    ofstream file;
    file.open(filename, std::ios_base::app);
    if(file.is_open())  {
        file << endl;
    }
    else {
        cout << "could not open file/n";
    }
    file.close();
}

void PrintToFileSweepPrepare(const string& filename, const string& sweepparameter, const vector<double>& Parameters)  { // at the beginning, add description of current sweep parameter and values of the sweep parameter to the file
    if(!PrintFile)
        return;
    ofstream file;
    file.open(filename, std::ios_base::app);
    if(file.is_open())  {
        cout << "opened file" << endl;
        cout << "sweep parameter is " << sweepparameter << endl;  // console output to indicate sweep parameter
        file << "sweep parameter is " << sweepparameter << endl;  // indicate which parameter sweep is done
        for (int i=0; i<Parameters.size(); ++i) { // print sweep parameter values
            file << Parameters[i] << " ";
        }
        file << endl;
        file << endl;
    }
    else {
        cout << "could not open file/n";
    }
    file.close();
}


void PrepareFileEnsemble(const string& filename)  { // indicate in file that spike times will be recorded
    if(!PrintFile)
        return;
    ofstream file;
    file.open(filename, std::ios_base::app);
    if(file.is_open())  {
        cout << "prepared file" << endl;
        file << "Record spike times" << endl;
    }
    else {
        cout << "could not open file/n";
    }
    file.close();
}

void ResetSystem (double& Volt, double& HalfVolt, const double& HalfVoltStartValue, double& VoltShiftT){  // reset system after firing of AP
    Volt = EquilibriumVoltage;  // reset voltage back to membrane resting potential
    HalfVolt = HalfVoltStartValue;  // set voltage at half maximum to start value
    VoltShiftT=0.;  // set shift in voltage at half maximum to zero
}


void MakeSingleTimeSeries (double& HalfVolStartVal){  // run the dynamics to generate a time series of the system (system is subject to feedback in the control parameter before the temperature is changed)
    vector<double> VoltageOverTime;  // vector to record voltage over time
    VoltageOverTime.reserve(NumberRecordings+1);  // reserve size according to the expected number of recordings
    vector<double> HalfVoltageOverTime;  // vector to record voltage at half maximum over time
    HalfVoltageOverTime.reserve(NumberRecordings+1);  // reserve size according to the expected number of recordings
    vector<double> VoltageShiftTemperatureOverTime;  // vector to record external shift in voltage at half maximum over time (e.g., due to change in temperature)
    VoltageShiftTemperatureOverTime.reserve(NumberRecordings+1);  // reserve size according to the expected number of recordings
    vector<double> Time;  // vector to record time points
    Time.reserve(NumberRecordings+1);  // reserve size according to the expected number of recordings
    double Voltage = EquilibriumVoltage;  // set voltage to membrane resting potential
    double HalfVoltage = HalfVolStartVal;  // set voltage at half maximum to starting value
    for (int i=0; i<SimulationStepsStartRecording; ++i) {  // perform simulation steps with feedback before recording; the number of channels in batch i%VoltageUpdatesPerChannelOpening is updated (that is, batches are updated one by one in a deterministic, periodically repeating pattern)
        OneSimulationStep(Voltage, HalfVoltage, VoltageShiftTemperature, ListNumberOpenChannelsBatches, i%VoltageUpdatesPerChannelOpening);
    }
    int TrackSimulationStep=0;  // start to track simulation steps after initial unrecorded simulation phase
    if constexpr (HalfVoltageFeebackAfterVoltageShift) {  // if the control parameter is subject to feedback after the temperature has been changed
        for (int i=SimulationStepsStartRecording; i<SimulationSteps; ++i) { // simulate system until total number of simulation steps has been reached
            OneSimulationStep(Voltage, HalfVoltage,VoltageShiftTemperature, ListNumberOpenChannelsBatches, i%VoltageUpdatesPerChannelOpening); // simulate one timestep of the system
            if (TrackSimulationStep%SimulationStepsBetweenRecording==0) {  // record state of the system every SimulationStepsBetweenRecording-th step
                VoltageOverTime.push_back(Voltage);  // record voltage
                HalfVoltageOverTime.push_back(HalfVoltage); // record voltage at half maximum
                VoltageShiftTemperatureOverTime.push_back(VoltageShiftTemperature); // record external shift of voltage
                Time.push_back(TrackSimulationStep);  // record time (in units of tau_c/number of batches)
            }
            if (i==SimStepShiftTemperature) {  // in the middle of the simulation, change the voltage at half maximum/the temperature (if NewTemperature is nonzero)
                VoltageShiftTemperature=NewTemperature;
            }
            TrackSimulationStep+=1;  // increase value of simulation step
        }
    }
    else {  // if the control parameter is not subject to feedback after the temperature has been changed; otherwise analogous to previous block
        for (int i=SimulationStepsStartRecording; i<SimulationSteps; ++i) {
            if(i<=SimStepShiftTemperature){
                OneSimulationStep(Voltage, HalfVoltage,VoltageShiftTemperature,ListNumberOpenChannelsBatches, i%VoltageUpdatesPerChannelOpening);
            }
            else{
                OneSimulationStepWithoutHalfVoltageFeedback(Voltage, HalfVoltage,VoltageShiftTemperature,ListNumberOpenChannelsBatches, i%VoltageUpdatesPerChannelOpening);
            }
            if (TrackSimulationStep%SimulationStepsBetweenRecording==0) {
                VoltageOverTime.push_back(Voltage);
                HalfVoltageOverTime.push_back(HalfVoltage);
                VoltageShiftTemperatureOverTime.push_back(VoltageShiftTemperature);
                Time.push_back(TrackSimulationStep);
            }
            if (i==SimStepShiftTemperature) {
                VoltageShiftTemperature=NewTemperature;
            }
            TrackSimulationStep+=1;
        }
    }
    PrintToFile(FileNameSingleTimeSeries, VoltageOverTime, HalfVoltageOverTime, VoltageShiftTemperatureOverTime, Time);  // print to file
    ResetSystem (Voltage, HalfVoltage, HalfVolStartVal, VoltageShiftTemperature);  // reset system

}

void MakeSingleTimeSeriesWithoutFeedbackBefore (double& HalfVolStartVal){ // run the dynamics to generate a time series of the system (system is not subject to feedback in the control parameter before the temperature is changed); otherwise analogous to MakeSingleTimeSeries
    vector<double> VoltageOverTime;
    VoltageOverTime.reserve(NumberRecordings+1);
    vector<double> HalfVoltageOverTime;
    HalfVoltageOverTime.reserve(NumberRecordings+1);
    vector<double> VoltageShiftTemperatureOverTime;
    VoltageShiftTemperatureOverTime.reserve(NumberRecordings+1);
    vector<double> Time;
    Time.reserve(NumberRecordings+1);
    double Voltage = EquilibriumVoltage;
    double HalfVoltage = HalfVolStartVal;
    for (int i=0; i<SimulationStepsStartRecording; ++i) {
        OneSimulationStepWithoutHalfVoltageFeedback(Voltage, HalfVoltage, VoltageShiftTemperature,ListNumberOpenChannelsBatches, i%VoltageUpdatesPerChannelOpening);
    }
    int TrackSimulationStep=0;
    if constexpr (HalfVoltageFeebackAfterVoltageShift) {
        for (int i=SimulationStepsStartRecording; i<SimulationSteps; ++i) {
            if(i<=SimStepShiftTemperature){
                OneSimulationStepWithoutHalfVoltageFeedback(Voltage, HalfVoltage,VoltageShiftTemperature,ListNumberOpenChannelsBatches, i%VoltageUpdatesPerChannelOpening);
            }
            else{
                OneSimulationStep(Voltage, HalfVoltage,VoltageShiftTemperature,ListNumberOpenChannelsBatches, i%VoltageUpdatesPerChannelOpening);
            }
            if (TrackSimulationStep%SimulationStepsBetweenRecording==0) {
                VoltageOverTime.push_back(Voltage);
                HalfVoltageOverTime.push_back(HalfVoltage);
                VoltageShiftTemperatureOverTime.push_back(VoltageShiftTemperature);
                Time.push_back(TrackSimulationStep);
            }
            if (i==SimStepShiftTemperature) {
                VoltageShiftTemperature=NewTemperature;
            }
            TrackSimulationStep+=1;
        }
    }
    else {
        for (int i=SimulationStepsStartRecording; i<SimulationSteps; ++i) {

            OneSimulationStepWithoutHalfVoltageFeedback(Voltage, HalfVoltage,VoltageShiftTemperature,ListNumberOpenChannelsBatches, i%VoltageUpdatesPerChannelOpening);
            
            if (TrackSimulationStep%SimulationStepsBetweenRecording==0) {
                VoltageOverTime.push_back(Voltage);
                HalfVoltageOverTime.push_back(HalfVoltage);
                VoltageShiftTemperatureOverTime.push_back(VoltageShiftTemperature);
                Time.push_back(TrackSimulationStep);
            }
            if (i==SimStepShiftTemperature) {
                VoltageShiftTemperature=NewTemperature;
            }
            TrackSimulationStep+=1;
        }
    }
    PrintToFile(FileNameSingleTimeSeries, VoltageOverTime, HalfVoltageOverTime, VoltageShiftTemperatureOverTime, Time);
    ResetSystem (Voltage, HalfVoltage, HalfVolStartVal, VoltageShiftTemperature);

}




void MakeSweepFixedHalfVoltage(double& helpepsilon, string& FileNameSweepAlpha){  // record spike times for different values of alpha (in the absence of feedback): in order to change alpha, the voltage at half maximum is offset by HalfVoltageOffsetSweep[SweepNumber] for sweep number SweepNumber
    PrepareFileEnsemble(FileNameSweepAlpha);  // prepare file
    PrintToFileSweepPrepare(FileNameSweepAlpha, "Half Voltage Offset", HalfVoltageOffsetSweep);
    double Vhalfbifurcation =  EquilibriumVoltage+VoltageShiftTemperature+MaxVoltageChange*(0.5-0.5*sqrt(1.-4.*helpepsilon)-helpepsilon*log((1.-sqrt(1.-4.*helpepsilon)-2.*helpepsilon)/2./helpepsilon));  // calculate value of the voltage at half maximum at the bifurcation; helpepsilon is identical to variable rho in the paper
    cout << "Half Voltage value at bifurcation is " << Vhalfbifurcation << endl;
    for (int SweepNumber=0; SweepNumber<HalfVoltageOffsetSweep.size(); ++SweepNumber) {  // run simulation for all values of the sweep parameter
        cout << "Sweep Number is " << SweepNumber << endl;  // track sweep number in console
        int RoughNumberSpikes = (int) 1.25*(SimulationSteps-SimulationStepsStartRecording)/2./SweepParameters[SweepNumber]+1;  // calculate an approximate upper bound for the number of APs fired during the simulation
        vector<double> SpikeTimesBeforeShift;  // initialize vector of spike times before a potential shift in the temperature
        vector<double> SpikeTimesAfterShift;  // initialize vector of spike times after a potential shift in the temperature
        double Voltage = EquilibriumVoltage;  // set voltage to membrane resting potential
        double HalfVoltageForThisSweep= Vhalfbifurcation + HalfVoltageOffsetSweep[SweepNumber];  // set half voltage to the value at the bifurcation offset by HalfVoltageOffsetSweep for the respective sweep
        double HalfVoltage = HalfVoltageForThisSweep;  // set half voltage to the value at the bifurcation offset by HalfVoltageOffsetSweep for the respective sweep
        cout << "Half voltage for this sweep is " << HalfVoltageForThisSweep << endl;  // print value of half voltage for this sweep
        for (int RunNumber=0; RunNumber<NumberRuns; ++RunNumber) {  // perform NumberRuns simulations
            SpikeTimesBeforeShift.reserve(RoughNumberSpikes);  // reserve space for the vector of spike times
            SpikeTimesAfterShift.reserve(RoughNumberSpikes);  // reserve space for the vector of spike times
            for (int i=0; i<SimulationStepsStartRecording; ++i) {  // simulate system before recording starts (no feedback); the number of channels in batch i%VoltageUpdatesPerChannelOpening is updated (that is, batches are updated one by one in a deterministic, periodically repeating pattern)
                OneSimulationStepWithoutHalfVoltageFeedback(Voltage, HalfVoltage, VoltageShiftTemperature,ListNumberOpenChannelsBatches, i%VoltageUpdatesPerChannelOpening);
            }
            int TrackSimulationStep=0;  // start to record
            for (int i=SimulationStepsStartRecording; i<SimulationSteps; ++i) {
                
                OneSimulationStepWithoutHalfVoltageFeedback(Voltage, HalfVoltage, VoltageShiftTemperature,ListNumberOpenChannelsBatches, i%VoltageUpdatesPerChannelOpening); // one simulation step (no feedback); the number of channels in batch i%VoltageUpdatesPerChannelOpening is updated (that is, batches are updated one by one in a deterministic, periodically repeating pattern)
                
                if (Voltage>=UpperThresholdVoltage&&i<SimStepShiftTemperature) {  // before temperature is shifted: if voltage has just passed the threshold for AP firing, record the time (i.e., record spike time); in units of simulation steps
                    SpikeTimesBeforeShift.push_back(TrackSimulationStep);
                }
                if (Voltage>=UpperThresholdVoltage&&i>SimStepShiftTemperature) { // after temperature is shifted: if voltage has just passed the threshold for AP firing, record the time (i.e., record spike time); in units of simulation steps
                    SpikeTimesAfterShift.push_back(TrackSimulationStep);
                }
                if (i==SimStepShiftTemperature) {  // in the middle of the simulation, update the temperature
                    VoltageShiftTemperature=NewTemperature;
                }
                TrackSimulationStep+=1;  // update simulation step
            }
            PrintToFileSweep(FileNameSweepAlpha, SpikeTimesBeforeShift, SpikeTimesAfterShift);  // print spike times to file
            SpikeTimesBeforeShift.clear();  // for case when simulation is run several times, clear spike times
            SpikeTimesAfterShift.clear();
            ResetSystem (Voltage, HalfVoltage, HalfVoltageForThisSweep, VoltageShiftTemperature);  // for case when simulation is run several times: reset system
        }
        PrintToFileBetweenSweeps(FileNameSweepAlpha);
    }
}

int main() {
    double epsilon=WidthSigmoidVoltage/MaxVoltageChange;  // epsilon is called rho in the paper (DeltaV c_mem)/(i_c tau_rest)
    cout << "constant epsilon is " << epsilon << endl;  // print value of epsilon
    cout << "Upper threshold value is " << UpperThresholdVoltage << endl;  // print value of threshold for AP firing
    if constexpr (SingleTimeSeriesSimulation){  // if single time series should be recorded
        double HalfVoltageStartValue = EquilibriumVoltage+VoltageShiftTemperature+MaxVoltageChange*(0.5-0.5*sqrt(1.-4.*epsilon)-epsilon*log((1.-sqrt(1.-4.*epsilon)-2.*epsilon)/2./epsilon)) + HalfVoltageOffset;  // set voltage at half maximum
        cout << "Half Voltage steady state at bifurcation + offset is " << HalfVoltageStartValue << endl;
        if constexpr(HalfVoltageFeebackBeforeVoltageShift){ // if voltage at half maximum should be subject to feedback before voltage (temperature) shift
            MakeSingleTimeSeries (HalfVoltageStartValue);  // generate time series
        }
        else{ // if voltage at half maximum should not be subject to feedback before voltage (temperature) shift
            MakeSingleTimeSeriesWithoutFeedbackBefore (HalfVoltageStartValue);  // generate time series
        }
    }
   if constexpr (FixedHalfVoltageSweep) {  // if sweep for spike times should be recorded
        string FileNameHelp1 = "alpha_";  // used to make file name for single time series
        for (int i=0; i<HalfVoltageOffsetSweep.size(); ++i) {
            FileNameHelp1 = FileNameHelp1 + to_string(double(HalfVoltageOffsetSweep[i])) + "_";
        }
        string FileNameHelp2 = "_N_" + to_string(int(NumberChannels));
        string FileNameHelpSweep = "_SimSteps2e"+to_string(int(log2(SimulationSteps)))+"_ChangeHalfVol"+ to_string(ChangeHalfVoltageReset)+"_tauve"+to_string(int(log10(VoltageTimeScale)))+"_tauc"+to_string(ChannelOpeningTime)+"_VoltUpdPerChannel"+to_string(int(VoltageUpdatesPerChannelOpening)) +"_FeedbackAfter"+to_string(int(HalfVoltageFeebackAfterVoltageShift));
        string FileNameSweepFixed = "sweep_" + FileNameHelp1 + "DT_" + to_string(NewTemperature)+ FileNameHelp2 + FileNameHelpSweep+"_sample"+to_string(NumberRuns) +".txt";
            
        MakeSweepFixedHalfVoltage(epsilon, FileNameSweepFixed);  // make sweep
    }

    

    
    return 0;
}



