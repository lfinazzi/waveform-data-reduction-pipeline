#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#define NUM_SEGMENTS                                                                          150000
#define SEGMENT_SIZE_IN_SAMPLES                                                                  100
#define TOTAL_REPORT_SIZE                                     NUM_SEGMENTS * SEGMENT_SIZE_IN_SAMPLES
#define SAMPLING_RATE                                                                          500E6   
#define PROCESSING_WINDOW                                                                        100
#define COINCIDENCE_WINDOW                                                                        50

#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <vector> 
#include <numeric>
#include <stdint.h>
#include <chrono>
#include <cstdlib>
#include <math.h>
#include "string.h"
#include <algorithm>
#include "spline.h"

#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TFitResult.h"
#include "TRandom.h"

using namespace std::chrono;

struct Segment {
    // variables loaded with raw data
    std::vector<signed char> samples0;
    std::vector<signed char> samples1;

    // variables determined after segment analysis
    signed char v_max[2];
    signed char v_min[2];
    signed char baseline[2];
    float baselineFine[2];
    signed char baseline_std[2];
    float trigs[2];
    float ToT[2];
    uint16_t energy[2];
};

/******************************
 *  FIRST STREAM FUNCTIONS
 ******************************/

// saves downscaled data from the first stream
void SaveSegmentData(std::vector<Segment> &segments, std::string save_path);

// generates segment time x-axis for plots and interpolations
std::vector<uint8_t> SegmentTime(uint8_t N, uint8_t step);

// loads raw data
void LoadBinaryFile(std::string path, std::pair<std::vector<signed char>, std::vector<signed char>>& buffer);

// Butterworth Low-Pass Filter of fourth order //
void lowPassButterworthFilter(std::vector<signed char> &indata, std::vector<double> &butterworthFilteredData, double CutOff, double SamplingRate);

// loads raw data into Segment structures
std::vector<Segment> SliceWaveforms(std::pair<std::vector<signed char>, std::vector<signed char>>&raw_data, std::vector<Segment> &segments, std::vector<int> params);

uint16_t GetEnergy(std::vector<signed char> &wf);

// calculates mean of a vector
float mean(std::vector<double> &v);

// calculates st. dev. of a vector
signed char stdev(std::vector<double> &v);

// calculates the accumulated sum of elements of a vector
short int sum(std::vector<double> &v);

// calculates the max value of a vector
signed char max(std::vector<double>& v);

// calculates the min value of a vector
signed char min(std::vector<double>& v);

// Function finds rough trigger position of wf and returns its index
uint8_t FindTriggerRough(std::vector<uint8_t> &time, std::vector<signed char> &signal, float trigger, std::string option);

// Function finds rough trigger position of wf and returns its index
uint8_t FindTriggerRough(std::vector<uint8_t> &time, std::vector<double> &signal, float trigger, std::string option);

// Linearly interpolates X and Y data
void LinearInterpolation(std::vector <uint8_t> &time, std::vector<signed char>& signal, uint8_t interpolation, std::vector<float> &time_interp, std::vector<float> signal_iterp);

// Function finds the fine trigger (around the previously found rought trigger with a wf interpolation 
float FindTriggerFine(std::vector<uint8_t> &time, std::vector<signed char> &signal, double trigger, std::string option);

// Function finds the fine trigger (around the previously found rought trigger with a wf interpolation 
float FindTriggerFineFit(std::vector<uint8_t> &time, std::vector<signed char> &signal, float trigger);

// Function to find the fine trigger with a cubic interpolation of wf
float FindTriggerFineCubic(std::vector<uint8_t> &time, std::vector<signed char> &signal, float trigger);

// Runs data analysis on segment 
int AnalyzeSegments(std::vector<Segment>& segments, std::vector<uint8_t>& seg_times, std::vector<int> params, std::string save_path, bool debugFlag);

// Runs raw data analysis
void RunPartialDataAnalysis(std::string path, std::string save_path, std::vector<int> params, int start_file, int num_files, bool debugFlag, bool useRun, int numOfRun);

// plots both channels of segment data
void PlotWaveforms(Segment &segment, float trigger);

#endif
