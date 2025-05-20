#include "functions.h"

void SaveSegmentData(std::vector<Segment> &segments, std::string save_path)
{
    FILE* pFile;
    pFile = fopen(save_path.c_str(), "wb");
    for (unsigned int i = 0; i < NUM_SEGMENTS; i++)
    {
        /* saved in little endian, same as Intel processor */
        // Writes 14B per segment acquired 
        fwrite(&segments[i].trigs, sizeof(float), 2, pFile);                // 8B
        fwrite(&segments[i].v_min, sizeof(signed char), 2, pFile);          // 2B
        fwrite(&segments[i].v_max, sizeof(signed char), 2, pFile);          // 2B
        fwrite(&segments[i].ToT, sizeof(float), 2, pFile);                  // 4B
        fwrite(&segments[i].baseline, sizeof(signed char), 2, pFile);       // 2B
        fwrite(&segments[i].baseline_std, sizeof(signed char), 2, pFile);   // 2B
        fwrite(&segments[i].energy, sizeof(uint16_t), 2, pFile);            // 4B
                                                                            // total: 24B (200B original) --> Reduction: 88%
    }
    fclose(pFile);
    return;
}

std::vector<uint8_t> SegmentTime(uint8_t N, uint8_t step)
{
    std::vector<uint8_t> v(N);
    for (unsigned int i = 0; i < N; i++)
    {
        v[i] = step * i;
    }
    return v;
}

void LoadBinaryFile(std::string path, std::pair<std::vector<signed char>, std::vector<signed char>>& buffer)
{
    FILE* pFile = fopen(path.c_str(), "rb");

    if (pFile) {
        for (int i = 0; i < TOTAL_REPORT_SIZE; i++) {
            fread(&buffer.first[i], sizeof(signed char), 1, pFile);
            fread(&buffer.second[i], sizeof(signed char), 1, pFile);
        }
    }

    fclose(pFile);

    return;
}

void lowPassButterworthFilter(std::vector<signed char> &indata, std::vector<double> &butterworthFilteredData, double CutOff, double SamplingRate) {
    const int arrSize = PROCESSING_WINDOW;
    double Dat2[arrSize + 4] = {}; // Array with 4 extra points front and back

    // Copy indata to Dat2
    for (long r = 0; r < arrSize; r++) {
        Dat2[2 + r] = indata[r];
    }
    Dat2[0] = indata[0];
    Dat2[1] = indata[0];
    Dat2[arrSize + 2] = indata[arrSize - 1];
    Dat2[arrSize + 3] = indata[arrSize - 1];
    double wc = tan(CutOff * M_PI / SamplingRate);
    double k1 = 1.414213562 * wc; // Sqrt(2) * wc
    double k2 = wc * wc;
    double a = k2 / (1 + k1 + k2);
    double b = 2 * a;
    double c = a;
    double k3 = b / k2;
    double d = -2 * a + k3;
    double e = 1 - (2 * a) - k3;

    // RECURSIVE TRIGGERS - ENABLE filter is performed (first, last points constant)
    double DatYt[arrSize + 4] = {};
    DatYt[0] = indata[0];
    DatYt[1] = indata[0];
    for (long s = 2; s < arrSize + 2; s++) {
        DatYt[s] = a * Dat2[s] + b * Dat2[s - 1] + c * Dat2[s - 2] + d * DatYt[s - 1] + e * DatYt[s - 2];
    }
    DatYt[arrSize + 2] = DatYt[arrSize + 1];
    DatYt[arrSize + 3] = DatYt[arrSize + 1];

    // FORWARD filter
    double DatZt[arrSize + 2] = {};
    DatZt[arrSize] = DatYt[arrSize + 2];
    DatZt[arrSize + 1] = DatYt[arrSize + 3];
    for (long t = -arrSize + 1; t <= 0; t++) {
        DatZt[-t] = a * DatYt[-t + 2] + b * DatYt[-t + 3] + c * DatYt[-t + 4] + d * DatZt[-t + 1] + e * DatZt[-t + 2];
    }

    // Calculated points copied for return
    for (long p = 0; p < arrSize; p++) {
        butterworthFilteredData[p] = DatZt[p];
    }
    return;
}

std::vector<Segment> SliceWaveforms(std::pair<std::vector<signed char>, std::vector<signed char>> &raw_data, std::vector<Segment> &segments)
{
    std::vector<int16_t> cut_data0(raw_data.first.begin(), raw_data.first.end());
    std::vector<int16_t> cut_data1(raw_data.second.begin(), raw_data.second.end());

    int segment_size = SEGMENT_SIZE_IN_SAMPLES;
    int preSamples = (int)(PROCESSING_WINDOW / 2);
    int postSamples = preSamples;

    for (int i = 0; i < NUM_SEGMENTS; i++)
    {
        segments[i].samples0.resize(preSamples + postSamples);
        segments[i].samples1.resize(preSamples + postSamples);
        //for (int j = (int)(SEGMENT_SIZE_IN_SAMPLES/2) - preSamples; j < (int)(SEGMENT_SIZE_IN_SAMPLES / 2) + postSamples; j++)
        for (int j = 0; j < preSamples + postSamples; j++)
        {
            segments[i].samples0[j] = cut_data0[i * segment_size + j + ((int)(segment_size / 2) - preSamples)];
            segments[i].samples1[j] = cut_data1[i * segment_size + j + ((int)(segment_size / 2) - preSamples)];
        }
    }   

    return segments;
}

uint16_t GetEnergy(std::vector<signed char> &wf)
{
    uint16_t energy = 0;
    for(uint i = 0; i < wf.size(); i++)
        energy += (uint16_t)wf[i] + 100;        // we sum 100 so that wf is always positive and a higher pulse count results in higher energy
    return energy;
}

float mean(std::vector<signed char> &v)
{
    float mean = 0;
    for (int i = 0; i < v.size(); i++)
        mean += (float)v[i];
    return mean / v.size();
}

signed char stdev(std::vector<signed char> &v)
{
    float mean_v = mean(v);
    float sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    float stdev = std::sqrt(sq_sum / v.size() - mean_v * mean_v);
    return (signed char)round(stdev);
}

short int sum(std::vector<signed char> &v)
{
    short int sum = 0;
    for (int i = 0; i < v.size(); i++)
        sum += abs(v[i]);
    return (short int)sum;
}

signed char max(std::vector<signed char>& v)
{
    signed char max = -128;
    for (int i = 0; i < v.size(); i++) {
        if (v[i] > max)
            max = v[i];
    }
    return (signed char)max;
}

signed char min(std::vector<signed char>& v)
{
    std::vector<signed char>::iterator result = std::min_element(v.begin(), v.end());
    return *result;
}

uint8_t FindTriggerRough(std::vector<uint8_t> &time, std::vector<signed char> &signal, float trigger, std::string option)
{
    if(option == "up")
    {
        for (unsigned int i = 1; i < signal.size(); i++) {
            if (signal[i] >= trigger)
                return i - 1;
        }
    }
    if(option == "down")
    {
        int count = 0;
        for (unsigned int i = signal.size() - 2; i > 0; i--) {
            if (signal[i] >= trigger)
                return i + 1;
        }
    }
    return 0;
}

uint8_t FindTriggerRough(std::vector<uint8_t> &time, std::vector<double> &signal, float trigger, std::string option)
{
    if(option == "up")
    {
        for (unsigned int i = 1; i < signal.size(); i++) {
            if (signal[i] >= trigger)
                return i - 1;
        }
    }
    if(option == "down")
    {
        int count = 0;
        for (unsigned int i = signal.size() - 2; i > 0; i--) {
            if (signal[i] >= trigger)
                return i + 1;
        }
    }
    return 0;
}

float FindTriggerFine(std::vector<uint8_t> &time, std::vector<signed char> &signal, double trigger, std::string option)
{

    std::vector<double> filtSignal(signal.size());
    lowPassButterworthFilter(signal, filtSignal, 60E6, SAMPLING_RATE);

    uint8_t trIndexRough = FindTriggerRough(time, filtSignal, trigger, option);
    if(trIndexRough == 0)
        return 0;
    
    uint8_t x1 = time[trIndexRough];
    uint8_t x2 = time[trIndexRough + 1];
    double y1 = filtSignal[trIndexRough];
    double y2 = filtSignal[trIndexRough + 1];

    double m = (double)(y2 - y1) / (double)(x2 - x1);
    double triggerFine = x1 + (trigger - y1) / m;

    //std::cout << std::fixed << std::setprecision(2) << "[" << static_cast<int>(y1) << ", " << static_cast<int>(y2) << ", " << static_cast<int>(y2 - y1) << ", " << std::to_string(m) << ", " << triggerFine << ", ";
    //std::cout << m << std::endl << ",";
    return (float)triggerFine;
    
}

// LS ANALYTICAL CALCULATION
float FindTriggerFineFit(std::vector<uint8_t> &time, std::vector<signed char> &signal, float trigger, std::string option)
{
    std::vector<double> filtSignal(signal.size());
    lowPassButterworthFilter(signal, filtSignal, 60E6, SAMPLING_RATE);

    uint8_t trIndexRough = FindTriggerRough(time, filtSignal, trigger, option);
    if(trIndexRough == 0)
        return 0;

    double sum_x = 0;
    double sum_x2 = 0;
    double sum_y = 0;
    double sum_xy = 0;
    int numPoints = 2;
    double triggerFine = 0;

    if(option == "up")
    {
        for(uint8_t i = 0; i < numPoints; i++){
            sum_x += time[trIndexRough + i];
            sum_x2 += pow(time[trIndexRough + i], 2);
            sum_y += filtSignal[trIndexRough + i];
            sum_xy += filtSignal[trIndexRough + i] * time[trIndexRough + i];
        }

        double m = ( numPoints * sum_xy - sum_x * sum_y ) / ( numPoints * sum_x2 - pow(sum_x, 2) );  
        double b = ( sum_y - m * sum_x ) / numPoints;

        triggerFine = (trigger - b) / m;
    }

    else if(option == "down")
    {
        for(uint8_t i = 0; i < numPoints; i--){
            sum_x += time[trIndexRough + i];
            sum_x2 += pow(time[trIndexRough + i], 2);
            sum_y += filtSignal[trIndexRough + i];
            sum_xy += filtSignal[trIndexRough + i] * time[trIndexRough + i];
        }

        double m = ( numPoints * sum_xy - sum_x * sum_y ) / ( numPoints * sum_x2 - pow(sum_x, 2) );  
        double b = ( sum_y - m * sum_x ) / numPoints;

        triggerFine = (trigger - b) / m;
    }


    //std::cout << m << std::endl << ",";
    return (float)triggerFine;
}

float FindTriggerFineCubic(std::vector<uint8_t> &time, std::vector<signed char> &signal, float trigger)
{
    std::vector<double> filtSignal(signal.size());
    lowPassButterworthFilter(signal, filtSignal, 60E6, SAMPLING_RATE);

    uint8_t trIndexRough = FindTriggerRough(time, filtSignal, trigger, "up");
    int numPoints = 10;
    if(trIndexRough == 0 || (trIndexRough + numPoints) > PROCESSING_WINDOW)
        return 0;


    std::vector<double> time_int(numPoints);
    std::vector<double> signal_int(numPoints);

    for(uint8_t i = 0; i < numPoints; i++){
        time_int[i] = time[trIndexRough + i];
        //std::cout << static_cast<int>(time[trIndexRough + i]) << " \n";
        signal_int[i] = filtSignal[trIndexRough + i];
    }

    tk::spline s(time_int, signal_int, tk::spline::cspline_hermite);

    std::vector<double> sols = s.solve(trigger);
    if(sols.size() == 0)
        return 0;
    else
        return (float)sols[0];
}

int AnalyzeSegments(std::vector<Segment>& segments, std::vector<uint8_t>& seg_times, std::vector<int> params, std::string save_path, bool debugFlag)
{
    double v_trigger = params[0];
    int baselineSamples = params[1];
    TRandom r(0);

    std::vector<signed char> baseline_samples0(baselineSamples);
    std::vector<signed char> baseline_samples1(baselineSamples);

    float trig0 = 0;
    float trig1 = 0;

    float trig0d = 0;
    float trig1d = 0;

    for(int i = 0; i < NUM_SEGMENTS; i++)
    {
        // calculate baseline (and baseline std) of both channels and save it in Segment structure
        for(int j = 0; j < baselineSamples; j++){
            baseline_samples0[j] = segments[i].samples0[j];
            baseline_samples1[j] = segments[i].samples1[j];
        }
        segments[i].baselineFine[0] = mean(baseline_samples0);
        segments[i].baselineFine[1] = mean(baseline_samples1);
        segments[i].baseline[0] = (signed char)segments[i].baselineFine[0];
        segments[i].baseline[1] = (signed char)segments[i].baselineFine[1];
        segments[i].baseline_std[0] = stdev(baseline_samples0);
        segments[i].baseline_std[1] = stdev(baseline_samples1);

        //double timeNoise = 2000E-12;

        // calculate relative trigger of both channels and save it in Segment structure
        trig0 = FindTriggerFine(seg_times, segments[i].samples0, segments[i].baselineFine[0] + v_trigger, "up");
        trig1 = FindTriggerFine(seg_times, segments[i].samples1, segments[i].baselineFine[1] + v_trigger, "up");

        trig0d = FindTriggerFine(seg_times, segments[i].samples0, segments[i].baselineFine[0] + v_trigger, "down");
        trig1d = FindTriggerFine(seg_times, segments[i].samples1, segments[i].baselineFine[1] + v_trigger, "down");
        
        //std::cout << std::fixed << std::setprecision(2) << segments[i].baselineFine[0] << ", " << segments[i].baselineFine[0] + v_trigger << "],\n";
        
        //std::cout << segments[i].baselineFine[0] << ", ";
        segments[i].trigs[0] = trig0;
        segments[i].trigs[1] = trig1;

        // Calculate min/max of both channels and save it in Segment structure
        segments[i].v_max[0] = max(segments[i].samples0);
        segments[i].v_max[1] = max(segments[i].samples1);
        segments[i].v_min[0] = min(segments[i].samples0);
        segments[i].v_min[1] = min(segments[i].samples1);
        
        // Calculate pulse area of both channels and save it in Segment structure
        segments[i].ToT[0] = trig0d - trig0;
        segments[i].ToT[1] = trig1d - trig1;

        segments[i].energy[0] = GetEnergy(segments[i].samples0);
        segments[i].energy[1] = GetEnergy(segments[i].samples1);

        // for debug purposes
        if (debugFlag)
        {
            PlotWaveforms(segments[i], v_trigger);
            std::cout << "segment: " << std::to_string(i) << std::endl; 
            std::cout << "\ntrig0: " << segments[i].trigs[0] << std::endl;
            std::cout << "\ntrig1: " << segments[i].trigs[1] << std::endl;
            std::cout << "\ndt: " << segments[i].trigs[1] - segments[i].trigs[0] << std::endl;
            std::cout << "vmax0: " << static_cast<int>(segments[i].v_max[0]) << std::endl;
            std::cout << "vmax1: " << static_cast<int>(segments[i].v_max[1]) << std::endl;
            std::cout << "baseline0: " << static_cast<int>(segments[i].baselineFine[0]) << std::endl;
            std::cout << "baseline1: " << static_cast<int>(segments[i].baselineFine[1]) << std::endl;
            std::cout << "baseline_std0: " << static_cast<int>(segments[i].baseline_std[0]) << std::endl;
            std::cout << "baseline_std1: " << static_cast<int>(segments[i].baseline_std[1]) << std::endl;
            std::cout << "tot0: " << segments[i].ToT[0] << std::endl;
            std::cout << "tot1: " << segments[i].ToT[1] << std::endl;
            std::cout << "energy0: " << segments[i].energy[0] << std::endl;
            std::cout << "energy1: " << segments[i].energy[1] << std::endl;
            return -1;
        }

    }

    // save Segment vector to file
    SaveSegmentData(segments, save_path);

    return 0;
}

void RunPartialDataAnalysis(std::string path, std::string save_path, std::vector<int> params, int start_file, int num_files, bool debugFlag, bool useRun, int numOfRun)
{
    if(!useRun)
        return;

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "running data pre-analysis...\n";
    std::vector<Segment> segments(NUM_SEGMENTS);
    std::vector<uint8_t> seg_times;

    std::pair<std::vector<signed char>, std::vector<signed char>> raw_data;
    raw_data.first.resize(TOTAL_REPORT_SIZE);
    raw_data.second.resize(TOTAL_REPORT_SIZE);

    for (int i = 0; i < 1; i++)
    {
        for (unsigned int j = start_file; j < num_files; j++)
        {
            LoadBinaryFile(path + "rapidBlock_" + std::to_string(j) + ".bin", raw_data);

            // could append many raw_data to a greater raw_data vector to save more segments per file and
            // it should work with minor modifications below (just saving more segments)

            if (raw_data.first.size() == raw_data.second.size() && raw_data.first.size() == NUM_SEGMENTS * SEGMENT_SIZE_IN_SAMPLES)
            {
                seg_times = SegmentTime(PROCESSING_WINDOW, 2);
                SliceWaveforms(raw_data, segments);
                if(AnalyzeSegments(segments, seg_times, params, save_path + "rapidBlock_preanalyzed_" + std::to_string(j) + ".bin", debugFlag) == -1)
                    return;
            }
            else
            {
                std::cout << "Data sanity check failed. Both channels have different number of samples of report size is incorrect in settings!\n";
                exit(1);
            }
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<seconds>(stop - start);
            if (j % 100 == 0)
                std::cout << "pre-analysis of " << std::to_string(j) << " files of " << std::to_string(num_files) << " completed in " << duration.count() << " seconds." << std::endl;

        }
        std::cout << "finished pre-analysis data saving (run " + std::to_string(numOfRun) + ", " + std::to_string(num_files) + " files)\n";
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    std::cout << "\npre-analysis finished in " << duration.count() << " seconds.\n";
    return;
}

void PlotWaveforms(Segment &segment, float trigger)
{
    std::vector<float> samples0(segment.samples0.begin(), segment.samples0.end());
    std::vector<float> samples1(segment.samples1.begin(), segment.samples1.end());

    std::vector<uint8_t> time = SegmentTime(PROCESSING_WINDOW, 2);
    std::vector<float> time0(time.begin(), time.end());

    TCanvas *c = new TCanvas();
    TGraph *g1 = new TGraph(time0.size(), &time0[0], &samples0[0]);
    TGraph *g2 = new TGraph(time0.size(), &time0[0], &samples1[0]);
    TF1 *l1 = new TF1("l1", "[0]", time0[0], time0[SEGMENT_SIZE_IN_SAMPLES - 1]);
    TF1 *l2 = new TF1("l2", "[0]", time0[0], time0[SEGMENT_SIZE_IN_SAMPLES - 1]);

    l1->SetParameter(0, segment.baselineFine[0] + trigger);
    l2->SetParameter(0, segment.baselineFine[1] + trigger);
    l1->SetLineColor(kBlack);
    l2->SetLineColor(kRed);

    TMultiGraph *mg = new TMultiGraph();
    g1->SetMarkerStyle(20);
    g2->SetMarkerStyle(20);
    g1->SetMarkerSize(1);
    g2->SetMarkerSize(1);
    g2->SetMarkerColor(kRed);
    mg->Add(g1);
    mg->Add(g2);

    c->cd();
    mg->Draw("AP");
    l1->Draw("SAME");
    l2->Draw("SAME");

    c->Update();
    c->Modified();

    return;
}

