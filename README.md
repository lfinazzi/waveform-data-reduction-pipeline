# waveform-data-reduction-pipeline
Pipeline software to reduce size of data acquired with picoscope 2406B for posterior analysis. This software was designed to go from 6 TB of waveform data acquired in coincidence to 600 GB of important waveform parameters for data analysis in a quantum optics experiment. Inputs are 2*N waveform events and outputs are 2*N set of important waveform parameters. This results in a data size reduction of 88%. 

REQUIREMENTS:

  1. ROOT library for corresponding platform (found in https://root.cern/install/)

INPUT FORMAT OF WAVEFORM DATA:

One file for both channels (waveforms) acquired in coincidence with picoscope 2406B. Each event is saved in an interleaved way. This means that CH1's first waveform sample (1B) is followed by CH2's first waveform sample (1B), followed by CH1's second sample (1B), followed by CH2's second sample (1B), and so on. Each event is saved following the previous one. The input file size in Bytes is 2 * NUM_SEGMENTS * SEGMENT_SIZE_IN_SAMPLES (input file size can be changed in functions.h).

PARAMETERS CALCULATED PER EVENT:

For each waveform, the following parameters are calculated:

  1. trigger time (float)
  2. waveform minimum (signed char)
  3. waveform maximum (signed char)
  4. Time-over-Threshold (float)
  5. baseline level (signed char)
  6. baseline standard deviation (signed char)
  7. energy (utin16_t)

These are saved in one output file and data is saved in an interleaved way (like the input waveform file), i.e. event1 CH1 trigger (4B), event1 CH2 trigger (4B), event1 CH1 waveform minimum (1B), event1 CH2 waveform minimum (1B), etc. 
