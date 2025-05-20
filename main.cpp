// data_analysis.cpp : Defines the entry point for the application.
//

#include "main.h"
#include "functions.h"
#include "TApplication.h"
#include "TF1.h"

int main(int argc, char **argv)
{

	TApplication *myApp = new TApplication("myApp", &argc, argv, 0, -1);
	gStyle->SetOptStat(111111);

	std::vector<int> params = { 30,					// trigger relative to baseline
							    8,					// baseline samples for baseline calculation
	};


	std::vector<std::string> paths = {"./raw_data/" 			// position of raw data acquired with picoscope 2406B
	};
	
	std::vector<std::string> save_paths = {"./reductions/",		// directory for reduced data
	};

	std::vector<int> numFiles = {150000 						// number of files in raw_data folder

	}; 	

	std::vector<bool> useRun = {true 							// use this run
	};

	
	int numRuns = 3;		// number of runs. Same as useRun.size()
	int start_file = 0;
	bool debugFlag = false;
	bool isPreAnalyzed = false;

	for(int i = 0; i < numRuns; i++){
		if(!isPreAnalyzed){
			RunPartialDataAnalysis(paths[i], save_paths[i], params, start_file, numFiles[i], debugFlag, useRun[i], i);
		}
	}
	
	//myApp->Run();

	return 0;
}
