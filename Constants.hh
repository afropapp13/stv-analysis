#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>

#include "TString.h"

namespace Constants {

	//----------------------------------------//

        constexpr double BIG_DOUBLE = 1e300;
        const std::string SAMPLE_NAME = "MicroBooNE_CC1MuNp_XSec_2D_PpCosp_nu_MC";
        const int LINE_THICKNESS = 2;
        const double TEXT_SIZE = 0.06;
        const int TEXT_FONT = 132;
        const int NDIVISIONS = 8;

        const TString PATH_XSEC = "/uboone/data/users/apapadop/xsec_analyzer/xsec/";
        const TString PATH_EVENTS = "/uboone/data/users/apapadop/xsec_analyzer/events/";
        const TString PATH_POSTPROCESS = "/uboone/data/users/apapadop/PeLEE_PostProcessing/";
        const TString WORK_AREA = "/uboone/app/users/apapadop/master_stv_analysis/";
      
	//----------------------------------------//		

}
#endif
