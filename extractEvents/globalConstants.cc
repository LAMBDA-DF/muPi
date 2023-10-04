#include "globalConstants.h"

int gHitMaxSize = 5000; //default hit max size
int gVerbosity = true;
const char *gBaseTNtupleVars[]  = {"runID", "ohdu", "nSat", "flag", "xMin", "xMax", "yMin", "yMax"};
const int   gNBaseTNtupleVars   = 8;
const char *gExtraTNtupleVars[] = {"E", "n", "xBary", "yBary", "xVar", "yVar"};
const int   gNExtraTNtupleVars  = 6;

const int kMaxLine = 10000;
const int kCompBL       = 1;
const int kCompBLWindow = 4;
const int kCompRowMedia = 8;
const int kCompNS       = 2;


const int kExtractedMask = -100000;
const int kEdgeFlag = 2;
const int kSatFlag = 4;
const int kBugFlag = 8;
const int kMerFlag = 16;
const double kSat  = 5e9;

const double kPrescan = 6;
const float kSatValue = 1e10;
const float kSatMargin = 0.9;
const float kBugEnergy = -50;

//colors
const char cyan[] = { 0x1b, '[', '1', ';', '3', '6', 'm', 0 };
const char magenta[] = { 0x1b, '[', '1', ';', '3', '5', 'm', 0 };
const char red[] = { 0x1b, '[', '1', ';', '3', '1', 'm', 0 };
const char green[] = { 0x1b, '[', '1', ';', '3', '2', 'm', 0 };
const char yellow[] = { 0x1b, '[', '1', ';', '3', '3', 'm', 0 };
const char blue[] = "\x1b[1;34m";

const char bold[] = "\x1b[1;39m";

const char whiteOnRed[]    = "\x1b[1;41m";
const char whiteOnGreen[]  = "\x1b[1;42m";
const char whiteOnPurple[] = "\x1b[1;45m";
const char whiteOnViolet[] = "\x1b[1;44m";
const char whiteOnBrown[]  = "\x1b[1;43m";
const char whiteOnGray[]   = "\x1b[1;47m";

const char normal[] = { 0x1b, '[', '0', ';', '3', '9', 'm', 0 };
