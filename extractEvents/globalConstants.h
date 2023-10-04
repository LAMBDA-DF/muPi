#ifndef _globalConstants_h_
#define _globalConstants_h_


extern int gHitMaxSize;
extern int gVerbosity;
extern const char *gBaseTNtupleVars[];
extern const int   gNBaseTNtupleVars;
extern const char *gExtraTNtupleVars[];
extern const int   gNExtraTNtupleVars;

inline double sq(double x){return x*x;};

extern const int kMaxLine;
extern const int kCompBL;
extern const int kCompBLWindow;
extern const int kCompRowMedia;
extern const int kCompNS;

extern const int kExtractedMask;
extern const int kEdgeFlag;
extern const int kSatFlag;
extern const int kBugFlag;
extern const int kMerFlag;
extern const double kSat;

extern const double kPrescan;
extern const float kSatValue;
extern const float kSatMargin;
extern const float kBugEnergy;

//colors
extern const char cyan[];
extern const char magenta[];
extern const char red[];
extern const char green[];
extern const char yellow[];
extern const char blue[];

extern const char bold[];

extern const char whiteOnRed[];
extern const char whiteOnGreen[];
extern const char whiteOnPurple[];
extern const char whiteOnViolet[];
extern const char whiteOnBrown[];
extern const char whiteOnGray[];

extern const char normal[];

#endif
