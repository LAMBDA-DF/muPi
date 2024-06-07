#ifndef _globalConstants_h_
#define _globalConstants_h_

#include <math.h> 
#include <map>

int gVerbosity = true;

inline double sq(double x){return x*x;};

const int    kMaxLine = 10000;

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


/*
const char black[] = "\x1b[0;30m";
const char blue[] = "\x1b[0;34m";
const char green[] = "\x1b[0;32m";
const char cyan[] = "\x1b[0;36m";

const char red[] = "\x1b[0;31m";
const char purple[] = "\x1b[0;35m";
const char brown[] = "\x1b[0;33m";
const char gray[] = "\x1b[0;37m";


    
const char darkGray[] = "\x1b[1;30m";
const char lightBlue[] = "\x1b[1;34m";
const char lightGreen[] = "\x1b[1;32m";
const char lightCyan[] = "\x1b[1;36m";
const char lightRed[] = "\x1b[1;31m";
const char lightPurple[] = "\x1b[1;35m";
const char yellow[] = "\x1b[1;33m";
const char white[] = "\x1b[1;30m";
*/

#endif
