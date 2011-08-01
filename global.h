#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include <utility>
#include <limits>
#include <cstring>
using std::pair;
// important: each node has at most 4 connected nets
enum DIRECTION{WEST, EAST, SOUTH, NORTH, BOTTOM, TOP, UNDEFINED};
enum NET_TYPE{RESISTOR, CURRENT, VOLTAGE};
const int NUM_NET_TYPE = 3;
enum LAYER_DIR{HR, VT, NA}; // NA means not available
enum RUN_MODE{IT,LU};
enum CIRCUIT_TYPE{WB, C4, UNKNOWN};

const size_t INFTY = std::numeric_limits<size_t>::max();
const int MAX_BUF = 256;
const int MAX_LAYER = 30;
const double VIRTUAL_RESISTER=1e-10;
const double FLOAT_ZERO = 1e-30;
const int OUTPUT_WIDTH_FLOAT = 15;
const int OUTPUT_WIDTH_INDEX = 3;
const int OUTPUT_WIDTH_STRING = 10;
// overlap_ratio should < 0.45, or will be divergence
const size_t X_BLOCKS = 3; // # of blocks along x axis
const size_t Y_BLOCKS = 3; // # of blocks along y axis
static size_t x_max = 0;
static size_t y_max = 0;
static size_t x_min = 0;
static size_t y_min = 0;
float overlap_ratio = 0.2;
typedef pair<size_t,size_t> SizeTPair;
class Triplet;
typedef Triplet Matrix;

//typedef vector<Net *> NetPtrVector;
#endif
