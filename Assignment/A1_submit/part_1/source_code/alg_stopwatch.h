/******************************************************************************
 *  File: alg_stopwatch.h
 * 
 *  A header file of a utility class for measuring the running time of an 
 *  algorithm. 
 * 
 *  Last modified by: Brian Rague
 *  Last modified on: Jan 1, 2023
 ******************************************************************************/

#ifndef _ADV_ALG_STOP_WATCH_H_
#define _ADV_ALG_STOP_WATCH_H_

#include <chrono>

class StopWatch {
private:
  std::chrono::high_resolution_clock::time_point start;

public:
  StopWatch();
  double elapsed_time();
  void reset();
};

#endif