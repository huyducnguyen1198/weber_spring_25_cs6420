/******************************************************************************
 *  File: alg_stopwatch.cpp
 * 
 *  An implementation of a utility class for measuring the running time of an 
 *  algorithm. 
 * 
 *  Last modified by: Brian Rague
 *  Last modified on: Jan 1, 2023
 ******************************************************************************/

#include "alg_stopwatch.h"

StopWatch::StopWatch(): start(std::chrono::high_resolution_clock::now()) {}

// Returns result in nanoseconds
double StopWatch::elapsed_time(){
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::nano> running_time = end - start;
  return running_time.count();
}

void StopWatch::reset(){
  start = std::chrono::high_resolution_clock::now();
}