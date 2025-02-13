/******************************************************************************
 *  File: alg_strings.h
 * 
 *  A header file defining string-matching classes. 
 *
 ******************************************************************************/

#ifndef _ADV_ALG_STRINGS_H_
#define _ADV_ALG_STRINGS_H_

#include <iostream>
#include <string>

/******************************************************************************
 *  Class: RabinKarp
 *  A class implementing the Rabin-Karp algorithm
 ******************************************************************************/
class RabinKarp {
private:
  std::string pat;
  long pat_hash;   // Pattern hash value
  int m;           // Pattern length
  long q;          // A large prime, small enough to avoid overflow
  int R;           // Radix
  long RM;         // R^(M-1) % Q

  long hash(const std::string &key, int m) const;
  bool check(const std::string &txt, int i) const;

  static long long_random_prime();

public:
  RabinKarp(const std::string &pat);
  int search(const std::string &txt) const;
};

#endif