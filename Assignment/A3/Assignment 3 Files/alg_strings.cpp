/******************************************************************************
 *  File: alg_strings.cpp
 * 
 *  An implementation of string-matching classes.
 *
 ******************************************************************************/

#include <iostream>
#include <stack>
#include <list>
#include <string>
#include <random>
#include <chrono>
#include "alg_strings.h"

/******************************************************************************
 *  Class: RabinKarp
 *  A class implementing the Rabin-Karp algorithm
 ******************************************************************************/
long RabinKarp::long_random_prime(){
  int min = 0, max = 49;
  long primes[] = {
    2094665479L, 1783990163L, 2094521287L, 2134397081L, 2126326253L, 
    1957216747L, 1436547389L, 1428780767L, 2075625529L, 1593123733L, 
    2132587157L, 1965562429L, 1164701777L, 1568991883L, 2130061793L, 
    1075370311L, 1711832929L, 2054631589L, 1587361861L, 1435348609L, 
    1332084959L, 1465215911L, 2088173753L, 1933073123L, 1319415599L, 
    1211741129L, 1487473783L, 1656920599L, 1817614213L, 1838911937L, 
    1697951429L, 1673793083L, 1971101663L, 1570547117L, 1869368041L, 
    1855484017L, 2057695543L, 1806695647L, 2082498797L, 2090345119L, 
    1349212999L, 1456810283L, 1271362889L, 1959057733L, 1073964823L, 
    1315871351L, 1308843649L, 1543027127L, 1230659387L, 1828780297L };
  
  static std::default_random_engine en;
  en.seed(std::chrono::system_clock::now().time_since_epoch().count());
  static std::uniform_int_distribution<int> dist{min, max}; 
    
  return primes[dist(en)];
}
long RabinKarp::hash(const std::string &key, int m) const {
  long h = 0;
  for (int j = 0; j < m; j++) {
    h = (R * h + key[j]) % q;
  }
  return h;
}

bool RabinKarp::check(const std::string &txt, int i) const {
  for (int j = 0; j < m; j++) {
    if (pat[j] != txt[i + j]) {
      return false;
    }
  }

  return true;
}

RabinKarp::RabinKarp(const std::string &pat):pat(pat), R(256) {
  m = pat.size();
  q = long_random_prime();

  RM = 1;
  for (int i = 1; i <= m-1; i++){
    RM = (R * RM) % q;
  }
  pat_hash = hash(pat, m);
}

int RabinKarp::search(const std::string &txt) const {
  int n = txt.length();
  if (n < m) return n;
  long txt_hash = hash(txt, m);

  // check for match at offset 0
  if ((pat_hash == txt_hash) && check(txt, 0))
    return 0;

  std::cout << "txt_hash: " << txt_hash << std::endl;
  // check for hash match; if hash match, check for exact match
  for (int i = m; i < n; i++) {
    // Remove leading digit, add trailing digit, check for match.
    txt_hash = (txt_hash + q - RM * txt[i - m] % q) % q;
    txt_hash = (txt_hash * R + txt[i]) % q;

	std::cout << "txt_hash: " << txt_hash << std::endl;
    // match
    int offset = i - m + 1;
    if ((pat_hash == txt_hash) && check(txt, offset)){
      return offset;
    }
  }

  return n; // no match
}													