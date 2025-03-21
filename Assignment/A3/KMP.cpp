#include <vector>
#include <iostream>
#include <string>

std::vector<int> compute_prefix_function(const std::string &pattern) {
  int m = pattern.size();

  std::vector<int> pi(m, 0);

  int k = 0;

  pi[0] = 0;

  for(int q = 1; q < m; q++) {
	while(k > 0 && pattern[k] != pattern[q]) {
	  	k = pi[k - 1]; // in psuedo code it is k = pi[k] because the array is 1-indexed
	}

	if (pattern[k] == pattern[q]) {
	  	k = k + 1;
	}

	pi[q] = k;
  }

  return pi;
}

void KMP(const std::string &text, const std::string &pattern) {
	int n = text.size();
	int m = pattern.size();

	std::vector<int> pi = compute_prefix_function(pattern);

	int q = 0;

	for(int i = 0; i < n; i++){
		while(q > 0 && pattern[q] != text[i]){
			q = pi[q - 1];
		}

		if(pattern[q] == text[i]){
			q = q + 1;
		}

		if(q == m){
			std::cout << "Pattern occurs with shift " << i - m + 1 << std::endl;
			q = pi[q - 1];
		}
	}
}

int main(int argc, char const *argv[]){
	std::string pattern = "";

	// get pattern from argv
	if (argc == 2) {
		pattern = argv[1];
	}else{
		std::cout << "Usage: " << argv[0] << " <pattern>" << std::endl;
		return 1;
	}

	std::vector<int> pi = compute_prefix_function(pattern);

	for(int i = 0; i < pi.size(); i++){
		std::cout << pi[i] << " ";
	}
	std::cout << std::endl;

	// KMP(text, pattern);
	// ComputePrefix ababaca ---> 0 0 1 2 3 0 1
	// ComputePrefix 01111000110101011010 ---> 0 0 0 0 0 1 1 1 2 3 1 2 1 2 1 2 3 1 2 1
	// ComputePrefix aaaaaaaabbbbbbbb ---> 0 1 2 3 4 5 6 7 0 0 0 0 0 0 0 0
	// ComputePrefix abab ---> 0 0 1 2

	return 0;
}