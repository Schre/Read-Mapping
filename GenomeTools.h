#include <iostream>
#include "SuffixTree.h"
#include "Alignment.h"
#include "Utilities.h"
#include <queue>

using std::queue;
class Engine{
 public:
  Engine(string &input){
    std::cout << "Initialize Read Mapping...\n";
    std::cout << "Building suffix tree...\n";
    Genome_ST = new SuffixTree(&input);
    std::cout << "Finished construction of suffix tree...\n";
  }
  
 private:
  SuffixTree * Genome_ST;
  queue<string> query_q;
}
