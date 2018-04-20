#include <iostream>
#include "SuffixTree.h"
#include <fstream>
#include <set>
#include <ctime>

using std::istream;
using std::set;

set<char> Alphabet;

void initAlphabet(char * path){
  std::ifstream ifile(path);
  if (ifile.fail()){
    std::cout << "Failed to open alphabet file: " << path << '\n';
  }
  string buffer;
  
  while(std::getline(ifile, buffer))
    for (int i = 0; i < buffer.size(); ++i)
      if (buffer[i] != '\n' && buffer[i] != EOF)
	Alphabet.insert(buffer[i]);
}
void showAlphabet(){
  std::cout << "[ ";
  for (auto c : Alphabet)
    std::cout << c << " ";
  std::cout << "]\n";
}
bool inAlphabet(char c){ 
  return (Alphabet.find(c) == Alphabet.end()) ? false : true;
}
int main(int argc, char * argv[]){

  if (argc < 3){
    std::cout << "Usage: ./suffixTree <filepath> <alphabetpath>\n";
    exit(0);
  }
 
  // try to open the file
  std::ifstream inputFile;
  inputFile.open(argv[1]);
  if(inputFile.fail()){
    std::cout << "Failed to open input file: " << argv[1] << '\n';
    exit(-1);
  }
  initAlphabet(argv[2]);
  //  showAlphabet();
  string s;
  string buffer = "";
  
  while (std::getline(inputFile, buffer)){
    for (int i = 0; i < buffer.length(); ++i)
      if (inAlphabet(buffer[i])) // <--- O(1)
	s += buffer[i];
    buffer.clear();
  }
  std::clock_t time;
  time = std::clock();
  SuffixTree * t = new SuffixTree(&s);
  time = std::clock() - time;
  //std::cout << "Time: " << time/(double)(CLOCKS_PER_SEC/1000000) << " micro seconds\n";

  std::vector<string> reads = {"a", "an", "ana", "na", "nana", "ban", "bana"};
  t->MapReads(reads, 0);
  return 1;
}
