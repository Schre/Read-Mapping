#pragma once
#include <stdio.h>
#include <iostream>
#include <string>
#include <math.h>

using std::string;

long int max2(long int arg1, long int arg2){
  return (arg1 > arg2) ? arg1 : arg2;
}
long int max3(long int arg1, long int arg2, long int arg3){
  return max2(max2(arg1,arg2),max2(arg2,arg3));
}
std::pair<int,int> OutputAlignment(string & s1, string & aligned, int outputLength, int & start, int m, int mi, int h, int g, int isLocal, bool print=false){
  std::pair<int,int> ret;
  int s1gaps = 0, alignedgaps = 0;
  int s1indx = 0, alignedindx = 0;
  int matches = 0, mismatches = 0, gaps = 0, opening_gaps = 0;
  
  while (1){
    int t1 = s1indx, t2 = alignedindx;
    
    if (start + outputLength >= s1.length()){
      outputLength = s1.length() - start; // remain
    }
    if (outputLength <= 0)
      break;
    for (int i = start; i < outputLength+start; ++i){
      if (s1[i] == '-'){
	if (i == 0 || s1[i-1] != '-')
	  opening_gaps++;
	s1gaps++;
      }
      else
	++s1indx;
    }
    if (print){
      printf("s1\t %d\t %s", t1+1, s1.substr(start, outputLength).c_str());
      printf("\t %lu\n", (unsigned long)(start + outputLength-s1gaps));
      printf("  \t  \t ");
    }
    for (int i = start; i < outputLength+start; ++i){
      if (s1[i] == aligned[i]){
	++ret.first;
	if (print)
	  printf("|");
	++matches;
      }
      else{
	if (print)
	  printf(" ");
	++mismatches;
      }
      if(aligned[i] == '-'){
	if (i == 0)
	  opening_gaps++;
	else if (aligned[i-1] != '-')
	  opening_gaps++;
	
	alignedgaps++;
      }
      else
	++alignedindx;
    }
    if (print){
      printf("\n");
      printf("s2\t %d\t %s\t %lu\n\n", t2+1,aligned.substr(start,outputLength).c_str(),(unsigned long)(start + outputLength - alignedgaps));
    }
    start+= outputLength;
  }
  gaps = s1gaps + alignedgaps;
  mismatches-= gaps;
  int totalLen = mismatches + matches + gaps;
  if (print)
    std::cout << std::endl;
  if (print){
  std::cout << "Report\n\n";
  if (!isLocal)
    std::cout << "Global ";
  else
    std::cout << "Local ";
  std::cout << "optimal score = " << (matches*m)+(mismatches*mi)+(gaps*g)+(opening_gaps*h) << "\n\n";
  std::cout << "Number of: matches = " << matches << ",mismatches = " << mismatches << ",gaps = " << gaps << ",opening gaps = " << opening_gaps << std::endl;
  printf("Identities = %d/%d (%d%%), Gaps = %d/%d = (%d%%)\n",matches,totalLen,(int)(round(100*(double)matches/(totalLen))), gaps, totalLen, (int)(round(100*(double)gaps/(totalLen))));
  }
  ret.second = matches + mismatches + gaps;
  if (print)
    std::cout << "Matched: " << ret.first << " Aligned Length: " << ret.second << '\n';
  return ret;
}
