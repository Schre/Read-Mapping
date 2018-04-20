#pragma once
#include <vector>
#include <string>
#include "Utilities.h"

#define GLOBAL false
#define LOCAL true
#define SUB 0
#define DEL 1
#define INS 2

using std::string;
using std::vector;
using std::cout;

struct DP_cell{
  int S;
  int D;
  int I;
};

typedef struct DP_cell dp_cell;

class GL_Alignment{
 public:
  GL_Alignment(int mv = 1, int miv = -2, int hv = -5, int gv = -2){
    this->match = mv;
    this->mismatch = miv;
    this->h = hv;
    this->g = gv;
  }
  std::vector<string> GlobalAlignment(string &s1, string &s2){
    return Align(s1,s2, false);
  }
  std::vector<string> LocalAlignment(string &s1, string &s2){
    return Align(s1,s2,true);
  }
 private:
  // PRIVATE METHODS
  int getScore(dp_cell & c){
  return max3(c.S,c.D,c.I);
}
  int Substitute(char c1, char c2){
    return (c1 == c2) ? match : mismatch;
  }
  dp_cell ** CreateCellArray(unsigned long m, unsigned long n){
    dp_cell **DP = (dp_cell **)malloc(sizeof(dp_cell *) *(m+1));
    for (int i = 0; i < m+1; ++i){DP[i] = (dp_cell *)malloc(sizeof(dp_cell)*(1+n));}
    return DP;
  }
  void InitArray(dp_cell **DP, unsigned long m, unsigned long n, bool flag = 0){
    if (flag)
      return;
    DP[0][0].S = 0;
    DP[0][0].I = DP[0][0].D = INT32_MIN;
    
    for (int i = 1; i < m+1; ++i){
      DP[i][0].I = DP[i][0].S = INT32_MIN;
      DP[i][0].D = h + g*i;
    }
    for (int j = 1; j < n+1; ++j){
      DP[0][j].D = DP[0][j].S = INT32_MIN;
      DP[0][j].I = h + g*j;
    }
  }
  int computeDir(int sval, int ival, int dval){
    int maxVal = max3(sval, ival, dval); // get first direction
    int dir;
    if (maxVal == sval)
      dir = SUB;
    else if (maxVal == ival)
      dir = INS;
    else
      dir = DEL;
    return dir;
  }
  std::vector<string> Backtrace (dp_cell ** DP, string &s1, string &s2, int *i, int * j, int f = GLOBAL){
    std::vector<string> ret;
    string aligneds1 = "";
    string aligneds2 = "";
    int val;
    bool flag = false;
    int numMatches = 0;
    /*    while (1){
      flag = false;
      if (f == LOCAL && getScore(DP[*i][*j]) == 0)
	break;
    */
      // determine direction obtained value from
        int dir = computeDir(DP[*i][*j].S,DP[*i][*j].I,DP[*i][*j].D);
    
    // start main loop
    while (1){
      if ((*i == 0 && *j == 0) || (f == LOCAL && max3(DP[*i][*j].S,DP[*i][*j].D,DP[*i][*j].I) == 0)){
	break;
      }
    switch (dir){
    case INS:
      // determine where to go next
      if (*i == 0)
	dir = INS;
      else
	dir = computeDir(DP[*i][*j-1].S+h+g, DP[*i][*j-1].I+g, DP[*i][*j-1].D+h+g);
 
      aligneds2.insert(aligneds2.begin(),s2[*j-1]);
      aligneds1.insert(aligneds1.begin(),'-');
      (*j)--;
      break;
    case SUB:
      dir = computeDir(DP[*i-1][*j-1].S, DP[*i-1][*j-1].I, DP[*i-1][*j-1].D);
      
      aligneds2.insert(aligneds2.begin(),s2[*j-1]);
      aligneds1.insert(aligneds1.begin(),s1[*i-1]);
      (*i)--; (*j)--;
      break;
    case DEL:
      if (*j == 0)
	dir = DEL;
      else
	dir = computeDir(DP[*i-1][*j].S+h+g, DP[*i-1][*j].I+h+g, DP[*i-1][*j].D+g);
      
      aligneds2.insert(aligneds2.begin(),'-');
      aligneds1.insert(aligneds1.begin(), s1[*i-1]);
      (*i)--;
      break;
    }
    }
    
    ret.push_back(aligneds1);
    ret.push_back(aligneds2);
    return ret;
  }
  
  vector<string> Align(string & s1, string & s2, bool flag); // if flag == 1 local align 
  //

  
  int match;
  int mismatch;
  int g;
  int h;
};

std::vector<string> GL_Alignment::Align(string & s1, string & s2, bool flag = 0){
  int m = s1.length();
  int n = s2.length();
  int option;
  std::cout << "KEY: \tMatch: " << match << "\tMismatch: "<< mismatch << "\tG: " << g << "\tH: " << h << std::endl; 
  // Set option to 0 for local, otherwise -inifinity
  (flag == 0) ? option = INT32_MIN : option = 0;
  dp_cell ** DP = CreateCellArray(m,n);
  InitArray(DP, m, n, flag);
  int row, col;
  unsigned int maxsofar = 0;
  // Begin Algorithm
  for (int i = 1; i < m+1; ++i){
    for (int j = 1; j < n+1; ++j){
      // Substitute score
      DP[i][j].S = max2(getScore(DP[i-1][j-1]) + Substitute(s1[i-1],s2[j-1]),option);
      // Insert score
      if (j == 1){
	DP[i][j].I = max2(DP[i][j-1].D+h+g, option);
      }
      else
	DP[i][j].I = max2(max3(DP[i][j-1].I+g,
			  DP[i][j-1].S+h+g,
			       DP[i][j-1].D+h+g),
			  option);
      // Delete score
      if (i == 1)
	DP[i][j].D = max2(DP[i-1][j].I+h+g,option);
      else
	DP[i][j].D = max2(max3(DP[i-1][j].D+g,
			  DP[i-1][j].S+h+g,
			       DP[i-1][j].I+h+g),
			  option);
      
      if (flag && getScore(DP[i][j]) != 0){
	if (maxsofar < getScore(DP[i][j])){
	    maxsofar = getScore(DP[i][j]);
	    row = i;
	    col = j;
	  }
      }
	//printf("%d\t ", getScore(DP[i][j]));
    }
      //std::cout << std::endl;
  }
  if (!flag){
    row = m;
    col = n;
  }
  std::cout << "SCORE: " << getScore(DP[row][col]) << '\n';
  std::vector<string> ret;
  char c = 0;
  ret = Backtrace(DP,s1,s2,&row,&col, flag);
  OutputAlignment(ret[0], ret[1], 60, row, col,  match, mismatch, h, g, flag);
  return ret;
}
