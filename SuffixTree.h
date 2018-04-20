/* TODO: Replace vector with map in node outgoing for easier access to indices */

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <map>

using std::string;
using std::pair;
using std::vector;
using std::map;

// Use map for easy lexicographical ordering and fast member access ( Theta(log n) )
// usage of Rope class motivated by wanting O(1) substr operation on strings and also for space efficiency

class Rope{
public:
	static string * str;

	Rope(){ 
		startIndx = 0;
		m_length = str->length();
	}
	Rope(size_t start, size_t len=str->length()){
		startIndx = start;
		m_length=len;
	}
	void print() const{
		for (size_t i = 0; i < m_length; ++i)
			putchar(this->at(i));
	}
	char at(size_t const index) const{
		if (index >= m_length || index < 0)
			return 0;
		return str->at(index + startIndx);
	}
	Rope substr(size_t const index, size_t const len) const{
		return Rope(index+startIndx, len);
	}
	string cpp_string() const{
		string ret = "";
		for (size_t i = 0; i < m_length; ++i)
			ret += this->at(i);
		return ret;
	}
	int start(){
		return startIndx;
	}
	int length(){
		return m_length;
	}
	bool operator==(Rope & rhs){
		if (this->startIndx == rhs.startIndx && this->m_length == rhs.m_length)
			return true;
		return false;
	}
	std::ostream& operator<<(std::ostream& os){
		this->print();
		return os;
	}
private:
	size_t startIndx;
	size_t m_length;
};
string * Rope::str = nullptr;

class Node{
public:
	Node(){
		m_parentEdge = std::make_pair(nullptr,Rope(0,0));
		m_sl = nullptr;
		start_leaf_index = end_leaf_index = -1;
	}
	map<char, pair<Node *, Rope>> m_outgoing; // Node, edge label
	Node * m_sl;
	std::pair<Node *, Rope> m_parentEdge;
	int start_leaf_index;
	int end_leaf_index;
	/* Note: start/end_leaf_index refers to the leafs underneath the current node
	   	 corresponding to the indexes in the leaflist. Useful for ReadMapping. If not
		 performing a ReadMap, remove these */
};

class SuffixTree{
public:
	SuffixTree(string * s){
		Rope::str = s; // initialize rope class static member str
		// Initialize root and add the first string
		Rope t;
		root = new Node();
		auto p = std::make_pair(root, Rope());
		root->m_parentEdge = p;
		Node * tempNode = new Node();
		root->m_outgoing['$'] = std::make_pair(tempNode,Rope(0,0));
		tempNode->m_parentEdge = p;
		//DoNaiveConstruction(); <--- Function is broken... it was too slow anyways.
		LinearTimeConstruction();
		TrimTree(root); // Delete superfluous nodes (artifacts of findpath)
		outputTree(root, "");
		PrepareST(root, 0);
		//		FindLoc("ss", 0);
		/*		std::cout << '\n';
		std::cout << "START BWT: ******************\n";
		for (int i = 0; i < BWT.size(); ++i)
		  std::cout << BWT[i] << '\n';
		std::cout << "END BWT:   ******************\n\n";
		getLongestSubStringLocs();*/			       		
	}

	/* BEGIN MAP READ FUNCTIONS */

	/* IMPORTANT! MUST CALL PrepareST() PRIOR TO CALLING THIS FUNCTION */
	void MapReads(std::vector<string> reads, int x){
	  int m = reads.size();
	  for (int i = 0; i < m; ++i){
	    /* 3 A */
	    int l = reads[i].length();
	    string & r_i = reads[i];
	    /* END 3 A */

	    /* 3 B */
	    auto starts = FindLoc(r_i, x);
	    std::cout << reads[i] << ':';
	    if (starts.size() == 0){std::cout << '\n'; // insufficient evidence that the read is a contained in our main sequence
	      continue;
	    }
	    for (int j : starts){
	      std::cout  << j << ';';
	    }
	    std::cout << '\n';
	    /* END 3 B */	   
	  }
	}

	/* TODO: Convert deepestNode into a list of potential deepestNode candidates. THEN, remove all nodes in that list that are less than the longest length found */
	std::vector<int> FindLoc(string s, int x){
	  // Exhaust characters in s
	  Node * curNode = root;
	  int strIndx = 0;
	  Node * deepestNode;
	  int deepestNodeLength = 0;
	  
	  auto outgoing = root->m_outgoing; // start from root
	  Node * prevNode = root;
	  int curLength = 0;
	  int prevIndx = 0;
	  int excess = 0;
	  
	  while(strIndx < s.size()){
	    int offset = 0;
	    if (outgoing.find(s[strIndx]) != outgoing.end()){ // if edge exists...
	      auto edge = outgoing[s[strIndx]];
	      prevNode = edge.first;
	      // has character... exhaust as many chars as possible
	      while (strIndx < s.length() && offset < edge.second.length() && s[strIndx] == edge.second.at(offset)){++strIndx; ++offset; ++curLength;}
	      if (curLength > deepestNodeLength){
		deepestNodeLength = curLength;
		deepestNode = prevNode;
		excess = edge.second.length() - offset;
	      }		

	      if (offset < edge.second.length()){ 
		// in middle
		// RESET
		curLength = 0; // reset curLength
		outgoing = root->m_outgoing; //  back to top of tree
		strIndx = prevIndx + 1; // try next index
		prevNode = root;
		++prevIndx;
	      }
	      else{
		outgoing = edge.first->m_outgoing;
	      }
	      continue;
	    }
	    // RESET
	    if (curLength > deepestNodeLength){
	      deepestNodeLength = curLength;
	      deepestNode = prevNode;
	    }
	    curLength = 0; // reset curLength
	    outgoing = root->m_outgoing; //  back to top of tree
	    strIndx = prevIndx + 1; // try next index
	    prevNode = root;
	    ++prevIndx;
	  }

	  //	  std::cout << "LCS: ";
	  //	  	  string lcs = "";
		  /*	   while (deepestNode != root && deepestNode){
			   lcs.insert(0,deepestNode->m_parentEdge.second.cpp_string());	      
			   deepestNode = deepestNode->m_parentEdge.first;
			   }*/
		  //std::cout << lcs;
	  
	  std::vector<int> ret;
	  //	  std::cout << " " << deepestNodeLength << '\n';
	  if (deepestNodeLength < x){
	    std::cout << "Insufficient number of characters matched to perform local alignment \n";
	    return ret; // empty
	  }
	  for (auto p : deepestNode->m_outgoing){
	    auto edge = p.second;
	    ret.push_back(edge.second.start() - deepestNodeLength - excess);
	  }
	  // NOTE: IF NOT OUTGOING EDGES, MUST BE A NON-REPEATED SUFFIX. THEN START LOCATION IS SequenceLength - deepestNodeLength
	  if (deepestNode->m_outgoing.size() == 0)
	    ret.push_back(Rope::str->length() - deepestNodeLength);
	   return ret;
	}
	/* END MAP READ FUNCTIONS */
	void getLongestSubStringLocs(){
	  auto res = getDeepestInternalNode(std::make_pair(root,0));
	  Node * deepestNode = res.first;
	  Node * deepestNodecpy = res.first;

	  std::cout << "Most repeated substring: ";
	  string longest = "";
	  while(deepestNode != root){
	    string t = longest;
	    longest = deepestNode->m_parentEdge.second.cpp_string() + t;
	    deepestNode = deepestNode->m_parentEdge.first;
	  }
	  int LCS_len = res.second;
		
	  std::cout << longest << '\n';
	  std::cout << "Length: " << LCS_len << '\n';
		
	  /* To do: Get all ends by doing for each edge in deepest node and taking their end indeces and doing (edge end - edge len - lcs len) */
	  //	LCS->print();
	  std::cout << "Starting positions:\n";
	  for (auto p : deepestNodecpy->m_outgoing){
	    auto edge = p.second;
	    std::cout << edge.second.start()-LCS_len + 1<< '\n';
	  }

	}

	// START MAY NOT BE NULL!
	void LinearTimeConstruction()
	{
		root->m_sl = root;
		Node * v = nullptr;
		Node * currentLeaf = nullptr;
		Node * u = root; // u internal node is initialized to be the root
		int alphaLen = 0;
		// add first string...
		Rope r(0,Rope::str->length());
		currentLeaf = FindPath(root, r)->m_parentEdge.first;
		for (int i = 1;i < Rope::str->length();++i){
		  Node * u_prime = u->m_parentEdge.first;
		  // CASE IA
		  if (u->m_sl != nullptr && u != root){
		    currentLeaf = HandleIA(u, alphaLen, i);
		  }
		  // CASE IB
		  else if (u->m_sl != nullptr && u == root){
		    currentLeaf = HandleIB(u, alphaLen, i);
		  }
		  // CASE IIA
		  else if (u->m_sl == nullptr && u_prime != root){
		    currentLeaf = HandleIIA(u, alphaLen, i);
		  }
		  // CASE IIB
		  else if (u->m_sl == nullptr && u_prime == root){
		    currentLeaf = HandleIIB(u, i);
		  }
		  if (currentLeaf == nullptr){
		    std::cout << "Leaf is nullptr... exit on iteration " << i << "out of " << Rope::str->length() << '\n';
		    //exit(1);
		    return;
		  }

		  // set u to internal node
		  u = currentLeaf->m_parentEdge.first;
		  // find the length of alpha (if current leaf's parent is root, the alpha is 0)
		  
		  alphaLen = (currentLeaf->m_parentEdge.first != root) ?
		    ((int)Rope::str->length() - (i+1)) - currentLeaf->m_parentEdge.second.length()
		    : 0;
		  if (currentLeaf->m_parentEdge.second.at(0) == '$')
		    ++alphaLen;
		  if (alphaLen < 0)
		    alphaLen += Rope::str->length();
		}
	}

	/* Functions for ReadMapping */
	int nextIndex = 0;
	void PrepareST(Node *& start, int slen){ // dfs procedure
	  // LEAF CASE
	  if (start->m_outgoing.size() == 0){ // Is a leaf node
	    int i = (int)(Rope::str->length()) - slen;
	    LeafList.push_back(i); // Leaf index
	    start->start_leaf_index = start->end_leaf_index = nextIndex++;
	    char c = (i==0) ? '$' : Rope::str->at(i-1);
	    if (c != '\n') // if string is from file, may be terminated by end line char
	      BWT.push_back(c);
	    //	    std::cout << "<" << start->start_leaf_index << "," << start->end_leaf_index << ">" << '\n';
	    return;
	  }
	   

	  bool firstIt = true;
	  start->start_leaf_index = start->end_leaf_index = nextIndex;
	  for (auto c : start->m_outgoing){
	    auto edge = c.second;
	    PrepareST(edge.first, slen + edge.second.length());
	    if (start->start_leaf_index > nextIndex - 1 || start->start_leaf_index == -1)
	      start->start_leaf_index = nextIndex - 1;
	    if (start->end_leaf_index < nextIndex - 1)
	      start->end_leaf_index = nextIndex - 1;
	  }  
	  /* std::cout << "<" << start->start_leaf_index << "," << start->end_leaf_index << ">";
	  	  for (int i = start->start_leaf_index; i <= start->end_leaf_index; ++i){
	    std::cout << LeafList[i]+1 << ",";
	    }
	  std::cout << '\n';
	  */
	}

 private:
	void TrimTree(Node * u){
	  for (auto c : u->m_outgoing){
	    auto edge = c.second;
	    if (edge.first->m_outgoing.size() == 1 && edge.first->m_outgoing.find('$') != edge.first->m_outgoing.end()){
	      // superfluous node
	      edge.first->m_outgoing.erase('$');
	      continue;
	    }
	    TrimTree(edge.first);
	  }
	}
	Node * HandleIA(Node * u, int alphaLen, int i){
		// 2.) Take SL(u) to v
		Node * v = u->m_sl;
		// 3.) Add suffix
		Rope r(i+alphaLen, Rope::str->length()-(i+alphaLen));
		return FindPath(v, r);

	}
	Node * HandleIB(Node * u, int alphaLen, int i){
		// 2.) Take SL(u) => root => v
		Node * v = u->m_sl;
		// 3.) Add Suffix
		Rope r(i+alphaLen, Rope::str->length()-(i+alphaLen));
		return FindPath(v, r);
	}
	Node * HandleIIA(Node * u, int alphaLen, int i){
	  // 2.) Go to u', and also record beta <-- edge label from u' to u
	  Node * u_prime = u->m_parentEdge.first;
	  Rope beta = u->m_parentEdge.second;
	  // 3.) Take Suffix Link
	  Node * v_prime = u_prime->m_sl;
	  // 4.) Go to v by tracing beta (end of alpha)
	  Node * v = NodeHop(v_prime, beta);
	  // 5.) Set sl(u) to v
	  u->m_sl = v;
	  // 6.) Find path
	  Rope r(i+alphaLen, Rope::str->length()-(i+alphaLen));
	  return FindPath(v, r);
	}
	Node * HandleIIB(Node * u,int i){
	  // 2.) Go to u', and also record beta <-- edge label from u' to u and beta'
	  Rope beta_prime(u->m_parentEdge.second);
	  beta_prime = beta_prime.substr(1, beta_prime.length() - 1);
	  Node * u_prime = u->m_parentEdge.first;
	  // 3.) Take Suffix Link
	  Node * v_prime = root;
	  // 4.) Go to v by tracing beta'
	  Node * v = NodeHop(root, beta_prime);
	  // 5.) Set sl(u) to v
	  u->m_sl = v;
	  // 6.) Find path
	  Rope r(i+beta_prime.length(), Rope::str->length()-(i+beta_prime.length()));
	  return FindPath(v, r);

	}
	Node * NodeHop(Node * node, Rope s){
	  int curIndex = 0;
	  int offset = 0;
	  while (curIndex < s.length()){
	    if (node->m_outgoing.find(s.at(curIndex)) != node->m_outgoing.end()){ // hop to node if there is an edge starting with s.at(0)

	      auto edge = node->m_outgoing[s.at(curIndex)];
	      if (curIndex + edge.second.length() <= s.length()){
		node = edge.first;
		curIndex += edge.second.length();
	      }
	      else{
		// create a new node to split string
		Node * newNode = SplitEdge(node, edge, s, offset, curIndex);
		return newNode;
	      }
	    }
	    else{
	      std::cout << "Fatal: In Node Hop\n";
	      exit(-1);
	    }
	    ++offset;
	  }
	  return node;
	}
	/*	void DoNaiveConstruction(){
	  root->m_sl = root;
	  Node * lastLeaf = nullptr;
	  for (int i = 0; i < s.length(); ++i){
	    Rope r(i, Rope::str->length()-i);
	    FindPath(root, r);
	  }
	  }*/
	Node * FindPath(Node * start, Rope s){ // Function creates or traces path
	  // Keep matching until cannot match characters, then diverge to create new suffix if need be
	  Node * it = nullptr; // it and indx will indicate which node it diverges and which index in that node
	  Node * ret = nullptr;
	  int lenMatch = 0;
	  int offset = 0;
	  Rope outString;

	  if (start->m_outgoing.find(s.at(0)) != start->m_outgoing.end()){
	    // found node to take
	    it = start->m_outgoing[s.at(0)].first;
	    outString = start->m_outgoing[s.at(0)].second;
	  }

	  if (it == nullptr){ // no out node to take.. make one and return it
	    if (s.length() != 0){
	      // Create node for label s
	      Node * temp = new Node();
	      // Set node's parent pointer to start
	      temp->m_parentEdge = std::make_pair(start,s);
	      auto p = std::make_pair(temp, s);
	      // scan to correct lexicographical ordering
	      start->m_outgoing[p.second.at(0)] = p; // add new node and edge with label s
	      ret = new Node();
	      // create backpointer to parent
	      ret->m_parentEdge = p;
	      p.first->m_outgoing['$'] = (std::make_pair(ret, Rope(s.start()+s.length(),0)));
	      ret = temp;
	    }
	    else{
	      ret = new Node();
	      // establish backpointer to parent
	      ret->m_parentEdge = std::make_pair(start,Rope(0,0));
	      start->m_outgoing['$'] = std::make_pair(ret, Rope(s.start(),0));
				
	    }

	  }
	  else{ // find longest common prefix in this outnode with string s
	    while (s.at(lenMatch) == outString.at(lenMatch) && lenMatch < (s.length())){lenMatch++;} // exhaust as many characters as possible
	    if(lenMatch < outString.length() && lenMatch < (s.length())){ // middle of outstring and not all characters exhausted
	      // Split outString
	      Rope firstSlice = s.substr(0, lenMatch); // matching slice
	      Rope secondSlice = outString.substr(lenMatch, outString.length() - lenMatch);

	      auto oldPair = start->m_outgoing[s.at(0)]; // save old pair, but set old pair's label to secondSlice
	      oldPair.second = secondSlice;
	      Node * temp = new Node();
	      // establish a backpointer to the parent
	      temp->m_parentEdge = std::make_pair(start,firstSlice); 
	      start->m_outgoing[firstSlice.at(0)].first = temp; 
	      start->m_outgoing[firstSlice.at(0)].second = firstSlice;
	      // set new node's parent to temp
	      oldPair.first->m_parentEdge = std::make_pair(temp, secondSlice);
	      temp->m_outgoing[secondSlice.at(0)] = (oldPair);
	      ret = FindPath(start->m_outgoing[s.at(0)].first, s.substr(lenMatch, s.length()-lenMatch));
	    }
	    else if (lenMatch == outString.length() && lenMatch < (s.length())){ // end of outstring and more characters to match
	      ret = FindPath(start->m_outgoing[s.at(0)].first, s.substr(lenMatch, s.length()-lenMatch));
	    }
	    else if (lenMatch < outString.length() && lenMatch == (s.length())){
	      // Split outString
	      Rope firstSlice = s.substr(0, lenMatch); // matching slice
	      Rope secondSlice = outString.substr(lenMatch, outString.length() - lenMatch);

	      auto oldPair = start->m_outgoing[firstSlice.at(0)]; // save old pair, but set old pair's label to secondSlice
	      oldPair.second = secondSlice;

	      Node * newNode = new Node();
	      Node * endNode = new Node();

	      endNode->m_parentEdge = std::make_pair(newNode, Rope(s.start() + s.length(),0));
	      newNode->m_parentEdge = std::make_pair(start, firstSlice);
	      oldPair.first->m_parentEdge = std::make_pair(newNode, secondSlice);

	      newNode->m_outgoing['$'] = (std::make_pair(endNode, Rope(s.start() + s.length(),0)));
	      start->m_outgoing[s.at(0)].first = newNode; 
	      start->m_outgoing[s.at(0)].second = firstSlice;
	      newNode->m_outgoing[secondSlice.at(0)] = (oldPair);
	      ret = endNode;
	    }
	    else {
	      if (it->m_outgoing.find('$') == it->m_outgoing.end()){
		Node * newNode = new Node();
		it->m_outgoing['$'] = std::make_pair(newNode,Rope(s.start() + s.length(),0));
		newNode->m_parentEdge = std::make_pair(it, Rope(s.start() + s.length(),0));
	      }
	      ret = it->m_outgoing['$'].first;
	    }
	  }
	  return ret;
	}

	Node * SplitEdge(Node * node, pair<Node*, Rope> edge, Rope & label, int edgeIndx, int startIndx){
	  Node * newNode = new Node();
	  Rope firstSlice = label.substr(startIndx, label.length()-startIndx); // remaining string			
	  Rope secondSlice = edge.second.substr(firstSlice.length(), edge.second.length()-firstSlice.length());
	  // now insert the new node
	  node->m_outgoing[firstSlice.at(0)].first = newNode;
	  node->m_outgoing[firstSlice.at(0)].second = firstSlice;
	  newNode->m_outgoing[secondSlice.at(0)] = (std::make_pair(edge.first, secondSlice));
	  // set parent pointers
	  newNode->m_parentEdge = std::make_pair(node, firstSlice);
	  edge.first->m_parentEdge = std::make_pair(newNode, secondSlice);
	  return newNode;
	}

	void outputTree(Node *& start, string s){ // dfs procedure
	  if (start->m_outgoing.size() == 0){ // Is a leaf node
	    std::cout << s << '\n';
	  }
	  for (auto c : start->m_outgoing){
	    auto edge = c.second;
	    outputTree(edge.first, s+edge.second.cpp_string());
	  }
	}
	std::pair<Node *, int> getDeepestInternalNode(std::pair<Node*,Rope> start, int depth=0){
	  std::pair<Node*, int> ret = std::make_pair(start.first, depth); 
	  
	  for (auto c : start.first->m_outgoing){
	    auto edge = c.second;
	    if (edge.first->m_outgoing.size() <= 1)
	      continue;
	    auto t = getDeepestInternalNode(edge, depth + edge.second.length());
	    if (t.second > ret.second){ // if internal node is deeper, choose this node
	      ret.first = t.first;
	      ret.second = t.second;
	    }
	  }
	  /*	  if (start.first->m_outgoing.find('$') != start.first->m_outgoing.end() && start.first->m_outgoing.size() == 1){ // if only outgoing edge is $ => at last level
	    // leaf node -> Get parent, if depth is deeper than ret then set ret to leaf's parent
	    auto parentEdge = start.first->m_parentEdge;
	    if (!ret.first  || (depth -  

	    }*/
	return ret;
	}
	void getBWT(Node *& start, int slen){ // dfs procedure
	  if (start->m_outgoing.size() == 0){ // Is a leaf node
	      int i = (int)(Rope::str->length()) - slen;
	      LeafList.push_back(i); // Leaf index
	      char c = (i==0) ? '$' : Rope::str->at(i-1);
	      if (c != '\n') // if string is from file, may be terminated by end line char
		BWT.push_back(c);
	    }

	  for (auto c : start->m_outgoing){
	    auto edge = c.second;	   
	    getBWT(edge.first, slen + edge.second.length());
	  }
	}
	std::vector<char> BWT;
	std::vector<int> LeafList; // Leaf indexes from left to right in ST
	Rope * LCS = nullptr;
	Node * root;
};
