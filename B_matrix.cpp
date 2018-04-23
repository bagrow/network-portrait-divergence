//  B_matrix.cpp
//  
//  Jim Bagrow 
//  Last Modified: 2007-04-23

#include <iostream>
#include <fstream>
#include <string>
#include <vector> // this code makes heavy use of the 
#include <queue>  // STL, since I suck at actual C
#include <map>

using namespace std;

vector<vector<int> > readEdgeList(string fN){
	// Graph is stored as a vector of vectors of ints.  The ith element of the 
	// vector 'graph' is a vector containing the neighbors of the ith node.
	// Nodes must be integers numbered SEQUENTIALLY from zero.  A gap in the
	// numbering will be interpreted as a node of degree zero. 
	// Input is text file containing a space-separated edgelist
	
	// read through file once to find maximum node:
	int maxid = -1;
	int s,t;
	ifstream fscan(fN.c_str(), ios::in);
	if (!fscan) {
	    cerr << "Unable to open input file.\n";
	    exit(1);   // call system to stop
	}
	while (fscan >> s >> t) {
		if (s > maxid)
			maxid = s;
		if (t > maxid)
			maxid = t;
	}
	fscan.close();
	
	// allocate space for graph
	vector<vector<int> > graph (maxid+1, vector<int>(0)); // vector of empty vectors
		
	// load up adjacency lists for graph:
	ifstream fin(fN.c_str(), ios::in);
	while (fin >> s >> t) {
		graph[s].push_back(t);
		graph[t].push_back(s);
	}
	fin.close();
	
	return graph;
}

vector<int> BFS( const vector<vector<int> >& graph, int start ) {
	vector<int> search;
	queue<int>  next;
	vector<int> seen(graph.size(), -1); // filled with negative ones!
	
	search.push_back(start);	
	seen[start] = 0; // seen tracks the depth the node was first seen at (-1 = unseen)
	int depth = 1;
	while ( !search.empty() ){
		for (vector<int>::const_iterator i = search.begin(); i != search.end(); ++i) {
			int u = *i;
			for (vector<int>::const_iterator j = graph[u].begin(); j != graph[u].end(); ++j) {
				int v = *j;
				if (seen[v] == -1) { // If v is unvisited.
					seen[v] = depth;
					next.push(v);
				}
			}
		}
		search.clear();
		while( !next.empty() ){
			search.push_back( next.front() );
			next.pop();
		}
		depth++;
	}
	return seen;
}

int main (int argc, char const* argv[]){
	// parse command line args, needs to be better:
	string inputFileName, outputFileName;
	if(argc == 1){
		cout << "\nUsage: ./B_matrix input_file output_file"                << endl;
		cout << "\ninput_file is an M x 2 edgelist of integers, where the"  << endl;
		cout << "     integers 0,1,...,N-1 represent the nodes of the "     << endl;
		cout << "     graph. Each line must have two integers separated by" << endl;
		cout << "     a space.  These are the edges of the graph."          << endl;
		cout << "\noutput_file is the file where the B matrix will be"      << endl;
		cout << "     saved."                                               << endl;
		return 0;
	}
	else if(argc == 3){
		inputFileName  = argv[1];
		outputFileName = argv[2];
	}
	else{
		cout << "Wrong number of arguments, aborting..." << endl;
		return 1;
	}
	
	// get the graph:
	vector<vector<int> > graph;
	graph = readEdgeList(inputFileName);
	if( graph.size() == 0 ){
		cout << "Graph has zero nodes, aborting... (check input file?)" << endl;
		return 1;
	}
	cout << "There are " << graph.size() << " nodes in the graph." << endl;
	
	// initialize B matrix, all entries are zero:
 	int num_rows = 500; // should be diameter+1, this is an ugly hack!
	int num_cols = graph.size();
	vector<vector<int> > shell_mat;
	for (int i = 0; i < num_rows; ++i)
		shell_mat.push_back( vector<int>(num_cols,0) ); // vector of num_col zeros
	
	// BFS from all starting nodes to build B matrix:
	cout << "Beginning BFS... " << endl;
	map<int,int> count_at_depth;
	map<int,int>::iterator iter;
	vector<int> seen;
	int depth, count, curr_max_depth, i, j;
	int max_depth = 0;
	int N = graph.size();
	for( unsigned int starting_node = 0; starting_node < graph.size(); starting_node += 1 ){
		seen = BFS(graph, starting_node);
		
		count_at_depth.clear();
		for( i = 0; i < seen.size(); i += 1 ) //  seen's keys must be sequential, since graph!
			count_at_depth[ seen[i] ]++; // = 0 if trying to increment new key
		
		curr_max_depth = 0; // iterate over map, incrementing matrix and getting max depth
		for( iter = count_at_depth.begin(); iter != count_at_depth.end(); iter++ ) {
			depth = iter->first; count = iter->second;
			if( depth != -1){ // if graph is disconnected, these show up
				shell_mat[depth][count]++;
				if( depth > curr_max_depth ) curr_max_depth = depth;
			}
		}
		if( curr_max_depth > max_depth ) max_depth = curr_max_depth;
	}
		
	// fill in zeroth column of B matrix:
	int curr_row_sum;
	for( i = 0; i < max_depth+1; i += 1 ){
		curr_row_sum = 0;
		for( j = 1; j < graph.size(); j += 1 )
			curr_row_sum += shell_mat[i][j];
		shell_mat[i][0] = N - curr_row_sum;
	}
	
	// save B-matrix to file:
	cout << "Saving B matrix..." << endl;
	ofstream fout(outputFileName.c_str(), ios::out);
	for( i = 0; i < max_depth+1; i += 1 ){
		for( j = 0; j < graph.size(); j += 1 )
			fout << shell_mat[i][j] << " ";
		fout << endl;
	}
	fout.close();
	cout << "Done." << endl;
	
	// no need to clean up memory since everything is vectors... I hope.
	
	return 0;
}
