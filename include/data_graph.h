#ifndef GRAPH_H_
#define GRAPH_H_

#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <unordered_map>
#include "util.h"

using namespace std;

//use for making binary
struct RawDataGraph {
	int max_vl_, max_el_;
	vector<int> vl_cnt_, el_cnt_;

	vector<vector<int>> vlabels_;
	vector<Edge> out_edges_; 
	vector<Edge> in_edges_; 

	vector<int> offset_;
	vector<int> label_;
	vector<int> adj_offset_;
	vector<int> adj_;

	vector<int> in_offset_;
	vector<int> in_label_;
	vector<int> in_adj_offset_;
	vector<int> in_adj_;

	vector<int> vl_offset_;
	vector<int> vl_;

	vector<vector<int>> vl_rel_;
#ifdef INTERSECTION
	vector<vector<int>> el_rel_;
	vector<vector<int>> in_el_rel_;
#else
	vector<vector<pair<int, int>>> el_rel_;
#endif
};

class DataGraph {
private:
	int vnum_, enum_, vl_num_, el_num_; 

	vector<int> vl_cnt_, el_cnt_;
	
	size_t encode_size_;

	const int* offset_; //data vertex id -> offset 
	const int* label_;
	const int* adj_offset_;
	const int* adj_;

	const int* in_offset_; //data vertex id -> offset 
	const int* in_label_;
	const int* in_adj_offset_;
	const int* in_adj_;

	const int* vl_offset_;
	const int* vl_; //not initialized in Alley

	const int* el_rel_offset_;
#ifdef INTERSECTION
	const int* el_rel_;
#else
	const pair<int, int>* el_rel_;
#endif

#ifdef INTERSECTION
	const int* in_el_rel_offset_;
	const int* in_el_rel_;
	int* label_bitmap_ = NULL;
#endif

	const int* vl_rel_offset_;
	const int* vl_rel_;
	
	RawDataGraph raw_;

	map<pair<int, int>, double> jvd_map_;
		
public:
	//build mode
	bool  HasBinary(const char*);
	void  ReadText(const char*);
	void  MakeBinary();
	size_t BinarySize();
	void  WriteBinary(const char*);
	void  ClearRawData();

	void  ReadBinary(const char*);
	range GetVertices(int);
	int   GetNumVertices();
	int   GetNumVertices(int);
#ifdef INTERSECTION
	range GetVertices(int, bool);
	int   GetNumVertices(int, bool);
	void  GroupELabels(vector<int>& l2g, vector<vector<int>>& g2l);
	void  GroupVLabels(vector<int>& l2g, vector<vector<int>>& g2l);
	inline bool PassFilter(int v, int& f) {
		return ((label_bitmap_[v] & f) == f);
	}
#endif
	int   GetNumEdges();
	int   GetNumEdges(int);
	int   GetNumVLabels(int = -1); 
	int   GetNumELabels(int = -1, bool = true); 
	range GetVLabels(int);
	range GetELabels(int, bool);
	int   GetELabelIndex(int, int, bool);
	range GetAdj(int, int, bool);
	int   GetAdjSize(int, int, bool);
	bool  HasEdge(int, int, int, bool);
	bool  HasVLabel(int, int);
	bool  HasELabel(int, int, bool);
	vector<int> GetRandomVertex(int);
	vector<int> GetVertex(int, int);
	vector<int> GetRandomEdge(int);
	vector<int> GetEdge(int, int);
	double GetJVD(int, int);
	vector<int> GetRandomEdge(int, int, bool);
};

#endif
