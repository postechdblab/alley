#ifndef QUERY_GRAPH_H_
#define QUERY_GRAPH_H_

#include <vector>
#include <string>
#include "util.h"

class QueryGraph {
private:
	int vnum_, enum_, vl_num_, el_num_; 
	vector<Edge> edge_;
	vector<vector<pair<int, int>>> adj_, in_adj_;
	vector<int> vl_;
	vector<int> bound_;
	vector<vector<AdjElem>> adj_elem_;

	friend class Pattern;
public:
	QueryGraph() { 
		vnum_ = enum_ = vl_num_ = el_num_ = 0;
		edge_.clear();
		adj_.clear();
		in_adj_.clear();
		vl_.clear();
		bound_.clear();
		adj_elem_.clear();
	}
	void ReadText(const char*);
    void ReadText(std::vector<std::string>&);
	inline int GetNumVertices() { return vnum_; }
	inline int GetNumEdges() { return enum_; }
	Edge GetEdge(int i) { return edge_[i]; }; 
	vector<pair<int, int>>& GetAdj(int, bool);
	vector<AdjElem>& GetAdj(int);
	int GetAdjSize(int);
	int GetELabel(int, int);
	inline int GetVLabel(int v) { return vl_[v]; }
	inline int GetBound(int v) { return bound_[v]; }

	bool IsStar();
	void DeleteEdge(int i); 
	void AddVertex();
	void AddEdge(Edge e); 
};

#endif
