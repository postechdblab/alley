#include <fstream>
#include <cassert>
#include <iostream>
#include <string>
#include "../include/query_graph.h"

void QueryGraph::ReadText(const char* fn) {
    vnum_ = enum_ = 0;
    edge_.clear();
    adj_.clear();
    in_adj_.clear();
    vl_.clear();
    bound_.clear();
	ifstream input(fn);
	string line;
	int max_vl = -1, max_el = -1;
    while (getline(input, line)) {
        auto tok = parse(line, " ");
        if (tok[0][0] == 'v') {
			int vl = stoi(tok[2]);
			max_vl = std::max(max_vl, vl);
			vl_.push_back(vl);
			int bound = stoi(tok[3]);
			bound_.push_back(bound);
			vnum_++;
        }
        if (tok[0][0] == 'e') {
            if (adj_.size() == 0) {
                adj_.resize(vnum_);
                in_adj_.resize(vnum_);
                adj_elem_.resize(vnum_);
            }
            int src = stoi(tok[1]);
            int dst = stoi(tok[2]);
            int el  = stoi(tok[3]);
            assert(src < vnum_);
            assert(dst < vnum_);
            max_el = std::max(max_el, el);
            edge_.emplace_back(src, dst, el);
            adj_[src].push_back(make_pair(dst, el));
            in_adj_[dst].push_back(make_pair(src, el));
            adj_elem_[src].emplace_back(dst, el, true);
            adj_elem_[dst].emplace_back(src, el, false);
            enum_++;
        }
    }
    if (adj_.size() == 0)
        adj_.resize(vnum_);
    if (in_adj_.size() == 0)
        in_adj_.resize(vnum_);
	vl_num_ = max_vl + 1;
	el_num_ = max_el + 1;
}

void QueryGraph::ReadText(std::vector<std::string> &text) {
    vnum_ = enum_ = 0;
    edge_.clear();
    adj_.clear();
    in_adj_.clear();
    vl_.clear();
    bound_.clear();
	string line;
	int max_vl = -1, max_el = -1;
    for (auto &line: text) {
        auto tok = parse(line, " ");
        if (tok[0][0] == 'v') {
			int vl = stoi(tok[2]);
			max_vl = std::max(max_vl, vl);
			vl_.push_back(vl);
			int bound = stoi(tok[3]);
			bound_.push_back(bound);
			vnum_++;
        }
        if (tok[0][0] == 'e') {
            if (adj_.size() == 0)
                adj_.resize(vnum_);
            if (in_adj_.size() == 0)
                in_adj_.resize(vnum_);
            int src = stoi(tok[1]);
            int dst = stoi(tok[2]);
            int el  = stoi(tok[3]);
            assert(src < vnum_);
            assert(dst < vnum_);
            max_el = std::max(max_el, el);
            edge_.emplace_back(src, dst, el);
            adj_[src].push_back(make_pair(dst, el));
            in_adj_[dst].push_back(make_pair(src, el));
            enum_++;
        }
    }
    if (adj_.size() == 0)
        adj_.resize(vnum_);
    if (in_adj_.size() == 0)
        in_adj_.resize(vnum_);
	vl_num_ = max_vl + 1;
	el_num_ = max_el + 1;
}

vector<pair<int, int>>& QueryGraph::GetAdj(int v, bool dir) { 
    return dir ? adj_[v] : in_adj_[v]; 
}

vector<AdjElem>& QueryGraph::GetAdj(int v) {
	return adj_elem_[v];
}

int QueryGraph::GetAdjSize(int v) {
	return adj_elem_[v].size();
}

//assume one edge label between any two query vertices
int QueryGraph::GetELabel(int u, int v) {
	for (auto& e : adj_[u]) {
		if (e.first == v)
			return e.second;
	}
	return -1;
}

bool QueryGraph::IsStar() {
	for (int i = 0; i < GetNumVertices(); i++) {
		if (adj_elem_[i].size() == enum_)
			return true;
	}
	return false;
}
	
void QueryGraph::DeleteEdge(int i) {
	Edge e = edge_[i];
	edge_.erase(edge_.begin() + i);
	enum_--;
	for (int j = 0; j < adj_[e.src].size(); j++) {
		if (adj_[e.src][j].first == e.dst && adj_[e.src][j].second == e.el) {
			adj_[e.src].erase(adj_[e.src].begin() + j);
			break;
		}
	}
	for (int j = 0; j < in_adj_[e.dst].size(); j++) {
		if (in_adj_[e.dst][j].first == e.src && in_adj_[e.dst][j].second == e.el) {
			in_adj_[e.dst].erase(in_adj_[e.dst].begin() + j);
			break;
		}
	}
	for (int j = 0; j < adj_elem_[e.src].size(); j++) {
		if (adj_elem_[e.src][j].id == e.dst && adj_elem_[e.src][j].label == e.el && adj_elem_[e.src][j].dir) {
			adj_elem_[e.src].erase(adj_elem_[e.src].begin() + j);
			break;
		}
	}
	for (int j = 0; j < adj_elem_[e.dst].size(); j++) {
		if (adj_elem_[e.dst][j].id == e.src && adj_elem_[e.dst][j].label == e.el && !adj_elem_[e.dst][j].dir) {
			adj_elem_[e.dst].erase(adj_elem_[e.dst].begin() + j);
			break;
		}
	}
}

void QueryGraph::AddVertex() {
	vl_.push_back(-1);
	bound_.push_back(-1);
	vnum_++;
	adj_.resize(adj_.size() + 1);
	in_adj_.resize(in_adj_.size() + 1);
	adj_elem_.resize(adj_elem_.size() + 1);
}

void QueryGraph::AddEdge(Edge e) {
	edge_.push_back(e);
	enum_++;
	adj_[e.src].push_back(make_pair(e.dst, e.el));
	in_adj_[e.dst].push_back(make_pair(e.src, e.el));
	adj_elem_[e.src].emplace_back(e.dst, e.el, true);
	adj_elem_[e.dst].emplace_back(e.src, e.el, false);
}
