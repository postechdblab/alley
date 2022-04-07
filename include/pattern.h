#ifndef PATTERN_H_
#define PATTERN_H_

#include <algorithm>
#include "query_graph.h"

#define PATTERN_EMPTY_ELABEL 0
#define PATTERN_EMPTY_VLABEL -1

typedef vector<int> VID_MAP;

struct Pattern {
	//only store outgoing edges
	vector<vector<pair<int, int>>> adj; //u -> <(order, v)>, sorted by order; then sorted by v (increasing order)
	vector<vector<pair<int, int>>> in_adj; //don't need to be sorted
	vector<int> vorder; 
	bool directed;

	Pattern() {
		directed=true;
	}

	Pattern(int order, bool isEdgeLabel=true) {
		directed = true;
		if(isEdgeLabel) {
			adj.resize(2);
			adj[0].emplace_back(order, 1);
			in_adj.resize(2);
			in_adj[1].emplace_back(order, 0);
			vorder.resize(2);
			vorder[0]=PATTERN_EMPTY_VLABEL;
			vorder[1]=PATTERN_EMPTY_VLABEL;
		}
		else {
			adj.resize(1);
			in_adj.resize(1);
			vorder.emplace_back(order);
		}
	}

	Pattern(const Pattern& p) : adj(p.adj), in_adj(p.in_adj), directed(p.directed), vorder(p.vorder) {}

	Pattern(QueryGraph& q, vector<int>& elabel2eorder, vector<int>& vlabel2vorder) {
		directed = true;
		adj.resize(q.GetNumVertices());
		in_adj.resize(q.GetNumVertices());
		vorder.resize(q.GetNumVertices(), -1);
		for (int u = 0; u < q.GetNumVertices(); u++) {
			if(vlabel2vorder.size() > 0) {
				vorder[u] = q.GetVLabel(u) == PATTERN_EMPTY_VLABEL? PATTERN_EMPTY_VLABEL : vlabel2vorder[ q.GetVLabel(u) ];
			}
			for (auto& e : q.GetAdj(u)) {
				int o = elabel2eorder.size() > 0? elabel2eorder[e.label] : e.label;
				if (e.dir)
					adj[u].emplace_back(o, e.id);
				else
					in_adj[u].emplace_back(o, e.id);
			}
		}
		sort();
	}

	Pattern(QueryGraph& q) {
		directed = true;
		adj.resize(q.GetNumVertices());
		in_adj.resize(q.GetNumVertices());
		vorder.resize(q.GetNumVertices(), -1);
		for (int u = 0; u < q.GetNumVertices(); u++) {
			vorder[u] = q.GetVLabel(u);
			for (auto& e : q.GetAdj(u)) {
				if (e.dir)
					adj[u].emplace_back(e.label, e.id);
				else
					in_adj[u].emplace_back(e.label, e.id);
			}
		}
		sort();
	}

	QueryGraph convertToQuery() {
		QueryGraph q;
		q.vnum_ = adj.size();
		q.enum_ = getNumEdges();
		q.vl_.resize(adj.size());
		q.bound_.resize(adj.size(), -1);
		q.adj_elem_.resize(adj.size());
		q.adj_.resize(adj.size());
		q.in_adj_.resize(adj.size());
		for (int u = 0; u < adj.size(); u++) {
			q.vl_[u]=vorder[u];
			for (auto& e : adj[u]) {
				int el = e.first;
				//id, label, dir
				q.adj_elem_[u].emplace_back(e.second, el, true); 
				//id, label
				q.adj_[u].emplace_back(e.second, el);
				q.edge_.emplace_back(u, e.second, el); 
			}
			for (auto& e : in_adj[u]) {
				int el = e.first; 
				//id, label, dir
				q.adj_elem_[u].emplace_back(e.second, el, false); 
				//id, label
				q.in_adj_[u].emplace_back(e.second, el);
			}
		}
		return q;
	}

	void convertVLabels(vector<int>& label2order) {
		for (int u = 0; u < getNumVertices(); u++) {
			assert(vorder.size() > u);
			
			vorder[u] = vorder[u] == PATTERN_EMPTY_VLABEL? PATTERN_EMPTY_VLABEL : label2order[ vorder[u] ];
		}
	}

	void convertELabels(vector<int>& label2order) {
		for (int u = 0; u < getNumVertices(); u++) {
			assert(adj.size() > u && in_adj.size() > u);
			for (auto& e : adj[u]) {
				assert(label2order.size() > e.first);
				e.first = label2order[e.first];
			}
			for (auto& e : in_adj[u]) {
				assert(label2order.size() > e.first);
				e.first = label2order[e.first];
			}
		}
		sort();
	}

	void removeAnEdge() {
		for (int u = 0; u < getNumVertices(); u++) {
			for (auto& e : adj[u]) {
				deleteEdge(u, e.second, e.first);
				return;
			}
			for (auto& e : in_adj[u]) {
				deleteEdge(e.second, u, e.first);
				return;
			}
		}
		assert(false);
	}

	bool isPath() {
		if(!isTree()) {
			return false;
		}
		for (int u = 0; u < getNumVertices(); u++) {
			int deg = getDegree(u);
			if (deg > 2)
				return false;
		}
		return true;
	}
	bool isTree() {
		return getNumVertices() == getNumEdges() + 1;
	}

	bool isDirected() {
		return directed;
	}

	bool ifPathGetLeaves(vector<int>& leaves) {
		for (int u = 0; u < getNumVertices(); u++) {
			int deg = getDegree(u);
			if (deg > 2)
				return false;
			else if (deg == 1)
				leaves.push_back(u);
		}
		if (leaves.empty())
			return false;
		assert(leaves.size() == 2);
		return true;
	}

	bool hasDirectedEdge(int u, int v) {
		assert(adj.size() > u);
		for (auto& e : adj[u])
			if (e.second == v)
				return true;
		return false;
	}

	bool hasUnacceptableCycleInNEC(vector< set<int> >& NEC, vector<int>& NEC_mapping, int cur_NEC_index, vector<set<int>>& decomposed_NEC) {
		assert(cur_NEC_index < NEC.size());
		if( NEC[cur_NEC_index].size() < 2 )
			return false;

		//Check cycle
		vector< vector<bool> > need_decompose(getNumVertices());
		for(int i=0; i<getNumVertices(); i++)
			need_decompose[i].resize(getNumVertices(), false);

		bool isNeedDecompose=false;
		for (auto& u: NEC[cur_NEC_index]) {
			//Perform DFS
			vector<int> stack;
			vector<bool> visited(getNumVertices(), false);
			for (int v = 0; v < getNumVertices(); v++)
				if (NEC[ NEC_mapping[v] ].size() == 1)
					visited[v]=true;
			stack.push_back(u);
			while(!stack.empty()) {
				int cur_u=stack.back();
				stack.pop_back();
				if(cur_u != u && NEC_mapping[cur_u] == cur_NEC_index) { //Find a problematic cycle
					need_decompose[u][cur_u]=true;
					isNeedDecompose=true;
				}
				for(auto& e: adj[cur_u]) {
					auto v=e.second;
					if (!visited[v]) {
						stack.push_back(v);
						visited[v]=true;
					}
				}
				for(auto& e: in_adj[cur_u]) {
					auto v=e.second;
					if (!visited[v]) {
						stack.push_back(v);
						visited[v]=true;
					}
				}
			}
		}
		if(!isNeedDecompose)
			return false;
		
		//Sequential greedy maximal independent set
		decomposed_NEC.clear();
		for(auto& u: NEC[cur_NEC_index]) {
			bool findNEC=false;
			for(int j=0; j< decomposed_NEC.size(); j++) {
				bool canBeAdded=true;
				for(auto& v: decomposed_NEC[j]) {
					if(need_decompose[u][v]) {
						canBeAdded=false;
						break;
					}
				}
				if(canBeAdded) {
					findNEC=true;
					decomposed_NEC[j].insert(u);
					break;
				}
			}
			if(!findNEC) {
				decomposed_NEC.push_back(set<int>({u}));
			}
		}
		return true;
		
	}

	map< std::pair<int, int>, int> getTopology(int v, bool outEdges, bool toSet=false, set<int> U=set<int>()) {
		map< std::pair<int, int>, int> topology;
		assert(v >= 0);
		topology.insert( make_pair( std::pair<int, int>(-1, vorder[v]), 1) );
		if(outEdges) {
			assert(adj.size() > v && v >= 0);
			for(auto u = adj[v].begin(); u != adj[v].end(); u++) {
				if(toSet && U.count(u->second) == 0) continue;
				auto it = topology.find( std::pair<int, int>(u->first, vorder[u->second]) );
				if(it == topology.end()) {
					topology.insert( make_pair( std::pair<int, int>(u->first, vorder[u->second]), 1) );
				}
				else {
					it->second=it->second+1;
				}
			}
		}
		else {
			assert(in_adj.size() > v);
			for(auto u = in_adj[v].begin(); u != in_adj[v].end(); u++) {
				if(toSet && U.count(u->second) == 0) continue;
				auto it = topology.find( std::pair<int, int>(u->first, vorder[u->second]) );
				if(it == topology.end()) {
					topology.insert( make_pair( std::pair<int, int>(u->first, vorder[u->second]), 1) );
				}
				else {
					it->second=it->second+1;
				}
			}
		}
		return topology;
	}

	map< std::pair<int, int>, set<int>> categorize(int v, bool outEdges) {
		map< std::pair<int, int>, set<int>> cat;
		if(outEdges) {
			assert(adj.size() > v && v >= 0);
			for (auto& edge : adj[v]) {
				pair<int, int> key(edge.first, vorder[edge.second]);
				cat[key].insert(edge.second);
			}
		}
		else {
			assert(in_adj.size() > v);
			for (auto& edge : in_adj[v]) {
				pair<int, int> key(edge.first, vorder[edge.second]);
				cat[key].insert(edge.second);
			}
		}
		return cat;
	}

	void deleteEdge(int u, int v, int o) {
		pair<int, int> e1(o, v);
		pair<int, int> e2(o, u);

		auto it = std::lower_bound(adj[u].begin(), adj[u].end(), e1);
		adj[u].erase(it);
		
		//fixed, originally searched for e1
		for(auto it = in_adj[v].begin(); it != in_adj[v].end(); it++) {
			if(*it == e2) {
				in_adj[v].erase(it);
				break;
			}
		}
	}

	size_t getNumEdges() {
		size_t ret = 0;
		for (auto& a : adj)
			ret += a.size();
		return ret;
	}

	void DFS(int u, vector<bool>& visited) {
		visited[u] = true;
		for (auto& e : adj[u]) {
			int v = e.second;
			if (!visited[v]) {
				DFS(v, visited);
			}
		}
		for (auto& e : in_adj[u]) {
			int v = e.second;
			if (!visited[v]) {
				DFS(v, visited);
			}
		}
	}

	bool DFSUtil(int u, vector<bool>& visited, vector<int>& deg) {
		visited[u] = true;
		for (auto& e : adj[u]) {
			int v = e.second;
			if (!visited[v]) {
				if (deg[u] < 2) {
					deg[v]--;
				}
				if (DFSUtil(v, visited, deg)) {
					deg[u]--;
				}
			}
		}
		for (auto& e : in_adj[u]) {
			int v = e.second;
			if (!visited[v]) {
				if (deg[u] < 2) {
					deg[v]--;
				}
				if (DFSUtil(v, visited, deg)) {
					deg[u]--;
				}
			}
		}
		return deg[u] < 2;
	}

	size_t getNumConEdges() {
		vector<bool> visited(getNumVertices(), false);
		vector<int>  deg(getNumVertices());

		int min_deg = MAX_INT;
		int s;

		for (int i = 0; i < getNumVertices(); i++) {
			deg[i] = getDegree(i);
			if (deg[i] < min_deg) {
				min_deg = deg[i];
				s = i;
			}
		}

		DFSUtil(s, visited, deg);

		size_t ret = 0;
		for (int u = 0; u < getNumVertices(); u++) {
			if (getDegree(u) == 1)
				ret++;
			else if (deg[u] >= 2) {
				for (auto& e : adj[u]) {
					if (deg[e.second] >= 2)
						ret++;
				}
			}
		}
		return ret;
	}

	//{(src, dst, order)}
	void getCoreEdges(vector<Edge>& edges) {
		vector<bool> visited(getNumVertices(), false);
		vector<int>  deg(getNumVertices());

		int min_deg = MAX_INT;
		int s;

		for (int i = 0; i < getNumVertices(); i++) {
			deg[i] = getDegree(i);
			if (deg[i] < min_deg) {
				min_deg = deg[i];
				s = i;
			}
		}

		DFSUtil(s, visited, deg);

		for (int u = 0; u < getNumVertices(); u++) {
			if (deg[u] >= 2) {
				for (auto& e : adj[u]) {
					if (deg[e.second] >= 2) {
						edges.emplace_back(u, e.second, e.first);
					}
				}
			}
		}
	}
	
	void getCoreLeafEdges(vector<Edge>& edges) {
		getCoreEdges(edges);

		for (int u = 0; u < getNumVertices(); u++) {
			if (getDegree(u) == 1) {
				if (adj[u].size())
					edges.emplace_back(u, adj[u][0].second, adj[u][0].first);
				else
					edges.emplace_back(in_adj[u][0].second, u, in_adj[u][0].first);
			}
		}
	}

	void getCoreDegree(vector<int>& deg) {
		vector<bool> visited(getNumVertices(), false);
		deg.resize(getNumVertices());

		int min_deg = MAX_INT;
		int s;

		for (int i = 0; i < getNumVertices(); i++) {
			deg[i] = getDegree(i);
			if (deg[i] < min_deg) {
				min_deg = deg[i];
				s = i;
			}
		}

		DFSUtil(s, visited, deg);
	}

	bool isConnected(int u, int v) {
		vector<bool> visited(getNumVertices(), false);
		DFS(u, visited);
		return visited[v];
	}

	void findAPs(vector<int>& matching_order, vector<bool>& ap, vector<vector<int>>& children) {
		assert(ap.size() == getNumVertices());
		vector<bool> visited(ap.size(), false);
		vector<int> parent(ap.size(), -1);
		vector<int> disc(ap.size());
		vector<int> low(ap.size());

		const int start = matching_order[0];
		APUtil(start, start, visited, disc, low, parent, ap, children);
		for(auto& e: adj[start])
			children[start].push_back(e.second);
		for(auto& e: in_adj[start])
			children[start].push_back(e.second);

		if (getNumVertices() == getNumEdges() + 1)
			for (int u = 0; u < getNumVertices(); u++) {
				if (getDegree(u) > 1)
					assert(ap[u]);
				else
					assert(!ap[u]);
			}

		for (int u = 0; u < getNumVertices(); u++)
			visited[u] = false;
		for (int u = 0; u < getNumVertices(); u++) {
			for (int c : children[u])
				visited[c] = true;
		}
		for (int i = 0; i < getNumVertices(); i++) {
			for (int c : children[ matching_order[i] ])
				visited[c] = false;
			children[matching_order[i]].clear();
			DFSAPUtil(matching_order[i], visited, children[matching_order[i]]);
		}

		for (int i = getNumVertices()-1; i > 0; i--) {
			int u = matching_order[i];
			//not ap, not root 
			if (!ap[u]) {
				assert(parent[u] != -1);
				children[parent[u]].insert(children[parent[u]].end(), children[u].begin(), children[u].end());
				children[u].clear();
			}
		}
	}

	void DFSAPUtil(const int u, vector<bool>& visited, vector<int>& children) {
		visited[u] = true;
		for (auto& e : adj[u]) {
			int v = e.second;
			if (!visited[v]) {
				children.push_back(v);
				DFSAPUtil(v, visited, children);
			}
		}
		for (auto& e : in_adj[u]) {
			int v = e.second;
			if (!visited[v]) {
				children.push_back(v);
				DFSAPUtil(v, visited, children);
			}
		}
	}

	void APUtil(const int start, const int u, vector<bool>& visited, vector<int>& disc, vector<int>& low, vector<int>& parent, vector<bool>& ap, vector<vector<int>>& children) {
		static int time = 0;
		int num_children = 0;
		visited[u] = true;
		disc[u] = low[u] = ++time;
		for (auto& e : adj[u]) {
			int v = e.second;
			if (!visited[v]) {
				num_children++;
				parent[v] = u;
				APUtil(start, v, visited, disc, low, parent, ap, children);

				//check if the subtree rooted at v has a connection to one of the ancestors of u
				low[u] = min(low[u], low[v]);
				//if u is loot and has two or more children
				if (parent[u] == -1 && num_children > 1)
					ap[u] = true;
				//if u is not root & low value of one of its child is more than discovery time of u
				if (parent[u] != -1 && low[v] >= disc[u]) {
					ap[u] = true;
					children[u].push_back(v);
				}
			}
			else if (v != parent[u]) {
				if(parent[u] != -1 && low[v] >= disc[u])
					children[u].push_back(v);
				low[u] = min(low[u], disc[v]);
			}
		}
		for (auto& e : in_adj[u]) {
			int v = e.second;
			if (!visited[v]) {
				num_children++;
				parent[v] = u;
				APUtil(start, v, visited, disc, low, parent, ap, children);

				//check if the subtree rooted at v has a connection to one of the ancestors of u
				low[u] = min(low[u], low[v]);
				//if u is loot and has two or more children
				if (parent[u] == -1 && num_children > 1)
					ap[u] = true;
				//if u is not root & low value of one of its child is more than discovery time of u
				if (parent[u] != -1 && low[v] >= disc[u]) {
					ap[u] = true;
					children[u].push_back(v);
				}
			}
			else if (v != parent[u]) {
				if(parent[u] != -1 && low[v] >= disc[u])
					children[u].push_back(v);
				low[u] = min(low[u], disc[v]);
			}
		}

		if (start == u)
			time = 0; //reset to prevent overflow
	}

	bool good() {
		assert(directed);
		assert(adj.size() == in_adj.size());
		assert(vorder.size() == adj.size());
		for (int u = 0; u < adj.size(); u++) {
			for (int ei = 0; ei < adj[u].size(); ei++) {
				if (ei > 0)
					if (adj[u][ei] <= adj[u][ei-1])
						return false;
				int v = adj[u][ei].second;
				bool found = false;
				for (auto& e : in_adj[v])
					if (e.second == u)
						found = true;
				if (!found)
					return false;
			}
		}
		int num_out_edge = 0, num_in_edge = 0;
		for (int u = 0; u < adj.size(); u++)
			num_out_edge += adj[u].size(); 
		for (int u = 0; u < in_adj.size(); u++)
			num_in_edge += in_adj[u].size();
		return num_out_edge == num_in_edge;
	}

	bool goodMapping(VID_MAP& mapping) {
		for(int v = 0; v < adj.size(); v++) {
			if (std::find(mapping.begin(), mapping.end(), v) == mapping.end())
				return false;
		}
		return true;
	}

	bool hasIncomingEdge(int a){
		assert(in_adj.size() > a);
		return in_adj[a].size() > 0;
	}

	size_t getNumVertices() const {
		assert(in_adj.size() == adj.size());
		return adj.size();
	}

	size_t getDegree(int u) {
		assert(in_adj.size() == adj.size());
		assert(in_adj.size() > u && adj.size() > u);
		return adj[u].size() + in_adj[u].size(); 
	}
	
	void addVertex(int vo=PATTERN_EMPTY_VLABEL) {
		assert(in_adj.size() == adj.size() && adj.size() == vorder.size());
		adj.resize(adj.size()+1);
		in_adj.resize(in_adj.size()+1);
		adj[adj.size()-1].clear();
		in_adj[in_adj.size()-1].clear();
		vorder.push_back(vo);
	}

	void addEdge(int u, int v, int o) {
		pair<int, int> e(o, v);
		assert(in_adj.size() == adj.size());
		assert(adj.size() > u);

		auto it = std::lower_bound(adj[u].begin(), adj[u].end(), e);
		adj[u].insert(it, e);
		
		assert(in_adj.size() > v);
		in_adj[v].emplace_back(o, u);
	}

	bool hasEdge(int u, int v) {
		for (auto& e : adj[u])
			if (e.second == v)
				return true;
		assert(in_adj.size() > u);
		for (auto& e : in_adj[u])
			if (e.second == v)
				return true;
		return false;
	}

	//label, dir
	pair<int, bool> getEdge(int u, int v) {
		for (auto& e : adj[u])
			if (e.second == v)
				return make_pair(e.first, true);
		for (auto& e : in_adj[u])
			if (e.second == v)
				return make_pair(e.first, false);
		assert(false);
	}

	void sort() {
		for (int v = 0; v < adj.size(); v++) {
			std::sort(adj[v].begin(), adj[v].end());
		}
	}

	void deleteVertex(int d) {
		for (int u = 0; u < getNumVertices(); u++) {
			for (auto& e : adj[u]) {
				if (e.second > d)
					e.second--;
			}
			for (auto& e : in_adj[u]) {
				if (e.second > d)
					e.second--;
			}
		}
		adj.erase(adj.begin() + d);
		in_adj.erase(in_adj.begin() + d);
		vorder.erase(vorder.begin() + d);
	}

	void deleteVertexWithoutRenumber(int d) {
		adj[d].clear();
		in_adj[d].clear();
		for (int u = 0; u < getNumVertices(); u++) {
			for (int i = 0; i < adj[u].size(); i++) {
				if (adj[u][i].second == d) {
					adj[u].erase(adj[u].begin() + i);
					i--;
				}
			}
			for (int i = 0; i < in_adj[u].size(); i++) {
				if (in_adj[u][i].second == d) {
					in_adj[u].erase(in_adj[u].begin() + i);
					i--;
				}
			}
		}
	}

	void reorder(VID_MAP& mapping) {
		assert(in_adj.size() == adj.size());
		assert(adj.size() == mapping.size());
		assert(goodMapping(mapping));

		vector<int> new_vorder(vorder.size());
		for (int v = 0; v < vorder.size(); v++) {
			new_vorder[ mapping[v] ] = vorder[v];
		}
		vorder = new_vorder;


		vector<vector<pair<int, int>>> new_adj(adj.size());

		for (int v = 0; v < adj.size(); v++) {
			new_adj[ mapping[v] ].clear();
			for (int ei = 0; ei < adj[v].size(); ei++) {
				new_adj[ mapping[v] ].emplace_back( adj[v][ei].first, mapping[adj[v][ei].second] );
			}
		}
		adj = new_adj;
		sort();

		vector<vector<pair<int, int>>> new_in_adj(in_adj.size());
		for (int v = 0; v < in_adj.size(); v++) {
			new_in_adj[ mapping[v] ].clear();
			for (int ei = 0; ei < in_adj[v].size(); ei++) {
				new_in_adj[ mapping[v] ].emplace_back( in_adj[v][ei].first, mapping[in_adj[v][ei].second] );
			}
		}
		in_adj = new_in_adj;
	}

	bool operator==(Pattern& other) {
		if (getNumVertices() != other.getNumVertices())
			return false;
		if (vorder.size() != other.vorder.size())
			return false;
		for (int u = 0; u < getNumVertices(); u++) {
			if (adj[u].size() != other.adj[u].size())
				return false;
			if (in_adj[u].size() != other.in_adj[u].size())
				return false;
			if (vorder[u] != other.vorder[u])
				return false;
			for (int ei = 0; ei < adj[u].size(); ei++) {
				if (adj[u][ei] != other.adj[u][ei])
					return false;
			}
		}
		return true;
	}
	
	friend ostream& operator<< (ostream&, Pattern&);
};

string encode_el(Pattern&, VID_MAP&);
string encode_vl(Pattern&, VID_MAP&);
string encode_evl(Pattern&, VID_MAP&);

#endif
