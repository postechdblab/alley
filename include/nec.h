#ifndef NEC_H_
#define NEC_H_

#include "pattern.h"

struct NEC_element {
	int label;
	int parent_id;
	int represent_node;

	NEC_element() {

	}

	NEC_element(int label, int parent_id, int represent_node) {
		this->label = label;
		this->parent_id = parent_id;
		this->represent_node = represent_node;
	}
};

struct NEC_set_array_element {
	int parent_id;
	int represent_node;
	int sum;
	NEC_set_array_element() {
	}
	NEC_set_array_element(int parent_id, int represent_node, int sum) {
		this->parent_id = parent_id;
		this->represent_node = represent_node;
		this->sum = sum;
	}
};

inline bool sort_by_NEC_label(NEC_element p1, NEC_element p2) {
	return (p1.label < p2.label);
}

const int MAX_QUERY_NODE = 400;

class NECMapper {
private:
	//don't need to set outside
	int NEC_map [MAX_QUERY_NODE];
	char visited_for_query[MAX_QUERY_NODE];
	vector<pair<int, int> > NEC_set_by_label_index; //include a redunant element to set the end
	int bin_query[MAX_QUERY_NODE];
	int pos_query[MAX_QUERY_NODE];
	int vert_query[MAX_QUERY_NODE];

	//need to set outside
	int cnt_node_query; //# query vertics
	int MAX_DEGREE_QUERY = 0; //max degree
	int core_number_query[MAX_QUERY_NODE]; //initially set to degrees
	int node_degree_query[MAX_QUERY_NODE]; //degrees
	int nodes_label_query[MAX_QUERY_NODE]; //labels 

	int query_nodes_array_info[MAX_QUERY_NODE + 1];
	int dfs_stack_query[MAX_QUERY_NODE];
	int residual_tree_match_seq[MAX_QUERY_NODE];
	pair <int,double> residual_tree_leaf_node [MAX_QUERY_NODE];
	vector <int> query_nodes_array;
	int * tree_node_parent;

	int cnt_unique_label;
	int * NEC_mapping;
	int * NEC_mapping_Actual;
	NEC_element * NEC_mapping_pair;
	NEC_set_array_element * NEC_set_array;

	NECMapper(DataGraph* g) {
		cnt_unique_label = g->GetNumVLabels();
		NEC_mapping = new int[(cnt_unique_label + 1) * MAX_QUERY_NODE]; //label and the parent node id
		memset(NEC_mapping, 0, sizeof(int) * (cnt_unique_label + 1) * MAX_QUERY_NODE);

		NEC_mapping_pair = new NEC_element[(cnt_unique_label + 1) * MAX_QUERY_NODE]; // <label, parent node id>
		NEC_set_array = new NEC_set_array_element[(cnt_unique_label + 1) * MAX_QUERY_NODE]; // <parent_id, sum>
	}

	~NECMapper() {
		delete [] NEC_mapping;
		delete [] NEC_mapping_pair;
		delete [] NEC_set_array;
	}

public:
	inline bool makeNEC(Pattern& p, VID_MAP& nec_map) {
		cnt_node_query = p.getNumVertices();
		MAX_DEGREE_QUERY = 0;
		for (int i = 0; i < cnt_node_query; i++) {
			int deg = p.getDegree(i);
			core_number_query[i] = node_degree_query[i] = deg;
			nodes_label_query[i] = p.vorder[i];
			if (deg > MAX_DEGREE_QUERY)
				MAX_DEGREE_QUERY = deg;
		}
		coreDecomposition_query(p);
		extractResidualStructures();
	}
private:

	inline void coreDecomposition_query(Pattern& p){

		//begin starting the core-decomposition, core number is the degree number
		int * bin = bin_query;	//	int bin [MAX_DEGREE_QUERY + 1];
		int * pos = pos_query;	//	int pos [cnt_node_query];
		int * vert = vert_query;//	int vert [cnt_node_query];

		memset(bin, 0, sizeof(int) * ( MAX_DEGREE_QUERY + 1) );

		for (int i = 0; i < cnt_node_query; i ++)
			bin[ core_number_query [i] ] ++;

		int start = 0;
		int num;

		for (int d = 0; d <= MAX_DEGREE_QUERY; d++){
			num = bin[d];
			bin[d] = start;
			start += num;
		}

		for (int i = 0; i < cnt_node_query; i++){
			pos[i] = bin[ core_number_query[i] ];
			vert[ pos[i] ] = i;
			bin[ core_number_query[i] ] ++;
		}

		for (int d = MAX_DEGREE_QUERY; d > 0; d --)
			bin[d] = bin[d-1];
		bin[0] = 0;

		for (int i = 0; i < cnt_node_query; i++){

			int v = vert[i];

			for (auto& e : p.adj[v]) {
				int u = e.second;

				if (core_number_query[u] > core_number_query[v]){

					int du = core_number_query[u];
					int pu = pos[u];

					int pw = bin[du];
					int w = vert[pw];

					if (u != w){	//if not the same node, switch the position of the two nodes.
						pos[u] = pw;
						pos[w] = pu;
						vert[pu] = w;
						vert[pw] = u;
					}

					bin[du] ++;
					core_number_query[u]--;
				}
			}
			for (auto& e : p.in_adj[v]) {
				int u = e.second;

				if (core_number_query[u] > core_number_query[v]){

					int du = core_number_query[u];
					int pu = pos[u];

					int pw = bin[du];
					int w = vert[pw];

					if (u != w){	//if not the same node, switch the position of the two nodes.
						pos[u] = pw;
						pos[w] = pu;
						vert[pu] = w;
						vert[pw] = u;
					}

					bin[du] ++;
					core_number_query[u]--;
				}
			}
		}
	}


	inline void extractResidualStructures(Pattern& p){
		memset(NEC_mapping, 0, sizeof(int) * (cnt_unique_label + 1) * MAX_QUERY_NODE);

		int residual_tree_match_seq_index = 0;
		int residual_tree_leaf_node_index = 0;
		int NEC_mapping_pair_index = 0;

		memset(NEC_map, -1, sizeof(int) * cnt_node_query );

		char * visited = visited_for_query; //indicate a node is visited or not
		memset(visited, 0, sizeof(char) * cnt_node_query);

		for (int i = 0; i < cnt_node_query; i++) {//for each node in the query

			if (core_number_query[i] < 2)//not in the two-core
				continue;

			//now i must be a 2-core node => next, we check whether i is a articulation node

			//for all of node i's children
			for (auto& e : p.adj[i]) { 
				int child = e.second; 

				if (core_number_query[child] < 2){ // the child node is not in the 2-core

					//two cases here, the NEC node or a residual tree

					if (node_degree_query[child] == 1){ //degree is one ==> NEC node
						//============ CASE ONE: ONE-DEGREE NODES => NEC nodes =====================

						int label = nodes_label_query[child]; 

						if (NEC_mapping[label * MAX_QUERY_NODE + i] == 0) {
							NEC_mapping_pair[NEC_mapping_pair_index ++] = NEC_element(label, i, child);// child is the representative node
							NEC_map [child] = child;//NEC map
							NEC_mapping_Actual[label * MAX_QUERY_NODE + i] = child;

#ifdef RESULT_ENUMERATION
							NEC_Node_array[child].node = child;
							NEC_Node_array[child].nextAddress = NULL;
#endif

						} else {
							NEC_map [child] = NEC_mapping_Actual[label * MAX_QUERY_NODE + i];//NEC map

#ifdef RESULT_ENUMERATION
							int rep =NEC_mapping_Actual[label * MAX_QUERY_NODE + i];
							NEC_Node_array[child].node = child;
							NEC_Node_array[child].nextAddress = NEC_Node_array[ rep ].nextAddress;
							NEC_Node_array[ rep ].nextAddress = &NEC_Node_array[child];
#endif
						}

						NEC_mapping[label * MAX_QUERY_NODE + i] ++; // the label with parent being i, nec_count ++

					} else {
						//============ CASE TWO: NORMAL CASE, THE QUERY TREE ================
						// extract the query tree for extra region candidate extraction, based on DFS
						// also give a DFS-based query sequence at the same time

						int * dfs_stack = dfs_stack_query;
						int dfs_stack_index = 0;

						visited[i] = 1; // this is the start node's parent node (a marked node)
						visited[child] = 1; // this is the start node

						dfs_stack[dfs_stack_index ++] = child;
						residual_tree_match_seq[residual_tree_match_seq_index ++] = child;

						tree_node_parent[child] = i;

						while (dfs_stack_index != 0) {

							int current_node = dfs_stack[dfs_stack_index - 1];
							dfs_stack_index--;

							int added_child = 0;

							for (int m = query_nodes_array_info[current_node]; m < query_nodes_array_info[current_node + 1]; m ++){

								int child_node = query_nodes_array[m];

								if (!visited[child_node]) {

									visited[child_node] = 1;

									//======== special treatment here: if a node is a leaf (degree being 1), then put it into nec node set
									if (node_degree_query[child_node] == 1){

										int label = nodes_label_query[child_node];

										if (NEC_mapping[label * MAX_QUERY_NODE + current_node] == 0) {
											NEC_mapping_pair[NEC_mapping_pair_index ++] = NEC_element(label, current_node, child_node);// child is the repesentive node
											NEC_map [child_node] = child_node;//NEC map
											NEC_mapping_Actual[label * MAX_QUERY_NODE + current_node] = child_node;
#ifdef RESULT_ENUMERATION
											NEC_Node_array[child_node].node = child_node;
											NEC_Node_array[child_node].nextAddress = NULL;
#endif
										}
										else{
											NEC_map [child_node] = NEC_mapping_Actual[label * MAX_QUERY_NODE + current_node];//NEC map
#ifdef RESULT_ENUMERATION
											int rep = NEC_mapping_Actual[label * MAX_QUERY_NODE + current_node];
											NEC_Node_array[child_node].node = child_node;
											NEC_Node_array[child_node].nextAddress = NEC_Node_array[ rep ].nextAddress;
											NEC_Node_array[ rep ].nextAddress = &NEC_Node_array[child_node];
#endif
										}
										NEC_mapping[label * MAX_QUERY_NODE + current_node]++; // the label with parent being i, nec_count ++
										continue;
									}
									//===========================================================
									tree_node_parent[child_node] = current_node;
									added_child ++;
									dfs_stack[dfs_stack_index ++] = child_node;
									residual_tree_match_seq[residual_tree_match_seq_index ++] = child_node;
								}

								if (added_child == 0)//this information is recorded for extracting the matching sequence for the tree matching sequence.
									residual_tree_leaf_node[residual_tree_leaf_node_index ++] = make_pair(current_node, 0);
							}
						}
					}
				}
			}
		}


		//================ construct the NEC set by label: each label is with a vector which contains many NECs with this label.=========
		sort(NEC_mapping_pair, NEC_mapping_pair + NEC_mapping_pair_index, sort_by_NEC_label);
		int last_label;
		int NEC_set_index = 0;
		NEC_set_by_label_index.clear();
		int sum;
		if (NEC_mapping_pair_index == 1){
			NEC_element & nec_ele = NEC_mapping_pair[0];
			int label = nec_ele.label;
			int parent_id = nec_ele.parent_id;
			int represent_child = nec_ele.represent_node;
			sum = NEC_mapping[label * MAX_QUERY_NODE + parent_id];
			NEC_mapping[label * MAX_QUERY_NODE + parent_id] = 0; //reset it back to 0
			NEC_set_by_label_index.push_back(make_pair(label, NEC_set_index));
			NEC_set_array[NEC_set_index ++] = NEC_set_array_element(parent_id, represent_child, sum);
			NEC_set_by_label_index.push_back(make_pair(-1, NEC_mapping_pair_index)); // redundant element to set the end
		} else {
			for (int i = 0; i < NEC_mapping_pair_index; i++) {

				NEC_element & nec_ele = NEC_mapping_pair[i];

				int label = nec_ele.label;
				int parent_id = nec_ele.parent_id;
				int represent_child = nec_ele.represent_node;
				sum = NEC_mapping[label * MAX_QUERY_NODE + parent_id];
				NEC_mapping[label * MAX_QUERY_NODE + parent_id] = 0; //reset it back to 0

				if (i == 0) {
					NEC_set_by_label_index.push_back(make_pair(label, NEC_set_index));
					NEC_set_array[NEC_set_index ++] = NEC_set_array_element(parent_id, represent_child, sum);
					last_label = label;
					continue;
				} else if (i == NEC_mapping_pair_index - 1) {
					if (label != last_label)
						NEC_set_by_label_index.push_back(make_pair(label, NEC_set_index));
					NEC_set_array[NEC_set_index ++] = NEC_set_array_element(parent_id, represent_child, sum);
					NEC_set_by_label_index.push_back(make_pair(-1, NEC_mapping_pair_index)); // redunant element to set the end
					continue;
				}

				if (label != last_label) {
					NEC_set_by_label_index.push_back(make_pair(label, NEC_set_index));
					last_label = label;
				}

				NEC_set_array[NEC_set_index ++] = NEC_set_array_element(parent_id, represent_child, sum);
			}
		}

#ifdef RESULT_ENUMERATION
		for (int i = 0; i < cnt_node_query; i++){
			if (node_degree_query[i] != 1)
				continue;
			NEC_Node * next = NEC_Node_array[i].nextAddress;
			if (next == NULL)
				continue;
			while (next != NULL)
				next = next->nextAddress;
		}

#endif


#ifdef	OUTPUT_EXTRA_INFO

		int sum_node = 0;

		if (NEC_mapping_pair_index != 0){
			for (int i = 0; i < NEC_set_by_label_index.size() - 1; i++) {
				int label = NEC_set_by_label_index[i].first;
				int start = NEC_set_by_label_index[i].second;
				int end = NEC_set_by_label_index[i + 1].second;

				for (int j = start; j < end; j++) {
					int parent_id = NEC_set_array[j].parent_id;
					int sum = NEC_set_array[j].sum;
					sum_node += sum;
					cerr << "label :" << label << " => parent id " << parent_id << " \t sum => " << sum
						<< "\t representative node is " << NEC_set_array[j].represent_node<< endl;
				}
			}
		}

		cerr << "NEC classes contained: " << NEC_mapping_pair_index << " classes with " << sum_node << " nodes " << endl;
		cerr << "Query trees with sum node: " << residual_tree_match_seq_index
			<< " and tree leaf index is " << residual_tree_leaf_node_index << endl;
#endif
	}
};
#endif
