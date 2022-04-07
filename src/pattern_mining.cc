#include <algorithm>
#include <queue>
#include "../include/pattern_mining.h"
#include "../include/intersection.h"

bool calculate_domains_context = false;

void MinerWrapper::BuildIndex(int maxL, double sampling_ratio, double failure_threshold, int num_bin, label_type lt) {
	((Alley*) estimator)->SetMode(true, maxL);
	if ( (lt == EDGE_LABEL && g->GetNumELabels() < num_bin/2+1) || (lt == VERTEX_LABEL && g->GetNumVLabels() < num_bin/2+1) ) {
		//NOTE THAT vertex label is only used when # ELabels = 1
		miners_.emplace_back(g, estimator, ALL_MINER, lt);
		miners_[0].Initialize(maxL, sampling_ratio, failure_threshold, true);
		miners_[0].GroupLabels(num_bin);
		if (lt == EDGE_LABEL)
			g->GroupELabels(miners_[0].getLabel2Order(), miners_[0].getOrder2Label());
		else if (lt == VERTEX_LABEL)
			g->GroupVLabels(miners_[0].getLabel2Order(), miners_[0].getOrder2Label());
		miners_[0].BuildIndex();
		miners_[0].Clear();
		cout << "built all index" << endl;
	}
	else {
		int initial_num_bin=num_bin;
		if(lt == EDGE_LABEL)
			initial_num_bin=min(g->GetNumELabels(), num_bin);
		else if(lt == VERTEX_LABEL)
			initial_num_bin=min(g->GetNumVLabels(), num_bin);
		else
			assert(false);

		miners_.emplace_back(g, estimator, PATH_MINER, lt);
		miners_.emplace_back(g, estimator, TREE_MINER, lt);
		miners_.emplace_back(g, estimator, GRAPH_MINER, lt);

		miners_[0].Initialize(maxL-1, sampling_ratio, failure_threshold, true);
		miners_[0].GroupLabels(initial_num_bin);
		if (lt == EDGE_LABEL)
			g->GroupELabels(miners_[0].getLabel2Order(), miners_[0].getOrder2Label());
		else if (lt == VERTEX_LABEL)
			g->GroupVLabels(miners_[0].getLabel2Order(), miners_[0].getOrder2Label());
		miners_[0].BuildIndex();
		miners_[0].Clear();
		cout << "built path index" << endl;

		PatternMining path_to_tree(g, NULL, PATH_MINER, lt);
		path_to_tree.Initialize(maxL, 0, 0, true);
		path_to_tree.GroupLabels(num_bin/2);
		path_to_tree.CopyIndex(miners_[0]);
        //assert(path_to_tree.unionTest());
		path_to_tree.UnionIndex();
		path_to_tree.Clear();
		cout << "copied path to tree" << endl;


		miners_[1].Initialize(maxL-1, sampling_ratio, failure_threshold, true);
		miners_[1].CopyOrder(path_to_tree);
		miners_[1].CopyKeys(path_to_tree);
		if (lt == EDGE_LABEL)
			g->GroupELabels(miners_[1].getLabel2Order(), miners_[1].getOrder2Label());
		else if (lt == VERTEX_LABEL)
			g->GroupVLabels(miners_[1].getLabel2Order(), miners_[1].getOrder2Label());
		miners_[1].BuildIndex(&path_to_tree);
		miners_[1].Clear();
		cout << "built tree index" << endl;

		PatternMining tree_to_graph(g, NULL, TREE_MINER, lt);
		tree_to_graph.Initialize(maxL, 0, 0, true);
		tree_to_graph.GroupLabels(num_bin/4);
		tree_to_graph.CopyIndex(path_to_tree);
		tree_to_graph.CopyIndex(miners_[1]);
		tree_to_graph.UnionIndex();
		tree_to_graph.Clear();
		cout << "copied tree to graph" << endl;

		miners_[2].Initialize(maxL, sampling_ratio, failure_threshold, true); 
		miners_[2].CopyOrder(tree_to_graph);
		miners_[2].CopyKeys(tree_to_graph);
		if (lt == EDGE_LABEL)
			g->GroupELabels(miners_[2].getLabel2Order(), miners_[2].getOrder2Label());
		else if (lt == VERTEX_LABEL)
			g->GroupVLabels(miners_[2].getLabel2Order(), miners_[2].getOrder2Label());
		miners_[2].BuildIndex(&tree_to_graph);
		miners_[2].Clear();
		cout << "built graph index" << endl;
	}
}

void MinerWrapper::Write(const char* fn) {
	FILE* tp = fopen(fn, "w");
	fprintf(tp, "%zu", miners_.size()); 
	for (auto& miner : miners_) {
		char type = miner.GetTypeInChar();
		char ltype = miner.GetLTypeInChar();
		fprintf(tp, "%c%c", type, ltype); 
		miner.Write((string(fn) + "." + type).c_str()); 
		//cout << "tangled count " << miner.review_cnt << endl;
		//cout << "tangled size " << miner.review_sum << endl;
		//cout << "tangled count " << miner.estimator_->review_cnt << endl;
		//cout << "tangled size " << miner.estimator_->review_sum << endl;
	}
	fclose(tp);

#ifdef PRINT_INTERSECTION
	cout << "S edge case 1, 2 count" << endl;
	cout << totalEdgeCase1 << ", " << totalEdgeCase2 << endl;
	cout << "S galloping stat: count, avg. |A|, avg. |B|, avg. |C|" << endl;
	cout << totalGalloping << ", " << (double)totalLenRare/totalGalloping << ", " << (double)totalLenFreq/totalGalloping << ", " << (double)totalCountGalloping/totalGalloping << endl;
	cout << "S shuffling stat: count, avg. |A|, avg. |B|, avg. |C|" << endl;
	cout << totalShuffling << ", " << (double)totalA/totalShuffling << ", " << (double)totalB/totalShuffling << ", " << (double)totalCountShuffling/totalShuffling << endl;

	cout << "CD edge case 1, 2 count" << endl;
	cout << CD_totalEdgeCase1 << ", " << CD_totalEdgeCase2 << endl;
	cout << "CD galloping stat: count, avg. |A|, avg. |B|, avg. |C|" << endl;
	cout << CD_totalGalloping << ", " << (double)CD_totalLenRare/CD_totalGalloping << ", " << (double)CD_totalLenFreq/CD_totalGalloping << ", " << (double)CD_totalCountGalloping/CD_totalGalloping << endl;
	cout << "CD shuffling stat: count, avg. |A|, avg. |B|, avg. |C|" << endl;
	cout << CD_totalShuffling << ", " << (double)CD_totalA/CD_totalShuffling << ", " << (double)CD_totalB/CD_totalShuffling << ", " << (double)CD_totalCountShuffling/CD_totalShuffling << endl;
#endif

#ifdef PRINT_SEARCH
	cout << "S search count, avg. |L|" << endl;
	cout << totalSearch << ", " << (double)totalSearchLen/totalSearch << endl;

	cout << "CD search count, avg. |L|" << endl;
	cout << CD_totalSearch << ", " << (double)CD_totalSearchLen/CD_totalSearch << endl;
#endif
}

void MinerWrapper::Read(const char* fn) {
	FILE* tp = fopen(fn, "r");
	size_t num_miners;
	char type, ltype;
	fscanf(tp, "%zu", &num_miners); 
	for (size_t i = 0; i < num_miners; i++) {
		fscanf(tp, "%c%c", &type, &ltype);
		miner_type t;
		label_type lt;
		switch(type) {
			case 'p': t=PATH_MINER; break;
			case 't': t=TREE_MINER; break;
			case 'g': t=GRAPH_MINER; break;
			case 'a': t=ALL_MINER; break;
			default: assert(false);
		}
		switch(ltype) {
			case 'e': lt=EDGE_LABEL; break;
			case 'v': lt=VERTEX_LABEL; break;
			default: assert(false);
		}
		miners_.emplace_back(g, estimator, t, lt);
		miners_.back().Read((string(fn) + "." + type).c_str()); 
	}
	fclose(tp);
}


ostream& operator<< (ostream& os, Pattern& pattern) {
	VID_MAP vid;
    string code = encode_el(pattern, vid);
    string code2 = encode_vl(pattern, vid);
	os << code << "[" << code2 << "]";
	return os;
}

ostream& operator<< (ostream& os, DAGNode& node) {
	os << "pattern: " << node.pattern << endl;
	return os;
}

void PatternMining::CopyKeys(PatternMining& other) {
	for (int L = 2; L <= other.maxL_; L++) {
		for (size_t i = 0; i < other.DAG_[L].size(); i++) {
			auto& original = other.DAG_[L][i];
            auto it = code2index_[L].find(original.code);
            if (it == code2index_[L].end()) {
                code2index_[L][original.code] = MAX_SIZE; 
            }
        }
    }
}

void PatternMining::CopyIndex(PatternMining& other) {
	//populate copied_DAG_ and copied_code2index_
    cout << "other maxL = " << other.maxL_ << endl;

	for (int L = 2; L <= other.maxL_; L++) {
        cout << "L = " << L << ": " << other.DAG_[L].size() << " original nodes" << endl;
        assert(other.DAG_.size() > L);
		for (size_t i = 0; i < other.DAG_[L].size(); i++) {
			auto& original = other.DAG_[L][i];
			Pattern p(original.pattern);
			//convert labels using label2order_ 
            if(ltype_ == EDGE_LABEL)
                p.convertELabels(label2order_);
            else if(ltype_ == VERTEX_LABEL)
                p.convertVLabels(label2order_);
            else
                assert(false);

			VID_MAP vid_map; 
			string code = this->encode(p, vid_map);

			size_t index;
			auto it = code2index_[L].find(code);
			if (it != code2index_[L].end()) {
				//if it already exists in DAG
				index = it->second;
                assert( index >= 0 && index < DAG_[L].size() );
			}
			else {
				code2index_[L][code] = index = DAG_[L].size();
				p.reorder(vid_map); //store canonicalized pattern
				DAG_[L].emplace_back(p);
				DAG_[L].back().code = code;
				DAG_[L].back().mat_info.resize(p.getNumVertices());
                for(auto& info: DAG_[L].back().mat_info)
                    info.clear();
			}

			//maintain connection from copied to original
			auto& copied = DAG_[L][index];
			for (int u = 0; u < p.getNumVertices(); u++) {
				int reu = vid_map[u];
                assert(reu >= 0 && reu < p.getNumVertices());
                assert(original.mat_info.size() > u);
                assert(original.mat_info[u].size() > 0);
                for (const MatInfo& info : original.mat_info[u]) {
                    assert(info.L <= other.maxL_);
                    auto r = info.pm->getMat(info.L, info.index, info.u);
                    assert(r.end - r.begin == info.size);
                }
				copied.mat_info[reu].insert(*(original.mat_info[u].begin()));
			}
		}
		if (L >= maxL_-1)
			break;
	}
}

bool PatternMining::unionTest() {
    auto& infos = DAG_[2][10].mat_info[2];
    for (auto& info : infos) {
        if (info.L == 2 && info.index == 49 && info.size == 3 && info.u == 1) {
            if (info.pm->getMat(2, 49, 1).begin[2] == 99029)
                return true;
        }
    }
    return false;
}

void PatternMining::UnionIndex() {
	for (int L = 2; L <= maxL_; L++) {
        cout << "L = " << L << ": " << DAG_[L].size() << " copied nodes" << endl;
        mat_old_[L].resize((L+1)*DAG_[L].size());
        mat_i_[L].resize((L+1)*DAG_[L].size());
        mat_o_[L].resize((L+1)*DAG_[L].size());
        mat_s_[L].resize((L+1)*DAG_[L].size());
        consideredEdge_[L].resize((L+1)*DAG_[L].size());
		//now iterate over copied_DAG_
		for (size_t i = 0; i < DAG_[L].size(); i++) {
			auto& copied = DAG_[L][i];
			for (int u = 0; u < copied.pattern.getNumVertices(); u++) {
                assert(u < copied.mat_info.size());
                assert(!copied.mat_info[u].empty());
				set<MatInfo> new_info;
				//merge domains (that point to original index)
				for (auto& info : copied.mat_info[u]) {
                    assert(info.L <= maxL_);
                    assert(info.pm != NULL);
                    assert(info.pm != this);
					assert(info.L < info.pm->DAG_.size());
					assert(info.index < info.pm->DAG_[info.L].size());
					auto& parent = info.pm->DAG_[info.L][info.index];
					//convert parent pattern in original index
					//to a pattern in this copied index 
					if (info.L < L) {
						Pattern p(parent.pattern);
                        if(ltype_ == EDGE_LABEL)
                            p.convertELabels(label2order_);
                        else if(ltype_ == VERTEX_LABEL)
                            p.convertVLabels(label2order_);
                        else
                            assert(false);

						VID_MAP vid_map;
						string code = this->encode(p, vid_map);
						//should be copied earlier
						assert(code2index_[info.L].find(code) != code2index_[info.L].end());
						size_t par_index = code2index_[info.L][code];
                        assert(par_index != MAX_SIZE);
                        assert(par_index >= 0 && par_index < DAG_[info.L].size());
						auto& par = DAG_[info.L][par_index];
                        assert(info.u >= 0 && info.u < vid_map.size());
						int par_u = vid_map[info.u];
                        assert(par_u >= 0 && par_u < par.mat_info.size());
						new_info.insert(*par.mat_info[par_u].begin());
					}
					else {
						assert(info.L == L);
                        //assert(info.pm->getMatIndex(info.L, info.index, info.u) == L);
					}
				}
				
				//all original domains are exact
				if (new_info.empty()) {
					//materialize
                    //vector<range> ranges;
                    set<range> ranges;
					for (auto& info : copied.mat_info[u]) {
                        assert(info.pm != this);
						range r = info.pm->getMat(info.L, info.index, info.u);
                        assert(r.end - r.begin == info.size);
						//ranges.push_back(r);
						ranges.insert(r);
					}
					int* arr;
					size_t size;
					//append the union to materialized_
					if (ranges.size() > 1) {
						int total_size = 0;
						for (auto& r : ranges)
							total_size += (r.end - r.begin);
                        if (total_size > g->GetNumVertices())
                            total_size = g->GetNumVertices();
						if (total_size > intersect_capacity_[u]) {
							intersect_capacity_[u] = std::max(total_size, 2 * intersect_capacity_[u]);
							int* b = intersect_[u];
							int* t = temp_[u];
							intersect_[u] = new int[intersect_capacity_[u]];
							temp_[u] = new int[intersect_capacity_[u]];
							delete [] b;
							delete [] t;
						}
						size = do_union(ranges, intersect_[u]); 
						arr = intersect_[u];
                        appendMat(L, i, u, arr, size);
					}
					else {
						size = ranges.begin()->end - ranges.begin()->begin;
						arr = const_cast<int*>(ranges.begin()->begin);
                        pointMat(L, i, u, arr, size);
					}
					//appendMat(L, i, u, arr, size);
                    copied.mat_info[u].clear();
                    copied.mat_info[u].emplace(this, L, i, u, size);
					for (auto& edge : copied.pattern.adj[u]) 
						setConsideredEdge(L, i, u, edge.first, true); 
					for (auto& edge : copied.pattern.in_adj[u]) 
						setConsideredEdge(L, i, u, edge.first, false); 
				}
				else {
					auto info = new_info.begin();
                    assert( info->L >= 0 && info->index >= 0 && info->L < DAG_.size() && info->index < DAG_[info->L].size() );
                    assert( mat_i_[L].size() > ( (info->L +1)*info->index + info->u ) );
                    setMatOld(L, i, u, false);
					setMatIndex(L, i, u, getMatIndex(info->L, info->index, info->u)); 
					setMatOffset(L, i, u, getMatOffset(info->L, info->index, info->u)); 
					setMatSize(L, i, u, getMatSize(info->L, info->index, info->u)); 

					copied.mat_info[u].clear();
					copied.mat_info[u].insert(*info);
					auto& par = DAG_[info->L][info->index].pattern;
					for (auto& edge : par.adj[info->u])
						setConsideredEdge(L, i, u, edge.first, true);
					for (auto& edge : par.in_adj[info->u])
						setConsideredEdge(L, i, u, edge.first, false);
				}
				assert(copied.mat_info[u].size() == 1);
			}
		}
	}
}

void PatternMining::CopyOrder(PatternMining& other) {
	order2label_ = other.order2label_;
	label2order_ = other.label2order_;
}

void PatternMining::Initialize(int maxL, double sampling_ratio, double failure_threshold, bool do_mining) {
    calculate_domains_time_ = 0;

	if (estimator_ != NULL) {
		((Alley*)estimator_)->SetData(g);
		((Alley*)estimator_)->ClearCache();
	}
    
	maxL_ = maxL;
	sampling_ratio_ = sampling_ratio;
    failure_threshold_ = failure_threshold;
    
    do_mining_ = do_mining;
    if (do_mining) {
        intersect_capacity_ = new int[maxL+1];
        intersect_ = new int*[maxL+1];
        temp_      = new int*[maxL+1];

        initial_intersect_capacity_ = new int[maxL+1];
        initial_intersect_ = new int*[maxL+1];
        initial_temp_      = new int*[maxL+1];

        for (int u = 0; u < maxL+1; u++) {
            intersect_capacity_[u] = INIT_INTER_CAPACITY;
            intersect_[u] = new int[INIT_INTER_CAPACITY];
            temp_[u]      = new int[INIT_INTER_CAPACITY];

            initial_intersect_capacity_[u] = INIT_INTER_CAPACITY;
            initial_intersect_[u] = new int[INIT_INTER_CAPACITY];
            initial_temp_[u]      = new int[INIT_INTER_CAPACITY];
        }

        marked_.clear();
        marked_.resize(maxL+1);
        for (int u = 0; u < maxL+1; u++) {
            marked_[u].resize(g->GetNumVertices(), false);
        }
    }

	DAG_.resize(maxL+1);
    for(auto& DAG_node: DAG_)
        DAG_node.clear();
	code2index_.resize(maxL+1);
    for(auto& c2i: code2index_)
        c2i.clear();
	
    mat_old_.resize(maxL+1);
    mat_i_.resize(maxL+1);
    mat_o_.resize(maxL+1);
    mat_s_.resize(maxL+1);
    consideredEdge_.resize(maxL+1);
    materialized_.resize(maxL+1);
}

void PatternMining::Clear() {
    if (do_mining_) {
        delete [] intersect_capacity_;

        delete [] initial_intersect_capacity_;
        for (int i = 0; i <= maxL_; i++) {
            delete [] intersect_[i];
            delete [] temp_[i];

            delete [] initial_intersect_[i];
            delete [] initial_temp_[i];
        }
        marked_.clear();
    }
}

void PatternMining::BuildIndex(PatternMining* other) {
    clock_t start, end;

    start = clock();
    //first, generate size-2 patterns (DAG_[2])
    if (ltype_ == VERTEX_LABEL)
        initPattern(-1);
    for (int i = 0; i < order2label_.size(); i++)
        initPattern(i);
    calculateDomains(2);
    end = clock();
    cout << "indexing time for size 2: " << (double)(end-start)/CLOCKS_PER_SEC << " secs for " << DAG_[2].size() << " nodes" << endl;
	//extend size-L patterns to size-L+1 patterns
	for (int L = 2; L < maxL_; L++) {
        sampling_cnt_ = 0;
        ((Alley*)estimator_)->calculate_domains_cnt_ = 0;
        ((Alley*)estimator_)->calculate_domains_time_ = 0;
        estimator_->initial_sample_space_time_ = 0;
        start = clock();
        int num_extend = 0;
		if (other != NULL) {
			for (size_t i = 0; i < other->DAG_[L].size(); i++) {
				num_extend++;
				extendPattern(L, i, other);
			}
		}
        for (size_t i = 0; i < DAG_[L].size(); i++) {
            num_extend++;
            extendPattern(L, i);
        }
        size_t prev_node_cnt = DAG_[L+1].size();
        cout << "to extend: " << num_extend << ", extended # nodes: " << DAG_[L+1].size() << endl;

		calculateDomains(L+1);

        end = clock();
        double total_time = (double)(end-start)/CLOCKS_PER_SEC;
        cout << "indexing time for size " << (L+1) << ": " << total_time << " secs for " << prev_node_cnt << " nodes pruned to " << DAG_[L+1].size() << endl;
        cout << "sampling for size " << (L+1) << ": " << (total_time - ((Alley*)estimator_)->calculate_domains_time_) << " secs for " << sampling_cnt_ << " nodes" << endl;
        cout << "initial sample space for size " << (L+1) << ": " << estimator_->initial_sample_space_time_ << " secs" << endl; 
        cout << "calculate domains time for size " << (L+1) << ": " << ((Alley*)estimator_)->calculate_domains_time_ << " secs for " << ((Alley*)estimator_)->calculate_domains_cnt_ << " nodes" << endl;
        cout << "tangled count: " << ((Alley*)estimator_)->review_cnt << ", tangled size: " << ((Alley*)estimator_)->review_sum << endl;
	}
}

void PatternMining::GroupVLabels(int num_bin) {
    label2order_.clear();
    order2label_.clear();

    label2order_.resize(g->GetNumVLabels());

    for (int vl = 0; vl < g->GetNumVLabels(); vl++)
        label2order_[vl] = vl;

    sort(label2order_.begin(), label2order_.end(), 
            [&](const int& a, const int& b) -> bool {
            int na = g->GetNumVertices(a); 
            int nb = g->GetNumVertices(b); 
            return na > nb;
            });

    //use all
    int o = 0;
    //equi-height
    int total_vertices = 0;
    for (int vl = 0; vl < g->GetNumVLabels(); vl++) {
        total_vertices += g->GetNumVertices(vl);
    }

    int size = 0;
    for (size_t i = 0; i < label2order_.size(); i++) {
        order2label_.resize(o+1);
        int threshold = total_vertices / num_bin;
        order2label_[o].push_back(label2order_[i]);
        size += g->GetNumVertices(label2order_[i]); 
        if (size > threshold) {
            total_vertices -= size;
            num_bin--;
            o++;
            size = 0;
        }
    }

    for (int i = 0; i < label2order_.size(); i++)
        label2order_[i] = -1;
    for (int i = 0; i < order2label_.size(); i++) {
        for (int l : order2label_[i]) {
            label2order_[l] = i;
        }
    }

    return;
}

void PatternMining::GroupELabels(int num_bin) {
    label2order_.clear();
    order2label_.clear();

    label2order_.resize(g->GetNumELabels());

    for (int el = 0; el < g->GetNumELabels(); el++)
        label2order_[el] = el;

    sort(label2order_.begin(), label2order_.end(), 
            [&](const int& a, const int& b) -> bool {
            int na = g->GetNumEdges(a); 
            int nb = g->GetNumEdges(b); 
            return na > nb;
            });

    //use all
    int o = 0;
    //equi-height
    int total_edges = g->GetNumEdges();
    int size = 0;
    for (size_t i = 0; i < label2order_.size(); i++) {
        order2label_.resize(o+1);
        int threshold = total_edges / num_bin;
        order2label_[o].push_back(label2order_[i]);
        size += g->GetNumEdges(label2order_[i]); 
        if (size > threshold) {
            total_edges -= size;
            num_bin--;
            o++;
            size = 0;
        }
    }
    for (int i = 0; i < label2order_.size(); i++)
        label2order_[i] = -1;
    for (int i = 0; i < order2label_.size(); i++) {
        for (int l : order2label_[i]) {
            label2order_[l] = i;
        }
    }
}

void PatternMining::initPattern(int order) {
	Pattern p(order, bool(ltype_ == EDGE_LABEL));

	if (type_ == GRAPH_MINER || type_ == ALL_MINER) {
#ifdef SELF_LOOP_PATTERN
        for (int i = order; i < (int) order2label_.size(); i++) {
            if (ltype_ == EDGE_LABEL) {
                Pattern p2(p);
                p2.addEdge(1, 0, i);
                addPatternToDAG(p2, 2);
            } else if(ltype_ == VERTEX_LABEL) {
                Pattern p2(p);
                p2.addVertex(i);
                p2.addEdge(0, 1, PATTERN_EMPTY_ELABEL); //
                p2.addEdge(1, 0, PATTERN_EMPTY_ELABEL); //DO NOT USE EDGE LABEL
                addPatternToDAG(p2, 2);
            } else {
                assert(false);
            }
        }
#endif
    }
    if (type_ == PATH_MINER || type_ == ALL_MINER) {
        if(ltype_ == EDGE_LABEL) {
            p.addVertex();
            for (int i = order; i < (int) order2label_.size(); i++) {
                for (int u = 0; u < 2; u++) {
                    {
                        Pattern p2(p);
                        p2.addEdge(u, 2, i);
                        addPatternToDAG(p2, 2);
                    }
                    {
                        Pattern p2(p);
                        p2.addEdge(2, u, i);
                        addPatternToDAG(p2, 2);
                    }
                }
            }
        } else if(ltype_ == VERTEX_LABEL) {
            for (int i = order; i < (int) order2label_.size(); i++) {
                for (int j = i; j < (int) order2label_.size(); j++) {
                    for (int u = 0; u < 3; u++) { //hubnode
                        int v1 = 0, v2= 2;
                        if(u == 0) v1 = 1;
                        if(u == 2) v2 = 1;

                        {
                            Pattern p2(p);
                            p2.addVertex(i);
                            p2.addVertex(j);
                            p2.addEdge(u, v1, PATTERN_EMPTY_ELABEL);
                            p2.addEdge(u, v2, PATTERN_EMPTY_ELABEL);
                            addPatternToDAG(p2, 2);
                        }

                        {
                            Pattern p2(p);
                            p2.addVertex(i);
                            p2.addVertex(j);
                            p2.addEdge(v1, u, PATTERN_EMPTY_ELABEL);
                            p2.addEdge(v2, u, PATTERN_EMPTY_ELABEL);
                            addPatternToDAG(p2, 2);

                        }

                        {
                            Pattern p2(p);
                            p2.addVertex(i);
                            p2.addVertex(j);
                            p2.addEdge(v1, u, PATTERN_EMPTY_ELABEL);
                            p2.addEdge(u, v2, PATTERN_EMPTY_ELABEL);
                            addPatternToDAG(p2, 2);

                        }

                        {
                            Pattern p2(p);
                            p2.addVertex(i);
                            p2.addVertex(j);
                            p2.addEdge(u, v1, PATTERN_EMPTY_ELABEL);
                            p2.addEdge(v2, u, PATTERN_EMPTY_ELABEL);
                            addPatternToDAG(p2, 2);

                        }
                    }
                }
            }
        } else {
            assert(false);
        }
    }
    if (type_ == TREE_MINER) {

    }
}

void PatternMining::addPatternToDAG(Pattern& cur_pattern, int L)
{
    VID_MAP vid_map; 
    string code = this->encode(cur_pattern, vid_map);

    auto it = code2index_[L].find(code);
    if (it != code2index_[L].end()) {
        //if it already exists in DAG
    }
    else {
        code2index_[L][code] = DAG_[L].size();
        cur_pattern.reorder(vid_map); //store canonicalized pattern
        DAG_[L].emplace_back(cur_pattern);
        DAG_[L].back().code = code;
        DAG_[L].back().mat_info.resize(cur_pattern.getNumVertices());
    }
}

//extend the ones that look hard (to reduce # encodes)
//then, even if a pattern is not found in the index, it is not empty
void PatternMining::extendPattern(int curL, size_t index, PatternMining* other) {
	auto& node = other == NULL ? DAG_[curL][index] : other->DAG_[curL][index];
	Pattern p(node.pattern);
	
    int num = p.getNumVertices();

	if (type_ == PATH_MINER) {
		vector<int> leaves;
		if (p.ifPathGetLeaves(leaves)) {
            if(ltype_ == EDGE_LABEL)
                p.addVertex();
			//path
			for (int leaf : leaves) {
                if (ltype_ == EDGE_LABEL) {
                    for (int i = 0; i < order2label_.size(); i++) {
                        {
                            Pattern p2(p);
                            p2.addEdge(leaf, num, i); 
                            addPatternToDAG(node, p2, curL+1);
                        }

                        {
                            Pattern p2(p);
                            p2.addEdge(num, leaf, i); 
                            addPatternToDAG(node, p2, curL+1);
                        }
                    }
                }
                else {
                    for (int i = -1; i < (int) order2label_.size(); i++) {
                        {
                            Pattern p2(p);
                            p2.addVertex(i);
                            p2.addEdge(leaf, num, PATTERN_EMPTY_ELABEL); 
                            addPatternToDAG(node, p2, curL+1);
                        }

                        {
                            Pattern p2(p);
                            p2.addVertex(i);
                            p2.addEdge(num, leaf, PATTERN_EMPTY_ELABEL); 
                            addPatternToDAG(node, p2, curL+1);
                        }
                    }
                }
			}
		}
	}
	else if (type_ == TREE_MINER) {
		vector<int> leaves;
		if (!p.ifPathGetLeaves(leaves)) {
			leaves.clear();
		}
        if (ltype_ == EDGE_LABEL)
            p.addVertex();
		//add to non-leaves
		//to avoid making paths
		for (int u = 0; u < num; u++) {
			if (find(leaves.begin(), leaves.end(), u) != leaves.end())
				continue;
            if (ltype_ == EDGE_LABEL) {
                for (int i = 0; i < order2label_.size(); i++) {
                    {
                        Pattern p2(p);
                        p2.addEdge(u, num, i);
                        addPatternToDAG(node, p2, curL+1);
                    }
                    {

                        Pattern p2(p);
                        p2.addEdge(num, u, i);
                        addPatternToDAG(node, p2, curL+1);
                    }
                }
            }
            else {
                for (int i = -1; i < (int) order2label_.size(); i++) {
                    {
                        Pattern p2(p);
                        p2.addVertex(i);
                        p2.addEdge(u, num, PATTERN_EMPTY_ELABEL); 
                        addPatternToDAG(node, p2, curL+1);
                    }
                    {
                        Pattern p2(p);
                        p2.addVertex(i);
                        p2.addEdge(num, u, PATTERN_EMPTY_ELABEL); 
                        addPatternToDAG(node, p2, curL+1);
                    }
                }
            }
		}
	}
	else if (type_ == GRAPH_MINER || type_ == ALL_MINER) {
		//add a cycle
		for (int u = 0; u < num; u++) {
			for (int v = u + 1; v < num; v++) {
				if (!p.hasDirectedEdge(u, v)) {
                    if(ltype_ == EDGE_LABEL) {
                        for (int i = 0; i < order2label_.size(); i++) {
                            Pattern p2(p);
                            p2.addEdge(u, v, i);
                            addPatternToDAG(node, p2, curL+1);
                        }
                    } else if(ltype_ == VERTEX_LABEL) {
                        Pattern p2(p);
                        p2.addEdge(u, v, PATTERN_EMPTY_ELABEL);
                        addPatternToDAG(node, p2, curL+1);
                    } else {
                        assert(false);
                    }
				}
				if (!p.hasDirectedEdge(v, u)) {
                    if(ltype_ == EDGE_LABEL) {
                        for (int i = 0; i < order2label_.size(); i++) {
                            Pattern p2(p);
                            p2.addEdge(v, u, i);
                            addPatternToDAG(node, p2, curL+1);
                        }
                    } else if(ltype_ == VERTEX_LABEL) {
                        Pattern p2(p);
                        p2.addEdge(v, u, PATTERN_EMPTY_ELABEL);
                        addPatternToDAG(node, p2, curL+1);
                    } else {
                        assert(false);
                    }
				}
			}
		}

		if (!p.isTree() || type_ == ALL_MINER) {
			//add a leaf vertex
            if(ltype_ == EDGE_LABEL)
                p.addVertex();
			for (int u = 0; u < num; u++) {
                if (ltype_ == EDGE_LABEL) {
                    for (int i = 0; i < order2label_.size(); i++) {
                        {
                            Pattern p2(p);
                            p2.addEdge(u, num, i);
                            addPatternToDAG(node, p2, curL+1);
                        }

                        {
                            Pattern p2(p);
                            p2.addEdge(num, u, i);
                            addPatternToDAG(node, p2, curL+1);
                        }
                    }
                }
                else {
                    for (int i = -1; i < (int) order2label_.size(); i++) {
                        {
                            Pattern p2(p);
                            p2.addVertex(i);
                            p2.addEdge(u, num, PATTERN_EMPTY_ELABEL);
                            addPatternToDAG(node, p2, curL+1);
                        }
                        {
                            Pattern p2(p);
                            p2.addVertex(i);
                            p2.addEdge(num, u, PATTERN_EMPTY_ELABEL);
                            addPatternToDAG(node, p2, curL+1);
                        }
                    }
                }
			}
		}
	}
	else
		assert(false);
}

void PatternMining::addPatternToDAG(DAGNode& parent, Pattern& cur_pattern, int L)
{
    VID_MAP vid_map; 
	string code = encode(cur_pattern, vid_map);

    size_t index;
    auto it = code2index_[L].find(code);
    if (it != code2index_[L].end()) {
        //if it already exists in DAG
        index = it->second;
    }
    else {
        index = code2index_[L][code] = DAG_[L].size();
        cur_pattern.reorder(vid_map); //store canonicalized pattern
        DAG_[L].emplace_back(cur_pattern);
        DAG_[L].back().code = code;
        DAG_[L].back().mat_info.resize(cur_pattern.getNumVertices());
    }

    //add parents' domains
    auto& node = DAG_[L][index];
    for (int u = 0; u < parent.pattern.getNumVertices(); u++) {
        int reu = vid_map[u];
        node.mat_info[reu].insert(*parent.mat_info[u].begin());
#ifdef DEBUG_UNION
        for (auto& info : node.mat_info[reu]) {
            range r = info.pm->getMat(info.L, info.index, info.u);
            assert(r.end - r.begin == info.size);
        }
#endif
    }
}

void PatternMining::calculateDomains(int L) {
    reserveMemory(L);
    //set this node's min_cand_size, sizes, and candidates
    if (L == 2) {
#ifdef PRINT_INTERSECTION
        calculate_domains_context = true;
#endif
        size_t i = 0;
        while (i != DAG_[L].size()) {
            string code = DAG_[L][i].code;
            //if empty, switch with back
            if (!calculateDomainsWedge(i)) {
                auto it = DAG_[L].begin() + i;
                *it = std::move(DAG_[L].back());
                DAG_[L].pop_back();
                //erase
                code2index_[L].erase(code);
            }
            else {
#ifdef DEBUG_INDEX
                vector<set<int>> gold_Dom;
                ((Alley*)estimator_)->calculateDomainsNaive(DAG_[L][i].pattern, gold_Dom);
                for (int u = 0; u < DAG_[L][i].pattern.getNumVertices(); u++) {
                    bool use_domain = true;
                    range r = GetMatQuery(L, i, u, use_domain);
                    if (use_domain) {
                        int j = 0;
                        for (int v : gold_Dom[u]) {
                            j = search(r.begin, j, r.end - r.begin, v);
                            assert(j != -1);
                            j++;
                        }
                    }
                }
#endif
                //update
                code2index_[L][code] = i;
                i++;
            }
        }
    }
    else {
        size_t index = 0;
        while (index != DAG_[L].size()) {
            auto& node = DAG_[L][index];
            string code = node.code; 
            bool empty = false;

            //check if the node has a ghost parent, i.e., not stored in the index
            if (hasEmptyParent(L, index)) {
                empty = true;
            }
            else {
                vector<MatInfo> infos, nec_infos;
                vector<range> domains, nec_domains;
                vector<star>  rels, nec_rels;
				vector<bool>  use_domain, nec_use_domain;

                vector<vector<int>> group;

                bool reduced = false;
                ///vector<vector<int>> Dom;

                int min_dom_size = MAX_INT;

                //calculate failure rate using domains
                //NEC here
                Pattern pattern(node.pattern); //copy
                VID_MAP nec_map;
                if (makeNEC(pattern, nec_map)) {
                    reduced = true;
                    assert(node.pattern.getNumVertices() > pattern.getNumVertices());

                    group.resize(pattern.getNumVertices());
                    for (int u = 0; u < node.pattern.getNumVertices(); u++)
                        group[nec_map[u]].push_back(u);

                    getMatInfo(node, group, nec_infos); 
                    min_dom_size = getDomains(pattern, nec_infos, nec_domains, nec_use_domain, nec_rels);

                    infos.resize(node.pattern.getNumVertices());
                    domains.resize(node.pattern.getNumVertices());
                    use_domain.resize(node.pattern.getNumVertices());
                    rels.resize(node.pattern.getNumVertices());

                    if (pattern.getNumVertices() == 2) {
                        setDefaultDomains(L, index, infos, rels);
                    }
                }
                else {
                    getMatInfo(node, infos); 
                    min_dom_size = getDomains(node.pattern, infos, domains, use_domain, rels);
                }

                if (min_dom_size < 10) {
                    setDefaultDomains(L, index, infos, rels);
                }
                else if (pattern.getNumVertices() != 2) {
					Alley* alley = (Alley*) estimator_;
                    if (reduced)
                        alley->SetDomains(nec_domains, nec_use_domain, nec_rels); 
                    else
                        alley->SetDomains(domains, use_domain, rels); 
                    //QueryGraph q = node.pattern.convertToQuery();
                    QueryGraph q = reduced ? pattern.convertToQuery() : node.pattern.convertToQuery();
#ifdef PRINT_INTERSECTION
                    calculate_domains_context = false;
#endif
                    estimator_->Run(g, &q, sampling_ratio_);
                    double r = estimator_->GetFailureRate();
                    sampling_cnt_++;

                    //if the rate > t
                    //calculate exact domains
                    if (r >= failure_threshold_) {
#ifdef PRINT_INTERSECTION
                        calculate_domains_context = true;
#endif
                        bool timeout = false;
                        bool print = false;
                        if (reduced && pattern.getNumVertices() > 2) {
                            //cout << "calculating reduced " << index << " " << node.code << endl;
                            empty = !(alley->calculateDomainsFast(pattern, nec_domains, nec_use_domain, nec_rels, print, timeout));
                        }
                        else {
                            //cout << "calculating " << index << " " << node.code << endl;
                            empty = !(alley->calculateDomainsFast(node.pattern, domains, use_domain, rels, print, timeout));
                        }
                        if (timeout) {
                            setDefaultDomains(L, index, infos, rels);
                        }
                        else if (!empty) {
                            //set node's mat_info and mat_ios
                            if (reduced) {
                                for (int nec = 0; nec < pattern.getNumVertices(); nec++) {
                                    int u = group[nec][0];
                                    assert(!alley->Dom_[nec].empty());
                                    appendMat(L, index, u, alley->Dom_[nec].data(), alley->Dom_[nec].size());
                                    setConsideredEdge(L, index, u, nec_rels[nec]); 
                                    size_t o = getMatOffset(L, index, u);
                                    size_t s = getMatSize(L, index, u);

                                    for (int i = 1; i < group[nec].size(); i++) {
                                        u = group[nec][i];
                                        setMatOld(L, index, u, false);
                                        setMatIndex(L, index, u, L);
                                        setMatOffset(L, index, u, o);
                                        setMatSize(L, index, u, s);
                                        setConsideredEdge(L, index, u, nec_rels[nec]); 
                                    }
                                }
                            }
                            else {
                                //append nec_Dom to materialized
                                for (int u = 0; u < node.pattern.getNumVertices(); u++) {
                                    assert(!alley->Dom_[u].empty());
                                    appendMat(L, index, u, alley->Dom_[u].data(), alley->Dom_[u].size()); 
                                    setConsideredEdge(L, index, u, rels[u]);
                                }
                            }
                            //set mat_info
                            for (int u = 0; u < node.pattern.getNumVertices(); u++) {
                                node.mat_info[u].clear();
                                node.mat_info[u].emplace(this, L, index, u, getMatSize(L, index, u));
#ifdef DEBUG_UNION
                                for (auto& info : node.mat_info[u]) {
                                    range r = info.pm->getMat(info.L, info.index, info.u);
                                    assert(r.end - r.begin == info.size);
                                }
#endif
                            }
                        }
                        //cout << "calculating done" << endl;
                    }
                    else {
                        setDefaultDomains(L, index, infos, rels);
                        assert(!empty);
                    }
                }
            }
            
            //if empty, switch with back
            if (empty) {
                auto it = DAG_[L].begin() + index;
                *it = std::move(DAG_[L].back());
                DAG_[L].pop_back();
                //erase
                code2index_[L].erase(code);
            }
            else {
                //update
                code2index_[L][code] = index;
                index++;
            }
        }
    }
}

bool PatternMining::hasEmptyParent(int L, size_t index) {
    auto& node = DAG_[L][index];

    vector<Edge> edges;
    node.pattern.getCoreLeafEdges(edges);

    //for each removal of core & leaf edge
    //BFS code
    for (auto& edge : edges) {
        //copy pattern
        Pattern p(node.pattern);
        //delete edge
        p.deleteEdge(edge.src, edge.dst, edge.el);
        //remove disconnected vertices
        if (p.getDegree(edge.src) == 0) {
            p.deleteVertex(edge.src);
        }
        else if (p.getDegree(edge.dst) == 0) {
            p.deleteVertex(edge.dst);
        }
        //encode
        VID_MAP vid_map; //node's vids (reformed) to parent's vids
        string code = encode(p, vid_map);

        //search this pattern
        auto it = code2index_[L-1].find(code);
        if (it == code2index_[L-1].end()) {
            return true;
        }
    }

    return false;
}

void PatternMining::getMatInfo(DAGNode& node, vector<vector<int>>& group, vector<MatInfo>& out_info) {
    for (int nec = 0; nec < group.size(); nec++) {
        MatInfo min_info(NULL, -1, MAX_SIZE, -1, MAX_SIZE);
        for (int u : group[nec]) {
            if (node.mat_info[u].empty()) {
                continue;
			}
            MatInfo info = *node.mat_info[u].begin();
            if (info.size < min_info.size) {
                min_info = info;
            }
        }
        assert(min_info.pm != NULL);
        assert(min_info.L >= 0);
        assert(min_info.u >= 0);
        assert(min_info.size > 0); 
        assert(min_info.size <= g->GetNumVertices());
        out_info.push_back(min_info);
    }
}

void PatternMining::getMatInfo(DAGNode& node, vector<MatInfo>& out_info) {
    for (int u = 0; u < node.pattern.getNumVertices(); u++) {
        if (node.mat_info[u].empty()) {
            out_info.emplace_back();
            continue;
        }
        MatInfo info = *node.mat_info[u].begin();
        assert(info.L > 1);
        range r = info.pm->getMat(info.L, info.index, info.u);
        assert(r.end - r.begin == info.size);
        assert(r.end - r.begin > 0);
        assert(info.L >= 0);
        assert(info.u >= 0);
        assert(info.size <= g->GetNumVertices());
        out_info.push_back(info);
    }
}

//set default search space used in both getLocalSearchSpace
//& getLocalSampleSpace in Alley
//returns min domain size
int PatternMining::getDomains(Pattern& p, vector<MatInfo>& infos, vector<range>& domains, vector<bool>& use_domain, vector<star>& rels) {
    int min_dom_size = MAX_INT;
    for (int u = 0; u < p.getNumVertices(); u++) {
        MatInfo info = infos[u];
        if (info.pm == NULL) {
			//always use vertex label! since edge label is always 0
			if (p.vorder[u] != -1) {
				domains.push_back(g->GetVertices(p.vorder[u]));
				rels.emplace_back();
			}
			else {
				size_t min_size = MAX_SIZE;
				int min_el;
				bool min_dir;
				for (auto& e : p.adj[u]) {
					size_t size = g->GetNumVertices(e.first, true);
					if (size < min_size) {
						min_size = size;
						min_el   = e.first;
						min_dir  = true;
					}
				}
				for (auto& e : p.in_adj[u]) {
					size_t size = g->GetNumVertices(e.first, false);
					if (size < min_size) {
						min_size = size;
						min_el   = e.first;
						min_dir  = false;
					}
				}

				domains.push_back(g->GetVertices(min_el, min_dir));
				int dom_size = g->GetNumVertices(min_el, min_dir);
				if (dom_size < min_dom_size)
					min_dom_size = dom_size;
				star s;
				s.emplace(min_el, min_dir);
				rels.push_back(s);
			}
            use_domain.push_back(false);
            continue;
        }
        assert(info.L > 1);
        range r = info.pm->getMat(info.L, info.index, info.u);
        assert(r.end - r.begin == info.size);
        assert(r.end - r.begin > 0);

        domains.push_back(r);
        if (info.size < min_dom_size)
            min_dom_size = info.size;

		//since we don't use filter now, use domain
		use_domain.push_back(true);
            
        star s;
        auto& ancestor = info.pm->DAG_[info.L][info.index].pattern;
        for (auto& e : ancestor.adj[info.u])
            s.emplace(e.first, true);
        for (auto& e : ancestor.in_adj[info.u])
            s.emplace(e.first, false);
        rels.push_back(s);
    }

    return min_dom_size;
}
	
void PatternMining::setDefaultDomains(int L, size_t index, vector<MatInfo>& infos, vector<star>& rels) {
    auto& node = DAG_[L][index];
    vector<MatInfo> to_insert;
    for (int u = 0; u < node.pattern.getNumVertices(); u++) {
        if (infos[u].pm == NULL) {
            size_t min_size = MAX_SIZE;
            int min_el;
            bool min_dir;
            for (auto& e : node.pattern.adj[u]) {
                size_t size = g->GetNumVertices(e.first, true);
                if (size < min_size) {
                    min_size = size;
                    min_el   = e.first;
                    min_dir  = true;
                }
            }
            for (auto& e : node.pattern.in_adj[u]) {
                size_t size = g->GetNumVertices(e.first, false);
                if (size < min_size) {
                    min_size = size;
                    min_el   = e.first;
                    min_dir  = false;
                }
            }
            setMatOld(L, index, u, 0);
            setMatIndex(L, index, u, 1);
            int vl = node.pattern.vorder[u];
            if (vl != -1 && g->GetNumVertices(vl) < min_size) {
                setMatOffset(L, index, u, vl); 
                setMatSize(L, index, u, 2);
                min_size = g->GetNumVertices(vl);
            }
            else {
                setMatOffset(L, index, u, min_el); 
                setMatSize(L, index, u, min_dir ? 1 : 0);
                setConsideredEdge(L, index, u, min_el, min_dir); 
            }
            to_insert.emplace_back(this, L, index, u, min_size);
            continue;
        }
        
        MatInfo& info = infos[u];
		if (info.pm != this)
			setMatOld(L, index, u, true);
        else
			setMatOld(L, index, u, false);
        setMatIndex(L, index, u, info.pm->getMatIndex(info.L, info.index, info.u));
        setMatOffset(L, index, u, info.pm->getMatOffset(info.L, info.index, info.u)); 
        setMatSize(L, index, u, info.pm->getMatSize(info.L, info.index, info.u)); 
        setConsideredEdge(L, index, u, rels[u]); 
        to_insert.push_back(info);
    }
    assert(to_insert.size() == node.pattern.getNumVertices());
    for (int u = 0; u < node.pattern.getNumVertices(); u++) {
        node.mat_info[u].clear();
        node.mat_info[u].insert(to_insert[u]);
#ifdef DEBUG_UNION
        for (auto& info : node.mat_info[u]) {
            range r = info.pm->getMat(info.L, info.index, info.u);
            assert(r.end - r.begin == info.size);
        }
#endif
    }
}

void PatternMining::reserveMemory(int L) {
    size_t mat_size = 0;
    mat_size = INIT_MAT_CAPACITY;
    mat_old_[L].resize((L+1) * DAG_[L].size()); 
    mat_i_[L].resize((L+1) * DAG_[L].size()); 
    mat_o_[L].resize((L+1) * DAG_[L].size()); 
    mat_s_[L].resize((L+1) * DAG_[L].size()); 
    consideredEdge_[L].resize((L+1) * DAG_[L].size());
    //enlarge on-the-fly
    materialized_[L].reserve(mat_size);
}

//calculate exact domains
//return whether the node has non-empty domains
//if true, set mat_ios and node.mat_info
bool PatternMining::calculateDomainsWedge(size_t index) {
    clearMatInfo(2, index);
    auto& node = DAG_[2][index];

    int center = -1;
    vector<int> uids;
    vector<int> els, vls;
    vector<bool> dirs;
    set<range> ranges;

    for (center = 0; center < node.pattern.getNumVertices(); center++) {
        if (node.pattern.getDegree(center) > 1) {
            for (auto& e : node.pattern.adj[center]) {
                uids.push_back(e.second);
                vls.push_back(node.pattern.vorder[e.second]);
                els.push_back(e.first);
                dirs.push_back(true);
                ranges.insert(g->GetVertices(els.back(), dirs.back()));
            }
            for (auto& e : node.pattern.in_adj[center]) {
                uids.push_back(e.second);
                vls.push_back(node.pattern.vorder[e.second]);
                els.push_back(e.first);
                dirs.push_back(false);
                ranges.insert(g->GetVertices(els.back(), dirs.back()));
            }
            if (node.pattern.vorder[center] != -1)
                ranges.insert(g->GetVertices(node.pattern.vorder[center]));
            break;
        }
    }

    if (ranges.size() == 1) {
        if (vls[0] == -1 && vls[1] == -1) {
            for (int u = 0; u < node.pattern.getNumVertices(); u++) {
                setMatOld(2, index, u, false);
                setMatIndex(2, index, u, 1);
                setMatOffset(2, index, u, els[0]); //edge label
                if (u == center) {
                    setMatSize(2, index, u, dirs[0] ? 1 : 0); //direction
                    setConsideredEdge(2, index, u, els[0], dirs[0]);
                }
                else {
                    setMatSize(2, index, u, dirs[0] ? 0 : 1); //direction
                    setConsideredEdge(2, index, u, els[0], !dirs[0]);
                }
            }

            //set mat_info
            for (int u = 0; u < node.pattern.getNumVertices(); u++) {
                node.mat_info[u].clear();
                node.mat_info[u].emplace(this, 2, index, u, u == center ? g->GetNumVertices(els[0], dirs[0]) : g->GetNumVertices(els[0], !dirs[0]));
#ifdef DEBUG_UNION
                for (auto& info : node.mat_info[u]) {
                    range r = info.pm->getMat(info.L, info.index, info.u);
                    assert(r.end - r.begin == info.size);
                }
#endif
            }

            return true;
        }
    }

	size_t to_reserve = ranges.begin()->end - ranges.begin()->begin + g->GetNumVertices(els[0], !dirs[0]) + g->GetNumVertices(els[1], !dirs[1]);
    size_t before = materialized_[2].size();
    materialized_[2].reserve(before + to_reserve);
	const int* cand = ranges.begin()->begin;
    int* arr;
    int size;
    
    if (ltype_ == EDGE_LABEL) {
        arr = materialized_[2].data() + before;
        size = intersect(ranges, arr, NULL);
    }
    else {
        int u = center;
        if (ranges.size() > 1) {
            if (to_reserve > intersect_capacity_[u]) {
                intersect_capacity_[u] = std::max((int)to_reserve, 2 * intersect_capacity_[u]);
                int* b = intersect_[u];
                int* t = temp_[u];
                intersect_[u] = new int[intersect_capacity_[u]];
                temp_[u] = new int[intersect_capacity_[u]];
                delete [] b;
                delete [] t;
            }
            size = intersect(ranges, intersect_[u], temp_[u]);
            arr = intersect_[u];
        }
        else {
            size = ranges.begin()->end - ranges.begin()->begin;
            arr = const_cast<int*>(ranges.begin()->begin);
        }
        //need to copy to materialized_[2] later
    }

    if (size == 0)
        return false;
    
    if (ltype_ == EDGE_LABEL) {
        materialized_[2].resize(before + size);
        setMatOffset(2, index, center, before);
        setMatSize(2, index, center, size);
    }
    //for VERTEX_LABEL, set later
    setConsideredEdge(2, index, center, els[0], dirs[0]); 
    setConsideredEdge(2, index, center, els[1], dirs[1]); 

    assert(uids.size() == 2);

    if (uids[0] == uids[1]) {
        assert(uids[0] != center);
        if (ltype_ == EDGE_LABEL)
            assert(els[0] != els[1] || dirs[0] != dirs[1]); 
        else {
            assert(vls[0] == vls[1]);
        }

        size_t s = materialized_[2].size(); 
        size_t e = s; 
        vector<int> center_ids;

        for (int i = 0; i < size; i++) {
            range r1 = g->GetAdj(arr[i], els[0], dirs[0]); 
            range r2 = g->GetAdj(arr[i], els[1], dirs[1]); 
            assert(r1.end != r1.begin);
            assert(r2.end != r2.begin);

            set<range> rs;
            rs.insert(r1);
            rs.insert(r2);
            
            if (ltype_ == EDGE_LABEL) {
                int other_size = intersect(rs, materialized_[2].data() + e, NULL); 
                e += other_size;
            }
            else {
                if (vls[0] != -1)
                    rs.insert(g->GetVertices(vls[0]));
                assert(rs.size() > 1);
                int u = uids[0];
                int min_size = rs.begin()->end - rs.begin()->begin;
                if (min_size > intersect_capacity_[u]) {
                    intersect_capacity_[u] = std::max(min_size, 2 * intersect_capacity_[u]);
                    int* b = intersect_[u];
                    int* t = temp_[u];
                    intersect_[u] = new int[intersect_capacity_[u]];
                    temp_[u] = new int[intersect_capacity_[u]];
                    delete [] b;
                    delete [] t;
                }
                int other_size = intersect(rs, intersect_[u], temp_[u]);
                if (other_size > 0) {
                    for (int i = 0; i < other_size; i++) {
                        int v = intersect_[u][i];
                        if (!marked_[u][v]) {
                            marked_[u][v] = true;
                            materialized_[2].push_back(v);
                            e++;
                        }
                    }
                    center_ids.push_back(arr[i]);
                }
            }
        }

        if (s == e) {
            assert(center_ids.empty());
            return false;
        }

        materialized_[2].resize(e);
        std::sort(materialized_[2].begin() + s, materialized_[2].end());
        if (ltype_ == EDGE_LABEL) {
            materialized_[2].erase(std::unique(materialized_[2].begin() + s, materialized_[2].end()));
            e = materialized_[2].size();
        }
        else {
            for (int i = s; i < e; i++)
                marked_[uids[0]][materialized_[2][i]] = false;
        }
        assert(e == materialized_[2].size());

        setMatOffset(2, index, uids[0], s);
        setMatSize(2, index, uids[0], e - s);

        //now set center's mat
        if (ltype_ == VERTEX_LABEL) {
            appendMat(2, index, center, center_ids.data(), center_ids.size()); 
        }
        setConsideredEdge(2, index, uids[0], els[0], !dirs[0]); 
        setConsideredEdge(2, index, uids[1], els[1], !dirs[1]); 
    }
    else {
        if (ltype_ == EDGE_LABEL) {
            for (int k = 0; k < 2; k++) {
                int u = uids[k];
                int el = els[k];
                bool dir = dirs[k];

                before = materialized_[2].size();
                for (int i = 0; i < size; i++) {
                    range r = g->GetAdj(arr[i], el, dir); 
                    assert(r.end != r.begin);
                    for (; r.begin != r.end; ++r.begin) {
                        int v = *r.begin;
                        if (!marked_[u][v]) {
                            materialized_[2].push_back(v);
                            marked_[u][v] = true;
                        }
                    }
                }

                std::sort(materialized_[2].begin() + before, materialized_[2].end());
                setMatOffset(2, index, u, before);
                setMatSize(2, index, u, materialized_[2].size() - before);
                setConsideredEdge(2, index, u, el, !dir);

                for (size_t o = before; o < materialized_[2].size(); o++) {
                    int v = materialized_[2][o];
                    assert(marked_[u][v]);
                    marked_[u][v] = false;
                }
            }
        }
        else {
            vector<int> center_ids;
            for (int i = 0; i < size; i++) {
                bool match = true;
                for (int k = 0; k < 2; k++) {
                    int u = uids[k];
                    if (vls[k] != -1) {
                        range r = g->GetAdj(arr[i], els[k], dirs[k]); 
                        range vr = g->GetVertices(vls[k]);
                        assert(r.end - r.begin < intersect_capacity_[u]);
                        int s;
                        if (r.end - r.begin < vr.end - vr.begin)
                            s = intersect(r, vr, intersect_[u]); 
                        else
                            s = intersect(vr, r, intersect_[u]); 
                        if (s == 0) {
                            match = false;
                            break;
                        }
                    }
                }

                if (match) {
                    center_ids.push_back(arr[i]);
                }
            }

            if (center_ids.empty())
                return false;
            
            //now set center's mat
            appendMat(2, index, center, center_ids.data(), center_ids.size());
            for (int k = 0; k < 2; k++) {
                before = materialized_[2].size();

                int u = uids[k];
                int el = els[k];
                bool dir = dirs[k];

                for (int c : center_ids) {
                    range r = g->GetAdj(c, el, dir); 
                    if (vls[k] != -1) {
                        range vr = g->GetVertices(vls[k]);
                        int s = intersect(r, vr, intersect_[u]); 
                        for (int i = 0; i < s; i++) {
                            int v = intersect_[u][i]; 
                            if (!marked_[u][v]) {
                                materialized_[2].push_back(v);
                                marked_[u][v] = true;
                            }
                        }
                    }
                    else {
                        int s = r.end - r.begin;
                        for (int i = 0; i < s; i++) {
                            int v = r.begin[i]; 
                            if (!marked_[u][v]) {
                                materialized_[2].push_back(v);
                                marked_[u][v] = true;
                            }
                        }
                    }
                }
                
                std::sort(materialized_[2].begin() + before, materialized_[2].end());
                setMatOffset(2, index, u, before);
                setMatSize(2, index, u, materialized_[2].size() - before);
                setConsideredEdge(2, index, u, el, !dir);
                for (size_t o = before; o < materialized_[2].size(); o++) {
                    int v = materialized_[2][o];
                    assert(marked_[u][v]);
                    marked_[u][v] = false;
                }
            }
        }
    }
    
    //set mat_info
    for (int u = 0; u < node.pattern.getNumVertices(); u++) {
        node.mat_info[u].clear();
        node.mat_info[u].emplace(this, 2, index, u, getMatSize(2, index, u));
#ifdef DEBUG_UNION
        for (auto& info : node.mat_info[u]) {
            range r = info.pm->getMat(info.L, info.index, info.u);
            assert(r.end - r.begin == info.size);
        }
#endif
    }

    return true;
}

bool makeNEC(Pattern& p, VID_MAP& NEC_mapping) {
    vector< set<int> > NEC;
    
    typedef map< pair<int, int>, int > TYPE_TOPOLOGY;
    NEC_mapping.resize(p.getNumVertices());

    // Make an initial unordered partition
    // First of all, grouped vertices having same edges
    // If two vertices in a group have an edge, separate them
    // For each vertex i,
    {
        vector< TYPE_TOPOLOGY > NEC_topology_out, NEC_topology_in;

        for (int i=0; i < p.getNumVertices(); i++) {
            TYPE_TOPOLOGY topology_out, topology_in;
            topology_out=p.getTopology(i, true);
            topology_in=p.getTopology(i, false);
            
            bool findNEC=false;
            for (int j=0; j < NEC.size(); j++) {
                assert(NEC.size() == NEC_topology_out.size() && NEC.size() == NEC_topology_in.size());
                if (NEC_topology_out[j] != topology_out || NEC_topology_in[j] != topology_in) {
                    continue;
                }
                bool canBeAdded=true;
                for (auto& v : NEC[j]) {
                    if (p.hasEdge(v, i) || p.hasEdge(i, v)) {
                        canBeAdded=false;
                        break;
                    }
                }
                if(canBeAdded) {
                    findNEC=true;
                    NEC[j].insert(i);
                    NEC_mapping[i]=j;
                    break;
                }
            }
            if (!findNEC) {
                NEC_mapping[i]=NEC.size();
                NEC.push_back( set<int>({i}) );
                NEC_topology_out.push_back( topology_out );
                NEC_topology_in.push_back( topology_in );
            }
        }
    }
    //Repeat calculate shattering
    while(1) {
        pair<int, int> shatter=std::pair<int, int>(-1, -1); //second shatters first
        vector< set<int> > new_NEC;
        bool isCycle=false;

        for(int i=0; i<NEC.size(); i++) {
            if(NEC[i].size() > 1) {
                //check if cycle only with other NECs whose size > 2
                if(p.hasUnacceptableCycleInNEC(NEC, NEC_mapping, i, new_NEC)) {
                    shatter.first=i;
                    isCycle=true;
                    break;
                }

                //calculate degree to Vj
                for(int j=0; j<NEC.size(); j++) {
                    if(j <= i && NEC[j].size() > 1) continue;

                    TYPE_TOPOLOGY topology_out_to_j, topology_in_to_j;
                    topology_out_to_j = p.getTopology(*(NEC[i].begin()), true, true, NEC[j]);
                    topology_in_to_j = p.getTopology(*(NEC[i].begin()), false, true, NEC[j]);
                    for(auto& v : NEC[i]) {
                        if(p.getTopology(v, true, true, NEC[j]) != topology_out_to_j || p.getTopology(v, false, true, NEC[j]) != topology_in_to_j) {
                            shatter=std::pair<int, int>(i,j);
                           break;
                       }
                    }
                    if(shatter.first != -1)
                        break;
                }
            }
            if(shatter.first != -1)
                break;
        }
        if(shatter.first == -1) {
            break;
        }

        //Decompose Vi
        if(!isCycle) { 
            //Grouping vertices having the same topology
            vector< TYPE_TOPOLOGY > NEC_topology_out, NEC_topology_in;
            for (auto& v : NEC[shatter.first]) {
                TYPE_TOPOLOGY topology_out, topology_in;
                topology_out=p.getTopology(v, true, true, NEC[shatter.second]);
                topology_in=p.getTopology(v, false, true, NEC[shatter.second]);
                
                bool findNEC=false;
                for (int j=0; j < new_NEC.size(); j++) {
                    assert(new_NEC.size() == NEC_topology_out.size() && new_NEC.size() == NEC_topology_in.size());
                    if (NEC_topology_out[j] != topology_out || NEC_topology_in[j] != topology_in) {
                        continue;
                    }
                    findNEC=true;
                    new_NEC[j].insert(v);
                    break;
                }
                if (!findNEC) {
                    new_NEC.push_back( set<int>( {v} ) );
                    NEC_topology_out.push_back( topology_out );
                    NEC_topology_in.push_back( topology_in );
                }
            }
        }

        //remove Vi and insert new NECs
        assert(new_NEC.size() > 1);
        //NEC mapping update
        for(int i=shatter.first+1; i<NEC.size(); i++) {
            for(auto& v : NEC[i]) {
                NEC_mapping[v]-=1;
            }
        }
        for(int i=0; i<new_NEC.size(); i++) {
            for(auto& v : new_NEC[i]) {
                NEC_mapping[v]=NEC.size()+i-1;
            }
        }

        NEC.erase(NEC.begin()+shatter.first);
        NEC.insert(NEC.end(), new_NEC.begin(), new_NEC.end());
    }

    if(NEC.size() == p.getNumVertices())
        return false;

    //make VID_MAP
    for(int i=0; i<NEC.size(); i++) {
        for(auto& v : NEC[i]) {
            NEC_mapping[v]=i;
        }
    }

    //Reshape p
    Pattern p_new;
    for(int i=0; i<NEC.size(); i++) {
        int v=*(NEC[i].begin());
        p_new.addVertex( p.vorder[v] );
    }
    for(int i=0; i<NEC.size(); i++) {
        int v=*(NEC[i].begin());
        for(auto& u : p.adj[v]) {
            if(!p_new.hasDirectedEdge( i, NEC_mapping[u.second] ) )
                p_new.addEdge( i , NEC_mapping[u.second] , u.first );
        }
    }
    p=p_new;

    return true;
}

//used in query mode
bool getDomainsRecursive(const vector<PatternMining*>& pms,  Pattern& q, const vector<int>& order2vid, int minMaxL, int maxMaxL, vector<range>& domains, vector<bool>& hasDomains, vector<star>& consideredElist, const int numVerticesWithDomains) 
{
    assert(order2vid.size() == q.getNumVertices());
    if(numVerticesWithDomains == 0) {

        bool replace = true;
        if (domains.empty()) {
            domains.resize(q.getNumVertices());
            for (auto& r : domains) {
                r.begin=NULL;
                r.end=NULL;
            }
            hasDomains.resize(q.getNumVertices(), false);
            consideredElist.resize(q.getNumVertices());
            replace = false;
        }
    }

    vector<int> vid2order(q.getNumVertices());
    vector<int> mindices(q.getNumVertices());
    for(int i=0; i<q.getNumVertices(); i++) {
        vid2order[ order2vid[i] ]=i;
    }

    int maxNumVinP=maxMaxL;
    if(pms.size() == 1)
        maxNumVinP+=1;
    // vid2order[ v ] % (MAX_INDEX_LEVEL+1) -> renumbered index
    auto order2patternvid = [maxNumVinP] (int order) {
        return order%maxNumVinP;
    };

    int cnt_getdom=0;
    int start_order=0;
    int end_order=std::min(start_order+maxNumVinP, (int)(q.getNumVertices()));
    Pattern p;
    vector<Edge> elist; // add vertex vi and edges (vi, X) (X, vi) to a pattern
    for(int cur_order=start_order; cur_order<end_order; cur_order++) {
        //Add edges to p
        p.addVertex(q.vorder[ order2vid[ cur_order ] ]);
        if(p.getNumVertices() == 1) continue;
        for(auto& v: q.adj[ order2vid[ cur_order ] ]) { //order[i] to v
            if(vid2order[ v.second ] >= start_order && vid2order[ v.second ] < cur_order) {
                assert(p.getNumVertices() > order2patternvid( vid2order[v.second] ));
                p.addEdge(order2patternvid(cur_order), order2patternvid( vid2order[v.second] ) , v.first);
                Edge e(order2patternvid(cur_order), order2patternvid( vid2order[v.second] ), v.first);
                elist.push_back(e);
            }
        }
        for(auto& v: q.in_adj[ order2vid[ cur_order ] ]) { //v to order[i]
            if(vid2order[ v.second ] >= start_order && vid2order[ v.second ] < cur_order) {
                p.addEdge(order2patternvid(vid2order[v.second]), order2patternvid(cur_order), v.first);
                Edge e(order2patternvid(vid2order[v.second]), order2patternvid(cur_order), v.first);
                elist.push_back(e);
            }
        }
        int total_edges=0;
        total_edges = elist.size();

        if(total_edges >= minMaxL || cur_order == end_order - 1) { // search index
            vector<vector<int>> considered_deleted_edges_stack;
            considered_deleted_edges_stack.push_back( vector<int>() );
            int cnt_stack=0;
            do { //search for each maximal connected subgrph
                vector<int> stack = considered_deleted_edges_stack.back();
                considered_deleted_edges_stack.pop_back();
                Pattern p2(p);
                for(int i=0; i<stack.size(); i++)
                    p2.deleteEdge( elist[ stack[i] ].src, elist[ stack[i] ].dst, elist[ stack[i] ].el );

                if(total_edges - stack.size() <= maxMaxL) {
                    //find an instance!
                    cnt_getdom++;
                    cnt_stack++;
                    int L = total_edges - stack.size();
                    size_t Matid;
                    for(int mi=0; mi < pms.size(); mi++) {
                        auto& pm=*(pms[mi]);
                        if ( mi > 0 ) {
                            if(pm.GetLTypeInChar() == 'e')
                                p2.convertELabels(pm.getLabel2Order());
                            else if(pm.GetLTypeInChar() == 'v')
                                p2.convertVLabels(pm.getLabel2Order());
                        }
                        if(pm.getMaxL() < L) continue;

                        VID_MAP mapping;
                        string code=pm.encode(p2, mapping);
                        if(pm.code2Index( L, code, Matid )) {
                            // patternvid2order: start_order + patternvid
                            // order2vid[ start_order + patternvid ]
                            for( int i=0; i<p2.getNumVertices(); i++ ) {
                                assert(p2.getDegree(i) != 0);
                                int vid=order2vid[ start_order + i ];
                                bool use_domain = true;
                                range r=pm.GetMatQuery( L, Matid, mapping[i], use_domain );
                                if (!use_domain) {
                                    continue;
                                }
                                hasDomains[ vid ]=true;
                                if(r.end - r.begin == 0) {
                                    assert(false);
                                    return false;
                                }
                                if(domains[ vid ].begin == NULL || domains[vid].end-domains[vid].begin > r.end-r.begin) {
                                    domains[ vid ]=r;
                                    mindices[ vid ]=mi;
                                    consideredElist[ vid ].clear();
                                    for(auto& u: p2.adj[i])
                                        consideredElist[ vid ].insert( make_pair( u.first, true ) );
                                    for(auto& u: p2.in_adj[i])
                                        consideredElist[ vid ].insert( make_pair( u.first, false ) );
                                    pm.GetConsideredElist( L, Matid, mapping[i], consideredElist[ vid ]); 
                                }
                            }
                            break;
                        }
                        else if( L > 1 ) {
                            if (pm.isAcceptablePattern(p2)) // We found a domain with size 0!!!!
                                return false;
                        }
                    }
                    if(cnt_stack > 10)
                        break;

                    continue;
                }
                else if(stack.size() > 0 && stack.back() == elist.size() - 1)
                    continue;

                vector<Edge> canBeDeleted;
                p2.getCoreEdges(canBeDeleted);
                for(auto& edge: canBeDeleted) {
                    int eid  = find(elist.begin(), elist.end(), edge) - elist.begin();
                    assert( eid < elist.size() );
                    if( stack.empty() || eid > stack.back() ) {
                        considered_deleted_edges_stack.push_back( stack );
                        considered_deleted_edges_stack.back().push_back( eid );
                    }
                }
            } while(!considered_deleted_edges_stack.empty());
        }
        if(cnt_getdom > 20)
            break;
    }
    start_order=end_order;
    int prunedVertices=maxNumVinP + numVerticesWithDomains;

    if(prunedVertices > (int)q.getNumVertices() - minMaxL) {
        //convert orders to labels
        for (int i = 0; i < consideredElist.size(); i++) {
            star s;
            for (auto it = consideredElist[i].begin(); it != consideredElist[i].end(); it++) {
                auto& pm=*(pms[mindices[i]]);
                if (pm.getOrder2LabelSize(it->first) == 1)
                    s.emplace(pm.getOrder2Label(it->first), it->second);
            }
            consideredElist[i] = s;
        }
        return true;
    }

    int numV = q.getNumVertices();
    vector<int> new_order2vid;
    vector<bool> visited(numV, false);

    priority_queue<pair<int, int>, vector< pair<int, int> >, greater<pair<int, int>> > pqueue;
    auto order2priority = [numV, maxNumVinP] (int order) {
        return (order+numV-maxNumVinP)%(numV);
    };
    pqueue.push( make_pair( order2priority(end_order), order2vid[end_order] ) );
    visited[ order2vid[end_order] ]=true;

    while(!pqueue.empty()) {
        auto u=pqueue.top();
        pqueue.pop();

        new_order2vid.push_back(u.second);

        for(auto& v: q.adj[u.second]) {
            if(! visited[v.second]) {
                pqueue.push( make_pair( order2priority( vid2order[v.second] ), v.second ) );
                visited[ v.second ] = true;
            }
        }
        for(auto& v: q.in_adj[u.second]) {
            if(!visited[v.second]) {
                pqueue.push( make_pair( order2priority( vid2order[v.second] ), v.second ) );
                visited[ v.second ] = true;
            }
        }
    }

    for(int u=0; u<q.getNumVertices(); u++) {
        if(!visited[u]) {
            new_order2vid.push_back(u);
            visited[u]=true;
        }
    }

    bool hasZeroDomain = getDomainsRecursive(pms, q, new_order2vid, minMaxL, maxMaxL, domains, hasDomains, consideredElist, prunedVertices);

    return hasZeroDomain;
}

//used in query mode
bool getDomains(vector<PatternMining*>& pms, Pattern& q, vector<int>& order2vid, int minMaxL, int maxMaxL, vector<range>& domains, vector<bool>& hasDomains, vector<star>& consideredElist) 
{
    assert(order2vid.size() == q.getNumVertices());
    vector<int> vid2order(q.getNumVertices());
    vector<int> mindices(q.getNumVertices());

    bool replace = true;
    if (domains.empty()) {
        domains.resize(q.getNumVertices());
        for (auto& r : domains) {
            r.begin=NULL;
            r.end=NULL;
        }
        hasDomains.resize(q.getNumVertices(), false);
        consideredElist.resize(q.getNumVertices());
        replace = false;
    }

    for(int i=0; i<q.getNumVertices(); i++) {
        vid2order[ order2vid[i] ]=i;
    }

    // vid2order[ v ] % (MAX_INDEX_LEVEL+1) -> renumbered index
    auto order2patternvid = [maxMaxL] (int order) {
        return order%(maxMaxL+1);
    };

    int cnt_getdom=0;
    int start_order=0;
    while(start_order < q.getNumVertices() ) {
        int end_order=std::min(start_order+(maxMaxL+1), (int)(q.getNumVertices()));
        Pattern p;
        vector<Edge> elist; // add vertex vi and edges (vi, X) (X, vi) to a pattern
        for(int cur_order=start_order; cur_order<end_order; cur_order++) {
            //Add edges to p
            p.addVertex(q.vorder[ order2vid[ cur_order ] ]);
            if(p.getNumVertices() == 1) continue;
            for(auto& v: q.adj[ order2vid[ cur_order ] ]) { //order[i] to v
                if(vid2order[ v.second ] >= start_order && vid2order[ v.second ] < cur_order) {
                    assert(p.getNumVertices() > order2patternvid( vid2order[v.second] ));
                    p.addEdge(order2patternvid(cur_order), order2patternvid( vid2order[v.second] ) , v.first);
                    Edge e(order2patternvid(cur_order), order2patternvid( vid2order[v.second] ), v.first);
                    elist.push_back(e);
                }
            }
            for(auto& v: q.in_adj[ order2vid[ cur_order ] ]) { //v to order[i]
                if(vid2order[ v.second ] >= start_order && vid2order[ v.second ] < cur_order) {
                    p.addEdge(order2patternvid(vid2order[v.second]), order2patternvid(cur_order), v.first);
                    Edge e(order2patternvid(vid2order[v.second]), order2patternvid(cur_order), v.first);
                    elist.push_back(e);
                }
            }
            int total_edges=0;
            total_edges = elist.size();

            if(total_edges >= minMaxL || cur_order == end_order - 1) { // search index
                vector<vector<int>> considered_deleted_edges_stack;
                considered_deleted_edges_stack.push_back( vector<int>() );
                int cnt_stack=0;
                do {
                    vector<int> stack = considered_deleted_edges_stack.back();
                    considered_deleted_edges_stack.pop_back();
                    Pattern p2(p);
                    for(int i=0; i<stack.size(); i++)
                        p2.deleteEdge( elist[ stack[i] ].src, elist[ stack[i] ].dst, elist[ stack[i] ].el );

                    if(total_edges - stack.size() <= maxMaxL) {
                        //find an instance!
                        cnt_getdom++;
                        cnt_stack++;
                        int L = total_edges - stack.size();
                        size_t Matid;
                        for(int mi=0; mi < pms.size(); mi++) {
                            auto& pm=*(pms[mi]);
                            if ( mi > 0 ) {
                                //p2.convertLabels(pm.getLabel2Order());
                                if(pm.GetLTypeInChar() == 'e')
                                    p2.convertELabels(pm.getLabel2Order());
                                else if(pm.GetLTypeInChar() == 'v')
                                    p2.convertVLabels(pm.getLabel2Order());
                            }
                            if(pm.getMaxL() < L) continue;

                            VID_MAP mapping;
                            string code=pm.encode(p2, mapping);
                            if(pm.code2Index( L, code, Matid )) {
                                // patternvid2order: start_order + patternvid
                                // order2vid[ start_order + patternvid ]
                                for( int i=0; i<p2.getNumVertices(); i++ ) {
                                    assert(p2.getDegree(i) != 0);
                                    int vid=order2vid[ start_order + i ];
                                    bool use_domain = true;
                                    range r=pm.GetMatQuery( L, Matid, mapping[i], use_domain );
                                    if (!use_domain) {
                                        continue;
                                    }
                                    hasDomains[ vid ]=true;
                                    if(r.end - r.begin == 0) {
                                        assert(false);
                                        return false;
                                    }
                                    if(domains[ vid ].begin == NULL || domains[vid].end-domains[vid].begin > r.end-r.begin) {
                                        domains[ vid ]=r;
                                        mindices[ vid ]=mi;
                                        consideredElist[ vid ].clear();
                                        for(auto& u: p2.adj[i])
                                            consideredElist[ vid ].insert( make_pair( u.first, true ) );
                                        for(auto& u: p2.in_adj[i])
                                            consideredElist[ vid ].insert( make_pair( u.first, false ) );
                                        pm.GetConsideredElist( L, Matid, mapping[i], consideredElist[ vid ]); 
                                    }
                                }
                                break;
                            }
                            else if( L > 1 ) {
                                if (pm.isAcceptablePattern(p2))
                                    return false;
                            }
                        }
                        if(cnt_stack > 15)
                            break;

                        continue;
                    }
                    else if(stack.size() > 0 && stack.back() == elist.size() - 1)
                        continue;

                    vector<Edge> canBeDeleted;
                    p2.getCoreEdges(canBeDeleted);
                    for(auto& edge: canBeDeleted) {
                        int eid  = find(elist.begin(), elist.end(), edge) - elist.begin();
                        assert( eid < elist.size() );
                        if( stack.empty() || eid > stack.back() ) {
                            considered_deleted_edges_stack.push_back( stack );
                            considered_deleted_edges_stack.back().push_back( eid );
                        }
                    }
                } while(!considered_deleted_edges_stack.empty());
            }
            if(cnt_getdom > 30)
                break;
        }
        start_order=end_order;
        break; 
    }

    //convert orders to labels
    for (int i = 0; i < consideredElist.size(); i++) {
        star s;
        for (auto it = consideredElist[i].begin(); it != consideredElist[i].end(); it++) {
            auto& pm=*(pms[mindices[i]]);
			if (pm.getOrder2LabelSize(it->first) == 1)
				s.emplace(pm.getOrder2Label(it->first), it->second);
		}
        consideredElist[i] = s;
    }

    return true;
}
