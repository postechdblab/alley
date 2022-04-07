#include <cassert>
#include <queue>
#include <algorithm>
#include <unordered_set>
#include <random>
#include <chrono>
#include "../include/intersection.h"
#include "../include/alley.h"

size_t init_totalEdgeCase1;
size_t init_totalEdgeCase2;
size_t init_totalLenRare;
size_t init_totalLenFreq;
size_t init_totalGalloping;
size_t init_totalShuffling;
size_t init_totalA;
size_t init_totalB;
size_t init_totalCountGalloping;
size_t init_totalCountShuffling;

/************************** EXTERNAL APIs **************************/

Alley::Alley(bool pm) {
    pattern_mining_ = pm;
    for (int i = 0; i < MAX_ORDER; i++) {
        intersect_capacity_[i] = INITIAL_CAPACITY;
        intersect_[i] = new int[INITIAL_CAPACITY];
        temp_[i] = new int[INITIAL_CAPACITY];
    }
}

Alley::~Alley() {
	for (int i = 0; i < MAX_ORDER; i++) {
		delete [] intersect_[i];
		delete [] temp_[i];
    }
    if (maxL_ != -1) {
        for (int u = 0; u < maxL_+1; u++) {
            delete [] initial_intersect_[u];
            delete [] initial_temp_[u];
        }

        delete [] initial_intersect_capacity_;
        delete [] initial_intersect_;
        delete [] initial_temp_;
    }
    ClearCache();
}

void Alley::PrepareSummaryStructure(DataGraph& g, double ratio) {
	if (pattern_mining_) {
		clock_t start, end;
		int maxL = stoi(getenv("GCARE_ALLEY_TPI_MAXL"));
        double failure_threshold = stod(getenv("GCARE_ALLEY_TPI_FAIL_THRESHOLD"));
        int num_group = stoi(getenv("GCARE_ALLEY_TPI_NUM_GROUP"));
		miner_ = new MinerWrapper(&g, this);
        start = clock();
        label_type lt = EDGE_LABEL;
        if(g.GetNumELabels() == 1 && g.GetNumVLabels() > 1)
            lt = VERTEX_LABEL;
		miner_->BuildIndex(maxL, ratio, failure_threshold, num_group, lt);
		end = clock();
		std::cout << "maxL " << maxL << " build index done in ... " << (double)(end-start)/CLOCKS_PER_SEC << " secs" << endl;
	}
}

void Alley::WriteSummary(const char* fn) {
	if (pattern_mining_) {
		miner_->Write(fn);
	}
}

void Alley::ReadSummary(const char* fn) {
	if (pattern_mining_) {
		miner_ = new MinerWrapper(); 
		miner_->Read(fn);
	}
}

void Alley::Init() {
#ifdef PRINT_INTERSECTION
    if (!in_build_)
        totalEdgeCase1 = totalEdgeCase2 = totalLenRare = totalLenFreq = totalGalloping = totalA = totalB = totalShuffling = totalCountGalloping = totalCountShuffling = 0;
#endif
    generator_.seed(rand());

    y1_ = y1_sum_ = sample_cnt_ = sample_size_ = last_ = 0;
    valid_ = true;

    plan_.clear();
    backward_.clear();

    rbi_.clear();
    rbi_.resize(q->GetNumVertices());

    initial_backward_.clear();
    sample_.clear();
    sample_space_size_.clear();
    num_adj_dvids_.clear();
    candidates_.clear();
    last_stopped_order_ = -1;
    plan_.resize(q->GetNumVertices(), -1);
    backward_.resize(q->GetNumVertices(), 0);
    initial_backward_.resize(q->GetNumVertices(), 0);
    sample_.resize(q->GetNumVertices(), -1);
    sample_space_size_.resize(q->GetNumVertices(), std::numeric_limits<int>::max());
    num_adj_dvids_.resize(q->GetNumVertices(), 0);

    reduced_rel_.clear();
    reduced_adj_.clear();
    reduced_rel_.resize(q->GetNumVertices());
    reduced_adj_.resize(q->GetNumVertices());

    //set sample_size_ and sample_cnt_ 
    if (in_build_) {
        recur_call_cnt_ = sample_cnt_ = sample_size_ = 0;
    }
    else {
        int sum = 0;
        int vl_cnt = 0;
        int e_cnt = q->GetNumEdges(); 
        for (int v = 0; v < q->GetNumVertices(); v++) {
            int vl = q->GetVLabel(v); 
            if (vl != -1) {
                sum += g->GetNumVertices(vl); 
                vl_cnt++;
            }
            for (auto& e : q->GetAdj(v, true)) 
                sum += g->GetNumEdges(e.second);
        }
        sample_size_ = sum / (vl_cnt + e_cnt); 
        sample_size_ = ceil(sample_size_ * sample_ratio);
        sample_cnt_ = sample_size_;
    }

    if (!validateBounds()) {
        valid_ = false;
        return;
    }

    getStarSelectivities();

    //set plan_
    determineOrder();

    if (pattern_mining_ && !in_build_) {
        domains_.clear();
        use_domain_.clear();
        dom_considered_rels_.clear();
        
        getDomainsStatic();
    }

    mergeIntersection();

    clock_t start, end;
    start = clock();
    if (in_build_)
        getLocalSampleSpaceBuild(0);
    else
        getLocalSampleSpaceQuery(0);
    end = clock();
    double t = (double)(end-start)/CLOCKS_PER_SEC;
    initial_sample_space_time_ += t;

    if (sample_space_size_[plan_[0]] == 0 || last_ == 0) {
        valid_ = false;
        return;
    }

    if (in_build_) {
		sum_y1_ = 0;
        //calculate max sampling count
        double upper_bound = 384;
        max_y1_cnt_ = ceil(upper_bound / (1 + upper_bound/sample_space_size_[plan_[0]])); 
        assert(max_y1_cnt_ <= upper_bound);
    }

#ifdef PRINT_INTERSECTION
    if (!in_build_) {
        init_totalA = totalA;
        init_totalB = totalB;
        init_totalLenRare = totalLenRare;
        init_totalLenFreq = totalLenFreq;
        init_totalEdgeCase1 = totalEdgeCase1;
        init_totalEdgeCase2 = totalEdgeCase2;
        init_totalGalloping = totalGalloping;
        init_totalShuffling = totalShuffling;
        init_totalCountGalloping = totalCountGalloping;
        init_totalCountShuffling = totalCountShuffling;

        totalEdgeCase1 = totalEdgeCase2 = totalLenRare = totalLenFreq = totalGalloping = totalA = totalB = totalShuffling = totalCountGalloping = totalCountShuffling = 0;
    }
#endif
}

int Alley::DecomposeQuery() {
    if (!valid_)
        return 0;
    return 1;
}

//set y1
bool Alley::GetSubstructure(int subquery_index) {
    if (in_build_) {
		int cur_y1_cnt_ = card_vec_.size();
        if (sample_cnt_ == 0 && cur_y1_cnt_ < max_y1_cnt_) {
            y1_ = recurBuildEntry(); 
			assert(y1_ == 0 || y1_ == 1);
			//stop if mean(y1_) converge
			if (cur_y1_cnt_ >= min_y1_cnt_) {
				double cur_mean = sum_y1_/cur_y1_cnt_;
				double new_mean = (sum_y1_+y1_)/(cur_y1_cnt_+1);
				if (std::abs(cur_mean - new_mean) < 0.05) {
                    sample_cnt_ = -1;
				}
			}
            sum_y1_ += y1_;
            return true;
        }
    }
    else {
        if (sample_cnt_ > 0) {
            y1_ = recurQuery(0); 
            return true;
        }
    }
    return false;
}

//after Run
double Alley::GetFailureRate() {
	if (card_vec_.empty())
		return 1.0;
	size_t size = card_vec_.size();
	return 1.0 - sum_y1_/card_vec_.size();
}

double Alley::EstCard(int subquery_index) {
    return y1_;
}

double Alley::AggCard() {
    if (card_vec_.size() == 0)
        return 0;
    double res = 0.0;
    for (double card : card_vec_) 
        res += card;
    return res / card_vec_.size();
}

double Alley::GetSelectivity() {
    if (!valid_) {
        if (last_ == 0)
            return sample_space_size_[plan_[0]];
        else
            return 0;
    }
    return 1;
}

void Alley::DumpEstimationDetails(std::fstream& fout) {
}

void Alley::PrintStat() {
#ifdef PRINT_INTERSECTION
    if (!in_build_) {
        cout << "init edge case 1, 2 count" << endl;
        cout << init_totalEdgeCase1 << ", " << init_totalEdgeCase2 << endl;
        cout << "init galloping stat: count, avg. |A|, avg. |B|, avg. |C|" << endl;
        cout << init_totalGalloping << ", " << (double)init_totalLenRare/init_totalGalloping << ", " << (double)init_totalLenFreq/init_totalGalloping << ", " << (double)init_totalCountGalloping/init_totalGalloping << endl;
        cout << "init shuffling stat: count, avg. |A|, avg. |B|, avg. |C|" << endl;
        cout << init_totalShuffling << ", " << (double)init_totalA/init_totalShuffling << ", " << (double)init_totalB/init_totalShuffling << ", " << (double)init_totalCountShuffling/init_totalShuffling << endl;

        cout << "edge case 1, 2 count" << endl;
        cout << totalEdgeCase1 << ", " << totalEdgeCase2 << endl;
        cout << "galloping stat: count, avg. |A|, avg. |B|, avg. |C|" << endl;
        cout << totalGalloping << ", " << (double)totalLenRare/totalGalloping << ", " << (double)totalLenFreq/totalGalloping << ", " << (double)totalCountGalloping/totalGalloping << endl;
        cout << "shuffling stat: count, avg. |A|, avg. |B|, avg. |C|" << endl;
        cout << totalShuffling << ", " << (double)totalA/totalShuffling << ", " << (double)totalB/totalShuffling << ", " << (double)totalCountShuffling/totalShuffling << endl;
    }
#endif
}


/************************** INTERNAL FUNCTIONS ************************/

bool Alley::validateBounds() {
    for (int i = 0; i < q->GetNumVertices(); i++) {
        int v = q->GetBound(i);
        if (v != -1) {
            for (auto& edge : q->GetAdj(i)) {
                int u = q->GetBound(edge.id);
                if (u != -1) {
                    if (!g->HasEdge(v, u, edge.label, edge.dir))
                        return false;
                    else 
                        sample_space_size_[edge.id] = 1;
                } else {
                    sample_space_size_[edge.id] = g->GetAdjSize(v, edge.label, edge.dir);
                    if (sample_space_size_[edge.id] == 0)
                        return false;
                }
                backward_[edge.id] |= (1 << i);
                initial_backward_[edge.id] |= (1 << i);
                num_adj_dvids_[edge.id]++;
            }
            sample_space_size_[i] = 1;
            sample_[i] = v;
            rbi_[i] = FIXED;
        } else {
            candidates_.push_back(i);
        }
    }

    last_ = candidates_.size() - 1;

    return true;
}

void Alley::getStarSelectivities() {
    star_sel_.clear();
    star_sel_.resize(q->GetNumVertices());
    edge_sel_.clear();
    edge_sel_.resize(q->GetNumVertices());

    for (int u : candidates_) {
        int min_sel = numeric_limits<int>::max();
        int min_sel_e = numeric_limits<int>::max();
        assert(q->GetBound(u) == -1);
        int vl = q->GetVLabel(u);
        if (vl != -1) {
            min_sel = g->GetNumVertices(vl);
        }
        for (auto& e : q->GetAdj(u)) {
            int sel = g->GetNumVertices(e.label, e.dir);
            int sel_e = g->GetNumEdges(e.label);
            if (sel < min_sel)
                min_sel = sel;
            if (sel_e < min_sel_e)
                min_sel_e = sel_e;
        }
        star_sel_[u] = min_sel;
        edge_sel_[u] = min_sel_e;
    }
}

void Alley::determineFirst() {
    //set plan_[0] and its neighbors' backwards
    int u = -1;
    int max_num = 0;
    bool has_dvid = false;
    for (auto i : candidates_) {
        if (num_adj_dvids_[i] > 0)
            has_dvid = true;
        if (num_adj_dvids_[i] > max_num) {
            max_num = num_adj_dvids_[i];
            u = i;
        }
    }
    if (u < 0) {
        int num_edges = q->GetNumEdges();
        if (num_edges > 1) {
            int num_vls = 0;
            for (auto i : candidates_)
                if (q->GetVLabel(i) != -1)
                    num_vls++;
            if (num_vls == 0) {
                for (auto i : candidates_) {
                    if (q->GetAdj(i).size() == num_edges) {
                        u = i;
                        break;
                    }
                }
            }
        }
        if (u < 0) {
            int min_sel = std::numeric_limits<int>::max(); 
            for (auto i : candidates_) {
                if (q->GetAdjSize(i) > 1 || q->GetVLabel(i) != -1) {
                    if (star_sel_[i] < min_sel) {
                        min_sel = star_sel_[i];
                        u = i;
                    }
                    else if (star_sel_[i] == min_sel) {
                        if (edge_sel_[i] < edge_sel_[u])
                            u = i;
                        else if (edge_sel_[i] == edge_sel_[u])
                            if (q->GetAdjSize(u) < q->GetAdjSize(i))
                                u = i;
                    }
                }
            }
        }
    }

    assert(u >= 0);
    plan_[0] = u;
    for (auto& edge : q->GetAdj(u)) {
        if (q->GetBound(edge.id) == -1) {
            backward_[edge.id] |= (1 << u);
        }
    }
}

void Alley::determineOrder(const int start) {
    vector<bool> visited(q->GetNumVertices(), false);

    int u;
    if (start == -1) {
        determineFirst();
        u = plan_[0];
    } else {
        plan_[0] = u = start;
        for (int i = 0; i < q->GetNumVertices(); i++)
            backward_[i] = initial_backward_[i];
        for (auto& edge : q->GetAdj(u))
            if (q->GetBound(edge.id) == -1)
                backward_[edge.id] |= (1 << u);
    }
    
    if (in_build_)
        use_domain_[u] = true;
    
    visited[u] = true;

    for (int i = 1; i < candidates_.size(); i++) {
        int min_u = -1;
        int max_cycle = 1;
        int min_sel = std::numeric_limits<int>::max();
        for (int u = 0; u < q->GetNumVertices(); u++) {
            if (!visited[u] && q->GetBound(u) == -1) {
                int cycle = __builtin_popcount(backward_[u]);
                if (cycle > max_cycle) {
                    min_u = u;
                    max_cycle = cycle;
                    min_sel = star_sel_[u];
                } else if (cycle == max_cycle) {
                    if (star_sel_[u] < min_sel) {
                        min_u = u;
                        min_sel = star_sel_[u];
                    }
                }
            }
        }
        assert(min_u != -1);

        plan_[i] = min_u;
        visited[min_u] = true;

        for (auto& edge : q->GetAdj(min_u)) {
            if (!visited[edge.id] && q->GetBound(edge.id) == -1) {
                backward_[edge.id] |= (1 << min_u);
            }
        }
    }

    makeRBIGraph();

    if (pattern_mining_ && !in_build_) {
        //has bound vertices
        if (candidates_.size() != q->GetNumVertices()) {
            int visited = 0;
            connected_plan_.clear();
            for (int i = 0; i < candidates_.size(); i++) {
                int u = plan_[i];
                assert(u != -1);
                visited |= (1 << u);
                bool connected = false;
                for (auto& edge : q->GetAdj(u)) {
                    if (visited & (1 << edge.id)) {
                        connected = true;
                        break;
                    }
                }
                if (connected) 
                    connected_plan_.push_back(u);
                else {
                    bool connected = false;
                    for (int v = 0; v < q->GetNumVertices(); v++) {
                        if (!(visited & (1 << v)) && q->GetBound(v) != -1) {
                            for (auto& edge : q->GetAdj(v)) {
                                if (visited & (1 << edge.id)) {
                                    connected = true;
                                    break;
                                }
                            }
                            if (connected) {
                                visited |= (1 << v);
                                connected_plan_.push_back(v);
                            }
                        }
                        if (connected)
                            break;
                    }
                    if (connected)
                        i--;
                }
            }
        }
        else
            connected_plan_ = plan_;
        assert(connected_plan_.size() == q->GetNumVertices());
    }
}


void Alley::getDomainsStatic() {
    if (!in_build_) {
        bool found = false;
        vector<PatternMining*> miners = miner_->GetMiners();
        int maxmaxL=0, minmaxL=MAX_INT;
        for (PatternMining* miner : miners) {
            maxmaxL=max(maxmaxL, miner->getMaxL());
            minmaxL=min(minmaxL, miner->getMaxL());
        }
        for (PatternMining* miner : miners) {
            Pattern pattern;
            vector<int> empty;

            if(miner->GetLTypeInChar() == 'e') {
                Pattern q2p(*q, miner->getLabel2Order(), empty);
                pattern=q2p;
            }
            else if(miner->GetLTypeInChar() == 'v') {
                Pattern q2p(*q, empty, miner->getLabel2Order());
                pattern=q2p;
            }
            else
                assert(false);
            VID_MAP nec_map;
            //nec
            if (makeNEC(pattern, nec_map)) {
                vector<vector<int>> nec2org(pattern.getNumVertices());
                vector<bool> visited(pattern.getNumVertices());
                vector<int> order;
                for (int o = 0; o < connected_plan_.size(); o++) {
                    int nec = nec_map[connected_plan_[o]];
                    if (!visited[nec]) {
                        order.push_back(nec);
                        visited[nec] = true;
                    }
                    nec2org[nec].push_back(connected_plan_[o]);
                }

                if (pattern.getNumVertices() > 2) {
                    vector<bool> nec_has_dom;
                    vector<range> nec_dom;
                    vector<set<pair<int, bool>>> nec_considered_rels;
                    if (getDomainsRecursive(miners, pattern, order, minmaxL, maxmaxL, nec_dom, nec_has_dom, nec_considered_rels, 0)) {
                        if (domains_.empty()) {
                            for (int u = 0; u < q->GetNumVertices(); u++) {
                                int nec = nec_map[u];
                                domains_.push_back(nec_dom[nec]);
                                use_domain_.push_back(nec_has_dom[nec]);
                                dom_considered_rels_.push_back(nec_considered_rels[nec]);
                            }
                        }
                        else {
                            for (int u = 0; u < q->GetNumVertices(); u++) {
                                int nec = nec_map[u];
                                if (nec_dom[nec].end - nec_dom[nec].begin < domains_[u].end - domains_[u].begin) {
                                    domains_[u] = nec_dom[nec];
                                    use_domain_[u] = nec_has_dom[nec];
                                    dom_considered_rels_[u] = nec_considered_rels[nec];
                                }
                            }
                        }
                        found = true;
                    }
                    else {
                    }
                }
            }
            else {
                if (!getDomainsRecursive(miners, pattern, connected_plan_, minmaxL, maxmaxL, domains_, use_domain_, dom_considered_rels_, 0)) {
                }
                else
                    found = true;
            }
            if (domains_.empty()) {
                domains_.resize(q->GetNumVertices());
                use_domain_.resize(q->GetNumVertices(), false);
                dom_considered_rels_.resize(q->GetNumVertices());
            }
            else {
            }
            assert(domains_.size() == q->GetNumVertices());
            assert(use_domain_.size() == q->GetNumVertices());
            assert(dom_considered_rels_.size() == q->GetNumVertices());
            break;
        }
    }
}

void Alley::mergeIntersectionOneVertex(int u) {
    reduced_rel_[u].clear();
    reduced_adj_[u].clear();

    map<pair<int, bool>, set<int>> adjs;
    set<pair<int, bool>> rels;

    for (auto& edge : q->GetAdj(u)) {
        if (backward_[u] & (1 << edge.id)) {
            adjs[make_pair(edge.label, edge.dir)].insert(edge.id);
        }
        else {
			if (pattern_mining_ && use_domain_[u]) {
				pair<int, bool> key(edge.label, edge.dir);
				auto it = dom_considered_rels_[u].find(key);
				if (it == dom_considered_rels_[u].end()) {
					rels.insert(make_pair(edge.label, edge.dir));
				}
				else {
				}
			}
			else
				rels.insert(make_pair(edge.label, edge.dir));
        }
    }
    for (auto& p : rels) {
        //if this rel is not subsumed by any adjs 
        if (adjs[p].size() == 0) {
            reduced_rel_[u].emplace(p.first, p.second);
        }
    }
    for (auto& p : adjs) {
        for (int id : p.second) {
            reduced_adj_[u].emplace_back(id, p.first.first, !p.first.second);
        }
    }
}

void Alley::mergeIntersection() {
    for (auto u : candidates_) {
        mergeIntersectionOneVertex(u);
    }
}

void Alley::makeRBIGraph() {
    rbi_[plan_[0]] = RED;
    for (int i = 1; i < candidates_.size(); i++) {
        int u = plan_[i];
        int num_red = 0;
        for (auto& edge : q->GetAdj(u)) {
            if (rbi_[edge.id] == RED || rbi_[edge.id] == FIXED) {
                num_red++;
            }
        }
        if (num_red == q->GetAdjSize(u)) {
            if (num_red == 1) 
                rbi_[u] = BLACK;
            else
                rbi_[u] = IVORY;
        } else {
            rbi_[u] = RED;
        }
    }
}

//return 0 or 1
double Alley::recurBuildEntry() {
    int u = plan_[0];
    int size = sample_space_size_[u];

    int index = std::uniform_int_distribution<>(0, size-1)(generator_);

    int v = sample_space_[0][index];
    sample_[u] = v;

    getLocalSampleSpaceBuild(1);
    cartesian_ = sample_space_size_[plan_[1]]; 
    double ret = recurBuild(1);

    sample_[u] = -1;
    return ret;
}

//return 0 or 1
double Alley::recurBuild(const int order) {
    recur_call_cnt_++;

    int u = plan_[order];
    int size = sample_space_size_[u];

    if (size == 0) {
        return 0;
    }

    if (rbi_[u] != RED) {
        if (size == 0 || last_ == order) {
            return size > 0 ? 1 : 0;
        }
        getLocalSampleSpaceBuild(order + 1);
        return recurBuild(order + 1);
    }

    double ret = 0.0;
    int skipped = 0;
    int branch = std::max(size >> 5, std::min(size, 8));

    unordered_set<int> indices;
    if (size == branch) {
        vector<int> temp; 
        for (int i = 0; i < size; i++)
            temp.push_back(i);
        std::random_shuffle(temp.begin(), temp.end());
        std::copy(temp.begin(), temp.end(), std::inserter(indices, indices.end())); 
    } else {
        for (int r = size - branch; r < size; r++) {
            int v = std::uniform_int_distribution<>(0, r)(generator_);
            if (!indices.insert(v).second)
                indices.insert(r);
        }
    }
    int i = 0;
    auto prev_cartesian = cartesian_;
    for (auto it = indices.begin(); it != indices.end(); ++it, ++i) {
        int v = sample_space_[order][*it];
        sample_[u] = v;

        getLocalSampleSpaceBuild(order + 1);
        cartesian_ = prev_cartesian * sample_space_size_[plan_[order+1]];
        if (cartesian_ > CARTESIAN_THRESHOLD) {
            cartesian_ = prev_cartesian;
            sample_[u] = -1;
            return 0;
        }
        if (recurBuild(order + 1) > 0) {
            cartesian_ = prev_cartesian;
            sample_[u] = -1;
            return 1;
        }
    }

    cartesian_ = prev_cartesian;
    sample_[u] = -1;
    return 0; 
}

double Alley::recurQuery(const int order) {
    if (sample_cnt_ <= 0)
        return 0.0;
    int u = plan_[order];
    int size = sample_space_size_[u];

    if (size == 0) {
        sample_cnt_--;
        return 0.0;
    }

    if (rbi_[u] != RED) {
        if (size == 0 || last_ == order) {
            sample_cnt_--;
            return size;
        }
        getLocalSampleSpaceQuery(order + 1);
        double ret = recurQuery(order + 1);
        return ret * size;
    }

    double ret = 0.0;
    int skipped = 0;
    int branch = std::max(size >> 5, std::min(size, 8));

    unordered_set<int> indices;
    if (size == branch) {
        vector<int> temp; 
        for (int i = 0; i < size; i++)
            temp.push_back(i);
        std::random_shuffle(temp.begin(), temp.end());
        std::copy(temp.begin(), temp.end(), std::inserter(indices, indices.end())); 
    } else {
        for (int r = size - branch; r < size; r++) {
            int v = std::uniform_int_distribution<>(0, r)(generator_);
            if (!indices.insert(v).second)
                indices.insert(r);
        }
    }
    int i = 0;

    for (auto it = indices.begin(); it != indices.end(); ++it, ++i) {
        if (sample_cnt_ <= 0) {
            skipped = branch - i;
            break;
        }
        int v = sample_space_[order][*it];
        sample_[u] = v;

        getLocalSampleSpaceQuery(order + 1);
        ret += recurQuery(order + 1);
    }

    sample_[u] = -1;
    return ret * size / (branch - skipped);
}

void Alley::getLocalSampleSpaceQuery(const int order) {
    //set sample_space_, sample_space_size_ for plan_[order]
    const int u = plan_[order];
    std::set<range> ranges;

    int vl = q->GetVLabel(u);
    if (vl != -1) {
        range r = g->GetVertices(vl);
        ranges.insert(r);
    }

    for (auto& edge : reduced_adj_[u]) {
        range r = g->GetAdj(sample_[edge.id], edge.label, edge.dir);
        ranges.insert(r);
    }

	if (pattern_mining_) {
		if (use_domain_[u]) {
			ranges.insert(domains_[u]);
        }
    }
	
    for (auto& edge : reduced_rel_[u]) {
        range r = g->GetVertices(edge.first, edge.second);
        ranges.insert(r);
    }

    if (ranges.size() > 1) {
        int min_size = ranges.begin()->end - ranges.begin()->begin;
        if (min_size > intersect_capacity_[order]) {
            intersect_capacity_[order] = std::max(min_size, 2 * intersect_capacity_[order]);
            int* b = intersect_[order];
            int* t = temp_[order];
            intersect_[order] = new int[intersect_capacity_[order]];
            temp_[order] = new int[intersect_capacity_[order]];
            delete [] b;
            delete [] t;
        }
        sample_space_size_[u] = intersect(ranges, intersect_[order], temp_[order]);
        sample_space_[order] = intersect_[order];
    } else {
        sample_space_size_[u] = ranges.begin()->end - ranges.begin()->begin;
        sample_space_[order] = const_cast<int*>(ranges.begin()->begin);
    }

    review_sum += sample_space_size_[u];
    review_cnt ++;
}

void Alley::getLocalSampleSpaceBuild(const int order) {
    //set sample_space_, sample_space_size_ for plan_[order]
    const int u = plan_[order];
    std::set<range> ranges;

    for (auto& edge : reduced_adj_[u]) {
        range r = g->GetAdj(sample_[edge.id], edge.label, edge.dir);
        ranges.insert(r);
    }

    //for VERTEX_LABEL, vertex label of u is always considered in domains
    assert(pattern_mining_);
    if (use_domain_[u]) {
        ranges.insert(domains_[u]);
        //if use_domain = true, ltype = VERTEX_LABEL, and label is set, the label should have been considered in domain 
    }
    else if (q->GetVLabel(u) != -1) {
        ranges.insert(g->GetVertices(q->GetVLabel(u)));
    }
	
    assert(!ranges.empty());

    for (auto& edge : reduced_rel_[u]) {
        range r = g->GetVertices(edge.first, edge.second);
        ranges.insert(r);
    }
    
    if (ranges.size() > 1) {
        int min_size = ranges.begin()->end - ranges.begin()->begin;
        bool use_cache = min_size > CACHE_THRESHOLD;
        if (use_cache) {
            //find in cache
            if (searchCache(ranges, &sample_space_[order], &sample_space_size_[u]))
                return;
        }

        if (min_size > intersect_capacity_[order]) {
            intersect_capacity_[order] = std::max(min_size, 2 * intersect_capacity_[order]);
            int* b = intersect_[order];
            int* t = temp_[order];
            intersect_[order] = new int[intersect_capacity_[order]];
            temp_[order] = new int[intersect_capacity_[order]];
            delete [] b;
            delete [] t;
        }
        sample_space_size_[u] = intersect(ranges, intersect_[order], temp_[order]);
        sample_space_[order] = intersect_[order];

        //store in cache
        if (use_cache) {
            storeCache(ranges, sample_space_[order], sample_space_size_[u]);
        }
    } else {
        sample_space_size_[u] = ranges.begin()->end - ranges.begin()->begin;
        sample_space_[order] = const_cast<int*>(ranges.begin()->begin);
    }
}

bool Alley::searchCache(set<range>& ranges, int** arr, int* size) {
    size_t hv = 0;
    for (auto& r : ranges) {
        boost::hash_combine(hv, r.begin);
        boost::hash_combine(hv, r.end);
    }
    auto it = intersect_cache_.find(hv);
    if (it != intersect_cache_.end()) {
        if (it->second.type == EMPTY) {
            *arr  = NULL;
            *size = 0;
        }
        else if (it->second.type == MIN) {
            *arr  = const_cast<int*>(ranges.begin()->begin);
            *size = ranges.begin()->end - ranges.begin()->begin; 
        }
        else {
            *arr  = it->second.arr;
            *size = it->second.size;
        }
        return true;
    }
    return false;
}

void Alley::storeCache(set<range>& ranges, const int* arr, int size) {
    size_t hv = 0;
    for (auto& r : ranges) {
        boost::hash_combine(hv, r.begin);
        boost::hash_combine(hv, r.end);
    }
    cache_val val;
    val.arr = NULL;
    if (size > 0) {
        int min_size = ranges.begin()->end - ranges.begin()->begin;
        if (min_size == size)
            val.type = MIN;
        else {
            assert(min_size > size);
            val.type = MAT;
            val.arr = new int[size];
            memcpy(val.arr, arr, sizeof(int) * size);
            val.size = size;
        }
    }
    else {
        assert(size == 0);
        val.type = EMPTY;
    }
    intersect_cache_[hv] = val;
}

/************************** MINING FUNCTIONS ************************/

void Alley::SetMode(bool b, int maxL) {
    //initialize
    if (b) {
        initial_search_space_.resize(maxL+1);
        initial_search_space_size_.resize(maxL+1);

        initial_intersect_capacity_ = new int[maxL+1];
        initial_intersect_ = new int*[maxL+1];
        initial_temp_      = new int*[maxL+1];

        for (int u = 0; u < maxL+1; u++) {
            initial_intersect_capacity_[u] = INITIAL_CAPACITY;
            initial_intersect_[u] = new int[INITIAL_CAPACITY];
            initial_temp_[u]      = new int[INITIAL_CAPACITY];
        }

        marked_.clear();
        pruned_.clear();
        Dom_.clear();
        Pr_.clear();
        Dom_.resize(maxL+1);
        Pr_.resize(maxL+1);

        for (int u = 0; u < maxL+1; u++) {
            Dom_[u].reserve(INITIAL_CAPACITY);
            Pr_[u].reserve(INITIAL_CAPACITY);
        }
    }
    in_build_ = b;
    maxL_ = maxL;
}

void Alley::SetData(DataGraph* gin) {
    g = gin;
    marked_.resize(g->GetNumVertices());
    pruned_.resize(g->GetNumVertices());
}

void Alley::ClearCache() {
    for (auto& p : intersect_cache_) {
        if (p.second.arr)
            delete [] p.second.arr;
    }
    intersect_cache_.clear();
    review_sum = review_cnt = 0;
}

void Alley::SetDomains(vector<range>& d, vector<bool>& u, vector<star>& r) {
    domains_ = d;
    //set true iff domain is not from the intersection for a star
    //but a more complex structure
    use_domain_ = u;
    dom_considered_rels_ = r;
}

bool Alley::calculateDomainsFast(Pattern& pattern, vector<range>& domains, vector<bool>& use_domain, vector<star>& rels, bool print, bool& fast_timeout) {
	bool timeout = false;

    clock_t start, end;
    start = clock();

    for (int u = 0; u < pattern.getNumVertices(); u++) {
        Dom_[u].clear();
        Pr_[u].clear();
        assert(Dom_[u].capacity() >= INITIAL_CAPACITY);
        assert(Pr_[u].capacity() >= INITIAL_CAPACITY);
    }
    bool res = treeCycleInternal(pattern, domains, use_domain, rels, 0, fast_timeout, print);

    end = clock();
    double t = (double)(end - start)/CLOCKS_PER_SEC;
    calculate_domains_time_ += t;
    calculate_domains_cnt_++;
    //cout << "calculating time " << t << endl;

    if (fast_timeout) {
        return true;
    }

    if (!res) {
        return false;
    }

    review_cnt++;
    for (int u = 0; u < pattern.getNumVertices(); u++) {
        review_sum += Dom_[u].size(); 
    }

    return true;
}

bool Alley::treeCycleInternal(Pattern& pattern, vector<range>& domains, vector<bool>& use_domain, vector<star>& rels, const int level, bool& timeout, bool print) {
    vector<int> deg;
    pattern.getCoreDegree(deg);

    if (level == 0) {
        assert(plan_.size() == pattern.getNumVertices());
        if (sample_space_size_[plan_[0]] == 0) {
            return false;
        }
        //Init() has been done
        //plan_, backward_, rbi_, sample_sapce_[0], 
        //sample_space_size_[plan_[0]],
        //reduced_rel_, reduced_adj_ are set already
        
        //need to set label_, subplans_, match_, CT_, ap_, to_clear_,
        //black_child_, is_black_marked_, is_black_calculated_, single_star_marked_, is_star_marked_
        
        for (int i = 0; i < plan_.size(); i++) {
            int u = plan_[i];
            set<range> ranges;

            //label, forward
            if (q->GetVLabel(u) != -1) {
                range r = g->GetVertices(q->GetVLabel(u));
                ranges.insert(r);
            }
            for (auto& edge : reduced_rel_[u]) {
                range r = g->GetVertices(edge.first, edge.second);
                ranges.insert(r);
            }
            for (auto& edge : reduced_adj_[u]) {
                range r = g->GetVertices(edge.label, !edge.dir);
                ranges.insert(r);
            }

            if (use_domain[u] || u == plan_[0]) {
                //label, dir
                for (auto& contained : rels[u]) {
                    reduced_rel_[u].erase(contained);
                    range r = g->GetVertices(contained.first, contained.second);
                    ranges.erase(r);
                }
                ranges.insert(domains[u]);
            }
            
            getInitialSearchSpace(u, ranges);

            if (initial_search_space_size_[u] == 0) {
                return false;
            }
            
            //now, no need to use domains[u] nor reduced_rel_[u]
        }

        std::sort(plan_.begin(), plan_.end(), [&](const int& a, const int& b) -> bool {
                if (deg[a] != deg[b])
                return deg[a] > deg[b];
                else if (pattern.getDegree(a) == 1 || pattern.getDegree(b) == 1)
                return pattern.getDegree(a) > pattern.getDegree(b);
                else {
                //prefer small space
                return initial_search_space_size_[a] < initial_search_space_size_[b];
                }
                });

        rank_.clear();
        rank_.resize(plan_.size());
        for (int r = 0; r < plan_.size(); r++) 
            rank_[plan_[r]] = r;
    }
    else {
        int min_u = -1;
        int min_rank = MAX_INT;
        int max_rank = MIN_INT;
        for (int u : subplans_[plan_[0]]) {
            if (rank_[u] < min_rank) {
                min_u = u;
                min_rank = rank_[u];
            }
            if (rank_[u] > max_rank) 
                max_rank = rank_[u];
        }
        rank_[plan_[0]] = max_rank + 1;
    }

    //reorder!
    plan_.clear();
    backward_.clear();
    forward_.clear();
    exact_backward_.clear();
    initial_backward_.clear();
    match_.clear();
    rbi_.clear();
    reduced_rel_.clear();
    reduced_adj_.clear();
    CT_.clear();

    plan_.resize(q->GetNumVertices(), -1);
    backward_.resize(q->GetNumVertices(), 0);
    forward_.resize(q->GetNumVertices());
    exact_backward_.resize(q->GetNumVertices());
    initial_backward_.resize(q->GetNumVertices(), 0);
    match_.resize(q->GetNumVertices(), -1);
    rbi_.resize(q->GetNumVertices());
    reduced_rel_.resize(q->GetNumVertices());
    reduced_adj_.resize(q->GetNumVertices());
    CT_.resize(q->GetNumVertices());

    //using rank_, set plan_, forward_, backward_, exact_backward_, rbi_, CT_
    determineOrderByRank();
    makeCoreFront();
    determineUseInitialSearchSpace();

    //set reduced_adj_ (used), reduced_rel_ (not used)
    mergeIntersection();

    //finally, set ap_, subplans_, to_clear_,
    //black_child_, is_black_marked_
    
    subplans_.clear();
    ap_.clear();
    to_clear_.clear();
    is_calculated_.clear();
    is_black_marked_.clear();

    search_space_.clear();
    space_size_.clear();

    subplans_.resize(plan_.size());
    ap_.resize(plan_.size(), false);
    to_clear_.resize(plan_.size());
    is_calculated_.resize(plan_.size(), false);
    is_black_marked_.resize(plan_.size(), false);
    
    search_space_.resize(plan_.size(), NULL);
    space_size_.resize(plan_.size(), 0);
    
    pattern.findAPs(plan_, ap_, subplans_);
        
	for (int u = 0; u < plan_.size(); u++) {
        if (ap_[u])
            assert(!subplans_[u].empty());
        else
            assert(u == plan_[0] || subplans_[u].empty());
    }

	//reorder subplans
    vector<int> vid2order(plan_.size(), -1); 
	for (int i = 0; i < plan_.size(); i++)
        vid2order[plan_[i]] = i;

	vector<char> visited(plan_.size(), 0);
	for (int u = 0; u < plan_.size(); u++) {
		std::sort(subplans_[u].begin(), subplans_[u].end(), 
				[&](const int& a, const int& b) -> bool {
				    return vid2order[a] < vid2order[b]; 
				});
		for (int v : subplans_[u]) {
			visited[v]++;
        }
	}
	assert(visited[plan_[0]] == 0);
	for (int i = 1; i < plan_.size(); i++)
		assert(visited[plan_[i]] == 1);

    rootplan_.clear();
    //if root is ap & core, should divide in half
    if (ap_[plan_[0]] && deg[plan_[0]] > 1) {
        Pattern p(pattern);
        p.deleteVertexWithoutRenumber(plan_[0]);
        int next = subplans_[plan_[0]][0];
        for (int i = 1; i < subplans_[plan_[0]].size(); i++) {
            //if next and u are disconnected & u is RED & u is ap, move u to rootplan_ 
            //if u is BLACK/IVORY, keep it in subplans_[root]
            //if u is RED but not ap, keep it in subplans_[root]
            int u = subplans_[plan_[0]][i];
            if (!p.isConnected(next, u) && rbi_[u] == RED && ap_[u]) {
                rootplan_.push_back(u);
                subplans_[plan_[0]].erase(subplans_[plan_[0]].begin() + i);
                i--;
            }
        }
    }

#ifdef PRINT_ALLEY
    printPlan();
    cout << "printing subplans ";
    for (auto& s : subplans_) {
        cout << "{";
        for (auto& v : s)
            cout << v << ",";
        cout << "}, "; 
    }
    cout << endl;
    cout << "printing rootplan ";
    for (auto& v : rootplan_)
        cout << v << ",";
    cout << endl;
#endif
    
    int rollback = -1;
    
    const int u = plan_[0];
    int size = space_size_[u] = initial_search_space_size_[u];
    const int* cand = search_space_[u] = initial_search_space_[u];
    assert(cand != NULL);
    assert(size > 0);

    clock_t start = clock();
    bool res = false;
    cartesian_ = 1;

    for (int i = 0; i < size; i++) {
        if (level > 1 && marked_.get(cand[i], u))
            continue;
        match_[u] = cand[i];

        for (int b : forward_[u])
            is_calculated_[b] = is_black_marked_[b] = false;
        if (treeCycleRecur(u, 0, rollback, domains, use_domain, print, level > 1, timeout)) {
            //due to rankUp
            if (!marked_.get(match_[u], u)) {
                marked_.set(match_[u], u);
                Dom_[u].push_back(match_[u]);
            }

            assert(rollback == -1);
            res = true;
            if (timeout)
                break;
        }
        else {
            if (!pruned_.get(match_[u], u)) {
                pruned_.set(match_[u], u);
                Pr_[u].push_back(match_[u]);
            }
        }
        if (rollback == u)
            rollback = -1;
        else 
            assert(rollback == -1);
        if (!to_clear_[u].empty()) {
            for (int c : to_clear_[u]) {
                for (int v : Pr_[c]) {
                    assert(pruned_.get(v, c));
                    pruned_.unset(v, c);
                }
                Pr_[c].clear();
            }
            to_clear_[u].clear();
        }
        if (!subplans_[u].empty()) {
            if (level <= 1) {
                clock_t end = clock(); //bottleneck
                if ((double)(end - start)/CLOCKS_PER_SEC > (5 * (level + 1))) {
                    match_[u] = -1;
                    assert(u == plan_[0]);

                    if (level == 0) {
                        if (treeCycleInternal(pattern, domains, use_domain, rels, level + 1, timeout, print))
                            res = true;
                    }
                    else {
                        //return timeout
                        timeout = true;
                        break;

                        size_t num_rotation = subplans_[u].size() + 1;
                        for (size_t i = 0; i < num_rotation; i++) {
                            printPlan();
                            if (treeCycleInternal(pattern, domains, use_domain, rels, level + 1, timeout, print))
                                res = true;
                        }
                    }
                    break;
                }
            }
        }
    }
    match_[u] = -1; 

    if (level == 0) {
        //clear marked_ and pruned_
        for (int u = 0; u < pattern.getNumVertices(); u++) {
            //sort
            std::sort(Dom_[u].begin(), Dom_[u].end());

            for (int v : Dom_[u]) {
                assert(marked_.get(v, u));
                marked_.unset(v, u);
            }
            for (int v : Pr_[u]) {
                assert(pruned_.get(v, u));
                pruned_.unset(v, u);
            }
        }
    }
        
#ifdef PRINT_ALLEY
    if (print)
        cout << "treeCycleInternal level " << level << " done" << endl;
#endif

    return res;
}

//if rankUp = true, find just any one match
bool Alley::treeCycleRecur(const int head, const int index, int& rollback, vector<range>& domains, vector<bool>& use_domain, bool print, bool find_one, bool& timeout) {
    if (index < subplans_[head].size()) {
        const int u = subplans_[head][index];
        
        if (is_calculated_[u]) {
            assert(space_size_[u] > 0);
        }
        else {
            getLocalSearchSpace(u);
            is_calculated_[u] = true;

            if (space_size_[u] == 0) {
                rollback = CT_[u][0];
                if (!pruned_.get(match_[rollback], rollback)) {
                    pruned_.set(match_[rollback], rollback);
                    Pr_[rollback].push_back(match_[rollback]);
                }
                if (CT_[u].size() > 1) {
                    to_clear_[CT_[u][1]].push_back(rollback);
                }
                return false;
            }
        }

        const int* cand = search_space_[u];
        int size = space_size_[u];

        assert(cand != NULL);
        if (rbi_[u] != RED) {
            return treeCycleRecur(head, index + 1, rollback, domains, use_domain, print, find_one, timeout);
        }

        cartesian_ *= space_size_[u];
        if (cartesian_ > CARTESIAN_THRESHOLD) {
            timeout = true;
            return true;
        }

        bool ret = false;
        for (int i = 0; i < size; i++) {
            if (pruned_.get(cand[i], u)) {
                continue;
            }
            match_[u] = cand[i];
            for (int b : forward_[u])
                is_calculated_[b] = is_black_marked_[b] = false;
            //found a match!
            if (treeCycleRecur(head, index + 1, rollback, domains, use_domain, print, find_one, timeout)) {
                if (timeout)
                    return true;
                if (find_one) {
                    match_[u] = -1;
                    if (!to_clear_[u].empty()) {
                        for (int c : to_clear_[u]) {
                            for (int v : Pr_[c]) {
                                assert(pruned_.get(v, c));
                                pruned_.unset(v, c);
                            }
                            Pr_[c].clear();
                        }
                        to_clear_[u].clear();
                    }
                    cartesian_ /= space_size_[u];
                    return true;
                }
                assert(rollback == -1);
                ret = true;
            }
            else {
                if (rollback != -1) {
                    if (rollback == u) {
                        rollback = -1;
                    }
                    else {
                        //skip this query vertex
                        match_[u] = -1;
                        if (!to_clear_[u].empty()) {
                            for (int c : to_clear_[u]) {
                                for (int v : Pr_[c]) {
                                    assert(pruned_.get(v, c));
                                    pruned_.unset(v, c);
                                }
                                Pr_[c].clear();
                            }
                            to_clear_[u].clear();
                        }
                        cartesian_ /= space_size_[u];
                        return false;
                    }
                }
            }
            if (!to_clear_[u].empty()) {
                for (int c : to_clear_[u]) {
                    for (int v : Pr_[c]) {
                        assert(pruned_.get(v, c));
                        pruned_.unset(v, c);
                    }
                    Pr_[c].clear();
                }
                to_clear_[u].clear();
            }
        }
        match_[u] = -1; 
        cartesian_ /= space_size_[u];
        return ret;
    }
    else {
        assert(index == subplans_[head].size());
        bool ret = true;
        //note that, match_[head] is matched to head for the first time
        //validate match_[head]
        if (head == plan_[0] && ap_[head]) {
            for (int u : rootplan_) {
                getLocalSearchSpace(u);
                is_calculated_[u] = true;

                if (space_size_[u] == 0) {
                    rollback = CT_[u][0];
                    if (!pruned_.get(match_[rollback], rollback)) {
                        pruned_.set(match_[rollback], rollback);
                        Pr_[rollback].push_back(match_[rollback]);
                    }
                    if (CT_[u].size() > 1) {
                        to_clear_[CT_[u][1]].push_back(rollback);
                    }
                    ret = false;
                    break;
                }
                
                const int* cand = search_space_[u];
                int size = space_size_[u];
                assert(cand != NULL);
                
                ret = false;
                for (int i = 0; i < size; i++) {
                    if (marked_.get(cand[i], u)) {
                        ret = true;
                        continue;
                    }
                    if (pruned_.get(cand[i], u)) {
                        continue;
                    }
                    match_[u] = cand[i];
                    for (int b : forward_[u])
                        is_calculated_[b] = is_black_marked_[b] = false;
                    unsigned long long prev_cartesian = cartesian_;
                    cartesian_ = 1;
                    //find all matches
                    if (treeCycleRecur(u, 0, rollback, domains, use_domain, print, false, timeout)) {
                        if (timeout)
                            return true;
                        assert(rollback == -1);
                        assert(!marked_.get(match_[u], u));
                        marked_.set(match_[u], u);
                        Dom_[u].push_back(match_[u]);
                        ret = true;
                    }
                    else {
                        //can be set true by rollback
                        if (!pruned_.get(match_[u], u)) {
                            pruned_.set(match_[u], u);
                            Pr_[u].push_back(match_[u]);
                        }

                        if (rollback == u) {
                            rollback = -1;
                        }
                        else
                            assert(rollback == -1);
                    }
                    cartesian_ = prev_cartesian;
                }
                match_[u] = -1; 
                if (!ret)
                    break;
            }
        }
        if (!ret)
            return ret;
        for (int u : subplans_[head]) {
            if (ap_[u] && rbi_[u] == RED) {

                assert(match_[u] != -1);
                if (marked_.get(match_[u], u)) {
                    assert(!pruned_.get(match_[u], u));
                    continue;
                }
                else if (pruned_.get(match_[u], u)) {
                    assert(!marked_.get(match_[u], u));
                    ret = false;
                    break;
                }
                else {
                    assert(!marked_.get(match_[u], u));
                    assert(!pruned_.get(match_[u], u));
                    for (int b : forward_[u])
                        is_calculated_[b] = is_black_marked_[b] = false;
                    unsigned long long prev_cartesian = cartesian_;
                    cartesian_ = 1;
                    //find all matches
                    if (treeCycleRecur(u, 0, rollback, domains, use_domain, print, false, timeout)) {
                        if (timeout)
                            return true;
                        assert(rollback == -1);
                        assert(!marked_.get(match_[u], u));
                        marked_.set(match_[u], u);
                        Dom_[u].push_back(match_[u]);
                        cartesian_ = prev_cartesian;
                    }
                    else {
                        //can be set true by rollback
                        if (!pruned_.get(match_[u], u)) {
                            pruned_.set(match_[u], u);
                            Pr_[u].push_back(match_[u]);
                        }

                        ret = false;
                        if (rollback == u) {
                            rollback = -1;
                        }
                        else
                            assert(rollback == -1);

                        cartesian_ = prev_cartesian;
                        break;
                    }
                }
            }
        }
        if (ret) {
            for (int u : subplans_[head]) {
                assert(u != head);
                assert(match_[head] != -1);
                if (rbi_[u] != RED) {
                    if (rbi_[u] == BLACK) {
                        if (is_black_marked_[u]) {
                            continue;
                        }
                        is_black_marked_[u] = true;
                    }

                    assert(subplans_[u].empty());
                    assert(match_[u] == -1);
                    assert(space_size_[u] > 0);
                    for (int i = 0; i < space_size_[u]; i++) {
                        int v = search_space_[u][i];
                        if (!marked_.get(v, u)) {
                            marked_.set(v, u);
                            Dom_[u].push_back(v);
                        }
                    }
                }
                else if (!ap_[u]) {
                    assert(match_[u] != -1);
                    if (!marked_.get(match_[u], u)) {
                        marked_.set(match_[u], u);
                        Dom_[u].push_back(match_[u]);
                    }
                }
                else {
                    assert(rbi_[u] == RED && ap_[u]);
                    assert(marked_.get(match_[u], u));
                }
            }
        }
        return ret;
    }
}

//set CT_
void Alley::determineOrderByRank() {
    vector<bool> visited(q->GetNumVertices(), false);

    int u;
    int min_rank = std::numeric_limits<int>::max();
    for (int i = 0; i < q->GetNumVertices(); i++) {
        if (min_rank > rank_[i]) {
            u = i;
            min_rank = rank_[i];
        }
    }
    plan_[0] = u;
    for (auto& edge : q->GetAdj(u)) {
        if (q->GetBound(edge.id) == -1) {
            if (std::find(CT_[edge.id].begin(), CT_[edge.id].end(), u) == CT_[edge.id].end())
                CT_[edge.id].push_back(u);
            backward_[edge.id] |= (1 << u);
            forward_[u].push_back(edge.id);
            exact_backward_[edge.id].emplace_back(u, edge.label, edge.dir); 
        }
    }
    
    visited[u] = true;
    CT_[u].push_back(-1);

    for (int i = 1; i < candidates_.size(); i++) {
        int min_u = -1;
        int max_cycle = 1;
        int min_rank = std::numeric_limits<int>::max();
        for (int u = 0; u < q->GetNumVertices(); u++) {
            if (!visited[u] && q->GetBound(u) == -1) {
                int cycle = __builtin_popcount(backward_[u]);
                if (cycle > max_cycle) {
                    min_u = u;
                    max_cycle = cycle;
                    min_rank = rank_[u]; 
                } else if (cycle == max_cycle) {
                    if (rank_[u] < min_rank) {
                        min_u = u;
                        min_rank = rank_[u];
                    }
                }
            }
        }
        assert(min_u != -1);

        plan_[i] = min_u;
        visited[min_u] = true;

        for (auto& edge : q->GetAdj(min_u)) {
            if (!visited[edge.id] && q->GetBound(edge.id) == -1) {
                if (std::find(CT_[edge.id].begin(), CT_[edge.id].end(), min_u) == CT_[edge.id].end())
                    CT_[edge.id].push_back(min_u);

                backward_[edge.id] |= (1 << min_u);
                forward_[min_u].push_back(edge.id);
                exact_backward_[edge.id].emplace_back(min_u, edge.label, edge.dir); 
            }
        }
    }

    for (u = 0; u < q->GetNumVertices(); u++) {
        std::reverse(CT_[u].begin(), CT_[u].end());
    }

    makeRBIGraph();
}


void Alley::makeCoreFront() {
    for (int i = 0; i < plan_.size(); i++) {
        int u = plan_[i];
        if (rbi_[u] == BLACK) {
            for (int j = i + 1; j < plan_.size(); j++) {
                //swap
                if (rbi_[plan_[j]] != BLACK) {
                    int temp = plan_[i];
                    plan_[i] = plan_[j];
                    plan_[j] = temp;
                    break;
                }

            }
            //not changed
            if (rbi_[plan_[i]] == BLACK) {
                for (int k = i; k < plan_.size(); k++)
                    assert(rbi_[plan_[k]] == BLACK);
                break;
            }
        }
    }
}

void Alley::determineUseInitialSearchSpace() {
    use_initial_search_space_.clear();
    use_initial_search_space_.resize(plan_.size(), false);
    //if IVORY, no need to use domain (selective enough)
    //if BLACK, no need to use domain (determined by parent)
    for (int u = 0; u < plan_.size(); u++)
        if (q->GetVLabel(u) != -1 || (rbi_[u] == RED && (!reduced_rel_[u].empty() || use_domain_[u])))
            use_initial_search_space_[u] = true;
}

void Alley::getInitialSearchSpace(const int u, set<range>& ranges) {
    assert(!ranges.empty());
    if (ranges.size() > 1) {
        int min_size = ranges.begin()->end - ranges.begin()->begin;
        //find in cache
        if (searchCache(ranges, &initial_search_space_[u], &initial_search_space_size_[u]))
            return;

        if (min_size > initial_intersect_capacity_[u]) {
            initial_intersect_capacity_[u] = std::max(min_size, 2 * initial_intersect_capacity_[u]);
            int* b = initial_intersect_[u];
            int* t = initial_temp_[u];
            initial_intersect_[u] = new int[initial_intersect_capacity_[u]];
            initial_temp_[u]      = new int[initial_intersect_capacity_[u]];
            delete [] b;
            delete [] t;
        }
        initial_search_space_size_[u] = intersect(ranges, initial_intersect_[u], initial_temp_[u]);
        initial_search_space_[u] = initial_intersect_[u];

        //store in cache
        storeCache(ranges, initial_search_space_[u], initial_search_space_size_[u]);
    } else {
        initial_search_space_size_[u] = ranges.begin()->end - ranges.begin()->begin;
        initial_search_space_[u] = const_cast<int*>(ranges.begin()->begin);
    }
    for (int i = 0; i < plan_.size(); i++)
        assert(initial_search_space_[u] != intersect_[i]); 
}

void Alley::getLocalSearchSpace(const int u) {
    set<range> ranges;
    for (auto& edge : reduced_adj_[u]) {
        range r = g->GetAdj(match_[edge.id], edge.label, edge.dir);
        ranges.insert(r);
    }
    assert(!ranges.empty());

    //initial_search_space_ points to initial_intersect_, cache, domain, or data graph 
    //if initial_intersect_, do not use cache
    if (use_initial_search_space_[u]) {
        range init;
        init.begin = initial_search_space_[u];
        init.end   = initial_search_space_[u] + initial_search_space_size_[u];
        assert(initial_search_space_size_[u] > 0);
        ranges.insert(init);
    }
    
    if (ranges.size() > 1) {
        int min_size = ranges.begin()->end - ranges.begin()->begin;
        bool use_cache = min_size > CACHE_THRESHOLD;
        if (use_cache) {
            if (use_initial_search_space_[u] && initial_search_space_[u] == initial_intersect_[u])
                use_cache = false;
            else {
                //find in cache
                if (searchCache(ranges, &search_space_[u], &space_size_[u]))
                    return;
            }
        }
        if (min_size > intersect_capacity_[u]) {
            intersect_capacity_[u] = std::max(min_size, 2 * intersect_capacity_[u]);
            int* b = intersect_[u];
            int* t = temp_[u];
            intersect_[u] = new int[intersect_capacity_[u]];
            temp_[u] = new int[intersect_capacity_[u]];
            delete [] b;
            delete [] t;
        }
        space_size_[u] = intersect(ranges, intersect_[u], temp_[u]);
        search_space_[u] = intersect_[u];

        //store in cache
        if (use_cache) {
            storeCache(ranges, search_space_[u], space_size_[u]);
        }
    } else {
        space_size_[u] = ranges.begin()->end - ranges.begin()->begin;
        search_space_[u] = const_cast<int*>(ranges.begin()->begin);
    }
}

void Alley::printPlan() {
    cout << "plan: ";
    for (auto u : plan_)
        cout << u << " ";
    cout << endl;
}



/************************** TEST FUNCTIONS ************************/

bool Alley::calculateDomainsNaive(Pattern& p, vector<set<int>>& Dom) {
    Dom.resize(p.getNumVertices());

    vector<int> match(p.getNumVertices(), -1);
    return naiveRecur(p, match, 0, Dom, clock());
}

bool Alley::naiveRecur(Pattern& p, vector<int>& match, int u, vector<set<int>>& Dom, clock_t start) {
    if (((double)(clock() - start)/CLOCKS_PER_SEC) > 5 * 60)
        return true;
    if (u == match.size()) {
        for (int i = 0; i < match.size(); i++) {
            Dom[i].insert(match[i]);
        }
        return true;
    }

    int size;
    const int* cand;

    //for each u in order
    set<range> ranges;
    if (p.vorder[u] != -1)
        ranges.insert(g->GetVertices(p.vorder[u]));
    for (auto& e : p.adj[u]) {
        //already visited
        if (e.second < u) {
            //reverse direction
            range r = g->GetAdj(match[e.second], e.first, false); 
            ranges.insert(r);
        }
        else {
            range r = g->GetVertices(e.first, true);
            ranges.insert(r);
        }
    }
    for (auto& e : p.in_adj[u]) {
        //already visited
        if (e.second < u) {
            //reverse direction
            range r = g->GetAdj(match[e.second], e.first, true); 
            ranges.insert(r);
        }
        else {
            range r = g->GetVertices(e.first, false);
            ranges.insert(r);
        }
    }
    if (ranges.size() > 1) {
        int min_size = ranges.begin()->end - ranges.begin()->begin;
        if (min_size > intersect_capacity_[u]) {
            intersect_capacity_[u] = std::max(min_size, 2 * intersect_capacity_[u]);
            int* b = intersect_[u];
            int* t = temp_[u];
            intersect_[u] = new int[intersect_capacity_[u]];
            temp_[u] = new int[intersect_capacity_[u]];
            delete [] b;
            delete [] t;
        }
        size = intersect(ranges, intersect_[u], temp_[u]);
        cand = intersect_[u];
    } else {
        size = ranges.begin()->end - ranges.begin()->begin;
        cand = ranges.begin()->begin;
    }
    for (int i = 0; i < size; i++) {
        match[u] = cand[i];
        if (naiveRecur(p, match, u+1, Dom, start))
            return true;
    }
    match[u] = -1;
    return false;
}

