#include <cassert>
#include <random>
#include "../include/wander_join.h"

void WanderJoin::PrepareSummaryStructure(DataGraph& g, double ratio) {

}

void WanderJoin::WriteSummary(const char* fn) {

}

void WanderJoin::ReadSummary(const char* fn) {

}

void WanderJoin::Init() {
    plans_generated_ = plan_chosen_ = false;
	walk_size_ = sample_size_ = sample_cnt_ = pos_ = inv_prob_ = 0;

    visited_.clear();
    join_from_.clear();
    join_to_.clear();
    node_to_v_.clear();
    counterparts_.clear();
    plans_.clear();
    counterpart_.clear();
    plan_.clear();
    success_cnt_.clear();
    est_.clear();
    num_idx_lookup_.clear();
    sampled_tuples_.clear();
    
    //set sample size
    int sum = 0;
    offset_ = 0;
    int e_cnt = q->GetNumEdges(); 
    for (int v = 0; v < q->GetNumVertices(); v++) {
        int vl = q->GetVLabel(v); 
        if (vl != -1) {
            sum += g->GetNumVertices(vl); 
            node_to_v_[offset_] = v;
            offset_++;
        }
        for (auto& e : q->GetAdj(v, true)) {
            sum += g->GetNumEdges(e.second);
        }
    }
    sample_size_ = sum / (offset_ + e_cnt); 
    sample_size_ = ceil(sample_size_ * sample_ratio);
    walk_size_ = offset_ + e_cnt;
    
    //set join from and to 
    for (int i = 0; i < walk_size_; i++) {
        for (int j = i + 1; j < walk_size_; j++) {
            //edges i and j can have join conditions on both query vertices, 
            //i.e., i and j makes a loop
            int c1[2] = {-1, -1}; 
            int c2[2] = {-1, -1};
            int n = 0;
            if (i < offset_) {
                int v = node_to_v_[i]; 
                int vl = q->GetVLabel(v);
                assert(vl != -1);
                if (j >= offset_) {
                    auto e2 = q->GetEdge(j - offset_);
                    if (v == e2.src) {
                        c1[n] = c2[n] = 0;
                        n++;
                    }
                    if (v == e2.dst) {
                        c1[n] = 0;
                        c2[n] = 1;
                        n++;
                    }
                } else {
                    int u = node_to_v_[j]; 
                    assert(u != v);
                }
                //do not consider j >= e_cnt since no q vertex has more than 1 label
                //so c1 and c2 remains -1
            } else {
                auto e1 = q->GetEdge(i - offset_);
                if (j >= offset_) {
                    auto e2 = q->GetEdge(j - offset_);
                    if (e1.src == e2.src) {
                        c1[n] = c2[n] = 0;
                        n++;
                    }
                    if (e1.dst == e2.src) {
                        c1[n] = 1;
                        c2[n] = 0;
                        n++;
                    }
                    if (e1.src == e2.dst) {
                        c1[n] = 0;
                        c2[n] = 1;
                        n++;
                    }
                    if (e1.dst == e2.dst) {
                        c1[n] = c2[n] = 1;
                        n++;
                    }
                } else {
                    int v = node_to_v_[j]; 
                    int vl = q->GetVLabel(v);
                    assert(vl != -1);
                    if (e1.src == v) {
                        c1[n] = c2[n] = 0;
                        n++;
                    }
                    if (e1.dst == v) {
                        c1[n] = 1;
                        c2[n] = 0;
                        n++;
                    }
                }
            }
            for (int m = 0; m < n; m++) { 
                if (c1[m] != -1 && c2[m] != -1) {
                    join_from_.push_back(make_pair(i, c1[m]));
                    join_to_.push_back(make_pair(j, c2[m]));
                    join_from_.push_back(make_pair(j, c2[m]));
                    join_to_.push_back(make_pair(i, c1[m]));
                }
            }
        }
    }

}

int WanderJoin::DecomposeQuery() {
    return 1;
}

//top function
void WanderJoin::generateWalkPlans() {
    for (int s = 0; s < walk_size_; s++) {
        visited_.clear();
        visited_.resize(walk_size_, false);
        visited_[s] = true;
        generateWalkPlans(0);
    }
}

void WanderJoin::generateWalkPlans(int cnt) {
    if (sample_cnt_ == 0) 
        return;
    if (cnt == walk_size_ - 1) {
        sample_cnt_--;
        assert(plan_.size() == walk_size_ - 1);
        assert(counterpart_.size() == walk_size_ - 1);
        plans_.push_back(plan_);
        counterparts_.push_back(counterpart_);
        return;
    }
    for (int i = 0; i < join_from_.size(); i++) {
        if (visited_[join_from_[i].first] && !visited_[join_to_[i].first]) {
            visited_[join_to_[i].first] = true;
            counterpart_.push_back(join_from_[i]);
            plan_.push_back(join_to_[i]);
            generateWalkPlans(cnt + 1);
            visited_[join_to_[i].first] = false;
            counterpart_.pop_back();
            plan_.pop_back();
        }
    }
}

#ifdef LEAF_SCANNING
//check bounds, non-tree edges
int WanderJoin::getNumGoodLeafTuples(vector<pair<int, int>>& plan, int node, int v, int c) {
    int ret = 0;
    vector<tuple<bool, int, int>> leaf_conditions;
    assert(sampled_tuples_.size() == plan.size());

    for (int i = 0; i < join_from_.size(); i++) {
		int n1 = join_from_[i].first;
		int c1 = join_from_[i].second;
		int n2 = join_to_[i].first;
		int c2 = join_to_[i].second;
		int pos1 = 0, pos2 = 0;
		for (int j = 0; j < plan.size(); j++) {
			if (plan[j].first == n1)
				pos1 = j + 1;
			if (plan[j].first == n2)
				pos2 = j + 1;
		}
        //keep and use it later
        if (pos1 == sampled_tuples_.size()) {
            assert(pos2 < sampled_tuples_.size());
            leaf_conditions.emplace_back(c1 == 0, pos2, c2);
        } else if (pos2 == sampled_tuples_.size()) {
            assert(pos1 < sampled_tuples_.size());
            leaf_conditions.emplace_back(c2 == 0, pos1, c1);
        } else { 
            assert(pos1 < sampled_tuples_.size());
            assert(pos2 < sampled_tuples_.size());
            if (sampled_tuples_[pos1][c1] != sampled_tuples_[pos2][c2])
                return 0;
        }
	}

	if (node >= offset_) {
		auto e = q->GetEdge(node - offset_);

        //don't need to scan
        if (leaf_conditions.empty()) {
            if (q->GetBound(e.src) < 0 && q->GetBound(e.dst) < 0)
                return g->GetAdjSize(v, e.el, c == 0);
            else if (q->GetBound(e.src) >= 0 && q->GetBound(e.dst) < 0) {
                if (c == 0) {
                    if (v == q->GetBound(e.src))
                        return g->GetAdjSize(v, e.el, true);
                }
                else {
                    if (g->HasEdge(q->GetBound(e.src), v, e.el, true))
                        return 1;
                }
            }
            else if (q->GetBound(e.src) < 0 && q->GetBound(e.dst) >= 0) {
                if (c == 1) {
                    if (v == q->GetBound(e.dst))
                        return g->GetAdjSize(v, e.el, false);
                }
                else {
                    if (g->HasEdge(v, q->GetBound(e.dst), e.el, true))
                        return 1;
                }
            }
            else {
                if (c == 0) {
                    if (v == q->GetBound(e.src) && g->HasEdge(v, q->GetBound(e.dst), e.el, true))
                        return 1;
                }
                else {
                    if (v == q->GetBound(e.dst) && g->HasEdge(v, q->GetBound(e.src), e.el, false))
                        return 1;
                }
            }
            return 0;
        }

        range r = g->GetAdj(v, e.el, c == 0);
        for (; r.begin != r.end; ++r.begin) {
            int temp = *r.begin;
            int src = (c == 0) ? v : temp;
            int dst = (c == 0) ? temp : v;
            if (q->GetBound(e.src) >= 0 && q->GetBound(e.src) != src)
                continue;
            if (q->GetBound(e.dst) >= 0 && q->GetBound(e.dst) != dst)
                continue;
            bool valid = true;
            for (auto& t : leaf_conditions) {
                //src
                if (get<0>(t)) {
                    if (src != sampled_tuples_[get<1>(t)][get<2>(t)]) {
                        valid = false;
                        break;
                    }
                }
                //dst
                else {
                    if (dst != sampled_tuples_[get<1>(t)][get<2>(t)]) {
                        valid = false;
                        break;
                    }
                }
            }
            if (valid)
                ret++;
        }
    }
    else {
        //don't need to check bounds or non-trees
        int u = node_to_v_[node];
        int vl = q->GetVLabel(u);
        if (g->HasVLabel(v, vl))
            ret = 1;
        else 
            ret = 0;
    }
    return ret;
}
#endif

//sample from the first node in the walk plan
//returns inverse probability
int WanderJoin::sampleTuple(int node) {
    int ret;
	if (node >= offset_) {
		auto e = q->GetEdge(node - offset_);
        vector<int> t;
#ifdef SELECTIVE_PREDICATE
        if (q->GetBound(e.src) >= 0) {
            t = g->GetRandomEdge(q->GetBound(e.src), e.el, true);
            ret = g->GetAdjSize(q->GetBound(e.src), e.el, true);
        } else if (q->GetBound(e.dst) >= 0) {
            t = g->GetRandomEdge(q->GetBound(e.dst), e.el, false);
            ret = g->GetAdjSize(q->GetBound(e.dst), e.el, false);
        }
        //fallback
        else {
#endif
            t = g->GetRandomEdge(e.el);
            ret = g->GetNumEdges(e.el);
#ifdef SELECTIVE_PREDICATE
        }
#endif
        sampled_tuples_.push_back(t);
    } else {
        int u = node_to_v_[node];
        int vl = q->GetVLabel(u);
#ifdef SELECTIVE_PREDICATE
        if (q->GetBound(u) >= 0) {
            if (g->HasVLabel(q->GetBound(u), vl)) { 
                sampled_tuples_.push_back(vector<int>(1, q->GetBound(u)));
                return 1;
            }
            else
                return 0;
        }
        //fallback
        else {
#endif
            auto t = g->GetRandomVertex(vl); 
            sampled_tuples_.push_back(t);
            ret = g->GetNumVertices(vl);
#ifdef SELECTIVE_PREDICATE
        }
#endif
	}
    review_sum += ret;
    review_cnt ++;
	return ret;
}

//sample from subsequent nodes, using graph adj. list
//returns inverse probability
int WanderJoin::sampleTuple(int node, int v, int c) {
    int ret;
    if (node >= offset_) {
        auto e = q->GetEdge(node - offset_);
#ifdef SELECTIVE_PREDICATE
        int bound = (c == 0) ? q->GetBound(e.src) : q->GetBound(e.dst);
        int other = (c == 0) ? q->GetBound(e.dst) : q->GetBound(e.src);
        if (bound >= 0 && bound != v)
            return 0;
        if (other >= 0) {
            if (g->HasEdge(v, other, e.el, c == 0)) {
                if (c == 0)
                    sampled_tuples_.push_back(vector<int>({v, other}));
                else
                    sampled_tuples_.push_back(vector<int>({other, v}));
                return 1;
            }
            else
                return 0;
        }
        //fallback
        else {
#endif
            auto t = g->GetRandomEdge(v, e.el, c == 0);
            if (t.size() > 0) {
                sampled_tuples_.push_back(t);
                ret = g->GetAdjSize(v, e.el, c == 0);
            } else
                ret = 0;
#ifdef SELECTIVE_PREDICATE
        }
#endif
    } else {
        int u = node_to_v_[node];
        int vl = q->GetVLabel(u);
#ifdef SELECTIVE_PREDICATE
        if (q->GetBound(u) >= 0) {
            if (q->GetBound(u) == v && g->HasVLabel(v, vl)) {
                sampled_tuples_.push_back(vector<int>(1, v));
                return 1;
            }
            else 
                return 0;
        }
        else { 
#endif
            if (g->HasVLabel(v, vl)) {
                sampled_tuples_.push_back(vector<int>(1, v));
                return 1;
            } else 
                return 0;
#ifdef SELECTIVE_PREDICATE
        }
#endif
    }
    review_sum += ret;
    review_cnt ++;
    return ret;
}

bool WanderJoin::checkBoundedVertices(int node, vector<int> t) {
#ifdef SELECTIVE_PREDICATE
    return true;
#endif
	if (node >= offset_) {
		auto e = q->GetEdge(node - offset_);
		if (q->GetBound(e.src) >= 0 && q->GetBound(e.src) != t[0])
			return false;
		if (q->GetBound(e.dst) >= 0 && q->GetBound(e.dst) != t[1])
			return false;
	} else {
		int u = node_to_v_[node];
		if (q->GetBound(u) >= 0 && q->GetBound(u) != t[0])
			return false;
	}
	return true;
}

bool WanderJoin::checkLabelStatistics(int node) { 
	if (node >= offset_) {
		auto e = q->GetEdge(node - offset_);
		if (g->GetNumEdges(e.el) == 0) {
			return false;
        }
	} else {
		int u = node_to_v_[node];
		int vl = q->GetVLabel(u);
		if (g->GetNumVertices(vl) == 0) {
			return false;
        }
	}
	return true;
}

bool WanderJoin::checkNonTreeEdges(vector<pair<int, int>>& plan) {
	for (int i = 0; i < join_from_.size(); i++) {
		int n1 = join_from_[i].first;
		int c1 = join_from_[i].second;
		int n2 = join_to_[i].first;
		int c2 = join_to_[i].second;
		int pos1 = 0, pos2 = 0;
		for (int j = 0; j < plan.size(); j++) {
			if (plan[j].first == n1)
				pos1 = j + 1;
			if (plan[j].first == n2)
				pos2 = j + 1;
		}
		if (sampled_tuples_[pos1][c1] != sampled_tuples_[pos2][c2])
			return false;
	}
	return true;
}

#ifdef TRIAL_WALKS
void WanderJoin::adjustZ() {
    int begin = z_opt_index_ == 0 ? 0 : z_opt_index_ + 1;
    int end   = plans_.size(); 
    int mid;
    while (begin != end) {
        mid = (begin + end) / 2;
        double cur_eq = (2 * y_info_.var / trial_info_[mid].var - 1) * y_info_.cnt; 
        if (cur_eq < trial_info_[mid].cnt) 
            begin = mid + 1;
        else 
            end = mid;
    }
    z_opt_index_ = std::min(begin, (int)plans_.size() - 1);
}
#endif

//generates all walk orders
//perform random walks
//returns si and P(si)
bool WanderJoin::GetSubstructure(int subquery_index) {
	if (!plans_generated_) {
		sample_cnt_ = sample_size_;
		//this consumes sample_cnt_, preventing generating many useless walk orders
		//(an optimization)
		generateWalkPlans();
        //restore sample_cnt_
		sample_cnt_ = sample_size_;
		success_cnt_.resize(plans_.size(), 0);
		est_.resize(plans_.size());
		num_idx_lookup_.resize(plans_.size());
#ifdef TRIAL_WALKS
		trial_info_.clear();
#endif
		plans_generated_ = true;
		if (plans_.size() == 0)
			return false;
	}
	if (sample_cnt_ <= 0)
		return false;
	while (sample_cnt_) {
		sample_cnt_--;

		sampled_tuples_.clear();
		valid_ = true;
		inv_prob_ = 1.0;
        int lookup = 0;

		int start_node = 0;
		if (counterparts_[pos_].size() > 0)
			start_node  = counterparts_[pos_][0].first;
		assert(start_node < walk_size_);
		if (!checkLabelStatistics(start_node))
			return false;
		//randomly sample an edge/vertex with edge/vertex label of start_node
		inv_prob_ *= sampleTuple(start_node);
        lookup++;
		if (inv_prob_ == 0 || !checkBoundedVertices(start_node, sampled_tuples_[0])) {
			if (!plan_chosen_) {
                est_[pos_].push_back(0);
                num_idx_lookup_[pos_].push_back(lookup);
				pos_ = (pos_ + 1) % plans_.size();
            }
            valid_ = false;
            return true;
		}

		for (int cur_order = 0; cur_order < plans_[pos_].size(); cur_order++) {
			int cur_node = plans_[pos_][cur_order].first;
			int prev_order = 0;
			for (int k = 0; k < cur_order; k++) 
				if (plans_[pos_][k].first == counterparts_[pos_][cur_order].first) 
					prev_order = k + 1;
			int c = counterparts_[pos_][cur_order].second;
			int v = sampled_tuples_[prev_order][c];
#ifdef LEAF_SCANNING
            //last table, return total # tuples that match non-tree edges and bounds
            //don't need to checkNonTreeEdges later
            if (cur_order == plans_[pos_].size()-1) {
                inv_prob_ *= getNumGoodLeafTuples(plans_[pos_], cur_node, v, plans_[pos_][cur_order].second); 
                lookup++;
                if (inv_prob_ == 0)
                    valid_ = false;
                break;
            }
            else {
#endif
                inv_prob_ *= sampleTuple(cur_node, v, plans_[pos_][cur_order].second); 
#ifdef LEAF_SCANNING
            }
#endif
            lookup++;
			if (inv_prob_ == 0) {
				valid_ = false;
				break;
			}
			if (!checkBoundedVertices(cur_node, sampled_tuples_.back())) {
				valid_ = false;
				break;
			}
		}

#ifdef TRIAL_WALKS
        //don't calculate variance since it requires O(y)

        if (plan_chosen_) {
            y_info_.mean = (y_info_.mean * y_info_.cnt + inv_prob_) / (y_info_.cnt + 1); 
            y_info_.cnt++;
            double cur_eq = (2 * y_info_.var / trial_info_[z_opt_index_].var - 1) * y_info_.cnt; 
            if (y_info_.cnt == 1 || cur_eq > trial_info_[z_opt_index_].cnt)
                adjustZ();
        }
#endif

#ifdef LEAF_SCANNING
		if (!valid_)
#else
		if (!valid_ || !checkNonTreeEdges(plans_[pos_])) 
#endif
        {
			if (!plan_chosen_) {
                est_[pos_].push_back(0);
                num_idx_lookup_[pos_].push_back(lookup);
				pos_ = (pos_ + 1) % plans_.size();
            }
            valid_ = false;
            return true;
        }
        success_cnt_[pos_]++;
        if (!plan_chosen_) {
            est_[pos_].push_back(inv_prob_);
            num_idx_lookup_[pos_].push_back(lookup);
            if (success_cnt_[pos_] >= 100) {
                int min_pos = -1;
                double min_val = std::numeric_limits<double>::max();
                for (int i = 0; i < plans_.size(); i++) {
                    if (success_cnt_[i] >= 50) {
                        double sum = 0;
                        for (double e : est_[i])
                            sum += e;
                        double mean = sum / est_[i].size();
                        double var = 0;
                        for (double e : est_[i])
                            var += (e - mean) * (e - mean);
                        //variance of estimates 
                        var = var / (est_[i].size() - 1); 

#ifdef TRIAL_WALKS
                        trial_info_.emplace_back(est_[i].size(), var, mean);
#endif

                        sum = 0;
                        for (double l : num_idx_lookup_[i])
                            sum += l;
                        //mean of index lookups
                        mean = sum / num_idx_lookup_[i].size();
                        if (var * mean < min_val) {
                            min_pos = i;
                            min_val = var * mean;
                        }
                    }
                    else {
#ifdef TRIAL_WALKS
                        if (est_[i].empty()) {
                            trial_info_.emplace_back(est_[i].size(), std::numeric_limits<double>::max(), 0);
                        }
                        else {
                            double sum = 0;
                            for (double e : est_[i])
                                sum += e;
                            double mean = sum / est_[i].size();
                            double var = 0;
                            if (est_[i].size() > 1) {
                                for (double e : est_[i])
                                    var += (e - mean) * (e - mean);
                                //variance of estimates 
                                var = var / (est_[i].size() - 1); 
                            }

                            trial_info_.emplace_back(est_[i].size(), var, mean);
                        }
#endif
                    }
                }
                assert(min_pos != -1);
                pos_ = min_pos;
                plan_chosen_ = true;

#ifdef TRIAL_WALKS
                assert(trial_info_.size() == plans_.size());
                //initially set to 0
                z_opt_index_ = 0;
                y_info_.cnt = 0;
                y_info_.var = trial_info_[pos_].var;
                y_info_.mean = 0;

                //sort plans by var
                std::sort(trial_info_.begin(), trial_info_.end(), 
                        [](const TrialInfo& a, const TrialInfo& b) -> bool { return a.var < b.var; });
                //calculate cumulative values
                //leave out the values for chosen plan
                for (int i = 1; i < plans_.size(); i++) {
                    auto& prev = trial_info_[i-1];
                    auto& cur  = trial_info_[i];
                    int c = cur.cnt;
                    cur.cnt = prev.cnt + cur.cnt;
                    assert(cur.cnt != 0);
                    cur.var = (prev.var * prev.cnt + cur.var * c) / cur.cnt;
                    cur.mean = (prev.mean * prev.cnt + cur.var * c) / cur.cnt;
                }
#endif
            }
            else 
                pos_ = (pos_ + 1) % plans_.size();
        }
		return true;
	}
	//used all sample_cnt_ before choosing an order
	if (sample_cnt_ == 0)
		return false;
}

//HT estimator,
//returns 1/P(si) if valid, 0 otherwise
double WanderJoin::EstCard(int subquery_index) {
    if (valid_) {
        return inv_prob_;
    }
    else
        return 0.0;
}

double WanderJoin::AggCard() {
#ifdef TRIAL_WALKS
    if (plan_chosen_) {
        return (y_info_.mean * y_info_.cnt + 
                trial_info_[z_opt_index_].mean * trial_info_[z_opt_index_].cnt)
            / (y_info_.cnt + trial_info_[z_opt_index_].cnt); 
    }
#endif
    double res = 0.0;
    if (card_vec_.size() == 0)
        return res;
    for (double card : card_vec_) res += card;
    return res / card_vec_.size();
}

double WanderJoin::GetSelectivity() {
    return 1;
}
