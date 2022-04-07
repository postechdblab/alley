#include <cassert>
#include <random>
#include <unordered_set>
#include "../include/index_sample.h"

void IndexSample::Init() {
    generator_.seed(rand());
    ibjs_started_ = false;

    WanderJoin::Init();
}

void IndexSample::appendTuple(int node, int i, Tuple& tuple) {
    if (node >= offset_) {
        auto e = q->GetEdge(node - offset_);
        auto edge = g->GetEdge(e.el, i);
        tuple.push_back(edge);
    }
    else {
        int u = node_to_v_[node];
        int vl = q->GetVLabel(u);
        assert(vl != -1);
        tuple.push_back(vector<int>({g->GetVertex(vl, i)}));
    }
}

void IndexSample::appendTuple(int node, int i, int c, int vid, Tuple& tuple) {
    if (node >= offset_) {
        auto e = q->GetEdge(node - offset_);
        int val = g->GetAdj(vid, e.el, c == 0).begin[i];
        tuple.push_back(c == 0 ? vector<int>({vid, val}) : vector<int>({val, vid}));
    }
    else {
        int u = node_to_v_[node];
        int vl = q->GetVLabel(u);
        assert(vl != -1);
        assert(g->HasVLabel(vid, vl));
        tuple.push_back(vector<int>({vid}));
    }
}

int IndexSample::getBoundedSpaceSize(int node, int c, int v) {
    if (node >= offset_) {
		auto e = q->GetEdge(node - offset_);
        return g->GetAdjSize(v, e.el, c == 0);
    }
    else {
        int vl = q->GetVLabel(node_to_v_[node]);
        if (g->HasVLabel(v, vl))
            return 1;
        else
            return 0;
    }
}

int IndexSample::getElementFromBoundedSpace(int node, int c, int v, int i) {
    if (node >= offset_) {
        auto e = q->GetEdge(node - offset_);
        return g->GetAdj(v, e.el, c == 0).begin[i];
    }
    else {
        assert(i == 0);
        return v;
    }
}

vector<Tuple> IndexSample::Iter(int p, vector<Tuple>& tuples) {
	vector<Tuple> ret;

	const int qt = plan_[p].first; 
	const int c = plan_[p].second;
	int pos = 0;
	for (int k = 0; k < p; k++)
		if (plan_[k].first == counterpart_[p].first)
			pos = k + 1;
	long long sum = 0;
	vector<long long> offs;

    int temp_size = 0;
	int cattr = 0, vidattr = 0;
	bool use_index = getBoundInfo(qt, cattr, vidattr, temp_size);
		
	if (use_index) {
        if (c == cattr) {
            for (int i = 0; i < tuples.size(); i++) {
                int v = tuples[i][pos][counterpart_[p].second];
                int cnt = 0;
                if (v == vidattr) {
                    cnt = getBoundedSpaceSize(qt, c, v); 
                }
                sum += cnt;
                offs.push_back(sum);
            }
        }
		else  {
			for (int i = 0; i < tuples.size(); i++) {
				int v = tuples[i][pos][counterpart_[p].second];
				int s = 0;
				int e = getBoundedSpaceSize(qt, c, v); 
				if (e > BINARY_THRESHOLD) {
					int mid;
					while (s < e) {
						mid = (s + e) / 2;
                        int val = getElementFromBoundedSpace(qt, c, v, mid);
						if (val < vidattr)
							s = mid + 1;
						else if (val > vidattr)
							e = mid;
						else {
							ret.push_back(tuples[i]);
							if (cattr == 1) ret.back().push_back(vector<int>({ v, val }));
							else ret.back().push_back(vector<int>({ val, v }));
							sum++;
                            offs.push_back(sum);
							break;
						}
					}
				} else {
					int node = qt;
					auto edge = q->GetEdge(node - offset_);
					range r = g->GetAdj(v, edge.el, c == 0);
					for (; r.begin != r.end; r.begin++) {
                        if (*r.begin == vidattr) {
                            ret.push_back(tuples[i]);
						    ret.back().push_back(cattr == 0 ? vector<int>({*r.begin, v}) : vector<int>({v, *r.begin}));
						    sum++;
                            offs.push_back(sum);
						    break;
                        }
					}          
                }
			}
			if (p == plan_.size() - 1) {
				inv_prob_ *= sum; 
            }

			return ret;
		}
	}

    if (offs.size() == 0) {
        for (int i = 0; i < tuples.size(); i++) {
            int v = tuples[i][pos][counterpart_[p].second];
            int cnt = getBoundedSpaceSize(qt, c, v);
            sum += cnt;
            offs.push_back(sum);
        }
    }

	inv_prob_ *= sum;

    assert(ret.empty());
    if (sum == 0) {
		return ret;
	}

	if (sum <= sample_cnt_) {
		for (int i = 0; i < tuples.size(); i++) {
			int v = tuples[i][pos][counterpart_[p].second];
			if (use_index && v != vidattr) {
                assert(c == cattr);
				continue;
            }
			int node = qt;
			if (node >= offset_) {
				auto e = q->GetEdge(node - offset_);
				range r = g->GetAdj(v, e.el, c == 0);
				for (; r.begin != r.end; r.begin++) {
					ret.push_back(tuples[i]);
					ret.back().push_back(c == 1 ? vector<int>({*r.begin, v}) : vector<int>({v, *r.begin}));
				}
			} else {
				int u = node_to_v_[node];
				int vl = q->GetVLabel(u);
                assert(vl != -1);
                if (g->HasVLabel(v, vl)) {
                    ret.push_back(tuples[i]);
                    ret.back().push_back(vector<int>({v}));
                }
            }
        }
        if (!use_index)
            assert(ret.size() == sum);
        assert(ret.size() <= sum);
	} else {
		std::unordered_set<long long> index;
		for (long long r = sum - sample_cnt_; r < sum; r++) {
			long long v = std::uniform_int_distribution<long long>(1, r)(generator_);
			if (!index.insert(v).second)
				index.insert(r);
		}
		assert(index.size() == sample_cnt_);
		std::set<long long> sorted(index.begin(), index.end());
		long long chosen = 0;
		for (auto j : sorted) {
			while (offs[chosen] <= j)
				chosen++;
			assert(chosen < tuples.size());
			int offset = chosen == 0 ? j : j - offs[chosen-1];
			assert(offset >= 0);
			int v = tuples[chosen][pos][counterpart_[p].second];
            if (use_index && v != vidattr) {
                assert(c == cattr);
				continue;
            }
			int node = qt;
			int val;
			if (node >= offset_) {
				auto e = q->GetEdge(node - offset_);
				val = g->GetAdj(v, e.el, c == 0).begin[offset];
				ret.push_back(tuples[chosen]);
				ret.back().push_back(c == 0 ? vector<int>({v, val}) : vector<int>({val, v}));
			} else {
				int u = node_to_v_[node];
				int vl = q->GetVLabel(u);
                assert(vl != -1);
                ret.push_back(tuples[chosen]);
                ret.back().push_back(vector<int>({v}));
			}
		}
        if (!use_index)
            assert(ret.size() == sample_cnt_);
        assert(ret.size() <= sample_cnt_);
	}

    if (p != plan_.size() - 1) {
        inv_prob_ /= ret.size();
    }

	return ret;
}

bool IndexSample::getBoundInfo(int node, int& c, int& v, int& size) {
    if (node >= offset_) {
		auto e = q->GetEdge(node - offset_);
		if (q->GetBound(e.src) >= 0) {
            c = 0;
			v = q->GetBound(e.src);
            assert(v != -1);
            size = g->GetAdjSize(v, e.el, true);
            return true;
		}
        if (q->GetBound(e.src) < 0 && q->GetBound(e.dst) >= 0) {
            c = 1;
			v = q->GetBound(e.dst);
            assert(v != -1);
            size = g->GetAdjSize(v, e.el, false);
            return true;
		}
	} else {
        assert(node < node_to_v_.size());
		int u = node_to_v_[node];
		if (q->GetBound(u) >= 0) {
            c = 0;
			v = q->GetBound(u);
            assert(v != -1);
			size = 1;
            if (q->GetVLabel(u) != -1)
                if (!g->HasVLabel(v, q->GetVLabel(u))) 
                    size = 0;
            return true;
		}
	}
    return false;
}

bool IndexSample::GetSubstructure(int subquery_index) {
    if (!plan_chosen_) {
        return WanderJoin::GetSubstructure(subquery_index);
    }
    else {
        if (ibjs_started_)
            return false;
    }
    assert(plan_chosen_);
    ibjs_started_ = true;
    valid_ = true;

    assert(plans_generated_);
    valid_ = true;

    counterpart_ = counterparts_[pos_];
    plan_ = plans_[pos_];

	int start_table = 0;
	if (counterpart_.size() > 0)
		start_table = counterpart_[0].first;
	int size = -1;
	if (start_table >= offset_) {
		auto e = q->GetEdge(start_table - offset_);
		size = g->GetNumEdges(e.el);
	} else {
        assert(start_table < node_to_v_.size());
		int u = node_to_v_[start_table];
		int vl = q->GetVLabel(u);
        assert(vl != -1);
		size = g->GetNumVertices(vl);
	}
	int c = -1, vidattr = -1;
	bool use_index = getBoundInfo(start_table, c, vidattr, size);
		
	vector<Tuple> tuples;
	if (size <= sample_cnt_) {
		tuples.resize(size);
		for (int i = 0; i < size; i++) {
			if (use_index) {
                appendTuple(start_table, i, c, vidattr, tuples[i]);
			}
			else {
                appendTuple(start_table, i, tuples[i]);
			}
		}
	} else {
		tuples.resize(sample_cnt_);
		std::unordered_set<int> start;
		for (int r = size - sample_cnt_; r < size; r++) {
			int v = std::uniform_int_distribution<>(1, r)(generator_);
			if (!start.insert(v).second)
				start.insert(r);
		}
		assert(start.size() == sample_cnt_);
		int i = 0;
		for (auto s : start) {
			if (use_index) {
                appendTuple(start_table, s, c, vidattr, tuples[i]); 
			}
			else {
                appendTuple(start_table, s, tuples[i]);
			}
			i++;
		}
	}

	inv_prob_ = tuples.size() == 0 ? 0 : (double)size/std::min(size, sample_cnt_);

	if (tuples.size() == 0) {
		return true;
	}

	for (int p = 0; p < plan_.size(); p++) {
		tuples = Iter(p, tuples);
		if (tuples.size() == 0) {
			return true;
		}
	}

    //check join conditions
    if (tuples.size() > 0 && tuples[0].size() == walk_size_) {
        int num_valid = 0;
        for (int i = 0; i < tuples.size(); i++) {
            if (checkNonTreeEdges(tuples[i]))
                num_valid++;
        }
        inv_prob_ *= (double)num_valid / tuples.size();
    }

	return true;
}

double IndexSample::AggCard() {
    double res = 0.0;
    if (ibjs_started_) {
        double ibjs_est = card_vec_.back(); 
        card_vec_.pop_back();
        double ibjs_weight = (double)sample_cnt_/sample_size_;
        res += ibjs_est * ibjs_weight;
    }
    double wj_est = WanderJoin::AggCard();
    double wj_weight = (double)(sample_size_ - sample_cnt_)/sample_size_;
    res += wj_est * wj_weight;
    return res;
}

bool IndexSample::checkNonTreeEdges(Tuple& tuple) {
    assert(tuple.size() == walk_size_);
    for (int i = 0; i < join_from_.size(); i++) {
		int n1 = join_from_[i].first;
		int c1 = join_from_[i].second;
		int n2 = join_to_[i].first;
		int c2 = join_to_[i].second;
		int pos1 = 0, pos2 = 0;
		for (int j = 0; j < plan_.size(); j++) {
			if (plan_[j].first == n1)
				pos1 = j + 1;
			if (plan_[j].first == n2)
				pos2 = j + 1;
		}
		if (tuple[pos1][c1] != tuple[pos2][c2])
			return false;
	}
	return true;
}
