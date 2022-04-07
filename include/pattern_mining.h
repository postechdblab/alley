#ifndef PATTERN_MINING_H_ 
#define PATTERN_MINING_H_ 

#include "query_graph.h"
#include "data_graph.h"
#include "pattern.h"
#include "estimator.h"
#include "alley.h"

#define MAX_NUM_GROUP 32
#define INIT_INTER_CAPACITY 1048576
#define INIT_NE_CAPACITY 1048576
#define INIT_MAT_CAPACITY 1048576

#ifdef PRINT_INTERSECTION
extern bool calculate_domains_context;

extern
size_t totalEdgeCase1;
extern
size_t totalEdgeCase2;
extern
size_t totalLenRare;
extern
size_t totalLenFreq;
extern
size_t totalGalloping;
extern
size_t totalA;
extern
size_t totalB;
extern
size_t totalShuffling;
extern
size_t totalCountGalloping;
extern
size_t totalCountShuffling;

extern
size_t CD_totalEdgeCase1;
extern
size_t CD_totalEdgeCase2;
extern
size_t CD_totalLenRare;
extern
size_t CD_totalLenFreq;
extern
size_t CD_totalGalloping;
extern
size_t CD_totalA;
extern
size_t CD_totalB;
extern
size_t CD_totalShuffling;
extern
size_t CD_totalCountGalloping;
extern
size_t CD_totalCountShuffling;
#endif

#ifdef PRINT_SEARCH
extern
size_t CD_totalSearch;
extern
size_t CD_totalSearchLen;

extern
size_t totalSearch;
extern
size_t totalSearchLen;
#endif

class Estimator;

class PatternMining;
bool makeNEC(Pattern&, VID_MAP&);

struct MatInfo {
	PatternMining* pm;
	int L;
	int u;
	size_t index;
	size_t size;

	MatInfo(PatternMining* _pm = NULL, int _L = -1, size_t _i = MAX_SIZE, int _u = -1, size_t _s = MAX_SIZE) : pm(_pm), L(_L), index(_i), u(_u), size(_s) {}

	bool operator< (const MatInfo& other) const {
		if (this->size != other.size) {
			return this->size < other.size;
		}
		else {
			if (this->L != other.L)
				return this->L > other.L;
			else {
				if (this->index != other.index)
					return this->index > other.index;
				return this->u < other.u;
			}
		}
	}
};

struct DAGNode {
	Pattern pattern; //canonicalized
	string code;
	// {u_p -> {L, i, u_q}}, ordered by size
	vector<set<MatInfo>> mat_info; 
    
    DAGNode(const Pattern& p) : pattern(p) {}
};

enum miner_type { PATH_MINER, TREE_MINER, GRAPH_MINER, ALL_MINER }; 
enum label_type { VERTEX_LABEL, EDGE_LABEL };

class PatternMining {
private:
	double calculate_domains_time_;

	miner_type type_;
	label_type ltype_;
	bool do_mining_;
	int maxL_;
	double sampling_ratio_;
	double failure_threshold_;

	DataGraph* g;

	vector<int> label2order_;
	vector<vector<int>> order2label_;
	
	vector<vector<DAGNode>> DAG_;
	vector<unordered_map<string, size_t>> code2index_;

	vector<LargeVector<uint64_t>> consideredEdge_; //edge label group size <= 32

	vector<LargeVector<bool>> mat_old_;
	vector<LargeVector<char>> mat_i_;
	vector<LargeVector<size_t>> mat_o_;
	vector<LargeVector<unsigned int>> mat_s_;
	vector<LargeVector<int>> materialized_; //candidates are located from materialized_[m] + o to materialized_[m] + o + s

	inline void appendMat(int L, size_t index, int u, const int* arr, size_t size) {
		setMatOld(L, index, u, false);
		setMatIndex(L, index, u, L);
		setMatOffset(L, index, u, materialized_[L].size());
		setMatSize(L, index, u, size);
		size_t o = materialized_[L].size();
		materialized_[L].resize(o + size);
		memcpy(materialized_[L].data() + o, arr, sizeof(int) * size);
	}

	inline void pointMat(int L, size_t index, int u, const int* arr, size_t size) {
		setMatOld(L, index, u, false);
		setMatIndex(L, index, u, 0);
		setMatOffset(L, index, u, (size_t)arr);
		setMatSize(L, index, u, size);
	}

	inline void setConsideredEdge(int L, size_t index, int u, star& Elist) {
		for(auto& e: Elist)
			setConsideredEdge(L, index, u, e.first, e.second);
	}

	inline void setConsideredEdge(int L, size_t index, int u, int o, bool dir) {
		//bitmap operator
		assert(o >= 0 && o < MAX_NUM_GROUP);
		consideredEdge_[L][(L+1)*index + u] |= (1 << (o + (dir ? 0 : MAX_NUM_GROUP)));
	}
	inline bool getConsideredEdge(int L, size_t index, int u, int o, bool dir) {
		//bitmap operator
		return consideredEdge_[L][(L+1)*index + u] & (1 << (o + (dir ? 0 : MAX_NUM_GROUP)));
	}

	inline bool getMatOld(int L, size_t index, int u) {
		return mat_old_[L][(L+1)*index + u];
	}

	inline void setMatOld(int L, size_t index, int u, bool old) {
		mat_old_[L][(L+1)*index + u] = old;
	}

	inline char getMatIndex(int L, size_t index, int u) {
		return mat_i_[L][(L+1)*index + u];
	}

	inline void setMatIndex(int L, size_t index, int u, char i) {
		mat_i_[L][(L+1)*index + u] = i;
	}

	inline size_t getMatOffset(int L, size_t index, int u) {
		return mat_o_[L][(L+1)*index + u];
	}

	inline void setMatOffset(int L, size_t index, int u, size_t o) {
		mat_o_[L][(L+1)*index + u] = o;
	}

	inline unsigned int getMatSize(int L, size_t index, int u) {
		return mat_s_[L][(L+1)*index + u];
	}

	inline void setMatSize(int L, size_t index, int u, unsigned int s) {
		mat_s_[L][(L+1)*index + u] = s;
	}

	inline void clearMatInfo(int L, size_t index) {
		memset(mat_old_[L].data() + (L+1)*index, 0, sizeof(bool) * (L+1));
		memset(mat_o_[L].data() + (L+1)*index, 0, sizeof(size_t) * (L+1));
		memset(mat_s_[L].data() + (L+1)*index, 0, sizeof(unsigned int) * (L+1));
		memset(consideredEdge_[L].data() + (L+1)*index, 0, sizeof(uint64_t) * (L+1));
		for (int u = 0; u < L+1; u++)
			setMatIndex(L, index, u, (char)L);
	}
	
	vector<vector<bool>> marked_;
	vector<vector<bool>> pruned_;
	//bitmaps
	//vector<int> marked_;
	//vector<int> pruned_;
	//vector<map<string, >> intersect_cache_;
	//vector<vector<bool>> passed_;

	vector<vector<int>> Dom_;
	vector<vector<int>> Pr_;

	int** intersect_;
	int** temp_;
	int*  intersect_capacity_;

	int** initial_intersect_;
	int** initial_temp_;
	int*  initial_intersect_capacity_;

	map<vector<AdjElem>, vector<int>> intersection_cache_;

	vector<int> plan_;
	vector<int> label_;
	vector<vector<int>> subplans_;
	vector<int> rank_;
    vector<RBI> rbi_;
    vector<vector<AdjElem>> backward_;
	vector<star> forward_;
	vector<int> match_;
	vector<vector<int>> CT_;
	vector<bool> ap_;
	vector<vector<int>> to_clear_;
	vector<const int*>  search_space_;
	vector<int> space_size_;

	vector<vector<int>> black_child_;
	//for black 
	vector<bool> is_black_marked_;
	vector<set<int>> is_single_star_marked_;
	//for ivory 
	vector<set<match_star>> is_star_marked_;
	
	void printPlan() {
		cout << "plan: ";
		for (auto u : plan_)
			cout << u << " ";
		cout << endl;
	}

	void initPattern(int order);
	void addPatternToDAG(Pattern&, int L);
	void extendPattern(int L, size_t index, PatternMining* other = NULL);
	void addPatternToDAG(DAGNode& parent, Pattern&, int L);
	void reserveMemory(int L);
	void calculateDomains(int L);
	bool calculateDomainsWedge(size_t index);

	bool hasEmptyParent(int L, size_t index);
	void getMatInfo(DAGNode&, vector<vector<int>>& group, vector<MatInfo>&);
	void getMatInfo(DAGNode&, vector<MatInfo>&);
	void setDefaultDomains(int L, size_t index, vector<MatInfo>& infos, vector<star>& rels);

public:
	Estimator* estimator_; 
	size_t sampling_cnt_ = 0;
	size_t calculate_domains_cnt_ = 0;
	size_t review_cnt = 0;
	size_t review_sum = 0;


	//function pointers -- for Vlabel or Elabel
	string (*encode)(Pattern&, VID_MAP&);
	void (PatternMining::*myGroupLabels)(int);
	
	PatternMining(DataGraph* _g = NULL, Estimator* _e = NULL, miner_type _t = ALL_MINER, label_type _lt = EDGE_LABEL) : type_(_t), g(_g), estimator_(_e), ltype_(_lt) {
	    //Initialize function pointers
		switch(ltype_) {
			case EDGE_LABEL:
				encode = &encode_el;
				myGroupLabels = &PatternMining::GroupELabels;
				break;
			case VERTEX_LABEL:
				encode = &encode_vl;
				myGroupLabels = &PatternMining::GroupVLabels;
				break;
			default:
				assert(false);
		}
	}

	void Initialize(int maxL, double sampling_ratio, double failure_threshold, bool do_mining);
	void Clear();
	void GroupLabels(int num_bin) {
		(this->*myGroupLabels)(num_bin);
	}
	void GroupELabels(int num_bin);
	void GroupVLabels(int num_bin);
	void CopyIndex(PatternMining& other);
	void CopyKeys(PatternMining& other);
	void CopyOrder(PatternMining& other);
	void UnionIndex();
	void BuildIndex(PatternMining* other = NULL);
	bool unionTest();
	
	char GetLTypeInChar() {
		switch(ltype_) {
			case EDGE_LABEL:
				return 'e';
			case VERTEX_LABEL:
				return 'v';
			default:
				assert(false);
		}
	}
	
	char GetTypeInChar() {
		switch(type_) {
			case PATH_MINER:
				return 'p';
			case TREE_MINER:
				return 't';
			case GRAPH_MINER:
				return 'g';
			case ALL_MINER:
				return 'a';
			default:
				assert(false);
		}
	}

	void Write(const char* fn) {
		FILE* tp = fopen(fn, "w");
		FILE* mp = fopen((string(fn) + ".mat").c_str(), "w");
		
		fprintf(tp, "%zu %zu\n", label2order_.size(), order2label_.size());
		for (size_t i = 0; i < order2label_.size(); i++) {
			fprintf(tp, "%zu ", order2label_[i].size());
			for (size_t j = 0; j < order2label_[i].size(); j++)
				fprintf(tp, "%d ", order2label_[i][j]);
		}
		fprintf(tp, "\n");

		fprintf(tp, "%d\n", maxL_);
		for (int L = 2; L <= maxL_; L++) {
			fprintf(tp, "%zu\n", DAG_[L].size()); 
			cout << "# nodes for L " << L << " = " << DAG_[L].size() << endl;
			for (size_t i = 0; i < DAG_[L].size(); i++) { 
				fprintf(tp, "%s %zu\n", DAG_[L][i].code.c_str(), i); 
			}
		}
		
		for (int L = 2; L <= maxL_; L++) {
			mat_old_[L].resize((L+1) * DAG_[L].size());
			mat_i_[L].resize((L+1) * DAG_[L].size());
			mat_o_[L].resize((L+1) * DAG_[L].size());
			mat_s_[L].resize((L+1) * DAG_[L].size());
			consideredEdge_[L].resize((L+1) * DAG_[L].size());
			assert(mat_old_[L].size() == mat_i_[L].size() && mat_i_[L].size() == mat_o_[L].size() && mat_o_[L].size() == mat_s_[L].size() && consideredEdge_[L].size() == mat_i_[L].size());
			cout << "mat ios size for L " << L << " = " << mat_i_[L].size() << ", mat size = " << materialized_[L].size() << endl;
			fprintf(tp, "%zu %zu\n", mat_i_[L].size(), materialized_[L].size());
			fwrite(mat_old_[L].data(), sizeof(bool), mat_old_[L].size(), mp);
			fwrite(mat_i_[L].data(), sizeof(char), mat_i_[L].size(), mp);
			fwrite(mat_o_[L].data(), sizeof(size_t), mat_o_[L].size(), mp);
			fwrite(mat_s_[L].data(), sizeof(unsigned int), mat_s_[L].size(), mp);
			fwrite(consideredEdge_[L].data(), sizeof(uint64_t), consideredEdge_[L].size(), mp);
			fwrite(materialized_[L].data(), sizeof(int), materialized_[L].size(), mp);
		}
		fclose(tp);
		fclose(mp);
	}

	void ReadTempFile(const char* fn) {
		//for quick debugging
	}
	void WriteTempFile(const char* fn) {
		//for quick debugging
	}

	void Read(const char* fn) {
		FILE* tp = fopen(fn, "r");
		FILE* mp = fopen((string(fn) + ".mat").c_str(), "r");
		
		size_t num_label, num_order;
		fscanf(tp, "%zu %zu", &num_label, &num_order); 
		label2order_.resize(num_label, -1);
		order2label_.resize(num_order);
		for (size_t i = 0; i < order2label_.size(); i++) {
			size_t order_size;
			fscanf(tp, "%zu", &order_size);
			order2label_[i].resize(order_size);
			for (size_t j = 0; j < order_size; j++)
				fscanf(tp, "%d", &order2label_[i][j]);
		}
		for (size_t i = 0; i < order2label_.size(); i++) {
			for (size_t j = 0; j < order2label_[i].size(); j++) {
				int l = order2label_[i][j];
				label2order_[l] = i;
			}
		}

		fscanf(tp, "%d", &maxL_); 
		code2index_.resize(maxL_+1);
		mat_old_.resize(maxL_+1);
		mat_i_.resize(maxL_+1);
		mat_o_.resize(maxL_+1);
		mat_s_.resize(maxL_+1);
		consideredEdge_.resize(maxL_+1);
		materialized_.resize(maxL_+1);
		for (int L = 2; L <= maxL_; L++) {
			size_t num_code;
			char key[1024];
			size_t value;
			fscanf(tp, "%zu", &num_code);
			for (size_t i = 0; i < num_code; i++) {
				fscanf(tp, "%s %zu", key, &value); 
				code2index_[L][key] = value;
			}
		}
		for (int L = 2; L <= maxL_; L++) {
			size_t ios_size, mat_size;
			fscanf(tp, "%zu %zu", &ios_size, &mat_size);
			mat_old_[L].resize(ios_size);
			mat_i_[L].resize(ios_size);
			mat_o_[L].resize(ios_size);
			mat_s_[L].resize(ios_size);
			consideredEdge_[L].resize(ios_size);
			fread(mat_old_[L].data(), sizeof(bool), ios_size, mp);
			fread(mat_i_[L].data(), sizeof(char), ios_size, mp);
			fread(mat_o_[L].data(), sizeof(size_t), ios_size, mp);
			fread(mat_s_[L].data(), sizeof(unsigned int), ios_size, mp);
			fread(consideredEdge_[L].data(), sizeof(uint64_t), ios_size, mp);
			materialized_[L].resize(mat_size);
			fread(materialized_[L].data(), sizeof(int), mat_size, mp);
			cout << "mat ios size for L " << L << " = " << ios_size << ", mat size = " << mat_size << endl;
		}
		fclose(tp);
		fclose(mp);
	}

	void SetDataGraph(DataGraph* _g) {
		g = _g;
	}

	int getMaxL() { return maxL_; }
	
	vector<int>& getLabel2Order() {
		return label2order_;
	}

	int getOrder2LabelSize(int order) { 
		return order == PATTERN_EMPTY_VLABEL && ltype_ == VERTEX_LABEL ? 1 : order2label_[order].size(); 
	}
	vector<vector<int>>& getOrder2Label() { return order2label_; }
	int getOrder2Label(int order) { 
		return order == PATTERN_EMPTY_VLABEL && ltype_ == VERTEX_LABEL ? PATTERN_EMPTY_VLABEL : order2label_[order][0]; 
	}

	bool code2Index(int L, string code, size_t& id) {
		assert (code2index_.size() > L && L >= 0);
		if(code2index_[L].find(code) == code2index_[L].end())
			return false;
		id=code2index_[L][code];
		return true;
	}

	//if index == 1, domain is not materialized, just refernce the data graph
	//also, offset = edge label and size = direction
	inline range getMat(int L, size_t index, int u) {
		assert(getMatOld(L, index, u) == false);
		assert(L > 1);
		char i = getMatIndex(L, index, u);
		//added for reducing memory
		if (i == 0) {
			range r;
			auto size = getMatSize(L, index, u);
			if (size == 0) {
				r.begin = r.end = NULL;
			}
			else {
				r.begin = (int *)getMatOffset(L, index, u);
				r.end = r.begin + size; 
			}
			return r;
		}
		if (i == 1) {
			//could be both vertex or edge label
			if (getMatSize(L, index, u) == 2) {
				int vl = getMatOffset(L, index, u);
				return g->GetVertices(vl);
			}
			else {
				int el = getMatOffset(L, index, u);
				bool dir = getMatSize(L, index, u) > 0;
				return g->GetVertices(el, dir);
			}
		}
		range r;
		size_t s = getMatSize(L, index, u);
		if (s == 0) {
			r.begin = r.end = NULL;
			return r;
		}
		size_t o = getMatOffset(L, index, u);
		r.begin = materialized_[i].data() + o;
		r.end   = materialized_[i].data() + o + s;
		return r;
	}
	inline void GetConsideredElist(int L, size_t index, int u, star& Elist) {
		//p need not to be reordered; we only use label
		star new_Elist;
		for(auto& e: Elist) {
			bool isConsideredEdge = getConsideredEdge(L, index, u, e.first, e.second);
			if(isConsideredEdge) {
				new_Elist.insert(e);
			}
		}
		Elist=new_Elist;
	}
	inline range GetMatQuery(int L, size_t index, int u, bool& use_domain) {
		assert(L > 1);
		range r;
		if (getMatOld(L, index, u)) {
			use_domain = false;
			r.begin = r.end = NULL;
			return r;
		}
		char i = getMatIndex(L, index, u);
		assert(i != 0);
		if (i == 1) {
			use_domain = false;
			r.begin = r.end = NULL;
			return r;
		}
		size_t s = getMatSize(L, index, u);
		if (s == 0) {
			r.begin = r.end = NULL;
			return r;
		}
		size_t o = getMatOffset(L, index, u);
		r.begin = materialized_[i].data() + o;
		r.end   = materialized_[i].data() + o + s;
		return r;
	}
	
	bool isAcceptablePattern(Pattern& p) {
		switch(type_) {
			case PATH_MINER:
				return p.isPath();
			case TREE_MINER:
				return !p.isPath() && p.isTree();
			case GRAPH_MINER:
				return !p.isTree();
			case ALL_MINER:
				return true;
			default:
				assert(0);
		}
	}

	int getDomains(Pattern&, vector<MatInfo>& infos, vector<range>& domains, vector<bool>& use_domain, vector<star>& rels);
};

class MinerWrapper {

	DataGraph* g;
	Estimator* estimator;

public:
	vector<PatternMining> miners_;

	MinerWrapper(DataGraph* _g = NULL, Estimator* _e = NULL) : g(_g), estimator(_e) {
	} 

	void BuildIndex(int maxL, double sampling_ratio, double failure_threshold, int num_bin, label_type lt=EDGE_LABEL); 

	void Write(const char* fn); 

	void Read(const char* fn); 

	vector<PatternMining*> GetMiners() {
		vector<PatternMining*> m;
		for (PatternMining& miner : miners_)
			m.push_back(&miner);
		return m;
	}
};
	
//for query mode in estimators
bool getDomainsRecursive(const vector<PatternMining*>& pms, Pattern& q, const vector<int>& order2vid, int minMaxL, int maxMaxL, vector<range>& domains, vector<bool>& hasDomains, vector<star>& consideredElist, const int numVerticesWithDomains);

bool getDomains(PatternMining& pm, Pattern& q, vector<int>& order2vid, int MAX_INDEX_LEVEL, vector<range>& domains, vector<bool>& hasDomains, vector<star>& consideredElist);
bool getDomains(vector<PatternMining*>& pms, Pattern& q, vector<int>& order2vid, int minMaxL, int maxMaxL, vector<range>& domains, vector<bool>& hasDomains, vector<star>& consideredElist);

bool getDomainsDynamic(PatternMining& pm, QueryGraph& q, vector<int>& order2vid, int MAX_INDEX_LEVEL, int failure_order, AdjElem& failed_edge, vector<range>& domains, vector<bool>& hasDomains, vector<star>& consideredElist);

#endif
