#ifndef ALLEY_H_ 
#define ALLEY_H_ 

#include <random>
#include "../include/estimator.h"
#include "../include/pattern_mining.h"

#define MAX_ORDER 32
#define CARTESIAN_THRESHOLD 100000000 

class Alley : public Estimator {
public:
	Alley(bool);
	~Alley();
	void PrepareSummaryStructure(DataGraph&, double); 
	void WriteSummary(const char*); 
	void ReadSummary(const char*); 
	void Init();
	int  DecomposeQuery(); 
	bool GetSubstructure(int); 
	double EstCard(int); 
	double AggCard();
	double GetSelectivity();
	void PrintStat();
	
	void DumpEstimationDetails(std::fstream& fout);
	double GetFailureRate();
	
	bool calculateDomainsNaive(Pattern& p, vector<set<int>>& Dom);
	bool calculateDomainsFast(Pattern&, vector<range>& domains, vector<bool>& use_domain, vector<star>& rels, bool debug, bool& timeout);
	bool treeCycleInternal(Pattern&, vector<range>& domains, vector<bool>& use_domain, vector<star>& rels, const int level, bool& timeout, bool print);
	bool treeCycleRecur(const int head, const int index, int& rollback, vector<range>& domains, vector<bool>& use_domain, bool print, bool find_one, bool& timeout);
	
	//for mining
	void SetMode(bool b, int maxL);
	void SetData(DataGraph*);
	void ClearCache();

	//in build mode only
	void SetDomains(vector<range>& d, vector<bool>& u, vector<star>& r);

	double calculate_domains_time_ = 0;
	size_t calculate_domains_cnt_ = 0;

	Bitmap single_star_marked_;
	Bitmap marked_;
	Bitmap pruned_;

	vector<int> SSM_;
	vector<vector<int>> Pr_;
	//returned to miner
	vector<vector<int>> Dom_;

private:
	bool validateBounds();
	void determineFirst(); 
	void determineOrder(const int = -1); 
	void determineOrderByRank(); 
	void makeCoreFront();
	void determineUseInitialSearchSpace();
	void getDomainsStatic();
	double recurBuildEntry();
	double recurBuild(const int);
	double recurQuery(const int);
	void getLocalSampleSpaceBuild(const int);
	void getLocalSampleSpaceQuery(const int);
	void getLocalSearchSpace(const int u);
	void getInitialSearchSpace(const int, set<range>&);
	void getStarSelectivities();
	void makeRBIGraph();
	void mergeIntersectionOneVertex(int);
	void mergeIntersection();
	void printPlan();

	bool searchCache(set<range>& ranges, int** arr, int* size);
	void storeCache(set<range>& ranges, const int* arr, int size);

	bool naiveRecur(Pattern& p, vector<int>& match, int u, vector<set<int>>& Dom, clock_t start);

	vector<int> star_sel_;
	vector<int> edge_sel_;
	vector<RBI> rbi_;
	bool pattern_mining_;
	double y1_, y1_sum_;
	int sample_cnt_;
	int sample_size_;
	bool valid_;
	int last_;

	mt19937 generator_;

	vector<int> plan_, connected_plan_;
	vector<int> backward_, initial_backward_;
	vector<int> sample_; 
	vector<int> sample_space_size_;
	vector<int> num_adj_dvids_;
	vector<int> candidates_;
	unsigned long long cartesian_;

	vector<vector<AdjElem>> reduced_adj_;
	vector<set<pair<int, int>>> reduced_rel_;
	vector<vector<bool>> to_rollback_rel_;

	int intersect_capacity_[MAX_ORDER];
	int* sample_space_[MAX_ORDER];
	int* intersect_[MAX_ORDER];
	int* temp_[MAX_ORDER];
	
	int last_stopped_order_;
	const int min_y1_cnt_ = 10;
	int max_y1_cnt_;
	double sum_y1_;
	int avg_sample_cnt_per_y1_;
	int recur_call_cnt_;

	//for mining
	MinerWrapper* miner_ = NULL;
	bool in_build_ = false;
	int maxL_ = -1;
	vector<range> domains_;
	vector<bool>  use_domain_;
	vector<set<pair<int, bool>>> dom_considered_rels_;
	
	vector<int*> initial_search_space_;
	vector<int>  initial_search_space_size_;

	vector<int*>  search_space_;
	vector<int> space_size_;

	int** initial_intersect_ = NULL;
	int** initial_temp_ = NULL;
	int*  initial_intersect_capacity_ = NULL;

	//intersection cache
	//can point to domains
	enum CACHE_TYPE { EMPTY, MIN, MAT };
	struct cache_val {
		CACHE_TYPE type;
		int size;
		int *arr;
	};
	unordered_map<size_t, cache_val> intersect_cache_; 
	
	vector<int> match_;
    vector<vector<AdjElem>> exact_backward_;
	vector<int> rank_;
	vector<int> rootplan_;
	vector<vector<int>> subplans_;
	vector<vector<int>> CT_;
	vector<bool> ap_;
	vector<bool> use_initial_search_space_;
	vector<vector<int>> to_clear_;

	//for any color
	vector<vector<int>> forward_;
	vector<bool> is_calculated_;
	//for black 
	vector<bool> is_black_marked_;
};

#endif
