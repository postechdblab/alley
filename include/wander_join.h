#ifndef WANDER_JOIN_H_
#define WANDER_JOIN_H_

#include <random>
#include "../include/estimator.h"

class WanderJoin : public Estimator {
public:
	//build mode
	void PrepareSummaryStructure(DataGraph&, double); 
	void WriteSummary(const char*); 

	//query mode
	void ReadSummary(const char*); 
	void Init();
	int  DecomposeQuery(); 
	bool GetSubstructure(int); 
	double EstCard(int); 
	double AggCard();
	double GetSelectivity();
	
protected:
	void generateWalkPlans();
	void generateWalkPlans(int);
#ifdef LEAF_SCANNING
	int getNumGoodLeafTuples(vector<pair<int, int>>&, int, int, int);
#endif
	bool checkBoundedVertices(int, vector<int>);
	bool checkLabelStatistics(int);
	bool checkNonTreeEdges(vector<pair<int, int>>&);
	int  sampleTuple(int);
	int  sampleTuple(int, int, int);

	int offset_; //# vertex labels in query
	bool plans_generated_, plan_chosen_;
	int walk_size_;
	int sample_size_, sample_cnt_;

	//each node represents a table in relational model
	//each edge (with an edge label) and a vertex (with a vertex label) corresponds to a node
	map<int, int> node_to_v_; //node id to query vertex id
	vector<bool> visited_;
	//(from_[i], to_[i]) is the i'th pair of joinable nodes
	vector<pair<int, int>> join_from_, join_to_; 

	int pos_;
	vector<pair<int, int>> counterpart_, plan_;
	vector<vector<pair<int, int>>> counterparts_, plans_;
	vector<int> success_cnt_;
	vector<vector<double>> est_;
	vector<vector<int>> num_idx_lookup_;

#ifdef TRIAL_WALKS
	struct TrialInfo {
		int cnt;
		double var;
		double mean;

		TrialInfo() {}
		TrialInfo(int c, double v, double m) : cnt(c), var(v), mean(m) {}
	};
	vector<TrialInfo> trial_info_;
	int z_opt_index_;
	TrialInfo y_info_;
	void adjustZ();
#endif

	vector<vector<int>> sampled_tuples_;
	bool valid_;
	double inv_prob_;
};

#endif

