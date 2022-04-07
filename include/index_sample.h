#ifndef INDEX_SAMPLE_H_
#define INDEX_SAMPLE_H_

#include <random>
#include "../include/wander_join.h"

typedef vector<vector<int>> Tuple;

class IndexSample : public WanderJoin {
public:
	void Init();
	bool GetSubstructure(int subquery_index);
	double AggCard();
	
private:
	void appendTuple(int node, int i, Tuple&);
	void appendTuple(int node, int i, int c, int v, Tuple&);
	int getBoundedSpaceSize(int node, int c, int v);
	int getElementFromBoundedSpace(int node, int c, int v, int i);
	bool getBoundInfo(int node, int& c, int& v, int& size);
	bool checkNonTreeEdges(Tuple&);
	vector<Tuple> Iter(int, vector<Tuple>&);

	mt19937 generator_;
	bool ibjs_started_;
};

#endif

