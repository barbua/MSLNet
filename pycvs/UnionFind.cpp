#pragma warning(disable:4786)

#include "UnionFind.h"

//
// UnionFind functions
//
void UnionFind::Construct(int N){
	_prev.assign(N,-1);
	_nSets=N;
}

int UnionFind::FindSet(int i){	
	//find the label of set containing i
	int j=i;
	while(_prev[j]!=-1){
	   j=_prev[j];
	}
	return j;
}

void UnionFind::SetUnion(int x,int y){	
	//updates P to reflect the union of the sets specified by x and y
	int i=FindSet(x);
	int j=FindSet(y);
	if (i>j){
		_prev[i]=j;
		_nSets--;
	}
	if (j>i){
		_prev[j]=i;
		_nSets--;
	}
}

void UnionFind::Reset(){
	Construct((int)_prev.size());
}

void UnionFind::GetSetOf( std::vector<int> &L, int idx){
	//find the elements of the set containing i and stores them in L
	//first element is the rep

	int rep,n=size();
	rep=FindSet(idx);
	L.clear();
	L.reserve(n);
 	L.push_back(rep);
	for (int i=0;i<n;i++){
		if ((i!=rep)&&(FindSet(i)==rep))
			L.push_back(i);
	}
}

int UnionFind::GetSetRep(int nComponentNumber){
	//find the repr of the k-th set
	int nComponentIdx,n=size();
	nComponentIdx=0;
	for (int i=0;i<n;i++){
		if (_prev[i]==-1){	
			if (nComponentIdx==nComponentNumber)
				return i;
			nComponentIdx++;
		}
	}
	return -1;
}
