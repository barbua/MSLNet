#pragma warning(disable:4786)
#ifndef _UNION_FIND_H
#define _UNION_FIND_H

#include <vector>

class UnionFind{
public:

	UnionFind(){_nSets=0;}

	void	Construct(int N);
	int		FindSet(int i);	//find the label of set containing i
	void	GetSetOf(std::vector<int> &L, int i);
	int		GetSetRep(int nComponentNumber);
	bool	IsRep(int i){return _prev[i]==-1;}
	int		NumSets(){return _nSets;}
	int		size(){return (int)_prev.size();}
	void	SetUnion(int i, int j);	//updates P to reflect the union of the sets specified by i and j
	void	Reset();
	int		GetPrev(int i){return _prev[i];}

protected:
	int		_nSets;
	std::vector<int> _prev;		//union-find data structure
};

#endif
