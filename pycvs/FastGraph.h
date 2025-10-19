#pragma warning(disable:4786)
#ifndef _FAST_GRAPH_H
#define _FAST_GRAPH_H

#include <vector>
#include <algorithm>
#include "VectorUtil.h"
#include "UnionFind.h"
#include "Edge.h"
#include "SimpleMatrix.h"


template<class Tp>
class FastGraph{
public:

	typedef  std::vector<std::pair<int, Tp> > EdgeList;
	
	std::vector<EdgeList> _edges;	//adjacency list with the Edges

	FastGraph(){}
	~FastGraph(){}
	
	int		size(){return (int)_edges.size();}
	virtual void MakeGraph(int nv);				//allocates memory for a graph with nv vertices
	void	Construct(SimpleMatrix<Tp> &adjMat);
	int		GetDegree(int i);
	void	InsertDirectedEdge(int i, int j, Tp w);
	void	InsertEdge(int i, int j, Tp w);
	Tp*		GetEdge(int i, int j);
	virtual void	GetNeighborNodes(int i,std::vector<int> &L);
	void	GetNbsAndWeights(std::vector<int> &nbs, std::vector<double> &wts, int i);
	void	GetWeightedEdges(std::vector<WeightedEdge<Tp> > &edge);
	void	GetWeightedEdges(std::vector<WeightedEdge<Tp> *> &edge);
	void	GetDirectedEdges(std::vector<WeightedEdge<Tp> > &edges);
	void	GrowConnectedComponent(std::vector<int> &out, int i, int radius);
	bool	Print(char * filename);
};

template<class Tp>
void ConnectedComponents(UnionFind &uf, std::vector<WeightedEdge<Tp> > &m, int nNodes){
	int i,n=(int)m.size();
	uf.Construct(nNodes);
	for (i=0;i<n;i++)
		uf.SetUnion(m[i]._i,m[i]._j);
}

template<class Tp>
void GetConnectedComponents(std::vector<std::vector<int> > &out, FastGraph<Tp> &graph){
	UnionFind uf;
	std::vector<WeightedEdge<int> > e;
	std::vector<int> idx;
	int i,n=(int)graph.size(),j,nn;
	graph.GetWeightedEdges(e);
	ConnectedComponents(uf,e,n);
	out.assign(uf.NumSets(),idx);
	idx.assign(n,-1);
	nn=0;
	for (i=0;i<n;i++){
		j=uf.FindSet(i);
		if (idx[j]==-1){
			idx[j]=nn;
			nn++;
		}
		else
			idx[i]=idx[j];
	}
	for (i=0;i<n;i++)
		out[idx[i]].push_back(i);
}

template<class Tp>
void MinimumSpanningTree(FastGraph<Tp> &out, std::vector<WeightedEdge<Tp> > &in, int nVertices){
	//finds MST for graph having edges in, stores it in out
	int u,v;
	typename std::vector< WeightedEdge<Tp> >::iterator ei;
	UnionFind uf;
	bool first=true;
	std::sort(in.begin(),in.end());
	out.MakeGraph(nVertices);
	uf.Construct(nVertices);
	for(ei=in.begin();ei!=in.end();++ei){
		u=(*ei)._i;v=(*ei)._j;
		if (uf.FindSet(u)!=uf.FindSet(v)){
			out.InsertEdge((*ei)._i,(*ei)._j,(*ei)._weight);
			uf.SetUnion(u,v);
		}
	}
}

template<class Graph>
void BuildNbGraph(Graph &g){
	int i,j,n,nn;
	std::vector<int> N;
	n=g.size();
	g.MakeGraph(n);
	for (i=0;i<n;i++){
		g.GetNeighbors(N,i);
		nn=(int)N.size();
		for (j=0;j<nn;j++){
			if (N[j]>i)
				g.InsertEdge(i,N[j],1);
		}
	}
}


//
// FastGraph functions
// 
template<class Tp>
void FastGraph<Tp>::MakeGraph(int nv){
	EdgeList e;
	_edges.assign(nv,e);
};

template<class Tp>
void FastGraph<Tp>::Construct(SimpleMatrix<Tp> &adjMat){
	int i,j,nv=adjMat.nx();
	EdgeList e;
	_edges.assign(nv,e);
	for (i=0;i<nv;i++)
	for (j=0;j<nv;j++)
		if(i!=j&&adjMat(i,j)!=0)
			InsertDirectedEdge(i,j,adjMat(i,j));
}

template<class Tp>
Tp *FastGraph<Tp>::GetEdge(int i, int j){
	//returns the weight of the Edge (i,j), 0 if there is no such Edge
	int ei=FindKey(_edges[i],j);
	if (ei>=0)
		return &_edges[i][ei].second;
	return 0;
}

template<class Tp>
int	FastGraph<Tp>::GetDegree(int i){
	//the number of Edges at vertex i
	return (int)_edges[i].size();
}


template<class Tp>
void FastGraph<Tp>::InsertDirectedEdge(int i, int j, Tp w){
	//inserts Edge without checking
	std::pair<int, Tp> ee;
	ee.first=j;ee.second=w;
	_edges[i].push_back(ee);
}

template<class Tp>
void FastGraph<Tp>::InsertEdge(int i, int j, Tp w){
	//inserts Edge without checking
	InsertDirectedEdge(i,j,w);
	InsertDirectedEdge(j,i,w);
}
template<class Tp>
void FastGraph<Tp>::GetNeighborNodes(int i,std::vector<int> &L){
	typename EdgeList::iterator ei;
	L.clear();
	for (ei=_edges[i].begin();ei!=_edges[i].end();++ei)
		L.push_back((*ei).first);	
}

template<class Tp>
void FastGraph<Tp>::GetNbsAndWeights(std::vector<int> &nbs, std::vector<double> &wts, int i){
	GetNeighborNodes(i,nbs);
	wts.assign(nbs.size(),1);
}

template<class Tp>
void FastGraph<Tp>::GetWeightedEdges(std::vector<WeightedEdge<Tp> > &edges){
	int i,j,n=size();
	WeightedEdge<Tp> edge;
	typename EdgeList::iterator ei;
	edges.clear();
	for (i=0;i<n;i++){
		for (ei=_edges[i].begin();ei!=_edges[i].end();++ei){
			j=(*ei).first;
			if (j>i){
				edge._i=i;
				edge._j=j;
				edge._weight=(*ei).second;
				edges.push_back(edge);
			}
		}
	}
}

template<class Tp>
void FastGraph<Tp>::GetWeightedEdges(std::vector<WeightedEdge<Tp> *> &edges){
	int i,j,n=size();
	WeightedEdge<Tp> *edge;
	typename EdgeList::iterator ei;
	edges.clear();
	for (i=0;i<n;i++){
		for (ei=_edges[i].begin();ei!=_edges[i].end();++ei){
			j=(*ei).first;
			if (j>i){
				edge=new WeightedEdge<double> (i,j,(*ei).second);
				edges.push_back(edge);
			}
		}
	}
}

template<class Tp>
void FastGraph<Tp>::GetDirectedEdges(std::vector<WeightedEdge<Tp> > &edges){
	int i,j,n=size();
	WeightedEdge<Tp> edge;
	typename EdgeList::iterator ei;
	edges.clear();
	for (i=0;i<n;i++){
		for (ei=_edges[i].begin();ei!=_edges[i].end();++ei){
			j=(*ei).first;
			edge._i=i;
			edge._j=j;
			edge._weight=(*ei).second;
			edges.push_back(edge);
		}
	}
}

template<class Tp>
void FastGraph<Tp>::GrowConnectedComponent(std::vector<int> &out, int i, int radius){
	static std::vector<int> visited;
	std::vector<int> N;
	int j,n=size(),l;
	if (visited.empty())
		visited.assign(n,-1);
	l=0;j=0;
	out.push_back(i);
	visited[i]=i;
	while (l<radius&&j<(int)out.size()){
		n=(int)out.size();
		for (;j<n;j++){// process the next batch
			GetNeighborNodes(out[j],N);
			for (int k=0;k<(int)N.size();k++)
				if (visited[N[k]]!=i){
					out.push_back(N[k]);
					visited[N[k]]=i;
				}
		}
		l++;
	}
}

template<class Tp>
bool FastGraph<Tp>::Print(char * filename){
	FILE *f=fopen(filename, "wt");
	if (f==0)
		return false;
	int i,j,n=(int)_edges.size();
	for (i=0;i<n;i++){
		fprintf(f,"%d:",i);
		for (j=0;j<(int)_edges[i].size();j++)
			fprintf(f," %d,%d",_edges[i][j].first,(int)_edges[i][j].second);
		fprintf(f,"\n");
	}
	fclose(f);
	return true;
}

#endif
