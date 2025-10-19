#pragma warning(disable:4786)
#ifndef _EDGE_H
#define _EDGE_H

class Edge{
public:
	Edge(){_i=0;_j=0;}
	Edge(int i, int j){_i=i;_j=j;}
	Edge(const Edge &e){_i=e._i;_j=e._j;}
	Edge	&operator=(const Edge &e){ 
		if(this!=&e){
			_i=e._i;_j=e._j;
		}
		return *this;
	}
	int _i,_j;
};

template<class Tp>
class WeightedEdge:public Edge{
public:
	WeightedEdge (){_i=0;_j=0;_weight=0;}
	WeightedEdge (int i, int j, Tp d):Edge(i,j){_weight=d;}
	bool	operator<(const WeightedEdge &e){return (_weight<e._weight);}
	WeightedEdge	&operator=(const WeightedEdge &e){ 
		if(this!=&e){
			_i=e._i;_j=e._j;_weight=e._weight;
		}
		return *this;
	}
	void SaveToLine(char * line){sprintf(line,"%d,%d,%g",_i,_j,_weight);}
	void ReadFromLine(const char * line){double w;sscanf(line,"%d,%d,%f",_i,_j,w);_weight=(Tp)w;}
	Tp _weight;
};
#endif

