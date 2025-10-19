#ifndef _PIXEL_GRAPH_H
#define _PIXEL_GRAPH_H
#pragma warning(disable:4786)

#include "SimpleMatrix.h"
#include "Point.h"
#include "FastGraph.h"

class PixelGraph:public FastGraph<int>{
public:
	PixelGraph		();

	void			Init(const std::vector<int> &pixx, const std::vector<int> &pixy, int nx, int ny);
	void			Init(std::vector<Pointi> &pix, int nx, int ny);
	template<class Tp>
	void	Init(SimpleMatrix<Tp> &map){	
		int x,y,nx=map.nx(),ny=map.ny(),n;
		_index.SetDimension(ny,nx);	_index.InitValue(-1);
		n=0;
		for (y=0;y<ny;y++)
		for (x=0;x<nx;x++)
			if (map(y,x)!=0){
				_index(y,x)=n;
				n++;
			}
		_points.clear(); 
		_points.resize(n);
		for (y=0;y<ny;y++)
		for (x=0;x<nx;x++){
			n=_index(y,x);
			if (n>=0) 
				_points[n].Set(x,y);
		}
	}

	void			BuildNbGraph(int nbType=8);
	virtual bool	IsNb(int x1,int y1,int x2,int y2,int nbtype);
	Pointi *		GetPoint(int i){return &_points[i];}
	int				GetIndex(int x, int y);
	void			GetNeighbors(int x, int y, std::vector<int> &N, int nbtype);
	void			GetNeighbors(Pointi *p, std::vector<int> &N, int nbtype);
	bool			CheckConsistent();
	int				GetOtherNb(int i1,std::vector<int> &N);

	// Curve Functions
	bool			GetCurveTowards(std::vector<int> &L, int i, int nbi);
	bool			GetCurveThrough(std::vector<int> &L, int i);
	bool			GetCurveThrough(std::vector<int> &L, int x, int y);
	void			GetAllCurves(std::vector<std::vector<int> > &L,int minlen);



	// union-find functions
	void			BuildNbSets(UnionFind &sets, int nbType=8);


	int				_nbType;
	SimpleMatrix<int>		_index;
	std::vector<Point<int> >  _points;
};

void AddPositiveRegions(SimpleMatrix<int> &out, SimpleMatrix<int> &m1, SimpleMatrix<int> &m2);
template<class Tp,class Tp2>
void ExtractCurves(std::vector<std::vector<Tp> > &out, SimpleMatrix<Tp2> &edge, int minSize, int nbType=6);
void MergeRegions(SimpleMatrix<int> &out, SimpleMatrix<int> &m1, SimpleMatrix<int> &m2);
void MergePositiveRegions(SimpleMatrix<int> &out, SimpleMatrix<int> &m1, SimpleMatrix<int> &m2);
void LabelRegions(SimpleMatrix<int> &out, UnionFind &U, int nx);

template<class Tp>
void Merge4NbSets(UnionFind &uf, SimpleMatrix<Tp> &m){
	int x,y,nx=m.nx(),ny=m.ny(),ynx;
	ynx=0;
	for (y=0;y<ny;y++){	
		for (x=1;x<nx;x++){
			if (m(y,x-1)==m(y,x))
				uf.SetUnion(ynx+x-1,ynx+x);
		}
		ynx+=nx;
	}
	ynx=nx;
	for (y=1;y<ny;y++){
		for (x=0;x<nx;x++){
			if (m(y-1,x)==m(y,x))
				uf.SetUnion(ynx+x-nx,ynx+x);
		}
		ynx+=nx;
	}
}

template<class Tp>
void Merge4NbSets(UnionFind &uf, SimpleMatrix<Tp> &m, float maxDiff){
	int x,y,nx=m.nx(),ny=m.ny(),ynx;
	ynx=0;
	for (y=0;y<ny;y++){	
		for (x=1;x<nx;x++){
			if (fabs((float)m(y,x-1)-m(y,x))<maxDiff)
				uf.SetUnion(ynx+x-1,ynx+x);
		}
		ynx+=nx;
	}
	ynx=nx;
	for (y=1;y<ny;y++){
		for (x=0;x<nx;x++){
			if (fabs((float)m(y-1,x)-m(y,x))<maxDiff)
				uf.SetUnion(ynx+x-nx,ynx+x);
		}
		ynx+=nx;
	}
}

template<class Tp>
void MergePos4NbSets(UnionFind &uf, SimpleMatrix<Tp> &m){
	int x,y,nx=m.nx(),ny=m.ny(),ynx;
	ynx=0;
	for (y=0;y<ny;y++){	
		for (x=1;x<nx;x++){
			if (m(y,x)>=0&&m(y,x-1)==m(y,x))
				uf.SetUnion(ynx+x-1,ynx+x);
		}
		ynx+=nx;
	}
	ynx=nx;
	for (y=1;y<ny;y++){
		for (x=0;x<nx;x++){
			if (m(y,x)>=0&&m(y-1,x)==m(y,x))
				uf.SetUnion(ynx+x-nx,ynx+x);
		}
		ynx+=nx;
	}
}

template<class Tp>
void Merge4NbSets(UnionFind &uf, SimpleMatrix<Tp> &m1, SimpleMatrix<Tp> &m2){
	int x,y,nx=m1.nx(),ny=m1.ny(),ynx;
	ynx=0;
	for (y=0;y<ny;y++){	
		for (x=1;x<nx;x++){
			if (m1(y,x-1)==m1(y,x)&&m2(y,x-1)==m2(y,x))
				uf.SetUnion(ynx+x-1,ynx+x);
		}
		ynx+=nx;
	}
	ynx=nx;
	for (y=1;y<ny;y++){
		for (x=0;x<nx;x++){
			if (m1(y-1,x)==m1(y,x)&&m2(y-1,x)==m2(y,x))
				uf.SetUnion(ynx+x-nx,ynx+x);
		}
		ynx+=nx;
	}
}

template<class Tp>
void LabelPositiveRegions(SimpleMatrix<int> &out, UnionFind &U, SimpleMatrix<Tp> &m){
	// label regions that are nonzero in m with the index of the set
	int i,j,nx=m.nx(),ny=m.ny(),n=nx*ny,setidx;
	out.SetDimension(ny,nx);out.InitValue(-1);
	setidx=0;
	for (i=0;i<n;i++)
		if (m[i]>=0){
			j=U.GetPrev(i);
			if (j<0){
				out[i]=setidx;
				setidx++;
			}
			else
				out[i]=out[j];
		}
}

template<class Tp>
void LabelPositiveRegions(SimpleMatrix<int> &out, UnionFind &U, SimpleMatrix<Tp> &m, int minSize){
	// label regions that are nonzero in m with the index of the set
	int i,j,nx=m.nx(),ny=m.ny(),n=nx*ny,setidx,lidx;
	std::vector<int> size;
	out.SetDimension(ny,nx);out.InitValue(-1);
	setidx=0;
	for (i=0;i<n;i++)
		if (m[i]>=0){
			j=U.GetPrev(i);
			if (j<0){
				out[i]=setidx;
				size.push_back(1);
				setidx++;
			}
			else{
				out[i]=out[j];
				size[out[i]]++;
			}
		}
	lidx=setidx=0;out.InitValue(-1);
	for (i=0;i<n;i++)
		if (m[i]>=0){
			j=U.GetPrev(i);
			if (j<0){
				if (size[setidx]>=minSize){
					out[i]=lidx;
					lidx++;
				}
				setidx++;
			}
			else
				out[i]=out[j];
		}
}

template<class Tp>
void LabelPositiveRegions(SimpleMatrix<int> &out, SimpleMatrix<Tp> &m){
	// label regions that are nonzero in m with the index of the set
	UnionFind U;
	U.Construct((int)m.size());
	MergePos4NbSets(U,m);
	LabelPositiveRegions(out,U,m);
}

template<class Tp,class Tp2>
void ExtractCurves(std::vector<std::vector<Tp> > &out, SimpleMatrix<Tp2> &edge, int minSize, int nbType){
	PixelGraph gr;
	std::vector<std::vector<int> > v;
	std::vector<Tp> cv;
	int i,c,n,nc;
	Tp p;
	gr.Init(edge);
	gr.BuildNbGraph(nbType);
//	gr.RemoveSpuriousPixels(100);
//	SaveImg("c:/tmp/cnt.bmp",(gr._index>-1)*255);
	gr.GetAllCurves(v,minSize);
	nc=(int)v.size();
	out.clear();
	out.reserve(nc);
	for (c=0;c<nc;c++){
		n=(int)v[c].size();
		if (n<minSize)
			continue;
		cv.resize(n);
		for (i=0;i<n;i++){
			p=*gr.GetPoint(v[c][i]);
			cv[i].Set(p);
		}
		out.push_back(cv);
	}
}

template<class Tp>
void ExtractBoundary(std::vector<Tp> &out, SimpleMatrix<int> &edge, int minSize){
	PixelGraph gr;
	Tp p;
	std::vector<std::vector<int> > v;
	int i,c,n,nc;
	gr.Init(edge);
	//gr.RemoveSpuriousPixels(100);
	//SaveImg("c:/tmp/cnt.bmp",(gr._index>-1)*255);
	gr.BuildNbGraph(8);
	gr.GetAllCurves(v,minSize);
	nc=(int)v.size();
	out.clear();
	out.reserve(nc*10);
	for (c=0;c<nc;c++){
		n=(int)v[c].size();
		if (n<minSize)
			continue;
		for (i=0;i<n;i++){
			p=*gr.GetPoint(v[c][i]);
			out.push_back(p);
		}
	}
}

template<class Tp>
void ExtractBoundary(std::vector<Tp> &out, SimpleMatrix<uchar> &map){
	SimpleMatrix<int> edge;
	ExtractBoundary(edge,map);
//	SaveImg("c:/tmp/cnt.bmp",edge*255);
	ExtractBoundary(out,edge);
}

#endif