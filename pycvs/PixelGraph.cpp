#include "PixelGraph.h"
using namespace std;

//
//	PixelGraph 
//

PixelGraph::PixelGraph(){
	_nbType=8;
}

void PixelGraph::Init(const std::vector<int> &x,const std::vector<int> &y, int nx, int ny){
	int i,n=(int)x.size();
	_points.resize(n);
	_index.SetDimension(ny,nx);
	_index.InitValue(-1);
	for (i=0;i<n;i++){
		_points[i].Set(x[i],y[i]);
		_index(y[i],x[i])=i;
	}
}

void PixelGraph::Init(std::vector<Pointi> &pix, int nx, int ny){
	int i,n=(int)pix.size();
	_points=pix;
	_index.SetDimension(ny,nx);
	_index.InitValue(-1);
	for (i=0;i<n;i++){
		_index(pix[i].y(),pix[i].x())=i;
	}
}

void PixelGraph::BuildNbGraph(int nbType){
	int i,j,n,nn;
	vector<int> N;
	n=(int)_points.size();
	_nbType=nbType;
	MakeGraph(n);
	for (i=0;i<n;i++){
		GetNeighbors(&_points[i],N,nbType);
		nn=(int)N.size();
		for (j=0;j<nn;j++){
			if (N[j]>i)
				InsertEdge(i,N[j],1);
		}
	}
}

int	PixelGraph::GetIndex(int x, int y){
	if (!_index.Valid (y,x))
		return -1;
	return _index(y,x);
}

bool PixelGraph::CheckConsistent(){
	int i,n=(int)_points.size();
	for (i=0;i<n;i++){
		if (_index(_points[i].y(),_points[i].x())!=i)
			return false;
	}
	SimpleMatrix<int> M;
	M=_index>-1;
	return Sum(M)==n;
}

void PixelGraph::BuildNbSets(UnionFind &sets,int nbType){
	int i,j,n,nx=_index.nx();
	Pointi *p,*q;
	n=size();
	_nbType=nbType;
	sets.Construct(n);
	for (i=0;i<n;i++){
		p=GetPoint(i);
		for (j=i+1;j<n;j++){
			q=GetPoint(j);
			if (IsNb(p->x(),p->y(),q->x(),q->y(),_nbType)){
				sets.SetUnion(i,j);
			}
		}
	}
	n=sets.NumSets();
}

void PixelGraph::GetNeighbors(Pointi *p, vector<int> &N, int _nbType){
	N.clear();
	if (p!=NULL)
		GetNeighbors(p->x(),p->y(),N,_nbType);
}

void PixelGraph::GetNeighbors(int x, int y, vector<int> &N,int _nbType){
	N.clear();
	if (x<0||y<0)
		return;
	if (x>0){
		if (y>0){
			int j=GetIndex(x-1,y-1);
			if (j>=0){
				if (IsNb(x,y,x-1,y-1,_nbType))
					N.push_back(j);
			}

		}
	}
	if (y>0){
		int j=GetIndex(x,y-1);
		if (j>=0)
			N.push_back(j);
		if (x<_index.nx()-1){
			int j=GetIndex(x+1,y-1);
			if (j>=0){
				if (IsNb(x,y,x+1,y-1,_nbType))
					N.push_back(j);
			}

		}
	}
	if (x>0){
		int j=GetIndex(x-1,y);
		if (j>=0)
			N.push_back(j);
	}
	if (x<_index.nx()-1){
		int j=GetIndex(x+1,y);
		if (j>=0)
			N.push_back(j);
	}
	if (y<_index.ny()-1){
		if (x>0){
			int j=GetIndex(x-1,y+1);
			if (j>=0){
				if (IsNb(x,y,x-1,y+1,_nbType))
					N.push_back(j);
			}

		}
		int j=GetIndex(x,y+1);
		if (j>=0)
			N.push_back(j);
	}
	if (x<_index.nx()-1){
		if (y<_index.ny()-1){
			int j=GetIndex(x+1,y+1);
			if (j>=0){
				if (IsNb(x,y,x+1,y+1,_nbType))
					N.push_back(j);
			}
		}
	}
}

bool PixelGraph::IsNb(int x1,int y1,int x2,int y2,int _nbType){
	int d1=abs(x1-x2),d2=abs(y1-y2);
	if ((d1>1)||(d2>1))				//if too far, not neighbors
		return false;
	if (_nbType==8)
		return true;
	if ((d1==0)||(d2==0))			//4 neighbors are neighbors
		return true;
	if (_nbType==4)
		return false;
	if ((GetIndex(x2,y1)<0)&&(GetIndex(x1,y2)<0)) //8 neighbors are neighbors if no 4 neighbors nearby
		return true;
	return false;
}

int	 PixelGraph::GetOtherNb(int i1,std::vector<int> &N){
	if (N.size()!=2)
		return -1;
	if (N[0]==i1)
		return N[1];
	return N[0];
}

bool PixelGraph::GetCurveThrough(vector<int> &C, int x, int y){
	vector<int> L;
	int i=GetIndex(x,y);
	if (i<0)
		return false;
	return GetCurveThrough(C,i);
}

bool PixelGraph::GetCurveTowards(vector<int> &C, int i, int nbi){
	int j,pj,oj;
	vector<int> N;
	C.clear();
	j=nbi;pj=i;
	C.push_back(i);
	while (j>=0){
		C.push_back(j);
		GetNeighborNodes(j,N);
		oj=j;
		j=GetOtherNb(pj,N);
		if (j==i)
			j=-1;
		pj=oj;
	}
	return true;
}
bool PixelGraph::GetCurveThrough(vector<int> &C, int i){
	int j,n1,n2;
	vector<int> C1,C2,N;
	GetNeighborNodes(i,N);
	if (N.size()>2||N.size()==0){
		C.clear();
		return false;
	}
	GetCurveTowards(C1,i,N[0]);
	if (N.size()==1)
		C=C1;
	else{
		GetCurveTowards(C2,i,N[1]);
		n1=(int)C1.size();n2=(int)C2.size();
		C.clear();
		C.resize(n1+n2-1);
		for (j=0;j<n1;j++)
			C[j]=C1[n1-j-1];
		for (j=1;j<n2;j++)
			C[j+n1-1]=C2[j];
	}
	return true;
}

void PixelGraph::GetAllCurves(std::vector<std::vector<int> > &out,int minlen){
	int i,n;
	vector<int> use,C;

	out.clear();
	n=(int)size();
	use.assign(size(),1);
	for (i=0;i<n;i++)
	if (use[i]==1&&GetDegree(i)<=2){
		GetCurveThrough(C,i);
		if (C.size()>=minlen)
			out.push_back(C);
		for (int j=0;j<(int)C.size();j++)
			use[C[j]]=0;
	}
}


void LabelRegions(SimpleMatrix<int> &out, UnionFind &U, int nx){
	// label regions that are nonzero in m with the index of the set
	int i,j,ny=U.size()/nx,n=nx*ny,setidx;
	out.SetDimension(ny,nx);out.InitValue(-1);
	setidx=0;
	for (i=0;i<n;i++){
		j=U.GetPrev(i);
		if (j<0){
			out[i]=setidx;
			setidx++;
		}
		else
			out[i]=out[j];
	}
}

void MergeRegions(SimpleMatrix<int> &out, SimpleMatrix<int> &m1, SimpleMatrix<int> &m2){
	UnionFind uf;
	uf.Construct((int)m1.size());
	Merge4NbSets(uf,m1,m2);
	LabelRegions(out,uf,m1.nx());
}

void MergePositiveRegions(SimpleMatrix<int> &out, SimpleMatrix<int> &m1, SimpleMatrix<int> &m2){
	int x,y,ynx,nx=m1.nx(),ny=m1.ny();
	int i,n=nx*ny;
	UnionFind uf;
	SimpleMatrix<int> bz;

	uf.Construct(nx*ny);
	bz.SetDimension(m1);bz.InitValue(1);
	for (i=0;i<n;i++)
		if (m1[i]<0&&m2[i]<0)
			bz[i]=-1;
	ynx=0;
	for (y=0;y<ny;y++){	
		for (x=1;x<nx;x++){
			if (bz(y,x)>0&&m1(y,x-1)==m1(y,x)&&m2(y,x-1)==m2(y,x))
				uf.SetUnion(ynx+x-1,ynx+x);
		}
		ynx+=nx;
	}
	ynx=nx;
	for (y=1;y<ny;y++){
		for (x=0;x<nx;x++){
			if (bz(y,x)>0&&m1(y-1,x)==m1(y,x)&&m2(y-1,x)==m2(y,x))
				uf.SetUnion(ynx+x-nx,ynx+x);
		}
		ynx+=nx;
	}
	
	LabelPositiveRegions(out,uf,bz);
}

void AddPositiveRegions(SimpleMatrix<int> &out, SimpleMatrix<int> &m1, SimpleMatrix<int> &m2){
	// keep regions of m1 and add regions of m2 where m1 is zero
	int x,y,ynx,nx=m1.nx(),ny=m1.ny();
	int i,n=nx*ny;
	UnionFind uf;
	SimpleMatrix<int> bz;

	bz.SetDimension(ny,nx);bz.InitValue(-1);
	uf.Construct(nx*ny);
	MergePos4NbSets(uf,m1);
	for (i=0;i<n;i++)
		if (m1[i]>=0||m2[i]>=0)
			bz[i]=1;
	ynx=0;
	for (y=0;y<ny;y++){	
		for (x=1;x<nx;x++){
			if (m2(y,x)>=0&&m1(y,x)<0&&m1(y,x-1)==m1(y,x)&&m2(y,x-1)==m2(y,x))
				uf.SetUnion(ynx+x-1,ynx+x);
		}
		ynx+=nx;
	}
	ynx=nx;
	for (y=1;y<ny;y++){
		for (x=0;x<nx;x++){
			if (m2(y,x)>=0&&m1(y,x)<0&&m1(y-1,x)==m1(y,x)&&m2(y-1,x)==m2(y,x))
				uf.SetUnion(ynx+x-nx,ynx+x);
		}
		ynx+=nx;
	}
	LabelPositiveRegions(out,uf,bz);
}
