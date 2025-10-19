#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <math.h>
#include <stdlib.h>
#include "PixelGraph.h"
#include <vector>

typedef std::vector < std::vector<int> > vvi;
typedef std::pair < std::vector<int>, std::vector<int> > pvi;
static std::vector<pvi> extract(const std::vector<int> &x, const std::vector<int> &y, int minSize, int nbtype) {


	PixelGraph gr;
	std::vector<std::vector<int> > v;
	//
	int i;

	int nx=GetMax(x)+5;
	int ny=GetMax(y)+5;

	gr.Init(x,y,nx,ny);
	gr.BuildNbGraph(nbtype);
	gr.GetAllCurves(v,minSize);
	//pybind11::array pt({2,2},0);
	int nc=(int)v.size();
	std::vector<pvi> pts(nc);
	int j=0;
	for (int c=0;c<nc;c++){
		int n=(int)v[c].size();
		if (n<minSize)
			continue;
		std::vector<int> xi(n),yi(n);
		for (i=0;i<n;i++){
			Point<int> *p=gr.GetPoint(v[c][i]);
			xi[i]=p->x();
			yi[i]=p->y();
		}		
		pvi xyi(xi, yi);
		pts[c] = xyi;
		j++;
	}
	return pts;
}

PYBIND11_MODULE(pycvs, m) {
    m.doc() = "usage: m = extractCurves(x,y,minSize,nbtype)"; // Optional module docstring
    
    m.def("extractCurves", &extract, "extracts curves");
}
