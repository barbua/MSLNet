import numpy as np
import torch
from scipy.optimize import linear_sum_assignment
from skimage.morphology import thin,medial_axis,skeletonize
import pycvs
from skimage.draw import line

def draw(m,p):
    nr,nc=m.shape
    p=torch.round(p).numpy().astype(np.int32)
    for i in range(1,p.shape[0]):
        r,c = line(p[i-1,0],p[i-1,1],p[i,0],p[i,1])
        r1=r[(r<nr)&(c<nc)]
        c1=c[(r<nr)&(c<nc)]
        m[r1,c1]=1
    return m

def KL(mu_q,invSigma_q,mu_p,Sigma_p):
    n=mu_p.shape[0]
    d=invSigma_q.shape[1]
    mu_p=mu_p.to(device)
    Sigma_p=Sigma_p.to(device)
    mu=mu_p-mu_q
    muSimu=(mu.unsqueeze(1)@invSigma_q@mu.unsqueeze(2)).squeeze()
    #score2=torch.sum(torch.diagonal(invSigma_q@Sigma_p, dim1=-1),dim=1)
    tr=torch.sum(invSigma_q.t()*Sigma_p,dim=[1,2])
    logdets=torch.logdet(Sigma_p)+torch.logdet(invSigma_q)
    return (muSimu+tr-logdets-d)/2
    
def PPCA(cov):
    u,d,vt=torch.linalg.svd(cov)
    si=d[:,1]
    si[si<.1]=.1
    d=d[:,0]-si
    P=vt[:,0,:]
    cov=P.unsqueeze(2)@(d.unsqueeze(1)*P).unsqueeze(1)+si.unsqueeze(1).unsqueeze(2)*torch.eye(2,device=P.device).unsqueeze(0)
    return cov,P,d,si
    
def PPCA(cov,la=1):
    u,d,vt=torch.linalg.svd(cov)
    d=d[:,0]
    P=vt[:,0,:]
    cov=P.unsqueeze(2)@(d.unsqueeze(1)*P).unsqueeze(1)+la*torch.eye(2,device=P.device).unsqueeze(0)
    return cov,P,d,la
    
    
def B_distance(mx1,Sigma1,mx2,Sigma2):
    mu=mx1-mx2
    Sigma=torch.linalg.inv((Sigma2+Sigma1)/2)
    #print(Sigma.shape,mu.shape)
    Simu=Sigma@mu.unsqueeze(2)
    #print(Simu.shape)
    dis= mu.unsqueeze(1)@Simu
    return dis.squeeze()/8

def scores(mx1,mx,a):
    n=mx.shape[0]
    nr,nc=a.shape
    s=torch.zeros(n,device=mx.device)
    l=torch.zeros(n,device=mx.device)
    x0=int(mx1[0].item())
    y0=int(mx1[1].item())
    for i in range(n):
        x1=int(mx[i,0].item())
        y1=int(mx[i,1].item())
        #print(x0,y0,x1,y1,nr,nc)
        ri,ci = line(x0,y0,x1,y1)
        sc=a[ci,ri]
        #print(sc)
        s[i]=torch.sum(sc)
        l[i]=ri.shape[0]
    return s,l


def validate(cl):
    n=cl.shape[0]
    clv=torch.zeros(n).long()-1
    for i in range(n):
        # if i%2==0:
        #     cvi=cv[i//2]
        j=cl[i]
        if (i<j)&(cl[j]==i):
            clv[i]=j
            clv[j]=i
    return clv   
        
def smooth_curves(cv,lam):
    spl = make_smoothing_spline(x, y, lam=lam)
    spl_n, u_n = make_splprep([xn, yn], s=0.1)

def draw_cv(m,p,b,P):
    nr,nc=m.shape
    p=torch.round(p).numpy().astype(np.int32)
    for i in range(1,p.shape[0]):
        r,c = line(p[i-1,0],p[i-1,1],p[i,0],p[i,1])
        r1=r[(r<nr)&(c<nc)]
        c1=c[(r<nr)&(c<nc)]
        m[r1,c1]=1
    return m

def draw_cvs(m,cvs):
    for cv in cvs:
        m=draw(m,cv.t())
    return m

def load_cvs(name):
    cv=loadmat(name)
    cv=cv['cv']
    cvs=[]
    for ci in cv:
        ci=torch.tensor(ci[0].astype(np.float32))-1
        cvs.append(ci)
        #print(ci)
        #plt.plot(ci[0,:],ci[1,:],'r')
    return cvs    

def AlignP(P,cv):
    n=len(cv)
    #print(n,P.shape)
    B=torch.zeros(2*n,2)
    for i in range(n):
        ci=cv[i]
        u=ci[:,0]-ci[:,2]
        q=torch.sum(P[2*i,:]*u)
        if q<0:
            P[2*i,:]*=-1
        u=ci[:,-1]-ci[:,-3]
        q=torch.sum(P[2*i+1,:]*u)
        if q<0:
            P[2*i+1,:]*=-1
        B[2*i]=ci[:,0]
        B[2*i+1]=ci[:,-1]
    return P,B

def find_matches(d,thr):
    n=d.shape[0]
    i,j= linear_sum_assignment(d)
    cl=torch.zeros(n).long()-1
    for i in range(n):
        if d[i,j[i]]<thr:
            cl[i]=j[i]
        #print(i,j[i],d[i,j[i]])
    return cl
    
def fitpoly(d0,d1,len):
    coef=torch.zeros(4);
    coef[1]=d0;
    coef[2]=-(2*d0+d1)/len;
    coef[3]=(d0+d1)/len**2;
    return coef

def fit_poly(p0,u0,p1,u1):
    dp=p1-p0
    len=torch.sum(dp**2)**0.5
    if len==0:
        return 0
    u=dp/len
    d0x=torch.sum(u*u0)
    d0y=u[0]*u0[1]-u[1]*u0[0]
    d0=d0y/d0x
    d1x=torch.sum(u*u1)
    d1y=u[0]*u1[1]-u[1]*u1[0]
    d1=d1y/d1x
    #print(d0x,d0y,d0,d1x,d1y,d1,len)
    return fitpoly(d0,d1,len)

def sumderiv22(coef,l):
    # int_0^l (der2^2)= int_0^l 4(b+3cx)^2=4(b+3cx)^3/9c|_0^l=4/9c*((b+3cl)^3-b^3)=4l((b+3cl)^2+b(b+3cl)+b^2)/3
    b,c=coef[2:4]
    return 4*l*(b**2+b*(b+3*c*l)+(b+3*c*l)**2)/3

def bdistmat(mx,P,Sigma,maxlen,rho):
	# Distance matrix using the Bhattacharyya distance
    n=mx.shape[0]
    d=torch.zeros(n,n)+10**5
    for j in range(n):
        mx0=mx[j,:]
        P0=P[j,:]/torch.sum(P[j,:]**2)**0.5
        dp=mx-mx0
        le=torch.sum(dp**2,1)**0.5
        #print(len)
        d[j,:]=B_distance(mx0,Sigma[j,:,:],mx,Sigma)
        for i in range(n):
            if (le[i]<=1)|(le[i]>maxlen):#*(lei+lej):
                d[j,i]=10**5
                continue
            p1=torch.sum(dp[i]*P0/le[i])
            p2=torch.sum(dp[i]*P[i,:]/le[i])
            if (p1<rho)|(p2>-rho): # at most 45 degrees
                d[j,i]=10**5
    return d

def distmat(mx,P,maxlen,rho):
	# Distance matrix using the continuity measure
    n=mx.shape[0]
    d=torch.zeros(n,n)+10**5
    for j in range(n):
        mx0=mx[j,:]
        P0=P[j,:]/torch.sum(P[j,:]**2)**0.5
        dp=mx-mx0
        le=torch.sum(dp**2,1)**0.5
        #print(len)
        for i in range(n):
            if (le[i]<=1)|(le[i]>maxlen):#*(lei+lej):
                continue
            p1=torch.sum(dp[i]*P0/le[i])
            p2=torch.sum(dp[i]*P[i,:]/le[i])
            if (p1>rho)&(p2<-rho): # at most 45 degrees
                c=fit_poly(mx0,P0,mx[i,:],P[i,:])
                d[i,j]=sumderiv22(c,le[i])
    return d
    
def compute_pcas(cv,nmin=10):
    n=len(cv)
    mx=torch.zeros(2*n,2)
    sx=torch.zeros(2*n,2,2)
    for i in range(n):
        ci=cv[i]
        ni=min(ci.shape[1],nmin)
        x=ci[:,:ni].t()
        mxi=torch.mean(x,dim=0)
        mx[2*i,:]=mxi
        x1=x-mxi
        sx[2*i,:,:]=x1.t()@x1
        x=ci[:,-ni:].t()
        mxi=torch.mean(x,dim=0)
        mx[2*i+1,:]=mxi
        x1=x-mxi
        sx[2*i+1,:,:]=x1.t()@x1
    return mx,sx

def polynomial(c,x):
    n=c.shape[0]
    xi=torch.ones(x.shape)
    sum=torch.zeros(x.shape)
    for i in range(n):
        sum+=c[i]*xi
        xi*=x
    return sum

def derivative(c,x):
    n=c.shape[0]
    xi=torch.ones(x.shape)
    sum=torch.zeros(x.shape)
    for i in range(1,n):
        sum+=i*c[i]*xi
        xi*=x
    return sum

def derivative2(c,x):
    n=c.shape[0]
    xi=torch.ones(x.shape)
    sum=torch.zeros(x.shape)
    for i in range(2,n):
        sum+=i*(i-1)*c[i]*xi
        xi*=x
    return sum

def get_cvpts(p0,d0,p1,d1):
    dp=p1-p0
    le=torch.sum(dp**2)**0.5
    c=fit_poly(p0,d0,p1,d1)
    x=torch.arange(0,le+0.9)
    y=polynomial(c,x)
    u=dp/le
    v=torch.tensor([-u[1],u[0]])
    return p0.view(1,-1)+u.view(1,-1)*x.view(-1,1)+v.view(1,-1)*y.view(-1,1)

def merge_cvs(v,cv,B,P,poly=0):
    n=v.shape[0]
    todo=torch.where(v==-1)[0].tolist()
    cvs=[]
#    allj=set(range(n))
    while len(todo)>0:
        j=todo[0]
        todo.remove(j)
        cvi=[]
        #js=[]
        j0=j
        while j>=0:
        #    js.append(j)
            if j%2==0:
                cvj=cv[j//2]
                bdj,edj=P[j,:],P[j+1,:] # directions
                j+=1
            else:
                cvj=torch.flip(cv[j//2],[1])
                j-=1
                bdj,edj=P[j+1,:],P[j,:] 
            #print(j,j//2,v[j].item(),bdj,edj)
            if len(cvi)==0:
                cvi.append(cvj)
                ei=cvj[:,-1]
                bdi,edi=bdj,edj
                #print(bi,ei)
            else:
                bj=cvj[:,0]
                d=torch.sum((bj-ei)**2)**0.5
                #print(d)
                if (d>5)&poly:
                    pi=get_cvpts(ei,edi,bj,bdj).t()
                    cvi.append(pi[:,1:-1])
                    #print(ei,edi,bj,bdj)
                    #print(pi)
                cvi.append(cvj)
                ei,edi=cvj[:,-1],edj
            j0=j
            j=v[j].item()
        todo.remove(j0)
#        allj=allj.difference(js)
        #print(js)
        cvi=torch.cat(cvi,dim=1)
        cvs.append(cvi)
    #print(allj)
    return cvs

def merge_curves(cv,nit,npts,dmax,rho,thr,poly):
    #plt.imshow(a,cmap='gray')
    for it in range(nit):
        if len(cv)==1:
            break
        mx,sx=compute_pcas(cv,npts)
        cov,P,d,si=PPCA(sx,10.)
        # if it==0:       
        #     plot_cvs(cv)
        #     plt.plot(mx[:,0],mx[:,1],'.')
        P,B=AlignP(P.cpu(),cv)
        #print(P)
        #print(B)
        d=distmat(mx,P,dmax[it],rho)
        #d=bdistmat(mx,P,cov,dmax[it],rho)
        #print(d)
        cl=find_matches(d,thr)
        # print(cl)
        v=validate(cl)
        cv=merge_cvs(v,cv,B,P,poly)
    return cv

def removeshortcurves(cv,minlen):
    cv1=[]
    for cvi in cv:
        if cvi.shape[1]>minlen:
            cv1.append(cvi)
    return cv1

def plot_cvs(cvs):
    for cvi in cvs:
        plt.plot(cvi[0,:].numpy(),cvi[1,:].numpy(), linewidth=0.5)

def save_plot(name,cv,im,dpi=600):
    a=np.ones((im.shape[0],im.shape[1],3))
    for i in range(3):
        a[:,:,i]=im
    plt.imshow(a)
    plot_cvs(cv)
    plt.axis('off')
    plt.savefig(name,dpi=dpi,pad_inches=0,bbox_inches='tight')
    plt.show()

def extractCurves(m,nbtype=6):
    a=thin(m)
    #a=medial_axis(m.numpy())
    #a=skeletonize(m.numpy())
    r,c=np.where(a>0)
    cvs=pycvs.extractCurves(c,r,3,nbtype)
    out=[]
    for cv in cvs:
        xi,yi=cv
        cvi=torch.stack((torch.tensor(xi),torch.tensor(yi))).float()
        out.append(cvi)
    return out
