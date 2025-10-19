import torch
import torch.nn as nn
import torch.nn.functional as Fn
from torchvision.models import resnet152, ResNet152_Weights
from time import time, strftime
from datetime import datetime

def get_y0(y,w):
    i,r,c=torch.where(y>0)
    r,c=r//w,c//w
    a=torch.zeros(y.shape[0],y.shape[1]//w,y.shape[2]//w,device=y.device)
    a[i,r,c]=1
    return a

def get_patches(y0,y,w=16):
    i1,r1,c1=torch.where(y0>0)
    n=r1.shape[0]
    c=torch.zeros(y.shape[0],w*w,y.shape[1]//w,y.shape[2]//w,device=y.device)
    for i in range(n):
        p=y[i1[i],r1[i]*w:(r1[i]+1)*w,c1[i]*w:(c1[i]+1)*w]
        c[i1[i],:,r1[i],c1[i]]=p.reshape(-1)
    return c

def get_allpatches(y,w=16):
    n,nr,nc=y.shape[0],y.shape[1]//w,y.shape[2]//w
    out=torch.zeros(n,w*w,nr,nc,device=y.device)
    for i in range(n):
        for r in range(nr):
            for c in range(nc): 
                p=y[i,r*w:(r+1)*w,c*w:(c+1)*w]
                out[i,:,r,c]=p.reshape(-1)
    return out

def set_nov_patches(xr,nx):
    n=xr.shape[0] # number of images
    nr=xr.shape[1] # number of rows
    nc=xr.shape[2] # number of cols

    xr=xr.reshape(xr.shape[0],xr.shape[1],xr.shape[2],nx,nx)
    x=torch.zeros(n,nr*nx,nc*nx,device=xr.device)
    for r in range(nr):
        for c in range(nc): # controls the patch
            x[:,r*nx:r*nx+nx,c*nx:c*nx+nx]=xr[:,r,c,:,:] 
    return x

def set_nov_patches2(xr,nx):
    n=xr.shape[0] # number of images
    nr=xr.shape[1] # number of rows
    nc=xr.shape[2] # number of cols

    xr=xr.reshape(xr.shape[0],xr.shape[1],xr.shape[2],nx,nx)
    x=torch.zeros(n,2,nr*nx,nc*nx,device=xr.device)
    for r in range(nr):
        for c in range(nc): # controls the patch
            x[:,1,r*nx:r*nx+nx,c*nx:c*nx+nx]=xr[:,r,c,:,:] 
    return x

class MSLNet(nn.Module):
    def __init__(self,device):
        super(MSLNet, self).__init__()
        self.net=resnet152(weights=ResNet152_Weights.DEFAULT)
        cv=self.net.conv1.weight.data[:,0,:,:].unsqueeze(1)
        self.net.conv1=nn.Conv2d(1,64,kernel_size=(7,7),stride=(2,2),padding=(3,3),bias=False).to(device)
        self.net.conv1.weight.data=cv
        self.net.fc0=nn.Conv2d(2048+1024,2,1).to(device) # coarse prediction
        self.net.fc=nn.Conv2d(2048+1024,256,1).to(device) # fine prediction
        
    def parameters(self):
        return self.net.parameters()
        
    def forward(self, x):
        x = self.net.conv1(x)
        x = self.net.bn1(x)                ### Base Block of ResNet18 : Input has 3 channels
        x = self.net.relu(x)
        x = self.net.maxpool(x)
    
        x = self.net.layer1(x)
        x = self.net.layer2(x)
        x = self.net.layer3(x)             ##### the four layers
        x1 = self.net.layer4(x)
        x2 = Fn.interpolate(x1,size=(x.shape[2],x.shape[3]),mode='bilinear')
        x=torch.cat((x,x2),dim=1)
        xr=self.net.fc(x)
        return self.net.fc0(x),xr

def merge(a):
	# make an rgb image with channels from the list of images a
    x=torch.zeros(a[0].shape[0],a[0].shape[1],3)
    for i in range(len(a)):
        x[:,:,i]=a[i]
    return x

def log2file(file,*args, also_print_to_console=True, add_timestamp=True):
    timestamp = int(time())
    dt_object = datetime.fromtimestamp(timestamp)

    if add_timestamp:
        args = (f"{dt_object}:", *args)

    successful = False
    max_attempts = 5
    ctr = 0
    while not successful and ctr < max_attempts:
        try:
            with open(file, 'a+') as f:
                for a in args:
                    f.write(str(a))
                    f.write(" ")
                f.write("\n")
            successful = True
        except IOError:
            print(f"{datetime.fromtimestamp(timestamp)}: failed to log: ", sys.exc_info())
            sleep(0.5)
            ctr += 1
    if also_print_to_console:
        print(*args)
