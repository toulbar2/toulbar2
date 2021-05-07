import pytoulbar2
import math, numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle
import torch
from torchvision import datasets, transforms
import itertools
import pandas as pd
import hashlib


##########################################################################
# Image output routines
##########################################################################
def fillImage(fig,axsf,g,ph,cs):
    for v,h in enumerate(g):
        axsf[v].set_axis_off()
        if h:
            if (ph[v]):
                if h != int(cs[v]):
                    h = int(cs[v])
                    mycmap = plt.get_cmap('Purples_r')
                elif ph[v] != int(cs[v]):
                    mycmap = plt.get_cmap('Greens_r')
                else:
                    mycmap = plt.get_cmap('Greys_r')
            else:
                if h != int(cs[v]):
                    mycmap = plt.get_cmap('Purples')
                else:
                    mycmap = plt.get_cmap('Greys')
            axsf[v].imshow(MNIST_image(cs,v,h),cmap=mycmap)
        else:
            axsf[v].imshow(np.zeros((28,28)))
    fig.tight_layout(pad=0.2,h_pad=0.2,w_pad=0.2) 

# Prepare figure with flat axis for easier access
mpl.rcParams['toolbar'] = 'None' 
plt.style.use('dark_background')
figs, axss = plt.subplots(9, 9,figsize=(5,5))
axssf = axss.flatten()

##########################################################################
# Loads MNIST images and outputs on every test set image 
##########################################################################
test_set = datasets.MNIST('./data', download = True, train=False, 
                          transform=transforms.Compose([
                                    transforms.ToTensor(),
                                    transforms.Normalize((0.1307,), (0.3081,))]))
data  = iter(torch.utils.data.DataLoader(test_set, batch_size=1, shuffle=False))
images = list(map(lambda x: x[0].reshape(28,28), data))
# Load MNIST outputs and image indices for every MNIST test digit
logits = pickle.load(open("MNIST_test_marginal", "rb"))
logits_len = list(map(lambda x: len(x), logits))
img_indces = pickle.load(open("MNIST_test_indices", "rb"))

def myhash(str):
    return int(hashlib.sha512(str.encode('utf-8')).hexdigest(), 16)

def MNIST_output(cg,p,val):
    h = myhash(cg+str(p))
    return logits[val][h % logits_len[val]]

def MNIST_image(cg,p,val):
    h = myhash(cg+str(p))
    return images[img_indces[val][h % logits_len[val]]]

##########################################################################
# Sudoku grids loading 
##########################################################################
# Load grid/solution pairs from the validation set of the RRN paper
valid = pd.read_csv("valid.csv.xz",sep=",", header=None).values
hints = valid[:][:,0]
sols = valid[:][:,1]
size = math.isqrt(len(sols[0]))
par = math.isqrt(size)

def MNIST_fails(lg):
    lf = []
    for i,cg in enumerate(lg):
        mygrid = [int(h) for h in cg]
        ok = True
        for v,h in enumerate(mygrid):
            if h and (np.argmin(MNIST_output(cg,v,h)) != h):
                ok = False
        if (not ok): lf.append(i)
    return lf

##########################################################################
# Auxiliary CFN functions for Sudoku 
##########################################################################
# Adds a clique of differences with violation "cost" on "varList"
def addCliqueAllDiff(theCFN, varList, cost):
    different = (cost*np.identity(size, dtype=np.int64)).flatten()
    for vp in itertools.combinations(varList,2):
        theCFN.AddFunction(vp,different)

# Sets the value of variable with index "vIdx" to "value" using a unary function
def setHint(theCFN,vIdx,value):
    costs = theCFN.GetUB()*np.ones(size, dtype = np.int64)
    costs[value-1] = 0
    theCFN.AddFunction([vIdx], costs)

# Add a MNIST minus log prob cost on vIdx. Uncalibrated yet decent
def setProbHint(theCFN,vIdx,mlp):
    theCFN.AddFunction([vIdx], mlp[1:])


##########################################################################
# Auxiliary CFN functions for Sudoku : 6 13
##########################################################################
CP_mode = False

grid_number = 13
cgrid = hints[grid_number]
csol = sols[grid_number]
grid = [int(h) for h in cgrid]

# list of row, column and cells variable indices
rows = [ [] for _ in range(size) ]
columns = [ [] for _ in range(size) ]
cells = [ [] for _ in range(size) ]


myCFN = pytoulbar2.CFN(1) if CP_mode else pytoulbar2.CFN(1000000,6)

# create variables and keep indices in row, columns and cells 
for i in range(size):
    for j in range(size):
        vIdx = myCFN.AddVariable("X"+str(i+1)+"."+str(j+1),range(1,size+1))
        columns[j].append(vIdx)
        rows[i].append(vIdx)
        cells[(i//par)*par+(j//par)].append(vIdx)

# add the clique constraints on rows, columns and cells
for scope in rows+columns+cells:
    addCliqueAllDiff(myCFN,scope, myCFN.GetUB())    

# assign/bias variables
pgrid = []
for v,h in enumerate(grid):
    if h:
        prediction = np.argmin(MNIST_output(csol,v,h))
        pgrid.append(prediction)
        if (prediction != h):
            row = v//size
            col = v % size
            print("Erreur MNIST on cell",row+1,col+1,", a", h,"has been predicted as", prediction)
        if (CP_mode):
            setHint(myCFN,v,prediction)
        else:
            setProbHint(myCFN,v,MNIST_output(csol,v,h))
    else:
        pgrid.append(0)

sol = myCFN.Solve()

if (sol):
    fillImage(figs,axssf,sol[0],pgrid,csol)
else:
    fillImage(figs,axssf,pgrid,pgrid,csol)
    print("No solution found")

plt.show()
