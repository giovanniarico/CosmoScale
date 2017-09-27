
import numpy as np
import matplotlib.pyplot as plt
import operator
import matplotlib.animation as animation
import time
from matplotlib.pyplot import cm




##################### 2D implementation ###############################
#defining the parameter space with a N-dimensional grid
OmegaMin, OmegaMax=0.2,0.4
SigmaMin, SigmaMax=0.7,0.9


#inserting random seeds
#seed=tuple((SigmaMin+(np.abs(SigmaMax-SigmaMin))*np.random.random_sample(), OmegaMin+np.random.random_sample()*(np.abs(OmegaMax-OmegaMin)) ))
seed=tuple((0.21,0.89))
max_move=np.minimum(np.abs(SigmaMax-SigmaMin),np.abs(OmegaMax-OmegaMin))

#defining a distance function to avoid over and under-flow.
def hypot(x,y):
    z=max(np.abs(x),np.abs(y))
    r=min(np.abs(x),np.abs(y))
    return (z*np.sqrt(1+r**2))

#defining the error function between two points
def ErrFunc(point1,point2):
    ErF=hypot(point2[0]-point1[0], point2[1]-point1[1])
    return (ErF)
#avoiding division per zero
def D0(quantity):
    if quantity<0.00001:
        quantity=0.0001
    return (quantity)
#PointinSpace: checking that the point is inside the parameter space
def PiS(point):
    if (point[0] < SigmaMin):
        point=list(point)
        point[0]= SigmaMin
        point=tuple(point)
    if (point[0] > SigmaMax):
        point=list(point)
        point[0]= SigmaMax
        point=tuple(point)
    if (point[1] < OmegaMin):
        point=list(point)
        point[1]= OmegaMin
        point=tuple(point)
    if (point[1] > OmegaMax):
        point=list(point)
        point[1]=OmegaMax
        point=tuple(point)
    return(point)

#defining a plot function for the random walk
def PlotRandWalk(path,title):
    path=list(path)
    row=np.array([x[0] for x in path])
    col=np.array([x[1] for x in path])
    N=len(row)
    data=np.append(row,col,axis=0)
    
    def update_line(num, paths, line):
        line.set_data(paths[..., :num])
        return line,
    
    
    fig = plt.figure()
    paths=np.reshape(data,(2,N))
    color=iter(cm.rainbow(np.linspace(0,1,N)))
    c=next(color)
    l, = plt.plot([], [], 'b-',marker='o',c=c)
    #plt.scatter(row,col)
    plt.title(title)
    plt.xlim(SigmaMin, SigmaMax)
    plt.ylim(OmegaMin, OmegaMax)
    plt.xlabel(r'$\sigma$')
    plt.ylabel(r'$\Omega_{m}$')
    line_ani = animation.FuncAnimation(fig, update_line, N, fargs=(paths, l),
    interval=10000/N, blit=False)
    #im_ani.save('animation.mp4', metadata={'artist':'Gio'})
    
    plt.show()


#defining two different functions that find the maximum error for a given point
def RandWalk (point,resolution): #1
    #ran_mov= tuple((max_move*np.random.random_sample(), max_move*np.random.random_sample())) #random within rmax
    step_size=max_move*0.001
    ran_mov= tuple((np.random.normal(0,step_size), np.random.normal(0,step_size))) #gaussian with mean rmax
    seed_move=PiS(tuple(map(operator.add, point,ran_mov)))
    nstep=0
    path=[point]
    while nstep<resolution:
        #moving randomly the random seeds. Using just the boundaries of the space?
        #ran_mov= tuple((max_move*np.random.random_sample(), max_move*np.random.random_sample())) #random within rmax
        ran_mov= tuple((np.random.normal(0,step_size), np.random.normal(0,step_size))) #gaussian with mean rma
        seed_move2=PiS(tuple(map(operator.add, seed_move,ran_mov)))
        #modify the size of the random walk depending the fraction of the error
        #max_move=max_move*ErrFunc(seed,seed_move1)/ErrFunc(seed,seed_move2)
        #probability function for accepting the movement
        Err1=D0(ErrFunc(point,seed_move))
        Err2=D0(ErrFunc(point,seed_move2))
        if np.random.random_sample() < np.exp((Err2**2-Err1**2)/(2*step_size**2)):
            seed_move=seed_move2
            Err1=Err2
            path.append(seed_move)
        
        
        nstep+=1
        #step_size=1/(10000000*Err1**3)
        step_size=1/(1000*Err1)
    #PlotRandWalk(path,'Maxima')
    return (Err1)


def OnlyBound (point, resolution): #2
    Omega=np.linspace(OmegaMin,OmegaMax,resolution)
    Sigma=np.linspace(SigmaMin,SigmaMax,resolution)
    bound=[]
    for i in Omega:
        bound.append(tuple((SigmaMin,i)))
        bound.append(tuple((SigmaMax,i)))
    for j in Sigma:
        bound.append(tuple((j,OmegaMin)))
        bound.append(tuple((j,OmegaMax)))
    k=0
    err0=ErrFunc(point, bound[0])
    MaxErrPt=bound[0]
    while k<len(bound)-1:
        if ErrFunc(point, bound[k+1]) > err0:
                   MaxErrPt=bound[k+1]
                   err0=ErrFunc(point, bound[k+1])
        k+=1
    return (err0)

step_size=max_move*0.001
ran_mov= tuple((np.random.normal(0,step_size), np.random.normal(0,step_size)))
point_move=PiS(tuple(map(operator.add, seed,ran_mov)))
nstep=0
path=[seed]
nob_time_begin=time.time()
while nstep<1000:
    #moving randomly the random seeds.
    #ran_mov= tuple((max_move*np.random.random_sample(), max_move*np.random.random_sample())) #random within rmax
    ran_mov= tuple((np.random.normal(0,step_size), np.random.normal(0,step_size)))
    point_move2=PiS(tuple(map(operator.add, point_move,ran_mov)))
    
    #probability function for accepting the movement
    Err1=D0(OnlyBound(point_move,1000))
    Err2=D0(OnlyBound(point_move2,1000))
    if np.random.random_sample() < np.exp((Err1**2-Err2**2)/(2*step_size**2)):
        point_move=point_move2
        Err1=Err2
        path.append(point_move)
    
    step_size=Err1**3
    nstep+=1

nob_time_end=time.time()
PlotRandWalk(path, 'Maximum with Boundaries Only')

print '***********RESULTS******************'
print 'Maxima found using just boundaries:'
print 'Minimum found: ', Err1, point_move
print 'Mean value of the path:', OnlyBound(np.mean(path,axis=0),1000),np.mean(path,axis=0)
print 'Execution time:', (nob_time_end-nob_time_begin), ' s'
print 'Number of accepted steps: ', len(path)

step_size=max_move*0.001
ran_mov= tuple((np.random.normal(0,step_size), np.random.normal(0,step_size)))
point_move=PiS(tuple(map(operator.add, seed,ran_mov)))
nstep=0
path=[seed]
rw_time_begin=time.time()

while nstep<1000:
    #moving randomly the random seeds.
    #ran_mov= tuple((max_move*np.random.random_sample(), max_move*np.random.random_sample())) #random within rmax
    ran_mov= tuple((np.random.normal(0,step_size), np.random.normal(0,step_size)))
    point_move2=PiS(tuple(map(operator.add, point_move,ran_mov)))
    
    #probability function for accepting the movement
    Err1=D0(RandWalk(point_move, 1000))
    Err2=D0(RandWalk(point_move2,1000))
    #if np.random.random_sample() < np.exp((Err1**2-Err2**2)/(2*step_size**2)):
    if np.random.random_sample() < np.exp(0.5*((Err1/step_size)**2)*(1-(Err2/Err1)**2)):
        point_move=point_move2
        path.append(point_move)
        Err1=Err2

    step_size=Err1**3
    nstep+=1
rw_time_end=time.time()
PlotRandWalk(path,'Maximum with Random Walk')

print 'Maxima found using Random Walk:'
print 'Minimum found: ', Err1, point_move
print 'Mean value of the path:', OnlyBound(np.mean(path,axis=0),1000),np.mean(path,axis=0)
print 'Execution time:', (rw_time_end-rw_time_begin), ' s'
print 'Number of accepted steps: ', len(path)
print 'Error Function in the center of the parameter space:', OnlyBound((0.8,0.3),1000), '(0.8,0.3)'
print 'Error Function in the corner of the parameter space:', OnlyBound((0.9,0.4),1000), '(0.9,0.4)'





