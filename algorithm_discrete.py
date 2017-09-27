
import numpy as np
import matplotlib.pyplot as plt
import operator
import matplotlib.animation as animation
import time





n_simulations=np.arange(0,10,1)
##################### 2D implementation ###############################
#defining the parameter space with a N-dimensional grid
OmegaMin, OmegaMax=0.2,0.4
SigmaMin, SigmaMax=0.7,0.9
parameter_space_resolution=100
Omega=np.linspace(OmegaMin,OmegaMax,parameter_space_resolution)
Sigma=np.linspace(SigmaMin,SigmaMax,parameter_space_resolution)


#defining an error field
#ErrorMatrix=np.random.rand(parameter_space_resolution,parameter_space_resolution)
i=0
ErrorMatrix=np.zeros((parameter_space_resolution,parameter_space_resolution))
while i<parameter_space_resolution:
	j=0
	while j<parameter_space_resolution:
		ErrorMatrix[i,j]=np.sqrt(Omega[j]**2+Sigma[i]**2)
		j=j+1
	i=i+1

#inserting random seeds
seed=tuple(np.random.random_integers(0,parameter_space_resolution-1,2))

#defining the Delta function between two points of the parameter space
def DELTA (point):
	i=0
	D=0
	while i<parameter_space_resolution:
		j=0
		while j<parameter_space_resolution:
			D=D+np.abs(ErrorMatrix[i,j]-ErrorMatrix[point])
			j=j+1
		i=i+1
	return D	 

print 'Delta function at the vertices of the parameter space:'
print DELTA((0,0))
print DELTA((0,parameter_space_resolution-1))
print DELTA((parameter_space_resolution-1,0))
print DELTA((parameter_space_resolution-1,parameter_space_resolution-1))

i=0
DE=[]

real_time_begin=time.time()
while i<parameter_space_resolution:
	j=0
	while j<parameter_space_resolution:
		DE.append(DELTA((i,j)))
		j=j+1
	i=i+1

real_time_end=time.time()
print '  Real minimum of the Delta field:'
print min(DE)
print 'execution time:'
print (real_time_end-real_time_begin), ' s'

i=0
epsilon=0.001
zeta=0.1
no_move=0
path=[]
Deltas=[]

guess_time_begin=time.time()
while no_move<20:
#moving randomly the random seeds, with a probability to accept the move
	#seed_move=tuple(seed+tuple(np.random.random_integers(-1,1,2))
	ran_mov= tuple(np.random.random_integers(-10,10,2))
	seed_move=tuple(map(operator.add, seed,ran_mov))
#boundaries for the random point
	if seed_move[0]>=parameter_space_resolution:
		seed_move=list(seed_move)
		seed_move[0]=parameter_space_resolution-1
		seed_move=tuple(seed_move)
	if seed_move[1]>=parameter_space_resolution:
		seed_move=list(seed_move)
		seed_move[1]=parameter_space_resolution-1
		seed_move=tuple(seed_move)
	if seed_move[0]<=-1:
		seed_move=list(seed_move)
		seed_move[0]=0
		seed_move=tuple(seed_move)
	if seed_move[1]<=-1:
		seed_move=list(seed_move)
		seed_move[1]=0
		seed_move=tuple(seed_move)
#probability function for accepting the movement
	if DELTA(seed_move) < DELTA(seed):
		seed=seed_move
		no_move=0
		path.append(seed)
		Deltas.append(DELTA(seed))
	elif DELTA(seed_move) < DELTA(seed)+epsilon:
		if np.random.rand() < 0.6:
			seed=seed_move
			no_move=0
			path.append(seed)
			Deltas.append(DELTA(seed))
	elif DELTA(seed_move) < DELTA(seed)+zeta:
		if np.random.rand() < 0.2:
			seed=seed_move
			no_move=0
			path.append(seed)
			Deltas.append(DELTA(seed))
	else:
	 no_move=no_move+1
	i=i+1

guess_time_end=time.time()


path=list(path)
row=np.array([Sigma[x[0]] for x in path])
col=np.array([Omega[x[1]] for x in path])
N=len(row)
data=np.append(row,col,axis=0)
print '   Minimum found:'
print Deltas[N-1]
print 'execution time:'
print (guess_time_end-guess_time_begin), ' s'

def update_line(num, paths, line):
    line.set_data(paths[..., :num])
    return line,

fig = plt.figure()
paths=np.reshape(data,(2,N))
l, = plt.plot([], [], 'b-')
plt.xlim(SigmaMin, SigmaMax)
plt.ylim(OmegaMin, OmegaMax)
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\Omega_{m}$')
plt.pcolor(Sigma, Omega, ErrorMatrix, cmap='RdBu')
plt.colorbar()
line_ani = animation.FuncAnimation(fig, update_line, 25, fargs=(paths, l),
                                   interval=50, blit=False)
im_ani.save('animation.mp4', metadata={'artist':'Gio'})


plt.plot(np.arange(0,N),Deltas)
plt.xlabel('step')
plt.ylabel(r'$\Delta$')
plt.save('rate.png')
plt.show()

