import numpy as np
import matplotlib.pyplot as plt
import math
from timeit import default_timer as timer
import time

h=0.2
l=0.15
w=1.843137254901961

nx=16
ny=21

dx=float(l/(nx-1))
dy=float(h/(ny-1))

beta=dx/dy
b2=math.pow(beta,2)
c2=2*(1+b2)

x=np.linspace(0,0.15,16)
y=np.linspace(0,0.2,21)
[Y,X]=np.meshgrid(y,x)


#-----------------------------------------------------------------------#
#----------------------Defining the error function----------------------#
#-----------------------------------------------------------------------#
def convergence_criteria(tn,tn1,nx,ny):
  error_sum=0
  error=(np.subtract(tn1,tn))
  for i in range(nx):
    for j in range(ny):
      error_sum+=abs(error[i][j])
  return(error_sum)    

#-----------------------------------------------------------------------#
#----------------------TDMA Solving Algorithm---------------------------#
#-----------------------------------------------------------------------#
def tdma(a, b, c, d):
    
    ne = len(d) # number of equations
    a1, b1, c1, d1 = map(np.array, (a, b, c, d)) # copying vectors
    for i in range(1, ne):
        m1 = a1[i-1]/b1[i-1]
        b1[i] = b1[i] - m1*c1[i-1] 
        d1[i] = d1[i] - m1*d1[i-1]
              
    x1 = b1
    x1[-1] = d1[-1]/b1[-1]

    for l in range(ne-2, -1, -1):
        x1[l] = (d1[l]-c1[l]*x1[l+1])/b1[l]

      
    return x1
#--------------------------------------------------------------------------#
#-----------------------Initializing the matrices--------------------------#
#--------------------------------------------------------------------------#

t_k=np.zeros((nx,ny))
t_k1=np.zeros((nx,ny))
temp=np.zeros((ny,nx))

#-----------------------------------------------------------------------------------------#
#--------------------------Setting up the Boundary conditions-----------------------------#
#-----------------------------------------------------------------------------------------#

for i in range(nx):
  for j in range(ny):
    if (j==0):
      t_k[i,j]=40
      t_k1[i,j]=40

temp=t_k1.T

# Defining the Neumann Boundary Condition
a=b2*np.ones(nx-1)
b=-(c2)*np.ones(nx)
c=b2*np.ones(nx-1)
d=np.ones(nx)

for i in range(nx):
	d[i]=-2*temp[ny-2][i]
temp[ny-1]=tdma(a,b,c,d)
t_k1=temp.T

# Restoring the Boundary conditions at wall    
for i in range(nx):
  for j in range(ny):
    if (j==0):
      t_k[i,j]=40
      t_k1[i,j]=40

# Defining the Neumann boundary condition
e=b2*np.ones(ny-1)
f=-(c2)*np.ones(ny)
g=b2*np.ones(ny-1)
h=np.ones(ny)

for j in range(ny):
	h[j]=-2*t_k1[nx-2][j]
t_k1[nx-1]=tdma(e,f,g,h)

#-----------------------------------------------------------------------------------------#
#----------------Performing the operation-------------------------------------------------#
#-----------------------------------------------------------------------------------------#

cc=0.01
iterate=0
ite=[]
err=[]

start= timer() # Time starts

while cc>=0.01:
    for i in range(1,nx-1):
      for j in range(1,ny-1):
      	t_k1[i][j] = (1-w)*t_k[i][j] + (w/c2)*((t_k[i+1][j] + t_k1[i-1][j] + b2*(t_k[i][j+1]+ t_k1[i][j-1])))
    cc= convergence_criteria(t_k,t_k1,nx,ny)
    
    # Declaring boundary conditions again
    for i in range(nx):
    	for j in range(ny):
    		if j==0:
    			t_k[i,j]=40
    			t_k1[i,j]=40
    temp=np.transpose(t_k1)
    
    a=b2*np.ones(nx-1)
    b=-(c2)*np.ones(nx)
    c=b2*np.ones(nx-1)
    d=np.ones(nx)
    for i in range(nx):
    	d[i]=-2*temp[ny-2][i]
    temp[ny-1]=tdma(a,b,c,d)
    t_k1=temp.T

    for i in range(nx):
    	for j in range(ny):
    		if (j==0):
    			t_k[i,j]=40
    			t_k1[i,j]=40

    e=b2*np.ones(ny-1)
    f=-(c2)*np.ones(ny)
    g=b2*np.ones(ny-1)
    h=np.ones(ny)
    for j in range(ny):
    	h[j]=-2*t_k1[nx-2][j]
    t_k1[nx-1]=tdma(e,f,g,h)
    for i in range(nx):
    	for j in range(ny):
    		if j==0:
    			t_k[i,j]=40
    			t_k1[i,j]=40
    			
    t_k=np.copy(t_k1)
    iterate+= 1
    err=np.concatenate((err,[cc]))
    ite=np.concatenate((ite,[iterate]))

end= timer() # time stops  


e=np.argmax(err)
print("\nTime taken to reach Steady State is (Reduced Domain):"+str(end-start)+"s")
print("Number of iterations to Steady state(Reduced Domain):",iterate)
print("Max error :",err[e])
plt.plot(ite,err)
plt.title("Error vs Iterations Plot\n Problem III")
plt.xlabel('iterations')
plt.ylabel('error')
plt.show()
  
plt.contourf(X,Y,t_k1,cmap='magma')
plt.colorbar()
plt.xlabel("Length in m") 
plt.ylabel("Width in m")
plt.title("Temperature distribution at Steady state\n Quadrant Problem III")
plt.show()

