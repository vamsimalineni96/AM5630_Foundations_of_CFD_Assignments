import numpy as np
import matplotlib.pyplot as plt
import math
from timeit import default_timer as timer

h=0.4
l=0.3
time=0

#nx=int(input("Enter the number of x gridpoints:"))
#ny=int(input("Enter the number of y gridpoints:"))
nx=31
ny=41

dx=float(l/(nx-1))
dy=float(h/(ny-1))
dt=0.1

beta=dx/dy
b=math.pow(beta,2)

x=np.linspace(0,l,nx)
y=np.linspace(0,h,ny)

[Y,X]=np.meshgrid(y,x)
#-----------------------------------------------------------------------#
#----------------------Defining the error function----------------------#
#-----------------------------------------------------------------------#
def convergence_criteria(tn,tn1,nx,ny):
  error_sum=0
  error=(np.subtract(tn1,tn))
  for j in range(ny):
    for i in range(nx):
      error_sum+=abs(error[j][i])
  return(error_sum)    

#-----------------------------------------------------------------------------------------#
#-----------------------Initializing the matrices-----------------------------------------#
#-----------------------------------------------------------------------------------------#

t_in=np.zeros((ny,nx))         # For plotting purpose (dummy)
t_k=np.zeros((ny,nx))          # Temperature matrix at kth iteration

#-----------------------------------------------------------------------------------------#
#--------------------Setting up the Boundary conditions-----------------------------------#
#-----------------------------------------------------------------------------------------#
# Defining the conditions at the corners to the initial t_n matrix :

for j in range(ny):
  for i in range(nx):
    if (i==0 and j==0):
      t_k[j,i]=20
    if(j==40 and i==0):
      t_k[j,i]=5 
    if (j==40 and i==30):
      t_k[j,i]=5 
    if (j==0 and i==30):
      t_k[j,i]=20    

# Defining the wall conditions :
for i in range(1,nx-1):
  t_k[0,i]=40
  t_k[ny-1,i]=10

#-----------------------------------------------------------------------------------------#
#----------------------------PSOR Implementation------------------------------------------#
#-----------------------------------------------------------------------------------------#


def psor(t_k,b,nx,ny,w):
  c=2*(1+b)
  cc=0.01
  iterate=0
  err=[]
  ite=[]
  t_k1=np.zeros((ny,nx),np.float)

  # Defining the corner points :  
  for j in range(ny):
      for i in range(nx):
        if (i==0 and j==0):
          t_k1[j,i]=20
        if(j==40 and i==0):
          t_k1[j,i]=5
        if (j==40 and i==30):
          t_k1[j,i]=5
        if (j==0 and i==30):
          t_k1[j,i]=20    

  # Defining the wall conditions :
  for i in range(1,nx-1):
    t_k1[0,i]=40
    t_k1[ny-1,i]=10

  start=timer()
  while (cc>=0.01):
    for j in range(1,ny-1):
      for i in range(1,nx-1):
        t_k1[j][i] = (1-w)*t_k[j][i] + (w/c)*(t_k[j+1][i] + t_k1[j-1][i] + b*(t_k[j][i+1]+ t_k1[j][i-1]))

    cc = convergence_criteria(t_k,t_k1,nx,ny)
    t_k = np.copy(t_k1)
    iterate+=1
    err=np.concatenate((err,[cc]))
    ite=np.concatenate((ite,[iterate]))

  end=timer()  
  e=np.argmax(err)
  print("\nTime taken to reach Steady State is :"+str(end-start)+"s")
  print("Number of iterations to Steady state:",iterate)
  print("Max error :",err[e])
  plt.plot(ite,err)
  plt.title("Error vs Iterations Plot\n PSOR Method")
  plt.xlabel('iterations')
  plt.ylabel('error')
  plt.show()
  
  #plt.imshow(t_k1,cmap='gnuplot2',interpolation='bilinear',origin='lower')
  plt.contourf(X,Y,np.transpose(t_k1),cmap='magma')
  plt.colorbar()
  plt.title("Temperature distribution at Steady State\n PSOR Method")
  plt.xlabel("x")
  plt.ylabel("y")
  plt.show()
  for i in range(1,nx-1):
    if(i%5==0):
      plt.plot(y,t_k1[:,i])
      plt.title("Temperature distribution at x="+str(0.01*i)+"m at Steady State\n PSOR Method",loc="left",fontstyle='italic')
      plt.xlabel("Width")
      plt.ylabel("Temperature in 0C")
      plt.show()
  xplotss=np.transpose(t_k1)
  for j in range(1,ny-1):
      if (j%5==0):
        plt.plot(x,xplotss[:,j])
        plt.title("Temperature distribution at y="+str(0.01*j)+"m at Steady State\n PSOR Method",loc="left",fontstyle='italic')
        plt.xlabel("Length")
        plt.ylabel("Temperature in 0C")
        plt.show()

  return iterate

#-----------------------------------------------------------------------------------------#
#----------------------------Optimum Relaxation Factor------------------------------------#
#-----------------------------------------------------------------------------------------#


def optimum_rf(t_k,b,nx,ny,n_points):
  a=np.linspace(1,2,n_points,endpoint=False)[1:]
  a=np.reshape(a,(len(a),1))
  res=np.zeros((len(a),1),np.int)

  for i in range(len(a)):
    res[i,0]=psor(t_k,b,nx,ny,a[i,0])

  
  plt.plot(a,res)
  plt.xlabel("Relaxation Factor")
  plt.ylabel("Number of Iterations")
  plt.title("Number of Iterations vs Relaxation Factor using "+str(n_points)+"points")
  plt.show()
  index=np.argmin(res)
  print("The Optimum Relaxation Factor is",a[index,0])

  return(a[index,0])

#-----------------------------------------------------------------------------------------#
# Optimum Relaxation factor found using 51 points
w_opti=1.843137254901961
psor(t_k,b,nx,ny,w_opti)
#optimum_rf(t_k,b,nx,ny,51)


