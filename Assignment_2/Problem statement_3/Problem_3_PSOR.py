import numpy as np
import matplotlib.pyplot as plt
import math
from timeit import default_timer as timer

h=0.4
l=0.3
time=0

nx=31
ny=41

dx=float(l/(nx-1))
dy=float(h/(ny-1))

beta=dx/dy
b2=math.pow(beta,2)


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
#-----------------------------------------------------------------------------------------#
#----------------------------PSOR_Method--------------------------------------------------#
#-----------------------------------------------------------------------------------------#
def psor(b2,nx,ny,w):
  c=2*(1+b2)
  cc=0.01
  iterate=0
  err=[]
  ite=[]
  x=np.linspace(0,0.3,31)
  y=np.linspace(0,0.4,41)
  [Y,X]=np.meshgrid(y,x)
  
  #-----------------------------------------------------------------------------------------#
  #-----------------------Initializing the matrices-----------------------------------------#
  #-----------------------------------------------------------------------------------------#

  t_k=np.zeros((nx,ny))
  t_k1=np.zeros((nx,ny))

  #-----------------------------------------------------------------------------------------#
  #--------------------------Setting up the Boundary conditions-----------------------------#
  #-----------------------------------------------------------------------------------------#
  for i in range(nx):
    for j in range(ny):
      if j==0:
        t_k[i,j]=40
        t_k1[i,j]=40
      elif j==40:
        t_k[i,j]=40
        t_k1[i,j]=40  
  #-----------------------------------------------------------------------------------------#
  #-----------------------------Implementing the scheme-------------------------------------#
  #-----------------------------------------------------------------------------------------#

  start =timer() # Time starts
  while cc>=0.01:
    for i in range(1,nx-1):
      for j in range(1,ny-1):
        t_k1[i][j] = (1-w)*t_k[i][j] + (w/c)*((t_k[i+1][j] + t_k1[i-1][j] + b2*(t_k[i][j+1]+ t_k1[i][j-1])))
    cc= convergence_criteria(t_k,t_k1,nx,ny)
    t_k=np.copy(t_k1)
    iterate+= 1
    err=np.concatenate((err,[cc]))
    ite=np.concatenate((ite,[iterate]))

  end = timer() # time stops  
  
  e=np.argmax(err)
  print("\nTime taken to reach Steady State (Complete Domain) :"+str(end-start)+"s")
  print("Number of iterations to Steady state(Complete Domain):",iterate)
  print("Max error :",err[e])
  plt.plot(ite,err)
  plt.title("Error vs Iterations Plot\n Problem III")
  plt.xlabel('iterations')
  plt.ylabel('error')
  plt.show()
  
  # Temperature Distribution at steady state shown as a contour plot
  plt.contourf(X,Y,t_k1,cmap='magma')
  plt.colorbar()
  plt.xlabel("Length in m") 
  plt.ylabel("Width in m")
  plt.title("Temperature distribution at Steady state\n Problem III")
  plt.show()

  return iterate

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

w_opti=1.843137254901961
psor(b2,nx,ny,w_opti)



