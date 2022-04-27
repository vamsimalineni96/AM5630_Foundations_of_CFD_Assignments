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
dt=0.1

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
#-----------------------Initializing the matrices-----------------------------------------#
#-----------------------------------------------------------------------------------------#

t_k=np.zeros((nx,ny))
t_k1=np.zeros((nx,ny))

#-----------------------------------------------------------------------------------------#
#--------------------------Setting up the Boundary conditions-----------------------------#
#-----------------------------------------------------------------------------------------#
# Defining the conditions at the corners 

for i in range(nx):
  for j in range(ny):
    if (i==0 and j==0):
      t_k[i,j]=20
      t_k1[i,j]=20
      
    if(i==30 and j==0):
      t_k[i,j]=20
      t_k1[i,j]=20  
    if (i==30 and j==40):
      t_k[i,j]=5
      t_k1[i,j]=5  
    if (i==0 and j==40):
      t_k[i,j]=5
      t_k1[i,j]=5    

# Defining the wall conditions :
for i in range(1,nx-1):
  t_k[i,0]=40
  t_k[i,ny-1]=10
  t_k1[i,0]=40
  t_k1[i,ny-1]=10


#-----------------------------------------------------------------------------------------#
#----------------------------Line_Gauss_Seidal_Method-------------------------------------#
#-----------------------------------------------------------------------------------------#
def line_gauss(t_k,t_k1,b2,nx,ny):
  c=2*(1+b2)
  cc=0.01
  iterate=0
  x=np.linspace(0,0.3,31)
  y=np.linspace(0,0.4,41)
  [Y,X]=np.meshgrid(y,x)
  temp=np.zeros((ny,nx))
  err=[]
  ite=[]

  start=timer()
  while(cc>=0.01):
    # Performing implicit operation for the inner points
    for j in range(1,ny-1):
        e=np.ones(nx-1)
        f=(-c)*(np.ones(nx))
        g=np.ones(nx-1)
        h=np.zeros(nx)

        for i in range(nx):
            if j==40:
                h[i]=(-b2)*t_k1[i][j-1]
            else:
                h[i]=(-b2)*(t_k[i][j+1]+t_k1[i][j-1])

        temp[j]=(tdma(e,f,g,h))  
    
    t_k1=np.transpose(temp)

    # Restoring the boundary conditions :
    for i in range(nx):
      for j in range(ny):
        if (i==0 and j==0):
          t_k1[i,j]=20
        if(i==30 and j==0):
          t_k1[i,j]=20 
        if (i==30 and j==40):
          t_k1[i,j]=5 
        if (i==0 and j==40):
         t_k1[i,j]=5
    # Restoring the wall conditions :
    for i in range(1,nx-1):
      t_k1[i,0]=40
      t_k1[i,ny-1]=10
    for j in range(1,ny-1):
      t_k1[0,j]=0
      t_k1[nx-1,j]=0     

    cc=convergence_criteria(t_k,t_k1,nx,ny)
    t_k=np.copy(t_k1)
    iterate+=1
    err=np.concatenate((err,[cc]))
    ite=np.concatenate((ite,[iterate]))
    end=timer()
  e=np.argmax(err)

  print("Time taken to reach Steady state is :"+str(end-start)+"s")
  print("Number of iterations to Steady state:",iterate)
  print("Max error :",err[e])

  plt.plot(ite,err)
  plt.title("Error vs Iterations Plot\n Line Gauss Seidal Method")
  plt.xlabel('iterations')
  plt.ylabel('error')
  plt.show()

  plt.contourf(X,Y,(t_k1),cmap='hot')
  plt.colorbar()
  plt.title("Temperature distribution at Steady State\n Line Gauss Seidal Method")
  plt.xlabel("Length in m")
  plt.ylabel("Height in m")
  plt.show()

  for i in range(1,nx-1):
    if(i%5==0):
      plt.plot(y,t_k1[i])
      plt.title("Temperature distribution at x="+str(0.01*i)+"m at Steady State\n Line Gauss Seidal Method",loc="left",fontstyle='italic')
      plt.xlabel("Width")
      plt.ylabel("Temperature in 0C")
      plt.show()

  xplotss=np.transpose(t_k1)
  for j in range(1,ny-1):
    if (j%5==0):
      plt.plot(x,xplotss[j])
      plt.title("Temperature distribution at y="+str(0.01*j)+"m at Steady State\n Line Gauss Seidal Method",loc="left",fontstyle='italic')
      plt.xlabel("Length")
      plt.ylabel("Temperature in 0C")
      plt.show()


#-----------------------------------------------------------------------------------------#

line_gauss(t_k,t_k1,b2,nx,ny)
      



