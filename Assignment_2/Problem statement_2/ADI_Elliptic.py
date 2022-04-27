import numpy as np
import matplotlib.pyplot as plt
import math
from timeit import default_timer as timer

h=0.4
l=0.3

nx=31
ny=41

dx=float(l/(nx-1))
dy=float(h/(ny-1))

dt=0.2

beta=dx/dy
b2=math.pow(beta,2)

alpha=11.234e-5
k=380

x=np.linspace(0,l,nx)
y=np.linspace(0,h,ny)

[Y,X]=np.meshgrid(y,x)
time_counter=0


#NOTE: ALWAYS PLOT THE CONTOUR USING THE TRANSPOSE OF t_n and so on...!!!

#-----------------------------------------------------------------------------------------#
#----------------------------------TDMA Solver--------------------------------------------#
#-----------------------------------------------------------------------------------------#
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
#------------------------------Defining the error function--------------------------------#
#-----------------------------------------------------------------------------------------#
def convergence_criteria(tn,tn1,nx,ny):
  error_sum=0
  error=(np.subtract(tn1,tn))
  for i in range(nx):
    for j in range(ny):
      error_sum+=abs(error[i][j])
  return(error_sum)    

#-----------------------------------------------------------------------------------------#
#-------------------------------Initializing the Matrices---------------------------------#
#-----------------------------------------------------------------------------------------#
t_n=np.zeros((nx,ny))
t_n12=np.zeros((nx,ny))
t_n1=np.zeros((nx,ny))
temp=np.zeros((ny,nx))


#-----------------------------------------------------------------------------------------#
#--------------------------Setting up the Boundary conditions-----------------------------#
#-----------------------------------------------------------------------------------------#
# Defining the conditions at the corners 

for i in range(nx):
  for j in range(ny):
    if (i==0 and j==0):
      t_n[i,j]=20
      t_n1[i,j]=20
      t_n12[i,j]=20
    if(i==30 and j==0):
      t_n[i,j]=20
      t_n12[i,j]=20 
      t_n1[i,j]=20  
    if (i==30 and j==40):
      t_n[i,j]=5
      t_n12[i,j]=5 
      t_n1[i,j]=5  
    if (i==0 and j==40):
      t_n[i,j]=5
      t_n12[i,j]=5
      t_n1[i,j]=5    

# Defining the wall conditions :
for i in range(1,nx-1):
  t_n[i,0]=40
  t_n[i,ny-1]=10
  t_n12[i,0]=40
  t_n12[i,ny-1]=10
  t_n1[i,0]=40
  t_n1[i,ny-1]=10

#---------------------------------------------------------------------------------#
#-------------------------------ADI Function--------------------------------------#
#---------------------------------------------------------------------------------#

def adi(t_n,t_n12,t_n1,b2,nx,ny):
  cc=0.01
  iterate=0
  start=timer()
  time_counter=0
  err=[]
  ite=[]

  while(cc>=0.01):

    # Explicit in Y direction and Implicit in X direction

    for j in range(ny):
      a=(np.ones(nx-1))
      b=-2*(1+b2)*(np.ones(nx))
      c=(np.ones(nx-1))
      d=np.zeros(nx)

      for i in range(nx):
        if j==0:
          d[i]=(-t_n[i][j+1])
        elif j==40:
          d[i]=(-t_n[i][j-1])
        else:
          d[i]=(-t_n[i][j+1]-t_n[i][j-1])

      temp[j]=np.transpose(tdma(a,b,c,d))
  
    t_n12=np.transpose(temp) 


    # Applying TDMA changes the Boundary conditions for Temperature matrix at 
    # n+1/2 th level, hence the boudnary temperature have to be reset 

    for i in range(nx):
      for j in range(ny):
        if (i==0 and j==0):
          t_n12[i,j]=20
        if(i==30 and j==0):
          t_n12[i,j]=20 
        if (i==30 and j==40):
          t_n12[i,j]=5 
        if (i==0 and j==40):
         t_n12[i,j]=5
    # Defining the wall conditions :
    for i in range(1,nx-1):
      t_n12[i,0]=40
      t_n12[i,ny-1]=10
    for j in range(1,ny-1):
      t_n12[0,j]=0
      t_n12[nx-1,j]=0      
  
    # Explicit in Y and Implicit in X ends here

    # Performing Implicit in Y direction and Explicit in X direction  

    for i in range(nx):
      e=(np.ones(ny-1))
      f=-2*(1+b2)*(np.ones(ny))
      g=(np.ones(ny-1))
      h=np.zeros(ny)

      for j in range(ny):
        if i==0:
          h[j]=(-t_n12[i+1][j])
        elif i==30:
          h[j]=(-t_n12[i-1][j])
        else:
          h[j]=(-t_n12[i-1][j]-t_n12[i+1][j])

      t_n1[i]=tdma(e,f,g,h)

    # Applying TDMA changes the Boundary conditions for Temperature matrix at 
    # n+1 th level, hence the boudnary temperature have to be reset 

    for i in range(nx):
      for j in range(ny):
        if (i==0 and j==0):
          t_n1[i,j]=20
        if(i==30 and j==0):
          t_n1[i,j]=20  
        if (i==30 and j==40):
          t_n1[i,j]=5  
        if (i==0 and j==40):
          t_n1[i,j]=5    

    # Defining the wall conditions :
    for i in range(1,nx-1):
      t_n1[i,0]=40
      t_n1[i,ny-1]=10
    for j in range(1,ny-1):
      t_n1[0,j]=0
      t_n1[nx-1,j]=0  

    # Implicit in Y and Explicit in X ends here      
    cc=convergence_criteria(t_n,t_n1,nx,ny)
    t_n=np.copy(t_n1)
    iterate+=1
    time_counter+=dt
    err=np.concatenate((err,[cc]))
    ite=np.concatenate((ite,[iterate]))
    
  
  end=timer()
  e=np.argmax(err)
  print("\nTime taken to reach Steady State is :"+str(end-start)+"s")
  print("Number of iterations to Steady state:",iterate)
  print("Max error :",err[e])
  
  plt.plot(ite,err)
  plt.title("Error vs Iterations Plot\n ADI (Elliptic Problem) Method")
  plt.xlabel('iterations')
  plt.ylabel('error')
  plt.show()
        

  t_n1=np.transpose(t_n1)
  plt.contourf(X,Y,np.transpose(t_n1),cmap='gnuplot2',origin='lower')
  plt.colorbar()
  plt.title("Temperature distribution at Steady State\n ADI (Elliptic Problem) Method")
  plt.show()

  for i in range(1,nx-1):
    if(i%5==0):
      plt.plot(y,t_n1[:,i])
      plt.title("Temperature distribution at x="+str(0.01*i)+"m at Steady State\n ADI Method",loc="left",fontstyle='italic')
      plt.xlabel("Width")
      plt.ylabel("Temperature in 0C")
      plt.show()

  xplotss=np.transpose(t_n1)
  for i in range(1,ny-1):
    if (i%5==0):
      plt.plot(x,xplotss[:,i])
      plt.title("Temperature distribution at y="+str(0.01*i)+"m at Steady State\n ADI Method",loc="left",fontstyle='italic')
      plt.xlabel("Length")
      plt.ylabel("Temperature in 0C")
      plt.show()


  return
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#

#Uncomment this to run the code :
adi(t_n,t_n12,t_n1,b2,nx,ny)
