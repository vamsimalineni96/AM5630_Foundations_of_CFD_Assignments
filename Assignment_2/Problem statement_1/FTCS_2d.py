import numpy as np
import matplotlib.pyplot as plt
import math
from timeit import default_timer as timer

h=0.4
l=0.3
iterate=0

nx=31
ny=41

dx=float(l/(nx-1))
dy=float(h/(ny-1))
dt=0.1

alpha=11.234e-5
k=380
ax=(alpha*dt)/(dx**2)
ay=(alpha*dt)/(dy**2)

x=np.linspace(0,0.3,31)
y=np.linspace(0,0.4,41)
[Y,X]=np.meshgrid(y,x)

#-----------------------------------------------------------------------#
#----------------------Defining the error function----------------------#
#-----------------------------------------------------------------------#
def convergence_criteria(tn,tn1,nx,ny):
  error_sum=0
  error=(np.subtract(tn1,tn))
  for i in range(ny):
    for j in range(nx):
      error_sum+=abs(error[i][j])
  return(error_sum)    

#-----------------------------------------------------------------------------------------#
#-----------------------Initializing the matrices-----------------------------------------#
#-----------------------------------------------------------------------------------------#

t_in=np.zeros((ny,nx))         # For plotting purpose (dummy)
t_n=np.zeros((ny,nx))          # Temperature matrix at nth level
t_n1=np.zeros((ny,nx))         # Tempearture matrix at n+1th level
t_inter=np.zeros((ny,nx))      # Intermediate Temperature matrix for ftcs calculation


#-----------------------------------------------------------------------------------------#
#--------------------Setting up the Boundary conditions-----------------------------------#
#-----------------------------------------------------------------------------------------#
# Defining the conditions at the corners to the initial t_n matrix :

for i in range(ny):
  for j in range(nx):
    if (i==0 and j==0):
      t_n[i,j]=20
    if(i==40 and j==0):
      t_n[i,j]=5 
    if (i==40 and j==30):
      t_n[i,j]=5 
    if (i==0 and j==30):
      t_n[i,j]=20    

# Defining the wall conditions :
for i in range(1,nx-1):
  t_n[0,i]=40
  t_n[ny-1,i]=10

t_in=t_n

plt.imshow(t_in,cmap='gnuplot2',origin='lower')
plt.colorbar()
plt.title("Initial Temperature distribution ",loc="left",fontstyle='italic')
plt.show()

#-----------------------------------------------------------------------------------------#
#-----------------------Implementing the FTCS Scheme--------------------------------------#
#-----------------------------------------------------------------------------------------#

cc=0.01
time_counter=0

start=timer()

while(cc>=0.01):
  for i in range(1,ny-1):
    for j in range(1,nx-1):
      t_inter[i,j]=ay*(t_n[i+1][j]-2*t_n[i][j]+t_n[i-1][j])+ax*(t_n[i][j+1]-2*t_n[i][j]+t_n[i][j-1])

  t_n1=np.add(t_n,t_inter)
  cc=convergence_criteria(t_n,t_n1,nx,ny)
  t_n=t_n1
  time_counter+=dt
  iterate+=1
  
  if(time_counter==9.99999999999998):
    for i in range(1,nx-1):
      if(i%5==0):
        plt.plot(y,t_n1[:,i])
        plt.title("Temperature distribution at x="+str(0.01*i)+"m at 10s",loc="left",fontstyle='italic')
        plt.xlabel("Width")
        plt.ylabel("Temperature in 0C")
        plt.show()

    xplot10s=np.transpose(t_n1)
    for i in range(1,ny-1):
      if (i%5==0):
        plt.plot(x,xplot10s[:,i])
        plt.title("Temperature distribution at y="+str(0.01*i)+"m at 10s",loc="left",fontstyle='italic')
        plt.xlabel("Length")
        plt.ylabel("Temperature in 0C")
        plt.show()
    plt.contourf(X,Y,np.transpose(t_n1),cmap='gnuplot2',origin='lower')
    plt.colorbar()
    plt.title("Temperature distribution at 10s")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()    

  if(time_counter==40.0000000000003):
    for i in range(1,nx-1):
      if(i%5==0):
        plt.plot(y,t_n1[:,i])
        plt.title("Temperature distribution at x="+str(0.01*i)+"m at 40s",loc="left",fontstyle='italic')
        plt.xlabel("Width")
        plt.ylabel("Temperature in 0C")
        plt.show()

    xplot40s=np.transpose(t_n1)
    for i in range(1,ny-1):
      if (i%5==0):
        plt.plot(x,xplot40s[:,i])
        plt.title("Temperature distribution at y="+str(0.01*i)+"m at 40s ",loc="left",fontstyle='italic')
        plt.xlabel("Length")
        plt.ylabel("Temperature in 0C")
        plt.show()
    plt.contourf(X,Y,np.transpose(t_n1),cmap='gnuplot2',origin='lower')
    plt.colorbar()
    plt.title("Temperature distribution at 40s")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()    
      


end=timer()    
print("\nTime taken to reach Steady State is :"+str(end-start)+"s")
print("Number of iterations to Steady state:",iterate)
# Steady State solution plotting 
for i in range(1,nx-1):
  if(i%5==0):
    plt.plot(y,t_n1[:,i])
    plt.title("Temperature distribution at x="+str(0.01*i)+"m at Steady State",loc="left",fontstyle='italic')
    plt.xlabel("Width")
    plt.ylabel("Temperature in 0C")
    plt.show()

xplotss=np.transpose(t_n1)
for i in range(1,ny-1):
  if (i%5==0):
    plt.plot(x,xplotss[:,i])
    plt.title("Temperature distribution at y="+str(0.01*i)+"m at Steady State",loc="left",fontstyle='italic')
    plt.xlabel("Length")
    plt.ylabel("Temperature in 0C")
    plt.show()

#plt.imshow(t_n1,cmap='gnuplot2',interpolation='bilinear',origin='lower')
plt.contourf(X,Y,np.transpose(t_n1),cmap='gnuplot2',origin='lower')
plt.colorbar()
plt.title("Temperature distribution Steady State\n FTCS Method")
plt.xlabel("Length in m")
plt.ylabel("Height in m")
plt.show()


