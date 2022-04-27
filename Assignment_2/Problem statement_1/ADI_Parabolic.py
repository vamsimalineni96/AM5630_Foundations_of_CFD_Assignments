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
dt=0.1

alpha=11.234e-5
k=380
ax=(alpha*dt)/(dx**2)
ay=(alpha*dt)/(dy**2)

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

def adi(t_n,t_n12,t_n1,nx,ny):
  cc=0.01
  iterate=0
  start=timer()
  time_counter=0

  while(cc>=0.01):

    # Explicit in Y direction and Implicit in X direction

    for j in range(ny):
      a=-ax*(np.ones(nx-1))
      b=(1+2*ax)*(np.ones(nx))
      c=-ax*(np.ones(nx-1))
      d=np.zeros(nx)

      for i in range(nx):
        if j==0:
          d[i]=(1-2*ay)*t_n[i][j]+ay*t_n[i][j+1]
        elif j==40:
          d[i]=(1-2*ay)*t_n[i][j]+ay*t_n[i][j-1]
        else:
          d[i]=(1-2*ay)*t_n[i][j]+ay*t_n[i][j+1]+ay*t_n[i][j-1]

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
      e=-ay*(np.ones(ny-1))
      f=(1+2*ay)*(np.ones(ny))
      g=-ay*(np.ones(ny-1))
      h=np.zeros(ny)

      for j in range(ny):
        if i==0:
          h[j]=(1-2*ax)*t_n12[i][j]+ax*t_n12[i+1][j]
        elif i==30:
          h[j]=(1-2*ax)*t_n12[i][j]+ay*t_n12[i-1][j]
        else:
          h[j]=(1-2*ax)*t_n12[i][j]+ay*t_n12[i-1][j]+ax*t_n12[i+1][j]

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
    
    if (time_counter == 9.99999999999998):
      for i in range(1,nx-1):
        if i%5 == 0:
          plt.plot(y,t_n1[i])
          plt.xlabel("Width in m")
          plt.ylabel("Temperature in °Celcius")
          plt.title("Temperature Distribution at t = 10 seconds\n"+"at x="+str(0.01*i)+"m  (ADI Method)")
          plt.show()
    
      S = np.transpose(t_n1)
      for j in range(1,ny-1):
        if j%5 == 0:
          plt.plot(x,S[j])
          plt.xlabel("Length in m")
          plt.ylabel("Temperature in °Celcius")
          plt.title("Temperature Distribution at t = 10 seconds\n"+"at y="+str(0.01*j)+"m  (ADI Method)")
          plt.show() 

      plt.contourf(X,Y,t_n1,cmap='gnuplot2')
      plt.colorbar()
      plt.xlabel("Length in m") 
      plt.ylabel("Width in m")
      plt.title("Temperature distribution at 10 seconds\n ADI Method")
      plt.show()
  
    elif (time_counter == 40.0000000000003):
      for i in range(1,nx-1):
        if  i%5== 0:
          plt.plot(y,t_n1[i])
          plt.xlabel("Width in m")
          plt.ylabel("Temperature in °Celcius")
          plt.title("Temperature Distribution at t = 40 seconds\n"+"at x="+str(0.01*i)+"m  (ADI Method)")
          plt.show() 

      S = np.transpose(t_n1)
      for j in range(1,ny-1):
        if j%5 == 0:
          plt.plot(x,S[j])
          plt.xlabel("Length in m")
          plt.ylabel("Temperature in °Celcius")
          plt.title("Temperature Distribution at t = 40 seconds\n"+"at y="+str(0.01*j)+"m  (ADI Method)")
          plt.show()
    
      plt.contourf(X,Y,t_n1,cmap='gnuplot2')
      plt.colorbar()
      plt.xlabel("Length in m") 
      plt.ylabel("Width in m")
      plt.title("Temperature distribution at 40 seconds\n ADI Method")
      plt.show()
  
  end=timer()

  print("\nTime taken to reach Steady State is :"+str(end-start)+"s")
  print("Number of iterations to Steady state:",iterate)
        
  
  plt.contourf(X,Y,(t_n1),cmap='gnuplot2',origin='lower')
  plt.xlabel("Length in m")
  plt.ylabel("Height in m")
  plt.colorbar()
  plt.title("Temperature Distribution at Steady State\n ADI Method")
  plt.show()

  for i in range(1,nx-1):
    if(i%5==0):
      plt.plot(y,t_n1[i,:])
      plt.title("Temperature distribution at x="+str(0.01*i)+"m at Steady State\n ADI Method (Parabolic)",loc="left",fontstyle='italic')
      plt.xlabel("Width")
      plt.ylabel("Temperature in 0C")
      plt.show()

  xplotss=np.transpose(t_n1)
  for i in range(1,ny-1):
    if (i%5==0):
      plt.plot(x,xplotss[i,:])
      plt.title("Temperature distribution at y="+str(0.01*i)+"m at Steady State\nADI Method (Parabolic)",loc="left",fontstyle='italic')
      plt.xlabel("Length")
      plt.ylabel("Temperature in 0C")
      plt.show()

  return t_n1
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#

#Uncomment this to run the code :
t1=adi(t_n,t_n12,t_n1,nx,ny)

#---------------------------------------------------------------------------------#
#------------------------Heat Flux Calculation------------------------------------#
#---------------------------------------------------------------------------------#

def heat_flux(t,k):
  t2=np.transpose(t)
  # Calculating heat flux at Bottom wall
  q_bw=np.subtract(t2[1],t2[0])*k/dy
  # Calculating heat flux at Top wall
  q_uw=np.subtract(t2[40],t2[39])*k/dy

  plt.plot(x,q_bw)
  plt.xlabel("Length in m")
  plt.ylabel("Heat flux in Watt per Metre_square")
  plt.title("Heat Transfer at Bottom Wall, y = 0")
  plt.show()
  
  plt.plot(x,q_uw)
  plt.xlabel("Length in m")
  plt.ylabel("Heat flux in Watt per Metre_square")
  plt.title("Heat Transfer at Bottom Wall, y = 40")
  plt.show()

#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#

#heat_flux(t1,k)  