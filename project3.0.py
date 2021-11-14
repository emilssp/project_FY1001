from numpy import *
import matplotlib.pyplot as plt

########<Set up constants>#########

g=9.81
R=0.54
r=0.03
mu_k=0.15
mu_s=0.2
c=1/2

#######<--The mass cancels out in the following calculations so it's ignored in the function formulas-->########
def Pot(g,theta1,theta2,R,r):
    return g*((R+r)*(cos(theta1)-(cos(theta2))))

def Kin(v,c): 
    return (1/2*(v**2))*(c+1)

def Work(mu_k,g,theta1,theta2,R,r,sluring): #calculates the work done by friction if sluring is 0 there is no work done 
    if sluring!=0: return mu_k*(g*cos(theta2)*(theta2-theta1)*(R+r)) #by friction so the function returns 0
    else: return 0

def E_tot(U,K,W):
    return U+K+W 

def v_2(E_max,P,W, c): #calculates v^2 by finding the kinetic energy and multiplying by 2/(1+c)
    return 2*(E_max-(P+W))/(1+c)

def time(theta0,theta1,R,r,v): #finds the time from between two angles theta0 and theta1
    return ((theta1-theta0)*(R+r))/v

def acc(v0,v1,t): #finds the acceleration from the equation v1=v0+a*t
    return (v1-v0)/t


#populate array with steps-number values for theta
steps=1000
theta_start=arcsin(r/R+r)
theta_end=pi/2 #Max value for theta
theta_arr=linspace(theta_start,theta_end,steps) #array with equally spaced thetas

###########starting conditions################
E_max=Pot(g,theta_arr[0],theta_arr[-1],R,r) #maximum energy when the objects is highest and not moving
E=[E_max] #total energy equals the maximum energy 
K=[0]
U=[E_max]
v=[0]
t=[0]
sluring=0
W=0
a=[]

for i in range(1,len(theta_arr)):
    if(i==1): #first time running the loop, assuming the objects starts to roll in the very beginning of the motion
        U.append(Pot(g,theta_arr[i],pi/2,R,r))
        v.append(sqrt(v_2(E_max,U[i],W,c)))
        t.append(time(theta_arr[i-1],theta_arr[i],R,r,v[i]))
        a.append(acc(v[i-1],v[i],t[i]))
        K.append(Kin((v[-2]+a[-1]*t[-1]),c))
        E.append(E_tot(U[i],K[i],W))
    
    max_fs=(mu_s*g*cos(theta_arr[i]))/c #max static friction acceleration
    U.append(Pot(g,theta_arr[i],pi/2,R,r))
    
    if a[i-1]>max_fs:    #if the acceleration is smaller that the maximum static friction the objects starts to slip
       if sluring==0:
           sluring=theta_arr[i] #variable "sluring" equals the angle the object starts to slip
           print('Starts slipping at an angle:',sluring)
    
    #append the speed, time, acceleration, Kinetic, Potential and Total energy in the system
    #add the work done by friction between theta[i-1] and theta[i] to the total work done by friction
    v.append(sqrt(v_2(E_max,U[i],W,c)))
    t.append(time(theta_arr[i-1],theta_arr[i],R,r,v[i]))
    a.append(acc(v[i-1],v[i],t[i]))
    
    W+=Work(mu_k,g,theta_arr[i-1],theta_arr[i],R,r,sluring)
    K.append(Kin(v[-2]+a[-1]*t[-1],c))
    E.append(E_tot(U[i],K[i],W))
    
    #######################################<End conditions>############################################
    if sluring!=0:    #if the objects started slipping at some point execute this
        if (sin(theta_arr[i])-mu_k*cos(theta_arr[i]))*g*cos(theta_arr[i])<=sin(theta_arr[i])*((v[i]**2)/(R+r)):
            print("Lost contact with the surface at angle:",theta_arr[i-1],"rad")
            break
    if sluring==0:    #if the objects was rolling all the way to the point of losing contact execute this
        if g*cos(theta_arr[i])<=g*(2*cos(theta_arr[0]))/(3+c):  
            print("Lost contact with the surface at angle:", theta_arr[i-1],"rad")
            break
    ###################################################################################################

#t is array with values for the time it takes to go from one theta to the next
#this loop creates an array with the time from the beginning until it reach reach angle
plotting_time=[]
for i in range(0,len(t)):
    plotting_time.append(sum(t[0:i])+t[i])

plt.figure(0)
plt.plot(plotting_time,E,color='black',label='Total energy in the system')
plt.plot(plotting_time,v,color='green',label='Velocity')
plt.plot(plotting_time[1:],a,color='yellow',label='Acceleration')
plt.legend(loc='upper left')
plt.title('Energy, velocity and acceleration over time')
plt.xlabel('Time')
plt.savefig('Fig')
plt.show()
