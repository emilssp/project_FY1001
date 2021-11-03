from numpy import *
import matplotlib.pyplot as plt

########<Set up constants>#########

g=9.81
R=0.54
r=0.03
mu_k=0.05
mu_s=0.2
c=1/2

#######<The mass cancels out in the following calculations so it's ignored in the function formulas>########
def Pot(g,theta1,theta2,R,r): #calculates potential energy (U=mgh/m)
    return g*((R+r)*(cos(theta1)-(cos(theta2))))

def Kin(E_max,P,W=0):#Kinetic energy equals total energy minus the change in the potential energy plus the work of friction 
    return E_max-(P+W)#Work of friction is 0 by default

def Work(mu_k,g,theta1,theta2,R,r): #work done by friction for a distance = arc between theta1 and theta2
    return mu_k*(g*cos(theta2)*(theta2-theta1)*(R+r))

def E_tot(E_max,g,theta2,R,r,W=0): #Energy in a system is conserved unless outer forces are applied. Total energy is the system's total
    return (Kin(E_max,Pot(g,theta2,pi/2,R,r),W)+Pot(g,theta2,pi/2,R,r))+W #kinetic and potential energy + the energy lost to friction W

def v_2(K,c): #calculates v^2 using the kinetic energy
    return 2*K/(1+c)

def time(theta0,theta1,R,r,v): # finds the time is the distance traveled from theta0(R+r) to theta1(R+r) over the speed v
    return ((theta1-theta0)*(R+r))/v

def acc(v0,v1,t): #finds the acceleration as the difference between the speeds v0 and v1 over time t
    return (v1-v0)/t

#populate array with values for theta
steps=100
theta_start=arcsin(r/R+r)
theta_end=pi/2 #Max value for theta
theta_arr=linspace(theta_start,theta_end,steps) #array with equally spaced thetas between theta_start and theta_end with steps-number elements

###########starting conditions################
E_max=Pot(g,theta_arr[0],theta_arr[-1],R,r)
E=[E_max]
K=[0]
U=[E_max]
v=[0]
t=[0]
sluring=0
W=0
w1=[]
a=[]

###############################################

for i in range(1,len(theta_arr)):
    if(i==1): #first time running the loop, assuming the objects starts to roll in the very begynning of the motion
        E.append(E_tot(E_max,g,theta_arr[i],R,r))
        U.append(Pot(g,theta_arr[i],pi/2,R,r))
        K.append(Kin(E_max,U[i]))
        v.append(sqrt(v_2(K[i],c)))
        t.append(time(theta_arr[i-1],theta_arr[i],R,r,v[i]))
        a.append(acc(v[i-1],v[i],t[i]))
    
    #speed, time and acceleration for distance between theta[i-1] and theta[i]
    max_fs=(mu_s*g*cos(theta_arr[i]))/c #maximum static friction acceleration
    
    if i!=1 and a[-1]<=max_fs:
        v.append(sqrt(v_2(K[i],c)))
        t.append(time(theta_arr[i-1],theta_arr[i],R,r,v[i]))
        a.append(acc(v[i-1],v[i],t[i]))
        
    elif i!=1 and a[-1]>max_fs:
        v.append(sqrt(v_2(K[i],c)))
        t.append(time(theta_arr[i-1],theta_arr[i],R,r,v[i]))
        a.append(g*(sin(theta_arr[i])-cos(theta_arr[i])*mu_k))   
    
    if a[i-1]<=max_fs:                                 #if the acceleration is smaller that the maximum static friction
        E.append(E_tot(E_max,g,theta_arr[i+1],R,r))  #the object is rolling
        U.append(Pot(g,theta_arr[i+1],pi/2,R,r))     #therefore the friction does not do any work
        K.append(Kin(E_max,U[i+1]))                  #Calculate the energies for next value of  theta
        
    if a[i-1]>max_fs:                                      #if the acceleration is smaller that the maximum static friction
        if sluring==0:
            sluring=theta_arr[i]
            print('Starts slipping at an angle:',sluring*180/pi)#the object is slipping. First iteration the object start slipping
        #w1.append(Work(mu_k,g,theta_arr[i],theta_arr[i+1],R,r))
        W+=Work(mu_k,g,theta_arr[i],theta_arr[i+1],R,r)  #store the angle of slipping in a variable
        E.append(E_tot(E_max,g,theta_arr[i],R,r,W))      #The friction does negative work
        U.append(Pot(g,theta_arr[i+1],pi/2,R,r))        
        K.append(Kin(E_max,U[i+1],W=W))                 
        
    ##################################<End conditions>############################################
    if sluring!=0:    #if the objects started slipping at some point execute this
        if (sin(theta_arr[i])-mu_k*cos(theta_arr[i]))*g*cos(theta_arr[i])<=sin(theta_arr[i])*((v[i]**2)/(R+r)):
            print("Lost contact with the surface at angle:",theta_arr[i-1]*180/pi)
            break
    if sluring==0:    #if the objects was rolling all the way to the point of losing contact execute this
        if g*cos(theta_arr[i])<=g*(2*cos(theta_arr[0]))/(3+c):  
            print("Lost contact with the surface at angle:", theta_arr[i-1]*180/pi)
            break
    ###################################################################################################
E.pop(-1)
plotting_time=[]
for i in range(0,len(t)):
    plotting_time.append(sum(t[0:i])+t[i])

print(a)
plt.figure(0)
plt.plot(plotting_time,E,color='black',label='Total energy in the system')
plt.plot(plotting_time,v,color='green',label='Velocity')
plt.plot(plotting_time[2:],a[1:],color='yellow',label='Acceleration')
plt.legend(loc='upper left')
plt.title('Energy, velocity and acceleration over time')
plt.xlabel('Time')
plt.show()