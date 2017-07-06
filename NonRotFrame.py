##
##Third attempt at the project, using a non rotating frame of reference.
##
##Calcuting the asteriods movement in Cartesian co-ordinates.
##

import numpy as np
from numpy import pi
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
import pdb
from matplotlib import animation



####Functions to define the derivatives for the ode solver####


##Calculates the x-direction acceleration
def xAccel(Y,g,rp):
##  define constants
    x = Y[0]
    xdot = Y[1]
    y = Y[2]
    ydot = Y[3]
    planetPosition = Y[4]
    r = (x**2 + y**2)**0.5
##  vectors going from asteriod to planet
    xtoplanet = (math.cos(planetPosition+pi/2)*rp - x)
    ytoplanet = (math.sin(planetPosition+pi/2)*rp - y)
    rtoplanet = (xtoplanet**2 + ytoplanet**2)**0.5
    
##  calculate star's and planet's acceleration in x direction using Newtonian Gravity formula and add them up
    star = -((g[0]/r**2)*x/r)
    planet = (g[1]/rtoplanet**2*xtoplanet/rtoplanet)
##    print(x,y,planetPosition,xtoplanet,ytoplanet,rtoplanet)
    accel = star + planet
    return accel


##Calculates the x-direction acceleration
def yAccel(Y,g,rp):
##  define constants
    x = Y[0]
    xdot = Y[1]
    y = Y[2]
    ydot = Y[3]
    planetPosition = Y[4]
    r = (x**2 + y**2)**0.5
##  vectors going from asteriod to planet
    xtoplanet = (math.cos(planetPosition+pi/2)*rp - x)
    ytoplanet = (math.sin(planetPosition+pi/2)*rp - y)
    rtoplanet = (xtoplanet**2 + ytoplanet**2)**0.5
    
##  calculate star's and planet's acceleration in y direction using Newtonian Gravity formula and add them up    
    star = -((g[0]/r**2)*y/r)
    planet = (g[1]/rtoplanet**2*ytoplanet/rtoplanet)
    
    accel = star + planet
    return accel



##returns derivatives in the order [xvel,xacc,yvel,yacc,omega] for the integrator
##the omega is to calculate the next position of the planet
def derivatives(Y,t,g,rp,omega):
    der = [Y[1], xAccel(Y,g,rp), Y[3], yAccel(Y,g,rp), omega]
##    print(der[1]**2+der[3]**2)
    return der


####################################################################################################################

####Setting up and running the ode solver####


##prints the options of planetary details
def planetOptions():
    print("Choose the planet's orbital details you wish to use")
    print("Enter 'd' for the default (Jupiter)")
    print("Enter 'm' for Mars")
    print("Enter 'v' for no planet")



##Jupiter's radius and mass
Mp = 0.001
rp = 5.2


##A loop to choose which planet's details to use
planetChoice = ' '
while(planetChoice!='d'):
    planetOptions()
    planetChoice = input('')
    print()

    if planetChoice == 'm':
        rp = 1.52
        Mp = 3*10**-7
        break
    if planetChoice == 'v':
        Mp = 0
        break



##Defining constants used throughout the program
Ms = 1
G = 4*(pi**2)
g = [G*Ms,G*Mp]
omega = (g[0]/(rp)**3)**0.5
no_of_orbits = 100
sampleperorbit = 10**2
no_of_animated_orbits = 3

##time in years
t = np.linspace(0,2*pi*no_of_orbits/omega, no_of_orbits*sampleperorbit)





##initial conditions for x
##pi/3 should be stable
init_r = rp
init_angle = pi/3
x0 = (init_r)*math.sin(init_angle)
y0 = (init_r)*math.cos(init_angle)
xdot0 = -omega*y0
ydot0 = omega*x0


##Y0 is an asteriod trailing behind the planet by an angle
##Y1 is an asteriod preceeding the planet by the same angle
Y0 = [x0,xdot0,y0,ydot0,0]
Y1 = [-x0,xdot0,y0,-ydot0,0]


##solving the ode for both initial conditions
##comment and uncomment them as required
solution1 = odeint(derivatives, Y0, t, args=(g,rp,omega))
solution2 = odeint(derivatives, Y1, t, args=(g,rp,omega))



######################################################################################################################

####below are various plots of the results####


##various graphs
##can be used with any number of orbits


r = (solution1[:,0]**2+solution1[:,2]**2)**0.5
speed = (solution1[:,1]**2+solution1[:,3]**2)**0.5
rtoplanet = ((np.cos(solution1[:,4]+pi/2)*rp - solution1[:,0])**2 + (np.sin(solution1[:,4]+pi/2)*rp - solution1[:,2])**2 )**0.5
vangle = np.arctan2(solution1[:,1],solution1[:,3])
theta = np.arctan2(solution1[:,0],solution1[:,2])
plt.figure(1)
plt.plot(t,r)
plt.title('Radius of orbit')
plt.figure(2)
plt.subplot(221)
plt.plot(t,0.5*speed**2+g[0]/r + g[1]/rtoplanet)
plt.title('Energy of asteriod')
plt.subplot(222)
plt.plot(t,speed*r*np.sin(theta-vangle))
plt.title('Angular momentum')
plt.figure(3)
plt.plot(t,rtoplanet)
plt.title('distance to planet')
plt.figure(4)
plt.plot(t,speed)
plt.title('Speed')
plt.show()





####animating the orbit into an mp4 file
####DO NOT USE WITH MORE THAN A FEW ORBITS



####setting up the objects and canvas for the animation
##fig = plt.figure()
##ax = plt.axes(xlim=(-10,10),ylim=(-10,10))
##ast1_line, = ax.plot([],[],'ko',label='Asteriod')
##ast2_line, = ax.plot([],[],'wo',label='Asteriod')
##planet_line, = ax.plot([],[],'ro',label='Planet')
##star_line, = ax.plot([],[],'y*',label='Star')
##ax.legend(loc='best')
##
##
####initialising the background to a blank space
##def init():
##    ast1_line.set_data([],[])
##    ast2_line.set_data([],[])
##    planet_line.set_data([],[])
##    star_line.set_data([],[])
##    return (ast1_line,ast2_line,planet_line,star_line)
##
##
####animation function
##def animate(i):
##    x = [solution1[i,0],solution2[i,0],rp*np.cos(solution1[i,4]+pi/2),0]
##    y = [solution1[i,2],solution2[i,2],rp*np.sin(solution1[i,4]+pi/2),0]
###    print(x,'\n',y,'\n\n')
##    ast1_line.set_data([x[0],x[0]+solution1[i,1]],[y[0],y[0]+solution1[i,3]])
##    ast2_line.set_data(x[1],y[1])
##    planet_line.set_data(x[2],y[2])
##    star_line.set_data(x[3],y[3])
##    return (ast1_line,ast2_line,planet_line,star_line)
##
##
####the animator, calls the animation function sequencally
####blit=True is faster (only changes are redrawn), but not supported on Mac OS
##anim = animation.FuncAnimation(fig,animate,init_func=init,frames=(min(no_of_orbits,no_of_animated_orbits)*sampleperorbit),interval=(t[1]*100),blit=False)
##
####change the first argument to change the saved file's name
##anim.save('trial_animation.mp4',fps=30,extra_args=['-vcodec','libx264'])
##
##plt.show()
##
##

