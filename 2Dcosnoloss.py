##Calculates the array pattern for a 1D array of cosine like elements
##
##theta = 0, and phi = 0 is the maximum of the array
##The first element is the LH most element

import numpy as np
import math
import pdb
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from math import pi

##gives the element pattern for each individual element
def elementFunction(theta,phi,n):
    return math.cos(theta)*math.cos(phi)

##gives the losses for each individual element
def elementLoss(theta,phi,n):
    return 1

##gives the array pattern
def arrayPattern(theta,phi,N,ycoords,zcoords,wavelength):
    distance = ( ycoords**2 + zcoords**2 )**0.5
    if theta == 0:
        perpdistance = distance
    else:
        b = math.atan2(-1 / (math.sin(theta)) , 1 / (math.cos(theta) * math.sin(phi)) )
        a = np.arctan2(zcoords,ycoords)
        perpdistance = distance * np.sin(a+b)
    
    elements = elementFunction(theta,phi,0) * elementLoss(theta,phi,0) * np.exp(1j*2*pi*perpdistance*math.sin(math.acos(math.cos(phi)*math.cos(theta)))/wavelength)
    return np.sum(elements)

############################################################################
##define the variables of the problem

freq = 1*10**12
##NOTE wavelength and coords must be in mm
wavelength = 3*(10**11) / freq

##desired angle in radians:
##yzplane is the direction of incident radian at an angle measured from the z axis
##xangle is the angle of inclination wrt to the x axis
xangle = 0.7
yzplane = 0.3
y=0
z=0
if xangle==0:
    y = wavelength
    z = wavelength
else:
    if yzplane==0:
        y=wavelength
    else:
        y = wavelength*math.sin(yzplane)/(math.sin(xangle))
    if yzplane==pi/2:
        z=wavelength
    else:
        z = wavelength*math.cos(yzplane)/(math.sin(xangle))


##uses 2D coords perpendicular to the theta=phi=0, theta is parallel to either y or z
dimensionsofarray = [5,5]
ycoords = np.empty((dimensionsofarray[0]*dimensionsofarray[1]))
zcoords = np.empty((dimensionsofarray[0]*dimensionsofarray[1]))

for i in range(dimensionsofarray[0]):
    for j in range(dimensionsofarray[1]):
        ycoords[dimensionsofarray[1]*i + j] = i*y
        zcoords[dimensionsofarray[1]*i + j] = j*z

N = ycoords.size

##############################################################################
##solve the problem

##print(arrayPattern(0,0,N,ycoords,zcoords,wavelength))

theta = np.linspace(-1,1,60)
phi = np.linspace(-1,1,60)

answer = np.empty((theta.size,phi.size))
answer2 = np.empty((theta.size,1))
coords = np.empty((theta.size*phi.size,3))

for i in range(theta.size):
    for k in range(phi.size):
        answer[i,k] = abs(arrayPattern(theta[i],phi[k],N,ycoords,zcoords,wavelength))
        coords[i*phi.size + k, 0] = math.cos(phi[k])*math.cos(theta[i])*abs(answer[i,k])
        coords[i*phi.size + k, 1] = math.sin(phi[k])*math.cos(theta[i])*abs(answer[i,k])
        coords[i*phi.size + k, 2] =- math.sin(theta[i])*abs(answer[i,k])

##for i in range(theta.size):
##    answer2[i] = arrayPattern(theta[i],0,N,ycoords,zcoords,wavelength)



fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(coords[:,0],coords[:,1], coords[:,2],'.')
plt.xlabel('x')
plt.ylabel('y')
##plt.zlabel('z')


##print(arrayPattern(0.2,0,N,ycoords,zcoords,wavelength))
plt.show()

