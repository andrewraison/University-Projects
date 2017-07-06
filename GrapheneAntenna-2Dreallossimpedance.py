##Calculates the array pattern for a 1D array of cosine like elements
##
##theta = 0, and phi = 0 is the maximum of each element
##The first element is the LH most element
##All variables that need to be changed between runs are in the first section

import numpy as np
import math
import pdb
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from math import pi
import csv

############################################################################

##Variables of this problem
##uses 2D coords perpendicular to the theta=phi=0:  theta is between the yz plane and the x plane
##                                                  phi is in the yz plane measured from the z axis
desiredPhi = pi/4
desiredTheta = 0.2
dimensionsofarray = [2,1]   ##[y,z]
freq = 1.5*10**12           ##in hz
##NOTE wavelength and coords must be in mm
wavelength = 3*(10**11) / freq
y=0.3
z=wavelength/2

##############################################################################

hbar = 1.055e-34
k = 1.381e-23
T = 300
e = 1.6e-19
tau = 1e-12
Z0 = 377
A = pi*(hbar**2)/(2*(e**2)*k*T)
ycoords = np.empty((dimensionsofarray[0]*dimensionsofarray[1]))
zcoords = np.empty((dimensionsofarray[0]*dimensionsofarray[1]))
N = ycoords.size
L = np.empty((dimensionsofarray[0]*dimensionsofarray[1]))
mu = np.empty((dimensionsofarray[0]*dimensionsofarray[1]))


##################################################################################

i = -1
f = open('2eFF.csv', newline='')
reader = csv.reader(f, quoting=csv.QUOTE_NONE)
rowcount = sum(1 for row in f)
elFunc = np.empty((rowcount-1))
theta1 = np.empty((rowcount-1))
phi1 = np.empty((rowcount-1))
theta = np.empty((rowcount-1))
phi = np.empty((rowcount-1))
f = open('2eFF.csv', newline='')
reader = csv.reader(f, quoting=csv.QUOTE_NONE)

for row in reader:
    if i!=-1:
        phi1[i] = row[0]
        theta1[i] = row[1]
        elFunc[i] = row[2]
        phi[i] = phi1[i]*pi/180
        theta[i] = theta1[i]*pi/180
    i=i+1

################################################################################
##Main functions of the problem

##gives the element pattern for each individual element
def elementFunction(theta,phi,n):
    return elFunc[i]


##gives the losses for each individual element
def elementLoss(theta,phi,n,gamma):
    return 1#np.real(gamma)


##gives the array pattern
def arrayPattern(theta,phi,N,ycoords,zcoords,wavelength,gamma,i):
    distance = ( ycoords**2 + zcoords**2 )**0.5
    if theta == 0:
        perpdistance = distance
    else:
        a = np.arctan2(ycoords,zcoords)
        perpdistance = distance * np.cos(a - phi)
    
    elements = elementFunction(theta,phi,i) * elementLoss(theta,phi,0,gamma) * np.exp(1j*(2*pi*perpdistance*math.sin(theta)/wavelength + np.imag(gamma)))
    return np.sum(elements)


############################################################################
##calculate the variables of the problem

for i in range(dimensionsofarray[0]):
    for j in range(dimensionsofarray[1]):
        ycoords[dimensionsofarray[1]*i + j] = i*y -y/2
        zcoords[dimensionsofarray[1]*i + j] = j*z


def Z(mu):
##  assume that mu is sufficently large that cosh(mu/2*k*T) = exp(mu/2*k*T)/2 to good approx
    return A*2*k*T/mu * (1/tau + 1j*2*pi*freq)


##Finds a new that creates the desired phase difference, reset changes the calculation by 2pi
def findMu(mU, newL, L,z,theta,reset):
    if reset:
        g = abs((np.sin(theta)*z))%1 * np.sign((np.sin(theta)*z))
    else:
        g = abs(-1 + abs(np.sin(theta)*z)) * np.sign(-1 + abs(np.sin(theta)*z)) * np.sign((np.sin(theta)*z))
    print(g)
    return 2*pi*freq*2*k*T*newL*A/(L*np.imag(Z(mU)) + Z0 * wavelength *2*pi* g)


##Note that muLimit depends strongly on L/wavelength;
##due to limitations of elements mu cannot be too small or large
muLimit = 0.0
for i in range(dimensionsofarray[0]):
    for j in range(dimensionsofarray[1]):
        current = dimensionsofarray[1]*i + j
        L[current] = 100*wavelength

        ##set the initial mu, depending on whether mu will increase or decrease alongthe array
        if i==0 and j==0:
            if desiredTheta<0 or desiredTheta==0:
                mu[0] = 2*e
                muLimit = (4*pi*freq*k*T*A*L[0]/(Z0*wavelength))/(-2*pi+4*pi*freq*k*T*A*L[0]/(Z0*wavelength*mu[0]))
            else:
                mu[0] = 5*e
                muLimit = (4*pi*freq*k*T*A*L[0]/(Z0*wavelength))/(2*pi+4*pi*freq*k*T*A*L[0]/(Z0*wavelength*mu[0]))
        else:

            ##define which element the phase is being calculated wrt
            previous=0
            if j==0:
                previous = dimensionsofarray[1]*(i-1)
            else:
                previous = current-1

                
            ##Calculate per dist from last point
            distancelp = ( (ycoords[current]-ycoords[previous])**2 + (zcoords[current] - zcoords[previous])**2 )**0.5
            if desiredTheta == 0:
                perpdistancelp = distancelp
            else:
                a = np.arctan2(ycoords[current]-ycoords[previous],zcoords[current]-zcoords[previous])
                perpdistancelp = distancelp * np.cos(a - desiredPhi)

            ##Find the next mu
            mu[current] = findMu(mu[previous],L[current],L[previous],perpdistancelp/wavelength,desiredTheta,True)

            ##If it is too large or small, it will be shifted back the equivalent of 2pi phase
            if (muLimit>mu[0] and mu[current]>muLimit) or (muLimit<mu[0] and muLimit>mu[current]):
                mu[current] = findMu(mu[previous],L[current],L[previous],perpdistancelp/wavelength,desiredTheta,False)
                   

gamma = Z(mu)*L/(Z0*wavelength)
##print(gamma)


##############################################################################
##solve the problem


##theta = np.linspace(-1,1,30)
##phi = np.linspace(0,pi,30)

answer = np.empty((theta.size))
coords = np.empty((theta.size,3))
coordstest=np.empty((theta.size,3))

for i in range(theta.size):
    answer[i] = abs(arrayPattern(theta[i],phi[i],N,ycoords,zcoords,wavelength,gamma,i))
    coords[i, 2] = math.cos(phi[i])*math.sin(theta[i])*abs(answer[i])
    coords[i, 1] = math.sin(phi[i])*math.sin(theta[i])*abs(answer[i])
    coords[i, 0] = math.cos(theta[i])*abs(answer[i])
##    coordstest[i, 0] = math.cos(phi[i])*math.sin(theta[i])*elFunc[i]
##    coordstest[i, 1] = math.sin(phi[i])*math.sin(theta[i])*elFunc[i]
##    coordstest[i, 2] = math.cos(theta[i])*elFunc[i]



fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(coords[:,0],coords[:,1], coords[:,2])
plt.xlabel('x')
plt.ylabel('y')
plt.show()

##fig = plt.figure()
##ax = fig.gca(projection='3d')
##ax.plot(coordstest[:,0],coordstest[:,1], coordstest[:,2],'.')
##plt.xlabel('x')
##plt.ylabel('y')
##plt.show()

f = open("\u03B8=%f \u03C6working=%f.txt" %(desiredTheta,desiredPhi),"w+")
f.write("Desired angles:\n\u03B8=%f\t\u03C6=%f\n\n" %(desiredTheta,desiredPhi))
f.write("Element no\ty coord/mm\tz coord/mm\tmu/e\n")
for i in range(N):
    f.write("%d\t\t%f\t%f\t%f\n" %(i,ycoords[i],zcoords[i],mu[i]/e))
f.write("\nSpacing (per wavelength):\ty=%f\tz=%f"%(y/wavelength,z/wavelength))


##l = np.linspace(1,100,100)
##leng = wavelength*l
##plt.plot(l,(4*pi*freq*k*T*A*leng/(Z0*wavelength))/(2*pi+4*pi*freq*k*T*A*leng/(Z0*wavelength*mu[0]))/e)
##plt.show()
