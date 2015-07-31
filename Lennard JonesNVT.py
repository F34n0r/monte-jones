import numpy as np
import math

#all numbers are in Lennard Jones units
n = 100.0                   #number of particles
rho = 0.5                   #particles density
l = (n/rho)** (1/3.0)       #box length
rc = l/2                    #cut-off radius
temperature = 2.0       

k = 10000                   #number of steps

acceptance = 0
particles = 0
dMove = l/20.0              #size of the move step

energyList = list()

def distributeParticles():
    "Distributes the particles randomly"
    
    
    global n, l, particles
    particles = l * np.random.rand(n,3)

 
def lattice():
    global n, l
    #for i in range(int(n)):
        
        
   
def potential(r2):
    "The Lennard-Jones potential. Distance squared is the parameter"
    #print "The distance is: " + str(math.sqrt(r2)) + ", and the potential is: " + str((1/r2)**6.0 - (1/r2)**3.0)
    return (1/r2)**6.0 - (1/r2)**3.0


def eTotal():
    "Calculates the total energy of the system directly. Expensive"
    global particles, rc, energy
    en = 0
    for i in range(int(n)):
        for j in range(i):
            r2 = abs((particles[i,0:3] - particles[j,0:3]))
            r2[0] = min(r2[0], l-r2[0])
            r2[1] = min(r2[1], l-r2[1])
            r2[2] = min(r2[2], l-r2[2])
            r2 = np.sum(r2**2)
            #print "Distance between " +str(i) + " and " + str(j) + " is " + str(math.sqrt(r))
            if r2 <= rc*rc:
                en += ((1/r2)**6) - ((1/r2)**3)
                #print en
            #else:
                #print str(j) + " is too far from " + str(i)
        #print "Energy of the " +str(i) + " particle is: " + str(en) 
    return 4*en        
   
   
   
def eChange(j,x,y,z):
    global n, particles, energy, rc
    tempE = 0
    prevE = 0
    for i in range(int(n)):
        r2 = abs((particles[i,0:3] - [x,y,z]))
        r2[0] = min(r2[0], l-r2[0])
        r2[1] = min(r2[1], l-r2[1])
        r2[2] = min(r2[2], l-r2[2])
        r2 = np.sum(r2**2)
        if (r2 <= rc*rc) and not(i == j):
            tempE += potential(r2)
                
    for i in range(int(n)): 
        r2 = abs((particles[i,0:3] - particles[j,0:3]))
        r2[0] = min(r2[0], l-r2[0])
        r2[1] = min(r2[1], l-r2[1])
        r2[2] = min(r2[2], l-r2[2])
        r2 = np.sum(r2**2)
        if (r2 <= rc*rc) and not(i == j):
                prevE += potential(r2)
                
    rChange = 4 *( tempE - prevE)
    return rChange
    
def moveParticle(pChoice):
    global n, particles, dMove, l, energy, acceptance, temperature
    #print pChoice
    #print dr
    #print particles[pChoice,1:4]
    [x,y,z] = particles[pChoice, 0:3] + dMove*(np.random.rand(1,3)[0]-[0.5,0.5,0.5]) 
    
    if x > l:
        x -= l
    elif x < 0:
        x += l
    if y > l:
        y -= l
    elif y < 0:
        y += l
    if z > l:
        z -= l
    elif z < 0:
        z += l
        
    dE = eChange(pChoice,x,y,z)
    #print "Energy change: " + str(dE)
    if dE < 0:
        #print "Energy is smaller than zero: " + str(dE)
        #print "Old coordinates: " + str(particles[pChoice,1:4])
        particles[pChoice,0:3] = [x,y,z]
        #print "New coordinates: " + str(particles[pChoice,1:4])
        acceptance += 1
        energy += dE
    else:
        choice = np.random.random()
        #print "Chance is: " + str(math.exp(-dE/temperature))
        if choice < math.exp(-dE/temperature):
            #print "Old coordinates: " + str(particles[pChoice,1:4])
            particles[pChoice,0:3] = [x,y,z]
            #print "New coordinates: " + str(particles[pChoice,1:4])
            acceptance += 1
            energy +=dE

         
   
def simulate(steps):
    global n, particles, dMove, energy
    for i in range(steps):
        pChoice = np.random.randint(n)
        moveParticle(pChoice)
        energyList.append(energy)
        print "New total energy is: " + str(energy)
    #print "Energy change is: " + str(change)
    print "New energy is: " + str(eTotal())    
    
distributeParticles()
energy = eTotal()
startEnergy = energy


print "Box size is: " + str(l)
print "Total energy is: " + str(energy)

for i in range(1):
    simulate(k)
    acceptance = acceptance/float(k)

    np.savetxt("nout"+str(temperature)+"try"+str(i),energyList)