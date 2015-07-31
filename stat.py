import numpy as np
from matplotlib import pyplot as plt
import math

d = 100000
term = 2000
data18 = np.loadtxt("output1.8")
data19 = np.loadtxt("output1.9")
data20 = np.loadtxt("output2.0")
data21 = np.loadtxt("output2.1")
data22 = np.loadtxt("output2.2")
data23 = np.loadtxt("output2.3")

data18 = data18[term:d]
data19 = data19[term:d]
data20 = data20[term:d]
data21 = data21[term:d]
data22 = data22[term:d]
data23 = data23[term:d]

#print data
#data = data[2000:100000]
#plt.hist( data, bins = 50)
#plt.ylabel("Frequency")
#plt.show()

m18 = np.mean(data18)
m19 = np.mean(data19)
m20 = np.mean(data20)
m21 = np.mean(data21)
m22 = np.mean(data22)
m23 = np.mean(data23)

print "T=1.8: " + str(np.mean(data18)) + str(np.std(data18))    
print "T=1.9: " + str(np.mean(data19)) + str(np.std(data19))
print "T=2.0: " + str(np.mean(data20)) + str(np.std(data20))
print "T=2.1: " + str(np.mean(data21)) + str(np.std(data21))
print "T=2.2: " + str(np.mean(data22)) + str(np.std(data22))
print "T=2.3: " + str(np.mean(data23)) + str(np.std(data23))


"""
binVariance = np.array()
for i in range(20):
    for j in range(len(data)):
        binVariance = np.append(binVariance, np.std(data[j*i,]))
        """

x = [1.8, 1.9, 2.0, 2.1, 2.2, 2.3]

a, b =  np.polyfit(x, [m18,m19,m20,m21,m22,m23],1)
y = b + np.multiply(a,x)
print "Specific heat (derivative): " + str(a)
data20sq = data20**2
meanSq = np.mean(data20sq)
sqMean = np.mean(data20)**2
specHeat = (meanSq-sqMean)/(4)
print "Specific heat (variance): " + str(specHeat)

plt.scatter(x, [m18,m19,m20,m21,m22,m23])
plt.plot([1.8, 1.9, 2.0, 2.1, 2.2, 2.3], y)
#plt.text(2, -293, 'Cv',)
plt.xlabel("Temperature")
plt.ylabel("U")
plt.show()