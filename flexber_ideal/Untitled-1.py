from scipy.special import *
import math
from matplotlib import pyplot



a=ellipk(0.45)
print(a)


a =0
for i in range(100000):
    theta=((math.pi/2)/100000*i)
    a = a + ((1/((1 - 0.45*(math.sin(theta))**2)**0.5)) * (1/100000))
print(a)

a = []
theta = []
for i in range(1000):
    theta.append((math.pi/2)/1000*i)
    a.append((1/((1 - 0.45*(math.sin(theta[i]))**2)**0.5)) * (1/1000))


pyplot.plot(theta, a)
pyplot.show()


