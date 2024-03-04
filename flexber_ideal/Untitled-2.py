from scipy.special import *
import math
from matplotlib import pyplot



a=ellipk(0.45)
#print(a)

a = 0
theta = 0
a1 = 0
a2 = math.pi/2
n = 10000000
d_theta = (a2 - a1) / n

for i in range(n):
    a = a + ((((1 - 0.45*(math.sin(theta))**2)**(-0.5)))*d_theta)
    theta = theta + d_theta

#print(a)

def F(a1,a2,p):
    n = 10000000
    d_theta = (a2 - a1) / n
    theta = 0
    a = 0
    
    for i in range(n):
        a = a + ((((1 - (p**2)*(math.sin(theta))**2)**(-0.5)))*d_theta)
        theta = theta + d_theta
    return a

print(F(0,math.pi/2,0.45**0.5))
x = F(0,math.pi/2,0.45**0.5)
print(x)



