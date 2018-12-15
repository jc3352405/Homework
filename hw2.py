from pylab import *
from numpy import polyfit

def some():             #define parameter
        global G
        G=6.67*10**(-11)
        global h
        h=2.7*10**15
        global c
        c=3.0*10**8
        global M
        M=1.989*10**30
        
some()


def Runge_Kutta(g, a, b, h1, initial):  #Runge-Kutta Method
        z = []
        y = []
        state = initial
        h2 = h1/2
        for t in arange(a, b+h1, h1):
                z.append(state[0])
                y.append(state[1])
                k1 = g(t, state)
                k2 = g(t+h2, state+h2*k1)
                k3 = g(t+h2, state+h2*k2)
                k4 = g(t+h1, state+h1*k3)
                state += h1/6 * (k1 + 2*k2 + 2*k3 + k4)
        return z,y


def f1(t, state):               #Eq.1 with relativistic term
        [u, d_u] = state
        y = d_u
        d_y = -u + 1.0 + 3*((G*M/c/h)**2)*(u**2)
        return array([y, d_y])


def f2(t, state):               #Eq.1 without relativistic term
        [u, d_u] = state
        y = d_u
        d_y = -u + 1.0
        return array([y, d_y])


# 1)


[z2,d2] = Runge_Kutta(f2, 0, 50.0, 0.01, [1.0+0.205630, 0])     #solve Eq.1 without relativistic term with different timesteps
[z3,d3] = Runge_Kutta(f2, 0, 50.0, 0.05, [1.0+0.205630, 0])
[z4,d4] = Runge_Kutta(f2, 0, 50.0, 0.1, [1.0+0.205630, 0])

u1=[]
u2=[]   
u3=[]                                   #Exact Newtonian result with different timesteps
for i in arange(0, 50.01, 0.01):
        u1.append(1.0 + 0.205630*cos(i))
for i in arange(0, 50.05, 0.05):
        u2.append(1.0 + 0.205630*cos(i))
for i in arange(0, 50.1, 0.1):
        u3.append(1.0 + 0.205630*cos(i))

phi1 =  [i for i in arange(0, 50.01, 0.01)]
phi2 =  [i for i in arange(0, 50.05, 0.05)]
phi3 =  [i for i in arange(0, 50.1, 0.1)]

plot(phi1,u1,phi1,z2,phi2,u2,phi2,z3,phi3,u3,phi3,z4)   #plot of the experimental value and theoretical value
xlabel("$\phi$")
ylabel('u/(GM/h^2)')
title('u vs $\phi$ with different timesteps')
show()


error1=[]
error2=[]
error3=[]
for i in range(5001):
        error1.append((z2[i]-u1[i])/u1[i])
for i in range(1001):
        error2.append((z3[i]-u2[i])/u2[i])
for i in range(501):
        error3.append((z4[i]-u3[i])/u3[i])
plot(phi1,error1,phi2,error2,phi3,error3)               #plot of the relative error
xlabel("$\phi$")
ylabel('relative error ')
title('relative error with different timesteps')
show()  



[z1,d1] = Runge_Kutta(f1, 0, 15.0, 0.001, [1.0+0.205630, 0])    #solve Eq.1 with relativistic term
phi4 =  [i for i in arange(0, 15.001, 0.001)]
plot(phi4, z1)
xlabel("$\phi$")
ylabel('u/(GM/h^2)')
show()

def turn(d):
        a=[]                    #find the turning point
        for i in range(15000):
                if d[i]<=0 and d[i+1]>=0:
                        a.append(i)
        return a
a=turn(d1)
x1=[phi4[a[0]-1],phi4[a[0]],phi4[a[0]+1]]
x2=[phi4[a[1]-1],phi4[a[1]],phi4[a[1]+1]]
y1=[z1[a[0]-1],z1[a[0]],z1[a[0]+1]]
y2=[z1[a[1]-1],z1[a[1]],z1[a[1]+1]]

p1=polyfit(x1,y1,2)             #fit a quadratic
fit1=poly1d(p1)
p2=polyfit(x2,y2,2)
fit2=poly1d(p2)
xfit=linspace(0,15,30)
plot(x1,y1,'.',xfit,fit1(xfit),'-',x2,y2,'.',xfit,fit2(xfit),'-')
xlabel("$\phi$")
ylabel('u/(GM/h^2)')
title('Quadratic fitting')
show()

def minimum(a):                 #find the minimum
        t=-a[1]/2/a[0]
        return t

delta_phi=minimum(p2)-minimum(p1)       #the difference in angle between one perihelion and the next
print 'The angle shift per period is',delta_phi-2*pi
print 'The angle shift per century is',(delta_phi-2*pi)*100*365/88*180/pi*3600,'"'


# 2)

def f3(t, state):       #New equation in a form analogous to Eq.1
        [u, d_u] = state
        y = d_u
        d_y = -u + u**0.05 + 3*((G*M/c/h)**2)*(u**2)
        return array([y, d_y])

[z5,d5] = Runge_Kutta(f3, 0, 15.0, 0.001, [1.0+0.205630, 0])
plot(phi4, z5)
xlabel("$\phi$")
ylabel('u/(GM/h^2)')
show()

b=turn(d5)
x3=[phi4[b[0]-1],phi4[b[0]],phi4[b[0]+1]]
x4=[phi4[b[1]-1],phi4[b[1]],phi4[b[1]+1]]
y3=[z5[b[0]-1],z5[b[0]],z5[b[0]+1]]
y4=[z5[b[1]-1],z5[b[1]],z5[b[1]+1]]

p3=polyfit(x3,y3,2)             #fit a quadratic
fit3=poly1d(p3)
p4=polyfit(x4,y4,2)
fit4=poly1d(p4)
plot(x3,y3,'.',xfit,fit3(xfit),'-',x4,y4,'.',xfit,fit4(xfit),'-')
xlabel("$\phi$")
ylabel('u/(GM/h^2)')
title('Quadratic fitting')
show()

delta_phi2=minimum(p4)-minimum(p3)      #the difference in angle between one perihelion and the next
print 'The angle shift under new gravitatioinal Newtonian force per period is',delta_phi2-2*pi
print 'The angle shift under new gravitatioinal Newtonian force per century is',(delta_phi2-2*pi)*100*365/88*180/pi,'degree'

# 3)

def Runge_Kutta2(g, a, b, h1, initial): #Runge-Kutta Method for stellar structure equation
        z = []
        y = []
        state = initial
        h2 = h1/2
        for t in arange(a, b+h1, h1):
                if state[0]>=0:
                        z.append(state[0])
                        y.append(state[1])
                        k1 = g(t, state)
                        k2 = g(t+h2, state+h2*k1)
                        k3 = g(t+h2, state+h2*k2)
                        k4 = g(t+h1, state+h1*k3)
                        state += h1/6 * (k1 + 2*k2 + 2*k3 + k4)
                        
        return z,y
#3)

k=3.0006*10**(-3)
a=G*M/c**2/1000

def f4(r, state):       #Rewrite Newtonian equation
        [p, m] = state
        ep=2.4216*(p**0.6)+2.8663*p
        d_p = -a*m*ep/r**2
        d_m = 4*k*pi*ep*r**2
        return array([d_p, d_m])

[P1,M1] = Runge_Kutta2(f4, 0.001, 15.0, 0.001, [0.01, 0])
r1 =  [i for i in arange(0.001, 15.001, 0.001)]

plot(r1, P1)
xlabel("radius")
ylabel('pressure')
show()

plot(r1, M1)
xlabel("radius")
ylabel('mass')
show()


def f5(r, state):       #Rewrite TOV equation
        [p, m] = state
        ep=2.4216*(p**0.6)+2.8663*p
        d_p = (-a*m*ep/r**2)*(1+p/ep)*(1+4*pi*k*r*p/m)/(1-2*a*m/r)
        d_m = 4*k*pi*ep*r**2
        return array([d_p, d_m])

[P2,M2] = Runge_Kutta2(f5, 0.001, 13.0, 0.001, [0.01, 0.00001])
r2 =  [i for i in arange(0.001, 13.001, 0.001)]
plot(r2, P2)
xlabel("radius")
ylabel('pressure')
show()

plot(r2, M2)
xlabel("radius")
ylabel('mass')
show()

R=[]
M=[]
for p in arange(0.01, 0.05, 0.001):
        print p
        [P3,M3] = Runge_Kutta2(f5, 0.001, 13.0, 0.001, [p, 0.00001])
        s=min(P3)
        R1=r2[P3.index(s)]
        M4=M3[P3.index(s)]
        R.append(R1)
        M.append(M4)

plot(R, M)
xlabel("radius")
ylabel('mass')
show()
