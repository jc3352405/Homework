from pylab import *
from math import factorial

#  1)
	
print 'Using Forward-Difference Algorithm'
def forward(f,x,h):
	a=[]
	b=[]
	while (h>10**-7):
		h/=2.0
		f2=(f(float32(x+h))-f(float32(x)))/float32(h)
		a.append(h)
		b.append(f2)
				
		
	return a,b
[h,f2] = forward(cos,0.1, 0.1)
error= abs((f2+sin(float32(0.1)))/-sin(float32(0.1)))
print 'step h is', h
print 'The first derivative by using Forward-Difference Algorithm is', f2
print 'The relative error is',error
loglog(h,error,basex=10)
grid(True)
xlabel('h')
ylabel('relative error')
title('log-log plot of derivatives of cos (0.1) using forward-difference algorithm')	
show()

[h,f2] = forward(cos,10, 1.0)
error= abs((f2+sin(float32(10.0)))/-sin(float32(10.0)))
print 'step h is', h
print 'The first derivative by using Forward-Difference Algorithm is', f2
print 'The relative error is',error
loglog(h,error,basex=10)
grid(True)
xlabel('h')
ylabel('relative error')
title('log-log plot of derivatives of cos (10.0) using forward-difference algorithm')	
show()

[h,f2] = forward(exp,0.1, 0.1)
error= abs((f2-exp(float32(0.1)))/exp(float32(0.1)))
print 'step h is', h
print 'The first derivative by using Forward-Difference Algorithm is', f2
print 'The relative error is',error
loglog(h,error,basex=10)
grid(True)
xlabel('h')
ylabel('relative error')
title('log-log plot of derivatives of exp (0.1) using forward-difference algorithm')	
show()

[h,f2] = forward(exp,10.0, 1.0)
error= abs((f2-exp(float32(10.0)))/exp(float32(10.0)))
print 'step h is', h
print 'The first derivative by using Forward-Difference Algorithm is', f2
print 'The relative error is',error
loglog(h,error,basex=10)
grid(True)
xlabel('h')
ylabel('relative error')
title('log-log plot of derivatives of exp (10.0) using forward-difference algorithm')	
show()

print 'Using Central-Difference Algorithm'
def central(f,x,h):
	a=[]
	b=[]
	while (h>10**-7):
		h/=2.0
		f2=(f(float32(x+h/2))-f(float32(x-h/2)))/float32(h)
		a.append(h)
		b.append(f2)
				
	return a,b

[h,f2] = central(cos,0.1, 0.1)
error= abs((f2+sin(float32(0.1)))/-sin(float32(0.1)))
print 'step h is', h
print 'The first derivative by using Central-Difference Algorithm is', f2
print 'The relative error is',error
loglog(h,error,basex=10)
grid(True)
xlabel('h')
ylabel('relative error')
title('log-log plot of derivatives of cos (0.1) using central-difference algorithm')	
show()

[h,f2] = central(cos,10.0, 1.0)
error= abs((f2+sin(float32(10.0)))/-sin(float32(10.0)))
print 'step h is', h
print 'The first derivative by using Central-Difference Algorithm is', f2
print 'The relative error is',error
loglog(h,error,basex=10)
grid(True)
xlabel('h')
ylabel('relative error')
title('log-log plot of derivatives of cos (10.0) using central-difference algorithm')	
show()

[h,f2] = central(exp,0.1, 0.1)
error= abs((f2-exp(float32(0.1)))/exp(float32(0.1)))
print 'step h is', h
print 'The first derivative by using Central-Difference Algorithm is', f2
print 'The relative error is',error
loglog(h,error,basex=10)
grid(True)
xlabel('h')
ylabel('relative error')
title('log-log plot of derivatives of exp (0.1) using central-difference algorithm')	
show()

[h,f2] = central(exp,10.0, 1.0)
error= abs((f2-exp(float32(10.0)))/exp(float32(10.0)))
print 'step h is', h
print 'The first derivative by using Central-Difference Algorithm is', f2
print 'The relative error is',error
loglog(h,error,basex=10)
grid(True)
xlabel('h')
ylabel('relative error')
title('log-log plot of derivatives of exp (10.0) using central-difference algorithm')	
show()

print 'Using Extrapolated-Difference Algorithm'
def extra(f,x,h):
	a=[]
	b=[]
	while (h>10**-7):
		f1=(f(float32(x+h/2))-f(float32(x-h/2)))/float32(h)
		f2=(f(float32(x+h/4))-f(float32(x-h/4)))/float32(h/2)
		h/=2.0
		f4=(4*f2-f1)/3
		a.append(h)
		b.append(f4)
		
		
	return a,b

[h,f2] = extra(cos,0.1, 1.0)
error= abs((f2+sin(float32(0.1)))/-sin(float32(0.1)))
print 'step h is', h
print 'The first derivative by using Extrapolated-Difference Algorithm is', f2
print 'The relative error is',error
loglog(h,error,basex=10)
grid(True)
xlabel('h')
ylabel('relative error')
title('log-log plot of derivatives of cos (0.1) using extrapolated-difference algorithm')	
show()

[h,f2] = extra(cos,10.0, 0.1)
error= abs((f2+sin(float32(10.0)))/-sin(float32(10.0)))
print 'step h is', h
print 'The first derivative by using Extrapolated-Difference Algorithm is', f2
print 'The relative error is',error
loglog(h,error,basex=10)
grid(True)
xlabel('h')
ylabel('relative error')
title('log-log plot of derivatives of cos (10.0) using extrapolated-difference algorithm')	
show()

[h,f2] = extra(exp,0.1, 1.0)
error= abs((f2-exp(float32(0.1)))/exp(float32(0.1)))
print 'step h is', h
print 'The first derivative by using Extrapolated-Difference Algorithm is', f2
print 'The relative error is',error
loglog(h,error,basex=10)
grid(True)
xlabel('h')
ylabel('relative error')
title('log-log plot of derivatives of exp (0.1) using extrapolated-difference algorithm')	
show()

[h,f2] = extra(exp,10.0, 0.1)
error= abs((f2-exp(float32(10.0)))/exp(float32(10.0)))
print 'step h is', h
print 'The first derivative by using Extrapolated-Difference Algorithm is', f2
print 'The relative error is',error
loglog(h,error,basex=10)
grid(True)
xlabel('h')
ylabel('relative error')
title('log-log plot of derivatives of exp (10.0) using extrapolated-difference algorithm')	
show()


#  2)

#    Trapezoid rule
 
def trapezoid(f,a,b,n):
	m=[]
	E=[]
	for j in range (1, n, 10):
		s = -(f(a) + f(b))/2.0
		for i in linspace(a,b,j+1):
			s += f(i)
		I= s * (b-a)/j
		error=abs((I-(e-1)/e)/((e-1)/e))
		m.append(j)
		E.append(error)
		print I

	return I,m,E
    
[Int1,N,error] = trapezoid(lambda t:exp(-t),(float32(0.0)), (float32(1.0)), 2000)
loglog(N,error,basex=10)
grid(True)
xlabel('N')
ylabel('relative error')
title('log-log plot of the integral of exponential using trapezoid rule')	
show()

#    Simpson rule

def simpson(f,a,b,n):
	m=[]
	E=[]
	for j in range (1, n, 5):
		s1 = 0.0
		s2 = 0.0
		x1 = a
		x2 = a
		h=(b - a)/j
		for i in range(1, j, 2):
			x1 += 2*h
			s1 += 4*f(x1)
		for i in range(2, j-1, 2):
			x2 += 2*h
			s2 += 2*f(x2)

		I = (s1 + s2 + exp(a) + exp(b))* h/3
		error=abs((I-(e-1)/e)/((e-1)/e))
		m.append(j)
		E.append(error)	
		print I
	return  I,m,E
 
def function(x): return x

[Int2,N,error] = simpson(lambda t:exp(-t),(float32(0.0)), (float32(1.0)), 200)
loglog(N,error,basex=10)
grid(True)
xlabel('N')
ylabel('relative error')
title('log-log plot of the integral of exponential using simpson rule')	
show()

def gaussNodes(m,tol=10e-9):
 
    def legendre(t,m):
	p=0
        p0 = 1.0; p1 = t
        for k in range(1,m):
            p = ((2.0*k + 1.0)*t*p1 - k*p0)/(1.0 + k )
            p0 = p1; p1 = p
        dp = m*(p0 - t*p1)/(1.0 - t**2)
        return p,dp
 
    A = zeros(m)   
    x = zeros(m)   
    nRoots = (m + 1)/2          # Number of non-neg. roots
    for i in range(nRoots):
        t = cos(pi*(i + 0.75)/(m + 0.5))  # Approx. root
        for j in range(30): 
            p,dp = legendre(t,m)          # Newton-Raphson
            dt = -p/dp; t = t + dt        # method         
            if abs(dt) < tol:
                x[i] = t; x[m-i-1] = -t
                A[i] = 2.0/(1.0 - t**2)/(dp**2) # Eq.(6.25)
                A[m-i-1] = A[i]
                break
    return x,A

def gaussQuad(f,a,b,N): 
    	c1 = (b + a)/2.0
    	c2 = (b - a)/2.0
	M=[]
	E=[]
	for y in range(1,N,10):
    		x,A = gaussNodes(y)
    		sum = 0.0
    		for i in range(len(x)):
        		sum = sum + A[i]*f(c1 + c2*x[i])
		I=c2*sum
    		error=abs((I-(e-1)/e)/((e-1)/e))
		M.append(y)
		E.append(error)	
		print I
		
	return I,M,E

[Int3,N,error] = gaussQuad(lambda t:exp(-t),(float32(0.0)), (float32(1.0)), 300)
loglog(N,error,basex=10)
grid(True)
xlabel('N')
ylabel('relative error')
title('log-log plot of the integral of exponential using Gaussian-Quadrature')	
show()


#  3)
#   a)

a=9301
c=49297
m=233280

def lcg():
        global xi        
        xi = (a*xi + c)%m
        return xi

def walk (nMax, N):
        global iseed
        trials =[]
       	
	for iseed in range(1,N+1):

		global xi
		xi = iseed
		x = range(nMax+1)
		x[0]= 0         #start at the origin

		for k in range(1, nMax+1):       
 
		        rnd = lcg()
		        if ((int(rnd*pi))%2 == 0):
		                x[k]=x[k-1] + 1
		        else:
		                x[k]=x[k-1] - 1
					              
		trials.append(x)   #save the result for each xi value	
		
	A=[]
	for i in range(1,nMax+1):
		B=[]
		for j in range(N):
			B.append(trials[j][i])
		A.append(B)
	return A

Xn=walk(500,1000)

def powerlist(list,p):
	K=[]
	for i in range(len(list)):
		K.append(list[i]**p)
	return K


def expect(m,n,p,list):
	A=[]
	for i in range (m):
		A.append(sum(powerlist(list[i],p))/float32(n)) 
	return A

Xn1=expect(500,1000.0,1,Xn)
Xn2=expect(500,1000.0,2,Xn)
Xn3=expect(500,1000.0,3,Xn)
Xn4=expect(500,1000.0,4,Xn)

Xn12=powerlist(Xn1,2)

M=[]
for i in range(500):
	M.append(Xn2[i]-Xn12[i])

plot(range(1,501),M)
xlabel('steps')
ylabel('<Xn^2>-<Xn>^2')
title('Plot of $\sigma$^2')
show()

Xn23=powerlist(Xn2,1.5)

N=[]
for i in range(500):
	N.append(Xn3[i]/Xn23[i])

plot(range(1,501),N)
xlabel('steps')
ylabel('<Xn^3>/$\sigma$^3')
title('Plot of s3')
show()

Xn22=powerlist(Xn2,2)

P=[]
for i in range(500):
	P.append(Xn4[i]/Xn22[i]-3)

plot(range(1,501),P)
xlabel('steps')
ylabel('<Xn^4>/$\sigma$^4-3')
title('Plot of s4')
show()
