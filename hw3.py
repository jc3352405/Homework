from pylab import *
import numpy as np
from math import *
from math import acos
from mpl_toolkits.mplot3d import axes3d, Axes3D


def powerlist(list,p):			# for powering the list
        K=[]
        for i in range(len(list)):
                K.append(list[i]**p)
        return K

def avgrms(list,n):			# for calculate average and rms
        avg=sum(list)/float(n)
        rms=sqrt(sum(powerlist(list,2))/float(n) -(avg)**2)
        return avg,rms

def fft1d(list,n,d,):			# 1d FFT and shift
        fftresult=np.fft.fft(list)
        power=np.abs(fftresult)**2
        fftfreqs=np.fft.fftfreq(n,d)
        x=np.argsort(fftfreqs)
        y=[]
        y2=[]
        for i in x:
                y.append(power[i])
                y2.append(fftresult[i])
                i=i+1
        return fftresult, power, fftfreqs, x, y, y2

def readdata(list,n1,n2,a,b):		# read the Gaussian data and coarse the grid
        phi=np.zeros(n1).reshape((a, a))
        phi2=np.zeros(n2).reshape((b, b))
        for i in range(a):
                phi[i]=list[i*a:i*a+a]
                if i < a/2:
                        phi2[i]=list[i*a:i*a+a:2]
        fft_phi=np.fft.fft2(phi)
        fft_phi2=np.fft.fft2(phi2)
        power_phi=np.abs(fft_phi)**2
        power_phi2=np.abs(fft_phi2)**2
        return power_phi,power_phi2


def circlesum(list,n1):			# for circle integral to get the log-log plot
        sump=range(n1)
        a=range(n1)
        for h in range(n1):
                for i in range(0,n1,1):
                        for j in range(0,n1,1):
                                if (i**2+j**2)**0.5 > float(h) and (i**2+j**2)**0.5 <= float(h+1):
                                        sump[h]+=(list[i+n1][j+n1])
                                        a[h]=sump[h]/(pi/4*((h+1)**2-h**2))
        return a

def CIC1d(list,t,s,d):
        disnear=range(t)
        disfar=range(t)
        x=range(t)
        for j in range(t):
                disnear[j]=[]
                disfar[j]=[]
                x[j]=[]
                for i in range(s):
                        p=int(list[j][i]/d)
                        r1=abs(p*d-list[j][i])
                        r2=abs((p+1)*d-list[j][i])
                        disnear[j].append(r1)
                        disfar[j].append(r2)
                        x[j].append(p)
        return disnear, disfar,x


x1,y1,p1=list(loadtxt('gaussian1.dat', skiprows=1, unpack=True))		#read in the files
x2,y2,p2=list(loadtxt('gaussian2.dat', skiprows=1, unpack=True))
x3,y3,p3=list(loadtxt('gaussian3.dat', skiprows=1, unpack=True))
temperature, humidity = list(loadtxt('TandH.dat', skiprows=1, unpack=True))

## problem 1
N=104

first_temp=temperature[0:N/2]; second_temp=temperature[N/2:N]
first_humi=humidity[0:N/2]; second_humi=humidity[N/2:N]

a1temp,r1temp=avgrms(first_temp,N/2)		#calculate the average and rms
a2temp,r2temp=avgrms(second_temp,N/2)
atemp,rtemp=avgrms(temperature,N)
a1humi,r1humi=avgrms(first_humi,N/2)
a2humi,r2humi=avgrms(second_humi,N/2)
ahumi,rhumi=avgrms(humidity,N)

print 'The average temperature of the first year is', a1temp, 'and the rms is', r1temp		#print the result
print 'The average temperature of the second year is',a2temp, 'and the rms is', r2temp
print 'The average humidity of the first year is',a1humi, 'and the rms is', r1humi
print 'The average humidity of the second year is',a2humi, 'and the rms is', r2humi
print 'The average temperature of two years is',atemp, 'and the rms is', rtemp
print 'The average humidity of two years is',ahumi, 'and the rms is', rhumi

temp=temperature-atemp			# minus average to get clear result
humi=humidity-ahumi

fftt,pst,freqs,xt,yt,yt2=fft1d(temp,N,1.0)		# FFT for temperature and humidity
ffth,psh,freqs,xt,yh,yh2=fft1d(humi,N,1.0)
     
plot(freqs[xt], yt, label="Temperature")      		# plot the power spectrum
plot(freqs[xt], yh, label="Humidity")
xlabel("Frequency(1/week)")
ylabel('Power')
legend(loc="upper left")
title('Power Spectrum of Temperature and Humidity')
plt.savefig("fig1.jpg")
show()

yt2[50]=0; yt2[54]=0			# take away the seasonal variation
yh2[50]=0; yh2[54]=0

temp2=np.fft.ifft(yt2)+atemp		# inverse FFT after taking away the seasonal variation
humi2=np.fft.ifft(yh2)+ahumi

temp3=[]
humi3=[]
for i in range(N):
        temp3.append(temp2[i].real)		# only get the real part
        humi3.append(humi2[i].real)
        i=i+1
        
plot(freqs[xt], np.abs(yt2)**2, label="Temperature")      	#plot the result
plot(freqs[xt], np.abs(yh2)**2, label="Humidity")
xlabel("Frequency(1/week)")
ylabel('Power')
legend(loc="upper left")
title('Power Spectrum of Temperature and Humidity without seasonal variation')
plt.savefig("fig2.jpg")
show()

t=[i for i in range(N)]
plot(t,temperature,label="with variation")
plot(t,temp3,label="without variation")
xlabel("Weeks")
ylabel('Temperature')
legend(loc="lower right")
title('Temperature value w/o Seasonal Variation')
plt.savefig("fig3.jpg")
show()

plot(t,humidity,label="with variation")
plot(t,humi3,label="without variation")
xlabel("Weeks")
ylabel('Humidity')
legend(loc="lower right")
title('Humidity value w/o Seasonal Variation')
plt.savefig("fig4.jpg")
show()

atemp2,rtemp2=avgrms(temp3,N)		# compare the fluctuation between temperature and humidity
ahumi2,rhumi2=avgrms(humi3,N)

print 'The average temperature without variation is',atemp2, 'and the rms is', rtemp2
print 'The average humidity without variation is',ahumi2, 'and the rms is', rhumi2
print 'Comparing normalized rms of temperature and humidity', (rtemp2/atemp2),'and',(rhumi2/ahumi2), 'so humidity has larger fluctuation.'

## Problem 2

n=100
power_phi1,power_phi4=readdata(p1,n**2,(n/2)**2,n,n/2)		# read in the data
power_phi2,power_phi5=readdata(p2,n**2,(n/2)**2,n,n/2)
power_phi3,power_phi6=readdata(p3,n**2,(n/2)**2,n,n/2)

freqphi=np.fft.fftfreq(n,d=0.01)		# FFT and shift
freqphi2=np.fft.fftfreq(n/2,d=0.02)
x=np.argsort(freqphi)
x2=np.argsort(freqphi2)
kx = np.linspace(-n/2, n/2, n)
kx2 = np.linspace(-n/4, n/4, n/2)
X, Y = np.meshgrid(kx, kx)
X2, Y2 = np.meshgrid(kx2, kx2)

orderpowerphi1=np.fft.fftshift(power_phi1)
orderpowerphi2=np.fft.fftshift(power_phi2)
orderpowerphi3=np.fft.fftshift(power_phi3)
orderpowerphi4=np.fft.fftshift(power_phi4)
orderpowerphi5=np.fft.fftshift(power_phi5)
orderpowerphi6=np.fft.fftshift(power_phi6)

fig = plt.figure(1)		# plot the power spectrum of Gaussian data
ax = Axes3D(fig)
ax.plot_surface(X,Y, orderpowerphi1)
ax.set_xlabel('kx (kf)')
ax.set_ylabel('ky (kf)')
ax.set_zlabel('Power')
ax.set_title('Gaussian data 1')
plt.savefig("fig5.jpg")
show()

fig = plt.figure(2)
ax = Axes3D(fig)
ax.plot_surface(X,Y, orderpowerphi2)
ax.set_xlabel('kx (kf)')
ax.set_ylabel('ky (kf)')
ax.set_zlabel('Power')
ax.set_title('Gaussian data 2')
plt.savefig("fig6.jpg")
show()

fig = plt.figure(3)
ax = Axes3D(fig)
ax.plot_surface(X,Y, orderpowerphi3)
ax.set_xlabel('kx (kf)')
ax.set_ylabel('ky (kf)')
ax.set_zlabel('Power')
ax.set_title('Gaussian data 3')
plt.savefig("fig7.jpg")
show()

sump1=circlesum(orderpowerphi1,n/2)		# circle integral
sump2=circlesum(orderpowerphi2,n/2)            
sump3=circlesum(orderpowerphi3,n/2)

kn=[ i for i in arange(1.0,float(n/2+1))]
logkn=[math.log(y) for y in kn]
logsump1=[math.log(y) for y in sump1]
logsump2=[math.log(y) for y in sump2]
logsump3=[math.log(y) for y in sump3]

plot(logkn,logsump1,label="Gaussian data 1")		# log-log plot
plot(logkn,logsump3,label="Gaussian data 3")
xlabel("log k (kf)")
ylabel('log Power')
legend(loc="upper left")
title('Log-Log plot of Power Spectrum')
plt.savefig("fig8.jpg")
show()

plot(logkn,logsump2,label="Gaussian data 2")
xlabel("log k (kf)")
ylabel('log Power')
legend(loc="upper left")
title('Log-Log plot of Power Spectrum')
plt.savefig("fig9.jpg")
show()

fig = plt.figure(4)			# plot the coarsening data
ax = Axes3D(fig)
ax.plot_surface(X2,Y2, orderpowerphi4)
ax.set_xlabel('kx (kf)')
ax.set_ylabel('ky (kf)')
ax.set_zlabel('Power')
ax.set_title('Gaussian data 1 Coarse')
plt.savefig("fig10.jpg")
show()

fig = plt.figure(5)
ax = Axes3D(fig)
ax.plot_surface(X2,Y2, orderpowerphi5)
ax.set_xlabel('kx (kf)')
ax.set_ylabel('ky (kf)')
ax.set_zlabel('Power')
ax.set_title('Gaussian data 2 Coarse')
plt.savefig("fig11.jpg")
show()

fig = plt.figure(6)
ax = Axes3D(fig)
ax.plot_surface(X2,Y2, orderpowerphi6)
ax.set_xlabel('kx (kf)')
ax.set_ylabel('ky (kf)')
ax.set_zlabel('Power')
ax.set_title('Gaussian data 3 Coarse')
plt.savefig("fig12.jpg")
show()
	
sump4=circlesum(orderpowerphi4,n/4)		# circle integral for coarsening data
sump5=circlesum(orderpowerphi5,n/4)            
sump6=circlesum(orderpowerphi6,n/4)

kn2=[ i for i in arange(1.0,float(n/4+1))]
logkn2=[math.log(y) for y in kn2]
logsump4=[math.log(y) for y in sump4]
logsump5=[math.log(y) for y in sump5]
logsump6=[math.log(y) for y in sump6]

plot(logkn2,logsump4,label="Gaussian data 1 Coarse")	#log-log plot
plot(logkn2,logsump6,label="Gaussian data 3 Coarse")
xlabel("log k (kf)")
ylabel('log Power')
legend(loc="lower left")
title('Log-Log plot of Power Spectrum')
plt.savefig("fig13.jpg")
show()

plot(logkn2,logsump5,label="Gaussian data 2 Coarse")
xlabel("log k (kf)")
ylabel('log Power')
legend(loc="lower left")
title('Log-Log plot of Power Spectrum')
plt.savefig("fig14.jpg")
show()

## Problem 3
Lc=250
step=99
alpha=1.5
r0=1
times=500
ngrid=401
gridsize=float(Lc)/float((ngrid-1))
kgrid=(ngrid+1)/2

xnext=np.zeros(times*(step+1)).reshape((times,step+1))
ynext=np.zeros(times*(step+1)).reshape((times,step+1))
znext=np.zeros(times*(step+1)).reshape((times,step+1))

for i in range(times):
        initial=[np.random.uniform(0,Lc),np.random.uniform(0,Lc),np.random.uniform(0,Lc)]       # create initial point
        particle=np.random.uniform(0,1,(step,3))        #create following steps
        r=r0*particle[:,0]**(-1/alpha)
        theta=[]
        for j in range(step):
                theta.append(acos(1-2*particle[j,1]))
        phi=particle[:,2]*2*pi

        xnext[i][0]=initial[0]  #read the position as list
        ynext[i][0]=initial[1]
        znext[i][0]=initial[2]
        
        for k in range(step):           # periodic condition
                xnext[i][k+1]=(xnext[i][k]+r[k]*sin(theta[k])*cos(phi[k]))%Lc
                ynext[i][k+1]=(ynext[i][k]+r[k]*sin(theta[k])*sin(phi[k]))%Lc
                znext[i][k+1]=(znext[i][k]+r[k]*cos(theta[k]))%Lc


fig = plt.figure(7)             # plot two cluster
ax = Axes3D(fig)
ax.plot(xnext[0], ynext[0], znext[0], label='particle 1')
ax.plot(xnext[1], ynext[1], znext[1], label='particle 2')
ax.legend()
ax.set_xlim3d(0, Lc)
ax.set_ylim3d(0, Lc)
ax.set_zlim3d(0, Lc)
ax.set_title('Random Walk')
plt.savefig("fig15.jpg")
show()


rx1,rx2,px=CIC1d(xnext,times,step+1,gridsize)           # Using CIC interpolation to get particle density
ry1,ry2,py=CIC1d(ynext,times,step+1,gridsize)
rz1,rz2,pz=CIC1d(znext,times,step+1,gridsize)


gridbox=np.zeros(ngrid**3).reshape((ngrid,ngrid,ngrid))
for j in range(times):
        for i in range(step+1):
                px0=px[j][i]
                py0=py[j][i]
                pz0=pz[j][i]
                
                gridbox[px0][py0][pz0]=(gridsize-rx1[j][i])*(gridsize-ry1[j][i])*(gridsize-rz1[j][i])/(gridsize**3)
                gridbox[px0+1][py0][pz0]=(gridsize-rx2[j][i])*(gridsize-ry1[j][i])*(gridsize-rz1[j][i])/(gridsize**3)
                gridbox[px0][py0+1][pz0]=(gridsize-rx1[j][i])*(gridsize-ry2[j][i])*(gridsize-rz1[j][i])/(gridsize**3)
                gridbox[px0+1][py0+1][pz0]=(gridsize-rx2[j][i])*(gridsize-ry2[j][i])*(gridsize-rz1[j][i])/(gridsize**3)
                gridbox[px0][py0][pz0+1]=(gridsize-rx1[j][i])*(gridsize-ry1[j][i])*(gridsize-rz2[j][i])/(gridsize**3)
                gridbox[px0][py0+1][pz0+1]=(gridsize-rx1[j][i])*(gridsize-ry2[j][i])*(gridsize-rz2[j][i])/(gridsize**3)
                gridbox[px0+1][py0][pz0+1]=(gridsize-rx2[j][i])*(gridsize-ry1[j][i])*(gridsize-rz2[j][i])/(gridsize**3)
                gridbox[px0+1][py0+1][pz0+1]=(gridsize-rx2[j][i])*(gridsize-ry2[j][i])*(gridsize-rz2[j][i])/(gridsize**3)


p=range(ngrid)                  # plot particle density vs distance
q=range(ngrid)
for i in range(ngrid):
        p[i]=[]
        for j in range(ngrid):
                p[i].append(sum(gridbox[i][j]))
for i in range(ngrid):
        q[i]=sum(p[i])
x=range(ngrid)
plot(x,q)
xlabel("distance")
ylabel('Particle density')
title('Particle density vs. distance')
plt.savefig("fig16.jpg")
show()
        
print 'starting fftn'           # fftn
FFT3= np.fft.fftn(gridbox)
SFFT3=np.fft.fftshift(FFT3)
PFFT3=np.abs(SFFT3)**2

print 'starting sphere integral'                # sphere integral for log log plot
n1=kgrid-1
n2=ngrid
sumPFFT=range(n1)
sumPFFT3=range(n1)
csumPFFT=range(n1)
csumPFFT3=range(n1)
for i in range(0,n1,1):
        for j in range(0,n1,1):
                for k in range(0,n1,1):
                        h=int((i**2+j**2+k**2)**0.5)
                        if h<=n1-1:
                                sumPFFT[h]+=(PFFT3[i+n1][j+n1][k+n1])
                                csumPFFT[h]+=(PFFT3[i+n1][j+n1][k+n1])/(sinc(float(i)/float(n2))**4*sinc(float(j)/float(n2))**4*sinc(float(k)/float(n2))**4)
for m in range(n1):
        sumPFFT3[m]=sumPFFT[m]/(pi/6*((m+1)**3-m**3))
        csumPFFT3[m]=csumPFFT[m]/(pi/6*((m+1)**3-m**3))


kn3=[ i for i in arange(1.0,float(kgrid))]
logkn3=[math.log(y) for y in kn3]
logsumPFFT3=[math.log(y) for y in sumPFFT3]
logcsumPFFT3=[math.log(y) for y in csumPFFT3]

plot(logkn3,logsumPFFT3,label="Before window correction")       # log log plot of powerspectrum and wavenumber
plot(logkn3,logcsumPFFT3,label="After window correction")
xlabel("log wavenumber")
ylabel('log Power')
legend(loc="lower left")
title('Power Spectrum for Rayleigh-Levy random walk')
plt.savefig("fig17.jpg")
show()

