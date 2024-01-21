import numpy
import matplotlib.pyplot

# 1

K = 100000
r = 0.4
t0 = 75
tf = 115
re = 101
x = numpy.linspace(t0,tf,re)
h = x[1]-x[0]

gom = [10]
for i in range(1,re):
    gom.append(gom[i-1]+h*r*gom[i-1]*numpy.log(K/gom[i-1]))

ver = [10]
for i in range(1,re):
    ver.append(ver[i-1]+h*r*ver[i-1]*(1-ver[i-1]/K))

matplotlib.pyplot.plot(x,gom)
matplotlib.pyplot.plot(x,ver)
matplotlib.pyplot.legend(["Gompertz","Verhulst"])
matplotlib.pyplot.xlabel("Time")
matplotlib.pyplot.ylabel("Number of cancer cells")
matplotlib.pyplot.title("Gompertz and Verhulst growth models")
matplotlib.pyplot.show()


# 2
N1 = [3]
N2 = [4]
e1,y1,h1 = (1.25,0.5,0.1)
e2,y2,h2 = (0.5,0.2,0.2)
re = 101
t0 = 0
tf = 30
x = numpy.linspace(t0,tf,re)
for i in range(1,re):
    N1.append(N1[i-1]+(e1-y1*(h1*N1[i-1]+h2*N2[i-1]))*N1[i-1]/10)
    N2.append(N2[i-1]+(e2-y2*(h1*N1[i-1]+h2*N2[i-1]))*N2[i-1]/10)

matplotlib.pyplot.plot(x,N1)
matplotlib.pyplot.plot(x,N2)
#matplotlib.pyplot.show()

N1 = [3]
N2 = [4]
e1,y1,h1 = (5,4,1)
e2,y2,h2 = (5,8,4)
re = 101
t0 = 0
tf = 30
x = numpy.linspace(t0,tf,re)
for i in range(1,re):
    N1.append(N1[i-1]+(e1-y1*(h1*N1[i-1]+h2*N2[i-1]))*N1[i-1]/1000)
    N2.append(N2[i-1]+(e2-y2*(h1*N1[i-1]+h2*N2[i-1]))*N2[i-1]/1000)

matplotlib.pyplot.plot(x,N1)
matplotlib.pyplot.plot(x,N2)
matplotlib.pyplot.legend(["N1a","N2a","N1b","N2b"])
matplotlib.pyplot.xlabel("Time")
matplotlib.pyplot.ylabel("Pop count")
matplotlib.pyplot.title("Competition model with different environmental values")
matplotlib.pyplot.show()

N1 = [4]
N2 = [8]
e1,g1,h1 = (0.8,1,0.3)
e2,g2,h2 = (0.4,0.5,0.4)
re = 1001
t0 = 0
tf = 30
x = numpy.linspace(t0,tf,re)
for i in range(1,re):
    N1.append(N1[i-1]+(e1-g1*(h1*N1[i-1]+h2*N2[i-1]))*N1[i-1]/10)
    N2.append(N2[i-1]+(e2-g2*(h1*N1[i-1]+h2*N2[i-1]))*N2[i-1]/10)

matplotlib.pyplot.plot(N1,N2)


N1 = [8]
N2 = [8]
e1,g1,h1 = (0.8,1,0.3)
e2,g2,h2 = (0.4,0.5,0.4)
re = 1001
t0 = 0
tf = 30
x = numpy.linspace(t0,tf,re)
for i in range(1,re):
    N1.append(N1[i-1]+(e1-g1*(h1*N1[i-1]+h2*N2[i-1]))*N1[i-1]/1000)
    N2.append(N2[i-1]+(e2-g2*(h1*N1[i-1]+h2*N2[i-1]))*N2[i-1]/1000)

matplotlib.pyplot.plot(N1,N2)


N1 = [12]
N2 = [8]
e1,g1,h1 = (0.8,1,0.3)
e2,g2,h2 = (0.4,0.5,0.4)
re = 1001
t0 = 0
tf = 30
x = numpy.linspace(t0,tf,re)
for i in range(1,re):
    N1.append(N1[i-1]+(e1-g1*(h1*N1[i-1]+h2*N2[i-1]))*N1[i-1]/1000)
    N2.append(N2[i-1]+(e2-g2*(h1*N1[i-1]+h2*N2[i-1]))*N2[i-1]/1000)

matplotlib.pyplot.plot(N1,N2)


x = numpy.linspace(0,13,20)
y = numpy.linspace(1,9,20)

X,Y = numpy.meshgrid(x,y)
dX = numpy.zeros(X.shape)
dY = numpy.zeros(Y.shape)

for i in range(X.shape[0]):
    for j in range(Y.shape[0]):
        dX[i,j] = (e1-g1*(h1*X[i,j]+h2*Y[i,j]))*X[i,j]
        dY[i,j] = (e2-g2*(h1*X[i,j]+h2*Y[i,j]))*Y[i,j]

matplotlib.pyplot.quiver(X,Y,dX,dY)




matplotlib.pyplot.legend(["N1,N2 = 4,8","N1,N2 = 8,8","N1,N2 = 12,8"])
matplotlib.pyplot.xlabel("N1 pop count")
matplotlib.pyplot.ylabel("N2 pop count")
matplotlib.pyplot.title("Phase graph for three different starting populations")

matplotlib.pyplot.show()
