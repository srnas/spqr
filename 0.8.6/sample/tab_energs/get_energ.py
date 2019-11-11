import numpy, sys
FILENAME=sys.argv[1]

xmin=-10
ymin=-10
zmin=-1.7
dx=0.25
dy=0.25
dz=0.2
NX=81
NY=81
NZ=18
xc=float(sys.argv[2])
yc=float(sys.argv[3])
zc=float(sys.argv[4])
shift=int(sys.argv[5])

NTOT=NX*NY*NZ*3
NLOC=NX*NY*NZ
ARRAY=[]

file=open(FILENAME,'r')
lin=file.readline()
lin=file.readline()
lin=file.readline()


for i in xrange(0,NTOT):
    lin=file.readline()
    #print lin.split()[0],i
    ARRAY.append(float(lin.split()[0]))
    xi=int((xc-xmin+0.5*dx)/dx)
    yi=int((yc-ymin+0.5*dy)/dy)
    zi=int((zc-zmin+0.5*dz)/dz)
    ind=zi+NZ*yi+NZ*NY*xi+shift*NLOC

print ARRAY[ind]
