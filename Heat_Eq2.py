# finite difference scheme to solve the heat equation with Neumann boundary conditions. 
# Specifically, the function should implement an explicit scheme which is forward difference in time and central difference
# in space.
# Ut = D(Uxx)
# Ux(0, t) = alpha, Ux(1, t) = beta
# U(x, 0) = f(x)

def heat_eq_matrix(dt, dx, N, f, D, alpha, beta):
    tfinal    = N*dt
    xvals     = np.linspace(0,1,int(1/dx)+1)
    tvals     = np.linspace(0,tfinal,N+1)

    tt        = len(tvals)
    xx        = len(xvals)

    u         = np.zeros((tt,xx))
    u[0,:]    = f(xvals) #1st row and all columns



    l1         =(D*dt)/(dx**2)
    l2         =1-((2*D*dt)/(dx**2))
    A          =np.zeros((xx,xx))
    A[0,0]     =l2
    A[0,1]     =2*l1
    A[-1,-1] =l2
    A[-1,-2] =2*l1
   
    c          =np.zeros((xx))
    c[0]       =-2*l1*dx*alpha
    c[-1]      =2*l1*dx*beta
    #print(c)

    for i in range(1,xx-1):
        A[i,i-1]  =l1
        A[i,i]    =l2
        A[i,i+1]  =l1
            
    for t in range(1,N+1):
        u[t,:]=A@u[t-1,:]+c
    return u 