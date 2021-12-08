# Finite difference scheme to solve the heat equation.
# Specifically, the function should implement an implicit scheme which is backward difference in time and central 
# difference in space.
# Ut = D(Uxx)
# U(0, t) = alpha, U(1, t) = beta
# U(x, 0) = f(x)


def heat_eq_matrix(dt, dx, N, f, D, alpha, beta):
    xvals     = np.linspace(0,1,int(1/dx)+1)
    tfinal    =N*dt
    tvals     = np.linspace(0,tfinal,N+1)
    xx        =len(xvals)
    tt        =len(tvals)
    A         =np.zeros((xx,xx))
    A[0,0]    = 1
    A[-1,-1]  = 1
    u         =np.zeros((tt,xx))
    u[0,:]    = f(xvals) #1st row and all columns
    u[:,0]    = alpha #all rows and 1st column
    u[:,-1]   = beta #all rows and last column
    for i in range(1,xx-1):
        A[i,i]   = 1+((2*D*dt)/(dx)**2) #lamda2
        A[i,i-1] = -(D*dt)/(dx)**2 #lamda1
        A[i,i+1] = -(D*dt)/(dx)**2 #lamda1
    for t in range(1,N+1):
        u[t,:]=LA.inv(A)@u[t-1,:]
    return u 