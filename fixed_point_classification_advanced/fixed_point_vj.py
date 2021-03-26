import matplotlib.pyplot as plt
import numpy as np

def partial_derivative_3p(x1,x2,function,h=0.001):
    return (function(x1+h,x2)-function(x1-h,x2))/(2*h),(function(x1,x2+h)-function(x1,x2-h))/(2*h)

def parameter_matrix(zeros,function1,function2,h=0.001,decimal=5):
    A=np.zeros((2,2))
    A[0,0],A[0,1]=df.partial_derivative_3p(zeros[0],zeros[1],function1,h)
    A[1,0],A[1,1]=df.partial_derivative_3p(zeros[0],zeros[1],function2,h)
    return np.round(A,decimal)
def trace_and_determinant(A):
    return (A[0,0]+A[1,1],A[0,0]*A[1,1]-A[0,1]*A[1,0])
    
def dx1_dt_dx2_dt_rk2(t_initial,t_final,x1_initial,x2_initial,dx1dt,dx2dt,dt=0.01):
    n=int((t_final-t_initial)/dt)
    if n<0:
        n=-n
        dt=-dt
    x1=np.zeros(n+1)
    x2=np.zeros(n+1)
    t=np.zeros(n+1)
    t[0]=t_initial
    x1[0]=x1_initial
    x2[0]=x2_initial
    for i in range(n):
        k1_x1,k1_x2=dx1dt(x1[i],x2[i],t[i]),dx2dt(x1[i],x2[i],t[i])
        k2_x1,k2_x2=dx1dt(x1[i]+dt/2*k1_x1,x2[i]+dt/2*k1_x2,t[i]+dt/2),dx2dt(x1[i]+dt/2*k1_x1,x2[i]+dt/2*k1_x2,t[i]+dt/2)
        x1[i+1]=x1[i]+dt*k2_x1
        x2[i+1]=x2[i]+dt*k2_x2
        t[i+1]=t[i]+dt
    return x1,x2,t

def near_trajectory(function1,function2,center,points=4,radius=0.1,t_final=1,dt=0.05):
    info={}
    for i in range(points):
        x1ini=radius*np.cos(2*np.pi*i/points)+center[0]
        x2ini=radius*np.sin(2*np.pi*i/points)+center[1]
        info['x1_'+str(i+1)],info['x2_'+str(i+1)],info['t_'+str(i+1)]=dx1_dt_dx2_dt_rk2(0,t_final,x1ini,x2ini,lambda x1,x2,t:function1(x1,x2),lambda x1,x2,t:function2(x1,x2),dt)
    return info

def unique(nlist):
    n=len(nlist)
    uni=[]
    for i in range(n):
        s=0
        for j in range(i+1,n):
            if nlist[i]==nlist[j]:
                s=1
        if s==0:
            uni.append(nlist[i])
    return uni

def newton_rapshon_2d(x1_start,x2_start,function1,function2,epsilon=1e-6,h=0.001,max_iteration=1000):
    temp1,temp2=2*epsilon,2*epsilon
    x1=x1_start
    x2=x2_start
    itr=0
    if function1(x1,x2)==0 and function2(x1,x2)==0:
        return x1,x2,itr
    x=np.array([x1,x2]).reshape(2,1)
    f=np.zeros((2,1))
    empty=np.zeros((2,2))
    
    while (abs(temp1)>epsilon or abs(temp2)>epsilon):
        f[0,0],f[1,0]=function1(x[0,0],x[1,0]),function2(x[0,0],x[1,0])
        empty[0,0],empty[0,1]=partial_derivative_3p(x[0,0],x[1,0],function1,h)
        empty[1,0],empty[1,1]=partial_derivative_3p(x[0,0],x[1,0],function2,h)
        x=x-np.dot(np.linalg.inv(empty),f)
        temp1,temp2=function1(x[0,0],x[1,0]),function2(x[0,0],x[1,0])
        if itr>max_iteration:
            print('newton rapshon 2d method may not converging')
            return x[0,0],x[1,0]
        itr=itr+1
    return x[0,0],x[1,0]

def non_linear_zeros_2d(function1,function2,box_size=10,search_time=20,decimal=3,h=0.001,epsilon=1e-6,max_iteration=100):
    t1=np.random.uniform(-box_size,box_size,search_time)
    t2=np.random.uniform(-box_size,box_size,search_time)
    zeros=[]
    for i in range(search_time):
        x1_start=t1[i]
        x2_start=t2[i]
        x1,x2=newton_rapshon_2d(x1_start,x2_start,function1,function2,epsilon,h,max_iteration)
        zeros.append([round(x1,decimal),round(x2,decimal)])
    unique_zeros=unique(zeros)
    return unique_zeros

def classify(trace,determinant,matrix=None,decimal=5):
    t=trace
    d=determinant
    disc=t**2-4*d
    if t==0 and d==0:
        return 'plane of fixed points(non isolated fixed points)'
    elif t==0 and d>0:
        return 'center'
    elif t>0 and disc<0:
        return 'unstable spiral'
    elif t<0 and disc<0:
        return 'stable spiral'
    elif disc==0:
        if type(matrix)!=type(None):
            eigenval,eigenvec=np.linalg.eig(matrix)
            eigenvec=np.round(eigenvec,decimal)
            e1,e2=eigenvec[:,0],eigenvec[:,1]
            if (np.array_equal(e1,e2) or np.array_equal(e1,-e2)):
                if t>0:
                    return 'unstable degenerate node'
                if t<0:
                    return 'stable degenerate node'
            else:
                if t>0:
                    return 'unstable star'
                if t<0:
                    return 'stable star'
        else:
            print('No matrix given to check conditions')
            return 'None'
    elif disc>0 and t>0 and d>0:
        return 'unstable node'
    elif disc>0 and t<0 and d>0:
        return 'stable node'
    elif d==0 and t>0:
        return 'unstable line of fixed points(non islolated fixed points)'
    elif d==0 and t<0:
        return 'stable line of fixed points(non islolated fixed points)'
    elif d<0:
        return 'saddle'
    else:
        return 'cannot classify'
    
def classify_fixed_points(function1,function2,box_size=50,search_time=20,points=12,radius=0.05,t_final=0.5,dt=0.01,arrow=2,detail=True,plot=True,decimal=5,h=0.001,epsilon=1e-6,max_iteration=500):
    zeros=non_linear_zeros_2d(function1,function2,box_size,search_time,decimal,h,epsilon,max_iteration)
    for val in zeros:
        x_zero,y_zero=val[0],val[1]
        A=parameter_matrix((x_zero,y_zero),function1,function2,h,decimal)
        trace,determinant=trace_and_determinant(A)
        name='Fixed point = '+'('+str(x_zero)+' , '+str(y_zero)+')'
        ans=classify(trace,determinant,A,decimal)
        print(name+' : '+ans)
        if detail:
            print('system matrix at '+'('+str(x_zero)+' , '+str(y_zero)+')'+'=')
            print(A)
        if plot:
            check=True
            if ans=='plane of fixed points(non isolated fixed points)' or ans=='cannot classify':
                check=False
            if check:
                center=(x_zero,y_zero)
                info=near_trajectory(function1,function2,center,points,radius,t_final,dt=dt)
                plt.figure()
                plt.scatter(center[0],center[1],marker='o',color='k', s=150)
                for i in range(points):
                    plt.plot(info['x1_'+str(i+1)],info['x2_'+str(i+1)],color='r')
                    plt.arrow(info['x1_'+str(i+1)][0],info['x2_'+str(i+1)][0],info['x1_'+str(i+1)][arrow]-info['x1_'+str(i+1)][0],info['x2_'+str(i+1)][arrow]-info['x2_'+str(i+1)][0],color='b',head_width=0.005,head_length=0.005)
                plt.title(name+' : '+ans)
                plt.show()
