import matplotlib.pyplot as plt
import numpy as np
import derivative as df
import matplotlib.image as mpimg
def parameter_matrix(zeros,function1,function2,h=0.001,decimal=5):
    A=np.zeros((2,2))
    A[0,0],A[0,1]=df.partial_derivative_3p(zeros[0],zeros[1],function1,h)
    A[1,0],A[1,1]=df.partial_derivative_3p(zeros[0],zeros[1],function2,h)
    return np.round(A,decimal)
def trace_and_determinant(A):
    return (A[0,0]+A[1,1],A[0,0]*A[1,1]-A[0,1]*A[1,0])

def classify(trace,determinant,matrix=None,detail=True,decimal=5,set_title='fixed point'):
    t=trace
    d=determinant
    disc=t**2-4*d
    if t==0 and d==0:
        return 'plane of fixed points(non isolated fixed points)'
    elif t==0 and d>0:
        if detail:
            img = mpimg.imread('center.gif')
            plt.figure()
            fig = plt.imshow(img)
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            plt.title(set_title)
            plt.show()
        return 'center'
    elif t>0 and disc<0:
        if detail:
            img = mpimg.imread('unstablespiral.PNG')
            plt.figure()
            fig = plt.imshow(img)
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            plt.title(set_title)
            plt.show()
        return 'unstable spiral'
    elif t<0 and disc<0:
        if detail:
            img = mpimg.imread('stablespiral.PNG')
            plt.figure()
            fig = plt.imshow(img)
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            plt.title(set_title)
            plt.show()
        return 'stable spiral'
    elif disc==0:
        if type(matrix)!=type(None):
            eigenval,eigenvec=np.linalg.eig(matrix)
            eigenvec=np.round(eigenvec,decimal)
            e1,e2=eigenvec[:,0],eigenvec[:,1]
            if (np.array_equal(e1,e2) or np.array_equal(e1,-e2)):
                if t>0:
                    if detail:
                        img = mpimg.imread('unstabledegnode.PNG')
                        plt.figure()
                        fig = plt.imshow(img)
                        fig.axes.get_xaxis().set_visible(False)
                        fig.axes.get_yaxis().set_visible(False)
                        plt.title(set_title)
                        plt.show()
                    return 'unstable degenerate node'
                if t<0:
                    if detail:
                        img = mpimg.imread('stabledegnode.PNG')
                        plt.figure()
                        fig = plt.imshow(img)
                        fig.axes.get_xaxis().set_visible(False)
                        fig.axes.get_yaxis().set_visible(False)
                        plt.title(set_title)
                        plt.show()
                    return 'stable degenerate node'
            else:
                if t>0:
                    if detail:
                        img = mpimg.imread('unstablestar.PNG')
                        plt.figure()
                        fig = plt.imshow(img)
                        fig.axes.get_xaxis().set_visible(False)
                        fig.axes.get_yaxis().set_visible(False)
                        plt.title(set_title)
                        plt.show()
                    return 'unstable star'
                if t<0:
                    if detail:
                        img = mpimg.imread('stablestar.PNG')
                        plt.figure()
                        fig = plt.imshow(img)
                        fig.axes.get_xaxis().set_visible(False)
                        fig.axes.get_yaxis().set_visible(False)
                        plt.title(set_title)
                        plt.show()
                    return 'stable star'
        else:
            print('No matrix given to check conditions')
            return None
    elif disc>0 and t>0 and d>0:
        if detail:
            img = mpimg.imread('unstablenode.PNG')
            plt.figure()
            fig = plt.imshow(img)
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            plt.title(set_title)
            plt.show()
        return 'unstable node'
    elif disc>0 and t<0 and d>0:
        if detail:
            img = mpimg.imread('stablenode.PNG')
            plt.figure()
            fig = plt.imshow(img)
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            plt.title(set_title)
            plt.show()
        return 'stable node'
    elif d==0 and t>0:
        if detail:
            img = mpimg.imread('unstableline.PNG')
            plt.figure()
            fig = plt.imshow(img)
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            plt.title(set_title)
            plt.show()
        return 'unstable line of fixed points(non islolated fixed points)'
    elif d==0 and t<0:
        if detail:
            img = mpimg.imread('stableline.PNG')
            plt.figure()
            fig = plt.imshow(img)
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            plt.title(set_title)
            plt.show()
        return 'stable line of fixed points(non islolated fixed points)'
    elif d<0 and t>0:
        if detail:
            img = mpimg.imread('saddle1.PNG')
            plt.figure()
            fig = plt.imshow(img)
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            plt.title(set_title)
            plt.show()
        return 'saddle kind1'
    elif d<0 and t<0:
        if detail:
            img = mpimg.imread('saddle2.PNG')
            plt.figure()
            fig = plt.imshow(img)
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            plt.title(set_title)
            plt.show()
        return 'saddle kind2'
    
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
        empty[0,0],empty[0,1]=df.partial_derivative_3p(x[0,0],x[1,0],function1,h)
        empty[1,0],empty[1,1]=df.partial_derivative_3p(x[0,0],x[1,0],function2,h)
        x=x-np.dot(np.linalg.inv(empty),f)
        temp1,temp2=function1(x[0,0],x[1,0]),function2(x[0,0],x[1,0])
        if itr>max_iteration:
            print('newton rapshon 2d method may not converging')
            return x[0,0],x[1,0]
        itr=itr+1
    return x[0,0],x[1,0]


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

