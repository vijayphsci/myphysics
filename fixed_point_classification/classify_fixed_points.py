import numpy as np
import fixedpointclassification as ls

def classify_fixed_point(zeros,function1,function2,detail=True,h=0.001,decimal=5):
    A=ls.parameter_matrix(zeros,function1,function2,h,decimal)
    trace,determinant=ls.trace_and_determinant(A)
    ans=ls.classify(trace,determinant,A,detail,decimal)
    return ans

def classify_fixed_point_automatic(function1,function2,box_size=50,search_time=20,detail=True,decimal=5,h=0.001,epsilon=1e-6,max_iteration=500):
    zeros=ls.non_linear_zeros_2d(function1,function2,box_size,search_time,decimal,h,epsilon,max_iteration)
    for val in zeros:
        x_zero,y_zero=val[0],val[1]
        A=ls.parameter_matrix((x_zero,y_zero),function1,function2,h,decimal)
        trace,determinant=ls.trace_and_determinant(A)
        name='Fixed point = '+'('+str(x_zero)+' , '+str(y_zero)+')'
        ans=ls.classify(trace,determinant,A,detail,decimal,set_title=name)
        print(name+' : '+ans)
        if detail:
            print('system matrix at '+'('+str(x_zero)+' , '+str(y_zero)+')'+'=')
            print(A)

