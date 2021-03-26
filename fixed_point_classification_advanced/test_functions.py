from fixed_point_classify_vj import classify_fixed_points
"""
dx1/dt = function1(x1,x2)

dx2/dt = function2(x1,x2)

"""
def function1(x1,x2):
    return  x1*(3-x1-2*x2)

def function2(x1,x2):
    return x2*(2-x1-x2)

box_size=5
search_time=50
detail=True
decimal=5
h=0.001
epsilon=1e-6
max_iteration=500
plot=True
points=12
radius=0.05
t_final=1
dt=0.05
arrow=2

classify_fixed_points(function1=lambda x1,x2:function1(x1,x2),function2=lambda x1,x2:function2(x1,x2),box_size=box_size,search_time=search_time,points=points,radius=radius,t_final=t_final,dt=dt,arrow=arrow,detail=detail,plot=plot,decimal=decimal,h=h,epsilon=epsilon,max_iteration=max_iteration)