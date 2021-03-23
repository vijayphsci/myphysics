import classify_fixed_points as cfp

def function1(x1,x2):
    return x2
def function2(x1,x2):
    return -x1-x2


box_size=50
search_time=5
detail=True
decimal=5
h=0.001
epsilon=1e-6
max_iteration=500



cfp.classify_fixed_point_automatatic(function1=lambda x1,x2:function1(x1,x2),function2=lambda x1,x2:function2(x1,x2),box_size=box_size,search_time=search_time,detail=detail,decimal=decimal,h=h,epsilon=epsilon,max_iteration=max_iteration)

