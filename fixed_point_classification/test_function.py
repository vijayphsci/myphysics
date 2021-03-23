import classify_fixed_points as cfp

def function1(x1,x2):
    return x1*(3-x1-2*x2)
def function2(x1,x2):
    return x2*(2-x1-x2)


box_size=5
search_time=50
detail=True
decimal=5
h=0.001
epsilon=1e-6
max_iteration=500



cfp.classify_fixed_point_automatic(function1=lambda x1,x2:function1(x1,x2),function2=lambda x1,x2:function2(x1,x2),box_size=box_size,search_time=search_time,detail=detail,decimal=decimal,h=h,epsilon=epsilon,max_iteration=max_iteration)
