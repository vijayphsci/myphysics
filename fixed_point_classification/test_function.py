import classify_fixed_points as cfp

def function1(x1,x2):
    return x2
def function2(x1,x2):
    return -x1

cfp.classify_fixed_point_automatatic(lambda x1,x2:function1(x1,x2),lambda x1,x2:function2(x1,x2),search_time=5)