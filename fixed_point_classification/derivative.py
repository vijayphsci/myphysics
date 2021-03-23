def derivative_2p(x,function,h=0.001):
    return (function(x+h)-function(x))/h
def derivative_3p(x,function,h=0.001):
    return (function(x+h)-function(x-h))/(2*h)
def derivative_5p(x,function,h=0.001):
    return (function(x-2*h)+8*(function(x+h)-function(x-h))-function(x+2*h))/(12*h)
def second_derivative_3p(x,function,h=0.001):
    return (function(x+h)+function(x-h)-2*function(x))/h**2
def second_derivative_5p(x,function,h=0.001):
    return (-function(x-2*h)+16*(function(x+h)+function(x-h))-function(x+2*h)-30*function(x))/(12*h**2)
def partial_derivative_2p(x1,x2,function,h=0.001):
    """
    function : lambda x1,x2=f(x1,x2)
    return df/dx1,df/dx2 at x1,x2
    """
    return (function(x1+h,x2)-function(x1,x2))/h,(function(x1,x2+h)-function(x1,x2))/h
def partial_derivative_3p(x1,x2,function,h=0.001):
    """
    function : lambda x1,x2=f(x1,x2)
    return df/dx1,df/dx2 at x1,x2
    """
    return (function(x1+h,x2)-function(x1-h,x2))/(2*h),(function(x1,x2+h)-function(x1,x2-h))/(2*h)



