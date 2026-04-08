import math

def Omega(a, b, d):
    
    arg = a*b/(math.sqrt((a**2 + d**2)*(b**2 + d**2)))

    return math.asin(arg)

def Surface(a, b):

    return a*b

def Surface_err(a, d_a, b, d_b):

    err_1 = b * d_a

    err_2 = a * d_b

    return math.sqrt(err_1**2 + err_2**2)

def Omega_err(a, d_a, b, d_b, d, d_d):

    O = Omega(a, b, d)

    err_1 = (math.tan(O)/(1+math.tan(O)**2))*(1 + 1/(1+(d/a)**2))*(d_a/a)

    err_2 = (math.tan(O)/(1+math.tan(O)**2))*(1 + 1/(1+(d/b)**2))*(d_b/b)

    err_3 = (math.tan(O)/(1+math.tan(O)**2)) * (d * d_d * (a**2 + b**2 + 2*d**2))/(d**4 + d**2 * (a**2 + b**2) + (a*b)**2)

    return math.sqrt(err_1**2 + err_2**2 + err_3**2)

print(Omega_err(58, 0.1, 21, 0.1, 3.8, 0.1), Omega_err(40, 0.1, 20, 0.1, 60.5, 0.1), Omega_err(50, 0.1, 40, 0.1, 41, 0.1))