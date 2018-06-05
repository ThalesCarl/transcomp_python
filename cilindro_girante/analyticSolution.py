def analyticSolution(x,y):
    re = 0.1 #[m]
    ri = 0.04 #[m]
    Ti = 250.0 #[ºC]
    Te = 30.0 #[ºC]
    
    import math
    
    r = math.sqrt(x*x+y*y)
    if r<ri or r>re:
        print("Raio inválido")
        T = math.nan
    else:
        T = Ti + (Te-Ti)*(math.log(r/ri)/math.log(re/ri))
    return T

