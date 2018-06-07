"""
Aluno: Thales Carl Lavoratti (151000656)
Funções auxiliares usadas no código do programa do cilindro girante
"""
def analyticSolution(x,y):
    re = 0.1 #[m]
    ri = 0.04 #[m]
    Ti = 250.0 #[ºC]
    Te = 30.0 #[ºC]
    
    import math as mt
    r = mt.sqrt(x*x+y*y)
    if r<ri or r>re:
        print("Raio inválido")
        T = mt.nan
    else:
        T = Ti + (Te-Ti)*(mt.log(r/ri)/mt.log(re/ri))
    return T

def uVelocity(omega,x,y):
    import math as mt
    r = mt.sqrt(x*x+y*y)
    theta = mt.atan(y/x)
    tangencialVelocity = omega*r
    return (-1)*mt.sin(theta)*tangencialVelocity

def vVelocity(omega,x,y):
    import math as mt
    r = mt.sqrt(x*x+y*y)
    theta = mt.atan(y/x)
    tangencialVelocity = omega*r
    return mt.cos(theta)*tangencialVelocity
    