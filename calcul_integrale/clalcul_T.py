import scipy.integrate as np
def calcul_T(T0,T1,p,Cp):
    """calcul de la fonction primitive de Cp(T), discrétisée avec un pas de 1/p
    sur l'intervalle [T0,T1], entiers, p entier"""

    R1=[0]
    R2=[0]#liste vide contenat les futurs couple (T,F(T))
    for i in range((T1-T0)*p):
        R1.append(T0+i*(1/p))
        R2.append(R2[-1]+np.quad(Cp,T0+i*(1/p),T0+(i+1)*(1/p))[0])
    R1.pop(0)
    R2.pop(0)
    return [R1,R2]
    
def f(T):
    a=2*sqrt(T)    
    
    return 3.5-0.000028*T+0.0000000224*(T**2)+(3000*3000*a)/(T*(a-1)*(a-1))

import matplotlib.pyplot as pypl
def trace(T0,T1,p,Cp):
    R1=calcul_T(T0,T1,p,Cp)[0]
    R2=calcul_T(T0,T1,p,Cp)[1]
    pypl.plot(R1,R2)
    pypl.show()