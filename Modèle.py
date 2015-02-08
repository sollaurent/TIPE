from calcul_T import *
import scipy
from scipy.integrate import quad


def gamma_temperature(fluide,alpha,a):
    """ Prend en entrée la nature du fluide : 0=air, 1=kerosene, 2= melange air-kerosene,
    la richesse alpha (valant 0 pour l'air et 1 pour le kerosene pur),
    la vitesse du son a.
    Renvoi gamma(T)"""
    if fluide==0:
        Cpr=lambda T: 3.5-0.000028*T+0.0000000224*T*T+(3000*3000*a)/(T*(a-1)*(a-1))
        gamma=lambda T: Cpr(T)/(Cpr(T)-1)
    if fluide==1:
        Cpr=lambda T: -0.0000018373*T*T+0.00801994*T+4.47659
        gamma=lambda T: Cpr(T)/(Cpr(T)-1)
    if fluide==2:
        Cpra=lambda T: 3.5-0.000028*T+0.0000000224*T*T+(3000*3000*a)/(T*(a-1)*(a-1))
        Cprk=lambda T: -0.0000018373*T*T+0.00801994*T+4.47659
        Cpr=lambda T: (Cpra(T)+alpha*Cprk(T))/(alpha+1)
        gamma=lambda T: Cpr(T)/(Cpr(T)-1)
    return gamma
    
def Cp_temperature(fluide,alpha,a):
    """ Prend en entrée la nature du fluide : 0=air, 1=kerosene, 2= melange air-kerosene,
    la richesse alpha (valant 0 pour l'air et 1 pour le kerosene pur),
    la vitesse du son a.
    Renvoi Cp(T)"""
    if fluide==0:
        Cpr=lambda T: 3.5-0.000028*T+0.0000000224*T*T+(3000*3000*a)/(T*(a-1)*(a-1))
    if fluide==1:
        Cpr=lambda T: -0.0000018373*T*T+0.00801994*T+4.47659
    if fluide==2:
        Cpra=lambda T: 3.5-0.000028*T+0.0000000224*T*T+(3000*3000*a)/(T*(a-1)*(a-1))
        Cprk=lambda T: -0.0000018373*T*T+0.00801994*T+4.47659
        Cpr=lambda T: (Cpra(T)+alpha*Cprk(T))/(alpha+1)
    return Cpr
    
def turboreacteur(T1,P1,ts,tcbp,tchp,tt,alpha,lamb,WA,WF,a):
    """Modélisation d'un turboréacteur
    Hypothèses : transformations isentropiques dans com  """
    
    #calcul des fonctions utilisés dans les blocs du turboréacteur
    f1=lambda T: gamma_temperature(0,0,a)(T)/(T*(gamma_temperature(0,0,a)(T)-1)) #fonction gamma/T*(gamma-1) de l'air
    F1=calcul_T(200,2500,5,f1) #calcul de la primitive de gamma/T*(gamma-1) pour l'air
    Finv1=[F1[1],F1[0]] #on inverse les deux listes pour determiner l'inverse

    f2=lambda T : 287*Cp_temperature(2,alpha,a)(T) #fonction 287*Cp pour mélange air/essence
    F2=calcul_T(200,2500,5,f2) #calcul de la primitive de 287*Cp pour mélange air/essence
    Finv2=[F2[1],F2[0]] #calcul de l'inverse
    
    def Finv1_calcul(K):
        """recherche de la valeur de Finv1 la plus proche de K"""
        assert K>Finv1[0][0] and K<Finv1[0][-1] #on vérifie que K est compris dans la plage connue
        for i in range(len(Finv1[0])): #dès que K passe au dessus d'une valeur de
        #Finv1[0), on lui associe la valeur correspondante, arrondie par défaut
            if K>Finv1[0][i]:
                return Finv1[1][i]

    def Finv2_calcul(K):
        """recherche de la valeur de Finv2 la plus proche de K"""
        assert K>Finv2[0][0] and K<Finv2[0][-1] #on vérifie que K appartient à la plage connue
        for i in range (len(Finv1[0])): #arrondi par défaut de la valeur de K par rapport à la plage connue
            if K>Finv2[0][i]:
                return Finv2[1][i]
              
            
    
    #Soufflante, obtenu par integration(Laplace adapté)
    T2=Finv1_calcul(F1_cacul(T1)+ln(ts))
    
    #Compresseur BP, obtenu par integration(Laplace adapté)
    T3=Finv1_calcul(F1_cacul(T2)+ln(tcbp))
        
    #Compresseur HP, obtenu par integration(Laplace adapté)
    T4=Finv1_calcul(F1_cacul(T3)+ln(tchp))
        
    #Chambre de combustion, 1er principe thermochimie
    DfCO2=394000
    DfH2O=280000
    DfN2=0
    DfO2=0
    DfCH4=88000
    avanct=WF/(0.012+4*0.001)
    Df=avanct(DfCO2+2*DfH2O-DfCH4)
    
    T5=Finv2_calcul(F2(T4)-Df)
    
    #Turbine HP, obtenu par equilibre HP
    T6=Finv2_calcul(F2(T5)+(F2(T4)-F2(T3))*WA/(WA+WF))
    
    #Turbine BP, obtenu par equilibre BP
    T7=Finv2_calcul(F2(T6)+(F2(T3)-F2(T1))*WA/(WA+WF)+(F2(T1)-F2(T2))*lamb*WA/(WA+WF))
    
    #Mélangeur, application 1er principe
    T8=Finv2_calcul(((WA+WF)*F2(T7)+lamb*WA*F2(T2))/(WF+WA+WA*lamb))
    
    #Tuyère, application 1er principe
    T9=Finv1_calcul(F1(T8)+ln(tt))
    C9=sqrt(2+F2(T9)-F2(T8))
    
    #Rendement
    Pcin=(WA+lamb*WA+WF)*C9*C9/2
    Pth=(WA+WF)*(F2(T5)-F2(T4))
    Rendement=Pcin/Pth
    
    return Rendement
    
    
    
    
    
    
    


    
    
"""    
#Soufflante, obtenu par integration(Laplace adapté)
T2=Finv1_calcul(F1(T1)+ln(ts))
    
#Compresseur BP, obtenu par integration(Laplace adapté)
T3=Finv1_calcul(F1(T2)+ln(tcbp))
        
#Compresseur HP, obtenu par integration(Laplace adapté)
T4=Finv1_calcul(F1(T3)+ln(tchp))
        
#Chambre de combustion, 1er principe thermochimie
DfCO2=394000
DfH2O=280000
DfN2=0
DfO2=0
DfCH4=88000
avanct=WF/(0.012+4*0.001)
Df=avanct(DfCO2+2*DfH2O-DfCH4)

T5=Finv2_calcul(F2(T4)-Df)

#Turbine HP, obtenu par equilibre HP
T6=Finv2_calcul(F2(T5)+(F2(T4)-F2(T3))*WA/(WA+WF))

#Turbine BP, obtenu par equilibre BP
T7=Finv2_calcul(F2(T6)+(F2(T3)-F2(T1))*WA/(WA+WF)+(F2(T1)-F2(T2))*lamb*WA/(WA+WF))

#Mélangeur, application 1er principe
T8=Finv2_calcul(((WA+WF)*F2(T7)+lamb*WA*F2(T2))/(WF+WA+WA*lamb))

#Tuyère, application 1er principe
T9=Finv1_calcul(F1(T8)+ln(tt))
C9=sqrt(2+F2(T9)-F2(T8))

#Rendement
Pcin=(WA+lamb*WA+WF)*C9*C9/2
Pth=(WA+WF)*(F2(T5)-F2(T4))
Rendement=Pcin/Pth
    
"""