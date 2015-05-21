from calcul_T import *
import scipy
from math import *
from scipy import interpolate


def gamma_temperature(fluide,alpha,a):
    """ Prend en entrée la nature du fluide : 0=air, 1=kerosene, 2= melange air-kerosene,
    la richesse alpha (valant 0 pour l'air et 1 pour le kerosene pur),
    la vitesse du son a.
    Renvoi gamma(T)"""
    if fluide==0:
        Cpr=lambda T: 3.5-0.000028*T+0.0000000224*T*T+(3090*3090*a)/(T*T*(a-1)*(a-1))
        gamma=lambda T: Cpr(T)/(Cpr(T)-1)
    if fluide==1:
        Cpr=lambda T: -0.0000018373*T*T+0.00801994*T+4.47659
        gamma=lambda T: Cpr(T)/(Cpr(T)-1)
    if fluide==2:
        Cpra=lambda T: 3.5-0.000028*T+0.0000000224*T*T+(3090*3090*a)/(T*T*(a-1)*(a-1))
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
        Cpr=lambda T: 3.5-0.000028*T+0.0000000224*T*T+(3090*3090*a)/(T*T*(a-1)*(a-1))
    if fluide==1:
        Cpr=lambda T: -0.0000018373*T*T+0.00801994*T+4.47659
    if fluide==2:
        Cpra=lambda T: 3.5-0.000028*T+0.0000000224*T*T+(3090*3090*a)/(T*T*(a-1)*(a-1))
        Cprk=lambda T: -0.0000018373*T*T+0.00801994*T+4.47659
        Cpr=lambda T: (Cpra(T)+alpha*Cprk(T))/(alpha+1)
    return Cpr
    
def turboreacteur(T1,P1,ts,tcbp,tchp,Tcomb,lamb,WA,VA,rs,rcbp,rchp,rtbp,rthp):
    """Modélisation d'un turboréacteur
    Hypothèses : transformations isentropiques dans compresseurs
    avec prise en compte des rendements   """
    
    #calcul des fonctions utilisés dans les blocs du turboréacteur
    fp1=lambda T: gamma_temperature(0,0,a)(T)/(T*(gamma_temperature(0,0,a)(T)-1)) #fonction gamma/T*(gamma-1) de l'air
    F1=calcul_T(200,2500,1,fp1) #calcul de la primitive de gamma/T*(gamma-1) pour l'air
    
    
    f1 = interpolate.interp1d(F1[0],F1[1])#fonction primitive
    finv1=interpolate.interp1d(F1[1],F1[0]) #reciproque

    fp2=lambda T : 287*Cp_temperature(2,alpha,a)(T) #fonction 287*Cp pour mélange air/essence
    F2=calcul_T(200,2500,1,fp2) #calcul de la primitive de 287*Cp pour mélange air/essence
    
    f2 = interpolate.interp1d(F2[0],F2[1])#fonction primitive
    finv2=interpolate.interp1d(F2[1],F2[0]) #fonction l'inverse            
            
    
    #Soufflante, obtenu par integration(Laplace adapté)
    T2=finv1(f1(T1)+log(ts)) #log correspond au logarithme néperien
    #Compresseur BP, obtenu par integration(Laplace adapté)
    T3=finv1(f1(T2)+log(tcbp))#log correspond au logarithme néperien
        
    #Compresseur HP, obtenu par integration(Laplace adapté)
    T4=finv1((f1(T3)+log(tchp)))#log correspond au logarithme néperien
        
    #Chambre de combustion, 1er principe thermochimie
    DfCO2=394000#DfCO2
    DfH2O=280000#DfH20
    DfCH4=88000
    DH=DfCO2-DfCH4-2*DfH2O #enthalpie de reaction
    MCH4(0.012+4*0.001)#masse molaire CH4
    T5=Tcomb
    
    alpha=MCH4/(WA*DH)*(f2(T5)-f2(T4))
    WF=alpha*WA

    print((T4))
    
    #T5=finv2(f2(T4)+Df)
    
    print(T5)
    
    #Turbine HP, obtenu par equilibre HP
    T6=finv2(f2(T5)-(f2(T4)-f2(T3))*WA/(WA+WF)/rchp/rthp)
    
    #Turbine BP, obtenu par equilibre BP
    T7=finv2(f2(T6)-(f2(T3)-f2(T1))*WA/(WA+WF)/rcbp/rtbp-(f2(T2)-f2(T1))*lamb*WA/(WA+WF)/rs/rtbp)
    
    #Mélangeur, application 1er principe
    T8=finv2(((WA+WF)*f2(T7)+lamb*WA*f2(T2))/(WF+WA+WA*lamb))

    
    #Tuyère, application 1er principe
    T9=finv1(f1(T8)+log(tt))
    C9=sqrt(2+f2(T9)-f2(T8))
    
    #Rendement
    Pcin=((WA+lamb*WA+WF)*C9*C9-(1+lamb)*WA*VA*VA)/2
    Pth=(WA+WF)*(f2(T5)-f2(T4))
    Rendement=Pcin/Pth
    
    return Rendement,alpha
    
    
