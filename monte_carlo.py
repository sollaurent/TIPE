

#fonction donnant la température et la pression en fonction de l'altitude    
ISA_temp = interpolate.interp1d([0,11000,20000,32000],[288,216.5,216.5,228.5])
ISA_P = interpolate.interp1d([0,11000,20000,32000],[101325,22632,5474.9,868.02])


def monte_carlo(N):
"Prend en entree le nombre d'essais N
"Retourne une liste contenant N listes de paramètres d'entrees
"Chaque liste correspond a une experience, les parametres etant choisis au hasard
    entrees=[]
    for i in range (N):
        experience=[]
        
        #Altitude
        z=randint(0,32000)
        
        TC = ISA_temp(z)
        PC= ISA_P(z)
        
        #Nombre de Mach
        M=random(0.5,2.5)
        #Temperature totale
        TT=TC*(1+0.2*M**2)
        #Pression totale
        PT=PC*(1+0.2*M**2)**3.5
        #Taux de compression
        ts=randint(1,10)
        tcbp=randint(1,10)
        tchp=randint(1,10)
        tt=randint(1,10)
        #Richesse
        alpha=random(0,1)
        #Coefficient de partage du flux
        lamb=random(0,1)
        #Flux
        WA=randint(100;1000)
        WF=randint(100,1000)
        #Vitesse
        VA=M*(1.4*237*TC)**0.5
        a=VA/M        
        
        experience.append(TT,PT,ts,tcbp,thbp,tt,alpha,lamb,WA,WF,a,turboreacteur(TT,PT,ts,tcbp,thbp,tt,alpha,lamb,WA,WF,a))
        entrees.append(experience)
    return entrees
    
    
T1,P1,ts,tcbp,tchp,tt,alpha,lamb,WA,WF,a
    
    
    

