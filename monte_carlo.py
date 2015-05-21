from Modèle import *
from random import *

#fonction donnant la température et la pression en fonction de l'altitude    
ISA_temp = interpolate.interp1d([0,11000,20000,32000],[288,216.5,216.5,228.5])
ISA_P = interpolate.interp1d([0,11000,20000,32000],[101325,22632,5474.9,868.02])


def monte_carlo(N):
    """Prend en entree le nombre d'essais N
    Retourne une liste contenant N listes de paramètres d'entrees
    Chaque liste correspond a une experience, les parametres etant choisis au 
    hasard"""
    entrees=[]
    resultat = open("resulat.txt","a")#crée fichier resultat
    for i in range (N):
        experience=[]
        
        #Altitude
        z=randint(0,15000)
        
        TC = ISA_temp(z)
        PC= ISA_P(z)
        
        #Nombre de Mach
        M=uniform(0,0.99)#loi uniforme sur 0 0.99
        #Temperature totale
        TT=TC*(1+0.2*M**2)
        #Pression totale
        PT=PC*(1+0.2*M**2)**3.5
        #Taux de compression
        ts=randint(0,50)#soufflante
        tcbp=randint(0,50)#compresseur bp
        tchp=randint(0,50)#compresseur hp
        #rendements
        rs=random(0.5,1)
        rcbp=random(0.5,1)
        rchp=random(0.5,1)
        rtbp=random(0.5,1)
        rthp=random(0.5,1)
        #Temperature chambre combustion desiree
        Tcomb=uniform(1700,2500)
        #Coefficient de partage du flux
        lamb=random()*100
        #Flux
        WA=randint(100,1000)
        #Vitesse
        VA=M*(1.4*237*TC)**0.5        
        
        experience=[TT,PT,ts,tcbp,tchp,Tcomb,lamb,WA,VA,rs,rcbp,rchp,rtbp,rthp,turboreacteur(TT,PT,ts,tcbp,tchp,Tcomb,lamb,WA,VA,rs,rcbp,rchp,rtbp,rthp)]
        entrees.append(experience)
        resultat.write(str(experience))#ajoute au fichier
    resultat.close()
    return entrees
