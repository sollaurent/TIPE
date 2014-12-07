def monte_carlo(N):
"Prend en entree le nombre d'essais N
"Retourne une liste contenant N listes de paramètres d'entrees
"Chaque liste correspond a une experience, les parametres etant choisis au hasard
    entrees=[]
    for i in range (N):
        experience=[]
        z=random.randint(-2000,50000)       #altitude
        M=random(0,2)                       #nb de Mach
        VC=random()                         #vitesse conventionnelle_domaine de test a preciser
        VA=random()                         #vitesse aerodynamique_domaine de test a preciser
        ISA=random()                        #ecart de temperature/atmosphere standard_domaine de test a preciser
        experience.append(z,M,VC,VA,ISA)
        entrees.append experience
    return entrees
    
    
    

