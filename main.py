"""
O codigo abaixo visa aplicacao da simulacao de Monte Carlo na resolucao de 
trelicas planas e espaciais de acordo com diferentes regimes de 
elastoplasticidade.
Substitua neste arquivo os dados da trelica a ser simulada.
Data: 27/11/2020
@autor: Gabriel de Carvalho Freitas
"""
from datetime import datetime
runtime = datetime.now() # cronometragem do tempo de execução
import numpy as np
import math
from scipy.stats import norm, lognorm
from nonlinear_solving import nonlinear_solving

# parâmetros nominais de material e geometria dos elementos
Diam_nom = 11.4 # Diâmetro nominal do tubo circular vazado. Unidade: cm
t_nom = 1.0 # Espessura do tubo. Unidade: cm
E_nom = 20000 # Módulo de elasticidade do material. Unidade: kN/cm^2
sigmay_nom = 24 # Tensão de escoamento do material. Unidade: kN/cm^2
P_nom = 390 # Carregamento imposto (força ou deslocamento). Unidade: kN

# parâmetros da simulação de monte carlo
MC_rep = 125000 # número de repetições da simulação de monte carlo
                
MC_Diam = norm.rvs(size=MC_rep,loc=1.0059*Diam_nom,scale=1.0059*.00181*Diam_nom)
MC_t = norm.rvs(size=MC_rep,loc=1.0069*t_nom,scale=1.0069*.0259*t_nom)
MC_E = norm.rvs(size=MC_rep,loc=1.04*E_nom,scale=1.04*.05*E_nom)
MC_P = norm.rvs(size=MC_rep,loc=1.21*P_nom,scale=1.21*.0405*P_nom)

p2 = math.sqrt(math.log(((.063*1.03*sigmay_nom)**2/(1.03*sigmay_nom)**2)+1))
p1 = math.log(1.03*sigmay_nom) - (p2**2)/2
MC_sigmay = lognorm.rvs(s=p2,size=MC_rep,loc=0,scale=math.exp(p1))

# model_type -> regime físico do material empregado:
# 0 = linear;
# 1 = elastoplástico com encruamento isotrópico linear;
# 2 = elastoplástico com encruamento segundo equação de ramberg-osgood.   
model_type = 2
    
# model_dimension -> define se a treliça é plana ou espacial
# 2 = treliça bidimensional;
# 3 = treliça tridimensional.   
model_dimension = 2
    
nrep = P_nom # número de passos de carga aplicados

if (model_type == 0): # mecanismo de "segurança", pois nrep deve ser sempre 1
    nrep = 1          # quando o modelo for linear
    
# nodes -> coordenadas dos nós (em cm): 1a linha = coordenada x;
# 2a linha = coordenada y; 3a linha (se houver) = coordenada z.  
nodes = np.array([[0,200,400,200], 
                  [0,0,0,200]], 
                 dtype=np.float64)  
    
# con -> conectividade das barras: 1a linha = nó inicial; 2a linha = nó final
con = np.array([[1,2,1,2,4], 
                [2,3,4,4,3]], 
               dtype=np.float64)
    
# cc_type -> tipo de condição de contorno: 0 = força prescrita;
# 1 = deslocamento prescrito; 1a linha = coordenada x;
# 2a linha = coordenada y; 3a linha (se houver) = coordenada z.    
cc_type = np.array([[1,0,0,0], 
                    [1,0,1,0]], 
                   dtype=np.float64)

sigma_u = 40 # tensão de ruptura do material (kN/cm²)
epsilon_y = 0.0012 # deformação correspondente à tensão de escoamento
epsilon_u = 0.2 # deformação correspondente à tensão de ruptura
K = 6563 # módulo de encruamento isotrópico (kN/cm²)
# obs: verificar coerência entre E, sigma_y e epsilon_y
n_nodes = len(nodes[0]) # cálculo do número de nós da estrutura
n_elm = len(con[0]) # cálculo do número de barras da estrutura
    
L = np.zeros(shape=(n_elm,1)) # declaração do vetor dos comprimentos das barras      
    
# cálculo do comprimento de cada elemento (cm) 
for i in range (n_elm):      
    for t in range (model_dimension):
        L[i] += (nodes[t][int(con[1][i])-1]-nodes[t][int(con[0][i])-1])**2
    L[i] = math.sqrt(L[i])   

# simulação de monte carlo
completion = 0
failure = 0

g1 = np.zeros(shape=(n_elm,1)) # vetor de probabilidade da 1a funcao de falha
g2 = np.zeros(shape=(n_elm,1)) # vetor de probabilidade da 2a funcao de falha
g3 = np.zeros(shape=(n_elm,1)) # vetor de probabilidade da 3a funcao de falha
g4 = np.zeros(shape=(n_elm,1)) # vetor de probabilidade da 4a funcao de falha
Pf = np.zeros(shape=(n_elm,1)) # vetor de probabilidade de falha de cada barra

for y in range (MC_rep): 

# valores iniciais de parâmetros simulados
    A_i = math.pi*MC_t[y]*(2*MC_Diam[y]-MC_t[y])/4
    sigma_y = MC_sigmay[y]  # tensão de escoamento do material (kN/cm²)
    
# E -> módulo de rigidez de cada barra (em ordem). Unidade: kN/cm^2.     
    E = np.array([MC_E[y],MC_E[y],MC_E[y],MC_E[y],MC_E[y]],
                 dtype=np.float64)    
# A -> área da seção transversal de cada elemento (em ordem). Unidade: cm^2.   
    A = np.array([A_i,A_i,A_i,A_i,A_i],
                 dtype=np.float64)
      
# n -> constante de encruamento de ramberg-osgood
    n = np.log(epsilon_u/epsilon_y)/np.log(sigma_u/sigma_y)
           
# cc_value -> valores das condições de contorno prescritas (kN ou cm).
# OBS: em caso de regime não linear, atribuir o deltaP ao(s) nó(s) desejado(s). 
    cc_value = np.array([[0,0,0,0], 
                         [0,-MC_P[y]/nrep,0,0]], 
                        dtype=np.float64)  
    
    prop = np.zeros(shape=(n_elm,1)) # declaração do vetor das rigidezes E*A 
                                     # das barras 

# cálculo da rigidez E*A de cada elemento
    for i in range (n_elm):
        prop[i] = E[i]*A[i]

    sigma_acum, D_acum, delta_L_acum = nonlinear_solving (model_type, 
                                                          model_dimension, 
                                                          nodes, cc_type, 
                                                          cc_value, con, E, A, 
                                                          sigma_y, epsilon_y, 
                                                          K, n, nrep, n_nodes, 
                                                          n_elm, prop, L)

# verificação de critérios de falha para cálculo de probabilidade de falha 
# de cada elemento

# I -> momento de inércia da seção transversal da peca
    I = (math.pi/64)*((MC_Diam[y]**4)-((MC_Diam[y]-MC_t[y])**4))
    for i in range (n_elm):
        if (D_acum[3] < -1.6):
            failure = 1
            g1[i] += 100/MC_rep
        if (sigma_acum[i]*A[i] <= -((math.pi**2)*E[i]*I/(L[i]**2))):
            failure = 1
            g2[i] += 100/MC_rep
        if (abs(sigma_acum[i]) > sigma_u):
            failure = 1
            g3[i] += 100/MC_rep
        if (abs(sigma_acum[i]) > sigma_y):
            failure = 1
            g4[i] += 100/MC_rep
        if (failure == 1):
            Pf[i] += 100/MC_rep
        failure = 0        
    
    if ((y/100) == int(y/100)):
        print (f'Completion percentage: {100*y/MC_rep}%')

runtime = datetime.now()-runtime
print (f'Total runtime: {runtime}')