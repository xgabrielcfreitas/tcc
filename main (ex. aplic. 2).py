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
from nonlinear_solving import nonlinear_solving

# parâmetros nominais de material e geometria dos elementos
Diam = 13.00 # Diâmetro nominal do tubo circular vazado. Unidade: cm
t = .63 # Espessura do tubo. Unidade: cm
E_nom = 21000 # Módulo de elasticidade do material. Unidade: kN/cm^2
sigma_y = 36 # Tensão de escoamento do material. Unidade: kN/cm^2
P = 2 # Carregamento imposto (força ou deslocamento). Unidade: kN ou cm

# parâmetros da simulação de monte carlo
MC_rep = 1 # número de repetições da simulação de monte carlo

# model_type -> regime físico do material empregado:
# 0 = linear;
# 1 = elastoplástico com encruamento isotrópico linear;
# 2 = elastoplástico com encruamento segundo equação de ramberg-osgood.   
model_type = 2
    
# model_dimension -> define se a treliça é plana ou espacial
# 2 = treliça bidimensional;
# 3 = treliça tridimensional.   
model_dimension = 3
    
nrep = 1000 # número de passos de carga aplicados

if (model_type == 0): # mecanismo de "segurança", pois nrep deve ser sempre 1
    nrep = 1          # quando o modelo for linear
    
# nodes -> coordenadas dos nós (em cm): 1a linha = coordenada x;
# 2a linha = coordenada y; 3a linha (se houver) = coordenada z.  

nodes = np.zeros(shape=(3,221))

# COORDENADAS X e Y
for i in range (11):
    for j in range (11):
        nodes[0][11*i+j] = 300*j
        nodes[1][11*j+i] = 300*j
for i in range (10):
    for j in range (10):
        nodes[0][121+10*i+j] = 150+300*j
        nodes[1][121+10*j+i] = 150+300*j        
#COORDENADAS Z
for i in range (100):
    nodes[2][121+i]=300
    
# con -> conectividade das barras: 1a linha = nó inicial; 2a linha = nó final

con = np.zeros(shape=(2,800))
# BARRAS EIXO X PLANO INFERIOR
aux = 0
for i in range (110):
    if (int(i/10) == i/10):
        aux += 1
    con[0][i] = i+aux
    con[1][i] = i+aux+1    
# BARRAS EIXO Y PLANO INFERIOR
for i in range (11):
    for j in range (11):
        con[0][110+10*i+j] = 11*j+i+1
        con[1][110+10*i+j] = 11*j+i+12
    con[0][220] = 0
    con[1][220] = 0
# BARRAS EIXO X PLANO SUPERIOR
aux = 0
for i in range (90):
    if (int(i/9) == i/9):
        aux += 1
    con[0][220+i] = 121+i+aux
    con[1][220+i] = 121+i+aux+1   
# BARRAS EIXO Y PLANO SUPERIOR
for i in range (10):
    for j in range (9):
        con[0][310+9*i+j] = 10*j+i+122
        con[1][310+9*i+j] = 10*j+i+132
    con[0][400] = 0
    con[1][400] = 0
# BARRAS DIAGONAIS 1 fileira
for i in range (400):
    con[0][400+i] = 122+int(i/4)
aux = 0
for i in range (100):
    if (int(i/10)==i/10):
        aux += 1
    con[1][400+4*i] = i + aux
    con[1][401+4*i] = i+1 + aux
    con[1][402+4*i] = i+11 + aux
    con[1][403+4*i] = i+12 + aux
   
# cc_type -> tipo de condição de contorno: 0 = força prescrita;
# 1 = deslocamento prescrito; 1a linha = coordenada x;
# 2a linha = coordenada y; 3a linha (se houver) = coordenada z.    
cc_type = np.zeros(shape=(3,221))

for i in range (3):
    for j in range (11):
        cc_type[i][j] = 1
        cc_type[i][j+110] = 1
    
sigma_u = 52 # tensão de ruptura do material (kN/cm²)
epsilon_y = sigma_y/E_nom # deformação correspondente à tensão de escoamento
epsilon_u = 0.2 # deformação correspondente à tensão de ruptura
K = 80.9 # módulo de encruamento isotrópico (kN/cm²)
# obs: verificar coerência entre E, sigma_y e epsilon_y
n_nodes = len(nodes[0]) # cálculo do número de nós da estrutura
n_elm = len(con[0]) # cálculo do número de barras da estrutura
    
L = np.zeros(shape=(n_elm,1)) # declaração do vetor dos comprimentos das barras      
    
# cálculo do comprimento de cada elemento (cm) 
for i in range (n_elm):      
    for r in range (model_dimension):
        L[i] += (nodes[r][int(con[1][i])-1]-nodes[r][int(con[0][i])-1])**2
    L[i] = math.sqrt(L[i])   

# simulação de monte carlo
for y in range (MC_rep): 

# valores iniciais de parâmetros simulados
    A_nom = math.pi*t*(2*Diam-t)/4
    
# E -> módulo de rigidez de cada barra (em ordem). Unidade: kN/cm^2.     
    E = np.ones(shape=(n_elm,1))
    E = E*E_nom
 
# A -> área da seção transversal de cada elemento (em ordem). Unidade: cm^2.   
    A = np.ones(shape=(n_elm,1))
    A = A*A_nom

# n -> constante de encruamento de ramberg-osgood
    n = np.log(epsilon_u/epsilon_y)/np.log(sigma_u/sigma_y)
           
# cc_value -> valores das condições de contorno prescritas (kN ou cm).
# OBS: em caso de regime não linear, atribuir o deltaP ao(s) nó(s) desejado(s). 
    cc_value = np.zeros(shape=(3,221))

    for i in range (100):
        cc_value[2][121+i] = -P/nrep

# cálculo da rigidez E*A de cada elemento
    prop = np.ones(shape=(n_elm,1))
    prop = prop*A_nom*E_nom

    sigma_acum, D_acum, delta_L_acum, E = nonlinear_solving (model_type, 
                                                          model_dimension, 
                                                          nodes, cc_type, 
                                                          cc_value, con, E, A, 
                                                          sigma_y, epsilon_y, 
                                                          K, n, nrep, n_nodes, 
                                                          n_elm, prop, L)

# verificação de critérios de falha para cálculo de probabilidade de falha 
# de cada elemento

# I -> momento de inércia da seção transversal da peca
    I = (math.pi/64)*((Diam**4)-((Diam-t)**4))
    for i in range (3*n_nodes):
        if (D_acum[i] < -12):
            print (f'deslocou demais na posição i = {i}')
    for i in range (n_elm):
        if (sigma_acum[i]*A[i] <= -((math.pi**2)*E_nom*I/(L[i]**2))):
            print (f'flambou na posição i = {i}')
        if (abs(sigma_acum[i]) > sigma_u):
            print (f'rompeu na posição i = {i}')
        if (abs(sigma_acum[i]) > sigma_y):
            print(f'escoou na posicao i = {i}')

maiordesloc = D_acum[3*n_nodes-1]
maiortensao = sigma_acum[n_elm-1]
menortensao = sigma_acum[n_elm-1]
for i in range (n_nodes*3-1):
    if (abs(D_acum[i])>abs(maiordesloc)):
        maiordesloc = D_acum[i]
for i in range (n_elm-1):
    if (sigma_acum[i]>maiortensao):
        maiortensao = sigma_acum[i]
    if (sigma_acum[i]<menortensao):
        menortensao = sigma_acum[i]

print (f'maior deslocamento: {maiordesloc}')
print (f'maior tensao: {maiortensao}')
print (f'menor tensao: {menortensao}')

runtime = datetime.now()-runtime
print (f'Total runtime: {runtime}')
