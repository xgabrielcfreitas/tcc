# -*- coding: utf-8 -*-
"""
O codigo abaixo visa a resolucao de trelicas planas e espaciais
de acordo com diferentes regimes de elastoplasticidade.
Data: 27/11/2020
@autor: Gabriel de Carvalho Freitas
"""

import numpy as np
from stiffness_matrix import stiffness_matrix
from direct_stiffness_solving import direct_stiffness_solving
from isotropic_stiffness_update import isotropic_stiffness_update
from RO_stiffness_update import RO_stiffness_update

def nonlinear_solving (model_type, model_dimension, nodes, cc_type, cc_value, 
                       con, E, A, sigma_y, epsilon_y, K, n, nrep, n_nodes, 
                       n_elm, prop, L):

# montagem do vetor de acoes e da matriz de rigidez da estrutura
    R, Rl, C, F = stiffness_matrix (L, prop, model_type, model_dimension, 
                                    n_nodes, n_elm, E, A, nodes, con, cc_type, 
                                    cc_value)

# resolucao da estrutura
    D, delta_L, N, sigma = direct_stiffness_solving (model_dimension, n_nodes, 
                                                     n_elm, prop, L, con, 
                                                     cc_type, cc_value, R, C, 
                                                     A, F)

# parâmetros de registro cumulativo de tensoes, deslocamentos e deformacoes
# para os modelos nao lineares
    sigma_acum = np.copy(sigma)
    D_acum = np.copy(D)
    delta_L_acum = np.copy(delta_L)
    E0 = np.copy(E)

# analise nao linear isotropica
    if (model_type == 1):
        for y in range (nrep-1):
            updated = 0 # verificador de encruamento de elementos

# verificacao de encruamento e atualizacao da rigidez dos elementos
            updated, E, prop, Rl, R = isotropic_stiffness_update (n_elm, 
                                                                  n_nodes, 
                                                                  model_dimension, 
                                                                  con, L, K, E, 
                                                                  E0, sigma_acum, 
                                                                  sigma, A, prop, 
                                                                  sigma_y, Rl, R, 
                                                                  updated)
        
# calculo dos deslocamentos e tensoes com a nova matriz de rigidez
            if (updated == 1):
                D, delta_L, N, sigma = direct_stiffness_solving (model_dimension, 
                                                                 n_nodes, n_elm, 
                                                                 prop, L, con, 
                                                                 cc_type, cc_value, 
                                                                 R, C, A, F)
            
# armazenamento das tensoes acumuladas para a proxima iteracao, assim como
# as variaveis de deslocamento e deformacao para fins de verificacao
            sigma_acum += np.copy(sigma)
            D_acum += np.copy(D)
            delta_L_acum += np.copy(delta_L)

# analise nao linear segundo a equacao de ramberg-osgood
    elif (model_type == 2):
        for y in range (nrep-1):
        
# atualizacao da rigidez dos elementos
            E, prop, Rl, R = RO_stiffness_update (n_elm, n_nodes, model_dimension, 
                                                  con, L, E0, sigma_acum, sigma, A, 
                                                  prop, sigma_y, Rl, n, epsilon_y)
        
# calculo dos deslocamentos e tensoes com a nova matriz de rigidez
            D, delta_L, N, sigma = direct_stiffness_solving (model_dimension, 
                                                             n_nodes, n_elm, prop, 
                                                             L, con, cc_type, 
                                                             cc_value, R, C, A, F)
        
# armazenamento das tensoes acumuladas para a proxima iteracao, assim como
# as variaveis de deslocamento e deformacao para fins de verificacao
            sigma_acum += np.copy(sigma)
            D_acum += np.copy(D)
            delta_L_acum += np.copy(delta_L)
            
    return (sigma_acum, D_acum, delta_L_acum)