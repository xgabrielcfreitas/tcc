import numpy as np

def direct_stiffness_solving (m_dim, n_nodes, n_elm, prop, L, con, 
                              cc_type, cc_value, R, C, A, F):

    D = [] # vetor de deslocamentos nodais globais (cm)
    delta_L = np.zeros(shape=(n_elm,1)) # vetor de deslocamentos nodais nos 
                                        # eixos dos elementos (cm)
    N = np.zeros(shape=(n_elm,1)) # vetor de esforços axiais (kN)
    sigma = np.zeros(shape=(n_elm,1)) # vetor de tensões axiais (kN/cm²)
    
# dl -> vetor de deformações axiais nas barras
    dl = np.zeros(shape=(n_elm,2*m_dim))

# nmg -> utilizado para atribuir valores prescritos de deslocamento
    nmg = 100000000000000 

# montagem do vetor de acoes globais com artifício numérico para atribuir
# valores prescritos de deslocamento
    for i in range (n_nodes):
        for j in range (m_dim):
            if (cc_type[j][i]==1):
                glib = m_dim*i+j
                R[glib][glib] = nmg

# cálculo dos deslocamentos nodais globais
    D = np.linalg.inv(R)@F
    
# extração dos deslocamentos nodais nos eixos das barras, cálculo das
# deformações e esforços axiais dos elementos
    for i in range (n_elm):
        dl[i][0:m_dim] = D[int(m_dim*(con[0][i]-1)):int(m_dim*con[0][i])]
        dl[i][m_dim:2*m_dim] = D[int(m_dim*(con[1][i]-1)):int(m_dim*con[1][i])]        
        for j in range (m_dim):
            delta_L[i] += np.copy(C[i][j]*(dl[i][m_dim+j]-dl[i][j]))
        N[i] = (prop[i]/L[i])*delta_L[i]
        sigma[i] = N[i]/A[i]
           
    return D, delta_L, N, sigma