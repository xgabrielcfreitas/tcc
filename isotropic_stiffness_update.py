import numpy as np

def isotropic_stiffness_update (n_elm, n_nodes, m_dim, con, L, K, E, E0, 
                                sigma_acum, sigma, A, prop, sigma_y, Rl, R,
                                updated):
    
# verificação de encruamento e atualização da matriz de rigidez do elemento
    for i in range (n_elm):
        if (E[i] == E0[i]):
            if (abs(sigma_acum[i]+sigma[i]) > abs(sigma_y)):
                Rl[i] = (L[i]/prop[i])*np.copy(Rl[i])
                E[i] = E[i]*K/(E[i]+K)
                prop[i] = E[i]*A[i]                
                Rl[i] = (prop[i]/L[i])*np.copy(Rl[i])
                updated = 1

# remontagem da matriz de rigidez da estrutura        
    if (updated == 1):
        R = np.zeros(shape=(m_dim*n_nodes,m_dim*n_nodes))
        for i in range (n_elm):
            
# variáveis auxiliares para alocação das submatrizes na matriz de rigidez
            nodes1 = int(m_dim*(con[0][i]-1))
            nodes2 = int(m_dim*con[0][i])
            nodes3 = int(m_dim*(con[1][i]-1))
            nodes4 = int(m_dim*con[1][i])

# alocação das submatrizes na matriz de rigidez global da estrutura    
            R[nodes1:nodes2,nodes1:nodes2] += Rl[i][0:m_dim,0:m_dim]
            R[nodes1:nodes2,nodes3:nodes4] += Rl[i][0:m_dim,m_dim:2*m_dim]
            R[nodes3:nodes4,nodes1:nodes2] += Rl[i][m_dim:2*m_dim,0:m_dim]
            R[nodes3:nodes4,nodes3:nodes4] += Rl[i][m_dim:2*m_dim,m_dim:2*m_dim]

    return updated, E, prop, Rl, R