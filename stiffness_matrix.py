import numpy as np

def stiffness_matrix (L, prop, model_type, m_dim, n_nodes, n_elm, E, A, nodes, 
                      con, cc_type, cc_value):
    
    Rl = np.zeros(shape=(n_elm, 2*m_dim,2*m_dim)) # vetor de matrizes de rigidez de cada elemento no sistema global
    C = np.zeros(shape=(n_elm,3)) # matriz que guarda os cossenos diretores de cada barra
    R = np.zeros(shape=(m_dim*n_nodes,m_dim*n_nodes)) # matriz de rigidez da estrutura
    F = np.zeros(shape=(m_dim*n_nodes)) # vetor de forcas globais (kN)
    
# nmg -> utilizado para atribuir valores prescritos de deslocamento
    nmg = 100000000000000
    
    for i in range (n_elm):

# cálculo dos cossenos diretores de cada elemento        
        for t in range (m_dim):
            C[i][t] += (nodes[t][int(con[1][i])-1]-nodes[t][int(con[0][i])-1])/L[i]

# cálculo da matriz de rigidez de cada elemento no sistema global de coordenadas            
        if m_dim == 2:
            Rl[i] = np.array([
                [C[i][0]**2,C[i][0]*C[i][1],-C[i][0]**2,-C[i][0]*C[i][1]],
                [C[i][0]*C[i][1],C[i][1]**2,-C[i][0]*C[i][1],-C[i][1]**2],
                [-C[i][0]**2,-C[i][0]*C[i][1],C[i][0]**2,C[i][0]*C[i][1]],
                [-C[i][0]*C[i][1],-C[i][1]**2,C[i][0]*C[i][1],C[i][1]**2]],
                dtype=np.float64)
        else:
            Rl[i] = np.array([
                [C[i][0]**2,C[i][0]*C[i][1],C[i][0]*C[i][2],-C[i][0]**2,-C[i][0]*C[i][1],-C[i][0]*C[i][2]],
                [C[i][0]*C[i][1],C[i][1]**2,C[i][1]*C[i][2],-C[i][0]*C[i][1],-C[i][1]**2,-C[i][1]*C[i][2]],
                [C[i][0]*C[i][2],C[i][1]*C[i][2],C[i][2]**2,-C[i][0]*C[i][2],-C[i][1]*C[i][2],-C[i][2]**2],
                [-C[i][0]**2,-C[i][0]*C[i][1],-C[i][0]*C[i][2],C[i][0]**2,C[i][0]*C[i][1],C[i][0]*C[i][2]],
                [-C[i][0]*C[i][1],-C[i][1]**2,-C[i][1]*C[i][2],C[i][0]*C[i][1],C[i][1]**2,C[i][1]*C[i][2]],
                [-C[i][0]*C[i][2],-C[i][1]*C[i][2],-C[i][2]**2,C[i][0]*C[i][2],C[i][1]*C[i][2],C[i][2]**2]],
                dtype=np.float64)

        Rl[i] = (prop[i]/L[i])*np.copy(Rl[i])
        
# variáveis auxiliares para alocação das submatrizes na matriz de rigidez
        nodes1 = int(m_dim*(con[0][i]-1))
        nodes2 = int(m_dim*con[0][i])
        nodes3 = int(m_dim*(con[1][i]-1))
        nodes4 = int(m_dim*con[1][i])

# montagem da matriz de rigidez global da estrutura    
        R[nodes1:nodes2,nodes1:nodes2] += Rl[i][0:m_dim,0:m_dim]
        R[nodes1:nodes2,nodes3:nodes4] += Rl[i][0:m_dim,m_dim:2*m_dim]
        R[nodes3:nodes4,nodes1:nodes2] += Rl[i][m_dim:2*m_dim,0:m_dim]
        R[nodes3:nodes4,nodes3:nodes4] += Rl[i][m_dim:2*m_dim,m_dim:2*m_dim]

# montagem do vetor de acoes globais        
        for i in range (n_nodes):
            for j in range (m_dim):
                if (cc_type[j][i]==0):
                    F[int(m_dim*i+j)] = cc_value[j][i]
                else:
                    glib = m_dim*i+j
                    F[glib] = nmg*cc_value[j][i]

    return R, Rl, C, F 