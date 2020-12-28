import numpy as np

def RO_stiffness_update (n_elm, n_nodes, m_dim, con, L, E0, sigma_acum, sigma, 
                         A, prop, sigma_y, Rl, n, epsilon_y):
    
# limpeza para posterior remontagem da matriz de rigidez da estrutura     
    R = np.zeros(shape=(m_dim*n_nodes,m_dim*n_nodes))
    
# preparo de variável para receber os modulos de rigidez dos elementos     
    E = np.zeros(shape=(n_elm))
    
# atualização da matriz de rigidez do elemento
    for i in range (n_elm):
        Rl[i] = (L[i]/prop[i])*np.copy(Rl[i])
# atualizacao do modulo de elasticidade utilizado em cada barra como sendo
# o modulo de elasticidade tangente no ponto atual de carregamento
# obs: esse metodo requer acrescimentos de carregamento suficientemente
# pequenos para nao afetar demasiadamente a qualidade dos resultados obtidos
        E[i] = np.copy(1/((1/E0[i])+epsilon_y*(n/sigma_y)*
                          ((abs(sigma_acum[i]/sigma_y))**(n-1))))
        prop[i] = E[i]*A[i]                
        Rl[i] = (prop[i]/L[i])*np.copy(Rl[i])
            
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
    
    return E, prop, Rl, R