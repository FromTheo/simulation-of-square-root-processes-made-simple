import numpy as np

from numba import njit, prange 

@njit(parallel=True, fastmath=True)
def simulate_hat_VUZ(V0, T, a, b, c, N, eta) : 
    M, n = eta.shape 
    dt = T/n 
    hat_V, hat_U, hat_Z = np.empty((M, n+1)), np.empty((M, n)), np.empty((M, n))

    term1 = (np.exp(b*dt) - 1) / b 
    term2 = a/b * (term1 - dt) 
    c_term1_inv = 1 / (c * term1) 
    c_term1_inv_sq = 1 / (c * term1)**2 
    for m in prange(M) : 
        hat_V[m, 0] = V0 
        for i in range(n) : 
            nmi = N[m, i] 
            y = nmi * nmi 
            alpha_i = hat_V[m, i] * term1 + term2 
            mu = alpha_i 
            lambd = alpha_i**2 * c_term1_inv_sq
            X = mu + 0.5 * mu**2 * y / lambd - 0.5 * mu / lambd * np.sqrt(4 * mu * lambd * y + mu**2 * y**2) 
            hat_U[m, i] = X if eta[m, i] <= mu / (mu + X) else mu**2 / X
            hat_Z[m, i] = c_term1_inv * (hat_U[m, i] - alpha_i)
            hat_V[m, i+1] = hat_V[m, i] + a*dt + b*hat_U[m, i] + c*hat_Z[m, i] 
    
    return hat_V, hat_U, hat_Z 

@njit(parallel=True, fastmath=True)
def simulate_hat_S(hat_U, hat_Z, S0, rho, N_) :         
    M, n = hat_U.shape
    logS = np.empty((M, n+1))
    logS0 = np.log(S0)
    sqrt_rho = np.sqrt(1-rho**2)
    for m in prange(M) : 
        logS[m, 0] = logS0 
        for i in range(n) : 
            logS[m, i+1] = logS[m, i] - 0.5 * hat_U[m, i] + rho * hat_Z[m, i] + sqrt_rho * np.sqrt(hat_U[m, i]) * N_[m, i]
    return np.exp(logS)

@njit(parallel=True) #faster than np.cumsum 
def get_UZ_trajectories(hat_U, hat_Z) :  
    M, n = hat_U.shape 
    U, Z = np.empty((M, n+1)), np.empty((M, n+1))
    for m in prange(M):
        U[m, 0] = 0.0
        Z[m, 0] = 0.0 
        for i in range(n):
            U[m, i+1] = U[m, i] + hat_U[m, i]
            Z[m, i+1] = Z[m, i] + hat_Z[m, i]
    return U, Z


