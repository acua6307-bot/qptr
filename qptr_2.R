h <- 6.62607015e-34         # J*s
hbar <- h/(2*pi)         # J*s
kB <- 1.380649e-23          # J*s
c_light <- 299792458.0      # m/s
N <- 2                  # Fock-space dimension (truncation)
nu <- 10e14              # frequency in Hz (example: THz)
omega <- 2*pi*nu        # angular frequency
Temp <- 300.0           # temperature in K  
# --- QuTiP operators ---
destroy <- function(N) {
  # Creamos una matriz de N x N llena de ceros
  mat <- matrix(0, nrow = N, ncol = N)
  # Llenamos la subdiagonal superior con sqrt(n)
  # n va de 1 hasta N-1
  for (i in 1:(N - 1)) {
    mat[i, i + 1] <- sqrt(i)
  }
  # Retornamos el objeto siguiendo la estructura de tu framework
  # Asumimos que Operator() es el constructor que ya definiste
  return(Operator(mat))
}
a <- destroy(N)
adag <- (a|>adjoint())
n_op <- adag[[a]]
Ham <- Tensor(hbar)*Tensor(omega)*(n_op + Tensor(0.5)*Identity(6))
# --- Thermal state ---
nbar <- 1/(exp(hbar*omega/(kB*Temp)) - 1.0)
thermal_dm <- function(N, nbar) {
  r <- nbar / (nbar + 1)
  probs <- (1 - r) * (r^(0:(N - 1)))
  # Renormalización: dividir por la suma de las probabilidades
  probs_norm <- probs / sum(probs)
  mat <- diag(probs_norm, nrow = N, ncol = N)
  return(Operator(mat))
}
rho_T <- thermal_dm(N,nbar)
rho_T[[n_op]]$data|>diag()|>sum()
rho_T[[Ham]]$data|>diag()|>sum()

T_test <- 300.0
nus <- 10^seq(10, 14, length.out = 80) #1e10 a 1e14 Hz
E_th <- numeric(length(nus))
E_class <- kB*T_test
for (i in seq_along(nus)){
  nu_i <- nus[i]
  omega_i = 2*pi*nu_i
  nbar_i = 1.0/(exp(hbar*omega_i/(kB*T_test)) - 1.0)
  E_th[i] <- hbar*omega_i*nbar_i
}

T_spec <- 3000.0 # K
nu_min <- 1e11
nu_max <- 3e15
nu_grid <- nu_grid <- 10^seq(log10(nu_min), log10(nu_max), length.out = 600)
x <- h*nu_grid/(kB*T_spec)
u_planck = (nu_grid**2) * (h*nu_grid) / (exp(x) - 1.0)
u_rj = (nu_grid**2) * (kB*T_spec)
# Normalizacion (solo para comparar formas)
u_planck_n = u_planck/max(u_planck)
u_rj_n = u_rj/max(u_planck)
