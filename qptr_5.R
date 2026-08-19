h <- 6.62607015e-34      # J s
c_e <- 2.99792458e8        # m/s
eC <- 1.602176634e-19     # C
m_e <- 9.1093837015e-31  # kg

lambda_c <- h / (m_e * c)
lambda_c

n_max <- 8
Z <- 1
n <- seq(1, n_max)
En_eV = -13.6*Z**2/n**2
En_J = En_eV * eC
print(paste("E_1 (eV) =", En_eV[1]))

library(tidyverse)
ggplot(data.frame(n = n, En_eV = En_eV), aes(x = n, y = En_eV)) + 
  geom_point() + 
  geom_line()

H_energy <- Operator(diag(En_J))
H_proj_energy <- Tensor(0) * H_energy
for(i in seq(En_J)){
  ket<-bitstring.numeric(rev(as.integer(intToBits(i-1)[1:(log(n_max,2))])));
  H_proj_energy <- H_proj_energy+Tensor(En_J[i])*ket[[adjoint(Tensor(1)*ket)]]
}

norm((H_energy-H_proj_energy)$data)
H_energy$data|>diag()
diag(H_energy$data/eC)
En_eV

rydberg_wavelength <- function(nf, ni, Z = 1, R = R_inf) {
  if (any(nf < 1) || any(ni < 1) || any(Z < 1)) {
    stop("nf, ni, and Z must be positive integers.")
  }
  if (any(ni <= nf)) {
    stop("Require ni > nf for a transition (ni -> nf).")
  }
  inv_lambda <- R * (Z^2) * (1 / nf^2 - 1 / ni^2)
  1 / inv_lambda
}

R_inf <- 10973731.568160 # m^-1
for (ni in 3:10) {
  lam_nm <- rydberg_wavelength(
    nf = 2,
    ni = ni,
    Z  = 1,
    R  = R_inf
  ) * 1e9
  cat(sprintf("n=%d -> 2: %.3f nm\n", ni, lam_nm))
}

# stick_spectrum <- function(nf, nmax, Z = 1, R = R_inf) {
#   lambda_nm <- numeric()
#   for (ni in (nf + 1):nmax) {
#     lambda_m <- rydberg_wavelength(
#       nf = nf,
#       ni = ni,
#       Z  = Z,
#       R  = R
#     )
#     lambda_nm <- c(
#       lambda_nm,
#       lambda_m * 1e9
#     )
#   }
#   lambda_nm
# }

stick_spectrum <- function(nf, nmax, Z = 1, R = R_inf) {
  ni <- (nf + 1):nmax
  rydberg_wavelength(
    nf = nf,
    ni = ni,
    Z = Z,
    R = R
  ) * 1e9
}

lam_nm <- stick_spectrum(
  nf = 2,
  nmax = n_max,
  Z = 1,
  R = R_inf
)

spec_df <- data.frame(
  wavelength = lam_nm,
  ymin = 0,
  ymax = 1
)

ggplot(data.frame(wavelength = lam_nm)) +
  geom_linerange(
    aes(
      x = wavelength,
      ymin = 0,
      ymax = 1
    )
  ) +
  labs(
    x = "Wavelength (nm)",
    y = "Intensity (arb.)",
    title = "Balmer series (stick spectrum), Z = 1"
  ) +
  theme_minimal()



spec_df <- tibble(
  Z = c(1, 2, 3),
  y0 = 1.2 * (0:2)
) %>%
  mutate(
    wavelength = map(
      Z,
      ~stick_spectrum(
        nf = 2,
        nmax = n_max,
        Z = .x,
        R = R_inf
      )
    )
  ) %>%
  unnest_longer(wavelength)

ggplot(
  spec_df,
  aes(
    x = wavelength,
    color = factor(Z)
  )
) +
  geom_segment(
    aes(
      xend = wavelength,
      y = y0,
      yend = y0 + 1
    ),
    linewidth = 0.7
  ) +
  labs(
    title = "Balmer series shift with Z (stick spectra)",
    x = "Wavelength (nm)",
    y = "Relative intensity (arb.)",
    color = "Z"
  ) +
  theme_minimal()


m_p <- 1.6726219259552e-27
M <- m_p
mu <- (m_e*M)/(m_e+M)
R_corr <- R_inf * (mu/m_e)
print(paste("mu/me =", mu/m_e))
print(paste("R_corr (1/m) =", R_corr))



ni_vals <- 3:n_max
lam_inf <- rydberg_wavelength(
  nf = 2,
  ni = ni_vals,
  Z = Z,
  R = R_inf
)
lam_corr <- rydberg_wavelength(
  nf = 2,
  ni = ni_vals,
  Z = Z,
  R = R_corr
)
delta_pm <- (lam_corr - lam_inf) * 1e12
for (i in seq_along(ni_vals)) {
  
  cat(
    sprintf(
      "n=%d -> 2: lambda_inf=%.6f nm, lambda_corr=%.6f nm, delta_lambda=%.3f pm\n",
      ni_vals[i],
      lam_inf[i] * 1e9,
      lam_corr[i] * 1e9,
      delta_pm[i]
    )
  )
}
df <- data.frame(
  ni = ni_vals,
  delta_pm = delta_pm
)
ggplot(df, aes(x = ni, y = delta_pm)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Reduced-mass correction for Balmer wavelengths (H)",
    x = expression("Initial level " * n[i]),
    y = expression(Delta * lambda ~ "(pm)")
  ) +
  theme_minimal()


n0 <- 5
rho0 <- (bitstring.numeric(rev(as.integer(intToBits(n0-1)[1:(log(n_max,2))]))))|>density()
gamma0 <- 1.0e8
c_ops <- list()
for (n_level in seq(2,n_max)){
  gamma_n <- gamma0 / (n_level**3)
  ket_down <- bitstring.numeric(rev(as.integer(intToBits(n_level-2)[1:(log(n_max,2))])))
  bra_up <- (Tensor(1)*bitstring.numeric(rev(as.integer(intToBits(n_level-1)[1:(log(n_max,2))]))))|>adjoint()
  c_ops[[n_level-1]] <- (Tensor(sqrt(gamma_n))*(ket_down*bra_up))
}

