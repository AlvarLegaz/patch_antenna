from math import acos, asin, cos, exp, log10, pi, sin, sqrt
from scipy.constants import c, pi
from scipy.special import sici
from scipy.integrate import quad
from scipy.special import j0

# === Parámetros de entrada: sustrato FR4 ===
Er = 4.2               # Permitividad relativa del sustrato
h = 1.6             # Altura del sustrato (mm) 
f0 = 1.41e9            # Frecuencia de resonancia (Hz)
c = 3e8                # Velocidad de la luz (m/s)
Z0 = 50              # Impedancia característica deseada linea alimentación(ohmios)

lambda_0 = (c / f0) * 1000  # Longitud de onda en el vacío (mm)


# === Cálculo del ancho para parche rectangular ===
W_rect = (c / (2 * f0)) * sqrt(2 / (Er + 1)) * 1000  # Ancho óptimo (mm) para parche rectangular

# === Cálculo refinado de permitividad efectiva con W definido ===
Ereff = (Er + 1) / 2 + ((Er - 1) / 2) * (1 + 12 * (h / (W_rect))) ** (-0.5)


# Cálculo de ΔL
delta_L_over_h = 0.412 * ((Ereff + 0.3) * (W_rect/h + 0.264)) / ((Ereff - 0.258) * (W_rect/h + 0.8))
delta_L = delta_L_over_h * h   # en mm

# === Recalculo del largo físico del parche con Delta L ===
L_eff = lambda_0 / (2 * sqrt(Ereff))  # Largo eléctrico efectivo actualizado
L = L_eff - 2 * delta_L           # Largo físico (mm)
W = W_rect                            # Antena rectangular

# === Resultados ===
print(f"\n--- Dimensiones parche---")

print(f"Parche rectangular resonante a {f0/1e9:.2f} GHz")

print(f"Permitividad efectiva (Ereff): {Ereff:.4f}")
print(f"Delta Lado del parche (L): {delta_L:.4f} mm")

print(f"Ancho del parche (W): {W:.4f} mm")
print(f"Lado del parche (L): {L:.4f} mm")


# === Cálculo impedancia entrada ===

# pasomos todo a metros
W= W / 1000
L= L / 1000

lambda_0 = c / f0       # Longitud de onda en el vacío (m)
k0 = 2 * pi / lambda_0  # Número de onda
X = k0 * (W)              # Parámetro X = k0 * W(m)

# === Aproximación de I1 usando la función seno integral ===
Si, _ = sici(X)  # scipy devuelve (Si(X), Ci(X))
I1 = -2 + cos(X) + X * Si + sin(X)/X

# === Cálculo de G1 ===
G1 = I1 / (120 * pi**2) 

# Función integrando para G12 (Balanis eq. 14-18.a)
def integrando(theta):
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    if abs(cos_theta) < 1e-6:
        return 0  # Evita división por cero
    argument = (k0 * W / 2) * cos_theta
    slot_term = (sin(argument) / cos_theta) ** 2
    bessel_term = j0(k0 * L * sin_theta)
    return slot_term * bessel_term * sin_theta**3

# Integración numérica
I12, _ = quad(integrando, 0, pi)
G12 = I12 / (120 * pi**2)

# === Impedancia de entrada estimada  ===
Zin_0 = 1 / (2 * (G1 + G12))

Zin = Z0 

y_0 = (L/pi) * acos(sqrt(Zin/Zin_0))/1000


print(f"\n--- Impedancia de entrada  parche---")
print(f"X = k₀W: {X:.4f}")
print(f"I₁ (aprox): {I1:.4f}")
print(f"Conductancia G₁: {G1:.4e} S")
print(f"Conductancia G12: {G12:.4e} S")
print(f"Impedancia de entrada Zₐ (borde): {Zin_0:.2f} Ω")
print(f"Recesion de alimentación y0 para Zin={Zin:.2f} Ω y0 = {y_0:.4f} mm")


# === Cálculo ancho linea alimentacion 50 ohm ===

def microstrip_width(Z0, Er, h):
    """Calcula el ancho W (en mm) de una línea microstrip sobre un sustrato dado."""
    A = Z0 / 60 * sqrt((Er + 1) / 2) + ((Er - 1) / (Er + 1)) * (0.23 + 0.11 / Er)
    W_h = (8 * exp(A)) / (exp(2 * A) - 2)
    
    # Para casos W/h > 2
    if W_h > 2:
        B = (377 * pi) / (2 * Z0 * sqrt(Er))
        W_h = (2 / pi) * (
            B - 1 - log(2 * B - 1) +
            ((Er - 1) / (2 * Er)) * (log(B - 1) + 0.39 - 0.61 / Er)
        )
    
    return W_h * h  # devuelve W en mm

W = microstrip_width(Z0, Er, h)
print(f"✅ Ancho de línea microstrip para 50 Ω sobre FR-4 (h = {h} mm): {W:.2f} mm")

