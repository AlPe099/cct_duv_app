import streamlit as st
import math

def planckian_xy(T):
    """Approximation polynomiale pour (x, y) du corps noir (4000..25000 K)."""
    invT = 1.0 / T
    invT2 = invT * invT
    invT3 = invT2 * invT

    x0 = -3.0258469e9 * invT3 + 2.1070379e6 * invT2 + 0.2226347e3 * invT + 0.240390
    y0 = -3.0 * (x0**2) + 2.87 * x0 - 0.275
    return (x0, y0)

def xy_to_uv(x, y):
    """Conversion (x, y) CIE 1931 -> (u', v') CIE 1976."""
    if y == 0:
        return (0.0, 0.0)
    X = x / y
    Z = (1.0 - x - y) / y
    denom = X + 15.0 + 3.0 * Z
    if denom == 0:
        return (0.0, 0.0)
    u = (4.0 * X) / denom
    v = (9.0 * 1.0) / denom
    return (u, v)

def uv_to_xy(u, v):
    """Conversion (u', v') CIE 1976 -> (x, y) CIE 1931."""
    denom = (6.0*u) - (16.0*v) + 12.0
    if denom == 0:
        return (0.0, 0.0)
    x = (9.0 * u) / denom
    y = (4.0 * v) / denom
    return (x, y)

def cct_duv_to_xy(T, Duv):
    """
    Calcule (x, y) à partir de T (Kelvin) et Duv.
    Méthode : on prend (x0, y0) à T et (x1, y1) à T+0.01 -> tangente -> applique Duv.
    """
    x0, y0 = planckian_xy(T)
    x1, y1 = planckian_xy(T + 0.01)

    u0, v0 = xy_to_uv(x0, y0)
    u1, v1 = xy_to_uv(x1, y1)

    du_tang = u1 - u0
    dv_tang = v1 - v0
    length_tang = math.sqrt(du_tang**2 + dv_tang**2)
    if length_tang < 1e-15:
        return (x0, y0)

    nx = -dv_tang / length_tang
    ny =  du_tang / length_tang

    u = u0 + Duv * nx
    v = v0 + Duv * ny

    return uv_to_xy(u, v)

# --- INTERFACE STREAMLIT ---
def main():
    st.title("CCT + Duv → (x, y)")

    T = st.number_input("Température de couleur (K)", min_value=1000.0, max_value=30000.0, value=6500.0)
    Duv = st.number_input("Duv", min_value=-0.1, max_value=0.1, value=0.0, format="%.5f")

    if st.button("Convertir"):
        x, y = cct_duv_to_xy(T, Duv)
        st.success(f"Résultat : x = {x:.4f}, y = {y:.4f}")

if __name__ == "__main__":
    main()
