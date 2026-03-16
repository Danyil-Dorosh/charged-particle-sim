"""
Charged particle motion in a static magnetic field.

Solves the Lorentz-force equation  m dv/dt = q (v × B)
using SciPy's RK4(5) integrator (solve_ivp with method='RK45').
"""

import numpy as np
from scipy.integrate import solve_ivp

# --------------- physical parameters ---------------
q = 1.602e-19      
m = 1.673e-27      
B0 = 0.1          
B_dir = np.array([0.0, 0.0, 1.0]) 

# --------------- initial conditions -----------------
v_perp = 1e5       
v_par  = 3e4        

r0 = np.array([0.0, 0.0, 0.0])                  
v0 = np.array([v_perp, 0.0, v_par])             
y0 = np.concatenate([r0, v0])                     

# --------------- time span --------------------------
omega_c = abs(q) * B0 / m         
T_c     = 2 * np.pi / omega_c     
t_span  = (0, 5 * T_c)
t_eval  = np.linspace(*t_span, 2000)


# --------------- equation of motion -----------------
def lorentz(t, y, q=q, m=m, B0=B0, B_dir=B_dir):
    """RHS of  dy/dt = [v, (q/m)(v × B)]."""
    v = y[3:6]
    B = B0 * B_dir
    a = (q / m) * np.cross(v, B)
    return np.concatenate([v, a])


# --------------- solve ------------------------------
sol = solve_ivp(lorentz, t_span, y0, method="RK45",
                t_eval=t_eval, rtol=1e-10, atol=1e-12)

x, y_coord, z = sol.y[0], sol.y[1], sol.y[2]


if __name__ == "__main__":
    print(f"Cyclotron frequency  ω_c = {omega_c:.4e} rad/s")
    print(f"Cyclotron period     T_c = {T_c:.4e} s")
    print(f"Larmor radius        r_L = {v_perp / omega_c:.4e} m")
    print(f"Integration points       = {sol.t.size}")
    print(f"Solver success           = {sol.success}")
