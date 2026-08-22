import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# --------------------------------
# Dynamical system
# --------------------------------
def oscillator(t, z, omega, delta):
    x, y = z

    dxdt = y
    dydt = -omega**2 * x - delta * y

    return [dxdt, dydt]


# --------------------------------
# Parameters
# --------------------------------
omega = 1.0

cases = [
    ("Critical damping", 2 * omega),
    ("Overdamping", 3 * omega),
    ("Underdamping", 1 * omega)
]


# Initial conditions
initial_conditions = [
    [1.0, 0.0],
    [0.0, 1.0],
    [-1.0, 0.5],
    [1.5, -1.0],
    [-1.5, -0.5]
]


# --------------------------------
# Phase-space grid
# --------------------------------
x_values = np.linspace(-2, 2, 20)
y_values = np.linspace(-2, 2, 20)

X, Y = np.meshgrid(x_values, y_values)


# --------------------------------
# Plot
# --------------------------------
fig, axes = plt.subplots(1, 3, figsize=(18, 5))


for ax, (title, delta) in zip(axes, cases):

    # Vector field
    U = Y
    V = -omega**2 * X - delta * Y

    # Normalize vectors so that arrows have similar lengths
    magnitude = np.sqrt(U**2 + V**2)

    U_norm = U / (magnitude + 1e-10)
    V_norm = V / (magnitude + 1e-10)

    ax.quiver(
        X, Y,
        U_norm, V_norm,
        alpha=0.5
    )

    # Trajectories
    for z0 in initial_conditions:

        sol = solve_ivp(
            lambda t, z: oscillator(t, z, omega, delta),
            (0, 20),
            z0,
            t_eval=np.linspace(0, 20, 1500),
            rtol=1e-9,
            atol=1e-9
        )

        x = sol.y[0]
        y = sol.y[1]

        ax.plot(x, y, lw=1.5)

        # Initial condition
        ax.plot(x[0], y[0], 'ko', markersize=4)

    # Fixed point
    ax.plot(0, 0, 'ko', markersize=6)

    ax.set_title(title)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y=\dot{x}$")

    ax.axhline(0, color='black', lw=0.5)
    ax.axvline(0, color='black', lw=0.5)

    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)

    ax.grid(alpha=0.2)


plt.tight_layout()

plt.savefig(
    "damped_oscillator.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
