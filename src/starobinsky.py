import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# ============================================================
# Parameters
# ============================================================

alpha = 0.01
a = np.sqrt(2.0 / 3.0)


# ============================================================
# Dynamical system
# ============================================================

def system(N, y):

    h, psi, chi = y

    h_prime = -4.0 * np.pi * h * chi**2

    psi_prime = chi

    potential_derivative = (
        2.0
        * alpha**2
        * a
        * np.exp(-a * psi)
        * (1.0 - np.exp(-a * psi))
    )

    chi_prime = (
        4.0 * np.pi * chi**3
        - 3.0 * chi
        - potential_derivative / h**2
    )

    return [h_prime, psi_prime, chi_prime]


# ============================================================
# Consistent initial Hubble parameter
# ============================================================

def initial_h(psi0, chi0):

    numerator = (
        (8.0 * np.pi / 3.0)
        * alpha**2
        * (1.0 - np.exp(-a * psi0))**2
    )

    denominator = (
        1.0
        - (4.0 * np.pi / 3.0) * chi0**2
    )

    if denominator <= 0:
        raise ValueError(
            "Initial chi is outside the physical region."
        )

    return np.sqrt(numerator / denominator)


# ============================================================
# Event: end of inflation
#
# epsilon = 4*pi*chi^2
#
# Inflation ends when epsilon = 1.
# ============================================================

def end_of_inflation(N, y):

    chi = y[2]

    epsilon = 4.0 * np.pi * chi**2

    return epsilon - 1.0


end_of_inflation.terminal = True
end_of_inflation.direction = 1


# ============================================================
# Initial conditions
# ============================================================

psi0 = 2.12

chi0_values = [
    -0.30,
    -0.20,
    -0.10,
     0.00,
     0.10,
     0.20,
     0.30
]


# ============================================================
# Plot
# ============================================================

plt.figure(figsize=(12, 10))


for chi0 in chi0_values:

    h0 = initial_h(psi0, chi0)

    y0 = [
        h0,
        psi0,
        chi0
    ]

    solution = solve_ivp(
        system,
        (0.0, 60.0),
        y0,
        method="RK45",
        max_step=0.05,
        rtol=1e-8,
        atol=1e-10,
        events=end_of_inflation
    )

    if not solution.success:

        print(
            f"Failed for chi0={chi0}: "
            f"{solution.message}"
        )

        continue

    h = solution.y[0]
    psi = solution.y[1]
    chi = solution.y[2]

    # --------------------------------------------------------
    # Epsilon
    # --------------------------------------------------------

    epsilon = 4.0 * np.pi * chi**2

    # --------------------------------------------------------
    # End of inflation
    # --------------------------------------------------------

    if len(solution.t_events[0]) > 0:

        N_end = solution.t_events[0][0]

        psi_end = solution.y_events[0][0][1]
        chi_end = solution.y_events[0][0][2]

        print(
            f"chi0 = {chi0:+.2f}   "
            f"N_end = {N_end:.4f}   "
            f"psi_end = {psi_end:.4f}   "
            f"chi_end = {chi_end:.4f}"
        )

        # Mark end of inflation
        plt.scatter(
            psi_end,
            chi_end,
            s=50,
            zorder=5
        )

    else:

        print(
            f"chi0 = {chi0:+.2f}: "
            "inflation did not end within the integration range."
        )

    # --------------------------------------------------------
    # Phase-space trajectory
    # --------------------------------------------------------

    plt.plot(
        psi,
        chi,
        linewidth=1.2,
        alpha=0.7,
        label=fr"$\chi_0={chi0:+.2f}$"
    )


# ============================================================
# Slow-roll attractor
# ============================================================

psi_sr = np.linspace(0.15, 3.0, 1000)

chi_sr = (
    -a
    / (4.0 * np.pi)
    * np.exp(-a * psi_sr)
    / (1.0 - np.exp(-a * psi_sr))
)

plt.plot(
    psi_sr,
    chi_sr,
    "--",
    linewidth=2.5,
    label="Slow-roll attractor"
)


# ============================================================
# Plot formatting
# ============================================================

plt.xlabel(
    r"$\psi = \phi/M_{\rm Pl}$",
    fontsize=18
)

plt.ylabel(
    r"$\chi = \psi' = d\psi/dN$",
    fontsize=18
)

plt.title(
    r"Starobinsky Model: Phase-Space Attractor",
    fontsize=20
)

plt.xlim(0.3, 2.5)
plt.ylim(-0.35, 0.15)

plt.axhline(
    0,
    linewidth=0.8,
    alpha=0.2
)

plt.axvline(
    0,
    linewidth=0.8,
    alpha=0.2
)
plt.grid(alpha=1)

plt.legend(
    fontsize=11
)

plt.tight_layout()

plt.savefig('phase_space_attractor')
plt.show()
# ============================================================
# Slow-roll approximation test
# ============================================================

plt.figure(figsize=(12, 8))

for chi0 in chi0_values:

    h0 = initial_h(psi0, chi0)

    y0 = [
        h0,
        psi0,
        chi0
    ]

    solution = solve_ivp(
        system,
        (0.0, 60.0),
        y0,
        method="RK45",
        max_step=0.05,
        rtol=1e-8,
        atol=1e-10,
        events=end_of_inflation
    )

    if not solution.success:
        print(
            f"Failed for chi0={chi0}: "
            f"{solution.message}"
        )
        continue

    psi = solution.y[1]
    chi = solution.y[2]

    # --------------------------------------------------------
    # Analytical slow-roll attractor evaluated at
    # the numerical psi values
    # --------------------------------------------------------

    chi_sr = (
        -a
        / (4.0 * np.pi)
        * np.exp(-a * psi)
        / (1.0 - np.exp(-a * psi))
    )

    # --------------------------------------------------------
    # Difference between numerical and slow-roll solutions
    # --------------------------------------------------------

    delta_chi = chi - chi_sr

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    plt.plot(
        solution.t,
        delta_chi,
        linewidth=1.5,
        label=fr"$\chi_0={chi0:+.2f}$"
    )


# ============================================================
# Formatting
# ============================================================

plt.axhline(
    0,
    linewidth=1.0,
    alpha=0.5
)

plt.xlabel(
    r"$N$",
    fontsize=18
)

plt.ylabel(
    r"$\Delta\chi=\chi_{\rm numerical}-\chi_{\rm SR}$",
    fontsize=18
)

plt.title(
    r"Deviation from the Slow-Roll Attractor",
    fontsize=20
)

plt.grid(alpha=0.2)

plt.legend(fontsize=11)

plt.tight_layout()

plt.savefig('Deviation_slow_roll')
plt.show()

# ============================================================
# Relative error of the slow-roll approximation
# ============================================================

plt.figure(figsize=(12, 8))

for chi0 in chi0_values:

    h0 = initial_h(psi0, chi0)

    y0 = [
        h0,
        psi0,
        chi0
    ]

    solution = solve_ivp(
        system,
        (0.0, 60.0),
        y0,
        method="RK45",
        max_step=0.05,
        rtol=1e-8,
        atol=1e-10,
        events=end_of_inflation
    )

    if not solution.success:
        print(
            f"Failed for chi0={chi0}: "
            f"{solution.message}"
        )
        continue

    psi = solution.y[1]
    chi = solution.y[2]

    # --------------------------------------------------------
    # Slow-roll attractor evaluated at the numerical psi
    # --------------------------------------------------------

    chi_sr = (
        -a
        / (4.0 * np.pi)
        * np.exp(-a * psi)
        / (1.0 - np.exp(-a * psi))
    )

    # --------------------------------------------------------
    # Relative error
    # --------------------------------------------------------

    relative_error = np.abs(
        (chi - chi_sr) / chi_sr
    )

    # Convert to percentage
    relative_error_percent = 100.0 * relative_error

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    plt.plot(
        solution.t,
        relative_error_percent,
        linewidth=1.5,
        label=fr"$\chi_0={chi0:+.2f}$"
    )


# ============================================================
# Formatting
# ============================================================

plt.xlabel(
    r"$N$",
    fontsize=18
)

plt.ylabel(
    r"Relative error in $\chi_{\rm SR}$ (\%)",
    fontsize=18
)

plt.title(
    r"Accuracy of the Slow-Roll Attractor",
    fontsize=20
)

# Log scale is useful because the error spans
# several orders of magnitude.
plt.yscale("log")

plt.grid(
    alpha=0.2,
    which="both"
)

plt.legend(
    fontsize=11
)

plt.tight_layout()

plt.savefig(
    "slow_roll_relative_error.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig('slow_roll_relative_error')
plt.show()

# ============================================================
# Slow-roll parameter epsilon
# ============================================================

plt.figure(figsize=(12, 8))

for chi0 in chi0_values:

    h0 = initial_h(psi0, chi0)

    y0 = [
        h0,
        psi0,
        chi0
    ]

    solution = solve_ivp(
        system,
        (0.0, 60.0),
        y0,
        method="RK45",
        max_step=0.05,
        rtol=1e-8,
        atol=1e-10,
        events=end_of_inflation
    )

    if not solution.success:
        print(
            f"Failed for chi0={chi0}: "
            f"{solution.message}"
        )
        continue

    chi = solution.y[2]

    # --------------------------------------------------------
    # First slow-roll parameter
    # --------------------------------------------------------

    epsilon = 4.0 * np.pi * chi**2

    # --------------------------------------------------------
    # Plot epsilon
    # --------------------------------------------------------

    plt.plot(
        solution.t,
        epsilon,
        linewidth=1.5,
        label=fr"$\chi_0={chi0:+.2f}$"
    )


# ============================================================
# End of inflation: epsilon = 1
# ============================================================

plt.axhline(
    1.0,
    linestyle="--",
    linewidth=1.5,
    alpha=0.7,
    label=r"$\epsilon=1$ (end of inflation)"
)


# ============================================================
# Formatting
# ============================================================

plt.xlabel(
    r"$N$",
    fontsize=18
)

plt.ylabel(
    r"$\epsilon = 4\pi\chi^2$",
    fontsize=18
)

plt.title(
    r"First Slow-Roll Parameter",
    fontsize=20
)

# Log scale makes both the slow-roll regime and the
# approach to epsilon = 1 easier to see.
plt.yscale("log")

plt.grid(
    alpha=0.2,
    which="both"
)

plt.legend(
    fontsize=11
)

plt.tight_layout()

plt.savefig(
    "epsilon_vs_N.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig('epsilon_vs_N')
plt.show()
