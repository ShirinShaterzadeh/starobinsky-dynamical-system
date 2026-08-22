# Dynamical Systems Analysis of Inflation

Numerical study of dynamical systems with an application to the inflationary dynamics of the Starobinsky model.

This project introduces phase-space analysis using the damped harmonic oscillator as a simple physical example, and then applies the same framework to investigate the inflationary attractor of the Starobinsky model.

## 1. Damped Harmonic Oscillator

As a simple physical example of a dynamical system, we consider the damped harmonic oscillator

$$
\ddot{x}+\delta\dot{x}+\omega^2x=0,
$$

where $$\omega$$ is the natural frequency and $$\delta$$ is the damping coefficient.

Introducing

$$
y=\dot{x},
$$

the second-order equation becomes the first-order dynamical system

$$
\begin{cases}
\dot{x}=y,\\
\dot{y}=-\omega^2x-\delta y.
\end{cases}
$$

The system has a fixed point at

$$
(x,y)=(0,0).
$$

The behavior of trajectories in phase space depends on the damping coefficient.

### Damping regimes

The eigenvalues of the system are determined by

$$
\lambda^2+\delta\lambda+\omega^2=0.
$$

Three qualitatively different regimes are considered.

#### Overdamping

For

$$
\delta>2\omega,
$$

the eigenvalues are distinct, real, and negative. The trajectories approach the fixed point without oscillating.

The fixed point is a **stable node**.

For the numerical example,

$$
\omega=1,\qquad \delta=3.
$$

#### Critical damping

For

$$
\delta=2\omega,
$$

the eigenvalues are repeated and negative. The system approaches the fixed point without oscillating.

The fixed point is a **degenerate stable node**.

For the numerical example,

$$
\omega=1,\qquad \delta=2.
$$

#### Underdamping

For

$$
\delta<2\omega,
$$

the eigenvalues are complex with negative real parts. The trajectories spiral toward the fixed point while their amplitude decreases.

The fixed point is a **stable spiral (spiral sink)**.

For the numerical example,

$$
\omega=1,\qquad \delta=1.
$$

### Numerical phase-space analysis

The three damping regimes are solved numerically using `scipy.integrate.solve_ivp`.

For each regime, several different initial conditions are used:

$$
(x_0,y_0)
=
(1,0),\,
(0,1),\,
(-1,0.5),\,
(1.5,-1),\,
(-1.5,-0.5).
$$

The trajectories are plotted in the phase space $$(x,y)$$.

Despite starting from different initial conditions, all trajectories approach the same fixed point. The geometry of the approach depends on the eigenvalue structure:

$$
\boxed{
\begin{array}{c}
\text{Overdamped} \rightarrow \text{stable node}\\
\text{Critical} \rightarrow \text{degenerate stable node}\\
\text{Underdamped} \rightarrow \text{stable spiral}
\end{array}
}
$$

This provides a simple example of how the eigenvalues of a dynamical system determine the local phase-space behavior and the type of attractor.
