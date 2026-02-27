# Reactor Kinetics Simulation
RP Activity 5: MATLAB-based transient simulator that solves the coupled non-linear ordinary differential equations for nuclear reactor point kinetics, lumped-parameter thermal-hydraulics, and a simplified 3-group decay heat model. The script models a dynamic reactor transient characterized by a continuous positive reactivity insertion (ramp) followed by an emergency reactor shutdown (SCRAM), actively accounting for Doppler and moderator temperature feedback mechanisms.

## Methodology

The simulation employs a state-space formulation, tracking 12 state variables over time: normalized neutron flux (1), delayed neutron precursors (6), fuel temperature (1), coolant temperature (1), and decay heat groups (3). The system of ordinary differential equations (ODEs) is stiff due to the varying time scales of prompt neutrons and delayed precursors, thus it is solved using MATLAB's `ode15s` solver.

Below are the governing physical equations used in the continuous ODE system and initial state formulation.

### 1. Neutronics (Point Kinetics with 6 Delayed Groups)

The neutron population is modeled using the standard Point Kinetics Equations (PKE). The flux $\phi(t)$ is normalized such that $\phi(0) = 1$, representing the relative power fraction.

**Neutron Flux Derivative:**


$$\frac{d\phi}{dt} = \left( \frac{\rho_{net}(t) - \beta}{\Lambda} \right) \phi(t) + \sum_{i=1}^6 \lambda_i C_i(t)$$

**Delayed Precursor Derivative (for groups $i = 1 \dots 6$):**


$$\frac{dC_i}{dt} = \frac{\beta_i}{\Lambda} \phi(t) - \lambda_i C_i(t)$$

**Initial Conditions (Steady-State):**


$$C_{i,0} = \frac{\beta_i}{\Lambda \cdot \lambda_i} \phi_0$$

*Where:*

* $\Lambda$ (`td`) = Prompt neutron lifetime / generation time
* $\beta_i$, $\lambda_i$ = Delayed neutron fractions and decay constants
* $\beta$ = Total delayed neutron fraction ($\sum \beta_i$)
* $C_i$ = Precursor concentrations

### 2. Reactivity and Feedback

Reactivity drives the transient and is influenced by both external control rod movements and internal temperature feedback (Doppler and moderator density effects).

**Net Reactivity:**


$$\rho_{net}(t) = \rho_{ext}(t) + \alpha_F (T_F(t) - T_{F,0}) + \alpha_C (T_C(t) - T_{C,0})$$

**External Reactivity Insertion ($\rho_{ext}$):**
The external reactivity simulates rod withdrawal (ramp) followed by a SCRAM (step drop).

$$
\rho_{ext}(t) =
\begin{cases}
5 \times 10^{-5} \cdot t & t < 20\,\text{s} \\
-5 \times 10^{-2} & t \ge 20\,\text{s}
\end{cases}
$$

*Where:*
* $\alpha_F$ = Fuel Doppler coefficient
* $\alpha_C$ = Coolant moderator coefficient
* $T_F$, $T_C$ = Fuel and Coolant temperatures

### 3. Decay Heat (3-Group Model)
The total reactor thermal power is the sum of instantaneous prompt fission power and delayed decay heat from fission products.

**Prompt Power Fraction:**

$$f_{prompt} = 1 - \sum_{j=1}^3 \gamma_j$$

**Decay Heat Derivative (for groups $j = 1 \dots 3$):**

$$\frac{dP_{d,j}}{dt} = \lambda_{h,j} (\gamma_j P_0 \phi(t)) - \lambda_{h,j} P_{d,j}(t)$$

**Initial Decay Heat (Steady-State):**

$$P_{d,j,0} = \gamma_j P_0 \phi_0$$

**Total Thermal Power:**

$$P_{total}(t) = P_0 \phi(t) f_{prompt} + \sum_{j=1}^3 P_{d,j}(t)$$

*Where:*
* $P_0$ = Nominal thermal power
* $\gamma_j$ = Fraction of total power for decay group $j$
* $\lambda_{h,j}$ = Decay constant for decay heat group $j$
* $P_{d,j}$ = Decay heat power contribution from group $j$

### 4. Thermal-Hydraulics (Lumped Parameter Model)
A zero-dimensional (0-D), two-node lumped parameter model calculates the average temperatures of the fuel and the coolant. 

**Fuel Temperature Derivative:**

$$\frac{dT_F}{dt} = \frac{(1 - f_c) P_{total}(t) - UA (T_F(t) - T_C(t))}{M_F C_{p,F}}$$

**Coolant Temperature Derivative:**
The coolant energy balance assumes a linear temperature rise across the core, meaning the exit temperature $T_{out} = 2T_C - T_{in}$.

$$\frac{dT_C}{dt} = \frac{f_c P_{total}(t) + UA (T_F(t) - T_C(t)) - 2 \dot{m} C_{p,C} (T_C(t) - T_{in})}{M_C C_{p,C}}$$

**Initial Conditions (Steady-State):**
Setting the time derivatives to zero yields the initial equilibrium temperatures:

$$T_{C,0} = T_{in} + \frac{P_0}{2 \dot{m} C_{p,C}}$$

$$T_{F,0} = T_{C,0} + \frac{(1 - f_c) P_0}{UA}$$

*Where:*
* $f_c$ = Fraction of heat deposited directly in the coolant (prompt gamma/neutron moderation)
* $UA$ = Overall heat transfer coefficient from fuel to coolant
* $M_F C_{p,F}$ = Heat capacity of the fuel
* $M_C C_{p,C}$ = Heat capacity of the coolant
* $\dot{m} C_{p,C}$ = Coolant mass flow heat capacity rate
* $T_{in}$ = Coolant inlet temperature
