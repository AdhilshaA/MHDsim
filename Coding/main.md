P464 - MHD Project
---
- Implemented by [Adhilsha A](https://github.com/AdhilshaA)
- Last updated on *14-03-2024*

This website follows the updates on my project on MHD. The project is about the `study of the Galactic Dynamos`. Each task in the project gets one step closer to the final goal. The project is divided into 3 tasks.
________

- [P464 - MHD Project](#p464---mhd-project)
- [TASK 1: Solving the diffusion equation](#task-1-solving-the-diffusion-equation)
  - [1.1 Introduction to Diffusion equation](#11-introduction-to-diffusion-equation)
  - [1.2 Equation forms](#12-equation-forms)
  - [1.3 Boundary conditions](#13-boundary-conditions)
  - [1.4 Implementation](#14-implementation)
    - [1.4.1 Spatial derivative](#141-spatial-derivative)
    - [1.4.2 RK4 time stepping](#142-rk4-time-stepping)
  - [1.5 Constants and parameters for simulation](#15-constants-and-parameters-for-simulation)
    - [1.5.1 Study of evolution of Magnetic field strength](#151-study-of-evolution-of-magnetic-field-strength)
    - [1.5.2 Study of evolution of pitch angle.](#152-study-of-evolution-of-pitch-angle)
    - [1.5.3 Study of Magnetic field evolution with different seed fields and different boundary conditions](#153-study-of-magnetic-field-evolution-with-different-seed-fields-and-different-boundary-conditions)
- [TASK 2: (TBD)](#task-2-tbd)
  - [2.1 Introduction and equations](#21-introduction-and-equations)
  - [2.3 Implementation](#23-implementation)
  - [2.4 Results](#24-results)
    - [2.4.1 Sample implementation](#241-sample-implementation)
    - [2.4.2 Decay rate at different spatial steps](#242-decay-rate-at-different-spatial-steps)
    - [2.4.3 Pitch angle evolution](#243-pitch-angle-evolution)
    - [2.4.4 Magnetic field evolution with different seed fields](#244-magnetic-field-evolution-with-different-seed-fields)
    - [2.4.5 Magnetic field evolution with different boundary conditions](#245-magnetic-field-evolution-with-different-boundary-conditions)
    - [2.4.6 Finding critical dynamo number](#246-finding-critical-dynamo-number)
- [Main Project](#main-project)
  - [Abstract](#abstract)
  - [1. Introduction](#1-introduction)
  - [2. Methods](#2-methods)
    - [2.1. Finite difference](#21-finite-difference)
    - [2.2. Time-stepping using RK4](#22-time-stepping-using-rk4)
    - [2.3. Implementation of BCs](#23-implementation-of-bcs)
    - [2.4. Calculation of Global growth rate](#24-calculation-of-global-growth-rate)
    - [2.5. Intuition behind Velocity profiles](#25-intuition-behind-velocity-profiles)
    - [2.6. Addition of Disk-flaring](#26-addition-of-disk-flaring)
  - [3. Simulation and results](#3-simulation-and-results)
    - [3.1. Effect of non-zero vertical outflow](#31-effect-of-non-zero-vertical-outflow)
    - [3.2. Effect of different profiles of vertical outflow](#32-effect-of-different-profiles-of-vertical-outflow)
    - [3.3. Disk flaring and vertical outflow](#33-disk-flaring-and-vertical-outflow)
  - [Conclusions](#conclusions)
- [References](#references)
________

## TASK 1: Solving the diffusion equation

### 1.1 Introduction to Diffusion equation

The first task is to solve the diffusion equation in 1D. The induction equation is given by:
$$ \frac{\partial \overline B_r}{\partial t} = - \frac{\overline V}{r} \frac{\partial}{\partial r}(r\overline B_r) -  \frac{\partial}{\partial z}(\overline V_z \overline B_r) -  \frac{\partial}{\partial z}(\alpha \overline B_{\phi}) + \eta_t \left[  \frac{\partial}{\partial r} \left( \frac{1}{r}  \frac{\partial}{\partial r}(r \overline B_r)\right) +  \frac{\partial ^2 \overline B_r}{\partial z^2} \right] $$

$$
\frac{\partial \overline B_{\phi}}{\partial t} = - q\Omega \overline B_{r} -  \frac{\partial}{\partial r}(\overline V_r \overline B_{\phi}) -  \frac{\partial}{\partial z}(\overline V_z \overline B_{\phi}) +  \frac{\partial}{\partial z}(\alpha \overline B_{r}) + \eta_t \left[  \frac{\partial}{\partial r} \left( \frac{1}{r}  \frac{\partial}{\partial r}(r \overline B_{\phi})\right) +  \frac{\partial ^2 \overline B_{\phi}}{\partial z^2} \right]
$$

where $\overline B_r$, $\overline B_{\phi}$ and $\overline B_z$ are the components of the magnetic field, $\overline V_r$, $\overline V_{\phi}$ and $\overline V_z$ are the components of the velocity field.

For this task, we solve the diffusion equation omitting the $\nabla \times (\overline V \times \overline B)$ and $\alpha$ terms. The equation we have to solve becomes:
$$
\frac{\partial \overline B_r}{\partial t} = \eta_t \left[  \frac{\partial}{\partial r} \left( \frac{1}{r}  \frac{\partial}{\partial r}(r \overline B_r)\right) +  \frac{\partial ^2 \overline B_r}{\partial z^2} \right] $$

$$
\frac{\partial \overline B_{\phi}}{\partial t} = \eta_t \left[  \frac{\partial}{\partial r} \left( \frac{1}{r}  \frac{\partial}{\partial r}(r \overline B_{\phi})\right) +  \frac{\partial ^2 \overline B_{\phi}}{\partial z^2} \right]
$$

Since both $\overline B_r$ and $\overline B_{\phi}$ have same form of equation, we generalize the form of the equation to:
$$
\frac{\partial \overline B}{\partial t} = \eta_t \left[  \frac{\partial}{\partial r} \left( \frac{1}{r}  \frac{\partial}{\partial r}(r \overline B)\right) +  \frac{\partial ^2 \overline B}{\partial z^2} \right]
$$

### 1.2 Equation forms

In case of $\overline B_r(z)$ and $\overline B_{\phi}(z)$, the equations become:
$$
\frac{\partial \overline B}{\partial t} = \eta_t \frac{\partial ^2 \overline B}{\partial z^2}
$$

In case of $\overline B_r(r)$ and $\overline B_{\phi}(r)$ under no-$z$ appproximation, the equations become:
$$
\frac{\partial \overline B}{\partial t} = \eta_t \left[  \frac{\partial}{\partial r} \left( \frac{1}{r}  \frac{\partial}{\partial r}(r \overline B)\right) -  \frac{\pi ^2}{4 h^2}\overline B \right] $$

On expanding the equation, we get:
$$
\frac{\partial \overline B}{\partial t} = \eta_t \left[ \frac{\partial ^2 \overline B}{\partial r^2} +  \frac{1}{r}  \frac{\partial \overline B}{\partial r} - \frac{\overline B}{r^2}  -  \frac{\pi ^2}{4 h^2}\overline B \right]
$$

These are the equations to be solved in this task.

### 1.3 Boundary conditions

We will be using the simplest boundary conditions for this task. The boundary conditions are:
$$
\overline B_r = \overline B_{\phi} = 0 \quad \text{at} \quad (z = -h, h)\text{ or }(r = 0,r_{max}) $$

We will be taking dimensionless units with $h = 1$, $\eta_t = 1$.

In $r$ implementation, we will be starting with $r$ = 0.01 since the value $r = 0$ is can possibly cause division by zero errors.

### 1.4 Implementation

For our implementation we will be using finite difference method for defining the spatial derivatives and RK4 for the time stepping / evolution. The code can be observed at [task1_code](task1_code.html).

#### 1.4.1 Spatial derivative

The spatial derivatives are calculated using the finite difference method. here, we are using 6th order finite difference method for the spatial derivatives. The 6th order finite difference method is given by:
$$
\frac{d f}{d x} = \frac{1}{60\delta x}\left( - f_{i-3} + 9f_{i-2} - 45f_{i-1} + 45f_{i+1} - 9f_{i+2} + f_{i+3} \right) $$
$$
\frac{d^2 f}{d x^2} = \frac{1}{180\delta x^2}\left( 2f_{i-3} - 27f_{i-2} + 270f_{i-1} - 490f_{i} + 270f_{i+1} - 27f_{i+2} + 2f_{i+3} \right) $$

Other orders can be used and is available for use in the code.

For the sake of obtaining smooth derivatives as well for the sake of maintaining the boundary conditions, we use `ghost zones`. By default, 'relative anti-symmetric' ghost zones are used where the ghost cell at $i$ distance away from boundary $b$ is given by:
$$ f_{b-i} = 2f_{b} - f_{b+i} \quad \text{at the start of the spatial domain} $$
$$ f_{b+i} = 2f_{b} - f_{b-i} \quad \text{at the end of the spatial domain} $$

where $i = 1, 2, 3, ..., n$ ($n$ is the half the order of the finite difference method).

Other ghost zones types, such as 'symmetric' and 'anti-symmetric' can also be used and is available for use in the code.

Using the 6th order finite difference method and 'relative anti-symmetric' ghost zones, a following magnetic field and the spatial derivatives are obtained (test case):

<div style="text-align:center;">
    <img src="figures\initial_condition_B.png" alt="Sample Field" />
</div>

<div style="text-align:center;">
    <div style="display:inline-block;">
        <img src="figures\first_spatial_derivative.png" alt="First Spatial Derivative" />
    </div>
    <div style="display:inline-block;">
        <img src="figures\second_spatial_derivative.png" alt="Second Spatial Derivative" />
    </div>
</div>


As you can see the spatial derivatives are smooth and the second derivative is zero at the boundaries as expected from the 'relative anti-symmetric' ghost zones. This works exceptionally well for the $z$ implementation. For the $r$ implementation, the ghost zones cannot help much in keeping the boundary conditions constant as the equation contains first and second derivatives as well as the function itself. Hence, the boundary conditions are not maintained well in the $r$ implementation. This is a known issue and is being worked on.

#### 1.4.2 RK4 time stepping

The time stepping is done using the RK4 method. Our evolution function can be given of form:
$$
\frac{\partial \overline B}{\partial t} = f(\overline B, t) $$
Then, the RK4 method evolution for a single time step can be given by:
$$
\kappa_1 = \delta t f(\overline B, t) $$
$$
\kappa_2 = \delta t f(\overline B + \frac{\kappa_1}{2}, t + \frac{\delta t}{2}) $$
$$
\kappa_3 = \delta t f(\overline B + \frac{\kappa_2}{2}, t + \frac{\delta t}{2}) $$
$$
\kappa_4 = \delta t f(\overline B + k_3, t + \delta t) $$

$$
\overline B(t + \delta t) = \overline B(t) + \frac{1}{6}(\kappa_1 + 2\kappa_2 + 2\kappa_3 + \kappa_4) $$


The thing to note here about the generalized form, with regards to the RK4 method, is that:
- The spatial derivatives do not change values with the change ($\overline B \rightarrow \overline B + \kappa/2$) in the RK4 method. This is due to the symmetry of the spatial derivative equations. Hence, the spatial derivatives are calculated only once at the start of the time step. Thus, it is not mentioned in the equation form.

Our method, as mentioned above, is to find the spatial derivatives using finite difference method at a certain time step and then use the RK4 method to evolve the system to the next time step. This is done repetitively to get the solution at all time steps.

### 1.5 Constants and parameters for simulation

For $z$ implementation, we will be using the following parameters:
- $z\ \epsilon\ (-1, 1)$ `# range of z`
- $t\ \epsilon\ [0, 1]$ `# range of time`
- $N_z = 80$ `# number of spatial steps`
- $N_t = 5000$ `# number of time steps`
- $order = 6$ `# order of the finite difference method`
- *ghostzone_type* = 'relative anti-symmetric' `# type of ghostzone`

For $r$ implementation, we will be using the following parameters:
- $r\ \epsilon\ (0.01, 8.01)$ `# range of r`
- $t\ \epsilon\ [0, 1]$ `# range of time`
- $N_r = 100$ `# number of spatial steps`
- $N_t = 4000$ `# number of time steps`
- $order = 6$ `# order of the finite difference method`
- *ghostzone_type* = 'relative anti-symmetric' `# type of ghostzone`

The seed field used for the z-implementation is:

$$ B_{\phi}(z) = 
3 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} \pi \right) + 1 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} 2\pi \right) + 2 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} 5\pi \right)
$$

whereas the seed field used for the r-implementation is:

$$ B_{\phi}(r) =
3 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} \pi \right) + 1 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} 2\pi \right) + 2 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} 5\pi \right)
$$

The results are:

<!-- ![image](figures\evolution_B_phi_z.png)
![image](figures\evolution_B_phi_r.png) -->
<div style="text-align:center;">
    <img src="figures\evolution_B_phi_z.png" alt="Evolution of B_phi_z" />
</div>

<div style="text-align:center;">
    <img src="figures\evolution_B_phi_r.png" alt="Evolution of B_phi_r" />
</div>


#### 1.5.1 Study of evolution of Magnetic field strength

For this part, I am setting $\overline B_r = 0$ and evolving $\overline B_{\phi}$ over time. The magnetic field strength, $|B_{\phi}|$, is given at different points in the spatial domain as a function of time. The results are (N.B. The decay rate is given in the legend off the plots):

<div style="text-align:center;">
    <img src="figures\B_phi_z_strength__timesteps.png" alt="Evolution of B_phi_z Strength Over Timesteps" />
</div>

<div style="text-align:center;">
    <img src="figures\B_phi_r_strength__timesteps.png" alt="Evolution of B_phi_r Strength Over Timesteps" />
</div>


from these results, taking the middle point in the spatial domain, the magnetic field strength, $|B_{\phi}|$, is given as a function of time in *log*-scale. Then it can be fitted to find the `decay rate` as well. The results are:

<!-- ![B_phi_z_decay](figures\B_phi_z_decay.png)
![B_phi_r_decay](figures\B_phi_r_decay.png) -->
<div style="text-align:center;">
    <img src="figures\B_phi_z_decay.png" alt="B_phi_z Decay" />
</div>

<div style="text-align:center;">
    <img src="figures\B_phi_r_decay.png" alt="B_phi_r Decay" />
</div>

#### 1.5.2 Study of evolution of pitch angle.

For this part, I am setting the seed $\overline B_r$ and $\overline B_{\phi}$ to non-zero values and evolving them over time. The pitch angle is given by:
$$p = \tan^{-1}\left(\frac{\overline B_{r}}{\overline B_{\phi}}\right)$$

The seed fields taken for this simulation are, for the z-implementation:
$$
\text{{B}}_{\phi}(z) = 3 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} \pi \right) + 1 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} 2\pi \right) + 2 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} 5\pi \right)
$$

$$
\text{{B}}_{r}(z) = 3 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} \pi \right) + 2.5 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} 2\pi \right) + 1 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} 5\pi \right)
$$
and for the r-implementation:
$$
\text{{B}}_{\phi}(r) = 3 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} \pi \right) + 1 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} 2\pi \right) + 2 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} 5\pi \right)
$$

$$
\text{{B}}_{r}(r) = 3 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} \pi \right) + 2.5 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} 2\pi \right) + 1 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} 5\pi \right)
$$

The time evolution of the fields are given below for reference:

<!-- ![Bphi_Br_z_r_evolution](figures\Bphi_Br_z_r_evolution.png) -->

<div style="text-align:center;">
    <img src="figures\Bphi_Br_z_r_evolution.png" alt="Bphi_Br_z_r_evolution" />
</div>


The results of the pitch angle evolution at certain specific spatial steps were obtained as well. The results:

<!-- ![pitch_angle_z_timeevolution_at_diff_spatialsteps](figures\pitch_angle_z_timeevolution_at_diff_spatialsteps.png)
![pitch_angle_r_timeevolution_at_diff_spatialsteps](figures\pitch_angle_r_timeevolution_at_diff_spatialsteps.png) -->

<div style="text-align:center;">
    <img src="figures\pitch_angle_z_timeevolution_at_diff_spatialsteps.png" alt="pitch_angle_z_timeevolution_at_diff_spatialsteps" />
</div>

<div style="text-align:center;">
    <img src="figures\pitch_angle_r_timeevolution_at_diff_spatialsteps.png" alt="pitch_angle_r_timeevolution_at_diff_spatialsteps" />
</div>

To get another perspective, the pitch angle variation over spatial domain at different time steps is also found. The results are:

<!-- ![pitch_angle_z_spatialevolution_at_diff_timesteps.png](figures\pitch_angle_z_spatialevolution_at_diff_timesteps.png)
![pitch_angle_r_spatialvariation_at_diff_timesteps](figures\pitch_angle_r_spatialevolution_at_diff_timesteps.png) -->

<div style="text-align:center;">
    <img src="figures\pitch_angle_z_spatialevolution_at_diff_timesteps.png" alt="pitch_angle_z_spatialevolution_at_diff_timesteps" />
</div>

<div style="text-align:center;">
    <img src="figures\pitch_angle_r_spatialevolution_at_diff_timesteps.png" alt="pitch_angle_r_spatialevolution_at_diff_timesteps" />
</div>

#### 1.5.3 Study of Magnetic field evolution with different seed fields and different boundary conditions

Here, I did two trials, first one with zero boundary conditions. Then, the second one with non-zero boundary conditions and different seed fields from the previous trial. 

The trial 1 seed field with zero boundary conditions is, for the z-implementation:
$$
B_{\phi}( z) = 3 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} \pi \right) + 1 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} 2\pi \right) + 2 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} 5\pi \right)
$$

$$
B_{\phi}( z) = 3 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} \pi \right) + 2.5 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} 2\pi \right) + 1 \sin \left( \frac{{z - z_i}}{{z_f - z_i}} 5\pi \right)
$$
and for the r-implementation:
$$
B_{\phi}( r) = 3 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} \pi \right) + 1 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} 2\pi \right) + 2 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} 5\pi \right)
$$

$$
B_{r}( z) = 3 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} \pi \right) + 2.5 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} 2\pi \right) + 1 \sin \left( \frac{{r - r_i}}{{r_f - r_i}} 5\pi \right)
$$
The results are:

<!-- ![Bphi_Br_z_r_evolution_trial1](figures\Bphi_Br_z_r_evolution_trial1.png) -->

<div style="text-align:center;">
    <img src="figures\Bphi_Br_z_r_evolution_trial1.png" alt="Bphi_Br_z_r_evolution_trial1" />
</div>

Here, with zero boundary conditions, the boundaries are maintained well in the z-implementation. The reason for this is 'relative anti-symmetric' ghost zones which keeps the $\frac{\partial^2 \overline B}{\partial z^@}$ zero at the boundaries.

However, the boundaries are not maintained well in the r-implementation, even though the it is hardly noticeable in the plots. This is because the magnetic field is decaying and the boundary condition is already zero. Hence, the boundary condition is 'seen' to be maintained well. But, the difference is noticeable in the next trial.

With the different seed fields used in the implementations, all magnetic field are seen to decay with time, as expected.

The trial 2 seed field with non-zero boundary conditions is, for the z-implementation:
$$
B_{\phi}( z) = 3 \cos \left( \frac{{z - z_i}}{{z_f - z_i}} \pi \right) + 1 \cos \left( \frac{{z - z_i}}{{z_f - z_i}} 2\pi \right) + 2 \cos \left( \frac{{z - z_i}}{{z_f - z_i}} 5\pi \right)
$$

$$
B_{r}(z) = 3 \cos \left( \frac{{z - z_i}}{{z_f - z_i}} \pi \right) + 2.5 \cos \left( \frac{{z - z_i}}{{z_f - z_i}} 2\pi \right) + 1 \cos \left( \frac{{z - z_i}}{{z_f - z_i}} 5\pi \right)
$$
and for the r-implementation:
$$
B_{\phi}( r) = 3 \cos \left( \frac{{r - r_i}}{{r_f - r_i}} \pi \right) + 1 \cos \left( \frac{{r - r_i}}{{r_f - r_i}} 2\pi \right) + 2 \cos \left( \frac{{r - r_i}}{{r_f - r_i}} 5\pi \right)
$$

$$
B_{r}( r) = 3 \cos \left( \frac{{r - r_i}}{{r_f - r_i}} \pi \right) + 2.5 \cos \left( \frac{{r - r_i}}{{r_f - r_i}} 2\pi \right) + 1 \cos \left( \frac{{r - r_i}}{{r_f - r_i}} 5\pi \right)
$$
The results are:

<!-- ![Bphi_Br_z_r_evolution_trial2](figures\Bphi_Br_z_r_evolution_trial2.png) -->

<div style="text-align:center;">
    <img src="figures\Bphi_Br_z_r_evolution_trial2.png" alt="Bphi_Br_z_r_evolution_trial2" />
</div>

Here in $z$ implementation, even the non-zero boundaries are maintained well as expected. 

As you can see in the $r$ implementation, the non-zero boundary conditions are not maintained at all. This is because the ghost zones cannot help much in keeping the boundary conditions constant as the equation contains first and second derivatives as well as the function itself. The 'relative anti-symmetric' ghost zones can only keep the second derivative zero at the boundaries. 

The decay with of magnetic field with time is evident here, even with different seed fields used in the implementations.

________

## TASK 2: (TBD)

### 2.1 Introduction and equations

From the original induction equation, we omit the terms involving $\overline V$ from both equations and the $\alpha$ term from the equation for $B_{\phi}$, in order to get a $\alpha\omega$ dynamo. Task 2 will be focused on solving these equations. THe equations are:

$$ \frac{\partial \overline B_r}{\partial t} = -  \frac{\partial}{\partial z}(\alpha \overline B_{\phi}) + \eta_t \left[  \frac{\partial}{\partial r} \left( \frac{1}{r}  \frac{\partial}{\partial r}(r \overline B_r)\right) +  \frac{\partial ^2 \overline B_r}{\partial z^2} \right] $$

$$ \frac{\partial \overline B_{\phi}}{\partial t} = - q\Omega \overline B_{r} + \eta_t \left[  \frac{\partial}{\partial r} \left( \frac{1}{r}  \frac{\partial}{\partial r}(r \overline B_{\phi})\right) +  \frac{\partial ^2 \overline B_{\phi}}{\partial z^2} \right] $$

By applying no-$z$ approximation and expanding, equations become:
$$ \frac{\partial \overline B_r}{\partial t} = -  \frac{2}{\pi h}(\alpha \overline B_{\phi}) + \eta_t \left[ \frac{\partial ^2 \overline B}{\partial r^2} +  \frac{1}{r}  \frac{\partial \overline B}{\partial r} - \frac{\overline B}{r^2} +  \frac{\pi^2 }{4h^2}\overline B_r \right] $$

$$ \frac{\partial \overline B_{\phi}}{\partial t} = - q\Omega \overline B_{r} + \eta_t \left[  \frac{\partial ^2 \overline B}{\partial r^2} +  \frac{1}{r}  \frac{\partial \overline B}{\partial r} - \frac{\overline B}{r^2} +  \frac{\pi^2 }{4h^2}\overline B_{\phi} \right] $$

### 2.3 Implementation

As done before, we will be using the dimensionless form of the equations. Hence it becomes:
$$ \frac{\partial \tilde B}{\partial \tilde t} = -  \frac{2}{\pi}(\alpha \tilde B_{\phi}) + \left[ \frac{\partial ^2 \tilde B}{\partial\tilde r^2} +  \frac{1}{\tilde r}  \frac{\partial \tilde B}{\partial\tilde r} - \frac{\tilde B}{\tilde r^2} +  \frac{\pi^2 }{4}\tilde B \right] $$

$$ \frac{\partial \tilde B}{\partial \tilde t} = - q\Omega \tilde B_{r} + \left[  \frac{\partial ^2 \tilde B}{\partial\tilde r^2} +  \frac{1}{\tilde r}  \frac{\partial \tilde B}{\partial\tilde r} - \frac{\tilde B}{\tilde r^2} +  \frac{\pi^2}{4}\tilde B \right] $$

We will use the following definitions in our implementation:

$$\alpha = \alpha_0(r)\hat\alpha(z);\ \ \ \hat \alpha = sin\left(\pi \tilde z\right);\ \ \ q=-\frac{r}{\Omega}\frac{\partial \Omega}{\partial r}; \ \ \ \Omega = \frac{\Omega_0}{\sqrt{1 + (r/r_\omega)^2}}$$
$$R_\alpha = \frac{\alpha_0h}{\eta_t};\ \ \ R_\omega = \frac{-q\Omega h^2}{\eta_t}$$

We will be keeping $\alpha$ a constant to properly study the dynamo, instead of zero at midplane.

The final equations are:
$$ \frac{\partial \tilde B}{\partial \tilde t} = -  \frac{2}{\pi}(R_\alpha\hat\alpha \tilde B_{\phi}) + \left[ \frac{\partial ^2 \tilde B}{\partial\tilde r^2} +  \frac{1}{\tilde r}  \frac{\partial \tilde B}{\partial\tilde r} - \frac{\tilde B}{\tilde r^2} +  \frac{\pi^2 }{4}\tilde B \right] $$

$$ \frac{\partial \tilde B}{\partial \tilde t} = R_\omega B_{r} + \left[  \frac{\partial ^2 \tilde B}{\partial\tilde r^2} +  \frac{1}{\tilde r}  \frac{\partial \tilde B}{\partial\tilde r} - \frac{\tilde B}{\tilde r^2} +  \frac{\pi^2}{4}\tilde B \right] $$

For solving these equations, we will use the same spatial derivative function using finite difference as given before in Section [1.4.1](#141-spatial-derivative). For time-stepping, we implement coupled RK4 which solved the coupled differential equations including every single spatial position.

Find the code at [task2_code](task2_code.html).

### 2.4 Results

#### 2.4.1 Sample implementation

To test the solver as well as obtain preliminary results, we ran the following simulation. THe results shows amplification as expected.

![image](figures\task2\Br_Bphi_time_evolution.png)
![image](figures\task2\B_strength_time_evolution.png)

#### 2.4.2 Decay rate at different spatial steps

To find the decay rate, we take the, we first plotted the magnetic field strength evolution over time in a log-plot, which shows straight lines as expected. Some local variation (as seen as one deviant curve in the figure below) could be due to the choice of the seed field. as it is not seen in every trial done.
![image](figures\task2\B_strength_time_evolution_at_different_r.png)

To find the decay rate, we take the slope of the log-plot of the magnetic field strength at the middle of the spatial domain, assuming exponential decay. The results are:
![image](figures\task2\B_strength_time_evolution_at_middle_r.png)

#### 2.4.3 Pitch angle evolution
The evolution of pitch angle was also tracked and is given below. The choice of seed field play a major role here since these simulations were of short time periods.
![image](figures\task2\pitch_angle_time_evolution.png)

#### 2.4.4 Magnetic field evolution with different seed fields

For studying the evolution with different seed fields, we ran two different simulations as shown below.
![image](figures\task2\Br_Bphi_time_evolution_seed_field1.png)
![image](figures\task2\Br_Bphi_time_evolution_seed_field2.png)
The magnetic strength evolution of the different seed fields are:
![image](figures\task2\B_strength_time_evolution_seed_field1.png)
![image](figures\task2\B_strength_time_evolution_seed_field2.png)


#### 2.4.5 Magnetic field evolution with different boundary conditions
To study the effect of boundary conditions, we ran the simulation with symmetric, anti-symmetric and relative anti-symmetric ghost zones, each of which ensures a different boundary condition. 
The symmetric boundary condition keeps the first derivative zero, the anti-symmetric boundary condition keeps the boundary at zero itself, whereas the relative anti-symmetric keeps the second derivative zero. We found that the anti-symmetric boundary conditions suits the current task and will be followed in the next section.

![image](figures\task2\Br_Bphi_time_evolution_ghost_zone_anti-symmetric.png)
![image](figures\task2\Br_Bphi_time_evolution_ghost_zone_symmetric.png)
![image](figures\task2\Br_Bphi_time_evolution_ghost_zone_relative_anti-symmetric.png)
The evolutions of magnetic field strength for each of the implementation above is:
![image](figures\task2\B_strength_time_evolution_ghost_zone_anti-symmetric.png)
![image](figures\task2\B_strength_time_evolution_ghost_zone_symmetric.png)
![image](figures\task2\B_strength_time_evolution_ghost_zone_relative_anti-symmetric.png)

#### 2.4.6 Finding critical dynamo number
To find the critical dynamo number at some time, we plot the local growth rate at that time and check at which spatial position does it cross zero (move from decay to growth). THe values used in this implementation is inspired from the numerical simulation in [(Chamandy et.al. 2014)](#references). It was used such that those $R_omega$ and $D_C$ values are present in our value range across $r$.


The variation of $R_\omega$ and Dynamo number $D$ across $r$ was obtained as:
![image](figures\task2\omega_Dynamo_number_vs_r_trial12.png)
The magnetic field strength evolution of this implementation turned out to be:
![image](figures\task2\B_strength_time_evolution_trial12.png)
The local growth rate and corresponding Critical Dynamo Number at unit time is marked in the figure below.
![image](figures\task2\gamma_vs_r_trial12_withDc.png)

The above method was done, drawing inspiration from SS21 Figure 11.1 [(Shukurov & Subramanian 2021)](#references) for our experimental setting. Reading into another paper [(Chamandy et.al. 2014)](#references), I believe it might have been better to run the simulation with a range of values for global dynamo number, each simulation extending to the stable regime. Then, the critical dynamo number can be found by checking the value of fynamo number at which the steady decay changes into steady growth. This required a large amount of time and thus could not be explored right now.

For doing a similar comparison, our simulation will be done for a longer time period till a global decay rate is obtained and will be compared with the theoretical growth rate predicted by  [(Chamandy et.al. 2014)](#references). The values there were given by:
$$D_C = - \frac{\pi^5}{32};\ \ \ \gamma = \frac{2}{\pi} t_d^{-1} (\sqrt{-D} - \sqrt{-D_C}) $$


With the theoretical values at $t_d = 1$ (unit diffusion time, since time is scaled) and the global decay rate obtained at a longer simulation, the below figure was plotted.

![image](figures\task2\gamma_vs_r_trial12_global.png)

## Main Project
________

### Abstract


Studying `galactic mean-field dynamos` involves exploring the intricate dynamics governing cosmic magnetic fields, which provide invaluable insights into the formation and evolution of galaxies. Our research focuses on understanding the `impact of vertical outflows` within these systems, as they greatly influence the distribution of gas and magnetic fields within galaxies and are crucial to the star formation rates and matter-energy transport processes. In this study, our objective was to investigate how vertical outflows, characterized by various velocity profiles, influence the evolution of galactic magnetic fields. 

Utilizing a combination of finite differencing and RK4 methods, we simulated the kinematic regime to observe the behavior of magnetic fields over time. Our findings demonstrate that the incorporation of velocity outflows, particularly stronger ones, `hinders the growth of magnetic fields`, and in some cases, even leads to decay. This highlights the significant role that vertical outflows play in `shaping the evolution of magnetic fields` in galactic dynamos, particularly emphasizing their interaction with radial terms.

These insights can not only deepen our understanding of the dynamics of galactic magnetic fields but also contribute to broader discussions surrounding galaxy formation. By studying the intricate interplay between various physical processes within cosmic structures, this project sheds light on aspects of galactic evolution that are closer to our comprehension of the universe's complex mechanisms.

### 1. Introduction

From the galactic mean field theory, we have magnetic field evolution equations that are given by the induction equation. The induction equation [(Chamandy 2024)](#references) is given by:

$$ \frac{\partial \overline B_r}{\partial t} = - \frac{\overline V_r}{r} \frac{\partial}{\partial r}(r\overline B_r) -  \frac{\partial}{\partial z}(\overline V_z \overline B_r) -  \frac{\partial}{\partial z}(\alpha \overline B_{\phi}) + \eta_t \left[  \frac{\partial}{\partial r} \left( \frac{1}{r}  \frac{\partial}{\partial r}(r \overline B_r)\right) +  \frac{\partial ^2 \overline B_r}{\partial z^2} \right] $$

$$
\frac{\partial \overline B_{\phi}}{\partial t} = - q\Omega \overline B_{r} -  \frac{\partial}{\partial r}(\overline V_r \overline B_{\phi}) -  \frac{\partial}{\partial z}(\overline V_z \overline B_{\phi}) +  \frac{\partial}{\partial z}(\alpha \overline B_{r}) + \eta_t \left[  \frac{\partial}{\partial r} \left( \frac{1}{r}  \frac{\partial}{\partial r}(r \overline B_{\phi})\right) +  \frac{\partial ^2 \overline B_{\phi}}{\partial z^2} \right]
$$

where $\overline B_r$, $\overline B_{\phi}$ and $\overline B_z$ are the components of the magnetic field, $\overline V_r$, $\overline V_{\phi}$ and $\overline V_z$ are the components of the velocity field, $\alpha$ is the dynamo alpha effect, $\Omega$ is the angular velocity of the galaxy, $q$ is the shear parameter, and $\eta_t$ is the turbulent magnetic diffusivity.

In this project, we will be focusing on the impact of `vertical outflows` on the evolution of galactic magnetic fields compared to task 2 where we we were studying the impact of $\alpha$-$\omega$ dynamo effects. The velocity field is given by $\overline V_z$ and terms involving $\overline V_r$ as well as $\alpha$ term from the equation for $B_{\phi}$ are omitted from the equations. Task 3 will be focused on solving these equations. THe equations are:

$$ \frac{\partial \overline B_r}{\partial t} = -  \frac{\partial}{\partial z}(\overline V_z \overline B_{r}) -  \frac{\partial}{\partial z}(\alpha \overline B_{\phi}) + \eta_t \left[  \frac{\partial}{\partial r} \left( \frac{1}{r}  \frac{\partial}{\partial r}(r \overline B_r)\right) +  \frac{\partial ^2 \overline B_r}{\partial z^2} \right] $$

$$ \frac{\partial \overline B_{\phi}}{\partial t} = -  \frac{\partial}{\partial z}(\overline V_z \overline B_{\phi}) - q\Omega \overline B_{r} + \eta_t \left[  \frac{\partial}{\partial r} \left( \frac{1}{r}  \frac{\partial}{\partial r}(r \overline B_{\phi})\right) +  \frac{\partial ^2 \overline B_{\phi}}{\partial z^2} \right] $$

By applying no-$z$ approximation and expanding, equations become:
$$ \frac{\partial \overline B_r}{\partial t} = -  \frac{2}{\pi h}(\overline V_z \overline B_{r}) -  \frac{2}{\pi h}(\alpha \overline B_{\phi}) + \eta_t \left[ \frac{\partial ^2 \overline B}{\partial r^2} +  \frac{1}{r}  \frac{\partial \overline B}{\partial r} - \frac{\overline B}{r^2} +  \frac{\pi^2 }{4h^2}\overline B_r \right] $$

$$ \frac{\partial \overline B_{\phi}}{\partial t} = -  \frac{2}{\pi h}(\overline V_z \overline B_{\phi}) - q\Omega \overline B_{r} + \eta_t \left[  \frac{\partial ^2 \overline B}{\partial r^2} +  \frac{1}{r}  \frac{\partial \overline B}{\partial r} - \frac{\overline B}{r^2} +  \frac{\pi^2 }{4h^2}\overline B_{\phi} \right] $$

We can convert these equations into dimensionless form and the same will be used in the simulations. Hence it becomes:
$$ \frac{\partial \tilde B_r}{\partial\tilde t} = -  \frac{2}{\pi\tilde h}(\tilde V_z \tilde B_{r}) -  \frac{2}{\pi\tilde h}(\tilde\alpha \tilde B_{\phi}) + \left[ \frac{\partial ^2 \tilde B}{\partial\tilde r^2} +  \frac{1}{\tilde r}  \frac{\partial \tilde B}{\partial\tilde r} - \frac{\tilde B}{\tilde r^2} +  \frac{\pi^2 }{4\tilde h^2}\tilde B_r \right] $$

$$ \frac{\partial \tilde B_{\phi}}{\partial\tilde t} = -  \frac{2}{\pi\tilde h}(\tilde V_z \tilde B_{\phi}) - q\tilde\Omega \tilde B_{r} + \left[  \frac{\partial ^2 \tilde B}{\partial\tilde r^2} +  \frac{1}{\tilde r}  \frac{\partial \tilde B}{\partial\tilde r} - \frac{\tilde B}{\tilde r^2} +  \frac{\pi^2 }{4\tilde h^2}\tilde B_{\phi} \right] $$

To make it dimensionless, we use the following definitions in our implementation [(Chamandy et.al. 2014)](#references):

$$\eta_t = \frac{1}{3}lu;\quad t_d = h^2/\eta_t;$$
$$\Omega = \frac{\Omega_0}{\sqrt{1 + (r/r_\omega)^2}};\quad q=-\frac{r}{\Omega}\frac{\partial \Omega}{\partial r};\quad \tilde \Omega = \frac{\Omega h^2}{\eta_t}$$
$$\alpha = \alpha_0(r)\hat\alpha(z);\quad \hat \alpha = sin\left(\pi \tilde z\right); \quad \alpha_0(r) = \frac{\Omega l^2}{h};\quad \tilde \alpha_0 = \frac{\alpha_0 h}{\eta_t}$$
$$R_\alpha = \tilde\alpha_0;\quad R_\omega = -q\tilde \Omega; \quad D = R_\alpha R_\omega$$
$$\tilde h = 1;\quad \tilde t_d = 1; \quad\tilde \eta_t = 1;\quad \tilde r = \frac{r}{h};\quad \tilde h = 1;\quad \tilde t_d = 1; \quad\tilde B = \frac{B}{B_0}; \quad \tilde V_z = \frac{V_z h}{\eta_t}$$

Here, $l$ is the correlation length of the turbulence, $u$ is the turbulent velocity, $h$ is the scale height of the disk, $r_\omega$ is the radial scale of the angular velocity profile, and $B_0$ is the initial magnetic field strength. The values of the parameter used ([(Chamandy et.al. 2014)](#references)) will be given in Section [3](#3-simulation-and-results).
We will be keeping $\hat\alpha = 1$ a constant instead of zero at midplane, for our study.
For disk flaring, we will use $\tilde h = h(r) / h$ and the same will affect the $\alpha_0$ term as well. This won't affect the conversion to dimensionless form, as we are using the same $h$ in the conversion.

The final equations are:
$$ \frac{\partial \tilde B_r}{\partial\tilde t} = -  \frac{2}{\pi\tilde h}(\tilde V_z \tilde B_{r}) -  \frac{2}{\pi\tilde h}(R_\alpha \tilde B_{\phi}) + \left[ \frac{\partial ^2 \tilde B}{\partial\tilde r^2} +  \frac{1}{\tilde r}  \frac{\partial \tilde B}{\partial\tilde r} - \frac{\tilde B}{\tilde r^2} +  \frac{\pi^2 }{4\tilde h^2}\tilde B_r \right] $$

$$ \frac{\partial \tilde B_{\phi}}{\partial\tilde t} = -  \frac{2}{\pi\tilde h}(\tilde V_z \tilde B_{\phi}) + R_\omega \tilde B_{r} + \left[  \frac{\partial ^2 \tilde B}{\partial\tilde r^2} +  \frac{1}{\tilde r}  \frac{\partial \tilde B}{\partial\tilde r} - \frac{\tilde B}{\tilde r^2} +  \frac{\pi^2 }{4\tilde h^2}\tilde B_{\phi} \right] $$

For solving these equations, we will use the same spatial derivative function using finite difference as given before in Section [2.1](#21-finite-difference). For time-stepping, we implement coupled RK4 which solved the coupled differential equations including every single spatial position.

Similar equations have been solved in the past by others [(Shukurov et.al. 2016)](#references) where though significantly different in terms of implementation and settings, the resutls can be compared in the kinematic regime. A result from the same is shown in Figure [1](#fig1) [(Shukurov et.al. 2016)](#references).

<center><img src="figures\task3\theory_pic.png" alt="image" id = 'fig1'></center>
<center><b>Figure 1</b>: The inset shows different levels of vertical outflow: absence (solid), weak (dashed), moderate (dotted) and strong (dot-dashed)</center><br>

Here, in the inset of Figure [1](#fig1), the magnetic field strength evolution in the kinematic regime can be seen for different vertical outflows. Stronger outflows are seen to hinder the growth of magnetic fields, and in some cases, even lead to decay. This is the main focus of our study.

Find the code at [task3_code](task3_code.html).

### 2. Methods
#### 2.1. Finite difference
The spatial derivatives are calculated using the finite difference method. here, we are using 6th order finite difference method for the spatial derivatives. The 6th order finite difference method is given by [(Bradenburg 2019)](#references) as follows:
$$
\frac{d f}{d x} = \frac{1}{60\delta x}\left( - f_{i-3} + 9f_{i-2} - 45f_{i-1} + 45f_{i+1} - 9f_{i+2} + f_{i+3} \right) $$
$$
\frac{d^2 f}{d x^2} = \frac{1}{180\delta x^2}\left( 2f_{i-3} - 27f_{i-2} + 270f_{i-1} - 490f_{i} + 270f_{i+1} - 27f_{i+2} + 2f_{i+3} \right) $$

Other orders can be used and is available for use in the code.

By using finite difference, we can calculate the spatial derivatives at every spatial point. This approximation is used to convert the partial differential equations into ordinary differential equations, which can be solved using time-stepping methods (RK4 in our case). This approximation method is called `method of lines`.

#### 2.2. Time-stepping using RK4
The time stepping is done using the RK4 method. Our evolution function can be given of form:
$$ \frac{\partial \overline B}{\partial t} = f(\overline B, t) $$
But, since our $\overline B$ has two components, we have two equations to solve.
$$ \frac{\partial \overline B_r}{\partial t} = f_r(\overline B_r, \overline B_{\phi}, t) $$
$$ \frac{\partial \overline B_{\phi}}{\partial t} = f_\phi(\overline B_r, \overline B_{\phi}, t) $$

Then, the RK4 method evolution for a single time step $\delta t$ from the current $\overline B(t)$  can be given by the coupled RK4 equations:
$$\kappa_{1r} = \delta t f_r(\overline B_r(t), \overline B_{\phi}(t), t) $$
$$\kappa_{1\phi} = \delta t f_\phi(\overline B_r(t), \overline B_{\phi}(t), t) $$
$$ \kappa_{2r} = \delta t f_r(\overline B_r(t) + \frac{\kappa_{1r}}{2}, \overline B_{\phi}(t) + \frac{\kappa_{1\phi}}{2}, t + \frac{\delta t}{2}) $$
$$ \kappa_{2\phi} = \delta t f_\phi(\overline B_r(t) + \frac{\kappa_{1r}}{2}, \overline B_{\phi}(t) + \frac{\kappa_{1\phi}}{2}, t + \frac{\delta t}{2}) $$
$$ \kappa_{3r} = \delta t f_r(\overline B_r(t) + \frac{\kappa_{2r}}{2}, \overline B_{\phi}(t) + \frac{\kappa_{2\phi}}{2}, t + \frac{\delta t}{2}) $$
$$ \kappa_{3\phi} = \delta t f_\phi(\overline B_r(t) + \frac{\kappa_{2r}}{2}, \overline B_{\phi}(t) + \frac{\kappa_{2\phi}}{2}, t + \frac{\delta t}{2}) $$
$$ \kappa_{4r} = \delta t f_r(\overline B_r(t) + \kappa_{3r}, \overline B_{\phi}(t) + \kappa_{3\phi}, t + \delta t) $$
$$ \kappa_{4\phi} = \delta t f_\phi(\overline B_r(t) + \kappa_{3r}, \overline B_{\phi}(t) + \kappa_{3\phi}, t + \delta t) $$
$$ \overline B_r(t + \delta t) = \overline B_r(t) + \frac{1}{6}(\kappa_{1r} + 2\kappa_{2r} + 2\kappa_{3r} + \kappa_{4r}) $$
$$ \overline B_{\phi}(t + \delta t) = \overline B_{\phi}(t) + \frac{1}{6}(\kappa_{1\phi} + 2\kappa_{2\phi} + 2\kappa_{3\phi} + \kappa_{4\phi}) $$


The thing to note here about the generalized form, with regards to the RK4 method, is that:
- The spatial derivatives do not change values with the change ($\overline B \rightarrow \overline B + \kappa/2$) in the RK4 method. This is due to the symmetry of the spatial derivative equations. Hence, the spatial derivatives are calculated only once at the start of the time step. Thus, it is not mentioned in the equation form.
- The order of calculation of $\kappa$ values is important. The order is $\kappa_{1}$, $\kappa_{2}$, $\kappa_{3}$, $\kappa_{4}$. Both components of $\overline B$ are calculated at each step before moving to the next step. This is important to maintain the integrity of the coupled equations.

#### 2.3. Implementation of BCs
For the sake of obtaining smooth derivatives as well for the sake of maintaining the boundary conditions, we use `ghost zones`. 
To implement ghost zones, we do the following:
- At the start and end of the spatial domain, we add ghost cells outside the domain of interest. We are free to choose the number of ghost cells to add and to set the values of the ghost cells.
- The number of ghost cells to add is half the order of the finite difference method. Hence, for 6th order finite difference method, we add 3 ghost cells at the start and end of the spatial domain. THis will be useful for calculating the derivatives at the boundaries smoothly.
- The values of the ghost cells are set according to the type of ghost zones used and the boundary conditions.

By default, '`anti-symmetric`' ghost zones are used where the ghost cell at $i$ distance away from boundary $b$ is given by:
$$ f_{b-i} = -f_{b+i} \quad \text{at the start of the spatial domain} $$
$$ f_{b+i} = -f_{b-i} \quad \text{at the end of the spatial domain} $$
This is used to implement the boundary conditions $B_r = 0$ and $B_{\phi} = 0$ at the boundaries.

`Symmetric` ghost zones are used where the ghost cell at $i$ distance away from boundary $b$ is given by:
$$ f_{b-i} = f_{b+i} \quad \text{at the start of the spatial domain} $$
$$ f_{b+i} = f_{b-i} \quad \text{at the end of the spatial domain} $$
This is used to implement the boundary conditions $\frac{\partial B_r}{\partial r} = 0$ and $\frac{\partial B_{\phi}}{\partial r} = 0$ at the boundaries.

`Relative anti-symmetric` ghost zones are used where the ghost cell at $i$ distance away from boundary $b$ is given by:
$$ f_{b-i} = 2f_{b} - f_{b+i} \quad \text{at the start of the spatial domain} $$
$$ f_{b+i} = 2f_{b} - f_{b-i} \quad \text{at the end of the spatial domain} $$
This is used to implement the boundary conditions $\frac{\partial^2 B_r}{\partial r^2} = 0$ and $\frac{\partial^2 B_{\phi}}{\partial r^2} = 0$ at the boundaries.

where $i = 1, 2, 3, ..., n$ ($n$ is the half the order of the finite difference method). These ghost zones are used to maintain the boundary conditions as well as to obtain smooth derivatives at the boundaries and are free to be changed in the code. By default, 'anti-symmetric' ghost zones are used to implement the boundary conditions $B_r = 0$ and $B_{\phi} = 0$ at the boundaries.

#### 2.4. Calculation of Global growth rate

In the calculation of the global growth rate and the iterative process for running simulations until the desired growth rate is achieved, several systematic steps are followed:

1. **Initial Evolution Period**: During the first 10 diffusion times, the magnetic field evolves without monitoring growth rates.

2. **Collection of Growth Rates**: At the end of each subsequent diffusion time, the growth rate is determined at every spatial point. This involves analyzing the evolution of magnetic field strength ($|B|$) over the last few hundred time steps at each spatial point and calculating the slope of the logarithm of $|B|$. This process is repeated for all spatial points, generating a distribution of growth rates ($\gamma$) across the radius ($r$).

3. **Mean Growth Rate Calculation**: The mean growth rate ($\overline{\gamma}$) across all radii is calculated from the obtained distribution of growth rates.

4. **Tolerance Check**: The absolute difference between the mean growth rate ($\overline{\gamma}$) and the individual growth rates ($\gamma$) is compared against a predefined tolerance threshold.

5. **Simulation Continuation**: The simulation continues until the absolute difference falls below the specified tolerance. Furthermore, this condition must be maintained for an additional 3 diffusion times to ensure stability.

By meticulously following these steps, we can ensure accurate determination of global growth rates in galactic dynamo simulations.

#### 2.5. Intuition behind Velocity profiles

Typical choice of the Velocity profile for vertical outflow as given in [(Chamandy et.al. 2014)](#references) is:
$$V_z = \frac{V_zz}{h}$$
Hence, at a certain $z$, the Vertical outflow will be a constant.

In addition to exploring different values of constant profiles, I was also thinking to explore radial variations of the Vertical outflow due to reasons discussed before in Section [1](#1-introduction).

Since the star formations are more in smaller radius, a radially decreasing profile seemed correct. In addition to this, if we consider the effect of higher gravity in smaller radius, we can also try a velocity profile which decreases till some radius $r_v$ and then increase due to the lower effect of gravity (henceforth mentioned as U-profile). To explore both using similar equation, I defined it in the following way:

$$ V_z = V_0 \times \left(1 + 0.5 \times \left(\left(\frac{r - r_{v}}{r_{v}}\right)^2 - 1\right)\right) $$

The last two profiles mentioned above are obtained by putting $r_v$ = $r_f$ (final radius) and $r_v$ = 10 $kpc$ respectively. The profiles are shown below in Figures [2](#fig2) and [3](#fig3).
<center><img src="figures\task3\trial31_V_z_vs_r.png" alt="image" id = 'fig2'></center>
<center><b>Figure 2</b>: Velocity profile that decreases till the end of spatial domain.</center><br>
<center><img src="figures\task3\trial32_V_z_vs_r.png" alt="image" id = 'fig3'></center>
<center><b>Figure 3</b>: Velocity U-profile that decreases till 10 kpc, then increases</center><br>
Note that this is just two sample profiles mimicking the shapes intuitive to the reasoning above.

#### 2.6. Addition of Disk-flaring

To include disk flaring, we changed applied the radial profile $h(r)$ on the dimensionless form and propagated it through the other affected parameters ($R_\alpha, D, etc/$) and evolution. The form of disk flaring used is:

$$ h(r) = h \times \sqrt{1 + \left(\frac{r}{r_h}\right)^2} $$

Here, $r_h$ controls the disk flaring and is set to $10 kpc$ as per [(Chamandy et.al. 2014)](#references)
The form of $h(r)$ is given below in Figure [4](#fig4).
<center><img src="figures\task3\trial25_h_vs_r.png" alt="image" id = 'fig4'></center>
<center><b>Figure 4</b>: Visualization of disk flaring</center><br>

Note that the scaling and making the parameters are kept unaffected by this change as we use $h$ and the corresponding $t_d$ for that purpose. See the implementation in [task_3-code](task3_code.html).

### 3. Simulation and results

All simulation paramters, trial settings, and the results are available in [Github](https://github.com/AdhilshaA/MHDsim_project) and [task_3-code](task3_code.html).

We used the following parameters and the forms discussed in Section [1](introduction) for the implementation. 

Practical parameters:
r_i = 0.01 kpc
r_f = 15.01 kpc
h = 0.4 kpc
omega_0 = 220.0000 km / (kpc s)
r_omega = 2.0 kpc
r_h = 10.0 kpc
B_0 = 1e-04 G
t_d = 0.4693 Gyr
eta_t = 0.3409 kpc2 / Gyr

Dimensionless parameters:
r_i = 0.0250
r_f = 37.5250
h = 1.0000
omega_0 = 105.6000
r_omega = 5.0000
r_h = 25.0000
B_0 = 1.0000
t_d = 1.0000
eta_t = 1.0000

The profiles of $R_\omega$, $R_\alpha$ and $D$ across $r$ are given below in Figures [5](#fig5) to [7](#fig7) for reference.
<center><img src="figures\task3\trial32_R_omega_vs_r.png" alt="image" id = 'fig5'></center>
<center><b>Figure 5</b>: Variation of R_omega with respect to r</center><br>
<center><img src="figures\task3\trial32_omega_vs_r.png" alt="image" id = 'fig6'></center>
<center><b>Figure 6</b>: Variation of R_alpha with respect to r</center><br>
<center><img src="figures\task3\trial32_Dynamo_number_vs_r.png" alt="image" id = 'fig7'></center>
<center><b>Figure 7</b>: Variation of Dynamo number with respect to r</center><br>

Note that the disk flaring is not used for studies in Sections [3.1](#31-effect-of-non-zero-vertical-outflow) and [3.2](#32-effect-of-different-profiles-of-vertical-outflow) for the sake of simplicity and to focus on the effect of vertical outflow. The disk flaring is used in the study in Section [3.3](#33-disk-flaring-and-vertical-outflow).
#### 3.1. Effect of non-zero vertical outflow

To study the effect of vertical outflow compared to the previous results from Task 2, we ran the simulations with $V_z$ = 0, 2, 4, 6 , 8, 10 $km/s$ for the modified equations and implementation in Section [1](#1-introduction). We will first focus on results from $V_z$ = 0 and $V_z$ = 4 $km/s$ (representing non-zero $V_z$) for comparison of different aspects of magnetic field evolution.

The results from $V_z$ = 0 and $V_z$ = 4 $km/s$ are shown below. Below, in Figures [8](#fig8) to [11](#fig11), are the final eigen modes and the magnetic field evolution of both simulations.
<center><img src="figures\task3\trial19_Final_eigen_modes.png" alt="image" id = 'fig8'></center>
<center><b>Figure 8</b>: Final eigen mode of magnetic field in the  absence of vertical outflow</center><br>
<center><img src="figures\task3\trial21_Final_eigen_modes.png" alt="image" id = 'fig9'></center>
<center><b>Figure 9</b>: Final eigen mode of magnetic field in the  presence of vertical outflow (V_z = 4 km/s)</center><br>
<center><img src="figures\task3\trial19_B_strength_evolution_slider_at_diff_r.gif" alt="image" id = 'fig10'></center>
<center><b>Figure 10</b>: Evolution of magnetic field strength in the  absence of vertical outflow</center><br>
<center><img src="figures\task3\trial21_B_strength_evolution_slider_at_diff_r.gif" alt="image" id = 'fig11'></center>
<center><b>Figure 11</b>: Evolution of magnetic field strength in the  presence of vertical outflow (V_z = 4 km/s)</center><br>

As you can see the final eigen modes are the same and is not affected by the vertical outflow here.

We can also look at the evolution of growth rate across $r$ for both simulations given below in Figures [12](#fig12) and [13](#fig13).
<center><img src="figures\task3\trial19_gamma_evolution.gif" alt="image" id = 'fig12'></center>
<center><b>Figure 12</b>: Evolution of growth rate in the  absence of vertical outflow</center><br>
<center><img src="figures\task3\trial21_gamma_evolution.gif" alt="image" id = 'fig13'></center>
<center><b>Figure 13</b>: Evolution of growth rate in the  presence of vertical outflow (V_z = 4 km/s)</center><br>

The global growth rate $\Gamma$ achieved by the simulation with Vertical outflow is lower than that without the vertical outflow. This is expected as the vertical outflow removes the magnetic field from the system that was growing. The pitch angle evolutions are not severely affected by the vertical outflow as seen below in Figures [14](#fig14) and [15](#fig15).
<center><img src="figures\task3\trial19_pitch_angle_vs_t_at_diff_r.png" alt="image" id = 'fig14'></center>
<center><b>Figure 14</b>: Evolution of pitch angle in the  absence of vertical outflow</center><br>
<center><img src="figures\task3\trial21_pitch_angle_vs_t_at_diff_r.png" alt="image" id = 'fig15'></center>
<center><b>Figure 15</b>: Evolution of pitch angle in the  presence of vertical outflow (V_z = 4 km/s)</center><br>

Looking at the general trend across all the above range of values for $V_z$, the important one is the change in the global growth rate $\Gamma$ obtained. The variation of $\Gamma$ across various values of $V_z$ is shown below in Figure [16](#fig16).
<center><img src="figures\task3\gamma_evolution_all_Vz2_new.gif" alt="image" id = 'fig16'></center>
<center><b>Figure 16</b>: Variation of global growth rate with respect to V_z</center><br>

As you can see, the global growth rate is negatively affected by the outflow, and the same is expected the outflow removes the magnetic field from the system that was growing. As $V_z$ further increases, the $\Gamma$ crosses zero and becomes negative, i.e. decaying magnetic field. This is also expected as the vertical outflow can remove magnetic faster than it evolves, causing decay.

Since, all other parameters are the same, we can see how the vertical outflow affected the growth of the magnetic field. The results can be compared to the result mentioned in [(Shukurov et.al. 2006)](#references) and discussed in Section [1](#1-introduction).

#### 3.2. Effect of different profiles of vertical outflow

Apart from constant velocity profile across $r$, we will also explore the decreasing and U-shaped profile effect on the magnetic field evolution.

The magnetic field evolution, pitch angle evolution and final eigen functions are given below for both simulations in Figures [17](#fig17) to [22](#fig22).
<center><img src="figures\task3\trial31_Final_eigen_modes.png" alt="image" id = 'fig17'></center>
<center><b>Figure 17</b>: Final eigen mode of magnetic field in the  decreasing vertical outflow profile</center><br>
<center><img src="figures\task3\trial32_Final_eigen_modes.png" alt="image" id = 'fig18'></center>
<center><b>Figure 18</b>: Final eigen mode of magnetic field in the  U-shaped vertical outflow profile</center><br>
<center><img src="figures\task3\trial31_B_strength_evolution_slider_at_diff_r.gif" alt="image" id = 'fig19'></center>
<center><b>Figure 19</b>: Evolution of magnetic field strength in the  decreasing vertical outflow profile</center><br>
<center><img src="figures\task3\trial32_B_strength_evolution_slider_at_diff_r.gif" alt="image" id = 'fig20'></center>
<center><b>Figure 20</b>: Evolution of magnetic field strength in the  U-shaped vertical outflow profile</center><br>
<center><img src="figures\task3\trial31_pitch_angle_vs_t_at_diff_r.png" alt="image" id = 'fig21'></center>
<center><b>Figure 21</b>: Evolution of pitch angle in the  decreasing vertical outflow profile</center><br>
<center><img src="figures\task3\trial32_pitch_angle_vs_t_at_diff_r.png" alt="image" id = 'fig22'></center>
<center><b>Figure 22</b>: Evolution of pitch angle in the  U-shaped vertical outflow profile</center><br>

Looking at their evolution of $\gamma$ given below, we can see that the decreasing profile has increased the global growth rate which can be explained looking at how the growth rates evolve and interact radially. See the animation below in Figure [23](#fig23).
<center><img src="figures\task3\gamma_evolution_diff_Vz_profiles.gif" alt="image" id = 'fig23'></center>
<center><b>Figure 23</b>: Evolution of growth rate in different Velocity profiles; constant (blue), decreasing (green) and U-shaped (red)</center><br>

As the growth at the radius 2 $kpc$ (maximum Dynamo number as shown in the beginning of Section [3](#3-simulation-and-results)) is influencing and spreading to other radial positions, the decreasing profile has the same effect as a decreased $V_z$ value which implies higher growth rate. This also implies how a different profile (actual profile as well) would interact different with magnetic field evolution.

#### 3.3. Disk flaring and vertical outflow

With the implementation of disk flaring as mentioned in Section [2.6](#26-addition-of-disk-flaring), we can see the interaction of the disk flaring with the vertical outflow. For this, we ran the simulation with disk-flaring and $V_z$ ranging from 0 to 10 $km/s$ in steps of 2 $km/s$. Looking at their evolution of growth rates given below, we can see how disk flaring affected the impact of vertical outflows in the system. See the animation below in Figure [24](#fig24).
<center><img src="figures\task3\gamma_evolution_all_Vz_disk_flaring_new.gif" alt="image" id = 'fig24'></center>
<center><b>Figure 24</b>: Evolution of growth rate in different Velocity profiles with disk flaring</center><br>

For higher vertical outflows, the growth rates are spread slowly in radius as the eariler interactions are affected by the disk flaring. We can compare these to the results from simulations without disk flaring to see the effect of disk flaring on the magnetic field evolution. In comparison to results shown in Figure [16](#fig16), the disk flaring has affected the growth rates and the spread of growth rates across the radius. See the animation below in Figure [25](#fig25).
<center><img src="figures\task3\gamma_evolution_all_Vz_disk_flaring_comparison.gif" alt="image" id = 'fig25'></center>
<center><b>Figure 25</b>: Comparison of evolution of growth rate in different Velocity profiles with and without disk flaring</center><br>

As you can see in Figure [25](#fig25), the radial interactions are affected by disk flaring, this shows that with additional terms dependent on $r$, the vertical outflows impact in various rates or in even more complicated manner which is further not explored here.

### Conclusions

In our comprehensive exploration of galactic mean-field dynamos, we have meticulously observed and drawn significant conclusions from the simulations and results obtained. Across the span of our investigations, several key observations have emerged, shedding light on the intricate dynamics of cosmic magnetic fields within galactic systems.

- Our analysis indicates that the presence of vertical outflows exerts a notable influence on the growth of magnetic fields, distinctly contrasting with scenarios where such outflows are absent. Notably, we have observed a linear hindrance to magnetic field growth in the presence of vertical outflows, with stronger outflows removing magnetic fields at a rate surpassing their growth, eventually leading to decay. This fundamental insight, clearly seen in Figure [16](#fig16), underscores the crucial role of vertical outflows in modulating magnetic field evolution within galaxies.

- Different profiles of vertical outflows, as discussed in Section [3.2](#32-effect-of-different-profiles-of-vertical-outflow), have revealed intriguing correlations between the evolution of growth rates and Dynamo number profiles (See Figure [23](#fig23)). Notably, we have observed that decreasing profiles exhibit effects akin to weaker velocity outflows, underscoring the intricate interplay between velocity profiles and magnetic field evolution. This highlights the nature of interactions within galactic dynamo systems, where variations in velocity profiles can show diverse responses in magnetic field evolution.

- Our inclusion of radially changing systems, as detailed in Section [3.3](#33-disk-flaring-and-vertical-outflow) and seen in Figure [25](#fig25), has shown the impact of radial variations on the influence of outflows, affecting convergence rates and growth dynamics. This underscores the complexity inherent in galactic dynamo systems, where radial variations can significantly modulate the effects of vertical outflows on magnetic field evolution.

The significance of our findings extends beyond these systems, as they hold profound implications for understanding real-world galactic systems. Traditionally, velocity outflow terms have been neglected in favor of simplifying system dynamics. However, our research underscores the critical relevance of considering these terms, particularly in light of their impact on shaping the magnetic field evolution and their relation to star formation rates and matter transport within and from the system.

Moving forward, our study suggests promising aspects to be explored for future research. Specifically, the formulation and exploration of various profiles of velocity outflows, accounting for effects such as star formation rates and disk flaring at different radii, hold considerable potential. Additionally, the possibility of coupling these considerations with the relaxation of various other assumptions and exploring system dynamics beyond the kinematic regime presents exciting opportunities for advancing our understanding of galactic dynamo theory.

In conclusion, our comprehensive investigation into the interactions of vertical outflows within galactic mean-field dynamos has yielded profound insights into the complex dynamics governing magnetic field evolution. Through meticulous analysis and simulation, we have illuminated the intricate interplay between velocity outflows and magnetic field dynamics, paving the way for further advancements in galactic dynamo theory and our understanding of cosmic structures.

## References

1. Brandenburg, A. 2019. Computational aspects of astrophysical MHD and turbulence. In *Advances in nonlinear dynamos* (CRC Press). [[link](https://arxiv.org/abs/astro-ph/0109497)]
2. Chamandy, L., Shukurov, A., Subramanian, K.  & Stoker, K. 2014. MNRAS 443. [[link](https://doi.org/10.1093/mnras/stu1274)]
3. Chamandy, L. 2024. P464 Class Notes. (NISER)
4. Shukurov, A., Sokoloff, D., Subramanian, K., & Brandenburg, A. 2006. Galactic dynamo and helicity losses through fountain flow. *Astronomy & Astrophysics*, 448(2), L33-L36. [[link](https://doi.org/10.1051/0004-6361:200600011)]
5. Shukurov, A., & Subramanian, K. 2021. Astrophysical Magnetic Fields: From Galaxies to the Early Universe. (Cambridge: Cambridge University Press). [[link](https://doi.org/10.1017/9781139046657)]