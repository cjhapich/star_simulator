
# This code has been adapted from the Fortran 95 program STATSTAR laid out by Carroll & Ostile (2007)
# Many of the equations were borrowed and translated from the Fortran code to Python 3.10
# Due to the use of a match function, a Python 3.10 interpreter is required to run this program.

# I translated this code for fun, I claim no ownership of the methods implemented in this code!

# Get packages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

plt.style.use('dark_background')

# All units MKS
G = 6.67 * 10 ** -11
kB = 1.38 * 10 ** -23
c = 2.998 * 10 ** 8
mH = 1.67 * 10 ** -27
sigma = 5.67 * 10 ** -8
a = 4 * sigma / c

# Opacity constants
Abf = 4.34 * 10 ** 21
Aff = 3.68 * 10 ** 18
Aes = 0.02

# Energy generation/transfer constants
App = 0.241
ACNO = 8.67 * 10 ** 20
AHe = 50.9
gamma = 5 / 3  # Assume ideal monatomic gases
gamma_ratio = gamma / (gamma - 1)

# Solar parameters
Msun = 1.989 * 10 ** 30
Rsun = 6.957 * 10 ** 8
Lsun = 3.86 * 10 ** 26


class Shell:
    def __init__(self, m, Teff, Leff, X, Z):
        self.Rs = math.sqrt(Leff / (4 * math.pi * sigma * Teff ** 4))
        self.r = self.Rs  # Starting radius
        self.dm = - m / 10 ** 10
        self.Ms = m  # Total mass
        self.m = m  # Starting mass
        self.X = X
        self.Y = 1 - X - Z
        self.Z = Z
        self.P = 0
        self.T = 0
        self.L = Leff
        self.Ls = Leff  # Referenced for progress
        self.mu = (2 * self.X + 0.75 * self.Y + 0.5 * self.Z) ** -1  # Total ionization assumed
        self.rho = 0
        self.epsilon = 0
        self.kappa = 0
        self.i = 0
        self.convective = False

        # Profile initializations
        self.mprofile = [m]
        self.rprofile = [self.Rs]
        self.rhoprofile = [self.rho]
        self.kappaprofile = [0]
        self.Tprofile = [0]
        self.Lprofile = [Leff]
        self.Pprofile = [0]
        self.epsilonprofile = [0]
        self.iprofile = [0]
        self.convectiveprofile = [True]

        # Zone variables
        self.surface = True
        self.near_surface = True
        self.center = False
        self.complete = False

# ZONE MANAGEMENT, HALTING PROCEDURES, & DATA MANIPULATION

    def zone_progress(self):
        """Updates step sizes and dictates the equations to use for steps"""
        if self.i == 10:
            self.surface = False
        elif self.i == 1000:
            self.near_surface = False
            self.dm = self.dm * 10 ** 3
        elif self.i == 10000:
            self.dm = self.dm * 10 ** 3

    def completeness(self):
        """Checks for completed models, either by failure or reaching the center"""
        if self.T < 0:
            self.complete = True
            print(
                "Negative temperature found: Terminating calculation\nPlease try again with different parameters\n")
            print("The last shell was number " + str(self.i))
            print("The data has been written to profile.csv")
        elif self.P < 0:
            self.complete = True
            print("Negative pressure found: Terminating calculation\nPlease try again with different parameters\n")
            print("The last shell was number " + str(self.i))
            print("The data has been written to profile.csv")
        elif self.rho < 0:
            self.complete = True
            print("Negative density found: Terminating calculation\nPlease try again with different parameters\n")
            print("The last shell was number " + str(self.i))
            print("The data has been written to profile.csv")
        elif self.L < 0:
            self.complete = True
            print(
                "Negative luminosity found: Terminating calculation\nPlease try again with different parameters\n")
            print("The last shell was number " + str(self.i))
            print("The data has been written to profile.csv")
        elif self.m < self.Ms * 0.05 and (self.r > self.Rs * 0.3 or self.L > self.Ls * 0.3):
            self.complete = True
            print("Too much mass and/or luminosity is generated in the center\n"
                  "Please try again with different parameters\n")
            print("The last shell was number " + str(self.i))
            print("The data has been written to profile.csv")
        elif self.m < self.Ms * 0.1 and self.r < self.Rs * 0.1 and self.L < self.Ls * 0.1:
            self.center = True
            print("Success! You made it close enough to the center before stuff got weird\n")
            print("The last shell was number " + str(self.i))
            print("The data has been written to profile.csv")
            # Add extrapolation to the center here? Otherwise don't complete and add new else to step

    def update_profiles(self):
        """Adds latest value to profiles"""
        self.iprofile.append(self.i)
        self.mprofile.append(self.m)
        self.rprofile.append(self.r)
        self.Pprofile.append(self.P)
        self.Tprofile.append(self.T)
        self.Lprofile.append(self.L)
        self.kappaprofile.append(self.kappa)
        self.epsilonprofile.append(self.epsilon)
        self.rhoprofile.append(self.rho)
        self.convectiveprofile.append(self.convective)

    def plot(self):
        """Built-in procedure to plot mass, temperature, pressure, density, luminosity, and nuclear energy profiles"""
        fig, ax = plt.subplots(3, 2)

        ax[0, 0].scatter(self.rprofile, [i / mass for i in self.mprofile], c='lightgrey')
        ax[0, 0].set_xlabel('Radius (m)');
        ax[0, 0].set_ylabel('Interior Mass (total)')

        ax[0, 1].scatter(self.rprofile, self.Tprofile, c='firebrick')
        ax[0, 1].set_xlabel('Radius (m)');
        ax[0, 1].set_ylabel('Temperature (K)')

        ax[1, 0].scatter(self.rprofile, self.Pprofile, c='olivedrab')
        ax[1, 0].set_xlabel('Radius (m)');
        ax[1, 0].set_ylabel('Pressure (Pa)')

        ax[1, 1].scatter(self.rprofile, self.rhoprofile, c='blueviolet')
        ax[1, 1].set_xlabel('Radius (m)');
        ax[1, 1].set_ylabel('Density (kg m^-3)')

        ax[2, 0].scatter(self.rprofile, self.Lprofile, c='gold')
        ax[2, 0].set_xlabel('Radius (m)');
        ax[2, 0].set_ylabel('Interior Luminosity (W)')

        ax[2, 1].scatter(self.rprofile, self.epsilonprofile, c='skyblue')
        ax[2, 1].set_xlabel('Radius (m)');
        ax[2, 1].set_ylabel('Energy Generation (W m^3 kg^-2)')

        fig.set_size_inches(9, 9)
        fig.suptitle("Stellar Profiles")
        fig.tight_layout()
        plt.show()

    def write_df(self):
        df = pd.DataFrame(list(zip(self.iprofile, self.mprofile, self.rprofile, self.rhoprofile, self.Pprofile,
                                   self.Tprofile, self.epsilonprofile, self.kappaprofile, self.convectiveprofile)),
                          columns=['i', 'Mass (kg)', 'Radius (m)', 'Density (kg m^-3)', 'Pressure (kg m^-1 s^-2)',
                                   'Temperature (K)', 'Energy Generation (W m^3 kg^-2)', 'Opacity (m^2 kg^-1)',
                                   'Convective Transport'])
        df.set_index('i', inplace=True)
        df.to_csv('profile.csv')

# STATE EQUATIONS

    def calc_rho(self, radiation):
        """Calculate density. If radiation is true, considers radiation pressure, otherwise it does not."""
        if radiation:
            return (self.P - self.T ** 4 * a / 3) * self.mu * mH / (kB * self.T)
        else:
            return self.P * self.mu * mH / (kB * self.T)

    def calc_mu(self):
        """Calculate mean molecular weight (assuming total ionization)"""
        return (2 * self.X + 0.75 * self.Y + 0.5 * self.Z) ** -1

    def calc_kappa(self):
        """Calculate opacity"""
        bf_guillotine = 0.708 * (self.rho * (1 + self.X)) ** 0.2
        kappa_bf = (Abf / bf_guillotine) * self.Z * (1 + self.X) * self.rho / self.T ** 3.5
        kappa_ff = Aff * (1 - self.Z) * (1 + self.X) * self.rho / self.T ** 3.5
        kappa_es = Aes * (1 + self.X)
        return kappa_bf + kappa_ff + kappa_es

    def calc_epsilon(self):
        """Calculate nuclear energy generation rate"""
        T6 = self.T / 10 ** 6
        T8 = self.T / 10 ** 8

        psi_pp = 1 + 1.412 * 10 ** 8 * (self.X ** -1 - 1) * math.e ** (-49.98 * T6 ** (-1/3))
        Cpp = 1 + 0.0123 * T6 ** (1/3) + 0.0109 * T6 ** (2/3) + 0.000938 * T6
        fpp = 1
        epsilon_pp = App * self.rho * self.X ** 2 * fpp * psi_pp * Cpp * T6 ** (-2/3) * math.e ** (-33.8 * T6 ** (-1/3))

        CCNO = 1 + 0.0027 * T6 ** (1/3) - 0.00778 * T6 ** (2/3) - 0.000149 * T6
        epsilon_CNO = ACNO * self.rho * self.X * (self.Z / 2) * CCNO * T6 ** (-2/3) * math.e ** (-152.28 * T6 ** (-1/3))

        fHe = 1
        epsilon_He = AHe * self.rho ** 2 * self.Y ** 3 * T8 ** -3 * fHe * math.e ** (-44.027 / T8)

        return epsilon_pp + epsilon_CNO + epsilon_He

    def convection(self):
        """Returns boolean for method of heat transfer. Always call before updating profiles of P and T for indexing."""
        P1 = self.Pprofile[-1]
        T1 = self.Tprofile[-1]
        dlnPdlnT = self.T * (self.P - P1) / (self.P * (self.T - T1))
        if dlnPdlnT <= gamma_ratio:
            return True
        else:
            return False

# SURFACE FUNCTIONS (initial density 0 produces trivial solution, steps in radius ok at low density)

    def PT_surf_rad(self, A, r):
        """Calculates pressure and temperature for a radiative surface"""
        T_s = G * self.Ms * self.mu * mH * (r ** -1 - self.Rs ** -1) / (4.25 * kB)
        P_s = math.sqrt(16 * math.pi * G * self.Ms * a * c * kB / (4.25 * 3 * self.L * A * self.mu * mH)) * T_s ** 4.25
        return P_s, T_s

    def PT_surf_conv(self, dP, dT, r):
        """Calculates pressure and temperature for a convective surface"""
        T_s = G * self.Ms * gamma_ratio ** -1 * self.mu * mH * (r ** -1 - self.Rs ** -1) / kB
        P_s = (T_s * dP / dT) ** gamma_ratio
        return P_s, T_s

    def surf_step(self):
        """Calculates pressure, temperature, radius, and radius steps for the first several zones"""
        Abf_s = Abf * self.Z * (1 + self.X) / 0.01
        Aff_s = Aff * (1 - self.Z) * (1 + self.X)
        A = Abf_s + Aff_s
        dr = self.Rs / 1000  # Need to step far enough into the star for differential equations to work
        R = self.r
        r = R - dr
        P_1, T_1 = self.PT_surf_rad(A, r)
        r -= dr
        P_2, T_2 = self.PT_surf_rad(A, r)
        dlnPdlnT = T_2 * (P_2 - P_1) / (P_2 * (T_2 - T_1))

        if dlnPdlnT <= gamma_ratio:  # Convection dominates, recalculate
            self.convective = True
            P, T = self.PT_surf_conv(P_2 - P_1, T_2 - T_1, R - dr)
            return P, T, self.r - dr, dr
        else:  # Radiation dominates, return values
            self.convective = False
            return P_1, T_1, self.r - dr, dr

# CENTRAL TEMPERATURE PROCEDURE

    def f(self, T):
        """Temperature-pressure relation. Set to zero to find roots of temperature in Newton-Raphson method."""
        return self.rho * kB * T / (self.mu * mH) + a * T ** 4 / 3 - self.P

    def f_prime(self, T):
        """Derivative of temperature-pressure relation for Newton-Raphson method"""
        return self.rho * kB / (self.mu * mH) + 4 * a * T ** 3 / 3

    def nr_solve(self):
        """Implementation of Newton-Raphson method. Accuracy set to a convergence of 1/1000th the target temperature."""
        T1 = self.T
        T2 = 0
        while np.abs(T2 - T1) > self.T / 1000:
            f = self.f(T1)
            f_prime = self.f_prime(T1)
            T2 = T1 - f / f_prime
        return T2

# DIFFERENTIAL EQUATIONS

    def r_derivative(self, m, r):
        """Holds the differential equation for radius. Also used to convert between differentials dm and dr."""
        return (4 * math.pi * self.rho * r ** 2) ** -1

    def P_derivative(self, m, P):
        """Holds the differential equation for pressure"""
        return self.r_derivative(m, self.r) * (-G * m * self.rho) / self.r ** 2

    def L_derivative(self, m, L):
        """Holds the differential equation for luminosity"""
        return self.r_derivative(m, self.r) * (4 * math.pi * self.rho * self.r ** 2 * self.calc_epsilon())

    def T_derivative(self, m, T):
        """Holds the differential equations for temperature"""
        dTdr = - 3 * self.calc_kappa() * self.rho * self.L / (4 * a * c * T ** 3 * 4 * math.pi * self.r ** 2)
        if self.convective:
            return - (1 - gamma ** -1) * self.mu * mH * G * self.m / (kB * self.r ** 2)
        else:
            return self.r_derivative(m, self.r) * dTdr

# RUNGE-KUTTA AND NEWTON-RAPHSON ALGORITHMS

    def rk_step(self, var):
        """Contains the Runge-Kutta step procedures for various parameters. Requires Python 3.10 to use match"""
        match var:
            case "r":
                k1 = self.r_derivative(self.m, self.r)
                k2 = self.r_derivative(self.m + 0.5 * self.dm, self.r + 0.5 * k1 * self.dm)
                k3 = self.r_derivative(self.m + 0.5 * self.dm, self.r + 0.5 * k2 * self.dm)
                k4 = self.r_derivative(self.m + self.dm, self.r + k3 * self.dm)
                r_new = self.r + (self.dm / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
                return r_new
            case "P":
                k1 = self.P_derivative(self.m, self.P)
                k2 = self.P_derivative(self.m + 0.5 * self.dm, self.P + 0.5 * k1 * self.dm)
                k3 = self.P_derivative(self.m + 0.5 * self.dm, self.P + 0.5 * k2 * self.dm)
                k4 = self.P_derivative(self.m + self.dm, self.P + k3 * self.dm)
                P_new = self.P + (self.dm / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
                return P_new
            case "L":
                k1 = self.L_derivative(self.m, self.L)
                k2 = self.L_derivative(self.m + 0.5 * self.dm, self.L + 0.5 * k1 * self.dm)
                k3 = self.L_derivative(self.m + 0.5 * self.dm, self.L + 0.5 * k2 * self.dm)
                k4 = self.L_derivative(self.m + self.dm, self.L + k3 * self.dm)
                L_new = self.L + (self.dm / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
                return L_new
            case "T":
                k1 = self.T_derivative(self.m, self.T)
                k2 = self.T_derivative(self.m + 0.5 * self.dm, self.T + 0.5 * k1 * self.dm)
                k3 = self.T_derivative(self.m + 0.5 * self.dm, self.T + 0.5 * k2 * self.dm)
                k4 = self.T_derivative(self.m + self.dm, self.T + k3 * self.dm)
                T_new = self.T + (self.dm / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
                return T_new
            case _:
                return 0

# STEP INSTRUCTIONS

    def step(self):
        """Instructions for taking a step at various regions"""
        if self.surface:  # Calculate surface values
            self.P, self.T, self.r, dr = self.surf_step()  # Calculate P and T values
            # Convection boolean updated in surf_step function

            # Update rho, mass, luminosity, mu, epsilon
            self.rho = self.calc_rho(False)
            self.epsilon = self.calc_epsilon()
            self.kappa = self.calc_kappa()
            self.m += 4 * math.pi * self.r ** 2 * self.rho * dr
            self.L = self.rk_step('L')
            self.i += 1

        elif self.near_surface:  # Small step size since surface is highly variable
            self.m += self.dm
            self.i += 1
            self.rho = self.calc_rho(False)  # Pressure not high enough to accommodate radiation
            self.epsilon = self.calc_epsilon()
            self.kappa = self.calc_kappa()
            self.r = self.rk_step('r')
            self.T = self.rk_step('T')
            self.L = self.rk_step('L')
            self.P = self.rk_step('P')
            self.convective = self.convection()

        elif self.center:  # Center is technically complete but needs calculated before self.complete turns True
            self.rho = 3 * self.m / (4 * math.pi * self.r ** 3)  # Estimate central density from last shell
            self.P = self.P + 2 * math.pi * G * self.rho ** 2 * self.r ** 2 / 3
            self.T = self.nr_solve()
            self.epsilon = self.L / self.m
            self.kappa = self.calc_kappa()
            self.convective = self.convection()

            # Definition values of center
            self.L = 0
            self.m = 0
            self.r = 0
            self.i += 1
            self.complete = True

        elif not self.complete:  # Allowed higher step size now, include radiation pressure
            self.m += self.dm
            self.i += 1
            self.rho = self.calc_rho(True)
            self.epsilon = self.calc_epsilon()
            self.kappa = self.calc_kappa()
            self.r = self.rk_step('r')
            self.T = self.rk_step('T')
            self.L = self.rk_step('L')
            self.P = self.rk_step('P')
            self.convective = self.convection()

# RUN INTEGRATION

    def run(self):
        """Contains the functions that need to be run at each step in the iteration"""
        while not self.complete:
            self.step()  # Compute next shell parameters
            self.zone_progress()  # Update step sizes if needed
            self.completeness()  # Check for errors or completion
            if not self.complete or self.center:
                self.update_profiles()  # Add profiles if they survive error check


# For the sun
'''
mass = Msun
Teff = 5776
Leff = Lsun
hydrogen = 0.73
metal = 0.02
'''

mass = float(input("Enter the mass (units Solar mass): ")) * Msun
Teff = float(input("Enter the effective temperature (K): "))
Leff = float(input("Enter the trial luminosity (units Solar luminosity): ")) * Lsun
hydrogen = float(input("Enter the hydrogen mass fraction: "))
metal = float(input("Enter the metal mass fraction: "))

star = Shell(mass, Teff, Leff, hydrogen, metal)
star.run()
star.plot()
star.write_df()
