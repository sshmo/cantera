"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.

Requires: cantera >= 2.5.0
"""

import cantera as ct
import numpy as np

class PorousFlame(ct.FlameBase):
    """A freely-propagating flat flame."""
    # __slots__ = ('inlet', 'flame', 'outlet')
    # _other = ('grid', 'velocity')

    def __init__(self, gas, grid=None, width=None):
        """
        A domain of type PorousFlow named 'flame' will be created to represent
        the flame and set to axisymmetric flow. The three domains comprising the stack
        are stored as ``self.inlet``, ``self.flame``, and ``self.outlet``.

        :param grid:
            A list of points to be used as the initial grid. Not recommended
            unless solving only on a fixed grid; Use the `width` parameter
            instead.
        :param width:
            Defines a grid on the interval [0, width] with internal points
            determined automatically by the solver.
        """
        self.inlet = ct.Inlet1D(name='reactants', phase=gas)
        self.outlet = ct.Outlet1D(name='products', phase=gas)
        
        # Create flame domain if not already instantiated by a child class
        self.flame = ct.PorousFlow(gas, name='flame')
        self.flame.set_axisymmetric_flow()

        if width is not None:
            grid = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]) * width
        
        super().__init__((self.inlet, self.flame, self.outlet), gas, grid)

        # Setting X needs to be deferred until linked to the flow domain
        
        self.inlet.T = gas.T
        self.inlet.X = gas.X

    def set_initial_guess(self, locs=[0.0, 0.3, 0.7, 1.0], data=None, group=None):
        """
        Set the initial guess for the solution. By default, the adiabatic flame
        temperature and equilibrium composition are computed for the inlet gas
        composition. Alternatively, a previously calculated result can be
        supplied as an initial guess via 'data' and 'key' inputs (see
        `FlameBase.set_initial_guess`).

        :param locs:
            A list of four locations to define the temperature and mass fraction
            profiles. Profiles rise linearly between the second and third
            location. Locations are given as a fraction of the entire domain
        """
        super(PorousFlame, self).set_initial_guess()
        # set fixed temperature

        self.gas.TPY = self.inlet.T, self.P, self.inlet.Y

        if not self.inlet.mdot:
            # nonzero initial guess increases likelihood of convergence
            self.inlet.mdot = 1.0 * self.gas.density

        Y0 = self.inlet.Y
        u0 = self.inlet.mdot / self.gas.density
        T0 = self.inlet.T

        # get adiabatic flame temperature and composition
        self.gas.equilibrate('HP')
        Teq = self.gas.T
        Yeq = self.gas.Y
        u1 = self.inlet.mdot / self.gas.density

        self.set_profile('velocity', locs, [u0, u0, u1, u1])
        self.set_profile('T', locs, [T0, T0, Teq, Teq])

        # Pick the location of the fixed temperature point, using an existing
        # point if a reasonable choice exists
        Tmid = 0.5 * T0 + 0.5 * Teq
        self.set_fixed_temperature = Tmid

        for n in range(self.gas.n_species):
            self.set_profile(self.gas.species_name(n),
                             locs, [Y0[n], Y0[n], Yeq[n], Yeq[n]])
    def solve(self, loglevel=1, refine_grid=True, auto=False):
        """
        Solve the problem.

        :param loglevel:
            integer flag controlling the amount of diagnostic output. Zero
            suppresses all output, and 5 produces very verbose output.
        :param refine_grid:
            if True, enable grid refinement.
        :param auto: if True, sequentially execute the different solution stages
            and attempt to automatically recover from errors. Attempts to first
            solve on the initial grid with energy enabled. If that does not
            succeed, a fixed-temperature solution will be tried followed by
            enabling the energy equation, and then with grid refinement enabled.
            If non-default tolerances have been specified or multicomponent
            transport is enabled, an additional solution using these options
            will be calculated.
        """
        if not auto:
            return super().solve(loglevel, refine_grid, auto)

        # Use a callback function to check that the domain is actually wide
        # enough to contain the flame after each steady-state solve. If the user
        # provided a callback, store this so it can called in addition to our
        # callback, and restored at the end.
        original_callback = self._steady_callback

        class DomainTooNarrow(Exception): pass

        def check_width(t):
            T = self.T
            x = self.grid
            mRef = (T[-1] - T[0]) / (x[-1] - x[0])
            mLeft = (T[1] - T[0]) / (x[1] - x[0]) / mRef
            mRight = (T[-3] - T[-1]) / (x[-3] - x[-1]) / mRef

            # The domain is considered too narrow if gradient at the left or
            # right edge is significant, compared to the average gradient across
            # the domain.
            if mLeft > 0.02 or mRight > 0.02:
                raise DomainTooNarrow()

            if original_callback:
                return original_callback(t)
            else:
                return 0.0

        self.set_steady_callback(check_width)

        for _ in range(100):
            try:
                super().solve(loglevel, refine_grid, auto)
                break
            except DomainTooNarrow:
                self.flame.grid *= 2
                if loglevel > 0:
                    print('Expanding domain to accommodate flame thickness. '
                          'New width: {} m'.format(
                          self.flame.grid[-1] - self.flame.grid[0]))
                if refine_grid:
                    self.refine(loglevel)

        self.set_steady_callback(original_callback)


# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
width = 0.2  # m
loglevel = 1  # amount of diagnostic output (0 to 8)
mdot = 0.4271
phi = 2.5
# Solution object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture

gas = ct.Solution('gri30.yaml')
gas.set_equivalence_ratio(phi, 'CH4', {'O2': 0.9, 'N2': 0.1 })
gas.TP = Tin, p

# Set up flame object
f = PorousFlame(gas, width=width)


f.flame.pore1 = 0.7
f.flame.pore2 = 0.7
f.flame.diam1 = 0.0015
f.flame.diam2 = 0.0015
f.flame.Omega1 = 0.8
f.flame.Omega2 = 0.8
f.flame.srho=510
f.flame.sCp=824
f.flame.zmid=0.035
f.flame.dzmid =0.002



f.set_refine_criteria(ratio=4, slope=0.4, curve=0.4)
f.inlet.mdot = mdot
f.show_solution()

#f.flame.set_steady_tolerances(rel = 1.0e-4, abs=1.0e-9)
#f.flame.set_transient_tolerances(rel = 1.0e-4, abs=1.0e-9)


# Solve with mixture-averaged transport model
f.energy_enabled = True
f.transport_model = 'Mix'
f.solve(loglevel=loglevel, auto=False)


# Solve with the energy equation enabled
f.save('porous_flame.xml', 'Mix',
           'solution with mixture transport')

# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('porous_flame.csv', quiet=False)
