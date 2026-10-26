# Two interacting scalar fields around a binary black hole

This example evolves two real scalar fields, `phi1` and `phi2`, minimally
coupled to gravity and interacting through their potential, in full GR around
a binary black hole. It is built on [Examples/BinaryBH](../BinaryBH) and keeps
all of its features: approximate boosted-BH or TwoPunctures initial data,
puncture tracking, the apparent horizon finder (when compiled with
`USE_AHFINDER`) and Weyl4 extraction.

## Physics setup

The matter action is

```
S = \int d^4x \sqrt{-g} [ - 1/2 (\nabla phi1)^2 - 1/2 (\nabla phi2)^2 - V(phi1, phi2) ]
```

with signature (-,+,+,+), so that each field obeys `Box phi_a = dV/dphi_a`.
The evolved momenta follow the `ScalarField` convention
`Pi_a = n^mu d_mu phi_a = (d_t phi_a - beta^i d_i phi_a) / lapse`.

The matter class is
[Source/Matter/TwoScalarFields.hpp](../../Source/Matter/TwoScalarFields.hpp).
It has the same interface as `ScalarField`, so it works as the `matter_t` of
`MatterCCZ4RHS`, `MatterConstraints`, `EMTensor` and `MatterWeyl4`. The level
class uses:

- `MatterCCZ4RHS` for the evolution;
- `MatterConstraints` for the constraint norms and plot files (`Ham`,
  `Mom1-3`);
- `MatterWeyl4` for the gravitational-wave extraction;
- `EMTensor` for the energy density `rho`, written to the plot files.

The data in `data/constraint_norms.dat` include the matter terms.

## The potential and how to change it

The default potential, in [Potential.hpp](./Potential.hpp), is

```
V = 1/2 mu1^2 phi1^2 + 1/2 mu2^2 phi2^2 + 1/2 lambda phi1^2 phi2^2
```

with parameters `scalar_mass1` (mu1), `scalar_mass2` (mu2) and `lambda`.
To use a different potential, edit `Potential::compute_potential`. It must set
`V`, `dV/dphi1` and `dV/dphi2`. If you add new parameters, add them to
`Potential::params_t` and read them in
[SimulationParameters.hpp](./SimulationParameters.hpp).

`SimulationParameters` warns when:

- an effective mass is not resolved on the coarsest level, using the same
  criterion as ScalarFieldBH: `mu_eff < 0.2 / coarsest_dx / dt_multiplier`,
  with `mu1_eff^2 = mu1^2 + lambda phi2_0^2` and
  `mu2_eff^2 = mu2^2 + lambda phi1_0^2`;
- `|K| * stop_time > 0.3`. The background expansion and the drift of the
  1+log lapse (`d_t lapse = -2 lapse K`) then matter, so `phi1_0` and
  `phi2_0` should be reduced;
- `lambda < 0`, since the potential is then unbounded below.

## Initial data

The vacuum binary BH metric is set first, using either `BinaryBH` or
TwoPunctures. The loop in
[InitialTwoFieldsData.hpp](./InitialTwoFieldsData.hpp) then runs separately.
The TwoPunctures loop cannot use SIMD, but this second loop can. It sets
uniform fields `phi1 = phi1_0`, `phi2 = phi2_0`, `Pi1 = Pi1_0` and
`Pi2 = Pi2_0` (the momenta default to 0), and sets

```
K = K_sign * sqrt(24 pi G rho),   rho = 1/2 Pi1_0^2 + 1/2 Pi2_0^2 + V(phi1_0, phi2_0).
```

Why this satisfies the constraints:

- Both vacuum data sets are conformally flat, with `K = 0` and a tracefree
  `A_ij`.
- The fields are uniform, so `S_i = 0` and `rho` is constant. `K` is then
  constant too.
- The momentum constraint `D_j A^ij - (2/3) D^i K = 8 pi G S^i` therefore
  reduces to the vacuum one.
- In the Hamiltonian constraint `R + (2/3) K^2 - A_ij A^ij = 16 pi G rho`, the
  term `(2/3) K^2` cancels `16 pi G rho` exactly. What remains is the vacuum
  constraint.

So TwoPunctures plus uniform fields satisfies the constraints to the accuracy
of the TwoPunctures solve. The approximate `BinaryBH` data has the same
constraint violation as the vacuum BinaryBH example.

`K_sign = -1` (the default) gives an expanding background. GRChombo's
convention is `K < 0` for expansion; see
[Examples/ScalarFieldCosmo/InitialK.hpp](../ScalarFieldCosmo/InitialK.hpp).
Note that [Examples/ScalarFieldBH](../ScalarFieldBH) uses `K > 0`.

## Boundary conditions

The uniform fields source an FRW-like background. Far from the BHs, `chi`,
`K`, `Theta`, the lapse and the fields all evolve in time. Sommerfeld
conditions with fixed asymptotic values would therefore violate the
constraints at the outer boundary. In a test with negligible BH masses, the
Hamiltonian constraint norm grew from 1e-10 to O(1) in one coarse step.

The params files therefore use mixed boundary conditions (`4`):

- `chi K Theta lapse phi1 Pi1 phi2 Pi2` are extrapolated, with
  `extrapolation_order = 1`;
- all other variables use Sommerfeld conditions.

With these, the cosmological test stays at round-off level.

## Restarting from an initial condition solver (GRTresna)

`postRestart` computes the constraints on every level after a restart and
prints their norms. A two-field output of an initial condition solver such as
[GRTresna](https://github.com/GRTLCollaboration/GRTresna) can then be used as
the `restart_file`. It needs variables named `phi1 Pi1 phi2 Pi2`, as in
Examples/ScalarFieldBH. The plan is to generate non-uniform field
configurations around the binary this way, instead of the uniform data used
here.

## Parameter files

- [params.txt](./params.txt): the production BinaryBH setup with
  `mu1 = 0.05`, `mu2 = 0.08`, `lambda = 0.1` and `phi1_0 = phi2_0 = 1e-4`.
  These are chosen so that the masses are resolved on the coarsest level and
  `|K| * stop_time ~ 0.13`.
- [params_two_punctures.txt](./params_two_punctures.txt): the same with
  TwoPunctures initial data.
- [params_very_cheap.txt](./params_very_cheap.txt): a small test run. It uses
  larger masses and amplitudes, so the matter has a visible effect over its
  short `stop_time`.

## Using TwoPunctures initial data

To build this example with [TwoPunctures](https://github.com/GRChombo/TwoPunctures)
initial data, set the environment variable `TWOPUNCTURES_SOURCE` to the Source
subdirectory of the local clone of the [TwoPunctures repository](https://github.com/GRChombo/TwoPunctures),
e.g.
```bash
export TWOPUNCTURES_SOURCE=/path/to/TwoPunctures/Source
```
Then, use `make` to build the `all-tp` target using a command such as
```bash
make all-tp -j 4
```
TwoPunctures requires the GNU Scientific Library (`-lgsl`). The created
executable will have the filename
```
Main_TwoFieldsBinaryBH_TP.<config>.ex
```
To stop using TwoPunctures, undefine the `TWOPUNCTURES_SOURCE` environment
variable:
```bash
unset TWOPUNCTURES_SOURCE
```
Note that the parameter names for TwoPunctures initial data differ from those
of the vanilla example: see
[params_two_punctures.txt](./params_two_punctures.txt).
