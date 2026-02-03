# tr-XPES
  This package calculates time-resolved X-ray Photoelectron Spectroscopy
(tr-XPES) signals.

  It is written in Python and Matlab.

  The quantum chemistry is expected to be done with NWChem.
 
  In the quantum chemistry, the cation excitations are described by
using an energy window such that we get transitions from the
corresponding core orbital to the valence range. These transitions are
equivalent to the photoemission from the neutral excited state.
Therefore we need to locate the correct excited state in the cation by
mapping onto the excited state of the nuetral molecule. Aditionally,
for the calculation of the Dyson orbitals, we need to swap the core
orbital index with the one of the orbital replaced in the main
transition of the neutral excited state (HOMO in the original case).

  Transition dipole moments and rotational averages are coded in
Matlab.

  The example corresponds to a single point computation of the HBQ
system.

  The notebook splits the computation in two steps:
  1- DO computation
  2- Transition moment computation

  From the transition moments the computation of the tr-XPES signal is
strighforward.

# Mayer Bond Order

  Functionality for calculating MBO is also included.