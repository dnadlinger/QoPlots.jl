QuantumOptics.jl plots
======================

Plotting functions for density matrices, etc. from the QuantumOptics.jl
package.

Supported plot types:
 - `plot_dm()` plots a density matrix in a 2D plot inspired by Hinton diagrams,
 	where symbol size represents magnitude and colour represents phase.
 - `plot_bloch()` plots an array of 1-qubit density matrices on the Bloch sphere.

**Note:** QoPlots.jl requrires a new-ish version of the PyCall Python interface package
with support for natural property syntax (JuliaPy/PyCall.jl#517). Installing the latest
version from the default repository (`add PyPlot`) should work without further
intervention.
