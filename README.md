# Computational-Electromagnetism

#########  Numerical study of 1D wave equation  ##########
• Write the appropriate code to study the one-dimensional wave equation partialdee^2 (u)/ partialdee(x^2)=(1/c2) partialdee^2 (u) / partialdee (t^2)
with an explicit and an implicit numerical scheme of second-order
accuracy.
• Assume that c =1.
• Assume a rectangular pulse which extends to 10 space steps, as the initial
condition, i.e., u_1^(-1) = 1, i = 2...11 or 0 elsewhere
• Use Dirichlet boundary conditions at the two ends of the propagation axis.
• Record and show snapshots of the travelling wave for the first 50 time steps, every 10 time-steps (i.e., produce figures of u_i^n for n = 20, 30, 40, 50, 60)
o for each of the two methods (explicit/implicit);
o for Δt = 0.9 Δx/c, Δt = Δx/c, and Δt = 1.1Δx/c.
• Comment the shape of the pulse in each of the six cases (2 methods × 3 time-steps) and try to explain it.


##########  Study of numerical dispersion in 2D FDTD  ##########
• Write the appropriate code to study dispersion in the 2D implementation of FDTD by employing Newton’s iterative method to estimate the numerical phase velocity (eq. 4.16a in A. Taflove and S.C Hagness, Computational Electrodynamics. The Finite-Difference Time-Domain Method, 2nd edition, Artech House Inc., Norwood, 2000).
• Create a graph showing the ratio of numerical to analytical phase velocity as a function of the propagation angle for a uniform grid (i.e., same grid step in both axes) and at least three different grid steps (λ/5, λ/10, λ/20).
• Create a graph showing the ratio of numerical to analytical phase velocity as a function of the propagation angle for a non-uniform grid with the grid step along the y-axis half of that in the x-axis and for at least two different grid steps (λ/5, λ/10) of the y-axis.
• Assume the propagation angle to be zero for the propagation along the x-axis.


##########  Total-field/Scattered-field implementation in 1D FDTD  ##########
• Write a 1D FDTD code that implements the total-field/scattered-field (TFSF) technique for a plane-wave source.
• At the boundaries of your 1D computational domain use Mur’s ABC of first order.
[Hint: You can use a TFSF only on the left side of the 1D domain, if you wish.]
• Use the code to show that the skin depth of a conductor is inversely proportional to (i) the frequency of the plane-wave impinging on it, and (ii) the conductivity of a conducting scatterer located in the total field region.
[Hint: Keep the plane-wave frequency (scatterer conductivity) constant and change the scatterer conductivity (or plane-wave frequency); note how the skin depth changes. Use at least three points to plot and fit your data for each of the two cases. Restrict the study to microwave frequencies.]
