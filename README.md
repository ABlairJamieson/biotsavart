# biotsavart

Simple code to calculate the magnetic field in Hyper-Kamkokande.

To build the code you will need root, but otherwise just type "make" at the command line.

Biotsavart calculation classes include:
- WireElement --  to represent a step along a wire
- typedef std::vector< WireElement > Wire -- to represent a wire
- BiotSavart -- to get magnetic field -- initialized with Wire and current

Classes specific to the calculation in HyperK include:
- HKEndCapCoils -- to build coils of wire on endcaps
- HKVerticalCoils -- to build circular coils to compensate B along z-axis
- HKHorizontalCoils -- to build box coils to compensate B along x-axis

Main code calls functions to:
- plot the magnetic field and save plots to rootfile
- do a chi^2 fit using minuit (commented out)
- do a markov-chain monte-carlo  (commented out)



