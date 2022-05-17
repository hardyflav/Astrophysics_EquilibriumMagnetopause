# EquilibriumMagnetopause

This repository contains instructions to find the position of an equilibrium magnetopause by solving the local pressure balance between external and internal pressure sources. It contains 3 sub-directories:

  - Dipole_Analytical: the pressure balance considered is (SW dynamic pressure) = (Dipolar planetary magnetic pressure); the numerical scheme makes use of an analytically-derived Jacobian matrix. The code is complete and operational.

  - DipoleCAN_Analytical: the pressure balance considered is (SW dynamic pressure) = (Dipolar planetary magnetic pressure + Hot plasma pressure + CAN Disk equatorial currents); the numerical scheme makes use of an analytically-derived Jacobian matrix. The code is NOT complete and NOT operational.

  - DipoleCAN_Numerical: the pressure balance considered is (SW dynamic pressure) = (Dipolar planetary magnetic pressure + Hot plasma pressure + CAN Disk equatorial currents); the numerical scheme makes use of an numerically-approximated Jacobian matrix. The code is complete and operational.

Related publications:
[Link to Google](https://www.google.com)

[Magnetopause Compressibility at Saturn With Internal Drivers](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GL086438)
Geophysical Research Letters
 
[The Magnetodisk Regions of Jupiter and Saturn](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/9781119815624.ch29)
American Geophysical Union

[A Self-Regulating Equilibrium Magnetopause Model with Applications to Saturn](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JA026751)
Journal of Geophysical Research
 
