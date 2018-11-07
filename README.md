# EquilibriumMagnetopause

This repository contains instructions to find the position of an equilibrium magnetopause by solving the local pressure balance between external and internal pressure sources. It contains 3 sub-directories:

  - Dipole_Analytical: the pressure balance considered is (SW dynamic pressure) = (Dipolar planetary magnetic pressure); the numerical scheme makes use of an analytically-derived Jacobian matrix. The code is complete and operational.

  - DipoleCAN_Analytical: the pressure balance considered is (SW dynamic pressure) = (Dipolar planetary magnetic pressure + Hot plasma pressure + CAN Disk equatorial currents); the numerical scheme makes use of an analytically-derived Jacobian matrix. The code is NOT complete and NOT operational.

  - DipoleCAN_Numerical: the pressure balance considered is (SW dynamic pressure) = (Dipolar planetary magnetic pressure + Hot plasma pressure + CAN Disk equatorial currents); the numerical scheme makes use of an numerically-approximated Jacobian matrix. The code is complete and operational.
