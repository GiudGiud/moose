#  External PETSc time stepper

This time stepper is used to query dt from PETSc TS to synchronize the time step size
between the wrapper app and the external PETSc solver.

The list of PETSc time steppers / integrators can be found on
[this summary](https://petsc.org/release/overview/integrator_table/#integrator-table).

This object expects the external PETSc solver to be using a PETSc TS (PETSc time stepper)
and retrieves the time step from this the TS.
