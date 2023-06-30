{
    const volScalarField rhoBlend
    (
        pos(alpha1 - dimensionedScalar("half", dimless, 0.5))
    );
    p = pd + (rhoBlend*rho1 + (1 - rhoBlend)*rho2)*gh;

    // Note: as p is constructed with "calculated" boundary conditions, there is
    // no need to call p.correctBoundaryConditions()
}
