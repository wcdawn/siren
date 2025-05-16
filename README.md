SIREN

Answering the Siren's call.

# Usage

To build.

```
git clone <url>
cd siren/src
make
```

Then, to run.

```
cd siren/cases/tworeg
../../src/siren.x ./twreg.inp
```

# Motivation

This is "to answer the siren's call" of investigating anisotropic scattering in hydrogeneous media.
It is my beginning of my response to Kord's "lost" hydrogen paper.

I realized that a lot of "theory" for anisotropic scattering is based on some hand-waving around the one-dimensional neutron transport equation, specifically the PN (spherical harmonics) form.
However, I haven't seen many places where these derivations have been demonstrated numerically.
What if we're missing something?
Or what if there's more that we can learn?

What SIREN is:
- A testbed for playing with the one-dimensional, multigroup PN equations.
- A jumping-off point for investigating neutron scattering in hydrogeneous media.

What SIREN is not:
- A general transport method. It works for one-dimensional, fixed meshes with constant spacing.
- Multiphysics. I've done it enough times, I don't want to.

# FAQ

1. Why Fortran?

I know Fortran. So do my coworkers. It allows for better collaboration.

1. Can I add XYX?

Yep. Make a pull request.

1. Will you add XYX?

Probably not.
