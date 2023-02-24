# LULEM

This is a module for the `LULEM` (LU decomposition with Large Expression Management) package. It contains the functions to solve systems of linear equations with large symbolic expressions. The module uses the LU decomposition to solve the system of equations. It is hopefully an improved version of the code provided in the Wenqin Zhou's PhD thesis *Symbolic Computation Techniques for Solving Large Expressions*.

We would like to thank *Jacques Carette* for providing the original code that we used to develop this module.

## Installation

To install the package you must have first installed Maple. Then copy the released MLA (Maple Library Archive) file `LULEM.mla` in the toolbox folder of Maple installation, which should be:

- OSX: `/Library/Frameworks/Maple.framework/Versions/Current/toolbox`;
- Windows: `C:/Programs/Maple/toolbox/`;
- Linux: `???` (if you managed to install Linux probably you know better than me where is the right folder 🫡).

If the toolbox folder does not exist, create it.

Then load the library in a Maple worksheet or document by typing:
```
> with(LULEM);
```
If the package is loaded without errors, it is done!

## Usage

If you want a full description of the `LULEM` package type:
```
> Describe(LULEM);
```
This command will generate a brief description of the module and all the procedures and other objects present in the `LULEM.mpl` file.

In case you are lazy as I am here is a simple worked example.

```
> restart;
> with(LULEM);

Generate a linear system of the type Ax=B
> d := 3;
> A := Matrix(d, d, symbol = a);
> B := Vector(d, symbol = b);

Set A[1,1] element equal to zero to force permutation (optional)
> A[1,1] := 0;

Perform the LU decomposition
> P, L, U, r := LUD(A, Q, VeilingStrategy, PivotStrategy, ZeroStrategy);

Check the LU decomposition
> simplify(SubsVeil(Q, P.A -~ L.U));

Let us see how many Q veiling variables we have
> LastUsed[Q];

Let us see the Q veiling variables
> ShowVeil(Q);

But now we want to forget about the Q veiling variables
> ForgetVeil(Q);

Now let us find the solution of the linear system
> sol := Solve(A, B, K, VeilingStrategy, PivotStrategy, ZeroStrategy);

Check the solution
> simplify(SubsVeil(K, A.sol-B)):

Let us see how many K veiling variables we have
> LastUsed[K];

Let us see the K veiling variables
> ShowVeil(K);

But now we want to forget about the K veiling variables
> ForgetVeil(K);
```

Notice that the available veiling, pivoting and zero detection strategies should be chosen between:
```
VeilingStrategy -> VeilingStrategy_n
                   VeilingStrategy_L
                   VeilingStrategy_Ls
                   VeilingStrategy_LB

PivotStrategy   -> PivotStrategy_Llength
                   PivotStrategy_Slength
                   PivotStrategy_Lindets
                   PivotStrategy_Sindets
                   PivotStrategy_numeric

ZeroStrategy    -> ZeroStrategy_length
                   ZeroStrategy_normalizer
```

## Authors

- *Davide Stocco* (maintainer) \
  Department of Industrial Engineering \
  University of Trento \
  email: davide.stocco@unitn.it

- *Matteo Larcher* (maintainer) \
  Department of Industrial Engineering \
  University of Trento \
  email: matteo.larcher@unitn.it

- *Enrico Bertolazzi* \
  Department of Industrial Engineering \
  University of Trento

- *Wenqin Zhou* (former affiliation) \
  Department of Applied Mathematics \
  University of Western Ontario

- *David J. Jeffrey* \
  Department of Applied Mathematics \
  University of Western Ontario

- *Jacques Carette* \
  Department of Computing and Software \
  McMaster University

- *Robert M. Corless* \
  Department of Applied Mathematics \
  University of Western Ontario

## References

```
@phdthesis{zhou2007symbolic,
  author = {Zhou, Wenqin},
  title = {Symbolic Computation Techniques for Solving Large Expression Problems from Mathematics
    and Engineering},
  year = {2007},
  school = {University of Western Ontario}
}
```

```
inproceedings{zhou2006hierarchical,
  title = {Hierarchical representations with signatures for large expression management},
  author = {Zhou, Wenqin and Carette, Jacques and Jeffrey, David J. and Monagan, Michael B.},
  booktitle = {Artificial Intelligence and Symbolic Computation: 8th International Conference,
    AISC 2006 Beijing, China, September 20-22, 2006 Proceedings 8},
  pages = {254--268},
  year = {2006},
  organization = {Springer}
}
```

```
@inproceedings{carette2006linear,
  title = {Linear algebra using Maple’s LargeExpressions package},
  author = {Carette, Jacques and Zhou, Wenqin and Jeffrey, David J. and Monagan, Michael B.},
  booktitle = {Proceedings of Maple Conference},
  pages = {14--25},
  year = {2006}
}
```

*Documented by Davide Stocco*
