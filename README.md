# LEM (Large Expressions Management)

This is a module for the `LEM` (Large Expressions Management) package. It contains the functions to deal with large expressions management. The module can be used to veil and unveil large expressions to avoid expression swell. The module is very similar to the `LargeExpressions` module provided by Maple. The main difference is that the 'LEM' module has some additional built-in functions to display, list and substitute veiled expressions.

The code in this repository is hopefully an improved version of the code provided in Wenqin Zhou's PhD thesis *Symbolic Computation Techniques for Solving Large Expressions*.

We would like to thank *Jacques Carette* for providing the original code that we used to develop this module.

## Installation

To install the package you must have first installed Maple. Then open the `PackAndGo.mw` file and use the `!!!` button to *execute the entire worksheet*.

Then test the library in a Maple worksheet or document by typing:
```
> with(LEM);
```
Alternatively, you can use one of the tests file provided in the `tests` folder. If the package is loaded without errors, it is done!

## Package description

If you want a full description of the `LEM` package type:
```
> Describe(LEM);
```
This command will generate a brief description of the module and all the procedures and other objects present in the `LEM.mpl` file, which will be (very) similar to the following code.
```
# Large Expressions Management module.
module LEM:

# Veil an expression <x> and return a label to it.
Veil(
  x::{anything},
  $)::{anything}

# Unveil the expression <x>.
UnVeil(
  x::{anything},
  $)::{anything}

# Unveil the expression <x> with internal permutation map.
UnVeilImap(
  x::{anything},
  $)::{anything}

# Return a list of the veiling variables labelled as <label>.
VeilUnorderedList(
  label::{symbol},
  $)::{anything}

# Return a list of the veiling variables labelled as <label>.
# If <label> is not given, return a list of all veiling variables.
VeilList(
  label::{symbol, list(symbol)} := _VeilLabels(),
  reverse_order::{boolean}      := false,
  $)::{list(anything)}

# Return the size of the table for symbol <label>.
VeilTableSize(
  label::{symbol},
  $)::{nonnegint}

# Return the table for symbol <label> and the permutation that sorts it.
VeilTableImap(
  label::{symbol},
  reverse_order::{boolean},
  $)::{table, list[nonnegint]}

# Append the veiled expression <x> to the veiling table with symbol
# <label>.
VeilTableAppend(
  label::{symbol},
  x::{anything},
  $)::{anything}

# Return a list of the veiling labels.
VeilLabels($)

# Substitute the reversed veiling variables of the veiling label <label> in
# the expression <x>. If <label> is not given, substitute the reversed veiling
# variables of all veiling labels.
VeilSubs(
  x::{anything},
  label::{symbol, list(symbol)} := _VeilLabels(),
  $)::{anything}

# Clear all the veiling variables of the veiling label <label>. If <label>
# is not given, clear all the veiling variables.
VeilForget(
  label::{symbol, list(symbol)} := _VeilLabels(),
  $)::{nothing}
```

## Worked example

In case you have no time to read the description and realize how it should or should not work, here is a simple worked example.

```
> restart:
> with(LEM);

Write some random polynomial
> p := randpoly([x,y,z], degree = 5, dense);
> LEM:-VeilList();

Veil the long expressions with veiling variable 'X'
> p_X := collect(p, x, LEM:-Veil[X]);

Veil the long expressions with veiling variable 'Y'
> p_Y := collect(p, y, LEM:-Veil[Y]);

Veil the long expressions with veiling variable 'Z'
> p_Z := collect(p, z, LEM:-Veil[Z]);

Get the list of veiling variables
> LEM:-VeilLabels();

Get the list of veiling variables 'X'
> LEM:-VeilList(X);

Get the list of veiling variables 'Y'
> LEM:-VeilList(Y);

Get the list of veiling variables 'Z'
> LEM:-VeilList(Z);

List all veiling variables
> LEM:-VeilList();

Substitute the veiling variables 'X' in the polynomial 'p_X'
> simplify(LEM:-VeilSubs(p_X, X) - p);

Substitute the veiling variables 'Y' in the polynomial 'p_Y'
> simplify(LEM:-VeilSubs(p_Y,Y) - p);

Substitute the veiling variables 'I' and 'J' in the polynomial 'p_X+p_Y'
> simplify(LEM:-VeilSubs(p_X + p_Y, [X,Y]) - 2*p);

Substitute all the veiling variables in the polynomial 'p_X+p_Y+p_Z'
> simplify(LEM:-VeilSubs(p_X + p_Y + p_Z) - 3*p);

Substitute all the veiling variables in the polynomials 'p_X' and 'p_Y'
> simplify(LEM:-VeilSubs(p_X, X) - LEM:-VeilSubs(p_Y, Y));

Get the 'X' veiling lists and number of veiling variables
> LEM:-VeilTableSize(X);
> LEM:-VeilList(X);
> LEM:-VeilUnorderedList(X);
> LEM:-VeilTableImap(X, reverse = true); # or reverse = false

Get the 'Y' veiling lists and number of veiling variables
> LEM:-VeilTableSize(Y);
> LEM:-VeilList(Y);
> LEM:-VeilUnorderedList(Y);
> LEM:-VeilTableImap(Y, reverse = true); # or reverse = false

Get the 'Z' veiling number of veiling variables, lists and table
> LEM:-VeilTableSize(Z);
> LEM:-VeilList(Z);
> LEM:-VeilUnorderedList(Z);
> LEM:-VeilTableImap(Z, reverse = true); # or reverse = false

Forget the veiling variables 'X'
> LEM:-VeilForget(X);
> LEM:-VeilList();

Forget the veiling variables 'Y'
> LEM:-VeilForget(Y);
> LEM:-VeilList();

Forget all the remaining veiling variables
> LEM:-VeilForget();

List all veiling variables
> LEM:-VeilList();

Applend a new veiling variable
> LEM:-VeilTableAppend(Q, 2*a*b*c);
> LEM:-VeilList();
```

## Authors

### Current version

- *Davide Stocco* \
  Department of Industrial Engineering \
  University of Trento \
  email: davide.stocco@unitn.it

- *Matteo Larcher* \
  Department of Industrial Engineering \
  University of Trento \
  email: matteo.larcher@unitn.it

- *Enrico Bertolazzi* \
  Department of Industrial Engineering \
  University of Trento
  email: enrico.bertolazzi@unitn.it

### Original code

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
