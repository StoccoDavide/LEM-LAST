# LEM (Large Expressions Management)

This is a module for the `LEM` (Large Expressions Management) package. It contains the functions to deal with large expressions management. The module can be used to veil and unveil large expressions to avoid expression swell.

The module is very similar to the 'LargeExpressions' module provided by Maple. The main difference is that the 'LEM' module has some additional built-in functions to display, list and substitute veiled expressions.

The code in this repository hopefully is an improved version of the code provided in Wenqin Zhou's PhD thesis *Symbolic Computation Techniques for Solving Large Expressions*.

We would like to thank *Jacques Carette* for providing the original code that we used to develop this module.

## Installation

To install the package you must have first installed Maple. Then follow the instructions in one of the following sections.

### Use the precompiled MLA file

Copy the latest [released](https://github.com/StoccoDavide/LEM/releases) MLA (Maple Library Archive) file `LEM.mla` in the toolbox folder of Maple installation, which should be:

- OSX: `/Library/Frameworks/Maple.framework/Versions/Current/toolbox`;
- Windows: `C:/Programs/Maple/toolbox/`;
- Linux: `???` (if you managed to install Linux probably you know better than me where is the right folder 🫡).

If the toolbox folder does not exist, create it.

Then load the library in a Maple worksheet or document by typing:
```
> with(LEM);
```
Alternatively, you can use the `test.mw` file provided in the repository. If the package is loaded without errors, it is done!

### Compile the library manually

This installation option is maintained only for OSX users. It is intended to be used by experienced Maple users that can edit and run the `CompileLibrary.mw` according to their OS.

## Package description

If you want a full description of the `LEM` package type:
```
> Describe(LEM);
```
This command will generate a brief description of the module and all the procedures and other objects present in the `LEM.mpl` file, which will be (very) similar to the following code.
```
# Large Expressions Management module module LEM:

# Veil an expression <x> and return a label to it.
Veil(
  x::{anything},
  $)::{anything}

# UnVeil the expression <x> up to <n> levels.
UnVeil(
  x::{anything},
  n::{infinity, nonnegint},
  $)::{anything}

# Return a list of the veiling labels.
VeilLabels(
  $)::{list(symbol)}

# Return a list of the veiling variables labelled as <label>.
# If <label> is not given, return a list of all veiling variables.
ListVeil(
  label::{symbol, list(symbol)} := VeilLabels(),
  $)::{list(anything)}

# Substitute the reversed veiling variables of the veiling label <label> in the expression <x>.
# If <label> is not given, substitute the reversed veiling variables of all veiling labels.
SubsVeil(
  x::{anything},
  label::{symbol, list(symbol)} := VeilLabels(),
  $)::{anything}

# Clear all the veiling variables of the veiling label <label>.
# If <label> is not given, clear all the veiling variables.
ForgetVeil(
  label::{symbol, list(symbol)} := VeilLabels(),
  $)::{nothing}

package LastUsed:
```

## Worked example

In case you have no time to read the description and realize how it should or should not work, here is a simple worked example.

```
> restart:
> with(LEM);

Write some random polynomial
> p := randpoly([x,y,z], degree = 5, dense);

Veil the long expressions with veiling variable 'X'
> p_X := collect(p, x, Veil[X]);

Veil the long expressions with veiling variable 'Y'
> p_Y := collect(p, y, Veil[Y]);

Veil the long expressions with veiling variable 'Z'
> p_Z := collect(p, z, Veil[Z]);

Get the list of veiling variables
> VeilLabels();

Get the list of veiling variables 'X'
> ListVeil(X);

Get the list of veiling variables 'Y'
> ListVeil(Y);

Get the list of veiling variables 'Z'
> ListVeil(Z);

List all veiling variables
> ListVeil();

Substitute the veiling variables 'X' in the polynomial 'p_X'
> SubsVeil(p_X, X);

Substitute the veiling variables 'Y' in the polynomial 'p_Y'
> SubsVeil(p_Y, Y);

Substitute the veiling variables 'I' and 'J' in the polynomial 'p_X+p_Y'
> SubsVeil(p_X + p_Y, [X,Y]);

Substitute all the veiling variables in the polynomial 'p_X+p_Y+p_Z'
> SubsVeil(p_X + p_Y + p_Z);

Get the number of 'X' veiling variables
> LastUsed[X];

Get the number of 'Y' veiling variables
> LastUsed[Y];

Get the number of 'Z' veiling variables
> LastUsed[Z];

Forget the veiling variables 'X'
> ForgetVeil(X);

Forget the veiling variables 'Y'
> ForgetVeil(Y);

Forget all the remaining veiling variables
> ForgetVeil();

List all veiling variables
> ListVeil();
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
