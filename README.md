# LEM (Large Expressions Management)

This is a module for the `LEM` (Large Expressions Management) package. It contains the functions to deal with large expressions management. The module can be used to veil and unveil large expressions to avoid expression swell. The module is very similar to the `LargeExpressions` module provided by Maple. The main difference is that the `LEM` module has some additional built-in functions to display, list and substitute veiled expressions.

The code in this repository is hopefully an improved version of the code provided in Wenqin Zhou's PhD thesis *Symbolic Computation Techniques for Solving Large Expressions*.

We would like to thank *Jacques Carette* for providing the original code that we used to develop this module.

## Installation

To install the package you must have first installed Maple. Then open the `PackAndGo.mw` file and use the `!!!` button to *execute the entire worksheet*.

Then test the library in a Maple worksheet or document by typing:
```
> LEM:-Info(LEM);
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
object LEM :: LEM:

# Print 'LEM' module information.
Info( )

# 'LEM' module load procedure.
ModuleLoad( )

# 'LEM' module unload procedure.
ModuleUnload( )

# Copy the objects <proto> into <self>.
ModuleCopy( _self::LEM, proto::LEM, $ )

# Set the veiling label to <label>.
SetVeilingLabel( _self::LEM, label::{string, symbol}, $ )

# Return the veiling label.
GetVeilingLabel( _self::LEM, $ ) :: symbol

# Clear the veiling table.
ClearVeilingTable( _self::LEM, $ )

# Check if the veiling strategy is verified, if true veil the expression
# <x> and return a label to it.
Veil( _self::LEM, x::anything,
      { force::boolean := false }, $ ) :: anything

# Compute the cost of the expression <x>.
ExpressionCost( _self::LEM, x::anything, $ ) :: integer

# Evaluate the veiling strategy for the expression <x>.
VeilingStrategy( _self::LEM, x::algebraic, $ ) :: boolean

# Set the veiling strategy parameters: maximum veiling cost <maxcost>,
# subscripts cost weight parameter <subscripts>, assignments cost weight
# parameter <assignments>, additions cost weight parameter <additions>,
# multiplications cost weight parameter <multiplications>, divisions cost
# weight parameter <divisions>, and functions cost weight parameter
# <functions>.
SetVeilingStrategyPars( _self::LEM,
                        { additions::nonnegint := 1,
                          assignments::nonnegint := 0,
                          divisions::nonnegint := 3,
                          functions::nonnegint := 2,
                          maxcost::nonnegint := 15,
                          multiplications::nonnegint := 2,
                          subscripts::nonnegint := 0 }, $ ) :: nothing

# Unveil the expression <x>.
UnVeil( _self::LEM, x::anything, $ ) :: anything

# Unveil the expression <x> with internal permutation map.
UnVeilImap( _self::LEM, x::anything, $ ) :: anything

# Return a list of the veiling labels.
VeilList( _self::LEM,
          { reverse::boolean := false }, $ ) :: list(anything)

# Return the unordered veiling list.
VeilUnorderedList( _self::LEM, $ ) :: list

# Return the size of the internal veiling table.
VeilTableSize( _self::LEM, $ ) :: nonnegint

# Return the veiling list and the permutation that sorts it.
VeilTableImap( _self::LEM,
                { reverse::boolean := false }, $ )
              :: {list(anything = anything), []}

# Append the veiled expression <x> to the veiling table.
VeilTableAppend( _self::LEM, x::anything, $ ) :: indexed

# Substitute the reversed veiling variables of the internal veiling table
# in the expression <x>.
VeilSubs( _self::LEM, x::anything, $ ) :: anything

# Clear all the veiling variables of the internal veiling table.
VeilForget( _self::LEM, $ )
```

## Usage

In case you have no time to read the description and realize how it should or should not work, here is a simple worked example.

### 🚧 Attention! 🚧

Maple object-oriented programming features have slightly changed in 2021, which online documentation states:

> As of Maple 2021, if the method has a formal parameter named `_self`, references to its object's local or exported variables may be written without prefixing them. That is, `_self:-variable` may be written as just `variable`. Maple will add the `self:-` prefix internally when the method is simplified.

> As of Maple 2021, a message-passing form of method call can be used, which will automatically pass the object as an argument if the method has a formal parameter named `_self`.

> Another way to invoke a method, similar to that used in other object-oriented languages, was introduced in Maple 2021. In this form, the object name is qualified by the method name, `object_name:-method_name( argument )` just as it can be in the function mechanism described above. However, the *object can be omitted from the argument sequence*.

For further information please refer to this [link](https://fr.maplesoft.com/support/help/Maple/view.aspx?path=object/methods).

### Maple < 2021

Here is a worked example of how to use the `LEM` module for Maple versions **before** 2021.

```
> restart;
> kernelopt(assertlevel = 2);

Print 'LEM' module informations
> LEM:-Info(LEM);

Write some random polynomial
> p := randpoly([x,y,z], degree = 5, dense);

Create some instaces of 'LEM'
> LEM_X := Object(LEM);
> LEM_X:-SetVeilingLabel(LEM_X, 'X');
> LEM_X:-GetVeilingLabel(LEM_X);
> LEM_Y := Object(LEM);
> LEM_Y:-SetVeilingLabel(LEM_Y, 'Y');
> LEM_Y:-GetVeilingLabel(LEM_Y);
> LEM_Z := Object(LEM);
> LEM_Z:-SetVeilingLabel(LEM_Z, 'Z');
> LEM_Z:-GetVeilingLabel(LEM_Z);

Veil the long expressions with veiling variable 'X'
> p_X := collect(p, x, i -> LEM_X:-Veil(LEM_X, i));

Veil the long expressions with veiling variable 'Y'
> p_Y := collect(p, y, i -> LEM_Y:-Veil(LEM_Y, i));

Veil the long expressions with veiling variable 'Z'
> p_Z := collect(p, z, i -> LEM_Z:-Veil(LEM_Z, i));

Get the list of veiling variables 'X'
> LEM_X:-VeilList(LEM_X);

Get the list of veiling variables 'Y'
> LEM_Y:-VeilList(LEM_Y);

Get the list of veiling variables 'Z'
> LEM_Z:-VeilList(LEM_Z);

Substitute the veiling variables 'X' in the polynomial 'p_X'
> simplify(LEM_X:-VeilSubs(LEM_X, p_X) - p);

Substitute the veiling variables 'Y' in the polynomial 'p_Y'
> simplify(LEM_Y:-VeilSubs(LEM_Y, p_Y) - p);

Substitute the veiling variables 'X' and 'Y' in the polynomial 'p_X+p_Y'
> simplify(LEM_X:-VeilSubs(LEM_X, LEM_Y:-VeilSubs(LEM_Y, p_X + p_Y)) - 2*p);

Substitute all the veiling variables in the polynomial 'p_X+p_Y+p_Z'
> simplify(LEM_X:-VeilSubs(LEM_X, LEM_Y:-VeilSubs(LEM_Y, LEM_Z:-VeilSubs(LEM_Z, p_X + p_Y + p_Z) - 3*p)));

Get the 'X' veiling lists and number of veiling variables
> LEM_X:-VeilTableSize(LEM_X);
> LEM_X:-VeilList(LEM_X);
> LEM_X:-VeilUnorderedList(LEM_X);
> LEM_X:-VeilTableImap(LEM_X, reverse = true);

Get the 'Y' veiling lists and number of veiling variables
> LEM_Y:-VeilTableSize(LEM_Y);
> LEM_Y:-VeilList(LEM_Y);
> LEM_Y:-VeilUnorderedList(LEM_Y);
> LEM_Y:-VeilTableImap(LEM_Y, reverse = true);

Get the 'Z' veiling number of veiling variables, lists and table
> LEM_Z:-VeilTableSize(LEM_Z);
> LEM_Z:-VeilList(LEM_Z);
> LEM_Z:-VeilUnorderedList(LEM_Z);
> LEM_Z:-VeilTableImap(LEM_Z, reverse = true);

Forget the veiling variables 'X'
> LEM_X:-VeilForget(LEM_X);
> LEM_X:-VeilList(LEM_X);

Forget the veiling variables 'Y'
> LEM_Y:-VeilForget(LEM_Y);
> LEM_Y:-VeilList(LEM_Y);

Forget the veiling variables 'Z'
> LEM_Z:-VeilForget(LEM_Z);
> LEM_Z:-VeilList(LEM_Z);

Append a new veiling variable with random veiling label generator
> LEM:-VeilTableAppend(LEM, 2*a*b*c);
> LEM:-VeilList(LEM);
> LEM:-GetVeilingLabel(LEM);

Try to change the defult (random) veiling label to 'A' (error)
> # LEM:-SetVeilingLabel(LEM, 'A');

Clear the veiling table
> L := LEM:-VeilList(LEM);
> LEM:-ClearVeilingTable(LEM);

Try to change the defult (random) veiling label to 'A' (success)
> LEM:-SetVeilingLabel(LEM, 'A');
```

For further information open the file `Test_Maple2021minus.mw` in the `tests` folder.

### Maple >= 2021


Here is a worked example of how to use the `LEM` module for Maple 2021 and the versions **after** 2021.

```
> restart;
> kernelopt(assertlevel = 2);

Print 'LEM' module informations
> LEM:-Info(LEM);

Write some random polynomial
> p := randpoly([x,y,z], degree = 5, dense);

Create some instaces of 'LEM'
> LEM_X := Object(LEM);
> LEM_X:-SetVeilingLabel('X');
> LEM_X:-GetVeilingLabel();
> LEM_Y := Object(LEM);
> LEM_Y:-SetVeilingLabel('Y');
> LEM_Y:-GetVeilingLabel();
> LEM_Z := Object(LEM);
> LEM_Z:-SetVeilingLabel('Z');
> LEM_Z:-GetVeilingLabel();

Veil the long expressions with veiling variable 'X'
> p_X := collect(p, x, i -> LEM_X:-Veil(i));

Veil the long expressions with veiling variable 'Y'
> p_Y := collect(p, y, i -> LEM_Y:-Veil(i));

Veil the long expressions with veiling variable 'Z'
> p_Z := collect(p, z, i -> LEM_Z:-Veil(i));

Get the list of veiling variables 'X'
> LEM_X:-VeilList();

Get the list of veiling variables 'Y'
> LEM_Y:-VeilList();

Get the list of veiling variables 'Z'
> LEM_Z:-VeilList();

Substitute the veiling variables 'X' in the polynomial 'p_X'
> simplify(LEM_X:-VeilSubs(p_X) - p);

Substitute the veiling variables 'Y' in the polynomial 'p_Y'
> simplify(LEM_Y:-VeilSubs(p_Y) - p);

Substitute the veiling variables 'X' and 'Y' in the polynomial 'p_X+p_Y'
> simplify(LEM_X:-VeilSubs(LEM_Y:-VeilSubs(p_X + p_Y)) - 2*p);

Substitute all the veiling variables in the polynomial 'p_X+p_Y+p_Z'
> simplify(LEM_X:-VeilSubs(LEM_Y:-VeilSubs(LEM_Z:-VeilSubs(p_X + p_Y + p_Z) - 3*p)));

Get the 'X' veiling lists and number of veiling variables
> LEM_X:-VeilTableSize();
> LEM_X:-VeilList();
> LEM_X:-VeilUnorderedList();
> LEM_X:-VeilTableImap(reverse = true);

Get the 'Y' veiling lists and number of veiling variables
> LEM_Y:-VeilTableSize();
> LEM_Y:-VeilList();
> LEM_Y:-VeilUnorderedList();
> LEM_Y:-VeilTableImap(LEM_Y, reverse = true);

Get the 'Z' veiling number of veiling variables, lists and table
> LEM_Z:-VeilTableSize();
> LEM_Z:-VeilList();
> LEM_Z:-VeilUnorderedList();
> LEM_Z:-VeilTableImap(reverse = true);

Forget the veiling variables 'X'
> LEM_X:-VeilForget();
> LEM_X:-VeilList();

Forget the veiling variables 'Y'
> LEM_Y:-VeilForget();
> LEM_Y:-VeilList();

Forget the veiling variables 'Z'
> LEM_Z:-VeilForget();
> LEM_Z:-VeilList();

Append a new veiling variable with random veiling label generator
> LEM:-VeilTableAppend(2*a*b*c);
> LEM:-VeilList();
> LEM:-GetVeilingLabel();

Try to change the defult (random) veiling label to 'A' (error)
> # LEM:-SetVeilingLabel('A');

Clear the veiling table
> L := LEM:-VeilList();
> LEM:-ClearVeilingTable();

Try to change the defult (random) veiling label to 'A' (success)
> LEM:-SetVeilingLabel('A');
```

For further information open the file `Test_Maple2021plus.mw` in the `tests` folder.

## Authors

### Current version authors:

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

### Inspired by the work of:

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
