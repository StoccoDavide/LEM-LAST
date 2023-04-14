# LEM (Large Expressions Management)

This is a module for the `LEM` (Large Expressions Management) module. It contains the functions to deal with large expressions management. The module can be used to veil and unveil large expressions to avoid expression swell.

The code in this repository is hopefully an improved version of the code provided in Wenqin Zhou's PhD thesis *Symbolic Computation Techniques for Solving Large Expressions*.

We would like to thank *Jacques Carette* for providing the original code that we used to develop this module.

## Installation

To install the module you must have first installed Maple. Then open the `PackAndGo.mw` file and use the `!!!` button to *execute the entire worksheet*.

Then test the module in a Maple worksheet or document by typing:

```
> LEM:-Info(LEM); # For Maple versions up to 2020
> LEM:-Info();    # For Maple versions starting from 2021
```

Alternatively, you can use one of the tests file provided in the `tests` folder. If the module is loaded without errors, it is done!

## Module description

If you want a full description of the `LEM` module type:

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

    # Enable the expression signature calculation.
    EnableSignature( _self::LEM, $ )

    # Disable the expression signature calculation.
    DisableSignature( _self::LEM, $ )

    # Return the veiling label.
    GetVeilingLabel( _self::LEM, $ ) :: symbol

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

    # Unveil the expression <x> with the veiled variables.
    Unveil( _self::LEM, x::anything, $ ) :: anything

    # Unveil the expression <x> with signature values.
    UnveilSig( _self::LEM, x::anything, $ ) :: anything

    # Unveil the expression <x> with veiling labels permutation map.
    UnveilImap( _self::LEM, x::anything, $ ) :: anything

    # Unveil the expression <x> with signature values permutation map.
    SigImap( _self::LEM, x::anything, $ ) :: anything

    # Return a list of the veiling labels.
    VeilList( _self::LEM,
              { reverse::boolean := false }, $ ) :: list(anything)

    # Return a list of the signature labels.
    SigList( _self::LEM,
             { reverse::boolean := false }, $ ) :: list(anything)

    # Return the unordered veiling list.
    VeilUnorderedList( _self::LEM, $ ) :: list

    # Return the unordered signature list.
    SigUnorderedList( _self::LEM, $ ) :: list

    # Return the size of the internal veiling table.
    VeilTableSize( _self::LEM, $ ) :: nonnegint

    # Return the size of the internal signature table.
    SigTableSize( _self::LEM, $ ) :: nonnegint

    # Return the veiling list and the permutation that sorts it.
    VeilTableImap( _self::LEM,
                   { reverse::boolean := false }, $ )
                 :: {list(anything = anything), []}

    # Return the signature list and the permutation that sorts it.
    SigTableImap( _self::LEM,
                  { reverse::boolean := false }, $ )
                :: {list(anything = anything), []}

    # Append the veiled expression <x> to the veiling table.
    TablesAppend( _self::LEM, x::anything, $ ) :: indexed

    # Substitute the reversed signature values of the internal signature table
    # in the expression <x>.
    SubsSig( _self::LEM, x::anything, $ ) :: anything

    # Substitute the reversed veiling variables of the internal veiling table
    # in the expression <x>.
    SubsVeil( _self::LEM, x::anything, $ ) :: anything

    # Clear all the veiling variables of the internal veiling table.
    ForgetVeil( _self::LEM, $ )

    # Compute the signature of the expression <x> modulo <p> (the default is
    # the internal signature value) also by using the internal signature table
    # of the veiled expressions.
    Signature( _self::LEM, x::anything, p::prime := _self:-m_SigValue, $ )
             :: nonnegint
```

## Usage

In case you have no time to read the description and realize how it should or should not work, refer to the files `Test_Maple2021minus.mw` and `Test_Maple2021plus.mw` in the `tests` folder.

### 🚧 Attention! 🚧

Maple object-oriented programming features have slightly changed in 2021, which online documentation states:

> As of Maple 2021, if the method has a formal parameter named `_self`, references to its object's local or exported variables may be written without prefixing them. That is, `_self:-variable` may be written as just `variable`. Maple will add the `self:-` prefix internally when the method is simplified.

> As of Maple 2021, a message-passing form of method call can be used, which will automatically pass the object as an argument if the method has a formal parameter named `_self`.

> Another way to invoke a method, similar to that used in other object-oriented languages, was introduced in Maple 2021. In this form, the object name is qualified by the method name, `object_name:-method_name( argument )` just as it can be in the function mechanism described above. However, the *object can be omitted from the argument sequence*.

For further information please refer to this [link](https://fr.maplesoft.com/support/help/Maple/view.aspx?path=object/methods).

# SIG (Expression Signature Module)

The module `SIG` (Signature) module in also contained in the `LEM` repository. It contains the functions and workarounds to calculate the signature of "almost" any kind of expression. The signature is a hashing value that can be used to calculate the similarity between two expressions, or to check how likely an expression is algebraically non-null (see Schwarz-Zippel Lemma).

## Module description

If you want a full description of the `SIG` module type:

```
> Describe(SIG);
```

This command will generate a brief description of the module and all the procedures and other objects present in the `SIG.mpl` file, which will be (very) similar to the following code.

```
# Expression Signature module.
# Compute the signature of an expression <expr> modulo a prime number <p>.
module SIG( expr::algebraic, p::prime, max_iter::posint := 10, $ ) :: algebraic

    # Print 'SIG' module information.
    Info( )

    # 'SIG' module load procedure.
    ModuleLoad( )

    # 'SIG' module unload procedure.
    ModuleUnload( )

    # Extract the arguments of a function <func> from an expression <expr>.
    ExtractArgs( expr::algebraic, func::name, $ ) :: list(algebraic)

    # Transform the absolute value function in an expression <expr>.
    AbsTransform( expr::algebraic, $ ) :: algebraic

    # Transform the logarithm function in an expression <expr>.
    LogTransform( expr::algebraic, $ ) :: algebraic

    # Transform the trigonometric functions in an expression <expr>.
    WeierstrassTransform( expr::algebraic, $ ) :: algebraic

    # Transform an expression <expr> using a maximum number of substitutions
    # recursions <max_iter>.
    Transform( expr::algebraic, max_iter::posint := 10, $ ) :: algebraic

    # Compute the signature of an expression <expr> modulo a prime number <p>.
    Signature( expr::algebraic, p::prime, $ ) :: nonnegint
```

## Authors

### Current version authors:

- *Davide Stocco*,
  Department of Industrial Engineering,
  University of Trento \
  email: davide.stocco@unitn.it

- *Matteo Larcher*,
  Department of Industrial Engineering,
  University of Trento \
  email: matteo.larcher@unitn.it

- *Enrico Bertolazzi*,
  Department of Industrial Engineering,
  University of Trento \
  email: enrico.bertolazzi@unitn.it

### Inspired by the work of:

- *Wenqin Zhou* (former affiliation),
  Department of Applied Mathematics,
  University of Western Ontario

- *David J. Jeffrey*,
  Department of Applied Mathematics,
  University of Western Ontario

- *Jacques Carette*,
  Department of Computing and Software,
  McMaster University

- *Robert M. Corless*,
  Department of Applied Mathematics,
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
@inproceedings{zhou2006hierarchical,
  title = {Hierarchical representations with signatures for large expression management},
  author = {Zhou, Wenqin and Carette, Jacques and Jeffrey, David J. and Monagan, Michael B.},
  booktitle = {Proceedings AISC 2006, LNCS 4120}
  editor = {Calmet, J. and Ida, T. and Wang, D.},
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

```
@article{zhou2008fraction,
  title = {Fraction-free matrix factors: new forms for LU and QR factors},
  author = {Zhou, Wenqin and Jeffrey, David J.},
  journal = {Frontiers of Computer Science in China},
  volume = {2},
  pages = {67--80},
  year = {2008},
  publisher = {Springer}
}
```

*Documented by Davide Stocco*