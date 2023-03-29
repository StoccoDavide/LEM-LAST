# LULEM (Linear Algebra Symbolic Toolbox)

This is a module for the `LAST` (Linear Algebra Symbolic Toolbox) module. It contains the functions to solve linear systems of equations with large symbolic expressions. The module uses symbolic full pivoting LU and QR decompositions to solve linear systems. The `LEM` (Large Expressions Management) module is used to avoid expression swell.

The code in this repository is hopefully an improved version of the code provided in Wenqin Zhou's PhD thesis *Symbolic Computation Techniques for Solving Large Expressions*.

We would like to thank *Jacques Carette* for providing the original code that we used to develop this module.

## Dependencies

The `LULEM` Maple module depends on the `LEM` (Large Expressions Management) module. This module contains functions to avoid expression swell and it is freely available at this [link](https://github.com/StoccoDavide/LEM).

## Installation

Firstly, install the dependency. The `LEM` module is freely available at this [link](https://github.com/StoccoDavide/LEM). For the installation of this dependency module, refer to the instructions present on the repository homepage.

To install the module you must have first installed Maple. Then open the `PackAndGo.mw` file and use the `!!!` button to *execute the entire worksheet*.

Then test the library in a Maple worksheet or document by typing:

```
> LAST:-Info(LAST); # Maple < 2021
> LAST:-Info();     # Maple >= 2021
```

Alternatively, you can use one of the tests file provided in the `tests` folder. If the module is loaded without errors, it is done!

## Module description

If you want a full description of the `LAST` module type:

```
> Describe(LAST);
```

This command will generate a brief description of the module and all the procedures and other objects present in the `LAST.mpl` file, which will be (very) similar to the following code.

```
# Linear Algebra Symbolic Toolbox module.
object LAST :: LAST:

    # Print 'LAST' module information.
    Info( )

    # 'LAST' module load procedure.
    ModuleLoad( )

    # 'LAST' module unload procedure.
    ModuleUnload( )

    # Copy the objects <proto> into <self>.
    ModuleCopy( _self::LAST, proto::LAST, $ )

    # Initialize the 'LEM' object with veiling label <label>.
    InitLEM( _self::LAST, label::{string, symbol} := NULL, $ )

    # Clear the 'LEM' object.
    ClearLEM( _self::LAST, $ )

    # Set the 'LEM' object <obj>.
    SetLEM( _self::LAST, obj::LEM, $ )

    # Get the 'LEM' object.
    GetLEM( _self::LAST, $ ) :: LEM

    # Enable the verbosity of the module.
    EnableVerboseMode( _self::LAST, $ )

    # Disable the verbosity of the module.
    DisableVerboseMode( _self::LAST, $ )

    # Set the time limit of the module to <x>.
    SetTimeLimit( _self::LAST, x::numeric, $ )

    # Get the time limit of the module.
    GetTimeLimit( _self::LAST, $ ) :: numeric

    # Clear the results of the last factorization.
    ClearResults( $ )

    # Get the results of the last factorization. If <field> is specified, only
    # the field Results['field'] is returned.
    GetResults( _self::LAST, field::string := "all", $ )

    # Solve the factorized linear system (LU)*x=b or (QR)*x=b.
    SolveLinearSystem( _self::LAST, b::Vector, $ )

    # Plot of non-zero values of the matrix <A>.
    Spy( _self::LAST, A::Matrix, $ ) :: anything

    # Plot of non-zero values of the matrices <A>, <L> and <U> with fill-in
    # values.
    SpyLU( _self::LAST, A::Matrix, L::Matrix, U::Matrix, $ ) :: anything

    # Compute the LU decomposition premutation matrices provided the rows pivot
    # vector <r> and the columns pivot vector <c>.
    PermutationMatrices( _self::LAST, r::Vector(nonnegint),
                         c::Vector(nonnegint), $ ) :: Matrix(nonnegint)

    # Compute the LU decomposition pivots vectors with minum degree provided
    # the step <k>, the temporary LU (NAG) matrix <M>, the rows permutation <r>
    # and the columns permutation <c>.
    Pivoting( _self::LAST, k::integer, M::Matrix, r::Vector(nonnegint),
              c::Vector(nonnegint), $ ) :: table

    # Compute the cost of the pivot <x>.
    PivotCost( _self::LAST, x::algebraic, $ ) :: integer

    # Get the degree matrices of the matrix <A>.
    GetDegrees( _self::LAST, A::Matrix, $ ) :: Matrix(nonnegint)

    # Set the strategy <str> for the minimum degree ordering.
    SetMinDegreeStrategy( _self::LAST, str::string := "product_1", $ )

    # Set the strategy <str> for the minimum degree ordering.
    DegreeCost( _self::LAST, val::table, $ ) :: integer

    # Compute the pivoting degree cost.
    DegreeCost_none( _self::LAST, val::table, $ ) :: integer

    # Compute the pivoting degree cost.
    DegreeCost_row( _self::LAST, val::table, $ ) :: integer

    # Compute the pivoting degree cost.
    DegreeCost_col( _self::LAST, val::table, $ ) :: integer

    # Compute the pivoting degree cost.
    DegreeCost_sum( _self::LAST, val::table, $ ) :: integer

    # Compute the pivoting degree cost.
    DegreeCost_prod( _self::LAST, val::table, $ ) :: integer

    # Compute the pivoting degree cost.
    DegreeCost_prod_1( _self::LAST, val::table, $ ) :: integer

    # Compute the pivoting degree cost.
    DegreeCost_min( _self::LAST, val::table, $ ) :: integer

    # Compute the pivoting degree cost.
    DegreeCost_max( _self::LAST, val::table, $ ) :: integer

    # Compute the pivoting strategy: given the current pivot <cur> and the next
    # pivot <val>, decide if to the next pivot is better than the current pivot
    # or not.
    PivotingCompare( _self::LAST, cur::table, val::table, $ ) :: boolean

    # Compute the LU decomposition of a square matrix <A>.
    LU( _self::LAST, A::Matrix, $ )

    # Solve the linear system Ax=b using LU decomposition provided the vector
    # <b>.
    LUsolve( _self::LAST, b::Vector, $ ) :: Vector

    # Compute the Givens QR decomposition of a square matrix <A>.
    QR( _self::LAST, A::Matrix, $ )

    # Solve the linear system Ax=b using QR decomposition provided the vector
    # <b>.
    QRsolve( _self::LAST, b::Vector, $ ) :: Vector
```

## Usage

In case you have no time to read the description and realize how it should or should not work, refer to the files `Test_00_LU.mw` and `Test_01_QR.mw` in the `tests` folder.

### 🚧 Attention! 🚧

Maple object-oriented programming features have slightly changed in 2021, which online documentation states:

> As of Maple 2021, if the method has a formal parameter named `_self`, references to its object's local or exported variables may be written without prefixing them. That is, `_self:-variable` may be written as just `variable`. Maple will add the `self:-` prefix internally when the method is simplified.

> As of Maple 2021, a message-passing form of method call can be used, which will automatically pass the object as an argument if the method has a formal parameter named `_self`.

> Another way to invoke a method, similar to that used in other object-oriented languages, was introduced in Maple 2021. In this form, the object name is qualified by the method name, `object_name:-method_name( argument )` just as it can be in the function mechanism described above. However, the *object can be omitted from the argument sequence*.

For further information please refer to this [link](https://fr.maplesoft.com/support/help/Maple/view.aspx?path=object/methods).

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
