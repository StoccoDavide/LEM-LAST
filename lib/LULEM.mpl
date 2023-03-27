# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                        _    _   _ _     _____ __  __                        #
#                       | |  | | | | |   | ____|  \/  |                       #
#                       | |  | | | | |   |  _| | |\/| |                       #
#                       | |__| |_| | |___| |___| |  | |                       #
#                       |_____\___/|_____|_____|_|  |_|                       #
#          LU and QR decomposition with Large Expressions Management          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Authors of the current version:
#  Davide Stocco     (University of Trento)
#  Matteo Larcher    (University of Trento)
#  Enrico Bertolazzi (University of Trento)
#
# Authors of the original code:
#   Wenqin Zhou       (University of Western Ontario) - Former affiliation
#   David J. Jeffrey  (University of Western Ontario)
#   Jacques Carette   (McMaster University)
#   Robert M. Corless (University of Western Ontario)
#
# License: BSD 3-Clause License
#
# This is a module for the 'LULEM' (LU and QR decomposition Large Expressions
# Management) package. It contains the functions to solve linear systems of
# equations with large symbolic expressions. The module uses a symbolic full
# pivoting LU decomposition to solve linear systems. The 'LEM' (Large Expressions
# Management) package is used to avoid expression swell. Moreover, it also
# provides a full symbolic QR decomposition.
#
# The following code is hopefully an improved version of the original version
# provided inthe following PhD thesis:
#
#   Wenqin Zhou, Symbolic Computation Techniques for Solveing Large Expressions
#   Problems from Mathematics and Engineering (2007), Faculty of Graduate Studies,
#   The University of Western Ontario London, Ontario, Canada.
#
# We would like to thank Jacques Carette for providing the original code that we
# have used to develop this module.

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

unprotect(LULEM);
LULEM := module()

  export  SetVerbosity,
          EnableVerbosity,
          DisableVerbosity,
          SetTimeLimit,
          PermutationMatrices,
          GetDegrees,
          Spy,
          SpyLU,
          SolveLinearSystem,
          SetMinDegreeStrategy,
          LU,
          QR,
          QR2,
          FFLU,
          FF2LU;

  local   ModuleLoad,
          ModuleUnload,
          Cost,
          PivotCost,
          Pivoting,
          LUsolve,
          QRsolve,
          QR2solve,
          FFLUsolve,
          InitLULEM,
          Verbose,
          TimeLimit,
          DegreeCost_fun,
          DegreeCost_none,
          DegreeCost_row,
          DegreeCost_col,
          DegreeCost_sum,
          DegreeCost_prod,
          DegreeCost_prod2,
          DegreeCost_min,
          DegreeCost_max,
          PivotingCompare;

  option  package,
          load   = ModuleLoad,
          unload = ModuleUnload;

  description "LU decomposition with LEM (Large Expressions Management).";

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #   __  __           _       _
  #  |  \/  | ___   __| |_   _| | ___
  #  | |\/| |/ _ \ / _` | | | | |/ _ \
  #  | |  | | (_) | (_| | |_| | |  __/
  #  |_|  |_|\___/ \__,_|\__,_|_|\___|
  #

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ModuleLoad := proc()

    description "'LULEM' module load procedure.";

    local i, lib_base_path;

    printf(
      "'LULEM' module version 1.0 - BSD 3-Clause License - Copyright (c) 2023\n"
      "Current version: D. Stocco, M. Larcher, E. Bertolazzi.\n"
      "Inspired by the code of: W. Zhou, D. J. Jeffrey, J. Carette and R. M. Corless.\n"
    );

    lib_base_path := null;
    for i in [libname] do
      if (StringTools:-Search("LULEM", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if (lib_base_path = null) then
      error "Cannot find 'LULEM' module" ;
    end if;

    LULEM:-InitLULEM();

    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ModuleUnload := proc()
    description "'LULEM' module unload procedure.";
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  InitLULEM := proc()

    description "Initialize 'LULEM' module internal variables.";

    LULEM:-Verbose        := false;
    LULEM:-TimeLimit      := 1;
    LULEM:-DegreeCost_fun := LULEM:-DegreeCost_prod;
    return NULL;
  end proc: # InitLULEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Code inspired by:
  # https://www.mapleprimes.com/questions/235996-Is-There-Any-Command-Or-Function-For

  GetDegrees := proc(
    A::{Matrix},
    $)::{Matrix(nonnegint)}, {Matrix(nonnegint)};

    description "Get the degree matrices of the matrix <A>.";

    local i, j, k, m, n, r, c, ro, co;

    m, n := LinearAlgebra:-Dimensions(A);
    ro := Vector[column](m, k -> 1);
    co := Vector[row](n, k -> 1);
    r  := Vector[column](
      [seq(rtable_scanblock(A, [i,..], ':-NonZeros'), i = 1..m)]
    ):
    c  := Vector[row](
      [seq(rtable_scanblock(A, [..,j], ':-NonZeros'), j = 1..n)]
    ):

    # Old code
    #r  := Vector[column](m);
    #c  := Vector[row](n);
    #for i from 1 to m do
    #  for j from 1 to n do
    #    if (A[i,j] <> 0) then
    #      r[i] := r[i]+1;
    #      c[j] := c[j]+1;
    #    end if;
    #  end do;
    #end do;

    return r.co, ro.c;
  end proc: # GetDegrees

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Spy := proc(
    A::{Matrix},
    $)::{anything};

    description "Plot of non-zero values of the matrix <A>.";

    return plots:-sparsematrixplot(A, 'matrixview');
  end proc: # Spy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SpyLU := proc(
    A::{Matrix},
    L::{Matrix},
    U::{Matrix},
    $)::{anything};

    description "Plot of non-zero values of the matrices <A>, <L> and <U> "
      "with fill-in values.";

    local A0, A1, A2;
    A0 := map(x->`if`(x = 0, 0, 1), A);
    A1 := map(x->`if`(x = 0, 0, 1), L+U);
    A2 := map(x->`if`(x = 1, 1, 0), A0+A1);

    return [plots:-sparsematrixplot(A0, 'matrixview', color="Blue"),
            plots:-sparsematrixplot(A2, 'matrixview', color="Red")];
  end proc: # SpyLU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PermutationMatrices := proc(
    r::{Vector(nonnegint)},
    c::{Vector(nonnegint)},
    $)::{Matrix(nonnegint)}, {Matrix(nonnegint)};

    description "Compute the LU decomposition premutation matrices provided "
      "the rows pivot vector <r> and the columns pivot vector <c>.";

    local m, n, i, P, Q;

    m := LinearAlgebra:-RowDimension(r):
    n := LinearAlgebra:-RowDimension(c):
    P := Matrix(m, m);
    Q := Matrix(n, n);
    for i from 1 to m by 1 do
      P[i, r[i]] := 1;
    end do;
    for i from 1 to n by 1 do
      Q[c[i], i] := 1;
    end do;
    return P, Q;
  end proc: # Permutatiomnatrices

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SetVerbosity := proc(
    x::{boolean},
    $)

    description "Set the verbosity of the package to <x>.";

    LULEM:-Verbose := x;
    return NULL;
  end proc: # SetVerbosity

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  EnableVerbosity := proc( $ )

    description "Enable the verbosity of the package.";

    LULEM:-Verbose := true;
    return NULL;
  end proc: # EnableVerbosity

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  DisableVerbosity := proc( $ )

    description "Disable the verbosity of the package.";

    LULEM:-Verbose := false;
    return NULL;
  end proc: # DisableVerbosity

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SetTimeLimit := proc(
    x::{numeric},
    $)

    description "Set the time limit of the package to <x>.";

    if (x < 0) then
      error "LULEM::SetTimeLimit(...): time limit must be a non-negative number.";
    end if;

    LULEM:-TimeLimit := x;
    return NULL;
  end proc: # SetTimeLimit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SolveLinearSystem := proc(
    T::{table},
    b::{Vector},
    V::{symbol},
    $)

    description "Solve the factorized linear system <T> * x = b using the "
      "method specified in <T>['method'] field and veiling label <V>.";

    if (T["method"] = "LU") then
      return LULEM:-LUsolve(T, b, V);
    elif (T["method"] = "FFLU") then
      return LULEM:-FFLUsolve(T, b, V);
    elif (T["method"] = "QR") then
      return LULEM:-QRsolve(T, b, V);
    elif (T["method"] = "QR2") then
      return LULEM:-QR2solve(T, b, V);
    end
  end proc: # SolveLinearSystem

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$include "./lib/LULEM_Pivoting.mpl"
$include "./lib/LULEM_LU.mpl"
$include "./lib/LULEM_FFLU.mpl"
$include "./lib/LULEM_QR.mpl"
$include "./lib/LULEM_QR2.mpl"

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # LULEM

# That's all folks!
