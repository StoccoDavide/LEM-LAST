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
# pivoting LU decomposition to solve linear systems. The `LEM` (Large Expressions
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
          SolveLinearSystem,
          VeilingStrategy,
          SetVeilingStrategyCost,
          SetPivotStrategy,
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
          VeilingStrategy_par,
          PivotStrategy_type,
          PivotingStrategy_Row,
          PivotingStrategy_Col,
          PivotingStrategy_Sum,
          PivotingStrategy_Prod,
          PivotingStrategy_Min,
          PivotingStrategy_Val;

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
      "Original code: W. Zhou, D. J. Jeffrey, J. Carette and R. M. Corless.\n"
    );

    lib_base_path := null;
    for i in [libname] do
      if (StringTools[Search]("LULEM", i) <> 0) then
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

    LULEM:-Verbose             := false;
    LULEM:-TimeLimit           := 1;
    LULEM:-VeilingStrategy_par := 15;
    LULEM:-PivotStrategy_type  := LULEM:-PivotingStrategy_Sum;
    return NULL;
  end proc: # InitLULEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  GetDegrees :=  proc(
    A::{Matrix},
    $)::{Matrix(nonnegint)};

    description "Get the degree of the matrix <A>.";

    local i, j, k, m, n, r, c, ro, co;

    m, n := LinearAlgebra[Dimensions](A);
    r  := Vector[column](m);
    c  := Vector[row](n);
    ro := Vector[column](m, k -> 1);
    co := Vector[row](n, k -> 1);
    for i from 1 to m do
      for j from 1 to n do
        if (A[i,j] <> 0) then
          r[i] := r[i]+1;
          c[j] := c[j]+1;
        end if;
      end do;
    end do;
    return r.co, ro.c;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Spy := proc(
    A::{Matrix},
    $)::{anything};

    description "Plot of non-zero values of the matrix <A>.";

    return plots:-sparsematrixplot(A, 'matrixview');
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PermutationMatrices := proc(
    r::{Vector(nonnegint)},
    c::{Vector(nonnegint)},
    $)

    description "Compute the LU decomposition premutation matrices provided "
      "the rows pivot vector <r> and the columns pivot vector <c>.";

    local m, n, i, P, Q;

    m := LinearAlgebra[RowDimension](r):
    n := LinearAlgebra[RowDimension](c):
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
    $)::{nothing};

    description "Set the verbosity of the package to <x>.";

    LULEM:-Verbose := x;
    return NULL;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  EnableVerbosity := proc(
    $)::{nothing};

    description "Enable the verbosity of the package.";

    LULEM:-Verbose := true;
    return NULL;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  DisableVerbosity := proc(
    $)::{nothing};

    description "Disable the verbosity of the package.";

    LULEM:-Verbose := false;
    return NULL;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SetTimeLimit := proc(
    x::{numeric},
    $)::{nothing};

    description "Set the time limit of the package to <x>.";

    if (x < 0) then
      error "LULEM::SetTimeLimit(...): time limit must be a non-negative number.";
    end if;

    LULEM:-TimeLimit := x;
    return NULL;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Cost := proc(
    x::{anything},
    $)::{integer};

    description "Compute the cost of the expression <x>.";

    local tmp;

    if not( type(x,'algebraic') or type(x,'list') ) then
      tmp := convert(x,'list');
    else
      tmp := x;
    end if;

    return subs(
      subscripts      = 0,
      assignments     = 0,
      additions       = 1,
      multiplications = 2,
      divisions       = 3,
      functions       = 2,
      codegen[cost](tmp)
    );
  end proc: # PivotCost

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SolveLinearSystem := proc(
    T::{table},
    b::{Vector},
    V::{symbol, function},
    $)

    description "Solve the factorized linear system <T> * x = b using the "
      "method specified in <T>['method'] field. The input <V> is the veiling "
      "strategy to be used.";

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

  #  __     __   _ _ _
  #  \ \   / /__(_) (_)_ __   __ _
  #   \ \ / / _ \ | | | '_ \ / _` |
  #    \ V /  __/ | | | | | | (_| |
  #     \_/ \___|_|_|_|_| |_|\__, |
  #                          |___/

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilingStrategy := proc(
    x::{algebraic},
    $)::{boolean};

    description "Comupte the veiling strategy value for the value <x>.";

    return evalb(LULEM:-Cost(x) > LULEM:-VeilingStrategy_par);
  end proc: # VeilingStrategy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SetVeilingStrategyCost := proc(
    c::{nonnegint},
    $)::{nothing};

    description "Set the veiling strategy parameter to <c>.";

    LULEM:-VeilingStrategy_par := c;
    return NULL;
  end proc: # SetVeilingStrategyCost

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$include "./lib/LULEM_Pivoting.mpl"
$include "./lib/LULEM_LU.mpl"
$include "./lib/LULEM_FFLU.mpl"
$include "./lib/LULEM_QR.mpl"
$include "./lib/LULEM_QR2.mpl"

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module:

# That's all folks!
