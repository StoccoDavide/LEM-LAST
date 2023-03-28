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

unprotect('LULEM');
LULEM := module()

  description "LU decomposition with LEM (Large Expressions Management).";

  option object;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local m_Verbose        := false;
  local m_TimeLimit      := 1;
  local m_DegreeCost_fun := NULL;
  local m_LEM            := Object(LEM);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Info::static := proc()

    description "Print 'LEM' module information.";

    printf(
      "+--------------------------------------------------------------------------+\n"
      "| 'LULEM' module version 1.0 - BSD 3-Clause License - Copyright (c) 2023   |\n"
      "| Current version authors:                                                 |\n"
      "|   D. Stocco, M. Larcher, E. Bertolazzi.                                  |\n"
      "| Inspired by the work of:                                                 |\n"
      "|   W. Zhou, D. J. Jeffrey, J. Carette and R. M. Corless.                  |\n"
      "+--------------------------------------------------------------------------+\n"
    );
    return NULL;
  end proc: # Info

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleLoad::static := proc()

    description "'LULEM' module load procedure.";

    local i, lib_base_path;

    lib_base_path := NULL;
    for i in [libname] do
      if (StringTools:-Search("LULEM", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if (lib_base_path = NULL) then
      error "Cannot find 'LULEM' module";
    end if;
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleUnload::static := proc()
    description "'LULEM' module unload procedure.";
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Code inspired by:
  # https://www.mapleprimes.com/questions/235996-Is-There-Any-Command-Or-Function-For

  export GetDegrees::static := proc(
    _self::LULEM,
    A::Matrix,
    $)::Matrix(nonnegint), Matrix(nonnegint);

    description "Get the degree matrices of the matrix <A>.";

    local i, j, k, m, n, r, c, ro, co;

    m, n := LinearAlgebra:-Dimensions(A);
    ro := Vector[column](m, k -> 1);
    co := Vector[row](n, k -> 1);
    r  := Vector[column](
      [seq(rtable_scanblock(A, [i,..], ':-NonZeros'), i = 1..m)]
    );
    c  := Vector[row](
      [seq(rtable_scanblock(A, [..,j], ':-NonZeros'), j = 1..n)]
    );
    return r.co, ro.c;
  end proc: # GetDegrees

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetLEM::static := proc(
    _self::LULEM,
    obj::LEM,
    $)

    description "Set the 'LEM' object <obj>.";

    _self:-m_LEM := obj;
    return NULL;
  end proc: # SetLEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetLEM::static := proc(
    _self::LULEM,
    $)::LEM;

    description "Get the 'LEM' object.";

    return _self:-m_LEM;
  end proc: # GetLEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Spy::static := proc(
    _self::LULEM,
    A::Matrix,
    $)::anything;

    description "Plot of non-zero values of the matrix <A>.";

    return plots:-sparsematrixplot(A, 'matrixview');
  end proc: # Spy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SpyLU::static := proc(
    _self::LULEM,
    A::Matrix,
    L::Matrix,
    U::Matrix,
    $)::anything;

    description "Plot of non-zero values of the matrices <A>, <L> and <U> "
      "with fill-in values.";

    local mat_0, mat_1, mat_2;
    mat_0 := map(x -> `if`(x = 0, 0, 1), A);
    mat_1 := map(x -> `if`(x = 0, 0, 1), L + U);
    mat_2 := map(x -> `if`(x = 1, 1, 0), mat_0 + mat_1);

    return [plots:-sparsematrixplot(mat_0, 'matrixview', color = "Black"),
            plots:-sparsematrixplot(mat_2, 'matrixview', color = "Red")];
  end proc: # SpyLU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export PermutationMatrices::static := proc(
    _self::LULEM,
    r::Vector(nonnegint),
    c::Vector(nonnegint),
    $)::Matrix(nonnegint), Matrix(nonnegint);

    description "Compute the LU decomposition premutation matrices provided "
      "the rows pivot vector <r> and the columns pivot vector <c>.";

    local m, n, i, P, Q;

    m := LinearAlgebra:-RowDimension(r);
    n := LinearAlgebra:-RowDimension(c);
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

  export SetVerbosity::static := proc(
    _self::LULEM,
    x::boolean,
    $)

    description "Set the verbosity of the package to <x>.";

    _self:-m_Verbose := x;
    return NULL;
  end proc: # SetVerbosity

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableVerbosity::static := proc(
    _self::LULEM,
    $)

    description "Enable the verbosity of the package.";

    _self:-m_Verbose := true;
    return NULL;
  end proc: # EnableVerbosity

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableVerbosity::static := proc(
    _self::LULEM,
    $)

    description "Disable the verbosity of the package.";

    _self:-m_Verbose := false;
    return NULL;
  end proc: # DisableVerbosity

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetTimeLimit::static := proc(
    _self::LULEM,
    x::numeric,
    $)

    description "Set the time limit of the package to <x>.";

    if (x < 0) then
      error "time limit must be a non-negative number.";
    end if;

    _self:-m_TimeLimit := x;
    return NULL;
  end proc: # SetTimeLimit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SolveLinearSystem::static := proc(
    _self::LULEM,
    T::table,
    b::Vector,
    $)

    description "Solve the factorized linear system <T> * x = b using the "
      "method specified in <T>['method'] field.";

    if (T["method"] = "LU") then
      return _self:-LUsolve(_self, T, b);
    elif (T["method"] = "FFLU") then
      return _self:-FFLUsolve(_self, T, b);
    elif (T["method"] = "QR") then
      return _self:-QRsolve(_self, T, b);
    elif (T["method"] = "QR2") then
      return _self:-QR2solve(_self, T, b);
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
