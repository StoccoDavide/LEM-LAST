# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           _        _    ____ _____                          #
#                          | |      / \  / ___|_   _|                         #
#                          | |     / _ \ \___ \ | |                           #
#                          | |___ / ___ \ ___) || |                           #
#                          |_____/_/   \_\____/ |_|                           #
#                       Linear Algebra Symbolic Toolbox                       #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco     (University of Trento)
#   Matteo Larcher    (University of Trento)
#   Enrico Bertolazzi (University of Trento)
#
# Inspired by the work of:
#   Wenqin Zhou       (University of Western Ontario) - Former affiliation
#   David J. Jeffrey  (University of Western Ontario)
#   Jacques Carette   (McMaster University)
#   Robert M. Corless (University of Western Ontario)
#
# License: BSD 3-Clause License
#
# This is a module for the 'LAST' (Linear Algebra Symbolic Toolbox) module. It
# contains the functions to solve linear systems of equations with large
# symbolic expressions. The module uses symbolic full pivoting LU and QR
# decompositions to solve linear systems. The 'LEM' (Large Expressions
# Management) module is used to avoid expression swell.
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

unprotect('LAST');
module LAST()

  description "Linear Algebra Symbolic Toolbox module.";

  option object;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local m_LEM               := NULL;
  local m_VerboseMode       := false;
  local m_WarningMode       := true;
  local m_TimeLimit         := 0.1;
  local m_MinDegreeStrategy := "product_1rc";
  local m_Results           := NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Info::static := proc()

    description "Print 'LAST' module information.";

    printf(
      "+--------------------------------------------------------------------------+\n"
      "| 'LAST' module version 1.0 - BSD 3-Clause License - Copyright (c) 2023    |\n"
      "| Current version authors:                                                 |\n"
      "|   D. Stocco, M. Larcher and E. Bertolazzi.                               |\n"
      "| Inspired by the work of:                                                 |\n"
      "|   W. Zhou, D. J. Jeffrey, J. Carette and R. M. Corless.                  |\n"
      "+--------------------------------------------------------------------------+\n"
    );
    return NULL;
  end proc: # Info

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleLoad::static := proc()

    description "'LAST' module load procedure.";

    local i, lib_base_path;

    lib_base_path := NULL;
    for i in [libname] do
      if (StringTools:-Search("LAST", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if (lib_base_path = NULL) then
      error("cannot find 'LAST' module.");
    end if;
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleUnload::static := proc()
    description "'LAST' module unload procedure.";
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleCopy::static := proc(
    _self::LAST,
    proto::LAST,
    $)

    description "Copy the objects <proto> into <self>.";

    _self:-m_LEM               := proto:-m_LEM;
    _self:-m_VerboseMode       := proto:-m_VerboseMode;
    _self:-m_WarningMode       := proto:-m_WarningMode;
    _self:-m_TimeLimit         := proto:-m_TimeLimit;
    _self:-m_MinDegreeStrategy := proto:-m_MinDegreeStrategy;
    _self:-m_Results           := proto:-m_Results;
  end proc: # ModuleCopy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export CheckInit::static := proc(
    _self::LAST,
    $)

    description "Check if the 'LAST' object is initialized.";

    if not type(_self:-m_LEM, LEM) then
      error("the 'LEM' object is not initialized, use 'LAST:-InitLEM(...)'' "
        "or other appropriate initialization methods first.");
    end if;
    return NULL;
  end proc: # CheckInit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export InitLEM::static := proc(
    _self::LAST,
    label::{symbol, string} := NULL,
    $)

    description "Initialize the 'LEM' object with veiling label <label>.";

    _self:-m_LEM := Object(LEM);
    if type(label, symbol) then
      _self:-m_LEM:-SetVeilingLabel(_self:-m_LEM, label);
    elif type(label, string) then
      _self:-m_LEM:-SetVeilingLabel(_self:-m_LEM, parse(label));
    end if;
    return NULL;
  end proc: # InitLEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ClearLEM::static := proc(
    _self::LAST,
    $)

    description "Clear the 'LEM' object.";

    _self:-m_LEM := NULL;
  end proc: # ClearLEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetLEM::static := proc(
    _self::LAST,
    obj::LEM,
    $)

    description "Set the 'LEM' object <obj>.";

    _self:-m_LEM := obj;
    return NULL;
  end proc: # SetLEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetLEM::static := proc(
    _self::LAST,
    $)::LEM;

    description "Get the 'LEM' object.";

    return _self:-m_LEM;
  end proc: # GetLEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableVerboseMode::static := proc(
    _self::LAST,
    $)

    description "Enable the verbosity of the module.";

    _self:-m_VerboseMode := true;
    _self:-m_LEM:-EnableVerboseMode(_self:-m_LEM);
    return NULL;
  end proc: # EnableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableVerboseMode::static := proc(
    _self::LAST,
    $)

    description "Disable the verbosity of the module.";

    _self:-m_VerboseMode := false;
    _self:-m_LEM:-DisableVerboseMode(_self:-m_LEM);
    return NULL;
  end proc: # DisableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetVerboseMode::static := proc(
    _self::LAST,
    mode::boolean,
    $)

    description "Set the verbosity of the module to <mode>.";

    _self:-m_VerboseMode := mode;
    _self:-m_LEM:-SetVerboseMode(_self:-m_LEM, mode);
    return NULL;
  end proc: # SetVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableWarningMode::static := proc(
    _self::LAST,
    $)

    description "Enable the warning mode of the module.";

    _self:-m_WarningMode := true;
    _self:-m_LEM:-EnableWarningMode(_self:-m_LEM);
    return NULL;
  end proc: # EnableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableWarningMode::static := proc(
    _self::LAST,
    $)

    description "Disable the warning mode of the module.";

    _self:-m_WarningMode := false;
    _self:-m_LEM:-DisableWarningMode(_self:-m_LEM);
    return NULL;
  end proc: # DisableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetWarningMode::static := proc(
    _self::LAST,
    mode::boolean,
    $)

    description "Set the warning mode of the module to <mode>.";

    _self:-m_WarningMode := mode;
    _self:-m_LEM:-SetWarningMode(_self:-m_LEM, mode);
    return NULL;
  end proc: # SetWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetTimeLimit::static := proc(
    _self::LAST,
    x::numeric,
    $)

    description "Set the time limit of the module to <x>.";

    if (x < 0) then
      error("time limit must be a non-negative number.");
    end if;

    _self:-m_TimeLimit := x;
    return NULL;
  end proc: # SetTimeLimit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetTimeLimit::static := proc(
    _self::LAST,
    $)::numeric;

    description "Get the time limit of the module.";

    return _self:-m_TimeLimit;
  end proc: # GetTimeLimit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export CheckResults::static := proc(
    _self::LAST,
    $)

    description "Check if the results of the last factorization are available.";

    if (numelems(_self:-m_Results) = 0) then
      error("the results of the last factorization are not available, use "
        "appropriate factorization and solution methods first.");
    end if;
    return NULL;
  end proc: # CheckResults

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetResults::static := proc(
    _self::LAST,
    field::string := "all",
    $)

    description "Get the results of the last factorization. If <field> is "
      "specified, only the field Results['field'] is returned.";

    # Check if the results are available
    _self:-CheckResults(_self);

    # Retrieve the results
    if (field = "all") then
      return _self:-m_Results;
    else
      return _self:-m_Results[field];
    end if;
  end proc: # GetResults

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ClearResults::static := proc( $ )

    description "Clear the results of the last factorization.";

    _self:-m_Results := table([]);
    return NULL;
  end proc: # ClearResults

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SolveLinearSystem::static := proc(
    _self::LAST,
    b::Vector,
    $)

    description "Solve the factorized linear system.";

    # Check if the LEM object is initialized
    _self:-CheckInit(_self);

    # Check if the results are available
    _self:-CheckResults(_self);

    # Solve the linear system
    if (_self:-m_Results["method"] = "LU") then
      return _self:-LUsolve(_self, b);
    elif (_self:-m_Results["method"] = "FFLU") then
      return _self:-FFLUsolve(_self, b);
    elif (_self:-m_Results["method"] = "QR") then
      return _self:-QRsolve(_self, b);
    elif (_self:-m_Results["method"] = "GJ") then
      return _self:-GJsolve(_self, b);
    else
      error("wrong or not available decomposition, use 'LAST:-LU()' or "
        "'LAST:-FFLU()' or 'LAST:-QR()' or 'LAST:-GJ()' first.");
    end if;

  end proc: # SolveLinearSystem

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ApplyLP::static := proc(
    _self::LAST,
    b::Vector,
    $)::Vector;

    description "Apply L^(-1)*P to the vector <b>.";

    # Check if the LEM object is initialized
    _self:-CheckInit(_self);

    # Check if the results are available
    _self:-CheckResults(_self);

    # Retrieve the results
    if (_self:-m_Results["method"] = "LU") then
      return _self:-LUapplyLP(_self, b);
    elif (_self:-m_Results["method"] = "FFLU") then
      return _self:-FFLUapplyLP(_self, b);
    else
      error("wrong or not available decomposition, use 'LAST:-LU()' or "
        "'LAST:-FFLU()' first.");
    end if;
  end proc: # GetUQT

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetUQT::static := proc(
    _self::LAST,
    $)::Matrix;

    description "Return the matrix U^T*Q.";

    # Check if the LEM object is initialized
    _self:-CheckInit(_self);

    # Check if the results are available
    _self:-CheckResults(_self);

    # Retrieve the results
    if (_self:-m_Results["method"] = "LU") then
      return _self:-LUgetUQT(_self);
    elif (_self:-m_Results["method"] = "FFLU") then
      return _self:-FFLUgetUQT(_self);
    else
      error("wrong or not available decomposition, use 'LAST:-LU()' or "
        "'LAST:-FFLU()' first.");
    end if;
  end proc: # GetUQT

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Spy::static := proc(
    _self::LAST,
    A::Matrix,
    $)::anything;

    description "Plot of non-zero values of the matrix <A>.";

    return plots:-sparsematrixplot(A, 'matrixview');
  end proc: # Spy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local GCD::static := proc(
    _self::LAST,
    expr::{list, Vector, Matrix},
    $)::algebraic;

    description "Compute the greatest common divisor of the elements of the "
      "expression <expr>.";

    local expr_tmp, num, den, i;

    if not type(expr, list) then
      expr_tmp := convert(expr, list);
    else
      expr_tmp := expr;
    end if;

    if (nops(expr_tmp) > 0) then
      num := numer(expr_tmp[1]);
      den := denom(expr_tmp[1]);
      if (nops(expr_tmp) > 1) then
        for i in expr_tmp[2..-1] do
          num := gcd(num, numer(i));
          den := gcd(den, denom(i));
        end do;
      end if;
      return num/den;
    else
      return 1;
    end if;
  end proc: # GCD

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SpyLU::static := proc(
    _self::LAST,
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
    _self::LAST,
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

$include "./lib/LAST_Pivoting.mpl"
$include "./lib/LAST_LU.mpl"
$include "./lib/LAST_FFLU.mpl"
$include "./lib/LAST_GJ.mpl"
$include "./lib/LAST_QR.mpl"

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # LAST

# That's all folks!
