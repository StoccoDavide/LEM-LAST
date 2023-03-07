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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

unprotect(LULEM);
LULEM := module()

  export  SetVerbosity,
          PermutationMatrices,
          SolveLinearSystem,
          VeilingStrategy,
          LU,
          QR,
          FFLU,
          FF2LU;

  local   ModuleLoad,
          ModuleUnload,
          PivotCost,
          DoPivoting,
          LUsolve,
          QRsolve,
          FFLUsolve,
          InitLULEM,
          Verbose,
          VeilingStrategy_par1,
          VeilingStrategy_par2;

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

    description "'LULEM' module load procedure";

    local i, lib_base_path;

    printf(
      "'LEM' module version 1.0, BSD 3-Clause License - Copyright (C) 2023\n"
      "D. Stocco, M. Larcher, E. Bertolazzi,\n"
      "W. Zhou, D. J. Jeffrey, J. Carette and R. M. Corless.\n"
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
    description "Module 'LULEM' module unload procedure";
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  InitLULEM := proc()
    description "Initialize 'LULEM' module internal variables";
    LULEM:-Verbose              := false;
    LULEM:-VeilingStrategy_par1 := 10;
    LULEM:-VeilingStrategy_par2 := 100;
    return NULL;
  end proc: # InitLULEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PermutationMatrices := proc( r::{Vector}, c::{Vector}, $ )

    description "Compute the LU decomposition premutation matrix provided the "
                "rows the pivot vector <r> and the columns the pivot vector <c>.";

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

  SetVerbosity := proc( x::{boolean}, $ )::{nothing};
    description "Set the verbosity of the package to <x>.";
    LULEM:-Verbose := x;
    return NULL;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SolveLinearSystem := proc( T::{table}, b::{Vector}, V::{symbol, function}, $)

    if T["method"] = "LU" then
      return LULEM:-LUsolve( T, b, V );
    elif T["method"] = "FFLU" then
      return LULEM:-FFLUsolve( T, b, V );
    elif T["method"] = "QR" then
      return LULEM:-QRsolve( T, b, V );
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

  VeilingStrategy := proc( x::{algebraic}, $ )::{boolean};
    description "Veiling strategy: ......";
    return evalb( (nops(indets(x))+1)*length(x) > LULEM:-VeilingStrategy_par1
                  or length(x) > LULEM:-VeilingStrategy_par2 );
  end proc: # VeilingStrategy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$include "./lib/LULEM_Pivoting.mpl"
$include "./lib/LULEM_LU.mpl"
$include "./lib/LULEM_FFLU.mpl"
$include "./lib/LULEM_QR.mpl"

end module:

# That's all folks!
