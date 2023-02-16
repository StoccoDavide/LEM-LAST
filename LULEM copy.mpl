# # # # # # # # # # # # # # # # # # # # # # # # # #
#          _    _   _ _     _____ __  __          #
#         | |  | | | | |   | ____|  \/  |         #
#         | |  | | | | |   |  _| | |\/| |         #
#         | |__| |_| | |___| |___| |  | |         #
#         |_____\___/|_____|_____|_|  |_|         #
# LU Decomposition for Large Expression Matrices  #
# # # # # # # # # # # # # # # # # # # # # # # # # #

# Authors: Davide Stocco and Matteo Larcher
# Date:    12/01/2023

# This is a module for the LULEM package. It contains the functions to solve
# systems of linear equations with large expressions. The module uses the LU
# decomposition of the matrix of the system. The module is hopefully a better
# version of the code provided in the following PhD thesis:
#
#   Wenqin Zhou, Symbolic Computation Techniques for Solveing Large Expressions
#   Problems from Mathematics and Engineering (2007), Faculty of Graduate Studies,
#   The University of Western Ontario London, Ontario, Canada.

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LULEM := module()

  # Exported variables
  export  SetModuleOptions,
          Veil,
          VeilDepth,
          Unveil,
          ShowVeil,
          SubsVeil,
          ForgetVeil,
          SquareLUwithpivoting,
          SolveSquareLUpivot,
          Solve,
          LEMStrategy_n,
          LEMStrategy_L,
          LEMStrategy_Ls,
          LEMStrategy_LB,
          PivotStrategy_Llength,
          PivotStrategy_Slength,
          PivotStrategy_Lindets,
          PivotStrategy_Sindets,
          PivotStrategy_numeric,
          ZeroStrategy_length,
          ZeroStrategy_normalizer,
          Info,
          License;

  # Local variables
  local   ModuleLoad,
          LastUsed,
          ModuleUnload,
          InitLULEM,
          Protect,
          NextLabel,
          Signature,
          UnVeilTable,
          lib_base_path;


  # Package options
  option  package,
          load   = ModuleLoad,
          unload = ModuleUnload;

  description "LU Decomposition for Large Expression Matrices";

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

    local i;

    # Display module init message
    printf(cat(
      "'LULEM' module version beta-0.0, ",
      "MIT License - Copyright (C) 2023, D. Stocco & M. Larcher, ",
      "University of Trento, Italy"
      ));

    # Library path
    lib_base_path := null;
    for i in [libname] do
      if (StringTools[Search]("LULEM", i) <> 0) then
        lib_base_path := i;
      end;
    end;
    if (lib_base_path = null) then
      error "cannot find 'LULEM' module" ;
    end:

    # Initialize the module variables
    InitLULEM();

    # Protect module keywords
    Protect();

    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ModuleUnload := proc()

    description "Module 'LULEM' module unload procedure";

    printf("Unloading 'LULEM'\n");
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  InitLULEM := proc()

    description "Initialize 'LULEM' module internal variables";

    # Define module variables
    UnVeilTable := table();
    LastUsed    := table();
  end proc: # InitLULEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Protect := proc()

    # Protect global variables
    #protect(
    #);

  end proc: # Protect

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #   _                    _
  #  | |    ___   ___ __ _| |
  #  | |   / _ \ / __/ _` | |
  #  | |__| (_) | (_| (_| | |
  #  |_____\___/ \___\__,_|_|
  #

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Signature := proc(
    f,           # TODO
    p::{posint}, # TODO
    A,           # TODO
    $)

    description "???"; # TODO

    local t, sig;
    option remember;

    if (f::'rational') then
      f mod p;
    elif (f::'symbol') then
      Signature(f, p, A) := rand() mod p;
    elif (f::'indexed') then
      if type(f, A[anything]) then
        Signature(Unveil(f, 1), p, A)
      else
        Signature(f, p, A) := rand() mod p;
      end if;
    elif (f::`+`) then
      add(Signature(t, p, A) mod p, t = f) mod p;
    elif (f::`*`) then
      sig := 1;
      for t in f do
        sig := sig * Signature(t, p, A) mod p;
      end do;
      sig;
    elif f::(anything^rational) then
      Signature(op(1, f), p, A) &^ op(2, f) mod p;
    elif f::(anything^polynom) then
      sig := numtheory['phi'](p);
      t := Signature(op(2, f), sig, A);
      Signature(op(1, f), p, A) &^ t mod p;
    else
      ERROR(
        "LULEM::Signature(...): expressions involving %1 not done.", f
        );
    end if;
  end: # Signature

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  NextLabel := proc(
    label::{symbol},
    vars::{set},
    $)

    description "Calculate the next veiling label for a given label <lapel> and "
      "a optional set of variables <vars>.";

    LastUsed[label] := LastUsed[label] + 1;
    label[LastUsed[label], `if`(nargs = 2, vars, NULL)]; # ???
  end proc: # NextLabel

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #   _____                       _
  #  | ____|_  ___ __   ___  _ __| |_
  #  |  _| \ \/ / '_ \ / _ \| '__| __|
  #  | |___ >  <| |_) | (_) | |  | |_
  #  |_____/_/\_\ .__/ \___/|_|   \__|
  #             |_|

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Veil := proc(
    x, # The expression to be veiled
    $)

    description "Veil an expression <x> and return a label to it.";

    local A, label;

    # label := op(procname);
    # A := NextLabel(label);
    # UnVeilTable[A] := x;
    # return A;

    label := op(procname);
    LastUsed[label] := LastUsed[label] + 1; # Inlining NextLabel for efficiency
    A := label[LastUsed[label]];
    UnVeil(A, 1) := x; # ???
    #UnVeilTable[A] := x; # ???
    return A;
  end proc: # Veil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilDepth := proc(
    x,           # The expression to be evaluated
    p::{posint}, # ???
    $)

    description "Return the veil depth of an expression <x> and return a label "
      "to it.";

    local A, s;

    s := Signature(x, p, op(procname));
    A := NextLabel(op(procname), op(3, s));
    UnVeilTable[A] := x;
    return A;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Unveil := proc(
    x,                                    # The expression to be Unveiled
    n::{nonnegint, identical (infinity)}, # The number of levels to Unveil
    $)

    description "Unveil the expression <x> up to <n> levels.";

    local out;
    option remember;

    if (nargs = 1) then
      return Unveil(x, 1);
    elif (nargs = 2) and (n = 0) then
      return x;
    elif x::atomic then
      return x;
    else
      return map(Unveil, x, n-1);
    end if;
  end proc: # Unveil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ShowVeil := proc(
    label::{symbol}, # The veiling label to be shown
    $)

    description "Show the veiling variables of the veiling label <label>.";

    local i;

    for i from 1 to LastUsed[label] do
      print(label[i] = UnVeilTable[label[i]]);
    end do;
  end proc: # ShowVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SubsVeil := proc(
    label::{symbol}, # The label of the veiling table to be shown
    x::{Vector},     # The expression to substitute
    $)

    description "Substitute the veiling variables of the veiling label <label> "
      "in the expression <x>.";

    local i, sub_vec;

    sub_vec := [seq(label[i] = UnVeilTable[label[i]], i = 1..LastUsed[label])]:
    return subs(op(ListTools[Reverse](sub_vec)), x);
  end proc: # SubsVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ForgetVeil := proc(
    label::{symbol}, # The veiling label to be forgotten
    $)

    description "Clear all the veiling variables of the veiling label <label>.";

    local i;

    subsop(4 = NULL, eval(Signature));
    LastUsed[label] := 0;
    UnVeilTable := table();
    return NULL;
  end proc: # ForgetVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SquareLUwithpivoting := proc(
    LU_NAG::{Matrix},              # Compact LU matrix (NAG-style form)
    Q::{symbol},                   # Symbol for the veiling
    Strategy_Veiling::{procedure}, # Veiling strategy
    Strategy_Pivots::{procedure},  # Pivoting strategy procedure
    Strategy_Zero::{procedure},    # Zero recognition strategy procedure
    $)

    description "Compute the LU decomposition of a square matrix provided the "
      "NAG-style LU matrix <LU_NAG> using the LEM strategy <Strategy_Veiling>, "
      "the veiling symbol <Q>, the pivoting strategy <Strategy_Pivots> and the "
      "zero recognition strategy <Strategy_Zero>.";

    local M, n, k, i, ii, j, m, L, U, mltip, l, r, kp, p, temp, normalize, row,
      flag, z;

    # Copy the input matrix
    M := copy(LU_NAG);
    n := LinearAlgebra[RowDimension](M):
    m := LinearAlgebra[ColumnDimension](M);

    # Check if the input matrix is square
    if (m <> n) then
      ERROR(
        "LULEM::SquareLUwithpivoting(...): the input matrix should be square."
        );
    end if;

    # L is lower matrix with l's as diagonal entries
    L := Matrix(LinearAlgebra[IdentityMatrix](n), shape = triangular[lower]);

    # U is upper matrix.
    U := Matrix(n, n, shape = triangular[upper]):

    # Create pivot vector
    r := Vector(n, k -> k);

    normalize := z -> `if`(Strategy_Veiling(z) > 0, Veil[Q](z), z);

    # Loop over columns and find the pivot element among the entries M[r[kp], kp],
    # M[r[kp+1], kp], ... M[r[n], kp] according to different pivoting strategies.
    # Interchange entries in pivot vector.

    row := 1;
    for kp from 1 to m while (row <= n) do
      p := row;
      flag := evalb(Strategy_Zero(M[r[row], kp]) = 0);
      for i from (row + 1) to n do
        # Once a pivot is found -- not "best"!
        if flag or Strategy_Pivots(M[r[i], kp], M[r[p], kp]) and
            not (Strategy_Zero(M[r[i], kp]) = 0) then
          p := i;
          break;
        end if;
      end do;

      # Only when the pivot is not equal to zero, the row elimination is performed
      # (the whole if statement). Else we will continue with the next column.

      if (Strategy_Zero(M[r[p], kp]) = 0) then

        WARNING(
          "LULEM::SquareLUwithpivoting(...): the matrix appears to be singular."
          );

      else

        if (p <> row) then
          (r[p], r[row]) := (r[row], r[p]);
        end if;

        userinfo(3, SquareLUwithpivoting, 'kp', kp, 'r', r);

        # Now do Gauss elimination steps to get new M and also keep L information
        # in new M. Packing everything into the M matrix during the computation of
        # the factors. The reason is that this code can be ported to restricted
        # memory environments, or equally, applied to very large matrices in a
        # large memory environment.

        for i from (row + 1) to n do
          mltip := normalize(M[r[i], kp]/M[r[row], kp]);
          M[r[i], kp] := mltip;
          for j from (kp + 1) to m do
            z := M[r[i], j] - mltip * M[r[row], j] ;
            M[r[i], j] := `if`(Strategy_Veiling(z) > 0, Veil[Q](z), z);
          end do;
        end do;

        userinfo(3, SquareLUwithpivoting, `M`, M);

      end if;
      row := row + 1;

    end do:

    userinfo(2, SquareLUwithpivoting, `r`, r, `M`, M);

    # Seperate new M into L and U
    for i from 1 to n do
      for j from 1 to m do
        if (i <= j) then
          U[i, j] := M[r[i], j];
        else
          L[i, j] := M[r[i], j];
        end if
      end do;
    end do;

    # Return the LU decomposition and the pivot vector
    return L, U, r;
  end proc: # SquareLUwithpivoting

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Solve the linear system Ax=b given the LU decomposition (PA=LU).
  # NOTE: The pivot vector r is returned from SquareLUwithpivoting function.
  SolveSquareLUpivot := proc(
    A::{Matrix},                   # Linear system matrix A
    r::{Vector},                   # Pivot vector
    b::{Vector},                   # Linear system vector b
    Q::{symbol},                   # Symbol for the veiling
    Strategy_Veiling::{procedure}, # Veiling strategy
    $)

    description "Solve the linear system Ax=b given the LU decomposition (PA=LU) "
      "privided the matrix <A>, the pivot vector <r>, the vector <b>, the veiling "
      "label <Q> and the veiling strategy <Strategy_Veiling> for the veiling.";

    local y, x, i, s, j, n, normalizer:

    userinfo(2, SolveSquareLUpivot, `b`, b, `r`, r);

    # Get number of rows in matrix A
    n := LinearAlgebra[RowDimension](A);

    # Create vector for solution of Ly=Pb
    y := Vector(n);

    # Create vector for solution of Ux=y
    x := Vector(n);

    normalizer := (y) -> `if`(Strategy_Veiling(y) > 0, Veil[Q](y), y);

    # Perform forward substitution to solve Ly=Pb
    userinfo(3, SolveSquareLUpivot,`n`, n, `y`, y);
    y[1] := normalizer(b[r[1]]);
    for i from 2 to n do
        y[i] := normalizer(b[r[i]]) - add(normalizer(A[r[i], j] * y[j]), j = 1..i-1);
    end do;

    # Perform backward substitution to solve Ux=y
    x[n] := normalizer(y[n]/A[r[n], n]);
    for i from n-1 to 1 by -1 do
      s := normalizer(y[i]) - add(normalizer(A[r[i],j] * x[j]), j = i+1..n);
      x[i] := normalizer(s/A[r[i], i]);
    end do;

    # Return solution vector x
    return x;
  end proc: # SolveSquareLUpivot

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Solve := proc(
    A::{Matrix},                   # Linear system matrix A
    b::{Vector},                   # Linear system vector b
    Q::{symbol},                   # Symbol for the veiling
    Strategy_Veiling::{procedure}, # Veiling strategy
    Strategy_Pivots::{procedure},  # Pivoting strategy procedure
    Strategy_Zero::{procedure},    # Zero recognition strategy procedure
    $)

    description "Solve the linear system Ax=b using LULEM algorithm, privided the "
      "matrix <A>, the vector <b>, the veiling symbol <Q>, the veiling strategy "
      "<Strategy_Veiling>, the pivoting strategy <Strategy_Pivots> and the zero "
      "recognition strategy <Strategy_Zero>.";

    local L, U, r, x, LU_NAG:

    # Get LU decomposition of A
    printf("LU decomposition...   ");
    L, U, r := SquareLUwithpivoting(
      A, Q, Strategy_Veiling, Strategy_Pivots, Strategy_Zero
      );
    printf("DONE\n");

    # Built the LU matrix (NAG-style)
    printf("LU_NAG matrix...      ");
    LU_NAG := L + U - Matrix(LinearAlgebra[RowDimension](L), shape = identity);
    printf("DONE\n");

    #print(LinearAlgebra[Diagonal](L), 'L' = L);
    #print(LinearAlgebra[Diagonal](U), 'U' = U);
    #print('r' = r);
    #print('LU_NAG' = LU_NAG);

    # Solve the linear system Ax=b given the LU decomposition (PA=LU).
    printf("SolveSquareLUpivot... ");
    x := SolveSquareLUpivot(LU_NAG, r, b, Q, Strategy_Veiling);
    printf("DONE\n");

    # Return solution vector x and the LU decomposition data
    return x, L, U, r;

  end proc: # Solve

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #  __     __   _ _ _
  #  \ \   / /__(_) (_)_ __   __ _
  #   \ \ / / _ \ | | | '_ \ / _` |
  #    \ V /  __/ | | | | | | (_| |
  #     \_/ \___|_|_|_|_| |_|\__, |
  #                          |___/

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  LEMStrategy_n := proc(
    x, # Expression to be analyzed
    $)::{integer};

    description "Veiling strategy: number of indeterminates in expression <x> "
      "minus 4.";

    return nops(indets(x)) - 4;
  end proc: # LEMStrategy_n

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  LEMStrategy_L:= proc(
    x, # Expression to be analyzed
    $)::{integer};

    description "Veiling strategy: length of expression <x> minus 50.";

    return length(x) - 50;
  end proc: # LEMStrategy_L

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  LEMStrategy_Ls:= proc(
    x, # Expression to be analyzed
    $)::{integer};

    description "Veiling strategy: length of expression <x> minus 120.";

    return length(x) - 120;
  end proc: # LEMStrategy_Ls

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  LEMStrategy_LB:= proc(
    x, # Expression to be analyzed
    $)::{integer};

    description "Veiling strategy: length of expression <x> minus 260.";

    return length(x) - 260;
  end proc: # LEMStrategy_LB

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #   ____  _            _   _
  #  |  _ \(_)_   _____ | |_(_)_ __   __ _
  #  | |_) | \ \ / / _ \| __| | '_ \ / _` |
  #  |  __/| |\ V / (_) | |_| | | | | (_| |
  #  |_|   |_| \_/ \___/ \__|_|_| |_|\__, |
  #                                  |___/

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotStrategy_Llength := proc(
    x::{algebraic}, # Expression to be analyzed
    y::{algebraic}, # Expression to be analyzed
    $)::{boolean};

    description "Pivoting strategy: choose the pivot with the largest length "
      "between expressions <x> and <y>.";

    return evalb(length(x) - length(y) > 0);
  end proc: # PivotStrategy_Llength

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotStrategy_Slength := proc(
    x::{algebraic}, # Expression to be analyzed
    y::{algebraic}, # Expression to be analyzed
    $)::{boolean};

    description "Pivoting strategy: choose the pivot with the smallest length "
      "between expressions <x> and <y>.";

    return evalb(length(x) - length(y) < 0);
  end proc: # PivotStrategy_Slength

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotStrategy_Lindets := proc (
    x::{algebraic}, # Expression to be analyzed
    y::{algebraic}, # Expression to be analyzed
    $)::{boolean};

    description "Pivoting strategy: choose the pivot with the largest number of "
      "indeterminates between expressions <x> and <y>.";

    return evalb((nops(indets(x)) - nops(indets(y))) > 0);
  end proc: # PivotStrategy_Lindets

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotStrategy_Sindets := proc (
    x::{algebraic}, # Expression to be analyzed
    y::{algebraic}, # Expression to be analyzed
    $)::{boolean};

    description "Pivoting strategy: choose the pivot with the smallest number of "
      "indeterminates between expressions <x> and <y>.";

    return evalb((nops(indets(x)) - nops(indets(y))) < 0);
  end proc: # PivotStrategy_Sindets

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotStrategy_numeric := proc (
    x::{algebraic}, # Expression to be analyzed
    y::{algebraic}, # Expression to be analyzed
    $)::{boolean};

    description "Pivoting strategy: choose the pivot with the largest numeric "
      "value between expressions <x> and <y>.";

    return evalb((abs(x) - abs(y)) > 0);
  end proc: # PivotStrategy_numeric

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #   _____
  #  |__  /___ _ __ ___
  #    / // _ \ '__/ _ \
  #   / /|  __/ | | (_) |
  #  /____\___|_|  \___/
  #

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ZeroStrategy_length := proc(
    x::{algebraic}, # Expression to be analyzed
    $)::{integer};

    description "Zero recognition strategy: length of expression <x>.";

    return length(x);
  end proc: # ZeroStrategy_length

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ZeroStrategy_normalizer := proc(
    x::{algebraic}, # Expression to be analyzed
    $) # FIXME: What is the return type of this procedure?

    description "Zero recognition strategy: normalizer of expression <x>.";

    return Normalizer(x)
  end proc: # ZeroStrategy_normalizer

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #   __  __ _
  #  |  \/  (_)___  ___
  #  | |\/| | / __|/ __|
  #  | |  | | \__ \ (__
  #  |_|  |_|_|___/\___|
  #

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Info := proc(
    $)::{nothing};

    description "Print the information about the package.";

    printf(cat(
      "This is a module for the LULEM package. It contains the functions to solve\n",
      "systems of linear equations with large expressions. The module uses the LU\n",
      "decomposition of the matrix of the system. The module is hopefully a better\n",
      "version of the code provided in the following PhD thesis:\n",
      "\n",
      "  Wenqin Zhou, Symbolic Computation Techniques for Solveing Large Expressions\n",
      "  Problems from Mathematics and Engineering (2007), Faculty of Graduate\n",
      "  Studies, The University of Western Ontario London, Ontario, Canada.\n"
    ));

    return NULL;
  end proc: # Info

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  License := proc(
    $)::{nothing};

    description "Print the license of the package.";

    printf(cat(
      "MIT License\n",
      "\n",
      "Copyright (c) 2023, Davide Stocco and Matteo Larcher\n",
      "\n",
      "Permission is hereby granted, free of charge, to any person obtaining a copy\n",
      "of this software and associated documentation files (the ""Software""), to deal\n",
      "in the Software without restriction, including without limitation the rights\n",
      "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n",
      "copies of the Software, and to permit persons to whom the Software is\n",
      "furnished to do so, subject to the following conditions:\n",
      "\n",
      "The above copyright notice and this permission notice shall be included in all\n",
      "copies or substantial portions of the Software.\n",
      "\n",
      "THE SOFTWARE IS PROVIDED ""AS IS"", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n",
      "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n",
      "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n",
      "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n",
      "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n",
      "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n",
      "SOFTWARE.\n"
    ));

    return NULL;
  end proc: # License

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module:

# That's all folks!