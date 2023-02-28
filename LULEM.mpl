# # # # # # # # # # # # # # # # # # # # # # # # # # # #
#            _    _   _ _     _____ __  __            #
#           | |  | | | | |   | ____|  \/  |           #
#           | |  | | | | |   |  _| | |\/| |           #
#           | |__| |_| | |___| |___| |  | |           #
#           |_____\___/|_____|_____|_|  |_|           #
#  LU Decomposition with Large Expression Management  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Authors: Davide Stocco, Matteo Larcher, Wenqin Zhou, Jacques Carette and
#          Robert M. Corless
# Date:    12/01/2023

# This is a module for the LULEM package. It contains the functions to solve
# systems of linear equations with large expressions management. The module uses
# the pivoting LU decomposition to solve the system. The module is hopefully an
# improved version of the code provided in the following PhD thesis:
#
#   Wenqin Zhou, Symbolic Computation Techniques for Solveing Large Expressions
#   Problems from Mathematics and Engineering (2007), Faculty of Graduate Studies,
#   The University of Western Ontario London, Ontario, Canada.
#
# We would like to thank Jacques Carette for providing the original code that we
# have used to develop this module.
#
# The module is distributed under the BSD 3-clause License. See the LICENSE file
# for more information.

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

unprotect(LULEM);
LULEM := module()

  # Exported variables
  export  SetModuleOptions,
          Veil,
          UnVeil,
          ShowVeil,
          ListVeil,
          SubsVeil,
          ForgetVeil,
          AssignData,
          SubsData,
          UnAssignData,
          PermutationMatrix,
          LUD,
          FFLUD,
          Solve,
          VeilingStrategy_n,
          VeilingStrategy_L,
          VeilingStrategy_Ls,
          VeilingStrategy_LB,
          PivotStrategy_Llength,
          PivotStrategy_Slength,
          PivotStrategy_Lindets,
          PivotStrategy_Sindets,
          PivotStrategy_numeric,
          ZeroStrategy_length,
          ZeroStrategy_normalizer,
          LastUsed;

  # Local variables
  local   ModuleLoad,
          auxiliary,
          SolvePivotingLU,
          ModuleUnload,
          InitLULEM,
          UnVeilTable,
          StoredData,
          lib_base_path;


  # Package options
  option  package,
          load   = ModuleLoad,
          unload = ModuleUnload;

  description "LU Decomposition for Large Expression Matrices";

  LastUsed := table('sparse');

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
      "'LULEM' module version 1.1, ",
      "BSD 3-Clause License - Copyright (C) 2023, D. Stocco, M. Larcher, ",
      "E. Bertolazzi, W. Zhou, D. J. Jeffrey, J. Carette and R. M. Corless.\n"
      ));

    # Library path
    lib_base_path := null;
    for i in [libname] do
      if (StringTools[Search]("LULEM", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if (lib_base_path = null) then
      error "Cannot find 'LULEM' module" ;
    end if;

    # Initialize the module variables
    InitLULEM();

    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ModuleUnload := proc()

    description "Module 'LULEM' module unload procedure";

    unprotect(LastUsed);
    LastUsed    := NULL;
    UnVeilTable := NULL;
    UnAssignData();

    printf("Unloading 'LULEM'\n");
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  InitLULEM := proc()

    description "Initialize 'LULEM' module internal variables";

    # Define module variables
    UnVeilTable  := table('sparse' = table('sparse' = (0 = 0)));
    StoredData := [];
  end proc: # InitLULEM

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

    local i, s, c, label;

	  label := `if`(procname::indexed, op(procname), '_V');

	  if (label <> eval(label, 2)) then
	    error "LULEM::Veil(...): label %a is assigned a value already, please save"
        " its contents and unassign it.", label;
	  end if;

    # Recognize zero if we can, so that we don't hide zeros.
    c := Normalizer(x);

    # Remove the integer content and sign so that we don't hide them either.
    i := icontent(c);

    # And we really mean sign, here, and not signum, because the interesting
    # case is when c may be a polynomial.
    try
      s := sign(c); # sign is weak
      # Don't do anything if we can tell that the coefficient is just a number
      # or a simple multiple of a name (simple or indexed)
      if (s*i = c) or type(s*c/i,'name') then
        return c
      end if;
    catch:
      s := 1;
    end try;
    # Only if there is something complicated to hide we do actually hide it and
    # return a label.
    return s * i * auxiliary(s*c/i, label);
  end proc; # Veil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UnVeil := proc(
    x,                        # The expression to be Unveiled
    n::{nonnegint, infinity}, # The number of levels to UnVeil
    $)

    description "UnVeil the expression <x> up to <n> levels.";

    local a, b, level, label;

    label := `if`(procname::indexed, op(procname), '_V');
    level := `if`(nargs<2, 1, min(LastUsed[label]+1, n));
    a := x;
    # Always do at least 1 unveiling
    b := eval(a, op(eval(UnVeilTable[label]))[2]);
    from 2 to level while not Testzero(a - b) do
      a := b;
      b := eval(a, op(eval(UnVeilTable[label]))[2]);
    end do;
    return b;
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auxiliary := proc(
    x,     # The expression to be veiled
    label, # The veiling label to be used
    $)

    description "Auxiliary procedure to scope LastUsed etc and use option "
      "remember to detect duplicates. There is no nontrivial storage duplication "
      "because objects are hashed and stored uniquely. All this costs is a pointer.";

    option remember;

    unprotect(LastUsed);
    LastUsed[label] := LastUsed[label] + 1;
    protect(LastUsed);
    if LastUsed[label] = 1 then
      UnVeilTable[label] := table('sparse' = (0 = 0));
      UnVeilTable[label][label[LastUsed[label]]] := x;
    else
      UnVeilTable[label][label[LastUsed[label]]] := x;
    end if;
    return label[LastUsed[label]]
  end proc: # auxiliary;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ShowVeil := proc(
    label::{symbol}, # The veiling label to be shown
    $)

    description "Show the veiling variables of the veiling label <label>.";

    print(<ListVeil(label)>);
  end proc: # ShowVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ListVeil := proc(
    label::{symbol}, # The veiling label to be shown
    $)

    description "Return a list of the veiling variables labelled as <label>.";
    return sort(op(eval(UnVeilTable[label]))[2]):
  end proc: # ListVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SubsVeil := proc(
    label::{symbol}, # The label of the veiling table to be shown
    x,               # The expression to substitute
    $)

    description "Substitute the reversed veiling variables of the veiling label "
      "<label> in the expression <x>.";

    return subs[eval](op(ListTools[Reverse](ListVeil(label))), x);
  end proc: # SubsVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ForgetVeil := proc(
    label::{symbol}, # The veiling label to be forgotten
    $)

    description "Clear all the veiling variables of the veiling label <label>.";

    local i;

    # subsop(4 = NULL, eval(Signature)); # FIXME: Signature only used here
    unprotect(LastUsed);
    LastUsed[label] := 0;
    protect(LastUsed);
    UnVeilTable[label] := evaln(UnVeilTable[label]);
    forget(auxiliary);
    return NULL;
  end proc: # ForgetVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  AssignData := proc(
    x::{list, set}, # The list to be assigned
    $)

    description "Assign the list <x> to the local variable <StoredData>.";

    StoredData := x;
    return NULL;
  end proc: # AssignData

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SubsData := proc(
    x::{anything}, # The expression to substitute
    $)
    local label, tmp_x, i;

    description "Substitute the local variable <StoredData> in the expression <x>.";

    label := lhs~(op(op(LastUsed))[2..-1]);
    tmp_x := copy(x);
    if (nops(StoredData) > 0) then
      for i in label do
        tmp_x := subs[eval](StoredData, subs[eval](op(ListTools[Reverse](ListVeil(i))), tmp_x));
      end do;
    end if;

    return tmp_x;
  end proc: # SubsData

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UnAssignData := proc(
    $)

    description "Unassign the local variable <StoredData>.";

    StoredData := [];
    return NULL;
  end proc: # UnAssignData

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PermutationMatrix := proc(
    r::{Vector}, # Pivot vector (NAG-style form)
    $)

    description "Compute the LU decomposition premutation matrix provided the "
      "NAG-style the pivot vector <r>.";

    local dim, out, i;

    dim := LinearAlgebra[RowDimension](r):
    out := Matrix(dim, dim);
    for i from 1 to dim by 1 do
      out[i, r[i]] := 1;
    end do;
    return out;
  end proc: # PermutationMatrix

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  LUD := proc(
    A::{Matrix},                                            # Linear system matrix A
    Q::{symbol},                                            # Symbol for the veiling
    Strategy_Veiling::{procedure} := VeilingStrategy_n,     # Veiling strategy
    Strategy_Pivots::{procedure}  := PivotStrategy_Slength, # Pivoting strategy procedure
    Strategy_Zero::{procedure}    := ZeroStrategy_length,   # Zero recognition strategy procedure
    $)

    description "Compute the LU decomposition of a square matrix <A> using the "
      "LEM strategy <Strategy_Veiling>, the veiling symbol <Q>, the pivoting "
      "strategy <Strategy_Pivots> and the zero recognition strategy <Strategy_Zero>.";

    local P, M, n, k, i, ii, j, m, L, U, mltip, l, r, kp, p, temp, normalize,
      row, flag, z;

    # Copy the input matrix
    M := copy(A);
    n := LinearAlgebra[RowDimension](M):
    m := LinearAlgebra[ColumnDimension](M);

    # Check if the input matrix is square
    if (m <> n) then
      ERROR(
        "LULEM::LUD(...): the input matrix should be square."
        );
    end if;

    # L is lower matrix with l's as diagonal entries
    L := Matrix(LinearAlgebra[IdentityMatrix](n), shape = triangular[lower]);

    # U is upper matrix
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
        if (flag or Strategy_Pivots(M[r[i], kp], M[r[p], kp])) and
            not (Strategy_Zero(M[r[i], kp]) = 0) then
          p := i;
          break; # Once a pivot is found -- not "best"!
        end if;
      end do;

      # Only when the pivot is not equal to zero, the row elimination is performed
      # (the whole if statement). Else we will continue with the next column.

      if (Strategy_Zero(M[r[p], kp]) = 0) then

        WARNING(
          "LULEM::LUD(...): the matrix appears to be singular."
          );

      else

        if (p <> row) then
          (r[p], r[row]) := (r[row], r[p]);
        end if;

        userinfo(3, LUD, 'kp', kp, 'r', r);

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

        userinfo(3, LUD, `M`, M);

      end if;
      row := row + 1;

    end do:

    userinfo(2, LUD, `r`, r, `M`, M);

    # Seperate new M into L and U
    for i from 1 to n do
      for j from 1 to m do
        if (i <= j) then
          U[i, j] := M[r[i], j];
        else # (i > j)
          L[i, j] := M[r[i], j];
        end if
      end do;
    end do;

    # Compute the LU decomposition premutation matrix
    P := PermutationMatrix(r);

    # Return the LU decomposition and the pivot vector
    if (_nresults = 4) then
      return P, L, U, r;
    else
      return P, L, U;
    end if;
  end proc: # LUD

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Solve the linear system Ax=b given the LU decomposition (PA=LU) in NAG-style.
  # NOTE: The pivot vector r is returned from LUD function.
  SolvePivotingLU := proc(
    LU_NAG::{Matrix},              # Compact LU matrix (NAG-style form)
    r::{Vector},                   # Pivot vector
    b::{Vector},                   # Linear system vector b
    Q::{symbol},                   # Symbol for the veiling
    Strategy_Veiling::{procedure}, # Veiling strategy
    $)

    description "Solve the linear system Ax=b given the LU decomposition (PA=LU) "
      "provided the NAG_style matrix <LU_NAG>, the pivot vector <r>, the vector <b>, "
      "the veiling label <Q> and the veiling strategy <Strategy_Veiling> for the veiling.";

    local y, x, i, s, j, n, normalizer:

    userinfo(2, SolvePivotingLU, `b`, b, `r`, r);

    # Get number of rows in matrix A
    n := LinearAlgebra[RowDimension](LU_NAG);

    # Create vector for solution of Ly=Pb
    y := Vector(n);

    # Create vector for solution of Ux=y
    x := Vector(n);

    normalizer := (y) -> `if`(Strategy_Veiling(y) > 0, Veil[Q](y), y);

    # Perform forward substitution to solve Ly=Pb
    userinfo(3, SolvePivotingLU,`n`, n, `y`, y);
    y[1] := normalizer(b[r[1]]);
    for i from 2 to n do
      y[i] := normalizer(b[r[i]]) - add(normalizer(LU_NAG[i, j] * y[j]), j = 1..i-1);
    end do;

    # Perform backward substitution to solve Ux=y
    x[n] := normalizer(y[n] / LU_NAG[n, n]);
    for i from n-1 to 1 by -1 do
      s := normalizer(y[i]) - add(normalizer(LU_NAG[i, j] * x[j]), j = i+1..n);
      x[i] := normalizer(s / LU_NAG[i, i]);
    end do;

    # Return solution vector x
    if (_nresults = 1) then
      return x;
    else
      return x, r;
    end if;
  end proc: # SolvePivotingLU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Solve := proc(
    A::{Matrix},                   # Linear system matrix A
    b::{Vector},                   # Linear system vector b
    Q::{symbol, function},         # Symbol for the veiling
    Strategy_Veiling::{procedure} := VeilingStrategy_n,    # Veiling strategy
    Strategy_Pivots::{procedure} := PivotStrategy_Slength, # Pivoting strategy procedure
    Strategy_Zero::{procedure} := ZeroStrategy_length,     # Zero recognition strategy procedure
    $)

    description "Solve the linear system Ax=b using LULEM algorithm, privided the "
      "matrix <A>, the vector <b>, the veiling symbol <Q>, the veiling strategy "
      "<Strategy_Veiling>, the pivoting strategy <Strategy_Pivots> and the zero "
      "recognition strategy <Strategy_Zero>.";

    local x, P, L, U, r, LU_NAG:

    # Get LU decomposition of A
    P, L, U, r := FFLUD(A, Q, Strategy_Veiling, Strategy_Pivots, Strategy_Zero);

    # Built the LU matrix (NAG-style)
    LU_NAG := L + U - Matrix(LinearAlgebra[RowDimension](L), shape = identity);

    # Solve the linear system Ax=b given the LU decomposition (PA=LU).
    x := SolvePivotingLU(LU_NAG, r, b, Q, Strategy_Veiling);

    # Return solution vector x and the LU decomposition data
    if (_nresults = 5) then
      return x, P, L, U, r;
    elif (_nresults = 4) then
      return x, P, L, U;
    else
      return x;
    end if;
  end proc: # Solve

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #  __     __   _ _ _
  #  \ \   / /__(_) (_)_ __   __ _
  #   \ \ / / _ \ | | | '_ \ / _` |
  #    \ V /  __/ | | | | | | (_| |
  #     \_/ \___|_|_|_|_| |_|\__, |
  #                          |___/

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilingStrategy_n := proc(
    x, # Expression to be analyzed
    $)::{integer};

    description "Veiling strategy: number of indeterminates in expression <x> "
      "minus 4.";

    return nops(indets(x)) - 4;
  end proc: # VeilingStrategy_n

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilingStrategy_L:= proc(
    x, # Expression to be analyzed
    $)::{integer};

    description "Veiling strategy: length of expression <x> minus 50.";

    return length(x) - 50;
  end proc: # VeilingStrategy_L

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilingStrategy_Ls:= proc(
    x, # Expression to be analyzed
    $)::{integer};

    description "Veiling strategy: length of expression <x> minus 120.";

    return length(x) - 120;
  end proc: # VeilingStrategy_Ls

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilingStrategy_LB:= proc(
    x, # Expression to be analyzed
    $)::{integer};

    description "Veiling strategy: length of expression <x> minus 260.";

    return length(x) - 260;
  end proc: # VeilingStrategy_LB

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

    return evalb(length(SubsData(x)) - length(SubsData(y)) > 0);
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

    return evalb(evalf((abs(x) - abs(y)) > 0));
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

    local tmp;

    tmp := SubsData(x);
    if evalb((evalf(abs(tmp)) = 0.0)) then
      return length(0);
    else
      return length(tmp);
    end if;
  end proc: # ZeroStrategy_length

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ZeroStrategy_normalizer := proc(
    x::{algebraic}, # Expression to be analyzed
    $) # FIXME: What is the return type of this procedure?

    description "Zero recognition strategy: normalizer of expression <x>.";

    return Normalizer(SubsData(x));
  end proc: # ZeroStrategy_normalizer

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module:

# That's all folks!