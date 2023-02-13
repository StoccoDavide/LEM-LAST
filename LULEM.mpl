# # # # # # # # # # # # # # # # # # #
#   _    _   _ _     _____ __  __   #
#  | |  | | | | |   | ____|  \/  |  #
#  | |  | | | | |   |  _| | |\/| |  #
#  | |__| |_| | |___| |___| |  | |  #
#  |_____\___/|_____|_____|_|  |_|  #
#                                   #
# # # # # # # # # # # # # # # # # # #

# Authors: Davide Stocco and Matteo Larcher
# Date:    12/01/2023

# This is a module for the LULEM package. It contains the functions to solve
# systems of linear equations with large expressions. The module used the LU
# decomposition of the matrix of the system. The module is hopefully a better
# version of the code provided in the following PhD thesis:
#
#   Wenqin Zhou, Symbolic Computation Techniques for Solveing Large Expressions
#   Problems from Mathematics and Engineering (2007), Faculty of Graduate Studies,
#   The University of Western Ontario London, Ontario, Canada.

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

module LULEM()

  # Exported variables
  export ModuleApply,
         UnVeil,
         Veil,
         ShowVeil,
         SubsVeil,
         ForgetVeil,
         LastUsed;

  # Local variables
  local Signature,
        UnVeilTable,
        NextLabel;

  # Exported tables
  LastUsed    := table('sparse');
  UnVeilTable := table():

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # NOTUSED
  NextLabel := proc(
    label::symbol,
    vars::set,
    $)

    LastUsed[label] := LastUsed[label] + 1;
    label[LastUsed[label], `if` (nargs = 2, vars, NULL)];
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Veil := proc(
    x, # The expression to be veiled
    $)

    description "Veil an expression <x> and return a label to it.";

    local A, label;

    label := op(procname);
    LastUsed[label] := LastUsed[label] + 1;
    A := label[LastUsed[label]];
    UnVeil(A, 1)   := x;
    UnVeilTable[A] := x;
    return A;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ShowVeil := proc(
    label::symbol, # The veiling label to be shown
    $)

    description "Show the veiling variables of the veiling label <label>.";

    local i;

    for i from 1 to LastUsed[label] do
      print(label[i] = UnVeilTable[label[i]]);
    end do;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SubsVeil := proc(
    label::symbol, # The label of the veiling table to be shown
    x,             # The expression to substitute
    $)

    description "Substitute the veiling variables of the veiling label <label> "
      "in the expression <x>.";

    local i, sub_vec;

    sub_vec := [seq(0, i = 1..LastUsed[label])]:
    for i from LastUsed[label] to 1 by -1 do
      sub_vec[i] := label[i] = UnVeil(label[i]);
    end do:

    return subs(op(ListTools[Reverse](sub_vec)), x);;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Signature := proc(
    f,         # TODO
    p::posint, # TODO
    A,         # TODO
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
        Signature(UnVeil(f, 1), p, A)
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
  end:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ModuleApply := Signature;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UnVeil := proc(
    x,                                    # The expression to be unveiled
    n::{nonnegint, identical (infinity)}, # The number of levels to unveil
    $)

    description "Unveil the expression <x> up to <n> levels.";

    option remember;

    if (nargs = 1) then
      return UnVeil(x, 1);
    elif (nargs = 2) and (n = 0) then
      return x;
    elif x::atomic then
      return x;
    else
      return map(UnVeil, x, n-1);
    end if;
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ForgetVeil := proc(
    label::symbol, # The veiling label to be forgotten
    $)

    description "Clear all the veiling variables of the veiling label <label>.";

    local i;

    subsop(4 = NULL, eval(Signature));
    LastUsed[label] := 0;
    UnVeilTable     := table();

    return NULL;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   _____           _ __  __           _       _
#  | ____|_ __   __| |  \/  | ___   __| |_   _| | ___
#  |  _| | '_ \ / _` | |\/| |/ _ \ / _` | | | | |/ _ \
#  | |___| | | | (_| | |  | | (_) | (_| | |_| | |  __/
#  |_____|_| |_|\__,_|_|  |_|\___/ \__,_|\__,_|_|\___|
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SquareLUwithpivoting := proc(
  LU_NAG::Matrix,              # Compact LU matrix (NAG-style form)
  Q::symbol,                   # Symbol for the veiling
  Strategy_Veiling::procedure, # Veiling strategy
  Strategy_Pivots::procedure,  # Pivoting strategy procedure
  Strategy_Zero::procedure,    # Zero recognition strategy procedure
  $)

  local A, n, k, i, ii, j, m, L, U, mltip, l, r, kp, p, temp, normalize, row,
    flag, z;

  description "Compute the LU decomposition of a square matrix provided the "
    "NAG-style LU matrix <LU_NAG> using the LEM strategy <Strategy_Veiling>, the "
    "veiling symbol <Q>, the pivoting strategy <Strategy_Pivots> and the zero "
    "recognition strategy <Strategy_Zero>.";

  # Copy the input matrix
  A := copy(LU_NAG);
  n := LinearAlgebra[RowDimension](A):
  m := LinearAlgebra[ColumnDimension](A);

  # Check if the input matrix is square
  if (m <> n) then
    ERROR(
      "LULEM::SquareLUwithpivoting(...): the input matrix should be square."
      );
  end if;

  # L is lower matrix with l's as diagonal entries
  L := Matrix(LinearAlgebra[IdentityMatrix](n), shape = triangular[lower]);

  # U is upper matrix.
  U := Matrix(n,n,shape = triangular[upper]):

  # Create pivot vector
  r := Vector(n, k -> k);

  normalize := z -> `if`(Strategy_Veiling(z) > 0, LULEM:-Veil[Q](z), z);

  # Loop over columns and find the pivot element among the entries
  # A[r[kp],kp], A[r[kp+1],kp],	A[r[n],kp]
  # according to different pivoting strategies.
  # Interchange entries in pivot vector.

  row := 1;
  for kp from 1 to m while row <= n do
    p := row;
    flag := evalb(Strategy_Zero(A[r[row],kp]) = 0);
    for i from row+1 to n do
        # Once a pivot is found -- not "best"!
        if (flag or Strategy_Pivots(A[r[i], kp], A[r[p],kp])) and
            not (Strategy_Zero(A[r[i], kp]) = 0) then
          p := i;
        break;
        end if;
    end do;

    # Only when the pivot is not equal to zero, the row elimination is performed
    # (the whole if statement). Else we will continue with the next column.

    if (Strategy_Zero(AEr[p]) = 0) then
        WARNING(
          "LULEM::SquareLUwithpivoting(...): the matrix appears to be singular."
          );
    else
      if (p <> row) then
        (r[p], r[row]) := (r[row], r[p]);
      end if;
      userinfo(3, SquareLUwithpivoting, 'kp', kp, 'r', r);

      # Now do Gauss elimination steps to get new A and also keep L information
      # in new A. Packing everything into the A matrix during the computation of
      # the factors. The reason is that this code can be ported to restricted
      # memory environments, or equally, applied to very large matrices in a
      # large memory environment.

      for i from row+1 to n do
        mltip := normalize(A[r[i], kp]/A[r[row], kp]);
        A[r[i] ,kp] := mltip;
        for j from kp+1 to m do
          z := A[r[i],j] - mltip * A[r[row], j] ;
          A[r[i], j] := `if`(Strategy_Veiling(z) > 0, LULEM:-Veil[Q](z), z);
        end do;
      end do;
      userinfo(3, SquareLUwithpivoting, `A`, A);
    end if;
    row := row + 1;
  end do:
  userinfo(2, SquareLUwithpivoting, `r`, r, `A`, A);

  # Seperate new A into L and U
  for i from 1 to n do
    for j from 1 to m do
      if i <= j then
        U[i, j] := A[r[i], j];
      else
        L[i, j] := A[r[i], j];
      end if
    end do;
  end do;

  # Return the LU decomposition and the pivot vector
  return L, U, r;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Solve the linear system Ax=b given the LU decomposition (PA=LU).
# NOTE: The pivot vector r is returned from SquareLUwithpivoting function.
SolveSquareLUpivot := proc(
  A::Matrix,                   # Linear system matrix A
  r::Vector,                   # Pivot vector
  b::Vector,                   # Linear system vector b
  Q::symbol,                   # Symbol for the veiling
  Strategy_Veiling::procedure, # Veiling strategy
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

  normalizer := y -> `if`(Strategy_Veiling(y) > 0, LULEM:-Veil[Q](y), y);

  # Perform forward substitution to solve Ly=Pb
  userinfo(3, SolveSquareLUpivot,`n`, n, `y`, y);
  y[1] := normalizer( b[r[1]] );
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
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LULEMsolve := proc(
  A::Matrix,                   # Linear system matrix A
  b::Vector,                   # Linear system vector b
  Q::symbol,                   # Symbol for the veiling
  Strategy_Veiling::procedure, # Veiling strategy
  Strategy_Pivots::procedure,  # Pivoting strategy procedure
  Strategy_Zero::procedure,    # Zero recognition strategy procedure
  $)

  local L, U, r, x:

  description "Solve the linear system Ax=b using LULEM algorithm, privided the "
    "matrix <A>, the vector <b>, the veiling symbol <Q>, the veiling strategy "
    "<Strategy_Veiling>, the pivoting strategy <Strategy_Pivots> and the zero "
    "recognition strategy <Strategy_Zero>.";

  # Get LU decomposition of A
  L, U, r := SquareLUwithpivoting(A, Q, Strategy_Veiling, Strategy_Pivots, Strategy_Zero);

  # Solve the linear system Ax=b given the LU decomposition (PA=LU).
  x := SolveSquareLUpivot(A, r, b, Q, Strategy_Veiling);

  # Return solution vector x
  return x;

end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#  __     __   _ _ _
#  \ \   / /__(_) (_)_ __   __ _
#   \ \ / / _ \ | | | '_ \ / _` |
#    \ V /  __/ | | | | | | (_| |
#     \_/ \___|_|_|_|_| |_|\__, |
#                          |___/

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LEMStrategy_n := proc(
  x, # Expression 1 to be analyzed
  $)

  description "Veiling strategy: number of indeterminates in x minus 4.";

  return nops(indets(x)) - 4;
end proc;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LEMStrategy_L:= proc(
  x, # Expression 1 to be analyzed
  $)

  description "Veiling strategy: length of x minus 50.";

  return length(x) - 50;
end proc;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LEMStrategy_Ls:= proc(
  x, # Expression 1 to be analyzed
  $)

  description "Veiling strategy: length of x minus 120.";

  return length(x) - 120;
end proc;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LEMStrategy_LB:= proc(
  x, # Expression 1 to be analyzed
  $)

  description "Veiling strategy: length of x minus 260.";

  return length(x) - 260;
end proc;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   ____  _            _   _
#  |  _ \(_)_   _____ | |_(_)_ __   __ _
#  | |_) | \ \ / / _ \| __| | '_ \ / _` |
#  |  __/| |\ V / (_) | |_| | | | | (_| |
#  |_|   |_| \_/ \___/ \__|_|_| |_|\__, |
#                                  |___/

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotStrategy_Llength := proc(
  x, # Expression 1 to be analyzed
  y, # Expression 2 to be analyzed
  $)

  description "Pivoting strategy: choose the pivot with the largest length";

  return evalb(length(x) - length(y) > 0);
end proc;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotStrategy_Slength := proc(
  x, # Expression 1 to be analyzed
  y, # Expression 2 to be analyzed
  $)

  description "Pivoting strategy: choose the pivot with the smallest length";

  return evalb(length(x) - length(y) < 0);
end proc;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotStrategy_Lindets := proc (
  x, # Expression 1 to be analyzed
  y, # Expression 2 to be analyzed
  $)

  description "Pivoting strategy: choose the pivot with the largest number of "
    "indeterminates";

  return evalb((nops(indets(x)) - nops(indets(y))) > 0);
end proc;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotStrategy_Sindets := proc (
  x, # Expression 1 to be analyzed
  y, # Expression 2 to be analyzed
  $)

  description "Pivoting strategy: choose the pivot with the smallest number of "
    "indeterminates";

  return evalb((nops(indets(x)) - nops(indets(y))) < 0);
end proc;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotStrategy_numeric := proc (
  x, # Expression 1 to be analyzed
  y, # Expression 2 to be analyzed
  $)

  description "Pivoting strategy: choose the pivot with the largest numeric "
    "value";

  return evalb((abs(x) - abs(y)) > 0);
end proc;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   _____
#  |__  /___ _ __ ___
#    / // _ \ '__/ _ \
#   / /|  __/ | | (_) |
#  /____\___|_|  \___/
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ZeroStrategy_length := proc(
  x, # Expression to be analyzed
  $)

  description "Zero recognition strategy: length";

  return length(x)
end proc;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ZeroStrategy_normalizer := proc(
  x, # Expression to be analyzed
  $)

  description "Zero recognition strategy: normalizer";

  return Normalizer(x)
end proc;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# That's all folks!
