# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                        _    _   _ _     _____ __  __                        #
#                       | |  | | | | |   | ____|  \/  |                       #
#                       | |  | | | | |   |  _| | |\/| |                       #
#                       | |__| |_| | |___| |___| |  | |                       #
#                       |_____\___/|_____|_____|_|  |_|                       #
#          LU and QR decomposition with Large Expressions Management          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Authors of the current version:
#  Davide Stocco (University of Trento)
#  Matteo Larcher (University of Trento)
#  Enrico Bertolazzi (University of Trento)
#
# Authors of the original code:
#   Wenqin Zhou (University of Western Ontario) - Former affiliation
#   David J. Jeffrey (University of Western Ontario)
#   Jacques Carette (McMaster University)
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

  export  AssignData,
          SubsData,
          ForgetData,
          SetVerbosity,
          Permutatiomnatrices,
          LU,
          LUPivoting,
          SolveLU,
          FFLU,
          FF2LU,
          SolveFFLU,
          QR,
          # TODO: SolveQR,
          # TODO: FFQR,
          # TODO: SolveFFQR,
          VeilingStrategy_n,
          VeilingStrategy_L,
          VeilingStrategy_Ls,
          VeilingStrategy_LB,
          PivotingStrategy_Llength,
          PivotingStrategy_Slength,
          PivotingStrategy_Lindets,
          PivotingStrategy_Sindets,
          PivotingStrategy_numeric,
          ZeroStrategy_length,
          ZeroStrategy_Dlength,
          ZeroStrategy_normalizer;

  local   ModuleLoad,
          ModuleUnload,
          auxiliary,
          SolvePivotingLU,
          InitLULEM,
          UnVeilTable,
          StoredData,
          lib_base_path,
          Verbose;

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

    local i;

    printf(cat(
      "'LULEM' module version 1.1, ",
      "BSD 3-Clause License - Copyright (C) 2023, D. Stocco, M. Larcher, ",
      "E. Bertolazzi, W. Zhou, D. J. Jeffrey, J. Carette and R. M. Corless.\n"
      ));

    lib_base_path := null;
    for i in [libname] do
      if (StringTools[Search]("LULEM", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if (lib_base_path = null) then
      error "Cannot find 'LULEM' module" ;
    end if;

    InitLULEM();

    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ModuleUnload := proc()

    description "Module 'LULEM' module unload procedure";

    ForgetData();
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  InitLULEM := proc()

    description "Initialize 'LULEM' module internal variables";

    StoredData := [];
    Verbose    := false;

    return NULL;
  end proc: # InitLULEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  AssignData := proc(
    x::{list, set}, # The list to be assigned
    $)::{nothing};

    description "Assign the data list <x> to the local variable <StoredData>.";

    StoredData := x;

    return NULL;
  end proc: # AssignData

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SubsData := proc(
    x::{anything},
    $)::{anything};

    description "Substitute the local variable <StoredData> in the expression <x>.";

    local out, i;

    out := copy(x);
    if (nops(StoredData) > 0) then
      for i in label do
        out := subs[eval](StoredData, out);
      end do;
    end if;

    return out;
  end proc: # SubsData

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ForgetData := proc(
    $)::{nothing};

    description "Unassign the local variable <StoredData>.";

    StoredData := [];
  end proc: # ForgetData

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Permutatiomnatrices := proc(
    r::{Vector},
    c::{Vector},
    $)

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

  LUPivoting := proc(
    k::{integer},
    M::{Matrix},
    V::{symbol},
    r::{Vector(integer)},
    c::{Vector(integer)},
    VeilingStrategy::{procedure},
    PivotingStrategy::{procedure},
    ZeroStrategy::{procedure},
    $)

    description "Compute the LU decomposition pivot vector provided the step <k>, "
      "the temporary LU (NAG) matrix <M>, the veiling symbol <V>, the rows the "
      "pivot vector <r>, the columns the pivot vector <c>, the veiling strategy "
      "<VeilingStrategy>, the pivoting strategy <PivotingStrategy> and the zero "
      "recognition strategy <ZeroStrategy>.";

    local Mij, Mkk, m, n, i, j, apply_veil, pivot_is_zero, Mij_is_zero, z, tmp;

    # Extract dimensions
    m, n := LinearAlgebra[Dimensions](M):

    # Check if to veil or not
    apply_veil := (z) -> `if`(VeilingStrategy(z) > 0, LEM[Veil][V](z), z);

    # Check if M[r[k],c[k]] = 0, if not true it is the pivot
    Mkk := M[k, k];
    try
      pivot_is_zero := evalb(ZeroStrategy(LEM[SubsVeil](Mkk, V)) = 0);
      #pivot_is_zero := evalb(ZeroStrategy(Normalizer(Mkk)) = 0);
    catch:
      print("Mkk: Division by 0 or numerical exception.\n");
      print(Mkk);
      pivot_is_zero := true;
    end try;

    # Iterate over the columns and rows
    for j from k to n do
      for i from k to m do

        # Look for a non-zero pivot
        Mij := M[i, j];
        try
          Mij_is_zero := evalb(ZeroStrategy(LEM[SubsVeil](Mij, V)) = 0);
          #Mij_is_zero := evalb(ZeroStrategy(Normalizer(Mij)) = 0);
        catch:
          if Verbose then
            print("Mij: Division by 0 or numerical exception.\n");
            print(Mij);
          end if;
          Mij_is_zero := true;
        end try;

        if not Mij_is_zero and (pivot_is_zero or PivotingStrategy(Mij, Mkk)) then
          # A better pivot is found
          pivot_is_zero := false;
          if (i <> k) then
            (r[i], r[k])    := (r[k], r[i]);
            M[[i,k], 1..-1] := M[[k,i], 1..-1];
          end if;
          if (j <> k) then
            (c[j], c[k])    := (c[k], c[j]);
            M[1..-1, [j,k]] := M[1..-1 ,[k,j]];
          end if;
          Mkk := Mij;
        end if;

      end do;
    end do;

    return pivot_is_zero, Mkk;
  end proc: # LUPivoting

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SetVerbosity := proc(
    x::{boolean},
    $)::{nothing};

    description "Set the verbosity of the package to <x>.";

    Verbose := x;

    return NULL;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  LU := proc(
    A::{Matrix},
    V::{symbol},
    VeilingStrategy::{procedure}  := VeilingStrategy_n,
    PivotingStrategy::{procedure} := PivotingStrategy_Slength,
    ZeroStrategy::{procedure}     := ZeroStrategy_length,
    $)::{table};

    description "Compute the LU decomposition of a square matrix <A> using the "
      "veiling strategy <VeilingStrategy>, the veiling symbol <V>, the pivoting "
      "strategy <PivotingStrategy> and the zero recognition strategy "
      "<ZeroStrategy>.";

    local M, L, U, Mkk, m, n, mn, k, rnk, r, c, apply_veil, pivot_is_zero, Mij_is_zero, tmp;

    m, n := LinearAlgebra[Dimensions](A):

    # Create pivot vector
    r := Vector(m, k -> k);
    c := Vector(n, k -> k);

    # Check if to veil or not
    apply_veil := (z) -> `if`(VeilingStrategy(z) > 0, LEM[Veil][V](z), z);

    # Perform Gaussian elimination
    M   := copy(A);
    mn  := min(m, n);
    rnk := mn;
    for k from 1 to (mn - 1) do
      if Verbose then
        printf("LULEM::LU(...): processing %d-th row.\n", k);
      end;
      pivot_is_zero, Mkk := LUPivoting(
        k, M, V, r, c, VeilingStrategy, PivotingStrategy, ZeroStrategy
      );

      if pivot_is_zero then
        rnk := k;
        if Verbose then
          WARNING("LULEM::LU(...): the matrix appears not full rank.");
        end;
        break;
      end if;

      if Verbose then
        print("LULEM::LU(...): pivot:", Mkk);
      end;

      # Shur complement
      tmp         := [k+1..-1];
      M[tmp, k]   := M[tmp, k] / Mkk;
      M[tmp, tmp] := apply_veil~(Normalizer~(M[tmp, tmp] - M[tmp, k].M[k, tmp]));
    end do:

    L := Matrix(M[1..m, 1..m], shape = triangular[lower, unit]);
    U := Matrix(M, shape = triangular[upper]);

    # Return the LU decomposition
    return table([
      "methos"   = "LU",
      "L"        = L,
      "U"        = U,
      "V"        = V,
      "r"        = r,
      "c"        = c,
      "rank"     = rnk,
      "L_length" = length(L),
      "U_length" = length(U),
      "V_length" = length(LEM[ListVeil](V))
    ]);
  end proc: # LU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SolveLU := proc(
    T::{table},
    b::{Vector},
    V::{symbol, function},
    VeilingStrategy::{procedure} := VeilingStrategy_n,
    $)

    description "Solve the linear system Ax=b using LU decomposition <T>, "
      "provided the vector <b>, the veiling symbol <V> and the veiling strategy "
      "<VeilingStrategy>.";

    local L, U, r, c, m, n, rnk, LU_NAG, apply_veil, y, i, j, x, s:

    # Extract the LU decomposition
    L := T["L"];
    U := T["U"];
    r := T["r"];
    c := T["c"];

    # Get linear system dimension
    m, n := LinearAlgebra[Dimensions](A);

    # Check if the linear system is consistent
    assert(
      LinearAlgebra[RowDimension](L) = LinearAlgebra[ColumnDimension](b),
      "LULEM::SolveLU(...): inconsistent linear system."
    );

    # Create vector for solution of Ly=b[r]
    y := Vector(m);

    # Create vector for solution of Ux=y[c]
    x := Vector(n);

    # Get permutation matrices
    PP, QQ := Permutatiomnatrices(r, c);

    # Create a normalizer function
    apply_veil := (y) -> `if`(VeilingStrategy(y) > 0, LEM[Veil][V](y), y);

    # Perform forward substitution to solve Ly=b[r]
    y[1] := apply_veil(b[r[1]]);
    for i from 2 to m do
      y[i] := apply_veil(b[r[i]]) - add(apply_veil(L[i, j] * y[j]), j = 1..i-1);
    end do;

    print(simplify(LEM[SubsVeil]([L.y, PP.b])));

    # Perform backward substitution to solve Ux[c]=y
    x[c[n]] := apply_veil(y[n] / U[n, n]);
    for i from (n - 1) to 1 by -1 do
      s := apply_veil(y[i]) - add(apply_veil(U[i, j] * x[c[j]]), j = i+1..n);
      x[c[i]] := apply_veil(s / U[i, i]);
    end do;

    print(simplify(LEM[SubsVeil]([U.x, QQ.y])));

    # Return outputs
    return x;
  end proc: # SolveLUD

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FFLU := proc(
    A::{Matrix},
    V::{symbol},
    VeilingStrategy::{procedure}  := VeilingStrategy_n,
    PivotingStrategy::{procedure} := PivotStrategy_Slength,
    ZeroStrategy::{procedure}     := ZeroStrategy_length,
    $)::{table};

    description "Compute the Fracton-Free LU decomposition of a square matrix "
      "<A> using the veiling strategy <VeilingStrategy>, the veiling symbol <V>, "
      "the pivoting strategy <PivotingStrategy> and the zero recognition strategy "
      "<ZeroStrategy>.";

    local SS, M, Mkk, m, n, mn, i, j, k, ri, rk, rnk, r, c, apply_veil,
      pivot_is_zero, Mij_is_zero, z, tmp;

    m, n := LinearAlgebra[Dimensions](A):

    # Create pivot vector
    r := Vector(m, k -> k);
    c := Vector(n, k -> k);

    # check if Veil or not
    apply_veil := z -> `if`(VeilingStrategy(z) > 0, LEM[Veil][V](z), z);

    # Gauss Elimination main loop
    M   := copy(A);
    mn  := min(m, n);
    rnk := mn;
    SS  := Vector(mn);
    for k from 1 to mn-1 do
      if Verbose then
        printf("LULEM::FFLU(...): processing %d-th row.\n", k)
      end;
      pivot_is_zero, Mkk := LUPivoting(
        k, M, V, r, c, VeilingStrategy, PivotingStrategy, ZeroStrategy
      );

      if pivot_is_zero then
        rnk := k;
        if Verbose then
          WARNING("LULEM::LU(...): the matrix appears not full rank.");
        end;
        break;
      end if;

      if Verbose then
        print("LULEM::FFLU(...): pivot:", Mkk);
      end;

      SS[k] := Mkk;
      # Scaled Shur complement
      tmp        := [k+1..-1];
      M[tmp,tmp] := apply_veil~(simplify(Mkk*M[tmp,tmp]-M[tmp,k].M[k,tmp],size));
      #M[tmp,tmp] := apply_veil~(Normalizer~(Mkk*M[tmp,tmp]-M[tmp,k].M[k,tmp]));
    end do:

    # Return the FFLU decomposition
    return table([
      "method"   = "FFLU",
      "M"        = M,
      "V"        = V,
      "S"        = SS,
      "r"        = r,
      "c"        = c,
      "rank"     = rnk,
      "M_length" = length(M),
      "S_length" = length(SS),
      "V_length" = length(LEM[ListVeil](V))
    ]);
  end proc: # FFLU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FF2LU := proc(
    T::{table},
    $)

    local M, SS, r, c, rk, n, m, L, DG, L_list, D_list, U, i, j, k;

    M  := T["M"];
    SS := T["S"];
    r  := T["r"];
    c  := T["c"];
    rk := T["rank"];

    m, n := LinearAlgebra[Dimensions](M):

    L_list := [];
    D_list := [];

    #
    #   +             + +         +      +            +       +             + +            +
    #   |  1          | | 1       |      | 1          |       | 1           | | 1          |
    #   | -x  1       | |   q     | ==>  | -x q       |  ==>  |   1/q       | | x  1       |
    #   | -x     1    | |     q   |      | -x    q    |       |      1/q    | | x     1    |
    #   | -x        1 | |       q |      | -x       q |       |         1/q | | x        1 |
    #   +             + +         +      +            +       +             + +            +
    #
    for k from 1 to rk do
      L            := Matrix( m, m, shape = triangular[lower, unit]);
      L[k+1..-1,k] := M[k+1..-1,k];
      L_list       := [ op(L_list), L ];
      DG           := Matrix( LinearAlgebra[IdentityMatrix](m), shape = diagonal );
      for j from k+1 to m do
        DG[j,j] := SS[k];
      end;
      D_list := [ op(D_list), DG ];
    end;
    return L_list, D_list, Matrix(M, shape = triangular[upper]);
  end proc: # FF2LU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SolveFFLU := proc(
    T::{table},
    b::{Vector},
    V::{symbol, function},
    VeilingStrategy::{procedure} := VeilingStrategy_n,
    $)

    description "Solve the linear system Ax=b using FFLU decomposition <T>, "
      "provided the vector <b>, the veiling symbol <V> and the veiling strategy "
      "<VeilingStrategy>.";

  end proc: # SolveFFLU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  QR := proc(
    A::{Matrix},
    V::{symbol},
    VeilingStrategy::{procedure}  := VeilingStrategy_n,
    PivotingStrategy::{procedure} := PivotStrategy_Slength,
    ZeroStrategy::{procedure}     := ZeroStrategy_length,
    $)::{table};

    description "Compute the Householder QR decomposition of a square matrix <A> "
      "using the veiling strategy <VeilingStrategy>, the veiling symbol <V>, "
      "the pivoting strategy <PivotingStrategy> and the zero recognition strategy "
      "<ZeroStrategy>.";

    local m, n, z, k, Q, R, M, norm_x, s, u1, w, tau;

    # Extract the dimensions of the matrix A
    m, n := LinearAlgebra[Dimensions](A):

    # Check if the matrix A valid
    assert(
      m >= n,
      "LULEM::QR(...): invalid matrix A(m,n) detected, got (m >= n)."
    );

    # Initialize some variables
    z := Vector(n, k -> 0);
    Q := Matrix(LinearAlgebra[IdentityMatrix](m)); # Orthogonal transformation so far
    R := copy(A);                                  # Transformed matrix so far

    # Compute the Householder QR decomposition with veiling
    for k from 1 to n do

      # Find H = I-tau*w*w’ to put zeros below R[j,j]
      norm_x := LinearAlgebra[Norm](R[k..-1,k], 2);
      s      := -sign(R[k, k]);
      u1     := R[k, k] - s*norm_x;
      w      := R[k..-1,k]/u1;
      w[1]   := 1;
      tau    := -s*u1/norm_x;

      # Update R = HR and Q = QH
      R[k..-1, 1..-1] := R[k..-1, 1..-1] - tau.w.LinearAlgebra[Transpose](w).R[k..-1,1..-1];
      Q[1..-1, k..-1] := Q[1..-1, k..-1] - tau.Q[1..-1, k..-1].w.LinearAlgebra[Transpose](w);

    end do;

    # Return the QR decomposition
    return table([
      "method"   = "QR",
      "Q"        = Q,
      "R"        = R,
      "V"        = V,
      "Q_length" = length(Q),
      "R_length" = length(R),
      "V_length" = length(LEM[ListVeil](V))
    ]);
  end proc: # QR

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
(*
  SolveQR := proc(
    T::{table},
    b::{Vector},
    V::{symbol, function},
    VeilingStrategy::{procedure} := VeilingStrategy_n,
    $)

    description "Solve the linear system Ax=b using QR decomposition <T>, "
      "provided the vector <b>, the veiling symbol <V> and the veiling strategy "
      "<VeilingStrategy>.";

  end proc: # SolveQR
*)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
(*
  FFQR := proc(
    A::{Matrix},
    V::{symbol},
    VeilingStrategy::{procedure}  := VeilingStrategy_n,
    PivotingStrategy::{procedure} := PivotStrategy_Slength,
    ZeroStrategy::{procedure}     := ZeroStrategy_length,
    $)::{table};

    description "Compute the Fracton-Free QR decomposition of a square matrix "
      "<A> using the veiling strategy <VeilingStrategy>, the veiling symbol <V>, "
      "the pivoting strategy <PivotingStrategy> and the zero recognition strategy "
      "<ZeroStrategy>.";

  end proc: # QR
*)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
(*
  SolveFFQR := proc(
    T::{table},
    b::{Vector},
    V::{symbol, function},
    VeilingStrategy::{procedure} := VeilingStrategy_n,
    $)

    description "Solve the linear system Ax=b using FFQR decomposition <T>, "
      "provided the vector <b>, the veiling symbol <V> and the veiling strategy "
      "<VeilingStrategy>.";

  end proc: # SolveFFQR
*)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #  __     __   _ _ _
  #  \ \   / /__(_) (_)_ __   __ _
  #   \ \ / / _ \ | | | '_ \ / _` |
  #    \ V /  __/ | | | | | | (_| |
  #     \_/ \___|_|_|_|_| |_|\__, |
  #                          |___/

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilingStrategy_n := proc(
    x::{algebraic},
    $)::{integer};

    description "Veiling strategy: number of indeterminates in expression <x> "
      "minus 4.";

    return nops(indets(x)) - 4;
  end proc: # VeilingStrategy_n

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilingStrategy_L := proc(
    x::{algebraic},
    $)::{integer};

    description "Veiling strategy: length of expression <x> minus 50.";

    return length(x) - 50;
  end proc: # VeilingStrategy_L

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilingStrategy_Ls := proc(
    x::{algebraic},
    $)::{integer};

    description "Veiling strategy: length of expression <x> minus 120.";

    return length(x) - 120;
  end proc: # VeilingStrategy_Ls

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilingStrategy_LB := proc(
    x::{algebraic},
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

  PivotingStrategy_Llength := proc(
    x::{algebraic},
    y::{algebraic},
    $)::{boolean};

    description "Pivoting strategy: choose the pivot with the largest length "
      "between expressions <x> and <y>.";

    return evalb(length(SubsData(x)) - length(SubsData(y)) > 0);
  end proc: # PivotingStrategy_Llength

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotingStrategy_Slength := proc(
    x::{algebraic},
    y::{algebraic},
    $)::{boolean};

    description "Pivoting strategy: choose the pivot with the smallest length "
      "between expressions <x> and <y>.";

    return evalb(length(x) - length(y) < 0);
  end proc: # PivotingStrategy_Slength

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotingStrategy_Lindets := proc (
    x::{algebraic},
    y::{algebraic},
    $)::{boolean};

    description "Pivoting strategy: choose the pivot with the largest number of "
      "indeterminates between expressions <x> and <y>.";

    return evalb((nops(indets(x)) - nops(indets(y))) > 0);
  end proc: # PivotingStrategy_Lindets

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotingStrategy_Sindets := proc (
    x::{algebraic},
    y::{algebraic},
    $)::{boolean};

    description "Pivoting strategy: choose the pivot with the smallest number of "
      "indeterminates between expressions <x> and <y>.";

    return evalb((nops(indets(x)) - nops(indets(y))) < 0);
  end proc: # PivotingStrategy_Sindets

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotingStrategy_numeric := proc (
    x::{algebraic},
    y::{algebraic},
    $)::{boolean};

    description "Pivoting strategy: choose the pivot with the largest numeric "
      "value between expressions <x> and <y>.";

    return evalb(evalf((abs(x) - abs(y)) > 0));
  end proc: # PivotingStrategy_numeric

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #   _____
  #  |__  /___ _ __ ___
  #    / // _ \ '__/ _ \
  #   / /|  __/ | | (_) |
  #  /____\___|_|  \___/
  #

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ZeroStrategy_length := proc(
    x::{algebraic},
    $)::{integer};

    description "Zero recognition strategy: length of expression <x>.";

    return length(x);
  end proc: # ZeroStrategy_length

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ZeroStrategy_Dlength := proc(
    x::{algebraic},
    $)::{integer};

    description "Zero recognition strategy: length of expression <x> substituted "
      "with assigned data.";

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
    x::{algebraic},
    $) # FIXME: What is the return type of this procedure?

    description "Zero recognition strategy: normalizer of expression <x>.";

    return Normalizer(SubsData(x));
  end proc: # ZeroStrategy_normalizer

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module:

# That's all folks!
