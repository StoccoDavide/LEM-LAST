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

  export  AssignData,
          SubsData,
          ForgetData,
          SetVerbosity,
          PermutationMatrices,
          SolveLinearSystem,
          LU,
          FFLU,
          FF2LU,
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
          PivotingStrategy_numeric;

  local   ModuleLoad,
          ModuleUnload,
          auxiliary,
          DoPivoting,
          SolveLU,
          SolveFractionalFreeLU,
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

  AssignData := proc( x::{list, set}, $ )::{nothing};
    # x::{list, set} The list to be assigned
    description "Assign the data list <x> to the local variable <StoredData>.";
    StoredData := x;
    return NULL;
  end proc: # AssignData

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SubsData := proc( x::{anything}, $ )::{anything};

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

  ForgetData := proc( $ )::{nothing};
    description "Unassign the local variable <StoredData>.";
    StoredData := [];
  end proc: # ForgetData

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

  DoPivoting := proc(
    k::{integer},
    M::{Matrix},
    V::{symbol},
    r::{Vector(integer)},
    c::{Vector(integer)},
    VeilingStrategy::{procedure},
    PivotingStrategy::{procedure},
    $)

    description "Compute the LU decomposition pivot vector provided the step <k>, "
                "the temporary LU (NAG) matrix <M>, the veiling symbol <V>, the rows the "
                "pivot vector <r>, the columns the pivot vector <c>, the veiling strategy "
                "<VeilingStrategy> and the pivoting strategy <PivotingStrategy>.";

    local Mij, Mkk, uMij, uMkk, LV, m, n, i, j, ii, jj,
          apply_veil, apply_unveil, max_size_expression, pivot_is_zero, Mij_is_zero, z, tmp;

    max_size_expression := 1000;

    # Extract dimensions
    m, n := LinearAlgebra[Dimensions](M):

    # Check if to veil or not
    apply_veil   := (z) -> `if`( VeilingStrategy(z), LEM[Veil][V](z), z);
    LV           := LEM[ListVeil](V,true);
    #apply_unveil := (z) -> Normalizer(subs[eval](op( LV ), z ));
    apply_unveil := (z) -> subs( op( LV ), z );

    # Check if M[r[k],c[k]] = 0, if not true it is the pivot
    Mkk := M[k,k];
    ii  := k;
    jj  := k;
    try
      Mkk           := Normalizer(Mkk);
      pivot_is_zero := evalb( Mkk = 0 );
      if not pivot_is_zero then
        uMkk := apply_unveil(Mkk);
        if length(uMkk) < max_size_expression then
          uMkk := Normalizer(eval(uMkk));
          pivot_is_zero := evalb( uMkk = 0 );
        end if;
      end;
    catch:
      print("Mkk: Division by 0 or numerical exception.\n");
      print(uMkk);
      pivot_is_zero := true;
    end try;

    # Iterate over the columns and rows
    for j from k to n do
      for i from k to m do
        # Look for a non-zero pivot
        Mij := M[i,j];
        try
          Mij         := Normalizer(Mij);
          Mij_is_zero := evalb( Mij = 0 );
          if not Mij_is_zero then
            uMij := apply_unveil(Mij);
            if length(uMij) < max_size_expression then
              # timelimit required because sometime Normalizer stuck
              uMij        := timelimit( 0.5, eval(Normalizer(eval(uMij))) );
              Mij_is_zero := evalb( uMij = 0 );
            end if;
          end;
        catch:
          if Verbose then
            print("Mij: Division by 0 or numerical exception.\n");
            print(Mij);
          end if;
          Mij_is_zero := true;
        end try;

        if not Mij_is_zero then
          # found non zero pivot, check if it is better
          if pivot_is_zero or PivotingStrategy(Mij,Mkk) then
            # A better pivot is found
            pivot_is_zero := false;
            Mkk := Mij;
            ii  := i;
            jj  := j;
            if nops(indets(Mkk)) = 0 then
              break 2;
            end if;
          end
        end if;

      end do;
    end do;

    if not pivot_is_zero then
      i := ii;
      j := jj;
      # A better pivot is found
      if (i <> k) then
        (r[i], r[k])    := (r[k], r[i]);
        M[[i,k], 1..-1] := M[[k,i], 1..-1];
      end if;
      if (j <> k) then
        (c[j], c[k])    := (c[k], c[j]);
        M[1..-1, [j,k]] := M[1..-1 ,[k,j]];
      end if;
    end if;

    return pivot_is_zero, Mkk;
  end proc: # DoPivoting

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SetVerbosity := proc( x::{boolean}, $ )::{nothing};
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
    $)::{table};

    description "Compute the LU decomposition of a square matrix <A> using the "
                "veiling strategy <VeilingStrategy>, the veiling symbol <V> and the pivoting "
                "strategy <PivotingStrategy>.";

    local M, L, U, Mkk, m, n, mn, k, rnk, r, c, apply_veil, pivot_is_zero, Mij_is_zero, tmp;

    m, n := LinearAlgebra[Dimensions](A):

    # Create pivot vector
    r := Vector(m, k -> k);
    c := Vector(n, k -> k);

    # Check if to veil or not
    apply_veil := (z) -> `if`(VeilingStrategy(z), LEM[Veil][V](z), z);

    # Perform Gaussian elimination
    M   := copy(A);
    mn  := min(m, n);
    rnk := mn;
    for k from 1 to (mn - 1) do
      if Verbose then
        printf(
          "LULEM::LU(...): processing %d-th row. Length %d:%d\n",
          k, length(convert(M,list)), length(LEM[ListVeil](V))
        );
      end;
      pivot_is_zero, Mkk := DoPivoting( k, M, V, r, c, VeilingStrategy, PivotingStrategy );

      if pivot_is_zero then
        rnk := k;
        if Verbose then
          WARNING( "LULEM::LU(...): the matrix appears not full rank." );
        end;
        break;
      end if;

      if Verbose then
        print( "LULEM::LU(...): pivot:", Mkk );
      end;

      # Shur complement
      tmp         := [k+1..-1];
      M[tmp, k]   := apply_veil~(Normalizer~(M[tmp, k] / Mkk));
      M[tmp, tmp] := apply_veil~(Normalizer~(M[tmp, tmp] - M[tmp, k].M[k, tmp]));
      #M[tmp, k]   := apply_veil~(simplify~(M[tmp, k] / Mkk));
      #M[tmp, tmp] := apply_veil~(simplify~(M[tmp, tmp] - M[tmp, k].M[k, tmp]));
    end do:

    L := Matrix(M[1..m, 1..m], shape = triangular[lower, unit]);
    U := Matrix(M, shape = triangular[upper]);

    # Return the LU decomposition
    return table([
      "method"   = "LU",
      "L"        = L,
      "U"        = U,
      "V"        = V,
      "r"        = r,
      "c"        = c,
      "rank"     = rnk,
      "L_length" = length(convert(L,list)),
      "U_length" = length(convert(U,list)),
      "V_length" = length(LEM[ListVeil](V))
    ]);
  end proc: # LU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FFLU := proc(
    A::{Matrix},
    V::{symbol},
    VeilingStrategy::{procedure}  := VeilingStrategy_n,
    PivotingStrategy::{procedure} := PivotStrategy_Slength,
    $)::{table};

    description "Compute the Fracton-Free LU decomposition of a square matrix "
                "<A> using the veiling strategy <VeilingStrategy>, the veiling symbol <V> and "
                "the pivoting strategy <PivotingStrategy>.";

    local SS, M, Mkk, m, n, mn, i, j, k, ri, rk, rnk, r, c, apply_veil,
          pivot_is_zero, Mij_is_zero, z, tmp, bot, top;

    m, n := LinearAlgebra[Dimensions](A):

    # Create pivot vector
    r := Vector(m, k -> k);
    c := Vector(n, k -> k);

    # check if Veil or not
    apply_veil := z -> `if`(VeilingStrategy(z), LEM[Veil][V](z), z);

    # Gauss Elimination main loop
    M   := copy(A);
    mn  := min(m, n);
    rnk := mn;
    SS  := Vector(mn);
    for k from 1 to mn-1 do
      if Verbose then
        printf(
          "LULEM::FFLU(...): processing %d-th row. Length %d:%d\n",
          k, length(convert(M,list)), length(LEM[ListVeil](V))
        );
      end;
      pivot_is_zero, Mkk := DoPivoting( k, M, V, r, c, VeilingStrategy, PivotingStrategy );

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

      top   := apply_veil(Normalizer(numer(Mkk)));
      bot   := apply_veil(Normalizer(denom(Mkk)));
      SS[k] := top;
      # Scaled Shur complement
      tmp        := [k+1..-1];
      M[tmp, k]  := apply_veil~(Normalizer~(M[tmp, k]*bot));
      M[tmp,tmp] := apply_veil~(Normalizer~( top*M[tmp,tmp] - M[tmp,k].M[k,tmp]) );
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
      "M_length" = length(convert(M,list)),
      "S_length" = length(convert(SS,list)),
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

  SolveLinearSystem := proc(
    T::{table},
    b::{Vector},
    V::{symbol, function},
    VeilingStrategy::{procedure} := VeilingStrategy_n,
    $)

    if T["method"] = "LU" then
      SolveLU( T, b, V, VeilingStrategy );
    elif T["method"] = "FFLU" then
      SolveFractionalFreeLU( T, b, V, VeilingStrategy );
    end
  end proc: # SolveLinearSystem

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

    local L, U, r, c, m, n, rnk, apply_veil, x, y, i, j, s:

    # Extract the LU decomposition
    L := T["L"];
    U := T["U"];
    r := T["r"];
    c := T["c"];

    # Get linear system dimension
    m, n := LinearAlgebra[Dimensions](L);

    # Check if the linear system is consistent
    assert(
      LinearAlgebra[RowDimension](L) = LinearAlgebra[ColumnDimension](b),
      "LULEM::SolveLU(...): inconsistent linear system."
    );

    # Create a normalizer function
    apply_veil := (y) -> `if`(VeilingStrategy(y), LEM[Veil][V](y), y);

    # apply permutation P
    x := b[convert(r,list)];

    # Perform forward substitution to solve Ly=b[r]
    for i from 2 to m do
      x[i] := apply_veil( x[i] - add(L[i,1..i-1]*~x[1..i-1]) );
    end do;

    # Perform backward substitution to solve Ux[c]=y
    x[n] := apply_veil(x[n]/U[n,n]);
    for i from n-1 to 1 by -1 do
      s    := apply_veil( x[i] - add(U[i,i+1..n] *~ x[i+1..n]) );
      x[i] := apply_veil( s / U[i,i] );
    end do;

    # apply permutation
    x := x[convert(c,list)];

    # Return outputs
    return x;
  end proc: # SolveLU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SolveFractionalFreeLU := proc(
    T::{table},
    b::{Vector},
    V::{symbol, function},
    VeilingStrategy::{procedure} := VeilingStrategy_n,
    $)

    description "Solve the linear system Ax=b using FFLU decomposition <T>, "
      "provided the vector <b>, the veiling symbol <V> and the veiling strategy "
      "<VeilingStrategy>.";

    local m, n, M, S, r, c, rk, x, i, s,  apply_veil;

    M  := T["M"];
    S  := T["S"];
    r  := T["r"];
    c  := T["c"];
    rk := T["rank"];

    # Get linear system dimension
    m, n := LinearAlgebra[Dimensions](M);

    # Check if the linear system is consistent
    assert( m = n, "LULEM::SolveFractionalFreeLU(...): rectangular linear system." );

    # Create a normalizer function
    apply_veil := (y) -> `if`(VeilingStrategy(y), LEM[Veil][V](y), y);

    # apply permutation P
    x := b[convert(r,list)];

    # Perform forward substitution to solve Ly=b[r]
    for i from 1 to rk-1 do
      x[i+1..-1] := S[i] * x[i+1..-1]; # apply D
      x[i+1..-1] := apply_veil~( x[i+1..-1] - x[i]*M[i+1..-1,i]);
    end do;

    # Perform backward substitution to solve Ux[c]=y
    x[n] := apply_veil(x[n]/M[n,n]);
    for i from n-1 to 1 by -1 do
      s    := apply_veil( x[i] - add(M[i,i+1..n] *~ x[i+1..n]) );
      x[i] := apply_veil( s / M[i,i] );
    end do;

    # apply permutation
    x := x[convert(c,list)];

    # Return outputs
    return x;
  end proc: # SolveFractionalFreeLU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  QR := proc(
    A::{Matrix},
    V::{symbol},
    VeilingStrategy::{procedure}  := VeilingStrategy_n,
    PivotingStrategy::{procedure} := PivotStrategy_Slength,
    $)::{table};

    description "Compute the Householder QR decomposition of a square matrix <A> "
                "using the veiling strategy <VeilingStrategy>, the veiling symbol <V> and "
                "the pivoting strategy <PivotingStrategy>.";

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
      "Q_length" = length(convert(Q,list)),
      "R_length" = length(convert(R,list)),
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
    $)::{table};

    description "Compute the Fracton-Free QR decomposition of a square matrix "
      "<A> using the veiling strategy <VeilingStrategy>, the veiling symbol <V> and "
      "the pivoting strategy <PivotingStrategy>.";

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

  VeilingStrategy_n := proc( x::{algebraic}, $ )::{boolean};
    description "Veiling strategy: number of indeterminates in expression <x> minus 4.";
    return evalb( nops(indets(x)) > 4 and  length(x) > 50);
  end proc: # VeilingStrategy_n

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilingStrategy_L := proc( x::{algebraic}, $ )::{boolean};
    description "Veiling strategy: length of expression <x> minus 50.";
    return evalb(length(x) > 50);
  end proc: # VeilingStrategy_L

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilingStrategy_Ls := proc( x::{algebraic}, $ )::{boolean};
    description "Veiling strategy: length of expression <x> minus 120.";
    return evalb(length(x) > 120);
  end proc: # VeilingStrategy_Ls

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilingStrategy_LB := proc( x::{algebraic}, $ )::{boolean};
    description "Veiling strategy: length of expression <x> minus 260.";
    return evalb(length(x) > 260);
  end proc: # VeilingStrategy_LB

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #   ____  _            _   _
  #  |  _ \(_)_   _____ | |_(_)_ __   __ _
  #  | |_) | \ \ / / _ \| __| | '_ \ / _` |
  #  |  __/| |\ V / (_) | |_| | | | | (_| |
  #  |_|   |_| \_/ \___/ \__|_|_| |_|\__, |
  #                                  |___/

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotingStrategy_Llength := proc( x::{algebraic}, y::{algebraic}, $ )::{boolean};
    description "Pivoting strategy: choose the pivot with the largest length "
                "between expressions <x> and <y>.";
    return evalb( length(SubsData(x)) > length(SubsData(y)) );
  end proc: # PivotingStrategy_Llength

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotingStrategy_Slength := proc( x::{algebraic}, y::{algebraic}, $ )::{boolean};
    description "Pivoting strategy: choose the pivot with the smallest length "
                "between expressions <x> and <y>.";
    return evalb( nops(indets(x)) < nops(indets(y)) );
  end proc: # PivotingStrategy_Slength

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotingStrategy_Lindets := proc ( x::{algebraic}, y::{algebraic}, $ )::{boolean};
    description "Pivoting strategy: choose the pivot with the largest number of "
                "indeterminates between expressions <x> and <y>.";
    return evalb( nops(indets(x)) > nops(indets(y)) );
  end proc: # PivotingStrategy_Lindets

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotingStrategy_Sindets := proc ( x::{algebraic}, y::{algebraic}, $ )::{boolean};
    description "Pivoting strategy: choose the pivot with the smallest number of "
                "indeterminates between expressions <x> and <y>.";
    return evalb( nops(indets(x)) < nops(indets(y)) );
  end proc: # PivotingStrategy_Sindets

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotingStrategy_numeric := proc ( x::{algebraic}, y::{algebraic}, $ )::{boolean};
    description "Pivoting strategy: choose the pivot with the largest numeric "
                "value between expressions <x> and <y>.";
    return evalb( evalf(abs(x)) > evalf(abs(y)) );
  end proc: # PivotingStrategy_numeric

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module:

# That's all folks!
