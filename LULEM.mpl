# # # # # # # # # # # # # # # # # # # # # # # # # # # #
#            _    _   _ _     _____ __  __            #
#           | |  | | | | |   | ____|  \/  |           #
#           | |  | | | | |   |  _| | |\/| |           #
#           | |__| |_| | |___| |___| |  | |           #
#           |_____\___/|_____|_____|_|  |_|           #
#  LU Decomposition with Large Expression Management  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Authors of the current version:
#  Davide Stocco (University of Trento)
#  Matteo Larcher (University of Trento)
#  Enrico Bertolazzi (University of Trento)
#
# Author of the original code:
#   Wenqin Zhou (University of Western Ontario) - Former affiliation
#   David J. Jeffrey (University of Western Ontario)
#   Jacques Carette (McMaster University)
#   Robert M. Corless (University of Western Ontario)
#
#
# This is a module for the 'LULEM' (LU Decomposition Large Expressions Management)
# package. It contains the functions to solve systems of linear equations with
# LEM (Large Wxpressions Management). The module uses the pivoting LU decomposition
# to solve the system.
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

LULEM := module()

  export  AssignData,
          SubsData,
          ForgetData,
          SetVerbosity,
          PermutationMatrices,
          LU,
          LUPivoting,
          SolveLU,
          FFLU,
          SolveFFLU,
          #TODO: QR,
          #TODO: SolveQR,
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

  uses LEM;

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
  end proc: # InitLULEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  AssignData := proc(
    x::{list, set}, # The list to be assigned
    $)::{nothing};

    description "Assign the data list <x> to the local variable <StoredData>.";

    StoredData := x;
  end proc: # AssignData

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SubsData := proc(
    x::{anything},
    $)::{anything};

    description "Substitute the local variable <StoredData> in the expression <x>.";

    local label, out, i;

    label := lhs~(op(op(LastUsed))[2..-1]);
    out := copy(x);
    if (nops(StoredData) > 0) then
      for i in label do
        out := subs[eval](StoredData,
          subs[eval](op(ListTools[Reverse](ListVeil(i))), out)
        );
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

  PermutationMatrices := proc(
    r::{Vector},
    c::{Vector},
    $)

    description "Compute the LU decomposition premutation matrix provided the "
                "rows the pivot vector <r> and the columns the pivot vector <c>.";

    local n, m, i, P, Q;

    n := LinearAlgebra[RowDimension](r):
    m := LinearAlgebra[RowDimension](c):
    P := Matrix(n, n);
    Q := Matrix(m, m);
    for i from 1 to n by 1 do
      P[i, r[i]] := 1;
    end do;
    for i from 1 to m by 1 do
      Q[c[i], i] := 1;
    end do;
    return P, Q;
  end proc: # PermutationMatrices

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  LUPivoting := proc(
    k::{integer},
    M::{Matrix},
    Q::{symbol},
    r::{Vector(integer)},
    c::{Vector(integer)},
    VeilingStrategy::{procedure},
    PivotingStrategy::{procedure},
    ZeroStrategy::{procedure},
    $)

    description "Compute the LU decomposition pivot vector provided the step <k>, "
      "the temporary LU (NAG) matrix <M>, the veiling symbol <Q>, the rows the "
      "pivot vector <r>, the columns the pivot vector <c>, the veiling strategy "
      "<VeilingStrategy>, the pivoting strategy <PivotingStrategy> and the zero "
      "recognition strategy <ZeroStrategy>.";

    local Mij, Mkk, n, m, i, j, apply_veil, pivot_is_zero, Mij_is_zero, z, tmp;

    # Extract dimensions
    n, m := LinearAlgebra[Dimensions](M):

    # Check if to veil or not
    apply_veil := (z) -> `if`(VeilingStrategy(z) > 0, Veil[Q](z), z);

    # Check if M[r[k],c[k]] = 0, if not true it is the pivot
    Mkk := M[k,k];
    try
      pivot_is_zero := evalb(ZeroStrategy(SubsVeil(Mkk, Q)) = 0);
      #pivot_is_zero := evalb(ZeroStrategy(Normalizer(Mkk)) = 0);
    catch:
      print("Mkk: Division by 0 or numerical exception.\n");
      print(Mkk);
      pivot_is_zero := true;
    end try;

    # Iterate over the columns and rows
    for j from k to m do
      for i from k to n do

        # Look for a non-zero pivot
        Mij := M[i,j];
        try
          Mij_is_zero := evalb(ZeroStrategy(SubsVeil(Mij, Q)) = 0);
          #Mij_is_zero := evalb(ZeroStrategy(Normalizer(Mij)) = 0);
        catch:
          print("Mij: Division by 0 or numerical exception.\n");
          print(Mij);
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
    Q::{symbol},
    VeilingStrategy::{procedure}  := VeilingStrategy_n,
    PivotingStrategy::{procedure} := PivotingStrategy_Slength,
    ZeroStrategy::{procedure}     := ZeroStrategy_length,
    $)

    description "Compute the LU decomposition of a square matrix <A> using the "
      "veiling strategy <VeilingStrategy>, the veiling symbol <Q>, the pivoting "
      "strategy <PivotingStrategy> and the zero recognition strategy "
      "<ZeroStrategy>.";

    local M, L, U, Mkk, n, m, nm, k, rnk, r, c, apply_veil, pivot_is_zero, Mij_is_zero, tmp;

    n, m := LinearAlgebra[Dimensions](A):

    # Create pivot vector
    r := Vector(n, k -> k);
    c := Vector(m, k -> k);

    # Check if to veil or not
    apply_veil := (z) -> `if`(VeilingStrategy(z) > 0, Veil[Q](z), z);

    # Perform Gaussian elimination
    M   := copy(A);
    nm  := min(n, m);
    rnk := nm;
    for k from 1 to (nm - 1) do
      if Verbose then
        printf("LULEM::LU(...): processing %d-th row.\n", k);
      end;
      pivot_is_zero, Mkk := LUPivoting(
        k, M, Q, r, c, VeilingStrategy, PivotingStrategy, ZeroStrategy
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

    L := Matrix(M[1..n, 1..n], shape = triangular[lower, unit]);
    U := Matrix(M, shape = triangular[upper]);

    # Return the LU decomposition and the pivot vector
    if (_nresults = 5) then
      return L, U, r, c, rnk;
    else
      return M, r, c, rnk;
    end if;
  end proc: # LU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SolveLU := proc(
    A::{Matrix},
    b::{Vector},
    Q::{symbol, function},
    VeilingStrategy::{procedure}  := VeilingStrategy_n,
    PivotingStrategy::{procedure} := PivotingStrategy_Slength,
    ZeroStrategy::{procedure}     := ZeroStrategy_length,
    $)

    description "Solve the linear system Ax=b using LU decomposition, provided "
      "the matrix <A>, the vector <b>, the veiling symbol <Q>, the veiling "
      "strategy <VeilingStrategy>, the pivoting strategy <PivotingStrategy> and "
      "the zero recognition strategy <ZeroStrategy>.";

    local n, m, M, r, c, rnk, LU_NAG, apply_veil, y, i, x, s:

    # Get linear system dimension
    n := LinearAlgebra[RowDimension](A);
    m := LinearAlgebra[ColumnDimension](A);

    # Check if matrix A is square
    assert(
      n = LinearAlgebra[ColumnDimension](b),
      "LULEM::SolveLU(...): input matrix A is not square."
    );

    # Check if the linear system is consistent
    assert(
      n = m,
      "LULEM::SolveLU(...): the linear system Ax=b is not consistent."
    );

    # Get LU decomposition of A (NAG-style)
    M, r, c, rnk := LU(A, Q, VeilingStrategy, PivotingStrategy, ZeroStrategy);

    # Create vector for solution of Ly=Pb
    y := Vector(n);

    # Create vector for solution of Ux=y
    x := Vector(n);

    # Create a normalizer function
    apply_veil := (y) -> `if`(VeilingStrategy(y) > 0, Veil[Q](y), y);

    # Perform forward substitution to solve Ly=Pb
    y[1] := apply_veil(b[r[1]]);
    for i from 2 to n do
      y[i] := apply_veil(b[r[i]]) - add(apply_veil(M[i, j] * y[j]), j = 1..i-1);
    end do;

    # Perform backward substitution to solve Ux=y
    x[c[n]] := apply_veil(y[n] / M[n, n]);
    for i from (n - 1) to 1 by -1 do
      s := apply_veil(y[i]) - add(apply_veil(M[i, j] * x[c[j]]), j = i+1..n);
      x[c[i]] := apply_veil(s / M[i, i]);
    end do;

    # Return outputs
    if (_nresults = 2) then
      return x, rnk;
    else
      return x;
    end if;
  end proc: # SolveLUD

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FFLU := proc(
    A::{Matrix},                                            # Linear system matrix A
    Q::{symbol},                                            # Symbol for the veiling
    VeilingStrategy::{procedure} := VeilingStrategy_n,     # Veiling strategy
    PivotingStrategy::{procedure}  := PivotStrategy_Slength, # Pivoting strategy procedure
    ZeroStrategy::{procedure}    := ZeroStrategy_length,   # Zero recognition strategy procedure
    $)

    description "Compute the LU decomposition of a square matrix <A> using the "
                "LEM strategy <VeilingStrategy>, the veiling symbol <Q>, the pivoting "
                "strategy <PivotingStrategy> and the zero recognition strategy <ZeroStrategy>.";

    local SS, M, Mkk, n, m, nm, i, j, k, ri, rk, rnk, r, c, apply_veil, pivot_is_zero, Mij_is_zero, z, tmp;

    n, m := LinearAlgebra[Dimensions](A):

    # Create pivot vector
    r := Vector(n, k -> k);
    c := Vector(m, k -> k);

    # check if Veil or not
    apply_veil := z -> `if`(VeilingStrategy(z) > 0, Veil[Q](z), z);

    # Gauss Elimination main loop
    M   := copy(A);  # make a working copy
    nm  := min(n,m);
    rnk := nm;
    SS  := Vector(nm);
    for k from 1 to nm-1 do
      if Verbose then
        printf("Process row N.%d\n",k);
      end;
      pivot_is_zero, Mkk := LUPivoting( k, M, Q, r, c, VeilingStrategy, PivotingStrategy, ZeroStrategy );

      if pivot_is_zero then
        rnk := k;
        if Verbose then
          WARNING( "LULEM::LUD(...): the matrix appears not full rank." );
        end;
        break;
      end if;

      if Verbose then
        print("PIVOT:",Mkk);
      end;

      SS[k] := Mkk;
      # Scaled Shur complement
      tmp        := [k+1..-1];
      M[tmp,tmp] := apply_veil~(Normalizer~(Mkk*M[tmp,tmp]-M[tmp,k].M[k,tmp]));
    end do:

    # Return the LU decomposition and the pivot vector
    return M, SS, r, c, rnk;
  end proc: # FFLU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

    #local tmp;
    #tmp := SubsData(x);
    #if evalb((evalf(abs(tmp)) = 0.0)) then
    #  return length(0);
    #else
    #  return length(tmp);
    #end if;
    return length(x);
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