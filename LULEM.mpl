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

LULEM := module()

  export  AssignData,
          SubsData,
          ForgetData,
          SetVerbose,
          PermutationMatrices,
          LUD,
          SolveLUD,
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
          verbose;

  uses LEM;

  option  package,
          load   = ModuleLoad,
          unload = ModuleUnload;

  description "LU and QR decomposition with LEM (Large Expressions Management).";

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
    verbose    := false;
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
  end proc: # PermutationMatrices

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SetVerbose := proc( value::boolean )
    verbose := value;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  LUD := proc(
    A::{Matrix},
    Q::{symbol},
    VeilingStrategy::{procedure}  := VeilingStrategy_n,
    PivotingStrategy::{procedure} := PivotingStrategy_Slength,
    ZeroStrategy::{procedure}     := ZeroStrategy_length,
    $)

    description "Compute the LU decomposition of a square matrix <A> using the "
      "LEM strategy <VeilingStrategy>, the veiling symbol <Q>, the pivoting "
      "strategy <PivotingStrategy> and the zero recognition strategy "
      "<ZeroStrategy>.";

     local LL, DD, UU, M, Mij, Mkk, n, m, nm, i, j, k, rnk, mltip, r, c, temp, normalize,
          pivot_is_zero, Mij_is_zero, z;

    # Copy the input matrix
    M := copy(A);
    n := LinearAlgebra[RowDimension](M):
    m := LinearAlgebra[ColumnDimension](M);

    # Create pivot vector
    r := Vector(n, k -> k);
    c := Vector(m, k -> k);

    # Check if Veil or not
    normalize := z -> `if`(VeilingStrategy(z) > 0, Veil[Q](z), z);

    # Gauss elimination main loop
    nm  := min(n,m);
    rnk := nm;
    for k from 1 to nm do
      if verbose then
        printf("Process row N.%d\n",k);
      end;
      # check if M[r[k],c[k]] = 0, if not true it is the pivot
      Mkk := M[r[k],c[k]];
      try
        pivot_is_zero := evalb(ZeroStrategy(SubsVeil(Mkk, Q)) = 0);
        #pivot_is_zero := evalb(ZeroStrategy(Normalizer(Mkk)) = 0);
      catch:
        print("divide by 0 or numerical exception\n");
        print(Mkk);
        pivot_is_zero := true;
      end try;
      # search for a nonzero pivot
      for j from k to m do
        for i from k to n do
          Mij := M[r[i],c[j]];
          try
            Mij_is_zero := evalb(ZeroStrategy(SubsVeil(Mij, Q)) = 0);
            #Mij_is_zero := evalb(ZeroStrategy(Normalizer(Mij)) = 0);
          catch:
            print("divide by 0 or numerical exception\n");
            print(Mij);
            Mij_is_zero := true;
          end try;
          if not Mij_is_zero then
            if pivot_is_zero then
              # found better pivot
              pivot_is_zero := false;
              (r[i], r[k]) := (r[k], r[i]);
              (c[j], c[k]) := (c[k], c[j]);
            elif PivotingStrategy(Mij, Mkk) then
              (r[i], r[k]) := (r[k], r[i]);
              (c[j], c[k]) := (c[k], c[j]);
              # break; # Once a pivot is found -- not "best"!
            end if;
          end if;
        end do;
      end do;

      # Only when the pivot is not equal to zero, the row elimination is performed
      # (the whole if statement). Else we will continue with the next column.

      if pivot_is_zero then
        rnk := k;
        if verbose then
          WARNING( "LULEM::LUD(...): the matrix appears not full rank." );
        end;
        break;
      end if;

      # Now do Gauss elimination steps to get new M and also keep L information
      # in new M. Packing everything into the M matrix during the computation of
      # the factors. The reason is that this code can be ported to restricted
      # memory environments, or equally, applied to very large matrices in a
      # large memory environment.

      Mkk := M[r[k],c[k]];
      if verbose then
        print("PIVOT:",Mkk);
      end;

      for i from k+1 to n do
        # per azzerare M[r[i],...]
        # M[r[i],...] = M[r[i],...] - mltip*M[r[k],...]
        mltip        := Normalizer(M[r[i],c[k]]/Mkk);
        M[r[i],c[k]] := mltip;
        for j from k+1 to m do
          z := Normalizer(M[r[i],c[j]] - mltip * M[r[k],c[j]]); # si puo mettere j = c[j]?
          M[r[i],c[j]] := normalize(z);
        end do;
      end do;
    end do;

    # L is lower matrix with l's as diagonal entries
    LL := Matrix(LinearAlgebra[IdentityMatrix](n), shape = triangular[lower]);

    # D is lower matrix with l's as diagonal entries
    DD := Matrix(n, n, shape = diagonal);

    # U is upper matrix.
    UU := Matrix(n, m, shape = triangular[upper]);

    # Seperate new M into L and U
    for i from 1 to n do
      for j from 1 to m do
        z := M[r[i],c[j]];
        if (i <= j) then
          UU[i,j] := z;
        else
          LL[i,j] := z;
        end if
      end do;
    end do;

    # Return the LU decomposition and the pivot vector
    return LL, UU, r, c, rnk;

  end proc: # LUD

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SolveLUD := proc(
    A::{Matrix},
    b::{Vector},
    Q::{symbol, function},
    VeilingStrategy::{procedure}  := VeilingStrategy_n,
    PivotingStrategy::{procedure} := PivotingStrategy_Slength,
    ZeroStrategy::{procedure}     := ZeroStrategy_length,
    $)

    description "Solve the linear system Ax=b using LULEM algorithm, privided the "
      "matrix <A>, the vector <b>, the veiling symbol <Q>, the veiling strategy "
      "<VeilingStrategy>, the pivoting strategy <PivotingStrategy> and the zero "
      "recognition strategy <ZeroStrategy>.";

    local n, m, LL, UU, r, c, rnk, LU_NAG, normalizer, y, i, x, s:

    # Get linear system dimension
    n := LinearAlgebra[RowDimension](A);
    m := LinearAlgebra[ColumnDimension](A);

    # Check if matrix A is square
    assert(
      n = LinearAlgebra[ColumnDimension](b),
      "LULEM::SolveLUD(...): input matrix A is not square."
    );

    # Check if the linear system is consistent
    assert(
      n = m,
      "LULEM::SolveLUD(...): the linear system Ax=b is not consistent."
    );

    # Get LUD decomposition of A
    LL, UU, r, c, rnk := LUD(A, Q, VeilingStrategy, PivotingStrategy, ZeroStrategy);

    # Built the LU matrix (NAG-style)
    LU_NAG := LL + UU - Matrix(n, n, shape = identity);

    # Create vector for solution of Ly=Pb
    y := Vector(n);

    # Create vector for solution of Ux=y
    x := Vector(n);

    # Create a normalizer function
    normalizer := (y) -> `if`(VeilingStrategy(y) > 0, Veil[Q](y), y);

    # Perform forward substitution to solve Ly=Pb
    y[1] := normalizer(b[r[1]]);
    for i from 2 to n do
      y[i] := normalizer(b[r[i]]) - add(normalizer(LU_NAG[i, j] * y[j]), j = 1..i-1);
    end do;

    # Perform backward substitution to solve Ux=y
    x[c[n]] := normalizer(y[n]/LU_NAG[n, n]);
    for i from n-1 to 1 by -1 do
      s := normalizer(y[i]) - add(normalizer(LU_NAG[i, j] * x[c[j]]), j = i+1..n);
      x[c[i]] := normalizer(s/LU_NAG[i, i]);
    end do;

    # Return solution vector x and the LU decomposition data
    if (_nresults = 6) then
      return x, L, U, r, c, rnk;
    elif (_nresults = 2) then
      return x, rnk;
    else
      return x;
    end if;
  end proc: # SolveLUD

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