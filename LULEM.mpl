# # # # # # # # # # # # # # # # # # # # # # # # # # # #
#            _    _   _ _     _____ __  __            #
#           | |  | | | | |   | ____|  \/  |           #
#           | |  | | | | |   |  _| | |\/| |           #
#           | |__| |_| | |___| |___| |  | |           #
#           |_____\___/|_____|_____|_|  |_|           #
#  LU Decomposition with Large Expression Management  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

# LULEM module version 1.1
# BSD 3-Clause License - Copyright (C) 2023
#
# Package created by
# W. Zhou, D. J. Jeffrey, J. Carette and R. M. Corless
#
# Extended version by
# D. Stocco, M. Larcher, E. Bertolazzi
#

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
  export SetModuleOptions,
         SetVerbose,
         Veil,
         UnVeil,
         ShowVeil,
         ListVeil,
         SubsVeil,
         ForgetVeil,
         PermutationMatrices,
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
  local ModuleLoad,
        auxiliary,
        LUPivoting,
        SolvePivotingLU,
        ModuleUnload,
        InitLULEM,
        Protect,
        ##NextLabel,
        ##Signature,
        UnVeilTable,
        lib_base_path,
        verbose;


  # Package options
  option package,
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
      end;
    end;
    if (lib_base_path = null) then
      error "Cannot find 'LULEM' module" ;
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

    unprotect(LastUsed);
    LastUsed := NULL;
    UnVeilTable := NULL;

    printf("Unloading 'LULEM'\n");
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  InitLULEM := proc()

    description "Initialize 'LULEM' module internal variables";

    # Define module variables
    UnVeilTable := table('sparse' = table('sparse' = (0 = 0)));
    verbose     := false;
  end proc: # InitLULEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Protect := proc()

    # Protect global variables
    #protect(
    #);

  end proc: # Protect


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SetVerbose := proc( value::boolean )
    verbose := value;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #   _                    _
  #  | |    ___   ___ __ _| |
  #  | |   / _ \ / __/ _` | |
  #  | |__| (_) | (_| (_| | |
  #  |_____\___/ \___\__,_|_|
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #   _____                       _
  #  | ____|_  ___ __   ___  _ __| |_
  #  |  _| \ \/ / '_ \ / _ \| '__| __|
  #  | |___ >  <| |_) | (_) | |  | |_
  #  |_____/_/\_\ .__/ \___/|_|   \__|
  #             |_|
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Veil := proc( x, $ ) # The expression to be veiled

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

  UnVeil := proc( x, n::{nonnegint, infinity}, $ )
    # x The expression to be Unveiled
    # n The number of levels to UnVeil

    description "UnVeil the expression <x> up to <n> levels.";

    local a, b, level, label;

    label := `if`(procname::indexed, op(procname), '_V');
    level := `if`(nargs<2, 1, min(LastUsed[label]+1, n));
    a     := x;
    b     := eval(a, op(eval(UnVeilTable[label]))[2]); # Always do at least 1 unveiling
    from 2 to level while not Testzero(a - b) do
      a := b;
      b := eval(a, op(eval(UnVeilTable[label]))[2]);
    end do;
    return b;
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auxiliary := proc( x, label, $ )
    # x,    The expression to be veiled
    # label The veiling label to be used

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

  ShowVeil := proc( label::{symbol}, $ )
    # label The veiling label to be shown
    description "Show the veiling variables of the veiling label <label>.";
    print(<ListVeil(label)>);
  end proc: # ShowVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ListVeil := proc( label::{symbol}, $ )
    # label The veiling label to be shown
    description "Return a list of the veiling variables labelled as <label>.";
    return sort(op(eval(UnVeilTable[label]))[2]):
  end proc: # ListVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SubsVeil := proc( label::{symbol}, x, $)
    # label The label of the veiling table to be shown
    # x     The expression to substitute
    description "Substitute the reversed veiling variables of the veiling label "
      "<label> in the expression <x>.";
    return subs(op(ListTools[Reverse](ListVeil(label))), x);
  end proc: # SubsVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ForgetVeil := proc( label::{symbol}, $ )
    # label The veiling label to be forgotten
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

  PermutationMatrices := proc( r::{Vector}, c::{Vector}, $ )
    # r row pivot vector
    # c column pivot vector
    description "Compute the LU decomposition premutation matrix provided the "
                "NAG-style the pivot vector <r>.";

    local nr, nc, P, Q, i;

    nr := LinearAlgebra[RowDimension](r):
    nc := LinearAlgebra[RowDimension](c):
    P  := Matrix(nr,nr);
    Q  := Matrix(nc,nc);
    for i from 1 to nr by 1 do
      P[i, r[i]] := 1;
    end do;
    for i from 1 to nc by 1 do
      Q[c[i],i] := 1;
    end do;
    return P, Q;
  end proc: # PermutationMatrices

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  LUPivoting := proc(
    k::integer,                    # k-th step
    M::{Matrix},                   # Linear system matrix A
    Q::{symbol},                   # Symbol for the veiling
    r::{Vector(integer)},          # Symbol for the veiling
    c::{Vector(integer)},          # Symbol for the veiling
    Strategy_Veiling::{procedure}, # Veiling strategy
    Strategy_Pivots::{procedure},  # Pivoting strategy procedure
    Strategy_Zero::{procedure},    # Zero recognition strategy procedure
    $)

    description "Compute the LU decomposition of a square matrix <A> using the "
                "LEM strategy <Strategy_Veiling>, the veiling symbol <Q>, the pivoting "
                "strategy <Strategy_Pivots> and the zero recognition strategy <Strategy_Zero>.";

    local Mij, Mkk, n, m, i, j, check_veil, pivot_is_zero, Mij_is_zero, z, tmp;

    n, m := LinearAlgebra[Dimensions](M):

    # check if Veil or not
    check_veil := z -> `if`(Strategy_Veiling(z) > 0, Veil[Q](z), z);

    # check if M[r[k],c[k]] = 0, if not true it is the pivot
    Mkk := M[k,k];
    try
      pivot_is_zero := evalb(Strategy_Zero(SubsVeil(Q,Mkk)) = 0);
      #pivot_is_zero := evalb(Strategy_Zero(Normalizer(Mkk)) = 0);
    catch:
      print("divide by 0 or numerical exception\n");
      print(Mkk);
      pivot_is_zero := true;
    end try;
    # search for a nonzero pivot
    for j from k to m do
      for i from k to n do
        Mij := M[i,j];
        try
          Mij_is_zero := evalb(Strategy_Zero(SubsVeil(Q,Mij)) = 0);
          #Mij_is_zero := evalb(Strategy_Zero(Normalizer(Mij)) = 0);
        catch:
          print("divide by 0 or numerical exception\n");
          print(Mij);
          Mij_is_zero := true;
        end try;
        if not Mij_is_zero and (pivot_is_zero or Strategy_Pivots( Mij, Mkk )) then
          # found better pivot
          pivot_is_zero := false;
          if i <> k then
            (r[i], r[k])   := (r[k], r[i]);
            M[[i,k],1..-1] := M[[k,i],1..-1];
          end if;
          if j <> k then
            (c[j], c[k])   := (c[k], c[j]);
            M[1..-1,[j,k]] := M[1..-1,[k,j]];
          end if;
          Mkk := Mij;
        end if;
      end do;
    end do;
    return pivot_is_zero, Mkk
  end proc: # LUD

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

    local LL, UU, M, Mkk, n, m, nm, i, j, k, ri, rk, rnk, r, c, check_veil, pivot_is_zero, Mij_is_zero, z, tmp;

    n, m := LinearAlgebra[Dimensions](A):

    # Create pivot vector
    r := Vector(n, k -> k);
    c := Vector(m, k -> k);

    # check if Veil or not
    check_veil := z -> `if`(Strategy_Veiling(z) > 0, Veil[Q](z), z);

    # Gauss Elimination main loop
    M   := copy(A);  # make a working copy
    nm  := min(n,m);
    rnk := nm;
    for k from 1 to nm-1 do
      if verbose then
        printf("Process row N.%d\n",k);
      end;
      pivot_is_zero, Mkk := LUPivoting( k, M, Q, r, c, Strategy_Veiling, Strategy_Pivots, Strategy_Zero );

      if pivot_is_zero then
        rnk := k;
        if verbose then
          WARNING( "LULEM::LUD(...): the matrix appears not full rank." );
        end;
        break;
      end if;

      if verbose then
        print("PIVOT:",Mkk);
      end;

      # Shur complement
      tmp        := [k+1..-1];
      M[tmp,k]   := M[tmp,k]/Mkk;
      M[tmp,tmp] := check_veil~(Normalizer~(M[tmp,tmp]-M[tmp,k].M[k,tmp]));
    end do:

    LL := Matrix(M[1..n,1..n], shape = triangular[lower,unit]);
    UU := Matrix(M, shape = triangular[upper]);

    # Return the LU decomposition and the pivot vector
    return LL, UU, r, c, rnk;
  end proc: # LUD

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FFLUD := proc(
    A::{Matrix},                                            # Linear system matrix A
    Q::{symbol},                                            # Symbol for the veiling
    Strategy_Veiling::{procedure} := VeilingStrategy_n,     # Veiling strategy
    Strategy_Pivots::{procedure}  := PivotStrategy_Slength, # Pivoting strategy procedure
    Strategy_Zero::{procedure}    := ZeroStrategy_length,   # Zero recognition strategy procedure
    $)

    description "Compute the LU decomposition of a square matrix <A> using the "
                "LEM strategy <Strategy_Veiling>, the veiling symbol <Q>, the pivoting "
                "strategy <Strategy_Pivots> and the zero recognition strategy <Strategy_Zero>.";

    local SS, M, Mkk, n, m, nm, i, j, k, ri, rk, rnk, r, c, check_veil, pivot_is_zero, Mij_is_zero, z, tmp;

    n, m := LinearAlgebra[Dimensions](A):

    # Create pivot vector
    r := Vector(n, k -> k);
    c := Vector(m, k -> k);

    # check if Veil or not
    check_veil := z -> `if`(Strategy_Veiling(z) > 0, Veil[Q](z), z);

    # Gauss Elimination main loop
    M   := copy(A);  # make a working copy
    nm  := min(n,m);
    rnk := nm;
    SS  := Vector(nm);
    for k from 1 to nm-1 do
      if verbose then
        printf("Process row N.%d\n",k);
      end;
      pivot_is_zero, Mkk := LUPivoting( k, M, Q, r, c, Strategy_Veiling, Strategy_Pivots, Strategy_Zero );

      if pivot_is_zero then
        rnk := k;
        if verbose then
          WARNING( "LULEM::LUD(...): the matrix appears not full rank." );
        end;
        break;
      end if;

      if verbose then
        print("PIVOT:",Mkk);
      end;

      SS[k] := Mkk;
      # Scaled Shur complement
      tmp        := [k+1..-1];
      M[tmp,tmp] := check_veil~(Normalizer~(Mkk*M[tmp,tmp]-M[tmp,k].M[k,tmp]));
    end do:

    # Return the LU decomposition and the pivot vector
    return M, SS, r, c, rnk;
  end proc: # FFLUD

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Solve the linear system Ax=b given the LU decomposition (PA=LU) in NAG-style.
  # NOTE: The pivot vector r is returned from LUD function.
  SolvePivotingLU := proc(
    LU_NAG::{Matrix},              # Compact LU matrix (NAG-style form)
    r::{Vector},                   # row exchange
    c::{Vector},                   # column exchange
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
    x[c[n]] := normalizer(y[n]/LU_NAG[n, n]);
    for i from n-1 to 1 by -1 do
      s := normalizer(y[i]) - add(normalizer(LU_NAG[i, j] * x[c[j]]), j = i+1..n);
      x[c[i]] := normalizer(s/LU_NAG[i, i]);
    end do;

    # Return solution vector x
    return x;
  end proc: # SolvePivotingLU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Solve := proc(
    A::{Matrix},                   # Linear system matrix A
    b::{Vector},                   # Linear system vector b
    Q::{symbol, function},         # Symbol for the veiling
    Strategy_Veiling::{procedure} := VeilingStrategy_n,    # Veiling strategy
    Strategy_Pivots::{procedure}  := PivotStrategy_Slength, # Pivoting strategy procedure
    Strategy_Zero::{procedure}    := ZeroStrategy_length,     # Zero recognition strategy procedure
    $)

    description "Solve the linear system Ax=b using LULEM algorithm, privided the "
      "matrix <A>, the vector <b>, the veiling symbol <Q>, the veiling strategy "
      "<Strategy_Veiling>, the pivoting strategy <Strategy_Pivots> and the zero "
      "recognition strategy <Strategy_Zero>.";

    local x, P, L, U, r, LU_NAG:

    # Get LU decomposition of A
    P, L, U, r := LUD(A, Q, Strategy_Veiling, Strategy_Pivots, Strategy_Zero);

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

    return evalb(length(x) - length(y) > 0);
  end proc: # PivotStrategy_Llength

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PivotStrategy_Slength := proc(
    x::{algebraic}, # Expression to be analyzed
    y::{algebraic}, # Expression to be analyzed
    $)::{boolean};

    description "Pivoting strategy: choose the pivot with the smallest length "
      "between expressions <x> and <y>.";

    return evalb(length(x) < length(y));
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

    return length(x);
  end proc: # ZeroStrategy_length

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ZeroStrategy_normalizer := proc(
    x::{algebraic}, # Expression to be analyzed
    $) # FIXME: What is the return type of this procedure?

    description "Zero recognition strategy: normalizer of expression <x>.";

    return Normalizer(x);
  end proc: # ZeroStrategy_normalizer

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module:

# That's all folks!
