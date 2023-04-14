# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   _____ _____ _    _   _
#  |  ___|  ___| |  | | | |
#  | |_  | |_  | |  | | | |
#  |  _| |  _| | |__| |_| |
#  |_|   |_|   |_____\___/
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export FFLU::static := proc(
  _self::LAST,
  A::Matrix,
  $)::table;

  description "Compute the Fracton-Free LU decomposition of a square matrix <A>.";

  local V, SS, M, pivot, pivot_list, m, n, mn, i, j, k, rk, rnk, r, c, tmp, bot,
    top;

    # Check if the LEM is initialized
  if not type(_self:-m_LEM, LEM) then
    error "LEM is not initialized (use LAST::InitLEM() first).";
    return table([]);
  end if;

  # Get the veiling label
  V := _self:-m_LEM:-GetVeilingLabel(_self:-m_LEM);

  # Sanity check
  if has(A, V) then
    error "veiling symbol %1 is already present in matrix coefficient.", V;
    return table([]);
  end if;

  # Forget the veilings
  _self:-m_LEM:-ForgetVeil(_self:-m_LEM);

  # Get matrix dimensions
  m, n := LinearAlgebra:-Dimensions(A):

  # Create pivot vector
  r := Vector(m, k -> k);
  c := Vector(n, k -> k);

  # Perform Gaussian elimination
  M          := copy(A);
  mn         := min(m, n);
  rnk        := mn;
  SS         := Vector(mn);
  pivot_list := [];
  for k from 1 to mn-1 do
    if _self:-m_VerboseMode then
      printf(
        "LAST::FFLU(...): processing %d-th row, cost = %d, veilings = %d.\n",
        k, _self:-m_LEM:-ExpressionCost(_self:-m_LEM, M),
        nops(_self:-m_LEM:-VeilList(_self:-m_LEM))
      );
    end;

    pivot := _self:-Pivoting(_self, k, M, r, c);
    if not pivot["is_zero"] then
      pivot_list := [op(pivot_list), pivot["value"]];
    end if;

    if pivot["is_zero"] then
      rnk := rnk - 1;
      if _self:-m_WarningMode then
        WARNING("LAST::LU(...): the matrix appears to be not full rank.");
      end;
      break;
    end if;

    if _self:-m_VerboseMode then
      printf(
        "LAST::FFLU(...): M[%d,%d] = %a, cost = %d, degree_r = %d, degree_c = %d.\n",
        k, k, pivot["value"], pivot["cost"], pivot["degree_r"], pivot["degree_c"]
      );
    end if;

    top   := _self:-m_LEM:-Veil(_self:-m_LEM, Normalizer(numer(pivot["value"])));
    bot   := _self:-m_LEM:-Veil(_self:-m_LEM, Normalizer(denom(pivot["value"])));
    SS[k] := top;

    # Scaled Shur complement
    tmp         := [k+1..-1];
    M[tmp, k]   := _self:-m_LEM:-Veil~(
      _self:-m_LEM, Normalizer~(M[tmp, k]*bot)
    );
    M[tmp, tmp] := _self:-m_LEM:-Veil~(
      _self:-m_LEM, Normalizer~(top*M[tmp,tmp] - M[tmp,k].M[k,tmp])
    );
  end do;

  # Return the FFLU decomposition
  return table([
    "method" = "FFLU",
    "M"      = M,
    "V"      = V,
    "S"      = SS,
    "r"      = r,
    "c"      = c,
    "rank"   = rnk,
    "pivots" = pivot_list,
    "M_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, M),
    "S_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, SS),
    "V_cost" = _self:-m_LEM:-ExpressionCost(
      _self:-m_LEM, _self:-m_LEM:-VeilList(_self:-m_LEM)
    ),
    "M_nnz"  = nops(op(2, M)),
    "A_nnz"  = nops(op(2, A))
  ]);
end proc: # FFLU

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export FF2LU::static := proc(
  _self::LAST,
  T::table,
  $)::list, list, Matrix;

  description "Compute the LU decomposition of a square matrix from its "
    "Fracton-Free LU decomposition <T>.";

  local M, SS, r, c, rk, n, m, L, DG, L_list, D_list, i, j, k;

  M      := T["M"];
  SS     := T["S"];
  r      := T["r"];
  c      := T["c"];
  rk     := T["rank"];
  m, n   := LinearAlgebra:-Dimensions(M):
  L_list := [];
  D_list := [];

  # +             + +         +      +            +      +             + +            +
  # |  1          | | 1       |      | 1          |      | 1           | | 1          |
  # | -x  1       | |   q     |  =>  | -x q       |  =>  |   1/q       | | x  1       |
  # | -x     1    | |     q   |      | -x    q    |      |      1/q    | | x     1    |
  # | -x        1 | |       q |      | -x       q |      |         1/q | | x        1 |
  # +             + +         +      +            +      +             + +            +

  for k from 1 to rk do
    L            := Matrix(m, m, shape = triangular[lower, unit]);
    L[k+1..-1,k] := M[k+1..-1, k];
    L_list       := [op(L_list), L];
    DG           := Matrix(LinearAlgebra:-IdentityMatrix(m), shape = diagonal);
    for j from k+1 to m do
      DG[j,j] := SS[k];
    end;
    D_list := [op(D_list), DG];
  end;
  return L_list, D_list, Matrix(M, shape = triangular[upper]);
end proc: # FF2LU

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export FFLUsolve::static := proc(
  _self::LAST,
  T::table,
  b::Vector,
  $)

  description "Solve the linear system Ax=b using FFLU decomposition <T> and "
    "provided the vector <b>.";

  local m, n, M, S, r, c, rnk, x, y, i, s;

  # Check if the LEM is initialized
  if not type(_self:-m_LEM, LEM) then
    error "LEM is not initialized (use LAST::InitLEM() first).";
    return NULL;
  end if;

  M   := T["M"];
  S   := T["S"];
  r   := T["r"];
  c   := T["c"];
  rnk := T["rank"];

  # Get linear system dimension
  m, n := LinearAlgebra:-Dimensions(M);

  # Check if the linear system is consistent
  if not m = n then
    error "only square system can be solved (got M = %1 x %2).", m, n;
    return NULL;
  end if;

  # Check if the linear system is consistent
  if not n = rnk then
    error "only full rank linear system can be solved (got rank = %1, expected "
      "rank = %2).", rnk, n;
    return NULL;
  end if;

  # apply permutation P
  x := b[convert(r, list)];

  # Perform forward substitution to solve Ly=b[r]
  for i from 2 to n do
    if _self:-m_VerboseMode then
      printf("LAST::FFLUsolve(...): forward substitution of %d-th row.\n", i);
    end if;
    x[i..-1] := S[i-1] * x[i..-1]; # Apply D
    x[i..-1] := _self:-m_LEM:-Veil~(_self:-m_LEM, x[i..-1] - x[i-1]*M[i..-1, i-1]);
  end do;

  # Perform backward substitution to solve Ux[c]=y
  if _self:-m_VerboseMode then
    printf("LAST::FFLUsolve(...): dividision by M[%d,%d].\n", n, n);
  end if;
  x[n] := _self:-m_LEM:-Veil(_self:-m_LEM, x[n]/M[n, n]);
  for i from n-1 to 1 by -1 do
    if _self:-m_VerboseMode then
      printf("LAST::FFLUsolve(...): backward substitution of %d-th column.\n", i);
    end if;
    s    := x[i] - add(M[i, i+1..n] *~ x[i+1..n]);
    x[i] := _self:-m_LEM:-Veil(_self:-m_LEM, s/M[i, i]);
  end do;

  # Apply inverse permutation Q
  y := Vector[column](n);
  for i from 1 to n do
    y[c[i]] := x[i];
  end do;

  # Return outputs
  return y;
end proc: # FFLUsolve

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
