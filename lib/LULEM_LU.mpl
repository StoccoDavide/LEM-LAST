# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   _    _   _
#  | |  | | | |
#  | |  | | | |
#  | |__| |_| |
#  |_____\___/
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export LU::static := proc(
  _self::LULEM,
  A::Matrix,
  $)::table;

  description "Compute the LU decomposition of a square matrix <A>.";

  local V, M, L, U, pivot, pivot_list, m, n, mn, k, rnk, r, c, apply_veil, tmp;

  print("here1", _self:-m_LEM);

  # Get the veiling label
  _self:-m_LEM:-Info(_self:-m_LEM);
  V := _self:-m_LEM:-GetVeilingLabel(_self:-m_LEM);

  print("hereV", V);

  # Sanity check
  if has(A, V) then
    error "veiling symbol %1 is already present in matrix coefficient.", V;
    return table([]);
  end if;

  print("here2");

  # Forget the veilings
  _self:-m_LEM:-VeilForget(_self:-m_LEM);
  print("here3");

  # Get matrix dimensions
  m, n := LinearAlgebra:-Dimensions(A):

  # Create pivot vector
  r := Vector(m, k -> k);
  c := Vector(n, k -> k);

  # Check if to veil or not
  apply_veil := (z) -> _self:-m_LEM:-Veil(_self:-m_LEM, z);

  # Perform Gaussian elimination
  M          := copy(A);
  mn         := min(m, n);
  rnk        := mn;
  pivot_list := [];
  for k from 1 to mn do
    if _self:-m_Verbose then
      printf(
        "LULEM::LU(...): processing %d-th row, cost = %d, veilings = %d.\n",
        k, _self:-m_LEM:-ExpressionCost(_self:-m_LEM, M),
        nops(_self:-m_LEM:-VeilList(_self:-m_LEM))
      );
    end if;

    pivot := _self:-Pivoting(_self, k, M, r, c);
    if not pivot["is_zero"] then
      pivot_list := [op(pivot_list), pivot["value"]];
    end if;

    if pivot["is_zero"] then
      rnk := rnk - 1;
      if _self:-m_Verbose then
        WARNING("LULEM::LU(...): the matrix appears to be not full rank.");
      end if;
      break;
    end if;

    if _self:-m_Verbose then
      printf(
        "LULEM::LU(...): M[%d,%d] = %a, cost = %d, degree_r = %d, degree_c = %d.\n",
        k, k, pivot["value"], pivot["cost"], pivot["degree_r"], pivot["degree_c"]
      );
    end if;

    # Shur complement
    tmp         := [k+1..-1];
    M[k, k]     := apply_veil(pivot["value"]);
    M[tmp, k]   := apply_veil~(Normalizer~(M[tmp, k])) / pivot["value"];
    M[k, tmp]   := apply_veil~(Normalizer~(M[k, tmp]));
    M[tmp, tmp] := apply_veil~(Normalizer~(M[tmp, tmp] - M[tmp, k].M[k, tmp]));
  end do;

  L := Matrix(M[1..m, 1..m], shape = triangular[lower, unit]);
  U := Matrix(M, shape = triangular[upper]);

  # Return the LU decomposition
  return table([
    "method" = "LU",
    "L"      = L,
    "U"      = U,
    "V"      = V,
    "r"      = r,
    "c"      = c,
    "rank"   = rnk,
    "pivots" = pivot_list,
    "L_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, L),
    "U_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, U),
    "V_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, _self:-m_LEM:-VeilList(_self:-m_LEM)),
    "L_nnz"  = nops(op(2, L)),
    "U_nnz"  = nops(op(2, U)),
    "A_nnz"  = nops(op(2, A))
  ]);
end proc: # LU

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export LUsolve::static := proc(
  _self::LULEM,
  T::table,
  b::Vector,
  $)::Vector;

  description "Solve the linear system Ax=b using LU decomposition <T> and "
    "provided the vector <b>.";

  local L, U, r, c, m, n, p, q, apply_veil, x, y, z, i, j, s, rnk:

  # Extract the LU decomposition
  L   := T["L"];
  U   := T["U"];
  r   := T["r"];
  c   := T["c"];
  rnk := T["rank"];

  # Get linear system dimension
  m, n := LinearAlgebra:-Dimensions(L);
  p, q := LinearAlgebra:-Dimensions(U);

  # Check if the linear system is consistent
  # sanity check
  if not ((m = n) and (p = q)) then
    error "only square system can be solved, got L = %d x %d, and U = %d x %d.",
      m, n, p, q;
    return table([]);
  end if;

  # Check if the linear system is consistent
  if not n = rnk then
    error "only full rank linear system can be solved (got rank = %1, expected "
      "rank = %2).", rnk, n;
    return table([]);
  end if;

  # Create a normalizer function
  apply_veil := (z) -> _self:-m_LEM:-Veil(_self:-m_LEM, z);

  # apply permutation P
  x := b[convert(r, list)];

  # Perform forward substitution to solve Ly=b[r]
  for i from 2 to m do
    if _self:-m_Verbose then
      printf("LULEM::LUsolve(...): forward substitution of %d-th row.\n", i);
    end if;
    x[i] := apply_veil(x[i] - add(L[i, 1..i-1] *~ x[1..i-1]));
  end do;

  # Perform backward substitution to solve Ux[c]=y
  if _self:-m_Verbose then
    printf("LULEM::LUsolve(...): dividision by U[%d,%d].\n", n, n);
  end if;
  x[n] := apply_veil(x[n] / U[n, n]);
  for i from n-1 to 1 by -1 do
    if _self:-m_Verbose then
      printf("LULEM::LUsolve(...): backward substitution of %d-th column.\n", i);
    end if;
    s    := apply_veil(x[i] - add(U[i, i+1..n] *~ x[i+1..n]));
    x[i] := apply_veil(s / U[i, i]);
  end do;

  # Apply inverse permutation Q
  y := Vector[column](n);
  for i from 1 to n do
    y[c[i]] := x[i];
  end do;

  # Return outputs
  return y;
end proc: # LUsolve

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
