# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   _    _   _
#  | |  | | | |
#  | |  | | | |
#  | |__| |_| |
#  |_____\___/
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export LU::static := proc(
  _self::LAST,
  A::Matrix,
  $)

  description "Compute the LU decomposition of a square matrix <A>.";

  local V, M, L, U, pivot, pivot_list, m, n, mn, k, rnk, r, c, tmp;

  # Check if the LEM is initialized
  if not type(_self:-m_LEM, LEM) then
    error "LEM is not initialized (use LAST::InitLEM() first).";
  end if;

  # Get the veiling label
  V := _self:-m_LEM:-GetVeilingLabel(_self:-m_LEM);

  # Sanity check
  if has(A, V) then
    error "veiling symbol %1 is already present in matrix coefficient.", V;
    return table([]);
  end if;

  # Forget the veilings
  _self:-m_LEM:-VeilForget(_self:-m_LEM);

  # Get matrix dimensions
  m, n := LinearAlgebra:-Dimensions(A):

  # Create pivot vector
  r := Vector(m, k -> k);
  c := Vector(n, k -> k);

  # Perform Gaussian elimination
  M          := copy(A);
  mn         := min(m, n);
  rnk        := mn;
  pivot_list := [];
  for k from 1 to mn do
    if _self:-m_VerboseMode then
      printf(
        "LAST::LU(...): processing %d-th row, cost = %d, veilings = %d.\n",
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
      if _self:-m_VerboseMode then
        WARNING("LAST::LU(...): the matrix appears to be not full rank.");
      end if;
      break;
    end if;

    if _self:-m_VerboseMode then
      printf(
        "LAST::LU(...): M[%d,%d] = %a, cost = %d, degree_r = %d, degree_c = %d.\n",
        k, k, pivot["value"], pivot["cost"], pivot["degree_r"], pivot["degree_c"]
      );
    end if;

    # Shur complement
    tmp         := [k+1..-1];
    M[k, k]     := _self:-m_LEM:-Veil(_self:-m_LEM, pivot["value"]);
    M[tmp, k]   := _self:-m_LEM:-Veil~(_self:-m_LEM, Normalizer~(M[tmp, k])) / pivot["value"];
    M[k, tmp]   := _self:-m_LEM:-Veil~(_self:-m_LEM, Normalizer~(M[k, tmp]));
    M[tmp, tmp] := _self:-m_LEM:-Veil~(_self:-m_LEM, Normalizer~(M[tmp, tmp] - M[tmp, k].M[k, tmp]));
  end do;

  L := Matrix(M[1..m, 1..m], shape = triangular[lower, unit]);
  U := Matrix(M, shape = triangular[upper]);

  # Store the LU decomposition
  _self:-m_Results := table([
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
  return NULL;
end proc: # LU

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export LUsolve::static := proc(
  _self::LAST,
  b::Vector,
  $)::Vector;

  description "Solve the linear system Ax=b using LU decomposition provided "
  "the vector <b>.";

  local L, U, r, c, m, n, p, q, x, y, i, s, rnk;

  # Check if the LEM is initialized
  if not type(_self:-m_LEM, LEM) then
    error "LEM is not initialized (use LAST::InitLEM() first).";
  end if;

  # Check if the LU decomposition is available
  if not (_self:-m_Results["method"] = "LU") then
    error "wrong or not available LU decomposition (use LAST::LU() first).";
  end if;

  # Extract the LU decomposition
  L   := _self:-m_Results["L"];
  U   := _self:-m_Results["U"];
  r   := _self:-m_Results["r"];
  c   := _self:-m_Results["c"];
  rnk := _self:-m_Results["rank"];

  # Get linear system dimension
  m, n := LinearAlgebra:-Dimensions(L);
  p, q := LinearAlgebra:-Dimensions(U);

  # Check if the linear system is consistent
  # sanity check
  if not ((m = n) and (p = q)) then
    error "only square system can be solved, got L = %d x %d, and U = %d x %d.",
      m, n, p, q;
  end if;

  # Check if the linear system is consistent
  if not n = rnk then
    error "only full rank linear system can be solved (got rank = %1, expected "
      "rank = %2).", rnk, n;
  end if;

  # apply permutation P
  x := b[convert(r, list)];

  # Perform forward substitution to solve Ly=b[r]
  for i from 2 to m do
    if _self:-m_VerboseMode then
      printf("LAST::LUsolve(...): forward substitution of %d-th row.\n", i);
    end if;
    x[i] := _self:-m_LEM:-Veil(_self:-m_LEM, x[i] - add(L[i, 1..i-1]*~x[1..i-1]));
  end do;

  # Perform backward substitution to solve Ux[c]=y
  if _self:-m_VerboseMode then
    printf("LAST::LUsolve(...): dividision by U[%d,%d].\n", n, n);
  end if;
  x[n] := _self:-m_LEM:-Veil(_self:-m_LEM, x[n]/U[n, n]);
  for i from n-1 to 1 by -1 do
    if _self:-m_VerboseMode then
      printf("LAST::LUsolve(...): backward substitution of %d-th column.\n", i);
    end if;
    s    := _self:-m_LEM:-Veil(_self:-m_LEM, x[i] - add(U[i, i+1..n]*~x[i+1..n]));
    x[i] := _self:-m_LEM:-Veil(_self:-m_LEM, s/U[i, i]);
  end do;

  # Apply inverse permutation Q
  y := Vector[column](n);
  for i from 1 to n do
    y[c[i]] := x[i];
  end do;

  print("here1");

  # Return outputs
  return y;
end proc: # LUsolve

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
