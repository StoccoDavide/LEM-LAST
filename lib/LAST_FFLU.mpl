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
  {
  warm_start::list({nonnegint, table}) := [0, table([])],
  veil_sanity_check::boolean           := true
  }, $)

  description "Compute the Fraction-Free LU (FFLU) decomposition of a square "
    "matrix <A> and check  if the veiling symbol is already present in the "
    "matrix coefficients.";

  local V, M, L, U, pivot, pivot_list, s_pivot, m, n, mn, k, j, rnk, r, c,
    last_GCD, tmp, tmp_try;

  # Check if the LEM object is initialized
  _self:-CheckInit(_self);

  # Get the veiling label
  V := _self:-m_LEM:-GetVeilingLabel(_self:-m_LEM);

  # Sanity check
  if veil_sanity_check and has(A, V) then
    error("veiling symbol %1 is already present in matrix coefficient.", V);
    return table([]);
  end if;

  # Get matrix dimensions
  m, n := LinearAlgebra:-Dimensions(A):

  # Create pivot vector
  r := Vector(m, k -> k);
  c := Vector(n, k -> k);

  M          := copy(A);
  mn         := min(m, n);
  rnk        := mn;
  pivot_list := [];

  # Perform Fraction-Free Gaussian elimination
  for k from 1 to mn do
    if _self:-m_VerboseMode then
      printf(
        "LAST:-FFLU(...): processing %d-th row, cost = %d, veilings = %d.\n",
        k, _self:-m_LEM:-ExpressionCost(_self:-m_LEM, M),
        nops(_self:-m_LEM:-VeilList(_self:-m_LEM))
      );
    end if;

    pivot := _self:-Pivoting(_self, k, M, r, c, parse("full_rows_degree") = false);

    if pivot["is_zero"] then
      rnk := k - 1;
      if _self:-m_VerboseMode then
        WARNING("LAST:-FFLU(...): the matrix appears to be not full rank.");
      end if;
      break;
    end if;

    if _self:-m_VerboseMode then
      printf(
        "LAST:-FFLU(...): M[%d,%d] = %a, cost = %d, degree_r = %d, degree_c = %d.\n",
        pivot["i"], pivot["j"], pivot["value"], pivot["cost"], pivot["degree_r"],
        pivot["degree_c"]
      );
    end if;
    tmp := [k+1..-1];

    # Scaled pivots (GCD)
    s_pivot := gcd~(M[tmp, k], pivot["value"]);
    tmp_try := M[tmp, k] /~ s_pivot;
    try
      M[tmp, k] := timelimit(_self:-m_TimeLimit, Normalizer~(tmp_try));
    catch "time expired":
      M[tmp, k] := tmp_try;
    end try;
    tmp_try := pivot["value"] /~ s_pivot;
    try
      s_pivot := timelimit(_self:-m_TimeLimit, Normalizer~(tmp_try));
    catch "time expired":
      s_pivot := tmp_try;
    end try;

    # Shur-like complement
    for j from k+1 to mn do
      M[j, tmp] := s_pivot[j-k] * M[j, tmp];
    end do;
    tmp_try := M[tmp, k].M[k, tmp] - M[tmp, tmp];
    try
      tmp_try := timelimit(_self:-m_TimeLimit, Normalizer~(tmp_try));
    catch:
      M[tmp, tmp] := tmp_try;
    end try;
    M[tmp, tmp] := _self:-m_LEM:-Veil~(_self:-m_LEM, tmp_try);

    # Save pivot list
    pivot_list := [op(pivot_list), s_pivot];
  end do;

  # Store the FFLU decomposition
  _self:-m_Results := table([
    "method" = "FFLU",
    "M"      = M,
    "V"      = V,
    "r"      = r,
    "c"      = c,
    "rank"   = rnk,
    "pivots" = pivot_list,
    "M_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, M),
    "V_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, _self:-m_LEM:-VeilList(_self:-m_LEM)),
    "M_nnz"  = nops(op(2, M)),
    "A_nnz"  = nops(op(2, A))
  ]);
  return NULL;
end proc: # FFLU

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export FFLUsolve::static := proc(
  _self::LAST,
  b::Vector,
  $)::Vector;

  description "Apply L^(-1)*P to the vector <b>.";

  local M, r, c, m, n, x, y, i, s, rnk, pivots, tmp;

  # Check if the LEM object is initialized
  _self:-CheckInit(_self);

  # Check if the results are available
  _self:-CheckResults(_self);

  # Check if the FFLU decomposition is available
  if not (_self:-m_Results["method"] = "FFLU") then
    error("wrong or not available FFLU decomposition (use 'LAST:-FFLU()' first).");
  end if;

  # Extract the FFLU decomposition
  M      := _self:-m_Results["M"];
  pivots := _self:-m_Results["pivots"];
  r      := _self:-m_Results["r"];
  c      := _self:-m_Results["c"];
  rnk    := _self:-m_Results["rank"];

  # Get linear system dimension
  m, n := LinearAlgebra:-Dimensions(M);

  # Check if the linear system is consistent (sanity check)
  if not (m = n) then
    error("only square system can be solved, got M = %1 x %2.", m, n);
  end if;

  # Check if the linear system is consistent
  if not (n = rnk) then
    error("only full rank linear system can be solved (got rank = %1, expected "
      "rank = %2).", rnk, n);
  end if;

  # Apply permutation P
  x := b[convert(r, list)];

  # Perform forward substitution to solve Ly=b[r]
  for i from 1 to m-1 do
    if _self:-m_VerboseMode then
      printf("LAST:-FFLUsolve(...): backward substitution of %d-th column.\n", i);
    end if;
    tmp    := [i+1..-1];
    x[tmp] := _self:-m_LEM:-Veil~(_self:-m_LEM, M[tmp, i]*x[i] - pivots[i]*~x[tmp]);
  end do;

  # Perform backward substitution to solve Ux[c]=y
  if _self:-m_VerboseMode then
    printf("LAST:-LUsolve(...): division by U[%d,%d].\n", n, n);
  end if;
  x[n] := _self:-m_LEM:-Veil(_self:-m_LEM, x[n]/M[n, n]);
  for i from n-1 to 1 by -1 do
    if _self:-m_VerboseMode then
      printf("LAST:-FFLUsolve(...): backward substitution of %d-th column.\n", i);
    end if;
    s    := _self:-m_LEM:-Veil(_self:-m_LEM, x[i] - add(M[i, i+1..n]*~x[i+1..n]));
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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export FFLUapplyLP::static := proc(
  _self::LAST,
  b::Vector,
  $)::Vector;

  description "Return the matrix U^T*Q.";

  local M, r, x, i, rnk, pivots, tmp;

  # Check if the LEM object is initialized
  _self:-CheckInit(_self);

  # Check if the results are available
  _self:-CheckResults(_self);

  # Check if the FFLU decomposition is available
  if not (_self:-m_Results["method"] = "FFLU") then
    error("wrong or not available FFLU decomposition (use 'LAST:-FFLU()' first).");
  end if;

  # Extract the FFLU decomposition
  M      := _self:-m_Results["M"];
  pivots := _self:-m_Results["pivots"];
  r      := _self:-m_Results["r"];
  rnk    := _self:-m_Results["rank"];

  # Apply permutation P
  x := b[convert(r, list)];

  # Perform forward substitution to solve Ly=b[r]
  for i from 1 to rnk do
    if _self:-m_VerboseMode then
      printf("LAST:-FFLUsolve(...): backward substitution of %d-th column.\n", i);
    end if;
    tmp    := [i+1..-1];
    x[tmp] := _self:-m_LEM:-Veil~(_self:-m_LEM, M[tmp, i]*x[i] - pivots[i]*~x[tmp]);
  end do;

  # Return outputs
  return x;
end proc: # FFLUapplyLP

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export FFLUgetUQT::static := proc(
  _self::LAST,
  $)::Matrix;

  description "Return the matrix U^T*Q.";

  local M, U, UT, c, n, i, rnk;

  # Check if the LEM object is initialized
  _self:-CheckInit(_self);

  # Check if the results are available
  _self:-CheckResults(_self);

  # Check if the FFLU decomposition is available
  if not (_self:-m_Results["method"] = "FFLU") then
    error("wrong or not available FFLU decomposition (use 'LAST:-FFLU()' first).");
  end if;

  # Extract the FFLU decomposition
  M   := _self:-m_Results["M"];
  U   := Matrix(M, shape = triangular[upper]);
  c   := _self:-m_Results["c"];
  rnk := _self:-m_Results["rank"];
  n   := LinearAlgebra:-ColumnDimension(U);
  UT  := Matrix(rnk, n);

  # Apply inverse permutation Q
  for i from 1 to n do
    UT[1..-1, c[i]] := U[1..rnk, i];
  end do;

  # Return outputs
  return UT;
end proc: # FFLUgetUQT

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
