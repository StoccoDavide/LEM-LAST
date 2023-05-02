# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    ____     _
#   / ___|   | |
#  | |  _ _  | |
#  | |_| | |_| |
#   \____|\___/
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GJ::static := proc(
  _self::LAST,
  A::Matrix,
  {
  veil_sanity_check::boolean := true
  }, $)

  description "Compute the Gauss-Jordan decomposition of a rectangular matrix "
    "<A> and check  if the veiling symbol is already present in the matrix "
    "coefficients.";

  local V, M, pivot, pivot_list, m, n, mn, k, rnk, r, c, tr, tc;

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

  # Perform Gaussian-Jordan elimination
  M          := copy(A);
  mn         := min(m, n);
  rnk        := mn;
  pivot_list := [];
  for k from 1 to mn do
    if _self:-m_VerboseMode then
      printf(
        "LAST:-GJ(...): processing %d-th row, cost = %d, veilings = %d.\n",
        k, _self:-m_LEM:-ExpressionCost(_self:-m_LEM, M),
        nops(_self:-m_LEM:-VeilList(_self:-m_LEM))
      );
    end if;

    pivot := _self:-Pivoting(
      _self, k, M, r, c, parse("full_rows_degree") = true
    );
    if not pivot["is_zero"] then
      pivot_list := [op(pivot_list), pivot["value"]];
    end if;

    if pivot["is_zero"] then
      rnk := rnk - 1;
      if _self:-m_VerboseMode then
        WARNING("LAST:-GJ(...): the matrix appears to be not full rank.");
      end if;
      break;
    end if;

    if _self:-m_VerboseMode then
      printf(
        "LAST:-GJ(...): M[%d,%d] = %a, cost = %d, degree_r = %d, degree_c = %d.\n",
        pivot["i"], pivot["j"], pivot["value"], pivot["cost"], pivot["degree_r"],
        pivot["degree_c"]
      );
    end if;

    # Multiply by D_k^(-1) and then by I-v*e_k^T
    if (k = 1) then
      tr := [2..-1];
    elif (k = m) then
      tr := [1..m-1];
    else
      tr := [1..k-1, k+1..-1];
    end;
    tc        := [k+1..-1];
    M[k,  k]  := _self:-m_LEM:-Veil(_self:-m_LEM, pivot["value"]);
    M[k,  tc] := _self:-m_LEM:-Veil~(_self:-m_LEM, Normalizer~(M[k, tc]/pivot["value"]));
    M[tr, k]  := _self:-m_LEM:-Veil~(_self:-m_LEM, Normalizer~(M[tr, k]));
    M[tr, tc] := _self:-m_LEM:-Veil~(_self:-m_LEM, Normalizer~(M[tr, tc]-M[tr, k].M[k,tc]));
  end do;

  # Store the LU decomposition
  _self:-m_Results := table([
    "method" = "GJ",
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
end proc: # LU

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GJsolve::static := proc(
  _self::LAST,
  b::Vector,
  $)::Vector;

  description "Solve the linear system Ax=b using GJ decomposition provided "
  "the vector <b>.";

  local M, r, c, m, n, p, q, x, y, k, tr, tmp, rnk;

  # Check if the LEM object is initialized
  _self:-CheckInit(_self);

  # Check if the results are available
  _self:-CheckResults(_self);

  # Check if the LU decomposition is available
  if not (_self:-m_Results["method"] = "GJ") then
    error("wrong or not available GJ decomposition (use 'LAST:-GJ()' first).");
  end if;

  # Extract the GJ decomposition
  M   := _self:-m_Results["M"];
  r   := _self:-m_Results["r"];
  c   := _self:-m_Results["c"];
  rnk := _self:-m_Results["rank"];

  # Get linear system dimension
  m, n := LinearAlgebra:-Dimensions(M);

  # Check if the linear system is consistent (sanity check)
  if not (m = n) then
    error("only square system can be solved, got M = %d x %d.", m, n);
  end if;

  # Check if the linear system is consistent
  if not (n = rnk) then
    error("only full rank linear system can be solved (got rank = %1, expected "
      "rank = %2).", rnk, n);
  end if;

  # Apply permutation P
  if _self:-m_VerboseMode then
    printf("LAST:-GJsolve(...): Apply permutation P\n");
  end if;
  x := b[convert(r, list)];

  # Apply multiplications in reverse order
  for k from 1 to n do
    if _self:-m_VerboseMode then
      printf("LAST:-GJsolve(...): step k = %d.\n", k);
    end if;
    # Multiply by D_k^(-1) and then by I-v*e_k^T
    # Try to simplify the pivot expression
    try
      tmp := x[k]/M[k,k];
      x[k] := timelimit(_self:-m_TimeLimit, Normalizer(tmp) );
    catch:
      if _self:-m_VerboseMode then
        printf("LAST:-GJsolve(...): step k = %d, timelimit(1) reached.\n", k);
      end if;
      x[k] := tmp;
    end try;
    #x[k] := _self:-m_LEM:-Veil(_self:-m_LEM, x[k] );
    if (k = 1) then
      tr := [2..-1];
    elif (k = m) then
      tr := [1..m-1];
    else
      tr := [1..k-1,k+1..-1];
    end if;
    try
      tmp := x[tr]-x[k]*M[tr, k];
      x[tr] := timelimit(_self:-m_TimeLimit, Normalizer~(tmp));
    catch:
      if _self:-m_VerboseMode then
        printf("LAST:-GJsolve(...): step k = %d, timelimit(2) reached.\n", k);
      end if;
      x[tr] := tmp;
    end try;
    x[tr] := _self:-m_LEM:-Veil~(_self:-m_LEM, x[tr]);
  end do;

  # Apply inverse permutation Q
  if _self:-m_VerboseMode then
    printf("LAST:-GJsolve(...): applying inverse permutation Q.\n");
  end if;
  y := Vector[column](n);
  for k from 1 to n do
    y[c[k]] := x[k];
  end do;

  # Return outputs
  return y;
end proc: # LUsolve

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
