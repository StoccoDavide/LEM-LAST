# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   ____             _
#  |  _ \ __ _ _ __ | | __
#  | |_) / _` | '_ \| |/ /
#  |  _ < (_| | | | |   <
#  |_| \_\__,_|_| |_|_|\_\
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Rank::static := proc(
  _self::LAST,
  A::Matrix, {
  reduced := false
  }, $)

  description "Compute the rank of a square matrix <A> by transforming it in "
    "row echelon form. If <reduced> is true, the matrix is transformed in "
    "reduced row echelon form.";

  local M, pivot, pivot_list, m, n, mn, i, k, rnk, r, c, tmp, tmp_LEM, veil_fnc;

  # Creare temporary LEM
  tmp_LEM  := _self:-m_LEM;
  veil_fnc := (lem, x) -> `if`(
    try
      timelimit(_self:-m_TimeLimit, evalb(lem:-Unveil(lem, x) = 0))
    catch:
      false
    end try,
    0,
    lem:-Veil(lem, x)
  );

  # Get matrix dimensions
  m, n := LinearAlgebra:-Dimensions(A):
  if _self:-m_VerboseMode then
    printf("LAST:-Rank(...): %d x %d matrix detected.\n", m, n);
  end if;

  # Create pivot vector
  r := Vector(m, k -> k);
  c := Vector(n, k -> k);

  M          := copy(A);
  mn         := min(m, n);
  rnk        := mn;
  pivot_list := [];

  # Transform in row echelon form
  for k from 1 to mn do
    if _self:-m_VerboseMode then
      printf(
        "LAST:-Rank(...): processing %d-th row, veilings = %d.\n",
        k,  nops(tmp_LEM:-VeilList(tmp_LEM))
      );
      #printf(
      #  "LAST:-Rank(...): processing %d-th row, cost = %d, veilings = %d.\n",
      #  k, tmp_LEM:-ExpressionCost(tmp_LEM, M), nops(tmp_LEM:-VeilList(tmp_LEM))
      #);
    end if;

    pivot := _self:-Pivoting(
      _self, k, M, r, c, parse("full_rows_degree") = false, parse("fast_pivoting") = true
    );
    if not pivot["is_zero"] then
      pivot_list := [op(pivot_list), pivot["value"]];
    end if;

    if pivot["is_zero"] then
      rnk := k - 1;
      if _self:-m_VerboseMode then
        WARNING("LAST:-Rank(...): the matrix appears to be not full rank.");
      end if;
      break;
    end if;

    if _self:-m_VerboseMode then
      printf(
        "LAST:-Rank(...): M[%d,%d] = %a, cost = %d, degree_r = %d, degree_c = %d.\n",
        pivot["i"], pivot["j"], pivot["value"], pivot["cost"], pivot["degree_r"],
        pivot["degree_c"]
      );
    end if;

    # Gaussian elimination
    tmp := [k..-1];
    if reduced then
      M[k, 1..-1] := veil_fnc~(tmp_LEM, Normalizer~(M[k, 1..-1]/pivot["value"]));
    end if;
    for i from k+1 to mn do
      if reduced then
        M[i, tmp] := veil_fnc~(tmp_LEM, Normalizer~(M[i, tmp] - M[i, k]*M[k, tmp]));
      else;
        M[i, tmp] := veil_fnc~(tmp_LEM, Normalizer~(M[i, tmp] - M[i, k]/M[k, k]*M[k, tmp]));
      end if;
    end do;
  end do;

  # Return the rank
  return rnk;
end proc: # Rank

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
