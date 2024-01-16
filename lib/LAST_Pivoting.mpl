# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   ____  _            _   _
#  |  _ \(_)_   _____ | |_(_)_ __   __ _
#  | |_) | \ \ / / _ \| __| | '_ \ / _` |
#  |  __/| |\ V / (_) | |_| | | | | (_| |
#  |_|   |_| \_/ \___/ \__|_|_| |_|\__, |
#                                  |___/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Pivoting::static := proc(
  _self::LAST,
  k::integer,
  M::Matrix,
  r::Vector(nonnegint),
  c::Vector(nonnegint),
  {
  full_rows_degree::boolean := false
  }, $)::table;

  description "Compute the LU decomposition pivots vectors with minum degree "
    "provided the step <k>, the temporary LU (NAG) matrix <M>, the rows "
    "permutation <r> and the columns permutation <c>.";

  local uMij, M_degree_R, M_degree_C, perm, perm_R, perm_C, m, n, i, j, ipos,
    Mij, pivot, pivot_list, pivot_cost, M_data, V, V_data;

  # Copy the matrix and veils, and substitute the data
  M_data := subs(op(_self:-m_StoredData), M);
  V      := _self:-m_LEM:-VeilList(_self:-m_LEM, parse("reverse") = true);
  V_data := subs(op(_self:-m_StoredData), V);

  # Check if the LEM object is initialized
  _self:-CheckInit(_self);

  # Extract matrix dimensions
  m, n := LinearAlgebra:-Dimensions(M):

  # Calculate the degree
  M_degree_R := Matrix(m, n);
  M_degree_C := Matrix(m, n);
  if full_rows_degree then
    M_degree_R[1..-1, k..n], M_degree_C[1..-1, k..n] :=
      _self:-GetDegrees(_self, M_data[1..-1, k..n]);
  else
    M_degree_R[k..m, k..n], M_degree_C[k..m, k..n] :=
      _self:-GetDegrees(_self, M_data[k..m, k..n]);
  end if;

  # Build a list (i,j,degree,cost) and sort it
  pivot      := table([]);
  Mij        := table([]);
  pivot_list := Vector((n-k+1)*(m-k+1));
  pivot_cost := Vector((n-k+1)*(m-k+1));
  ipos       := 1;
  for j from k to n do
    for i from k to m do
      # Look for a non-zero pivot
      Mij["i"]           := i;
      Mij["j"]           := j;
      Mij["degree_r"]    := M_degree_R[i, j];
      Mij["degree_c"]    := M_degree_C[i, j];
      Mij["degree_cost"] := _self:-DegreeCost(_self, Mij);
      Mij["is_zero"]     := false;
      pivot_list[ipos]   := copy(Mij);
      pivot_cost[ipos]   := Mij["degree_cost"];
      ipos               := ipos + 1;
    end do;
  end do;

  # Sort the pivot list by degree cost
  perm := sort(pivot_cost, output = 'permutation');

  # Iterate over the columns and rows using estimated increasing cost o pivot
  pivot["is_zero"] := true;
  for ipos in perm do
    Mij          := copy(pivot_list[ipos]);
    i            := Mij["i"];
    j            := Mij["j"];
    Mij["value"] := M[i, j];
    try
      Mij["sig"] := _self:-m_LEM:-Signature(_self:-m_LEM, Mij["value"]);
      Mij["is_zero"] := evalb(Mij["value"] = 0 or Mij["sig"] = 0);
    catch "division by zero":
      if _self:-m_WarningMode then
        WARNING(
          "LAST:-Pivoting(...): division by zero, assumed 'M[%1,%2] = 0'.",
          i, j
        );
      end if;
      Mij["is_zero"] := true;
    catch:
      if _self:-m_WarningMode then
        WARNING(
          "LAST:-Pivoting(...): something went wrong in '%1', assumed 'M[%2,%3] <> 0'.",
          lastexception, i, j
        );
      end if;
      Mij["is_zero"] := false;
    end try;

    # if zero skip
    if Mij["is_zero"] then continue; end if;

    # Look for a non-zero pivot
    # Pre-check if next pivot cannot improve search
    if not pivot["is_zero"] and (Mij["degree_cost"] > pivot["degree_cost"]) then
      # If the degree cost is higher than degree cost of pivot cannot improve
      break;
    end if;

    # Try to simplify the pivot expression
    try
      Mij["value"] := timelimit(_self:-m_TimeLimit, Normalizer(Mij["value"]));
      Mij["is_zero"] := evalb(Mij["value"] = 0);
      if not Mij["is_zero"] then
        uMij := timelimit(_self:-m_TimeLimit, subs(op(_self:-m_StoredData), op(V_data), Mij["value"]));
        # Recalculate cost and value of the pivot
        Mij["cost"], Mij["numeric_value"] := _self:-PivotCost(_self, uMij);
        # Time limit required because sometimes normalizer get stuck
        uMij := timelimit(_self:-m_TimeLimit, eval(Normalizer(uMij)));
        Mij["is_zero"] := evalb(uMij = 0);
      end if;
    catch "time expired":
      if _self:-m_WarningMode then
        WARNING(
          "LAST:-Pivoting(...): time expired, assumed 'M[%1,%2] <> 0'.",
          i, j
        );
      end if;
      Mij["is_zero"] := false;
      uMij := subs(op(_self:-m_StoredData), Mij["value"]);
      # Recalculate cost and value of the pivot
      Mij["cost"], Mij["numeric_value"] := _self:-PivotCost(_self, uMij);
    catch "division by zero":
      if _self:-m_WarningMode then
        WARNING(
          "LAST:-Pivoting(...): division by zero, assumed 'M[%1,%2] = 0'.",
          i, j
        );
      end if;
      Mij["is_zero"] := true;
    catch:
      WARNING(
        "LAST:-Pivoting(...): something went wrong in '%1', assumed 'M[%2,%3] <> 0'.",
        lastexception, i, j
      );
      if _self:-m_VerboseMode then
        print("catch", Mij["value"] , Mij["sig"], Mij["is_zero"]);
      end if;
      Mij["is_zero"] := false;
    end try;

    if Mij["is_zero"] then
      M[i, j] := 0;
    else
      # Found a non-zero pivot, check if it is better
      if pivot["is_zero"] then
        # First non-zero pivot found
        pivot := copy(Mij);
      elif _self:-PivotingCompare(_self, pivot, Mij) then
        # A better pivot is found
        pivot := copy(Mij);
      end if;
    end if;

  end do;

  # Swap rows and columns
  if not pivot["is_zero"] then
    i := pivot["i"];
    j := pivot["j"];
    if (i <> k) then
      (r[i], r[k])     := (r[k], r[i]);
      M[[i, k], 1..-1] := M[[k, i], 1..-1];
    end if;
    if (j <> k) then
      (c[j], c[k])     := (c[k], c[j]);
      M[1..-1, [j, k]] := M[1..-1, [k, j]];
    end if;
  end if;

  return pivot;
end proc: # Pivoting

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PivotCost::static := proc(
  _self::LAST,
  x::algebraic,
  $)::integer, algebraic;

  description "Compute the cost of the pivot <x>.";

  if type(x, numeric) then
    if evalb(x = 0) then
      return 0, 0;
    else
      return 1, abs(x);
    end if;
  end if;
  return _self:-m_LEM:-ExpressionCost(_self:-m_LEM, x), infinity;
end proc: # PivotCost

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   __  __ _       ____
#  |  \/  (_)_ __ |  _ \  ___  __ _ _ __ ___  ___
#  | |\/| | | '_ \| | | |/ _ \/ _` | '__/ _ \/ _ \
#  | |  | | | | | | |_| |  __/ (_| | | |  __/  __/
#  |_|  |_|_|_| |_|____/ \___|\__, |_|  \___|\___|
#                             |___/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Code inspired by:
# https://www.mapleprimes.com/questions/235996-Is-There-Any-Command-Or-Function-For

export GetDegrees::static := proc(
  _self::LAST,
  A::Matrix,
  $)::Matrix(nonnegint), Matrix(nonnegint);

  description "Get the degree matrices of the matrix <A>.";

  local i, j, k, m, n, r, c, ro, co;

  m, n := LinearAlgebra:-Dimensions(A);
  ro := Vector[column](m, k -> 1);
  co := Vector[row](n, k -> 1);
  r  := Vector[column](
    [seq(rtable_scanblock(A, [i,..], ':-NonZeros'), i = 1..m)]
  );
  c  := Vector[row](
    [seq(rtable_scanblock(A, [..,j], ':-NonZeros'), j = 1..n)]
  );
  return r.co, ro.c;
end proc: # GetDegrees

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export SetMinDegreeStrategy::static := proc(
  _self::LAST,
  str::string := "product_1rc",
  $)

  description "Set the strategy <str> for the minimum degree ordering.";

  if not str in ["none", "row", "column", "sum", "product", "product_1rc",
    "product_1r", "product_1c"] then
    error("unknown minimum degree strategy %1.", str);
  else
    _self:-m_MinDegreeStrategy := str;
  end if;
  return NULL;
end proc: # SetMinDegreeStrategy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DegreeCost::static := proc(
  _self::LAST,
  val::table,
  $)::integer;

  description "Set the strategy <val> for the minimum degree ordering.";

  if (_self:-m_MinDegreeStrategy = "none") then
    return _self:-DegreeCost_none(_self, val);
  elif (_self:-m_MinDegreeStrategy = "row") then
    return _self:-DegreeCost_row(_self, val);
  elif (_self:-m_MinDegreeStrategy = "column") then
    return _self:-DegreeCost_col(_self, val);
  elif (_self:-m_MinDegreeStrategy = "sum") then
    return _self:-DegreeCost_sum(_self, val);
  elif (_self:-m_MinDegreeStrategy = "product") then
    return _self:-DegreeCost_prod(_self, val);
  elif (_self:-m_MinDegreeStrategy = "product_1rc") then
    return _self:-DegreeCost_prod_1rc(_self, val);
  elif (_self:-m_MinDegreeStrategy = "product_1r") then
    return _self:-DegreeCost_prod_1r(_self, val);
  elif (_self:-m_MinDegreeStrategy = "product_1c") then
    return _self:-DegreeCost_prod_1c(_self, val);
  elif (_self:-m_MinDegreeStrategy = "minimum") then
    return _self:-DegreeCost_min(_self, val);
  elif (_self:-m_MinDegreeStrategy = "maximum") then
    return _self:-DegreeCost_max(_self, val);
  else
    error("unknown minimum degree strategy %1.", val);
    return NULL;
  end if;
end proc: # SetMinDegreeStrategy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DegreeCost_none::static := proc(
  _self::LAST,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return 0;
end proc: # DegreeCost_none

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DegreeCost_row::static := proc(
  _self::LAST,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return val["degree_r"];
end proc: # DegreeCost_row

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DegreeCost_col::static := proc(
  _self::LAST,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return val["degree_c"];
end proc: # DegreeCost_col

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DegreeCost_sum::static := proc(
  _self::LAST,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return val["degree_c"] + val["degree_r"];
end proc: # DegreeCost_sum

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DegreeCost_prod::static := proc(
  _self::LAST,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return val["degree_c"] * val["degree_r"];
end proc: # DegreeCost_prod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DegreeCost_prod_1rc::static := proc(
  _self::LAST,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return max(val["degree_c"]-1, 0) * val["degree_r"] +
         max(val["degree_r"]-1, 0) * val["degree_c"];
end proc: # DegreeCost_prod_1rc

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DegreeCost_prod_1r::static := proc(
  _self::LAST,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return max(val["degree_c"], 0) * max(val["degree_r"]-1, 0);
end proc: # DegreeCost_prod_1r

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DegreeCost_prod_1c::static := proc(
  _self::LAST,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return max(val["degree_c"]-1, 0) * max(val["degree_r"], 0);
end proc: # DegreeCost_prod_1c

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DegreeCost_min::static := proc(
  _self::LAST,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return min(val["degree_c"], val["degree_r"]);
end proc: # DegreeCost_min

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DegreeCost_max::static := proc(
  _self::LAST,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return max(val["degree_c"], val["degree_r"]);
end proc: # DegreeCost_max

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   ____  _            _   _
#  |  _ \(_)_   _____ | |_(_)_ __   __ _
#  | |_) | \ \ / / _ \| __| | '_ \ / _` |
#  |  __/| |\ V / (_) | |_| | | | | (_| |
#  |_|   |_| \_/ \___/ \__|_|_| |_|\__, |
#                                  |___/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PivotingCompare::static := proc(
  _self::LAST,
  cur::table,
  val::table,
  $)::boolean;

  description "Compute the pivoting strategy: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current pivot or not.";

  if (val["numeric_value"] = infinity) and (cur["numeric_value"] = infinity) then
    # Both are expressions: use 'cost' to compare them
    return evalb(val["cost"] < cur["cost"]);
  elif (val["numeric_value"] = infinity) then
    return false;
  elif (cur["numeric_value"] = infinity) then
    return true;
  else
    return evalb(val["numeric_value"] > cur["numeric_value"]);
  end if;
end proc: # PivotingCompare

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
