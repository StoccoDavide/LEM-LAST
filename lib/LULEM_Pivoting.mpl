# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   ____  _            _   _
#  |  _ \(_)_   _____ | |_(_)_ __   __ _
#  | |_) | \ \ / / _ \| __| | '_ \ / _` |
#  |  __/| |\ V / (_) | |_| | | | | (_| |
#  |_|   |_| \_/ \___/ \__|_|_| |_|\__, |
#                                  |___/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Pivoting::static := proc(
  _self::LULEM,
  k::integer,
  M::Matrix,
  r::Vector(nonnegint),
  c::Vector(nonnegint),
  $)::table;

  description "Compute the LU decomposition pivots vectors with minum degree "
    "provided the step <k>, the temporary LU (NAG) matrix <M>, the rows "
    "permutation <r> and the columns permutation <c>.";

  local uMij, M_degree_R, M_degree_C, perm, perm_R, perm_C, m, n, i, j, ipos,
    apply_unveil, z, Mij, pivot, pivot_list, pivot_cost;

  # Check if to veil or not to veil
  apply_unveil := (z) -> _self:-m_LEM:-UnVeil(_self:-m_LEM, z);

  # Check degree cost function
  if (_self:-m_DegreeCost_fun = NULL) then
    _self:-SetMinDegreeStrategy(_self);
  end if;

  # Extract matrix dimensions
  m, n := LinearAlgebra:-Dimensions(M):

  # Calculate the degree
  M_degree_R := Matrix(m, n);
  M_degree_C := Matrix(m, n);
  M_degree_R[k..m, k..n], M_degree_C[k..m, k..n] := _self:-GetDegrees(_self, M[k..m, k..n]);

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
      Mij["degree_cost"] := _self:-m_DegreeCost_fun(_self, Mij);
      pivot_list[ipos]   := copy(Mij);
      pivot_cost[ipos]   := Mij["degree_cost"];
      ipos               := ipos+1;
    end do;
  end do;

  # Sort the pivot list by degree cost
  perm := sort( pivot_cost, output = 'permutation' );

  # Iterate over the columns and rows using estimated increasing cost o pivot
  pivot["is_zero"] := true;
  for ipos in perm do
    Mij          := copy(pivot_list[ipos]);
    i            := Mij["i"];
    j            := Mij["j"];
    Mij["value"] := M[i,j];

    # Look for a non-zero pivot
    # Pre-check if next pivot cannot improve search
    if not pivot["is_zero"] and Mij["degree_cost"] > pivot["degree_cost"] then
      # If the degree cost is higher than degree cost of pivot cannot improve
      break;
    end if;

    # Try to simplify the pivot expression
    try
      Mij["value"]   := Normalizer(Mij["value"]);
      Mij["is_zero"] := evalb(Mij["value"] = 0);
      if not Mij["is_zero"] then
        uMij := timelimit(_self:-m_TimeLimit, apply_unveil(Mij["value"]));
        # Recalculate cost and value of the pivot
        Mij["cost"], Mij["numeric_value"] := _self:-PivotCost(_self, uMij);
        # Time limit required because sometimes Normalizer get stuck
        uMij := timelimit(_self:-m_TimeLimit, eval(Normalizer(uMij)));
        Mij["is_zero"] := evalb(uMij = 0);
      end if;
    catch "time expired":
      printf("LULEM::Pivoting(...): simplify(Mij) failed, assumed <> 0.\n");
      Mij["is_zero"] := false;
    catch:
      printf("LULEM::Pivoting(...): Mij division by 0 or other exception.\n");
      if _self:-m_Verbose then
        print(Mij["value"]);
      end if;
      Mij["is_zero"] := true;
    end try;

    if Mij["is_zero"] then
      M[i,j] := 0;
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
  _self::LULEM,
  x::algebraic,
  $)::integer, algebraic;

  description "Compute the cost of the pivot <x>.";

  if type(x, integer) or type(x, float) then
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

export SetMinDegreeStrategy::static := proc(
  _self::LULEM,
  str::string := "prod_1",
  $)

  description "Set the strategy <str> for the minimum degree ordering.";

  if (str = "none") then
    _self:-m_DegreeCost_fun := _self:-DegreeCost_none;
  elif (str = "row") then
    _self:-m_DegreeCost_fun := _self:-DegreeCost_row;
  elif (str = "col") then
    _self:-m_DegreeCost_fun := _self:-DegreeCost_col;
  elif (str = "sum") then
    _self:-m_DegreeCost_fun := _self:-DegreeCost_sum;
  elif (str = "prod") then
    _self:-m_DegreeCost_fun := _self:-DegreeCost_prod;
  elif (str = "prod_1") then
    _self:-m_DegreeCost_fun := _self:-DegreeCost_prod_1;
  elif (str = "min") then
    _self:-m_DegreeCost_fun := _self:-DegreeCost_min;
  elif (str = "max") then
    _self:-m_DegreeCost_fun := _self:-DegreeCost_max;
  else
    error "unknown minimum degree strategy %1.", str;
  end if;
  return NULL;
end proc: # SetMinDegreeStrategy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

local DegreeCost_none::static := proc(
  _self::LULEM,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return 0;
end proc: # DegreeCost_none

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

local DegreeCost_row::static := proc(
  _self::LULEM,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return val["degree_r"];
end proc: # DegreeCost_row

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

local DegreeCost_col::static := proc(
  _self::LULEM,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return val["degree_c"];
end proc: # DegreeCost_col

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

local DegreeCost_sum::static := proc(
  _self::LULEM,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return val["degree_c"] + val["degree_r"];
end proc: # DegreeCost_sum

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

local DegreeCost_prod::static := proc(
  _self::LULEM,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return val["degree_c"] * val["degree_r"];
end proc: # DegreeCost_prod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

local DegreeCost_prod_1::static := proc(
  _self::LULEM,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return max(val["degree_c"]-1, 0) * max(val["degree_r"]-1, 0);
end proc: # DegreeCost_prod_1

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

local DegreeCost_min::static := proc(
  _self::LULEM,
  val::table,
  $)::integer;

  description "Compute the pivoting degree cost.";

  return min(val["degree_c"], val["degree_r"]);
end proc: # DegreeCost_min

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

local DegreeCost_max::static := proc(
  _self::LULEM,
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

local PivotingCompare::static := proc(
  _self::LULEM,
  cur::table,
  val::table,
  $)::boolean;

  description "Compute the pivoting strategy: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current pivot or not.";

  if (val["numeric_value"] < cur["numeric_value"]) then
    return true;
  elif (val["numeric_value"] > cur["numeric_value"]) then
    return false;
  elif (val["numeric_value"] = infinity) then
    return evalb(val["cost"] < cur["cost"]);
  else
    return evalb(val["numeric_value"] > cur["numeric_value"]);
  end if;
end proc: # PivotingCompare

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
