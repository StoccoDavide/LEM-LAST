# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   ____  _            _   _
#  |  _ \(_)_   _____ | |_(_)_ __   __ _
#  | |_) | \ \ / / _ \| __| | '_ \ / _` |
#  |  __/| |\ V / (_) | |_| | | | | (_| |
#  |_|   |_| \_/ \___/ \__|_|_| |_|\__, |
#                                  |___/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Pivoting := proc(
  k::{integer},
  M::{Matrix},
  V::{symbol},
  r::{Vector(nonnegint)},
  c::{Vector(nonnegint)},
  $)::table;

  description "Compute the LU decomposition pivots vectors provided the step "
    "<k>, the temporary LU (NAG) matrix <M>, the veiling symbol <V>, the rows "
    "the pivot vector <r>, the columns the pivot vector <c>.";

  local uMij, M_degree_R, M_degree_C, m, n, i, j, apply_unveil, z, Mij, pivot;

  # Extract matrix dimensions
  m, n := LinearAlgebra:-Dimensions(M):

  # Check if to veil or not to veil
  apply_unveil := (z) -> LEM:-UnVeil[V](z);

  # Calculate the degree
  M_degree_R := Matrix(m,n);
  M_degree_C := Matrix(m,n);
  M_degree_R[k..m, k..n], M_degree_C[k..m, k..n] := LULEM:-GetDegrees(M[k..m, k..n]);

  pivot := table([]);
  Mij   := table([]);

  # Iterate over the columns and rows
  pivot["is_zero"] := true;
  for j from k to n do
    for i from k to m do
      # Look for a non-zero pivot
      Mij["value"]    := M[i, j];
      Mij["i"]        := i;
      Mij["j"]        := j;
      Mij["degree_r"] := M_degree_R[i, j];
      Mij["degree_c"] := M_degree_C[i, j];
      Mij["cost"], Mij["numeric_value"] := LULEM:-PivotCost(Mij);
      try
        Mij["value"]   := Normalizer(Mij["value"]);
        Mij["is_zero"] := evalb(Mij["value"] = 0);
        if not Mij["is_zero"] then
          uMij := apply_unveil(Mij["value"]);
          # Recalculate cost and value of the pivot
          Mij["cost"], Mij["numeric_value"] := LULEM:-PivotCost(uMij);
          # Time limit required because sometimes Normalizer get stuck
          uMij := timelimit(TimeLimit, eval(Normalizer(uMij)));
          Mij["is_zero"] := evalb(uMij = 0);
        end if;
      catch "time expired":
        printf("LULEM::Pivoting(...): simplify(Mij) failed, assumed <> 0.\n");
        Mij["is_zero"] := false;
      catch:
        printf("LULEM::Pivoting(...): Mij division by 0 or other exception.\n");
        if LULEM:-Verbose then
          print(eval(Mij));
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
        elif LULEM:-PivotStrategy_type(pivot, Mij) then
          # A better pivot is found
          pivot := copy(Mij);
        end if;
      end if;
    end do;
  end do;

  if not pivot["is_zero"] then
    i := pivot["i"];
    j := pivot["j"];
    if (i <> k) then
      (r[i], r[k])    := (r[k], r[i]);
      M[[i,k], 1..-1] := M[[k,i], 1..-1];
    end if;
    if (j <> k) then
      (c[j], c[k])    := (c[k], c[j]);
      M[1..-1, [j,k]] := M[1..-1 ,[k,j]];
    end if;
  end if;

  return pivot;
end proc: # Pivoting

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotCost := proc(
  x::{algebraic},
  $)::{integer};

  description "Compute the cost of the pivot <x>.";

  if type(x, 'integer') or type(x, 'float') then
    if evalb(x = 0) then
      return 0, 0;
    else
      return 1, abs(x);
    end if;
  end if;
  return LULEM:-Cost(x), infinity;
end proc: # PivotCost

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SetPivotStrategy := proc(
  str::{string},
  $)::{nothing};
  if (str = "Row") then
    LULEM:-PivotStrategy_type := LULEM:-PivotingStrategy_Row;
  elif (str = "Col") then
    LULEM:-PivotStrategy_type := LULEM:-PivotingStrategy_Col;
  elif (str = "Sum") then
    LULEM:-PivotStrategy_type := LULEM:-PivotingStrategy_Sum;
  elif (str = "Prod") then
    LULEM:-PivotStrategy_type := LULEM:-PivotingStrategy_Prod;
  elif (str = "Min") then
    LULEM:-PivotStrategy_type := LULEM:-PivotingStrategy_Min;
  elif (str = "Val") then
    LULEM:-PivotStrategy_type := LULEM:-PivotingStrategy_Val;
  else
    error "unknown pivoting strategy %s.\n", str;
  end if;
  return NULL;
end proc: # SetVeilingStrategy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_Row := proc(
  cur::table,
  val::table,
  $)::{boolean};

  description "Compute the pivoting strategy 1: given the current pivot <cur> "
    "and the next pivot <new>, decide if to the new pivot is better than the "
    "current  pivot or not.";

  if (val["degree_r"] < cur["degree_r"]) then
    return true;
  elif (val["degree_r"] > cur["degree_r"]) then
    return false;
  elif (val["numeric_value"] < cur["numeric_value"]) then
    return true;
  elif (val["numeric_value"] > cur["numeric_value"]) then
    return false;
  elif (val["numeric_value"] = infinity) then
    # All equal and symbolic
    return evalb(val["cost"] < cur["cost"]);
  end if;
  # All equal and finite
  return evalb(val["numeric_value"] > cur["numeric_value"]);
end proc: # PivotingStrategy_Row

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_Col := proc(
  cur::table,
  val::table,
  $)::{boolean};

  description "Compute the pivoting strategy 1: given the current pivot <cur> "
    "and the next pivot <new>, decide if to the new pivot is better than the "
    "current  pivot or not.";

  if val["degree_c"] < cur["degree_c"] then
    return true;
  elif val["degree_c"] > cur["degree_c"] then
    return false
  elif val["numeric_value"] < cur["numeric_value"] then
    return true;
  elif val["numeric_value"] > cur["numeric_value"] then
    return false;
  elif val["numeric_value"] = infinity then
    # All equal and symbolic
    return val["cost"] < cur["cost"];
  end if;
  # All equal and finite
  return val["numeric_value"] > cur["numeric_value"];

end proc: # PivotingStrategy_Col

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_Sum := proc(
  cur::table,
  val::table,
  $)::{boolean};

  description "Compute the pivoting strategy 1: given the current pivot <cur> "
    "and the next pivot <new>, decide if to the new pivot is better than the "
    "current  pivot or not.";

  local dc, dv;

  dc := cur["degree_c"] + cur["degree_r"];
  dv := val["degree_c"] + val["degree_r"];

  if (dv < dc) then
    return true;
  elif (dc > dv) then
    return false
  elif (val["numeric_value"] < cur["numeric_value"]) then
    return true;
  elif (val["numeric_value"] > cur["numeric_value"]) then
    return false;
  elif (val["numeric_value"] = infinity) then
    # all equals and symbolic
    return evalb(val["cost"] < cur["cost"]);
  end if;
  # all equals and finite
  return evalb(val["numeric_value"] > cur["numeric_value"]);
end proc: # PivotingStrategy_Sum

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_Prod := proc(
  cur::table,
  val::table,
  $)::{boolean};

  description "Compute the pivoting strategy 2: given the current pivot <cur> "
    "and the next pivot <new>, decide if to the new pivot is better than the "
    "current  pivot or not.";

  local dc, dv;

  dc := cur["degree_c"] * cur["degree_r"];
  dv := val["degree_c"] * val["degree_r"];

  if (dv < dc) then
    return true;
  elif (dc > dv) then
    return false
  elif (val["numeric_value"] < cur["numeric_value"]) then
    return true;
  elif (val["numeric_value"] > cur["numeric_value"]) then
    return false;
  elif (val["numeric_value"] = infinity) then
    # All equal and symbolic
    return evalb(val["cost"] < cur["cost"]);
  end if;
  # All equal and finite
  return evalb(val["numeric_value"] > cur["numeric_value"]);
end proc: # PivotingStrategy_Prod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_Min := proc(
  cur::table,
  val::table,
  $)::{boolean};

  description "Compute the pivoting strategy 3: given the current pivot <cur> "
    "and the next pivot <new>, decide if to the new pivot is better than the "
    "current  pivot or not.";

  local dc, dv;

  dc := min(cur["degree_c"], cur["degree_r"]);
  dv := min(val["degree_c"], val["degree_r"]);

  if (dv < dc) then
    return true;
  elif (dc > dv) then
    return false
  elif (val["numeric_value"] < cur["numeric_value"]) then
    return true;
  elif (val["numeric_value"] > cur["numeric_value"]) then
    return false;
  elif (val["numeric_value"] = infinity) then
    # all equals and symbolic
    return evalb(val["cost"] < cur["cost"]);
  end if;
  # all equals and finite
  return evalb(val["numeric_value"] > cur["numeric_value"]);
end proc: # PivotingStrategy_Min

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_Val := proc(
  cur::table,
  val::table,
  $)::{boolean};

  description "Compute the pivoting strategy 4: given the current pivot <cur> "
    "and the next pivot <new>, decide if to the new pivot is better than the "
    "current  pivot or not.";

  if (val["numeric_value"] < cur["numeric_value"]) then
    return true;
  elif (val["numeric_value"] > cur["numeric_value"]) then
    return false;
  elif (val["numeric_value"] = infinity) then
    # All equal and symbolic
    return evalb(val["cost"] < cur["cost"]);
  end if;
  # All equal and finite
  return evalb(val["numeric_value"] > cur["numeric_value"]);
end proc: # PivotingStrategy_Val

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
