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
  $)::{table};

  description "Compute the LU decomposition pivots vectors with minum degree "
    "provided the step <k>, the temporary LU (NAG) matrix <M>, the veiling "
    "symbol <V>, the rows the pivot vector <r>, the columns the pivot vector "
    "<c>.";

  local uMij, M_degree_R, M_degree_C, perm_R, perm_C, m, n, i, j, apply_unveil,
    z, Mij, pivot;

  # Extract matrix dimensions
  m, n := LinearAlgebra:-Dimensions(M):

  # Check if to veil or not to veil
  apply_unveil := (z) -> LEM:-UnVeil[V](z);

  # Calculate the degree
  M_degree_R := Matrix(m, n);
  M_degree_C := Matrix(m, n);
  M_degree_R[k..m, k..n], M_degree_C[k..m, k..n] := LULEM:-GetDegrees(M[k..m, k..n]);
  perm_R := (k-1) +~ sort(M_degree_R[k..m, k], `<`, output = 'permutation');
  perm_C := (k-1) +~ sort(M_degree_C[k, k..n], `<`, output = 'permutation');

  pivot := table([]);
  Mij   := table([]);

  # Iterate over the columns and rows
  pivot["is_zero"] := true;
  for j in perm_C do
    for i in perm_R do

      # Look for a non-zero pivot
      Mij["value"]    := M[i, j];
      Mij["i"]        := i;
      Mij["j"]        := j;
      Mij["degree_r"] := M_degree_R[i, j];
      Mij["degree_c"] := M_degree_C[i, j];
      Mij["cost"], Mij["numeric_value"] := LULEM:-PivotCost(Mij);

      # Pre-check if the pivot is better than the previous one
      if not pivot["is_zero"] and not LULEM:-MinDegreeStrategy_fun(pivot, Mij) then
        break;
      end if;

      # Try to simplify the pivot expression
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
          print(Mij["value"]);
        end if;
        Mij["is_zero"] := true;
      end try;

      if Mij["is_zero"] then
        M[i, j] := 0;
      else
        # Found a non-zero pivot, check if it is better
        if pivot["is_zero"] then
          # First non-zero pivot found
          pivot := copy(Mij);
        elif LULEM:-PivotingStrategy_fun(pivot, Mij) then
          # A better pivot is found
          pivot := copy(Mij);
        end if;
      end if;
    end do;
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
#   __  __ _       ____
#  |  \/  (_)_ __ |  _ \  ___  __ _ _ __ ___  ___
#  | |\/| | | '_ \| | | |/ _ \/ _` | '__/ _ \/ _ \
#  | |  | | | | | | |_| |  __/ (_| | | |  __/  __/
#  |_|  |_|_|_| |_|____/ \___|\__, |_|  \___|\___|
#                             |___/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SetMinDegreeStrategy := proc(
  str::{string},
  $)::{nothing};
  if (str = "none") then
    LULEM:-MinDegreeStrategy_fun := LULEM:-MinDegreeStrategy_none;
  elif (str = "row") then
    LULEM:-MinDegreeStrategy_fun := LULEM:-MinDegreeStrategy_row;
  elif (str = "col") then
    LULEM:-MinDegreeStrategy_fun := LULEM:-MinDegreeStrategy_col;
  elif (str = "sum") then
    LULEM:-MinDegreeStrategy_fun := LULEM:-MinDegreeStrategy_sum;
  elif (str = "prod") then
    LULEM:-MinDegreeStrategy_fun := LULEM:-MinDegreeStrategy_prod;
  elif (str = "min") then
    LULEM:-MinDegreeStrategy_fun := LULEM:-MinDegreeStrategy_min;
  else
    error "unknown minimum degree strategy %s.\n", str;
  end if;
  return NULL;
end proc: # SetMinDegreeStrategy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MinDegreeStrategy_none := proc(
  cur::{table},
  val::{table},
  $)::{boolean};

  description "Compute the pivoting precheck: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current  pivot or not.";

  return true;
end proc: # MinDegreeStrategy_row

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MinDegreeStrategy_row := proc(
  cur::{table},
  val::{table},
  $)::{boolean};

  description "Compute the pivoting precheck: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current  pivot or not.";

  if (val["degree_r"] > cur["degree_r"]) then
    return false;
  else
    return true;
  end if;
end proc: # MinDegreeStrategy_row

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MinDegreeStrategy_col := proc(
  cur::{table},
  val::{table},
  $)::{boolean};

  description "Compute the pivoting precheck: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current  pivot or not.";

  if (val["degree_c"] > cur["degree_c"]) then
    return false;
  else
    return true;
  end if;
end proc: # MinDegreeStrategy_col

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MinDegreeStrategy_sum := proc(
  cur::{table},
  val::{table},
  $)::{boolean};

  description "Compute the pivoting precheck: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current  pivot or not.";

  local degree_cur, degree_val;

  degree_cur := cur["degree_c"] + cur["degree_r"];
  degree_val := val["degree_c"] + val["degree_r"];

  if (degree_val > degree_cur) then
    return false;
  else
    return true;
  end if;
end proc: # MinDegreeStrategy_sum

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MinDegreeStrategy_prod := proc(
  cur::{table},
  val::{table},
  $)::{boolean};

  description "Compute the pivoting precheck: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current  pivot or not.";

  local degree_cur, degree_val;

  degree_cur := cur["degree_c"] * cur["degree_r"];
  degree_val := val["degree_c"] * val["degree_r"];

  if (degree_val > degree_cur) then
    return false;
  else
    return true;
  end if;
end proc: # MinDegreeStrategy_prod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MinDegreeStrategy_min := proc(
  cur::{table},
  val::{table},
  $)::{boolean};

  description "Compute the pivoting precheck: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current  pivot or not.";

  local degree_cur, degree_val;

  degree_cur := min(cur["degree_c"], cur["degree_r"]);
  degree_val := min(val["degree_c"], val["degree_r"]);

  if (degree_val > degree_cur) then
    return false;
  else
    return true;
  end if;
end proc: # MinDegreeStrategy_min

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   ____  _            _   _
#  |  _ \(_)_   _____ | |_(_)_ __   __ _
#  | |_) | \ \ / / _ \| __| | '_ \ / _` |
#  |  __/| |\ V / (_) | |_| | | | | (_| |
#  |_|   |_| \_/ \___/ \__|_|_| |_|\__, |
#                                  |___/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SetPivotingStrategy := proc(
  str::{string},
  $)::{nothing};
  if (str = "val") then
    LULEM:-PivotingStrategy_fun := LULEM:-PivotingStrategy_val;
  elif (str = "row") then
    LULEM:-PivotingStrategy_fun := LULEM:-PivotingStrategy_row;
  elif (str = "col") then
    LULEM:-PivotingStrategy_fun := LULEM:-PivotingStrategy_col;
  elif (str = "sum") then
    LULEM:-PivotingStrategy_fun := LULEM:-PivotingStrategy_sum;
  elif (str = "prod") then
    LULEM:-PivotingStrategy_fun := LULEM:-PivotingStrategy_prod;
  elif (str = "min") then
    LULEM:-PivotingStrategy_fun := LULEM:-PivotingStrategy_min;
  else
    error "unknown pivoting strategy %s.\n", str;
  end if;
  return NULL;
end proc: # SetPivotingStrategy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_val := proc(
  cur::{table},
  val::{table},
  $)::{boolean};

  description "Compute the pivoting strategy: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current  pivot or not.";

  if (val["numeric_value"] < cur["numeric_value"]) then
    return true;
  elif (val["numeric_value"] > cur["numeric_value"]) then
    return false;
  elif (val["numeric_value"] = infinity) then
    return evalb(val["cost"] < cur["cost"]);
  else
    return evalb(val["numeric_value"] > cur["numeric_value"]);
  end if;
end proc: # PivotingStrategy_val

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_row := proc(
  cur::{table},
  val::{table},
  $)::{boolean};

  description "Compute the pivoting strategy: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current  pivot or not.";

  if (val["degree_r"] < cur["degree_r"]) then
    return true;
  elif (val["degree_r"] > cur["degree_r"]) then
    return false;
  else
    return LULEM:-PivotingStrategy_val(cur, val);
  end if;
end proc: # PivotingStrategy_row

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_col := proc(
  cur::{table},
  val::{table},
  $)::{boolean};

  description "Compute the pivoting strategy: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current  pivot or not.";

  if val["degree_c"] < cur["degree_c"] then
    return true;
  elif val["degree_c"] > cur["degree_c"] then
    return false;
  else
    return LULEM:-PivotingStrategy_val(cur, val);
  end if;
end proc: # PivotingStrategy_col

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_sum := proc(
  cur::{table},
  val::{table},
  $)::{boolean};

  description "Compute the pivoting strategy: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current  pivot or not.";

  local degree_cur, degree_val;

  degree_cur := cur["degree_c"] + cur["degree_r"];
  degree_val := val["degree_c"] + val["degree_r"];

  if (degree_val < degree_cur) then
    return true;
  elif (degree_cur > degree_val) then
    return false;
  else
    return LULEM:-PivotingStrategy_val(cur, val);
  end if;
end proc: # PivotingStrategy_sum

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_prod := proc(
  cur::{table},
  val::{table},
  $)::{boolean};

  description "Compute the pivoting strategy: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current  pivot or not.";

  local degree_cur, degree_val;

  degree_cur := cur["degree_c"] * cur["degree_r"];
  degree_val := val["degree_c"] * val["degree_r"];

  if (degree_val < degree_cur) then
    return true;
  elif (degree_cur > degree_val) then
    return false;
  else
    return LULEM:-PivotingStrategy_val(cur, val);
  end if;
end proc: # PivotingStrategy_prod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_min := proc(
  cur::{table},
  val::{table},
  $)::{boolean};

  description "Compute the pivoting strategy: given the current pivot <cur> "
    "and the next pivot <val>, decide if to the next pivot is better than the "
    "current  pivot or not.";

  local degree_cur, degree_val;

  degree_cur := min(cur["degree_c"], cur["degree_r"]);
  degree_val := min(val["degree_c"], val["degree_r"]);

  if (degree_val < degree_cur) then
    return true;
  elif (degree_cur > degree_val) then
    return false;
  else
    return LULEM:-PivotingStrategy_val(cur, val);
  end if;
end proc: # PivotingStrategy_min

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
