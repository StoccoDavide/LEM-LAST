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
  $)

  description "Compute the LU decomposition pivots vectors provided the step "
    "<k>, the temporary LU (NAG) matrix <M>, the veiling symbol <V>, the rows "
    "the pivot vector <r>, the columns the pivot vector <c>.";

  local Mij, uMij, Mij_is_zero, M_degree, m, n, i, j, ii, jj, apply_unveil, z,
        pivot_is_zero, Mij_cost, Mij_value, Mij_degree, Mij_degree_r, Mij_degree_c,
        pivot_cost, pivot_value, pivot_degree, pivot_degree_r, pivot_degree_c;

  # Extract matrix dimensions
  m, n := LinearAlgebra:-Dimensions(M):

  # Check if to veil or not to veil
  apply_unveil := (z) -> LEM:-UnVeil[V](z);

  # Calculate the degree
  M_degree             := Matrix(m, n);
  M_degree[k..m, k..n] := LULEM:-GetDegree(M[k..m, k..n]);

  # Iterate over the columns and rows
  pivot_is_zero := true;
  pivot_cost    := 0;
  pivot_value   := 0;
  for jj from k to n do
    for ii from k to m do
      # Look for a non-zero pivot
      Mij                 := M[ii, jj];
      Mij_degree          := M_degree[ii, jj];
      Mij_degree_r        := convert(M_degree[ii, 1..-1], `+`);
      Mij_degree_c        := convert(M_degree[1..-1, jj], `+`);
      Mij_cost, Mij_value := LULEM:-PivotCost(Mij);
      try
        Mij         := Normalizer(Mij);
        Mij_is_zero := evalb(Mij = 0);
        if not Mij_is_zero then
          uMij := apply_unveil(Mij);
          # Recalculate cost and value of the pivot
          Mij_cost, Mij_value := LULEM:-PivotCost(uMij);
          # Time limit required because sometimes Normalizer get stuck
          uMij        := timelimit(TimeLimit, eval(Normalizer(uMij)));
          Mij_is_zero := evalb(uMij = 0);
        end if;
      catch "time expired":
        printf("LULEM::Pivoting(...): simplify(Mij) failed, assumed <> 0.\n");
        Mij_is_zero := false;
      catch:
        printf("LULEM::Pivoting(...): Mij division by 0 or other exception.\n");
        if LULEM:-Verbose then
          print(Mij);
        end if;
        Mij_is_zero := true;
      end try;

      if Mij_is_zero then
        M[ii, jj] := 0;
      else
        # Found a non-zero pivot, check if it is better
        if pivot_is_zero then
          # First non-zero pivot found
          pivot_is_zero  := false;
          pivot_cost     := Mij_cost;
          pivot_value    := Mij_value;
          pivot_degree   := Mij_degree;
          pivot_degree_r := Mij_degree_r;
          pivot_degree_c := Mij_degree_c;
          i              := ii;
          j              := jj;
        elif PivotingStrategy_1r(
            pivot_value, pivot_degree, pivot_degree_r, pivot_degree_c, pivot_cost,
            Mij_value,   Mij_degree,   Mij_degree_r,   Mij_degree_c,   Mij_cost
          ) then
          # A better pivot is found
          pivot_cost     := Mij_cost;
          pivot_value    := Mij_value;
          pivot_degree   := Mij_degree;
          pivot_degree_r := Mij_degree_r;
          pivot_degree_c := Mij_degree_c;
          i              := ii;
          j              := jj;
        end if;
      end if;
    end do;
  end do;

  if not pivot_is_zero then
    if (i <> k) then
      (r[i], r[k])    := (r[k], r[i]);
      M[[i,k], 1..-1] := M[[k,i], 1..-1];
    end if;
    if (j <> k) then
      (c[j], c[k])    := (c[k], c[j]);
      M[1..-1, [j,k]] := M[1..-1 ,[k,j]];
    end if;
  end if;

  return pivot_is_zero, M[k, k], pivot_cost, pivot_degree;
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
  elif type(x, 'symbol') then
    return 2, infinity;
  #elif type(x, 'algebraic') then
  #  return infinity, 1;
  end if;
  return 2 + length(x), infinity;
end proc: # PivotCost

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_1r := proc(
  cur_value::{algebraic, numeric},
  cur_degree::{nonnegint},
  cur_degree_r::{nonnegint},
  cur_degree_c::{nonnegint},
  cur_cost::{nonnegint},
  new_value::{algebraic, numeric},
  new_degree::{nonnegint},
  new_degree_r::{nonnegint},
  new_degree_c::{nonnegint},
  new_cost::{nonnegint},
  $)::{boolean};

  description "Compute the pivoting strategy 1: given the current pivot <cur> "
    "and the next pivot <new>, decide if to the new pivot is better than the "
    "current  pivot or not.";

  return
    (new_degree < cur_degree) or
    (new_degree = cur_degree and new_degree_r < cur_degree_r) or
    (new_degree = cur_degree and new_cost < cur_cost) or
    (new_degree = cur_degree and new_cost = cur_cost and new_value > cur_value);

end proc: # PivotingStrategy_1r

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_1c := proc(
  cur_value::{algebraic, numeric},
  cur_degree::{nonnegint},
  cur_degree_r::{nonnegint},
  cur_degree_c::{nonnegint},
  cur_cost::{nonnegint},
  new_value::{algebraic, numeric},
  new_degree::{nonnegint},
  new_degree_r::{nonnegint},
  new_degree_c::{nonnegint},
  new_cost::{nonnegint},
  $)::{boolean};

  description "Compute the pivoting strategy 1: given the current pivot <cur> "
    "and the next pivot <new>, decide if to the new pivot is better than the "
    "current  pivot or not.";

  return
    (new_degree < cur_degree) or
    (new_degree = cur_degree and new_degree_c < cur_degree_c) or
    (new_degree = cur_degree and new_cost < cur_cost) or
    (new_degree = cur_degree and new_cost = cur_cost and new_value > cur_value);

end proc: # PivotingStrategy_1c

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_1rc := proc(
  cur_value::{algebraic, numeric},
  cur_degree::{nonnegint},
  cur_degree_r::{nonnegint},
  cur_degree_c::{nonnegint},
  cur_cost::{nonnegint},
  new_value::{algebraic, numeric},
  new_degree::{nonnegint},
  new_degree_r::{nonnegint},
  new_degree_c::{nonnegint},
  new_cost::{nonnegint},
  $)::{boolean};

  description "Compute the pivoting strategy 1: given the current pivot <cur> "
    "and the next pivot <new>, decide if to the new pivot is better than the "
    "current  pivot or not.";

  return
    (new_degree < cur_degree) or
    (new_degree = cur_degree and new_degree_r + new_degree_c < cur_degree_r + cur_degree_c) or
    (new_degree = cur_degree and new_cost < cur_cost) or
    (new_degree = cur_degree and new_cost = cur_cost and new_value > cur_value);

end proc: # PivotingStrategy_1rc

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_2 := proc(
  cur_value::{algebraic, numeric},
  cur_degree::{nonnegint},
  cur_degree_r::{nonnegint},
  cur_degree_c::{nonnegint},
  cur_cost::{nonnegint},
  new_value::{algebraic, numeric},
  new_degree::{nonnegint},
  new_degree_r::{nonnegint},
  new_degree_c::{nonnegint},
  new_cost::{nonnegint},
  $)::{boolean};

  description "Compute the pivoting strategy 2: given the current pivot <cur> "
    "and the next pivot <new>, decide if to the new pivot is better than the "
    "current  pivot or not.";

  return
    (new_degree < cur_degree) or
    (new_degree = cur_degree and new_cost < cur_cost) or
    (new_degree = cur_degree and new_cost = cur_cost and new_value > cur_value);

end proc: # PivotingStrategy_2

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_3 := proc(
  cur_value::{algebraic, numeric},
  cur_degree::{nonnegint},
  cur_degree_r::{nonnegint},
  cur_degree_c::{nonnegint},
  cur_cost::{nonnegint},
  new_value::{algebraic, numeric},
  new_degree::{nonnegint},
  new_degree_r::{nonnegint},
  new_degree_c::{nonnegint},
  new_cost::{nonnegint},
  $)::{boolean};

  description "Compute the pivoting strategy 3: given the current pivot <cur> "
    "and the next pivot <new>, decide if to the new pivot is better than the "
    "current  pivot or not.";

  return
    (new_cost < cur_cost) or
    (new_cost = cur_cost and new_value > cur_value);

end proc: # PivotingStrategy_3

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_4 := proc(
  cur_value::{algebraic, numeric},
  cur_degree::{nonnegint},
  cur_degree_r::{nonnegint},
  cur_degree_c::{nonnegint},
  cur_cost::{nonnegint},
  new_value::{algebraic, numeric},
  new_degree::{nonnegint},
  new_degree_r::{nonnegint},
  new_degree_c::{nonnegint},
  new_cost::{nonnegint},
  $)::{boolean};

  description "Compute the pivoting strategy 4: given the current pivot <cur> "
    "and the next pivot <new>, decide if to the new pivot is better than the "
    "current  pivot or not.";

  return
    evalb(new_value > cur_value);

end proc: # PivotingStrategy_4

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
