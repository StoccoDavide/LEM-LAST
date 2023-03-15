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

  local Mij, uMij, Mij_is_zero, Mij_cost, Mij_value, Mij_degree, M_degree,
        pivot_is_zero, pivot_cost, pivot_value, pivot_degree,
        m, n, i, j, ii, jj, apply_unveil, z;

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
          pivot_is_zero := false;
          pivot_cost    := Mij_cost;
          pivot_value   := Mij_value;
          pivot_degree  := Mij_degree;
          i             := ii;
          j             := jj;
        elif (Mij_degree < pivot_degree) or
             (Mij_degree = pivot_degree and Mij_cost < pivot_cost) or
             (Mij_degree = pivot_degree and Mij_cost = pivot_cost and Mij_value > pivot_value) then
          # A better pivot is found
          pivot_cost   := Mij_cost;
          pivot_value  := Mij_value;
          pivot_degree := Mij_degree;
          i            := ii;
          j            := jj;
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
