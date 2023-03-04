# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   ____  _            _   _
#  |  _ \(_)_   _____ | |_(_)_ __   __ _
#  | |_) | \ \ / / _ \| __| | '_ \ / _` |
#  |  __/| |\ V / (_) | |_| | | | | (_| |
#  |_|   |_| \_/ \___/ \__|_|_| |_|\__, |
#                                  |___/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

DoPivoting := proc(
  k::{integer},
  M::{Matrix},
  V::{symbol},
  r::{Vector(integer)},
  c::{Vector(integer)},
  VeilingStrategy::{procedure},
  $)

  description "Compute the LU decomposition pivot vector provided the step <k>, "
              "the temporary LU (NAG) matrix <M>, the veiling symbol <V>, "
              "the rows the pivot vector <r>, the columns the pivot vector <c> "
              "and the veiling strategy <VeilingStrategy>.";

  local Mij, uMij, Mij_is_zero, Mij_cost, pivot_is_zero, pivot_cost,
        LV, m, n, i, j, ii, jj,
        apply_veil, apply_unveil, z, tmp;

  # Extract dimensions
  m, n := LinearAlgebra[Dimensions](M):

  # Check if to veil or not to veil
  apply_veil   := (z) -> `if`( VeilingStrategy(z), LEM:-Veil[V](z), z);
  LV           := LEM:-VeilList(V,true);
  apply_unveil := (z) -> subs( op( LV ), z );

  pivot_is_zero := true;
  pivot_cost    := 0;

  # Iterate over the columns and rows
  #for jj from k to n do
  #  for ii from k to m do
  for ii from k to m do
    for jj from k to n do
      # Look for a non-zero pivot
      Mij      := M[ii,jj];
      Mij_cost := LULEM:-PivotCost(Mij);
      try
        Mij         := Normalizer(Mij);
        Mij_cost    := LULEM:-PivotCost(Mij); # recalculate, might be better
        Mij_is_zero := evalb( Mij = 0 );
        if not Mij_is_zero then
          uMij := apply_unveil(Mij);
          # timelimit required because sometime Normalizer stuck
          uMij        := timelimit( 0.5, eval(Normalizer(eval(uMij))) );
          Mij_is_zero := evalb( uMij = 0 );
        end;
      catch:
        if LULEM:-Verbose then
          print("Mij: Division by 0 or numerical exception.\n");
          print(Mij);
        end if;
        Mij_is_zero := true;
      end try;

      if not Mij_is_zero then
        # found non zero pivot, check if it is better
        if pivot_is_zero then
          # fist non zero pivot!
          pivot_is_zero := false;
          pivot_cost    := Mij_cost;
          i             := ii;
          j             := jj;
        #elif Mij_cost = pivot_cost then
        elif Mij_cost < pivot_cost then
          # A better pivot is found
          pivot_cost := Mij_cost;
          i          := ii;
          j          := jj;
        end
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

  return pivot_is_zero, M[k,k];

end proc: # DoPivoting

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotCost := proc( x::{algebraic}, $ )::{integer};
  return nops(indets(x));
end proc: # PivotCost

(*
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_Llength := proc( x::{algebraic}, y::{algebraic}, $ )::{boolean};
  description "Pivoting strategy: choose the pivot with the largest length "
              "between expressions <x> and <y>.";
  return evalb( length(x) > length(y) );
end proc: # PivotingStrategy_Llength

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_Slength := proc( x::{algebraic}, y::{algebraic}, $ )::{boolean};
  description "Pivoting strategy: choose the pivot with the smallest length "
              "between expressions <x> and <y>.";
  return evalb( nops(indets(x)) < nops(indets(y)) );
end proc: # PivotingStrategy_Slength

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_Lindets := proc ( x::{algebraic}, y::{algebraic}, $ )::{boolean};
  description "Pivoting strategy: choose the pivot with the largest number of "
              "indeterminates between expressions <x> and <y>.";
  return evalb( nops(indets(x)) > nops(indets(y)) );
end proc: # PivotingStrategy_Lindets

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_Sindets := proc ( x::{algebraic}, y::{algebraic}, $ )::{boolean};
  description "Pivoting strategy: choose the pivot with the smallest number of "
              "indeterminates between expressions <x> and <y>.";
  return evalb( nops(indets(x)) < nops(indets(y)) );
end proc: # PivotingStrategy_Sindets

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotingStrategy_numeric := proc ( x::{algebraic}, y::{algebraic}, $ )::{boolean};
  description "Pivoting strategy: choose the pivot with the largest numeric "
              "value between expressions <x> and <y>.";
  return evalb( evalf(abs(x)) > evalf(abs(y)) );
end proc: # PivotingStrategy_numeric
*)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
