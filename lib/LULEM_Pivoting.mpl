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
  $)

  description "Compute the LU decomposition pivot vector provided the step <k>, "
              "the temporary LU (NAG) matrix <M>, the veiling symbol <V>, "
              "the rows the pivot vector <r>, the columns the pivot vector <c> "
              "and the veiling strategy <VeilingStrategy>.";

  local Mij, uMij, Mij_is_zero, Mij_cost, Mij_value,
        pivot_is_zero, pivot_cost, pivot_value,
        LV, m, n, i, j, ii, jj, n_found,
        apply_unveil, z, tmp;

  # Extract dimensions
  m, n := LinearAlgebra[Dimensions](M):

  # Check if to veil or not to veil
  LV           := LEM:-VeilList(V,true);
  apply_unveil := (z) -> subs[eval](op( LV ), x );

  pivot_is_zero := true;
  pivot_cost    := 0;
  pivot_value   := 0;

  # Iterate over the columns and rows
  n_found := 0;
  for jj from k to n do
    for ii from k to m do
      # Look for a non-zero pivot
      Mij                 := M[ii,jj];
      Mij_cost, Mij_value := LULEM:-PivotCost(Mij);
      try
        Mij         := Normalizer(Mij);
        Mij_is_zero := evalb( Mij = 0 );
        if not Mij_is_zero then
          uMij                := apply_unveil(Mij);
          Mij_cost, Mij_value := LULEM:-PivotCost(uMij); # recalculate
          # timelimit required because sometime Normalizer stuck
          uMij        := timelimit( 0.5, eval(Normalizer(uMij)) );
          Mij_is_zero := evalb( uMij = 0 );
        end;
      catch:
        print("Mij: Division by 0 or another exception.\n");
        if LULEM:-Verbose then
          #print(lastexception);
          print(Mij);
        end if;
        Mij_is_zero := true;
      end try;

      if Mij_is_zero then
        M[ii,jj] := 0;
      else
        # found non zero pivot, check if it is better
        if pivot_is_zero then
          # fist non zero pivot!
          pivot_is_zero := false;
          pivot_cost    := Mij_cost;
          pivot_value   := Mij_value;
          i             := ii;
          j             := jj;
        #elif Mij_cost = pivot_cost then
        elif Mij_cost < pivot_cost or ( Mij_cost = pivot_cost and Mij_value > pivot_value ) then
          # A better pivot is found
          pivot_cost  := Mij_cost;
          pivot_value := Mij_value;
          i          := ii;
          j          := jj;
        end
      end if;
    end do;
    # check if found pivot
    #if not pivot_is_zero then
    #  if pivot_cost = 0 then
    #    break;
    #  end;
    #  if pivot_cost = 1 and n_found > 0 then
    #    break;
    #  end;
    #  if pivot_cost > 1 and n_found > 1 then
    #    break;
    #  end;
    #  n_found := n_found+1;
    #end if;
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

  return pivot_is_zero, M[k,k], pivot_cost;

end proc: # DoPivoting

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PivotCost := proc( x::{algebraic}, $ )::{integer};
  if type(x,'integer') or type(x,'float') then
    if evalb( x = 0 ) then
      return 0, 0;
    else
      return 1, abs(x);
    end if;
  elif type(x,'symbol') then
    return 2, infinity;
  #elif type(x,'algebraic') then
  #  return infinity, 1;
  end if;
  return 2+length(x), infinity;
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
