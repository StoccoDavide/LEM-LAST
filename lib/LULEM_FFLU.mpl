# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   _____ _____ _    _   _
#  |  ___|  ___| |  | | | |
#  | |_  | |_  | |  | | | |
#  |  _| |  _| | |__| |_| |
#  |_|   |_|   |_____\___/
#
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FFLU := proc(
  A::{Matrix},
  V::{symbol},
  VeilingStrategy::{procedure} := VeilingStrategy_n,
  $)::{table};

  description "Compute the Fracton-Free LU decomposition of a square matrix "
              "<A> using the veiling strategy <VeilingStrategy> and the veiling symbol <V>.";

  local SS, M, Mkk, m, n, mn, i, j, k, ri, rk, rnk, r, c, apply_veil,
        pivot_is_zero, Mij_is_zero, z, tmp, bot, top;

  m, n := LinearAlgebra[Dimensions](A):

  # Create pivot vector
  r := Vector(m, k -> k);
  c := Vector(n, k -> k);

  # check if Veil or not
  apply_veil := z -> `if`(VeilingStrategy(z), LEM:-Veil[V](z), z);

  # Gauss Elimination main loop
  M   := copy(A);
  mn  := min(m, n);
  rnk := mn;
  SS  := Vector(mn);
  for k from 1 to mn-1 do
    if LULEM:-Verbose then
      printf(
        "LULEM::FFLU(...): processing %d-th row. Length %d:%d\n",
        k, length(convert(M,list)), length(LEM:-VeilList(V))
      );
    end;

    pivot_is_zero, Mkk := DoPivoting( k, M, V, r, c, VeilingStrategy );

    if pivot_is_zero then
      rnk := k;
      if LULEM:-Verbose then
        WARNING("LULEM::LU(...): the matrix appears not full rank.");
      end;
      break;
    end if;
    if LULEM:-Verbose then
      print("LULEM::FFLU(...): pivot:", Mkk);
    end;
    top   := apply_veil(Normalizer(numer(Mkk)));
    bot   := apply_veil(Normalizer(denom(Mkk)));
    SS[k] := top;

    # Scaled Shur complement
    tmp        := [k+1..-1];
    M[tmp, k]  := apply_veil~(Normalizer~(M[tmp, k]*bot));
    M[tmp,tmp] := apply_veil~(Normalizer~( top*M[tmp,tmp] - M[tmp,k].M[k,tmp]) );
  end do:

  # Return the FFLU decomposition
  return table([
    "method"   = "FFLU",
    "M"        = M,
    "V"        = V,
    "S"        = SS,
    "r"        = r,
    "c"        = c,
    "rank"     = rnk,
    "M_length" = length(convert(M,list)),
    "S_length" = length(convert(SS,list)),
    "V_length" = length(LEM:-VeilList(V))
  ]);
end proc: # FFLU

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
FF2LU := proc( T::{table}, $ )
  local M, SS, r, c, rk, n, m, L, DG, L_list, D_list, U, i, j, k;

  M  := T["M"];
  SS := T["S"];
  r  := T["r"];
  c  := T["c"];
  rk := T["rank"];
  m, n := LinearAlgebra[Dimensions](M):
  L_list := [];
  D_list := [];
  #
  #   +             + +         +      +            +       +             + +            +
  #   |  1          | | 1       |      | 1          |       | 1           | | 1          |
  #   | -x  1       | |   q     | ==>  | -x q       |  ==>  |   1/q       | | x  1       |
  #   | -x     1    | |     q   |      | -x    q    |       |      1/q    | | x     1    |
  #   | -x        1 | |       q |      | -x       q |       |         1/q | | x        1 |
  #   +             + +         +      +            +       +             + +            +
  #
  for k from 1 to rk do
    L            := Matrix( m, m, shape = triangular[lower, unit]);
    L[k+1..-1,k] := M[k+1..-1,k];
    L_list       := [ op(L_list), L ];
    DG           := Matrix( LinearAlgebra[IdentityMatrix](m), shape = diagonal );
    for j from k+1 to m do
      DG[j,j] := SS[k];
    end;
    D_list := [ op(D_list), DG ];
  end;
  return L_list, D_list, Matrix(M, shape = triangular[upper]);
end proc: # FF2LU

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SolveFractionalFreeLU := proc(
  T::{table},
  b::{Vector},
  V::{symbol, function},
  VeilingStrategy::{procedure} := VeilingStrategy_n,
  $)

  description "Solve the linear system Ax=b using FFLU decomposition <T>, "
              "provided the vector <b>, the veiling symbol <V> and the "
              "veiling strategy <VeilingStrategy>.";

  local m, n, M, S, r, c, rk, x, i, s,  apply_veil;

  M  := T["M"];
  S  := T["S"];
  r  := T["r"];
  c  := T["c"];
  rk := T["rank"];

  # Get linear system dimension
  m, n := LinearAlgebra[Dimensions](M);

  # Check if the linear system is consistent
  assert( m = n, "LULEM::SolveFractionalFreeLU(...): rectangular linear system." );

  # Create a normalizer function
  apply_veil := (y) -> `if`(VeilingStrategy(y), LEM:-Veil[V](y), y);

  # apply permutation P
  x := b[convert(r,list)];

  # Perform forward substitution to solve Ly=b[r]
  for i from 1 to rk-1 do
    x[i+1..-1] := S[i] * x[i+1..-1]; # apply D
    x[i+1..-1] := apply_veil~( x[i+1..-1] - x[i]*M[i+1..-1,i]);
  end do;

  # Perform backward substitution to solve Ux[c]=y
  print(M[n-1..n,n-1..n]);
  x[n] := apply_veil(x[n]/M[n,n]);
  for i from n-1 to 1 by -1 do
    s    := apply_veil( x[i] - add(M[i,i+1..n] *~ x[i+1..n]) );
    x[i] := apply_veil( s / M[i,i] );
  end do;

  # apply permutation
  x := x[convert(c,list)];
  # Return outputs
  return x;
end proc: # SolveFractionalFreeLU

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
