# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   _    _   _
#  | |  | | | |
#  | |  | | | |
#  | |__| |_| |
#  |_____\___/
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LU := proc(
  A::{Matrix},
  V::{symbol},
  $)::{table};

  description "Compute the LU decomposition of a square matrix <A> using the "
    "veiling strategy <VeilingStrategy> and the veiling symbol <V>.";

  local M, L, U, pivot, pivot_list, m, n, mn, k, rnk, r, c, apply_veil, tmp;

  # Forget the veilings
  LEM:-VeilForget(V);

  # Get matrix dimensions
  m, n := LinearAlgebra[Dimensions](A):

  # Create pivot vector
  r := Vector(m, k -> k);
  c := Vector(n, k -> k);

  # Check if to veil or not
  apply_veil := (z) -> `if`(LULEM:-VeilingStrategy(z), LEM:-Veil[V](z), z);

  # Perform Gaussian elimination
  M          := copy(A);
  mn         := min(m, n);
  rnk        := mn;
  pivot_list := [];
  for k from 1 to mn do
    if LULEM:-Verbose then
      printf(
        "LULEM::LU(...): processing %d-th row. Length %d: %d.\n",
        k, length(convert(M,list)), length(LEM:-VeilList(V))
      );
    end if;

    pivot      := LULEM:-Pivoting( k, M, V, r, c );
    pivot_list := [op(pivot_list), pivot["value"]];

    if pivot["is_zero"] then
      rnk := k;
      if LULEM:-Verbose then
        WARNING("LULEM::LU(...): the matrix appears to be not full rank.");
      end if;
      break;
    end if;

    if LULEM:-Verbose then
      printf(
        "LULEM::LU(...): M[%d,%d] = %a, cost = %d, degree_r = %d, degree_c = %d.\n",
        k, k, pivot["value"], pivot["cost"], pivot["degree_r"], pivot["degree_c"]
      );
    end if;

    # Shur complement
    tmp         := [k+1..-1];
    M[k, k]     := apply_veil(pivot["value"]);
    M[tmp, k]   := apply_veil~(Normalizer~(M[tmp, k])) / pivot["value"];
    M[k, tmp]   := apply_veil~(Normalizer~(M[k, tmp]));
    M[tmp, tmp] := apply_veil~(Normalizer~(M[tmp, tmp] - M[tmp, k].M[k, tmp]));
  end do;

  L := Matrix(M[1..m, 1..m], shape = triangular[lower, unit]);
  U := Matrix(M, shape = triangular[upper]);

  # Return the LU decomposition
  return table([
    "method"   = "LU",
    "L"        = L,
    "U"        = U,
    "V"        = V,
    "r"        = r,
    "c"        = c,
    "rank"     = rnk,
    "pivots"   = pivot_list,
    "L_cost"   = LULEM:-Cost(L),
    "U_cost"   = LULEM:-Cost(U),
    "V_cost"   = LULEM:-Cost(LEM:-VeilList(V)),
    "L_nnz"    = nops(op(2,L)),
    "U_nnz"    = nops(op(2,U)),
    "A_nnz"    = nops(op(2,A))
  ]);
end proc: # LU

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LUsolve := proc(
  T::{table},
  b::{Vector},
  V::{symbol, function},
  $)

  description "Solve the linear system Ax=b using LU decomposition <T>, "
    "provided the vector <b>, the veiling symbol <V> and the veiling strategy "
    "<VeilingStrategy>.";

  local L, U, r, c, m, n, p, q, apply_veil, x, y, i, j, s, rnk:

  # Extract the LU decomposition
  L   := T["L"];
  U   := T["U"];
  r   := T["r"];
  c   := T["c"];
  rnk := T["rank"];

  # Get linear system dimension
  m, n := LinearAlgebra[Dimensions](L);
  p, q := LinearAlgebra[Dimensions](U);

  # Check if the linear system is consistent
  assert(
    (m = n) and (p = q),
    "LULEM::LUsolve(...): only square system can be solved.\n"
    "L is %d x %d, U is %d x %d\n", m, n, p, q
  );

  # Check if the linear system is consistent
  assert(
    n = rnk,
    "LULEM::LUsolve(...): only full rank linear system can be solved.\n"
    "Rank is %d expected %d.\n", rnk, n
  );

  # Create a normalizer function
  apply_veil := (y) -> `if`(LULEM:-VeilingStrategy(y), LEM:-Veil[V](y), y);

  # apply permutation P
  x := b[convert(r, list)];

  # Perform forward substitution to solve Ly=b[r]
  for i from 2 to m do
    if LULEM:-Verbose then
      printf("LULEM::LUsolve(...): forward substitution of %d-th row.\n", i);
    end if;
    x[i] := apply_veil(x[i] - add(L[i, 1..i-1] *~ x[1..i-1]));
  end do;

  # Perform backward substitution to solve Ux[c]=y
  if LULEM:-Verbose then
    printf("LULEM::LUsolve(...): dividision by U[%d,%d].\n", n, n);
  end if;
  x[n] := apply_veil(x[n] / U[n, n]);
  for i from n-1 to 1 by -1 do
    if LULEM:-Verbose then
      printf("LULEM::LUsolve(...): backward substitution of %d-th column.\n", i);
    end if;
    s    := apply_veil(x[i] - add(U[i, i+1..n] *~ x[i+1..n]));
    x[i] := apply_veil(s / U[i, i]);
  end do;

  # Apply inverse permutation Q
  y := Vector[column](n);
  for i from 1 to n do
    y[c[i]] := x[i];
  end do;

  # Return outputs
  return y;
end proc: # LUsolve

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
