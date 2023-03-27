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
    "veiling symbol <V>.";

  local M, L, U, pivot, pivot_list, m, n, mn, k, rnk, r, c, apply_veil, tmp;

  # sanity check
  if has(A, V) then
    error "veiling symbol %1 is already present in matrix coefficient.", V;
    return table([]);
  end if;

  # Forget the veilings
  LEM:-VeilForget(V);

  # Get matrix dimensions
  m, n := LinearAlgebra:-Dimensions(A):

  # Create pivot vector
  r := Vector(m, k -> k);
  c := Vector(n, k -> k);

  # Check if to veil or not
  apply_veil := (z) -> LEM:-Veil[V](z);

  # Perform Gaussian elimination
  M          := copy(A);
  mn         := min(m, n);
  rnk        := mn;
  pivot_list := [];
  for k from 1 to mn do
    if LULEM:-Verbose then
      printf(
        "LULEM::LU(...): processing %d-th row, cost = %d, veilings = %d.\n",
        k, LEM:-ExpressionCost(M), nops(LEM:-VeilList(V))
      );
    end if;

    pivot := LULEM:-Pivoting(k, M, V, r, c);
    if not pivot["is_zero"] then
      pivot_list := [op(pivot_list), pivot["value"]];
    end if;

    if pivot["is_zero"] then
      rnk := rnk - 1;
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
    "method" = "LU",
    "L"      = L,
    "U"      = U,
    "V"      = V,
    "r"      = r,
    "c"      = c,
    "rank"   = rnk,
    "pivots" = pivot_list,
    "L_cost" = LEM:-ExpressionCost(L),
    "U_cost" = LEM:-ExpressionCost(U),
    "V_cost" = LEM:-ExpressionCost(LEM:-VeilList(V)),
    "L_nnz"  = nops(op(2, L)),
    "U_nnz"  = nops(op(2, U)),
    "A_nnz"  = nops(op(2, A))
  ]);
end proc: # LU

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LUsolve := proc(
  T::{table},
  b::{Vector},
  V::{symbol},
  $)::{Vector};

  description "Solve the linear system Ax=b using LU decomposition <T>, "
    "provided the vector <b> and the veiling symbol <V>.";

  local L, U, r, c, m, n, p, q, apply_veil, x, y, z, i, j, s, rnk:

  # Extract the LU decomposition
  L   := T["L"];
  U   := T["U"];
  r   := T["r"];
  c   := T["c"];
  rnk := T["rank"];

  # Get linear system dimension
  m, n := LinearAlgebra:-Dimensions(L);
  p, q := LinearAlgebra:-Dimensions(U);

  # Check if the linear system is consistent
  # sanity check
  if not ((m = n) and (p = q)) then
    error "only square system can be solved, got L = %d x %d, and U = %d x %d.",
      m, n, p, q;
    return table([]);
  end if;

  # Check if the linear system is consistent
  if not n = rnk then
    error "only full rank linear system can be solved (got rank = %1, expected "
      "rank = %2).", rnk, n;
    return table([]);
  end if;

  # Create a normalizer function
  apply_veil := (z) -> LEM:-Veil[V](z);

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
