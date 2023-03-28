# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    ___  ____  ____
#   / _ \|  _ \|___ \
#  | | | | |_) | __) |
#  | |_| |  _ < / __/
#   \__\_\_| \_\_____|
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export QR2::static := proc(
  A::Matrix,
  $)::table;

  description "Compute the Givens QR decomposition of a square matrix <A>.";

  local V, m, n, Q, R, k, j, a, b, c, r, Rk, Rj, l, C, apply_veil;

  # Get the veiling label
  V := _self:-m_LEM:-GetVeilingLabel(_self:-m_LEM);

  # Sanity check
  if has(A, V) then
    error "veiling symbol %1 is already present in matrix coefficient.", V;
  end if;

  # Clear the veiling list
  _self:-m_LEM:-VeilForget(_self:-m_LEM);

  # Check if to veil or not
  apply_veil := (z) -> _self:-m_LEM:-Veil(_self:-m_LEM, z);

  # Extract the dimensions of the matrix A
  m, n := LinearAlgebra:-Dimensions(A):

  # Check if the matrix A valid
  if m < n then
    error "invalid matrix size A(m = %1, n = %2) detected (m >= n).", m, n;
  end if;

  # Initialize some variables
  Q := [];      # Orthogonal transformation as a list of Given rotations
  R := copy(A); # Transformed matrix so far

  # Create pivot vector
  c := Vector(n, k -> k);

  # Compute the Householder QR decomposition with veiling
  for k from 1 to m-1 do
    for l from k to n do
      if _self:-m_Verbose then
        printf(
        "LULEM::QR2(...): processing %d-th row, cost = %d, veilings = %d.\n",
        k, _self:-m_LEM:-ExpressionCost(_self:-m_LEM, R),
        length(_self:-m_LEM:-VeilList(_self:-m_LEM))
        );
      end if;
      if (l > k) then
        if _self:-m_Verbose then
          printf("LULEM::QR2(...): swap with colum %d\n", l+1);
        end if;
        C          := R[1..-1,k]; a    := c[k];
        R[1..-1,k] := R[1..-1,l]; c[k] := c[l];
        R[1..-1,l] := C;          c[l] := a;
      end if;
      for j from k+1 to n do
        # Check if expressions are equal to 0 (very costly)
        a := Normalizer(_self:-m_LEM:-UnVeil(_self:-m_LEM, R[k, k]));
        b := Normalizer(_self:-m_LEM:-UnVeil(_self:-m_LEM, R[j, k]));
        if (b = 0) then
          R[j,k] := 0;
          if (a = 0) then
            R[k,k] := 0;
          end if;
        else
          if (a = 0) then
            R[k,k] := 0;
            # simple case do a swap
            Q           := [op(Q), [k, j, 0, 0]];
            Rk          := R[k, k..-1];
            Rj          := R[j, k..-1];
            R[k, k..-1] := Rj;
            R[j, k..-1] := Rk;
          else
            a := R[k, k];
            b := R[j, k];
            r := apply_veil(Normalizer(a^2 + b^2));
            Q := [op(Q), [k, j, a, b]];

            Rk := R[k, k+1..-1];
            Rj := R[j, k+1..-1];

            R[k, k] := r; R[k, k+1..-1] := apply_veil~(Normalizer~(a*Rk + b*Rj));
            R[j, k] := 0; R[j, k+1..-1] := apply_veil~(Normalizer~(a*Rj - b*Rk));
          end if;
        end if;
      end do;
      if not (R[k, k] = 0) then
        break;
      end if;
    end do;
    if R[k,k] = 0 then
      printf("R[%a,%a] = 0", k, k);
      break;
    end if;
  end do;

  # Return the QR decomposition
  return table([
    "method" = "QR2",
    "Q"      = Q,
    "R"      = R,
    "V"      = V,
    "c"      = c,
    "Q_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, Q),
    "R_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, R),
    "V_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, _self:-m_LEM:-VeilList(_self:-m_LEM)),
    "Q_nnz"  = nops(op(2, Q)),
    "A_nnz"  = nops(op(2, A))
  ]);
end proc: # QR2

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export QR2solve::static := proc(
  T::table,
  xb::Vector,
  $)::Vector;

  description "Solve the linear system Ax=b using QR decomposition <T> and "
    "provided the vector <b>.";

  local Q, R, m, n, i, j, k, z, c, s, a, b, x, z1, z2, r, apply_veil;

  # Check if to veil or not
  apply_veil := (z) -> _self:-m_LEM:-Veil(_self:-m_LEM, z);

  # apply Q^T a rhs
  Q  := T["Q"];
  R  := T["R"];

  # Extract the dimensions of the matrix R
  m, n := LinearAlgebra:-Dimensions(R):
  x    := xb;

  for i from 1 to nops(Q) do
    k  := Q[i][1];
    j  := Q[i][2];
    c  := Q[i][3];
    s  := Q[i][4];
    if c = 0 then
      b    := x[k];
      x[k] := x[j];
      x[j] := b;
    else
      a    := apply_veil~(Normalizer~(c*x[k] + s*x[j]));
      b    := apply_veil~(Normalizer~(c*x[j] - s*x[k]));
      x[k] := a;
      x[j] := b;
    end if;
  end do;
  # solve R^(-1)

  x[n] := apply_veil(x[n] /R[n, n]);
  for i from n-1 to 1 by -1 do
    if _self:-m_Verbose then
      printf("LULEM::QR2solve, backward %d\n",i);
    end if;
    s    := apply_veil(x[i] - add(R[i, i+1..n] *~ x[i+1..n]));
    x[i] := apply_veil(s / R[i, i]);
  end do;

  # Return outputs
  return x;
end proc: # QR2solve

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
