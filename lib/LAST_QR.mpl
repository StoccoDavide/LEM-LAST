# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    ___  ____
#   / _ \|  _ \
#  | | | | |_) |
#  | |_| |  _ <
#   \__\_\_| \_\
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export QR::static := proc(
  _self::LAST,
  A::Matrix,
  {
  veil_sanity_check::boolean := true
  }, $)

  description "Compute the Givens QR decomposition of a square matrix <A>.";

  local m, n, Q, R, DG, k, j, a, b, z1, z2, r, Rk, Rj;

  # Check if the LEM object is initialized
  _self:-CheckInit(_self);

  # Sanity check
  if veil_sanity_check and has(A, V) then
    error("veiling symbol %1 is already present in matrix coefficient.", V);
  end if;

  # Clear the veiling list
  _self:-m_LEM:-ForgetVeil(_self:-m_LEM);

  # Extract the dimensions of the matrix A
  m, n := LinearAlgebra:-Dimensions(A):

  # Check if the matrix A valid
  if (m < n) then
    error("invalid matrix size A(m = %1, n = %2) detected (m >= n).", m, n);
  end if;

  # Initialize some variables
  Q  := [];      # Orthogonal transformation as a list of Given rotations
  R  := copy(A); # Transformed matrix so far
  DG := Vector[column](m, k -> 1);

  # Compute the Householder QR decomposition with veiling
  if _self:-m_VerboseMode then
    printf("LAST:-QR(...): %d x %d matrix detected.\n", m, n);
  end if;
  for k from 1 to m-1 do
    if _self:-m_VerboseMode then
      printf(
        "LAST:-QR(...): processing %d-th row, cost = %d, veilings = %d.\n",
        k, _self:-m_LEM:-ExpressionCost(_self:-m_LEM, R),
        length(_self:-m_LEM:-VeilList(_self:-m_LEM))
      );
    end if;
    for j from k+1 to n do
      a := R[k, k];
      b := R[j, k];
      if not (b = 0) then
        z1 := DG[k];
        z2 := DG[j];
        if (a = 0) then
          # Simple case: swap
          Q          := [op(Q), [k, j, 0, 0, 0, 0, 0]];
          Rk         := R[k, k..-1];
          Rj         := R[j, k..-1];
          R[k,k..-1] := Rj;
          R[j,k..-1] := Rk;
          DG[k]      := z2;
          DG[j]      := z1;
        else
          r  := _self:-m_LEM:-Veil(_self:-m_LEM, Normalizer~(z2*a^2 + z1*b^2));
          Q  := [op(Q), [k, j, a, b, z1, z2, r]];
          Rk := R[k, k+1..-1];
          Rj := R[j, k+1..-1];

          R[k, k+1..-1] := _self:-m_LEM:-Veil~(
            _self:-m_LEM, Normalizer~((a*z2)*Rk + (b*z1)*Rj)
          );
          R[j, k+1..-1] := _self:-m_LEM:-Veil~(
            _self:-m_LEM, Normalizer~(a*Rj - b*Rk)
          );

          R[k,k] := r;
          R[j,k] := 0;
          #DG[k]  := _self:-m_LEM:-Veil(_self:-m_LEM, Normalizer~(z1*z2*r));
          DG[k]  := z1*z2*r;
          DG[j]  := r;
        end if;
      end if;
    end do;
    if (R[k,k] = 0) then
      error("R[%1, %2] = 0 detected.", k, k);
    end if;
  end do;

  # Store the QR decomposition
  _self:-m_Results := table([
    "method" = "QR",
    "D"      = DG,
    "Q"      = Q,
    "R"      = R,
    "V"      = V,
    "Q_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, Q),
    "D_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, DG),
    "R_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, R),
    "V_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, LEM:-VeilList(_self:-m_LEM)),
    "Q_nnz"  = nops(op(2, Q)),
    "D_nnz"  = nops(op(2, DG)),
    "A_nnz"  = nops(op(2, A))
  ]);
  return NULL;
end proc: # QR

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export QRsolve::static := proc(
  _self::LAST,
  b::Vector,
  $)::Vector;

  description "Solve the linear system Ax=b using QR decomposition provided "
    "the vector <b>.";

  local Q, R, DG, m, n, i, j, k, c, s, a, d, x, z1, z2;

  # Check if the LEM object is initialized
  _self:-CheckInit(_self);

  # Check if the results are available
  _self:-CheckResults(_self);

  # Check if the QR decomposition is available
  if not (_self:-m_Results["method"] = "QR") then
    error("wrong or not available QR decomposition (use 'LAST:-QR()'' first).");
  end if;

  # apply Q^T a rhs
  Q  := _self:-m_Results["Q"];
  R  := _self:-m_Results["R"];
  DG := _self:-m_Results["D"];

  # Extract the dimensions of the matrix R
  m, n := LinearAlgebra:-Dimensions(R):
  x    := b;

  for i from 1 to nops(Q) do
    k  := Q[i][1];
    j  := Q[i][2];
    a  := Q[i][3];
    if (a = 0) then
      d    := x[k];
      x[k] := x[j];
      x[j] := d;
    else
      d  := Q[i][4];
      z1 := Q[i][5];
      z2 := Q[i][6];
      c  := _self:-m_LEM:-Veil~(_self:-m_LEM, Normalizer~((a*z2)*x[k] + (d*z1)*x[j]));
      s  := _self:-m_LEM:-Veil~(_self:-m_LEM, Normalizer~(a*x[j] - d*x[k]));
      x[k] := c;
      x[j] := s;
    end if;
  end do;

  # Solve R^(-1)
  x[n] := _self:-m_LEM:-Veil(_self:-m_LEM, x[n]/R[n, n]);
  for i from n-1 to 1 by -1 do
    if _self:-m_VerboseMode then
      printf("LAST:-QRsolve(...): backward substitution of %d-th row.\n", i);
    end if;
    s    := _self:-m_LEM:-Veil(_self:-m_LEM, x[i] - add(R[i, i+1..n] *~ x[i+1..n]));
    x[i] := _self:-m_LEM:-Veil(_self:-m_LEM, s/R[i, i]);
  end do;

  # Return outputs
  return x;
end proc: # QRsolve

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
