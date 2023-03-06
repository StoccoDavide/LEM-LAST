  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  QR := proc(
    A::{Matrix},
    V::{symbol},
    VeilingStrategy::{procedure}  := VeilingStrategy_n,
    $)::{table};

    description "Compute the Householder QR decomposition of a square matrix <A> "
                "using the veiling strategy <VeilingStrategy> and the veiling symbol <V>.";

    local m, n, z, k, Q, R, M, norm_x, s, u1, w, tau;

    LEM:-VeilForget(V);

    # Extract the dimensions of the matrix A
    m, n := LinearAlgebra[Dimensions](A):

    # Check if the matrix A valid
    assert(
      m >= n,
      "LULEM::QR(...): invalid matrix A(m,n) detected, got (m >= n)."
    );

    # Initialize some variables
    z := Vector(n, k -> 0);
    Q := Matrix(LinearAlgebra[IdentityMatrix](m)); # Orthogonal transformation so far
    R := copy(A);                                  # Transformed matrix so far

    # Compute the Householder QR decomposition with veiling
    for k from 1 to n do

      # Find H = I-tau*w*w’ to put zeros below R[j,j]
      norm_x := LinearAlgebra[Norm](R[k..-1,k], 2);
      s      := -sign(R[k, k]);
      u1     := R[k, k] - s*norm_x;
      w      := R[k..-1,k]/u1;
      w[1]   := 1;
      tau    := -s*u1/norm_x;

      # Update R = HR and Q = QH
      R[k..-1, 1..-1] := R[k..-1, 1..-1] - tau.w.LinearAlgebra[Transpose](w).R[k..-1,1..-1];
      Q[1..-1, k..-1] := Q[1..-1, k..-1] - tau.Q[1..-1, k..-1].w.LinearAlgebra[Transpose](w);

    end do;

    # Return the QR decomposition
    return table([
      "method"   = "QR",
      "Q"        = Q,
      "R"        = R,
      "V"        = V,
      "Q_length" = length(convert(Q,list)),
      "R_length" = length(convert(R,list)),
      "V_length" = length(LEM:-VeilList(V))
    ]);
  end proc: # QR

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  (*
  SolveQR := proc(
    T::{table},
    b::{Vector},
    V::{symbol, function},
    VeilingStrategy::{procedure} := VeilingStrategy_n,
    $)

    description "Solve the linear system Ax=b using QR decomposition <T>, "
      "provided the vector <b>, the veiling symbol <V> and the veiling strategy "
      "<VeilingStrategy>.";

  end proc: # SolveQR
  *)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  (*
  FFQR := proc(
    A::{Matrix},
    V::{symbol},
    VeilingStrategy::{procedure}  := VeilingStrategy_n,
    $)::{table};

    description "Compute the Fracton-Free QR decomposition of a square matrix "
      "<A> using the veiling strategy <VeilingStrategy> and the veiling symbol <V>.";

  end proc: # QR
  *)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  (*
  SolveFFQR := proc(
    T::{table},
    b::{Vector},
    V::{symbol, function},
    VeilingStrategy::{procedure} := VeilingStrategy_n,
    $)

    description "Solve the linear system Ax=b using FFQR decomposition <T>, "
      "provided the vector <b>, the veiling symbol <V> and the veiling strategy "
      "<VeilingStrategy>.";

  end proc: # SolveFFQR
  *)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
