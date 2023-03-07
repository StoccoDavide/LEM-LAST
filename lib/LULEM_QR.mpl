# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

QR := proc(
  A::{Matrix},
  V::{symbol},
  VeilingStrategy::{procedure}  := VeilingStrategy_n,
  $)::{table};

  description "Compute the Givens QR decomposition of a square matrix <A> "
              "using the veiling strategy <VeilingStrategy> and the veiling symbol <V>.";

  local m, n, z, a, b, c, s, r, i, j, k, Q, R, apply_veil;

  LEM:-VeilForget(V);

  # Check if to veil or not
  apply_veil := (z) -> `if`(VeilingStrategy(z), LEM:-Veil[V](z), z);

  # Extract the dimensions of the matrix A
  m, n := LinearAlgebra[Dimensions](A):

  # Check if the matrix A valid
  assert(
    m >= n,
    "LULEM::QR(...): invalid matrix A(m,n) detected, got (m >= n)."
  );

  # Initialize some variables
  Q := [];      # Orthogonal transformation as a list of Given rotations
  R := copy(A); # Transformed matrix so far

  # Compute the Householder QR decomposition with veiling
  for k from 1 to m-1 do
    if LULEM:-Verbose then
      printf(
        "LULEM::QR(...): processing %d-th colum. Length %d\n",
        k, length(LEM:-VeilList(V))
      );
    end;
    for j from k+1 to n do
      a := R[k,k];
      b := R[j,k];
      if not b = 0 then
        r := apply_veil(sqrt(a^2+b^2));
        c := apply_veil(a/r);
        s := apply_veil(-b/r);

        Q := [op(Q),[c,s]];

        R[k,k] := r;
        R[j,k] := 0;
        for i from k+1 to m do
          a      := c*R[k,i]-s*R[j,i];
          b      := s*R[k,i]+c*R[j,i];
          R[k,i] := apply_veil(Normalizer(a));
          R[j,i] := apply_veil(Normalizer(b));
        end do;
      end if;
    end do;
    if R[k,k] = 0 then
      error "R[%a,%a] = 0", k, k;
    end if;
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
