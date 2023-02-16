# Signature utilities
module Sig()
    export ModuleApply, UnVeil, Veil, LastUsed, forgetVeil;
    local SIG, UnVeilTable, NextLabel;

	LastUsed := table('sparse');
    UnVeilTable := table():

    NextLabel := proc(label::symbol,vars::set)
	    LastUsed[label] := LastUsed[label] + 1;
		label[LastUsed[label],`if`(nargs=2, vars, NULL)];
	end proc;

    Veil := proc(x,p) local A, label;
        label := op(procname);
        # s := SIG(x, p, label);
        # inlining NextLabel for efficiency
	    LastUsed[label] := LastUsed[label] + 1;
		A := label[LastUsed[label]];
        UnVeil(A,1) := x;
        A;
    end:
#    VeilDep := proc(x,p) local A, s;
#        s := SIG(x, p, op(procname));
#	    A := NextLabel(op(procname),op(3,s));
#        UnVeilTable[A] := x;
#        A;
#    end:
    SIG := proc(f,p::posint,A)
        local t, sig;
        option remember;
        if f::'rational' then f mod p;
        elif f::'symbol' then SIG(f,p,A) := rand() mod p;
        elif f::'indexed' then
            if type(f, A[anything]) then
                SIG(UnVeil(f,1),p,A)
            else
                SIG(f,p,A) := rand() mod p;
            end if;
        elif f::`+` then
             # is ok since this does not grow too much
             add(SIG(t,p,A) mod p, t=f) mod p;
             # sig := 0;
             # for t in f do
             #     T := SIG(t,p,A); sig := sig+T[2] mod p; od;
             # sig;
        elif f::`*` then
             # mul(SIG(t,p,A) mod p, t=f) mod p;
             sig := 1;
             for t in f do
                 sig := sig*SIG(t,p,A) mod p;
             od;
             sig;
        elif f::(anything^rational) then
             SIG(op(1,f),p,A) &^ op(2,f) mod p;
        # elif type(f,Matrix) then
        #      map(SIG,f,p,A);
        else error "expressions involving %1 not done", f;
        fi;
    end:
    ModuleApply := SIG;
    UnVeil := proc(x,n::{nonnegint,identical(infinity)})
        option remember;
        if nargs=1 then UnVeil(x,1)
        elif (nargs=2 and n=0) then x
#        elif assigned(UnVeilTable[x]) then
#            if n=1 then # don't bother recursing!
#                UnVeilTable[x];
#            else
#                UnVeil(UnVeilTable[x],n-1)
#            end if;
        elif x::atomic then x
        else map(UnVeil, x, n-1)
        end if;
    end proc;
    # Signature := proc(V) option inline; V[2] end:
    # Val := proc(V) option inline; V[1] end:
    forgetVeil := proc(label::symbol)
        local i;
        subsop(4=NULL, eval(SIG));
        LastUsed[label] := 0;
        UnVeilTable := table();
    end proc:
end module:

# This procedure is to computer n by m matrix LU decomposation. (Now we just think about n= m square system)
#This LU with LEM  and pivoting with different strategies.
# also about singular case and zero recognition strategies.
#In general, any square matrix A, singular or nonsingular, has a factorization PA=LU, where P is a permutation matrix,
#L is  lower matrix with 1's as diagonal entries,  U is upper matrix.

ZeroStrategy_Signature := proc(f,p,A) Sig(f,p,A) end:

SquareLUwithpivoting := proc(xA::Matrix, Strategy_LEM::procedure, Q::symbol, Prim::posint, Strategy_Pivots::procedure, Strategy_Zero::procedure)

# Strategy_LEM is for all the strategy for LEM, and Q is the sequence of vars used by veil.
#Strategy_Pivots is for all the strategy for choosing pivots.

local A, n, k, i, ii, j, m, L, U, mltip, l, r, kp, p, temp, normalize,
      row, flag, det, z, t;
A := copy(xA);
n := LinearAlgebra[RowDimension](A):
m := LinearAlgebra[ColumnDimension](A);
if m<>n then error "The input matrix should be square!" end if;
L := Matrix(n,n,shape=triangular[lower]);  # L is lower matrix
U := Matrix(n,n,shape=triangular[upper]):  # U is upper matrix.

r := Vector(n, k->k );                    # create pivot vector
normalize := z ->`if`(Strategy_LEM(z) > 0, Sig:-Veil[Q](z,Prim), z);

row := 1;
det := 1;
for kp from 1 to m while row <= n do             # loop over columns
    # We find the pivot element among the entries
    # A[r[kp],kp], A[r[kp+1],kp], .... A[r[n],kp]
    # according to different pivoting strategies.
    # Interchange entries in pivot vector/
    p := row;
    flag := evalb(Strategy_Zero(A[r[row],kp], Prim, Q) = 0);
    for i from row+1 to n do
        if (flag or Strategy_Pivots( A[r[i],kp], A[r[p],kp] )) and not
           (Strategy_Zero(A[r[i],kp], Prim, Q)=0) then
            p := i;
            break;  # once we've found a pivot -- not ``best'' !
        end if;
    end do;

    # only when the pivot is not equal to zero, we will do the row
    # elimination (the whole if statement). else we will continue with
    # next colomn.
    if Strategy_Zero(A[r[p],kp], Prim, Q)=0 then
        WARNING("the matrix appears to be singular.");
    else
        if (p <> row) then
            (r[p], r[row]) := (r[row], r[p]);
        end if;

        userinfo( 3, SquareLUwithpivoting, ` kp `, kp, ` r `, r );

        # Now do Guass elimination steps to get new A
        # Since it is fraction-free, we must keep track of ``diagonal''
        # of L explicitly?
        for i from row+1 to n do
            for j from kp+1 to m do
                z := A[r[i],j]*A[r[row],kp] - A[r[row],j]*A[r[i],kp];
                t := normal(z/det);
                A[r[i], j] := `if`(Strategy_LEM(t) > 0,
                    Sig:-Veil[Q](t,Prim), t);
            end do;
        end do;
        det := A[r[row],kp];
        userinfo( 3, SquareLUwithpivoting, ` A `, A );
        row := row + 1;
    end if;
end do:
userinfo( 2, SquareLUwithpivoting, ` r `, r , ` A `, A);

# Seperate new A into L and U
for j from 1 to m do
    for i from 1 to n do
        if i < j then
            U[i, j] := A[r[i], j];
        elif i=j then
            U[i, i] := A[r[i], i];
            L[i, i] := `if`(i=1,1, 1/A[r[i-1],i-1]);
        else # i>j
            t := mul(A[r[z],z],z=1..j);
            t := `if`(t=0,1,t);
            L[i, j] := A[r[i], j]/t;
        end if;
  end do;
end do;

L,U, r;
end proc:

# Solve the linear system Ax=b given the LU decomposition ( PA=LU ) in a returned from SquareLUwithpivoting
# The pivot vector r returned from SquareLUwithpivoting and the right hand side vector b.

SolveSquareLUpivot := proc( A, r, b, Strategy_LEM, P, Q)
local y, x, i, s, j, n, normalizer;

   userinfo( 2, SolveSquareLUpivot, ` b `, b, ` r `, r );

   n := LinearAlgebra[RowDimension](A);       # get number of rows in matrix A
   y := Vector(n);                            # create vector for solution of Ly=Pb
   x := Vector(n);                            # create vector for solution of Ux=y

   normalizer := y -> `if`( Strategy_LEM(y) > 0, Sig:-Veil[Q](y,P), y);
   # do forward substitution to solve Ly=Pb

   userinfo( 3, SolveSquareLUpivot, ` n `, n, ` y `, y );

   y[1] := normalizer( b[r[1]] );
   for i from 2 to n do
       y[i] := normalizer( b[r[i]] ) - add(normalizer( A[r[i],j]*y[j] ), j=1..i-1) ;
   od;

   # do backward substitution to solve Ux=y
   x[n] := normalizer( y[n]/A[r[n],n] );
   for i from n-1 to 1 by -1 do
      s := normalizer( y[i] ) - add(normalizer( A[r[i],j]*x[j] ), j=i+1..n);
      x[i] := normalizer( s/ A[r[i],i] );
   od;

   x;                                          # return solution vector
end proc:


# nvar is the max number of indets in one expression which will not be replaced by K[i]
#nops(indets(x))>3 then
#This means we just want to subs the expression more than 2 variables.
#So it will decrease the subs number but the computation time will increase!!!
LEMStrategy_n := proc(x) nops(indets(x)) - 4; end proc;

#n is the max length which will not be replaced by K[i].
#When we choose the parameter value: n=12, the number of substitution is the same as ChadLU code.
#When we choose the n=8, the number of substitution is larger than ChadLU code, but not too much!
#When we choose the parameter value: n>12, the number of substitution is still the same as ChadLU code.
LEMStrategy_L:= proc(x) length(x) - 50; end proc;

#n is the max length which will not be replaced by K[i].
#When we choose the parameter value: n=12, the number of substitution is the same as ChadLU code.
#When we choose the n=8, the number of substitution is larger than ChadLU code, but not too much!
#When we choose the parameter value: n>12, the number of substitution is still the same as ChadLU code.
LEMStrategy_LB:= proc(x) length(x) - 260 end proc;

#n is the max length which will not be replaced by K[i].
#When we choose the parameter value: n=12, the number of substitution is the same as ChadLU code.
#When we choose the n=8, the number of substitution is larger than ChadLU code, but not too much!
#When we choose the parameter value: n>12, the number of substitution is still the same as ChadLU code.
LEMStrategy_Ls:= proc(x) length(x) - 120 end proc;

# x , y are the entries from a symbolic matrix
# the output is the length difference of them.
# this strategy is to choose the entry with the longest length as pivot.
PivotStrategy_Llength := proc( x , y ) evalb( length(x) - length(y) > 0 ); end proc;

# x , y are the entries from a symbolic matrix
# the output is the difference of their length.
# this strategy is to choose the entry with the shortest length as pivot.
PivotStrategy_Slength := proc( x , y ) evalb( length(x) - length(y) < 0 ); end proc;

# x, y are the entries from a symbolic matrix
# the output is the difference of the independent vars
# this strategy is to choose the entry with the biggest indets number as pivot
PivotStrategy_Lindets := proc ( x, y) evalb ( (nops(indets(x)) - nops(indets(y))) > 0) end proc;

# x, y are the entries from a symbolic matrix
# the output is the difference of the independent vars
# this strategy is to choose the entry with the smallest indets number as pivot
PivotStrategy_Sindets := proc ( x, y) evalb ( (nops(indets(x)) - nops(indets(y))) < 0) end proc;

# x, y are the entries from a symbolic matrix
# the output is the difference of the independent vars
# this strategy is to choose the entry with the larger number as pivot
PivotStrategy_numeric := proc ( x, y) evalb ( (abs(x) - abs(y) ) > 0); end proc;

# x is the element needed to be recognized if it is zero or not
# the output is the length of x
ZeroStrategy_length := proc( x, zero, K_LEM ) length(x) end proc;

# x is the element needed to be test if it is zero or not
# the output is the normalized form of  x
ZeroStrategy_normalizer := proc( x, zero, K_LEM) Normalizer(x) end proc;
