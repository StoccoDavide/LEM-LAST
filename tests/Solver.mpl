
read "LinearSystem.mpl";
with(LinearAlgebra):
Vars := [seq(x||i,i=1..44)];
Eqns := GenerateEquations(A,Vars,B);
Eqns := {op(Eqns)}:
Vars := {op(Vars)}:
infolevel[solve] := 2;
Sols := solve(Eqns,Vars):
Sols := factor(Sols):


Structure := proc(EE,Y)
local E,T,N,e,f,S,U;
T := table();
U := table();
N := 1;
S := {};
E := factor(EE);
for e in E do
   if rhs(e)=0 then next fi;
   for f in factors(rhs(e))[2] do
       if assigned( T[f[1]] ) then next fi;
       if assigned( T[-f[1]] ) then print(HELLO) fi;
       T[f[1]] := true; # seen it
       U[Y[N]] := f[1]; # store it
       S := S union {f[1] = Y[N]};
       N := N+1;
   od;
od:

subs( S, E ), U;

end:


SSols,Z := Structure(Sols,z):
Zvals := op(2,eval(Z)):  # convert table to list of equations

numterms := proc(f) if f=0 then 0 elif type(f,`+`) then nops(f) else 1 fi end:
for var in map(op,[indices(Z)]) do printf("#%a = %d\n",var, numterms(Z[var])) od;

for xsol in SSols do printf("%a\n",xsol) od;
for zval in Zvals do printf("\n%a\n",zval); od;

