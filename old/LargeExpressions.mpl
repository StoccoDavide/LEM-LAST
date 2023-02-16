#
# LEM --- A module to help with Large Expression Management
# (c) Robert M. Corless 2001
#
# Description:
# Tools for veiling complicated coefficients in a polynomial:
#
#   collect( p, [vars], distributed, Veil[C] );
# to look at what you've done:
#   seq( C[i]=Unveil(C[i]), i=1..lastUsed );
#
#   C = name of symbol to use instead of complicated coefficients
#       (so the answer comes out like C[1]*x + C[2]*x^3 + ... etc)
#
#
# Example:
#
#   p := randpoly( [x,y,z], dense, degree=5 );
#          p := -85-61*x+85*y-99*z-55*x^4*y-37*x^4*z
#               -53*y*z+97*x^3*y^2+79*x^3*y+56*x^3*z^2
#               +49*x^3*z+23*y^2*z+88*y*z^2+57*x^2*y^3
#               +45*x^2*y^2+92*x^2*y+43*x^2*z^3-62*x^2*z^2
#               +77*x^2*z-y^3*z-86*y^2*z^2-50*y*z^3+54*x*y^4
#               +99*x*y^3-12*x*y^2-62*x*y+x*z^4-47*x*z^3
#               -91*x*z^2-47*x*z-58*y^4*z+53*y^3*z^2-85*x^5
#               -35*x^4-84*y^2+72*z^2+50*x^3*y*z-59*x^2*y^2*z
#               -8*x^2*y*z^2-93*x^2*y*z-5*x*y^3*z+63*x^3+94*y^3
#               +17*z^3+66*x^2-90*y^4+78*z^4-61*x*y^2*z^2-50*x*y^2*z
#               -18*x*y*z^3+31*x*y*z^2-26*x*y*z+41*y^5+49*z^5+83*y^2*z^3+19*y*z^4
#
#   compact := collect( p, [x,y], distributed, Veil[K] );
#          compact := K[1]+K[2]*x+K[3]*y-55*x^4*y+97*x^3*y^2+K[4]*x^3*y
#                     +57*x^2*y^3-K[5]*x^2*y^2-K[6]*x^2*y+54*x*y^4-K[7]*x*y^3
#                     -K[8]*x*y^2-K[9]*x*y-85*x^5-K[10]*x^4+K[11]*y^2
#                     +7*K[12]*x^3+K[13]*y^3+K[14]*x^2-2*K[15]*y^4+41*y^5
#
# which is obviously much shorter.  The data are still available, e.g.
# K[1] = Unveil( K[1] ) = -85+72*z^2+78*z^4-99*z+17*z^3+49*z^5.
#
# Verify:
#   zero := Unveil( compact, infinity ) - p;
#              zero := big ugly expression
#   normal( zero ); # returns zero.
#
#
# LEM --- A procedure to help with Large Expression Management
# (c) Robert M. Corless 2001
# Acknowledgement: David Jeffrey, Michael Monagan, Jacques Carette


#
# A procedure to help with Large Expression Management
# (c) Robert M. Corless 2001
# Acknowledgement: David Jeffrey, Michael Monagan, Jacques Carette

unprotect('LargeExpressions'):
module LargeExpressions()
      export Veil, Unveil, LastUsed;
      local auxiliary, labelledValues;
      global _V; # default name to use for veiling.
      option package;

      LastUsed := table('sparse'); # Begin with nothing recorded

      # We use a table to store the labelled values.
      labelledValues := table():

      Veil := proc( coefficient )
         local i, s, c, label;
         description "hide expressions behind labels.";

	 label := `if`(procname::indexed, op(procname), '_V');

	  if label<>eval(label,2) then
	      error "label %a is assigned a value already, please save its contents and unassign it", label;
	  end if;

         # Recognize zero if we can, so that we don't hide zeros.
         c := Normalizer( coefficient );
         # Remove the integer content and sign so that we don't hide them either.
         i := icontent( c );
         # And we really mean sign, here, and not signum,
         # because the interesting case is when c may be a polynomial.
         try
            s := sign( c ); # sign is weak
            # Don't do anything if we can tell that
            # the coefficient is just a number or a simple
            # multiple of a name (simple or indexed)
            if s*i=c or type(s*c/i,'name') then return c end if;
         catch:
            s := 1;
         end try;
         # Only if there is something complicated to hide
         # do we actually hide anything.
         s*i*auxiliary( s*c/i, label );
      end proc;

      Unveil := proc( c, ilevel::{nonnegint,infinity} )
         local a, b, i, level, label;
         description "reveal expressions hidden behind labels.";
	 label := `if`(procname::indexed, op(procname), '_V');
         level := `if`( nargs<2, 1, min(LastUsed[label]+1,ilevel) );
         a := c;
         # Always do at least 1
         b := eval( a, [seq(label[i]=labelledValues[label,i],i=1..LastUsed[label])] );
         from 2 to level while not Testzero(a-b) do
            a := b;
            b := eval( a, [seq(label[i]=labelledValues[label,i],i=1..LastUsed[label])] );
         end do;
         b;
      end proc;

      # Scope LastUsed etc and use option remember to detect duplicates.
      # There is no nontrivial storage duplication because objects are hashed and
      # stored uniquely.  All this costs is a pointer.
      auxiliary := proc( c, label )
         option remember;
         unprotect(LastUsed);
         LastUsed[label] := LastUsed[label] + 1;
         protect(LastUsed);
         labelledValues[ label, LastUsed[label] ] := c;
         label[ LastUsed[label] ]
      end proc:

end module:

#savelib(LargeExpressions):
