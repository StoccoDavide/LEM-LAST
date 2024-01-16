# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                ____ ___ ____                                #
#                               / ___|_ _/ ___|                               #
#                               \___ \| | |  _                                #
#                                ___) | | |_| |                               #
#                               |____/___\____|                               #
#                         Expression Signature Module                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco     (University of Trento)
#   Enrico Bertolazzi (University of Trento)
#
# License: BSD 3-Clause License
#
# This is a module for the 'SIG' (Signature) module. It contains the functions
# and workarounds to calculate the signature of "almost" any kind of expression.
# The signature is a hashing value that can be used to calculate the similarity
# between two expressions, or to check how likely an expression is algebraically
# non-null (see Schwarz-Zippel Lemma).

SIG := module()

  description "Expression Signature module.";

  option package,
         load   = ModuleLoad,
         unload = ModuleUnload;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Info := proc()

    description "Print 'SIG' module information.";

    printf(
      "+--------------------------------------------------------------------------+\n"
      "| 'SIG' module version 1.0 - BSD 3-Clause License - Copyright (c) 2023     |\n"
      "| Current version authors:                                                 |\n"
      "|   D. Stocco and E. Bertolazzi.                                           |\n"
      "+--------------------------------------------------------------------------+\n"
    );
    return NULL;
  end proc: # Info

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleLoad := proc()
    description "'SIG' module load procedure.";
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleUnload := proc()
    description "'SIG' module unload procedure.";
    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleApply := proc(
    expr::algebraic,
    p::prime,
    {
    verbose::boolean := false
    }, $)::anything;

    description "Compute the signature of an expression <expr> modulo a prime "
      "number <p>.";

    return SIG:-Signature(expr, p, parse("verbose") = verbose);
  end proc: # ModuleApply

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Signature := proc(
    expr::algebraic,
    p::prime,
    {
    verbose::boolean := false
    }, $)::anything;

    description "Compute the signature of an expression <expr> modulo a prime "
      "number <p>.";

    local out;

    try
      out := signature(expr, p);
    catch:
      if verbose then
        WARNING(
          "LEM:-Signature(...): expression signature not defined, assumed '1'."
        );
      end if;
      out := 1; # FIXME: is 1 the best choice?
    end try;

    return out;
  end proc: # Signature

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # SIG

# That's all folks!
