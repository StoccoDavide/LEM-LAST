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
    max_iter::posint := 10,
    $)::algebraic;

    description "Compute the signature of an expression <expr> modulo a prime "
      "number <p>.";

    return SIG:-Signature(SIG:-Transform(expr, max_iter), p);
  end proc: # ModuleApply

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ExtractArgs := proc(
    expr::algebraic,
    func::name,
    $)::list(algebraic);

    description "Extract the arguments of a function <func> from an expression "
      "<expr>.";

    return [op(map(op, indets(expr, ':-specfunc(func)')))];
  end proc: # ExtractArgs

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export AbsTransform := proc(
    expr::algebraic,
    $)::algebraic;

    description "Transform the absolute value function in an expression <expr>.";

    local new_abs, arg_abs;

    new_abs := u -> sqrt(u);
    arg_abs := SIG:-ExtractArgs(expr, 'abs');

    return subs[eval](
      op(abs~(arg_abs) =~ new_abs~(arg_abs)),
      expr
    );
  end proc: # AbsTransform

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export LogTransform := proc(
    expr::algebraic,
    $)::algebraic;

    description "Transform the logarithm function in an expression <expr>.";

    local new_log, arg_log;

    new_log := u -> u;
    arg_log := SIG:-ExtractArgs(expr, 'ln');

    return subs[eval](
      op(abs~(arg_log) =~ new_abs~(arg_log)),
      expr
    );
  end proc: # LogTransform

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export WeierstrassTransform := proc(
    expr::algebraic,
    $)::algebraic;

    description "Transform the trigonometric functions in an expression <expr>.";

    local new_cos, new_sin, new_tan, arg_cos, arg_sin, arg_tan;

    new_cos := u -> (1-u^2)/(1+u^2);
    new_sin := u -> (2*u)/(1+u^2);
    new_tan := u -> new_sin(u)/new_cos(u);
    arg_cos := SIG:-ExtractArgs(expr, 'cos');
    arg_sin := SIG:-ExtractArgs(expr, 'sin');
    arg_tan := SIG:-ExtractArgs(expr, 'tan');

    return subs[eval](
      op(cos~(arg_cos) =~ new_cos~(arg_cos)),
      op(sin~(arg_sin) =~ new_sin~(arg_sin)),
      op(tan~(arg_tan) =~ new_tan~(arg_tan)),
      expr
    );
  end proc: # WeierstrassTransform

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Transform := proc(
    expr::algebraic,
    max_iter::posint := 10,
    $)::algebraic;

    description "Transform an expression <expr> using a maximum number of "
      "substitutions recursions <max_iter>.";

    local out, i;

    # Copy the input expression
    out := copy(expr);

    # Indexed type elements
    if hastype(out, indexed) then
      error("the input expression has indexed type elements.");
    end if;

    # Absolute value function elements
    for i from 1 to max_iter while has(out, {'abs'}) do
      out := SIG:-AbsTransform(out);
    end do;
    if has(out, {'abs'}) then
      error("the input expression has absolute value function elements.");
    end if;

    # Logarithm function elements
    for i from 1 to max_iter while has(out, {'ln'}) do
      out := SIG:-LogTransform(out);
    end do;
    if has(out, {'ln'}) then
      error("the input expression has logarithm function elements.");
    end if;

    # Trigonometric function elements
    for i from 1 to max_iter while has(out, {'sin', 'cos', 'tan'}) do
      out := SIG:-WeierstrassTransform(out);
    end do;
    if has(out, {'sin', 'cos', 'tan'}) then
      error("the input expression has trigonometric function elements.");
    end if;

    return out;
  end proc: # Transform

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Signature := proc(
    expr::algebraic,
    p::prime,
    $)::integer;

    description "Compute the signature of an expression <expr> modulo a prime "
      "number <p>.";

    local out;

    try
      out := signature(expr, p);
    catch:
      WARNING("expression signature not defined, assumed '0'.");
      out := 0;
    end try;

    return out;
  end proc: # Signature

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # SIG

# That's all folks!
