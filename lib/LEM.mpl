# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                              _     _____ __  __                             #
#                             | |   | ____|  \/  |                            #
#                             | |   |  _| | |\/| |                            #
#                             | |___| |___| |  | |                            #
#                             |_____|_____|_|  |_|                            #
#                         Large Expressions Management                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Authors of the current version:
#  Davide Stocco     (University of Trento)
#  Matteo Larcher    (University of Trento)
#  Enrico Bertolazzi (University of Trento)
#
# Authors of the original code:
#   Wenqin Zhou       (University of Western Ontario) - Former affiliation
#   David J. Jeffrey  (University of Western Ontario)
#   Jacques Carette   (McMaster University)
#   Robert M. Corless (University of Western Ontario)
#
# License: BSD 3-Clause License
#
# This is a module for the 'LEM' (Large Expressions Management) package.
# It contains the functions to deal with large expressions management.
# The module can be used to veil and unveil large expressions to avoid
# expression swell.
#
# The module is very similar to the 'LargeExpression' module provided by Maple.
# The main difference is that the 'LEM' module can handle more than one veiling
# label at a time. Moreover, it has some built-in functions to display, list and
# substitute the veiled expressions.
#
# The following code is hopefully an improved version of the original version
# provided inthe following PhD thesis:
#
#   Wenqin Zhou, Symbolic Computation Techniques for Solveing Large Expressions
#   Problems from Mathematics and Engineering (2007), Faculty of Graduate Studies,
#   The University of Western Ontario London, Ontario, Canada.
#
# We would like to thank Jacques Carette for providing the original code that
# have been used to develop this module.

LEM := module()

  export  Veil,
          UnVeil,
          UnVeilImap,
          VeilUnorderedList,
          VeilList,
          VeilTableSize,
          VeilTableImap,
          VeilTableAppend,
          VeilLabels,
          VeilSubs,
          VeilForget,
          ExpressionCost,
          VeilingStrategy,
          SetVeilingStrategyPars;

  local   ModuleLoad,
          ModuleUnload,
          UnVeilTables,
          UnVeilLabels,
          Auxiliary,
          VeilingStrategy_maxcost,
          VeilingStrategy_subscripts,
          VeilingStrategy_assignments,
          VeilingStrategy_additions,
          VeilingStrategy_multiplications,
          VeilingStrategy_divisions,
          VeilingStrategy_functions;

  option  package,
          load   = ModuleLoad,
          unload = ModuleUnload;

  description "Large Expressions Management module.";

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ModuleLoad := proc()

    description "'LEM' module load procedure.";

    local i, lib_base_path;

    printf(
      "'LEM' module version 1.0 - BSD 3-Clause License - Copyright (c) 2023\n"
      "Current version: D. Stocco, M. Larcher, E. Bertolazzi.\n"
      "Original code: W. Zhou, D. J. Jeffrey, J. Carette and R. M. Corless.\n"
    );

    lib_base_path := null;
    for i in [libname] do
      if (StringTools[Search]("LEM", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if (lib_base_path = null) then
      error "Cannot find 'LEM' module";
    end if;

    # Initialize internal variables
    LEM:-UnVeilTables := table([]);
    LEM:-SetVeilingStrategyPars();

    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ModuleUnload := proc()

    description "'LEM' module unload procedure.";

    LEM:-UnVeilTables := NULL;

    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Veil := proc(
    x::{anything},
    {force::{boolean} := false},
    $)::{anything};

    description "Check if the veiling strategy is verified and veil an expression <x> and return a label to it.";

    local label, i, s, c;

    # Retreive the label variable
	  label := `if`(procname::{indexed}, op(procname), '_V');

    # Check if label is already assigned
	  if (label <> eval(label, 2)) then
	    error "LEM::Veil(...): label %a is already assigned, please save its "
        "contents and unassign it.", label;
	  end if;

    # Recognize zero if we can, so that we don't hide zeros.
    c := Normalizer(x);

    # Check the veiling strategy
    if not(force) and not(LEM:-VeilingStrategy(c)) then
      return c;
    end if;

    # Remove the integer content and sign so that we don't hide them either.
    i := icontent(c);

    # And we really mean sign, here, and not signum, because the interesting
    # case is when c may be a polynomial.
    try
      s := sign(c); # sign is weak
      # Don't do anything if we can tell that the coefficient is just a number
      # or a simple multiple of a name (simple or indexed)
      if (s*i = c) or type(s*c/i, indexed(integer)) or type(s*c/i, name) then
        return c;
      end if;
    catch:
      s := 1;
    end try;
    # Only if there is something complicated to hide we do actually hide it and
    # return a label.
    return s * i * LEM:-VeilTableAppend(label, s*c/i);
  end proc; # Veil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ExpressionCost := proc(
    x::{anything},
    $)::{integer};

    description "Compute the cost of the expression <x>.";

    local tmp;

    if not(type(x, algebraic) or type(x, list)) then
      tmp := convert(x, list);
    else
      tmp := x;
    end if;

    return subs(
      'subscripts'      = 0,
      'assignments'     = 0,
      'additions'       = 1,
      'multiplications' = 2,
      'divisions'       = 3,
      'functions'       = 2,
      codegen:-cost(tmp)
    );
  end proc: # ExpressionCost

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilingStrategy := proc(
    x::{algebraic},
    $)::{boolean};

    description "Comupte the veiling strategy value for the value <x>.";

    return evalb(LEM:-ExpressionCost(x) > LEM:-VeilingStrategy_maxcost);
  end proc: # VeilingStrategy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SetVeilingStrategyPars := proc({
    maxcost::{nonnegint}         := 15,
    subscripts::{nonnegint}      := 0,
    assignments::{nonnegint}     := 0,
    additions::{nonnegint}       := 1,
    multiplications::{nonnegint} := 2,
    divisions::{nonnegint}       := 3,
    functions::{nonnegint}       := 2
    }, $)::{nothing};

    description "Set the veiling strategy parameters: maximum veiling cost "
      "<maxcost>, subscripts cost weight parameter <subscripts>, assignments "
      "cost weight parameter <assignments>, additions cost weight parameter "
      "<additions>, multiplications cost weight parameter <multiplications>, "
      "divisions cost weight parameter <divisions>, and functions cost weight "
      "parameter <functions>.";

    LEM:-VeilingStrategy_maxcost         := maxcost;
    LEM:-VeilingStrategy_subscripts      := subscripts;
    LEM:-VeilingStrategy_assignments     := assignments;
    LEM:-VeilingStrategy_additions       := additions;
    LEM:-VeilingStrategy_multiplications := multiplications;
    LEM:-VeilingStrategy_divisions       := divisions;
    LEM:-VeilingStrategy_functions       := functions;
    return NULL;
  end proc: # SetVeilingStrategyPars

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UnVeil := proc(
    x::{anything},
    $)::{anything};

    description "Unveil the expression <x>.";

    local label;

    label := `if`(procname::{indexed}, op(procname), '_V');
    return eval['recurse'](x,LEM:-VeilUnorderedList(label));
  end proc: # UnVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UnVeilImap := proc(
    x::{anything},
    $)::{anything};

    description "Unveil the expression <x> with internal permutation map.";

    local label, T, perm, a, b, k;

	  label   := `if`(procname::{indexed}, op(procname), '_V');
    T, perm := LEM:-VeilTableImap(label, parse("reverse") = true);
    b       := copy(x);
    for k from 1 to nops(perm) do
      a := b;
      b := subs[eval](T[perm[k]], a);
    end do;
    return b;
  end proc: # UnVeilImap

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilList := proc(
    label::{symbol, list(symbol)} := VeilLabels(),
    {reverse::boolean := false},
    $)::list(anything);

    description "Return a list of the veiling variables labelled as <label>. "
      "If <label> is not given, return a list of all veiling variables.";

    local T, perm;

    if type(label, list) then
      return map(x -> op(LEM:-VeilList(x, parse("reverse") = reverse)), label);
    else
      T, perm := LEM:-VeilTableImap(label, parse("reverse") = reverse);
      return T[perm]:
    end if;
  end proc: # VeilList

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilUnorderedList := proc(
    label::{symbol},
    $)::{anything};

    description "Return a list of the veiling variables labelled as <label>.";

    if type(LEM:-UnVeilTables[label], table) then
      return op(eval(LEM:-UnVeilTables[label]));
    else
      return [];
    end if;
  end proc: # VeilUnorderedList

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilTableSize := proc(
    label::{symbol},
    $)::{nonnegint};

    description "Return the size of the table for symbol <label>.";

    return numelems(LEM:-VeilUnorderedList(label));
  end proc: # VeilTableSize

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilTableImap := proc(
    label::{symbol},
    {reverse::{boolean} := false},
    $)::{[], list(anything=anything)}, {[], list(nonnegint)};

    description "Return the table for symbol <label> and the permutation that "
      "sorts it.";

    local T, a, b, comparator;

    T := LEM:-VeilUnorderedList(label);
    if reverse then
      comparator := (a, b) -> evalb(op(1, lhs(a)) > op(1, lhs(b)));
    else
      comparator := (a, b) -> evalb(op(1, lhs(a)) < op(1, lhs(b)));
    end if;
    return T, sort(T, comparator, output = 'permutation');
  end proc: # VeilTableImap

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilTableAppend := proc(
    label::{symbol},
    x::{anything},
    $)::{anything};

    description "Append the veiled expression <x> to the veiling table with "
      "symbol <label>.";

    local k;

    k := 1;
    if type(LEM:-UnVeilTables[label], table) then
      k := numelems(op(eval(LEM:-UnVeilTables[label]))) + 1;
      LEM:-UnVeilTables[label][label[k]] := x;
    else
      LEM:-UnVeilTables[label] := table([label[1] = x]);
    end if;
    return label[k];
  end proc: # VeilTableAppend

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilLabels := proc( $ )::{list(symbol)};

    description "Return a list of the veiling labels.";

    return [indices(LEM:-UnVeilTables, 'nolist')];
  end proc: # VeilLabels

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilSubs := proc(
    x::{anything},
    label::{symbol, list(symbol)} := VeilLabels(),
    $)::{anything};

    description "Substitute the reversed veiling variables of the veiling label "
      "<label> in the expression <x>. If <label> is not given, substitute the "
      "reversed veiling variables of all veiling labels.";

    return subs[eval](op(LEM:-VeilList(label, parse("reverse") = true)), x);
  end proc: # VeilSubs

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilForget := proc(
    label::{symbol, list(symbol)} := VeilLabels(),
    $)::{nothing};

    description "Clear all the veiling variables of the veiling label <label>. "
      "If <label> is not given, clear all the veiling variables.";

    if type(label, list) then
      map(x -> LEM:-VeilForget(x), label);
    else
      LEM:-UnVeilLabels[label] := evaln(LEM:-UnVeilLabels[label]);
      LEM:-UnVeilTables[label] := evaln(LEM:-UnVeilTables[label]);
    end if;
    return NULL;
  end proc: # VeilForget

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # LEM

# That's all folks!
