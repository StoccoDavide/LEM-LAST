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
# The main difference is that the 'LEM' module has some built-in functions to
# collect, display, list and substitute the veiled expressions.
#
# The following code is hopefully an improved version of the original version
# provided inthe following PhD thesis:
#
#   Wenqin Zhou, Symbolic Computation Techniques for Solveing Large Expressions
#   Problems from Mathematics and Engineering (2007), Faculty of Graduate Studies,
#   The University of Western Ontario London, Ontario, Canada.
#
# We would like to thank Jacques Carette for providing the original code that
# have inspired this module.

unprotect('LEM');
module LEM()

  description "Large Expressions Management module.";

  option object;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local m_VeilingLabel                    := parse(
    cat("V_", StringTools:-Random(5, alnum))
  );
  local m_UnVeilTable                     := table([]);
  local m_VeilingStrategy_maxcost         := 15;
  local m_VeilingStrategy_subscripts      := 0;
  local m_VeilingStrategy_assignments     := 0;
  local m_VeilingStrategy_additions       := 1;
  local m_VeilingStrategy_multiplications := 2;
  local m_VeilingStrategy_divisions       := 3;
  local m_VeilingStrategy_functions       := 2;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Info::static := proc()

    description "Print 'LEM' module information.";

    printf(
      "+--------------------------------------------------------------------------+\n"
      "| 'LEM' module version 1.0 - BSD 3-Clause License - Copyright (c) 2023     |\n"
      "| Current version authors:                                                 |\n"
      "|   D. Stocco, M. Larcher, E. Bertolazzi.                                  |\n"
      "| Inspired by the work of:                                                 |\n"
      "|   W. Zhou, D. J. Jeffrey, J. Carette and R. M. Corless.                  |\n"
      "+--------------------------------------------------------------------------+\n"
    );
    return NULL;
  end proc: # Info

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleLoad::static := proc()

    description "'LEM' module load procedure.";

    local i, lib_base_path;

    lib_base_path := NULL;
    for i in [libname] do
      if (StringTools:-Search("LEM", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if (lib_base_path = NULL) then
      error "Cannot find 'LEM' module";
    end if;
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleUnload::static := proc()
    description "'LEM' module unload procedure.";
    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleCopy::static := proc(
    _self::LEM,
    proto::LEM,
    $)

    description "Copy the objects <proto> into <self>.";

    _self:-m_VeilingLabel                    := proto:-m_VeilingLabel;
    _self:-m_UnVeilTable                     := copy(proto:-m_UnVeilTable);
    _self:-m_VeilingStrategy_maxcost         := proto:-m_VeilingStrategy_maxcost;
    _self:-m_VeilingStrategy_subscripts      := proto:-m_VeilingStrategy_subscripts;
    _self:-m_VeilingStrategy_assignments     := proto:-m_VeilingStrategy_assignments;
    _self:-m_VeilingStrategy_additions       := proto:-m_VeilingStrategy_additions;
    _self:-m_VeilingStrategy_multiplications := proto:-m_VeilingStrategy_multiplications;
    _self:-m_VeilingStrategy_divisions       := proto:-m_VeilingStrategy_divisions;
    _self:-m_VeilingStrategy_functions       := proto:-m_VeilingStrategy_functions;
  end proc: # ModuleCopy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetVeilingLabel::static := proc(
    _self::LEM,
    label::{symbol, string},
    $)

    description "Set the veiling label to <label>.";

    if( _self:-VeilTableSize(_self) > 0) then
      error "the veiling veiling table is not empty, save the list if necessary "
        "and clear it before changing veiling label.";
      return NULL;
    end if;

    if type(label, string) then
      _self:-m_VeilingLabel := parse(label);
    else
      _self:-m_VeilingLabel := label;
    end if;
    return NULL;
  end proc: # SetVeilingLabel

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetVeilingLabel::static := proc(
    _self::LEM,
    $)::symbol;

    description "Return the veiling label.";

    return _self:-m_VeilingLabel;
  end proc: # GetVeilingLabel

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ClearVeilingTable::static := proc(
    _self::LEM,
    $)

    description "Clear the veiling table.";

    _self:-m_UnVeilTable := table([]);
    return NULL;
  end proc: # ClearVeilingTable

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Veil::static := proc(
    _self::LEM,
    x::anything,
    {force::boolean := false},
    $)::anything;

    description "Check if the veiling strategy is verified and veil an "
      "expression <x> and return a label to it.";

    local i, s, c;

    # Check if label is already assigned
	  if (_self:-m_VeilingLabel <> eval(_self:-m_VeilingLabel, 2)) then
	    error "LEM::Veil(...): label %a is already assigned, please save its "
        "contents and unassign it.", _self:-m_VeilingLabel;
	  end if;

    # Recognize zero if we can, so that we don't hide zeros.
    c := Normalizer(x);

    # Check the veiling strategy
    if not(force) and not(_self:-VeilingStrategy(_self, c)) then
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
    return s * i * _self:-VeilTableAppend(_self, s*c/i);
  end proc; # Veil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ExpressionCost::static := proc(
    _self::LEM,
    x::anything,
    $)::integer;

    description "Compute the cost of the expression <x>.";

    local tmp;

    if not(type(x, algebraic) or type(x, list)) then
      tmp := convert(x, list);
    else
      tmp := x;
    end if;

    return subs(
      'subscripts'      = _self:-m_VeilingStrategy_subscripts,
      'assignments'     = _self:-m_VeilingStrategy_assignments,
      'additions'       = _self:-m_VeilingStrategy_additions,
      'multiplications' = _self:-m_VeilingStrategy_multiplications,
      'divisions'       = _self:-m_VeilingStrategy_divisions,
      'functions'       = _self:-m_VeilingStrategy_functions,
      codegen:-cost(tmp)
    );
  end proc: # ExpressionCost

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilingStrategy::static := proc(
    _self::LEM,
    x::{algebraic},
    $)::{boolean};

    description "Comupte the veiling strategy value for the value <x>.";

    return evalb(_self:-ExpressionCost(_self, x) > _self:-m_VeilingStrategy_maxcost);
  end proc: # VeilingStrategy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetVeilingStrategyPars::static := proc(
    _self::LEM,
    {
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

    _self:-m_VeilingStrategy_maxcost         := maxcost;
    _self:-m_VeilingStrategy_subscripts      := subscripts;
    _self:-m_VeilingStrategy_assignments     := assignments;
    _self:-m_VeilingStrategy_additions       := additions;
    _self:-m_VeilingStrategy_multiplications := multiplications;
    _self:-m_VeilingStrategy_divisions       := divisions;
    _self:-m_VeilingStrategy_functions       := functions;
    return NULL;
  end proc: # SetVeilingStrategyPars

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export UnVeil::static := proc(
    _self::LEM,
    x::anything,
    $)::anything;

    description "Unveil the expression <x>.";

    return eval['recurse'](x, _self:-VeilTable(_self));
  end proc: # UnVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export UnVeilImap::static := proc(
    _self::LEM,
    x::anything,
    $)::anything;

    description "Unveil the expression <x> with internal permutation map.";

    local T, perm, a, b, k;

    T, perm := _self:-VeilTableImap(_self, parse("reverse") = true);
    b       := copy(x);
    for k from 1 to nops(perm) do
      a := b;
      b := subs[eval](T[perm[k]], a);
    end do;
    return b;
  end proc: # UnVeilImap

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilList::static := proc(
    _self::LEM,
    {reverse::boolean := false},
    $)::list(anything);

    description "Return a list of the veiling labels.";

    local T, perm;

    T, perm := _self:-VeilTableImap(_self, parse("reverse") = reverse);
    return T[perm]:
  end proc: # VeilList

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilUnorderedList::static := proc(
    _self::LEM,
    $)::list;

    description "Return the unordered veiling list.";

    if type(_self:-m_UnVeilTable, table) then
      return op(eval(_self:-m_UnVeilTable));
    else
      return [];
    end if;
  end proc: # VeilUnorderedList

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilTableSize::static := proc(
    _self::LEM,
    $)::nonnegint;

    description "Return the size of the internal veiling list.";

    return numelems(op(eval(_self:-m_UnVeilTable)));
  end proc: # VeilTableSize

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilTableImap::static := proc(
    _self::LEM,
    {reverse::boolean := false},
    $)::{[], list(anything = anything)}, {[], list(nonnegint)};

    description "Return the veiling list and the permutation that sorts it.";

    local T, a, b, comparator;

    T := _self:-VeilUnorderedList(_self);
    if reverse then
      comparator := (a, b) -> evalb(op(1, lhs(a)) > op(1, lhs(b)));
    else
      comparator := (a, b) -> evalb(op(1, lhs(a)) < op(1, lhs(b)));
    end if;
    return T, sort(T, comparator, output = 'permutation');
  end proc: # VeilTableImap

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilTableAppend::static := proc(
    _self::LEM,
    x::anything,
    $)::indexed;

    description "Append the veiled expression <x> to the veiling list.";

    local k;

    k := 1;
    if type(_self:-m_UnVeilTable, table) then
      k := numelems(op(eval(_self:-m_UnVeilTable))) + 1;
      _self:-m_UnVeilTable[_self:-m_VeilingLabel[k]] := x;
    else
      _self:-m_UnVeilTable := table([_self:-m_VeilingLabel[1] = x]);
    end if;
    return _self:-m_VeilingLabel[k];
  end proc: # VeilTableAppend

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilSubs::static := proc(
    _self::LEM,
    x::anything,
    $)::anything;

    description "Substitute the reversed veiling variables of the internal "
      "veiling label in the expression <x>.";

    return subs[eval](op(_self:-VeilList(_self, parse("reverse") = true)), x);
  end proc: # VeilSubs

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilForget::static := proc(
    _self::LEM,
    $)

    description "Clear all the veiling variables of the internal veiling label.";

    _self:-m_UnVeilTable := evaln(_self:-m_UnVeilTable);
    return NULL;
  end proc: # VeilForget

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # LEM

# That's all folks!
