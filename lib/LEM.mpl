# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                              _     _____ __  __                             #
#                             | |   | ____|  \/  |                            #
#                             | |   |  _| | |\/| |                            #
#                             | |___| |___| |  | |                            #
#                             |_____|_____|_|  |_|                            #
#                         Large Expressions Management                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco     (University of Trento)
#   Matteo Larcher    (University of Trento)
#   Enrico Bertolazzi (University of Trento)
#
# Inspired by the work of:
#   Wenqin Zhou       (University of Western Ontario) - Former affiliation
#   David J. Jeffrey  (University of Western Ontario)
#   Jacques Carette   (McMaster University)
#   Robert M. Corless (University of Western Ontario)
#
# License: BSD 3-Clause License
#
# This is a module for the 'LEM' (Large Expressions Management) module.
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

  # General
  local m_VerboseMode := false;
  local m_WarningMode := true;

  # Veiling
  local m_VeilingLabel := parse(cat("V_", StringTools:-Random(5, 'alnum')));
  local m_UnveilTable                     := table([]);
  local m_VeilingDependency               := [];
  local m_VeilingStrategy_maxcost         := 15;
  local m_VeilingStrategy_subscripts      := 2;
  local m_VeilingStrategy_assignments     := 0;
  local m_VeilingStrategy_additions       := 1;
  local m_VeilingStrategy_multiplications := 2;
  local m_VeilingStrategy_divisions       := 3;
  local m_VeilingStrategy_functions       := 2;

  # Signature
  local m_SigMode  := true;
  local m_SigTable := table([]);
  local m_SigValue := 1000000007;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Info::static := proc()

    description "Print 'LEM' module information.";

    printf(
      "+--------------------------------------------------------------------------+\n"
      "| 'LEM' module version 1.0 - BSD 3-Clause License - Copyright (c) 2023     |\n"
      "| Current version authors:                                                 |\n"
      "|   D. Stocco, M. Larcher and E. Bertolazzi.                               |\n"
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
      error "cannot find 'LEM' module.";
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

    # General
    _self:-m_VerboseMode := proto:-m_VerboseMode;
    _self:-m_WarningMode := proto:-m_WarningMode;

    # Veiling
    _self:-m_VeilingLabel                    := proto:-m_VeilingLabel;
    _self:-m_UnveilTable                     := copy(proto:-m_UnveilTable);
    _self:-m_VeilingDependency               := proto:-m_VeilingDependency;
    _self:-m_VeilingStrategy_maxcost         := proto:-m_VeilingStrategy_maxcost;
    _self:-m_VeilingStrategy_subscripts      := proto:-m_VeilingStrategy_subscripts;
    _self:-m_VeilingStrategy_assignments     := proto:-m_VeilingStrategy_assignments;
    _self:-m_VeilingStrategy_additions       := proto:-m_VeilingStrategy_additions;
    _self:-m_VeilingStrategy_multiplications := proto:-m_VeilingStrategy_multiplications;
    _self:-m_VeilingStrategy_divisions       := proto:-m_VeilingStrategy_divisions;
    _self:-m_VeilingStrategy_functions       := proto:-m_VeilingStrategy_functions;

    # Signature
    _self:-m_SigMode  := proto:-m_SigMode;
    _self:-m_SigTable := copy(proto:-m_SigTable);
    _self:-m_SigValue := proto:-m_SigValue;
  end proc: # ModuleCopy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableVerboseMode::static := proc(
    _self::LEM,
    $)

    description "Enable the verbosity of the module.";

    _self:-m_VerboseMode := true;
    return NULL;
  end proc: # EnableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableVerboseMode::static := proc(
    _self::LEM,
    $)

    description "Disable the verbosity of the module.";

    _self:-m_VerboseMode := false;
    return NULL;
  end proc: # DisableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetVerboseMode::static := proc(
    _self::LEM,
    mode::boolean,
    $)

    description "Set the verbosity of the module to <mode>.";

    _self:-m_VerboseMode := mode;
    return NULL;
  end proc: # SetVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableWarningMode::static := proc(
    _self::LEM,
    $)

    description "Enable the warning mode of the module.";

    _self:-m_WarningMode := true;
    return NULL;
  end proc: # EnableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableWarningMode::static := proc(
    _self::LEM,
    $)

    description "Disable the warning mode of the module.";

    _self:-m_WarningMode := false;
    return NULL;
  end proc: # DisableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetWarningMode::static := proc(
    _self::LEM,
    mode::boolean,
    $)

    description "Set the warning mode of the module to <mode>.";

    _self:-m_WarningMode := mode;
    return NULL;
  end proc: # SetWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetVeilingLabel::static := proc(
    _self::LEM,
    label::{symbol, string},
    $)

    description "Set the veiling label to <label>.";

    if (_self:-VeilTableSize(_self) > 0) then
      error "the veiling table is not empty, save the list if necessary and "
        "clear it before changing veiling label.";
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

  export SetVeilingDependency::static := proc(
    _self::LEM,
    dependency::list,
    $)

    description "Set the veiling dependency to <dependency>.";

    if (_self:-VeilTableSize(_self) > 0) then
      error "the veiling table is not empty, save the list if necessary and "
        "clear it before changing veiling dependency.";
      return NULL;
    end if;

    _self:-m_VeilingDependency := dependency;
    return NULL;
  end proc: # SetVeilingDependency

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetVeilingDependency::static := proc(
    _self::LEM,
    $)::list;

    description "Return the veiling dependency.";

    return _self:-m_VeilingDependency;
  end proc: # GetVeilingDependency

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ClearVeilingDependency::static := proc(
    _self::LEM,
    $)

    description "Clear the veiling dependency.";

    _self:-m_VeilingDependency := {};
    return NULL;
  end proc: # ClearVeilingDependency

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableSignatureMode::static := proc(
    _self::LEM,
    $)

    description "Enable the expression signature calculation.";

    if (_self:-VeilTableSize(_self) > 0) then
      error "the signature table is not empty, save the list if necessary and "
        "clear it before enabling signature calculation.";
      return NULL;
    end if;

    _self:-m_SigMode := true;
    return NULL;
  end proc: # EnableSignatureMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableSignatureMode::static := proc(
    _self::LEM,
    $)

    description "Disable the expression signature calculation.";

    if (_self:-VeilTableSize(_self) > 0) then
      error "the signature table is not empty, save the list if necessary and "
        "clear it before disabling signature calculation.";
      return NULL;
    end if;

    _self:-m_SigMode := false;
    return NULL;
  end proc: # DisableSignatureMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetSignatureMode::static := proc(
    _self::LEM,
    mode::boolean,
    $)

    description "Set the expression signature calculation of the module to "
      "<mode>.";

    if not mode and (_self:-VeilTableSize(_self) > 0) then
      error "the signature table is not empty, save the list if necessary and "
        "clear it before disabling signature calculation.";
      return NULL;
    end if;

    _self:-m_SigMode := mode;
    return NULL;
  end proc: # SetSignatureMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetSignatureValue := proc(
    _self::LEM,
    p::prime,
    $)

    description "Set the signature value to <p>.";

    _self:-m_SigValue := p;
    return NULL;
  end proc: # SetSignatureValue

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetSignatureValue := proc(
    _self::LEM,
    $)::prime;

    description "Return the signature value.";

    return _self:-m_SigValue;
  end proc: # GetSignatureValue

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetVeilingLabel::static := proc(
    _self::LEM,
    $)::symbol;

    description "Return the veiling label.";

    return _self:-m_VeilingLabel;
  end proc: # GetVeilingLabel

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Veil::static := proc(
    _self::LEM,
    x::anything,
    {
    depend::list   := [],
    force::boolean := false
    },
    $)::anything;

    description "Check if the veiling strategy is verified, if true veil the "
      "expression <x> and return a label to it.";

    local i, s, c, deps;
    # Check if label is already assigned
	  if (_self:-m_VeilingLabel <> eval(_self:-m_VeilingLabel, 2)) then
	    error "label %a is already assigned, please save its contents and "
        "unassign it.", _self:-m_VeilingLabel;
	  end if;

    # Recognize zero if we can, so that we don't hide zeros.
    c := Normalizer(x);

    # Check the veiling strategy
    if not force and not _self:-VeilingStrategy(_self, c) then
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
    if (nops(_self:-m_VeilingDependency) > 0) then
      convert(_self:-m_VeilingDependency, set) intersect indets(s*c/i);
      deps := remove[flatten](j -> evalb(j in %), _self:-m_VeilingDependency);
      if (nops(deps) > 0) then
        return s * i * _self:-TablesAppend(_self, s*c/i)(op(deps));
      end if;
    end if;
    return s * i * _self:-TablesAppend(_self, s*c/i);
  end proc; # Veil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ExpressionCost::static := proc(
    _self::LEM,
    x::anything,
    $)::integer;

    description "Compute the cost of the expression <x>.";

    local tmp;

    if not (type(x, algebraic) or type(x, list)) then
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
    x::algebraic,
    $)::boolean;

    description "Evaluate the veiling strategy for the expression <x>.";

    return evalb(_self:-ExpressionCost(_self, x) > _self:-m_VeilingStrategy_maxcost);
  end proc: # VeilingStrategy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetVeilingStrategyPars::static := proc(
    _self::LEM,
    {
    maxcost::nonnegint         := 15,
    subscripts::nonnegint      := 0,
    assignments::nonnegint     := 0,
    additions::nonnegint       := 1,
    multiplications::nonnegint := 2,
    divisions::nonnegint       := 3,
    functions::nonnegint       := 2
    }, $)

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

  export Unveil::static := proc(
    _self::LEM,
    x::anything,
    $)::anything;

    description "Unveil the expression <x> with the veiled variables.";

    return eval['recurse'](x, _self:-VeilList(_self));
  end proc: # Unveil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export UnveilSig::static := proc(
    _self::LEM,
    x::anything,
    $)::anything;

    description "Unveil the expression <x> with signature values.";

    return eval['recurse'](x, _self:-SigList(_self));
  end proc: # UnveilSig

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export UnveilImap::static := proc(
    _self::LEM,
    x::anything,
    $)::anything;

    description "Unveil the expression <x> with veiling labels permutation "
      "map.";

    local T, perm, a, b, k;

    T, perm := _self:-VeilTableImap(_self, parse("reverse") = true);
    b       := copy(x);
    for k from 1 to nops(perm) do
      a := b;
      b := subs[eval](T[perm[k]], a);
    end do;
    return b;
  end proc: # UnveilImap

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SigImap::static := proc(
    _self::LEM,
    x::anything,
    $)::anything;

    description "Unveil the expression <x> with signature values permutation "
      "map.";

    local T, perm, a, b, k;

    T, perm := _self:-SigTableImap(_self, parse("reverse") = true);
    b       := copy(x);
    for k from 1 to nops(perm) do
      a := b;
      b := subs[eval](T[perm[k]], a);
    end do;
    return b;
  end proc: # SigImap

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilDependencyList::static := proc(
    _self::LEM,
    {
    reverse::boolean := false
    }, $)::list(anything);

    description "Return a list of the veiling labels dependency substitution.";

    local vars, veil, i, deps;

    # Extract variables and veils
    vars := _self:-m_VeilingDependency;
    veil := _self:-VeilList(_self, parse("reverse") = reverse, parse("dependency") = false);

    # Iterate over veils to find dependencies
    for i from 1 to nops(veil) do
      convert(_self:-m_VeilingDependency, set) intersect indets(rhs(veil[i]));
      deps := remove[flatten](j -> evalb(j in %), _self:-m_VeilingDependency);
      veil[i]  := lhs(veil[i]) = `if`(
        nops(deps) > 0, lhs(veil[i])(op(deps)), lhs(veil[i])
      );
    end do;
    return veil;
  end proc: # VeilDependencyList

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilList::static := proc(
    _self::LEM,
    {
    reverse::boolean    := false,
    dependency::boolean := false
    }, $)::list(anything);

    description "Return a list of the veiling labels.";

    local T, perm;

    T, perm := _self:-VeilTableImap(_self, parse("reverse") = reverse);
    if (dependency) then
      T := subs(op(
        _self:-VeilDependencyList(_self, parse("reverse") = reverse)
        ), lhs~(T)
      ) =~ rhs~(T);
    end if;
    return T[perm]:
  end proc: # VeilList

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SigList::static := proc(
    _self::LEM,
    {reverse::boolean := false},
    $)::list(anything);

    description "Return a list of the signature labels.";

    local T, perm;

    T, perm := _self:-SigTableImap(_self, parse("reverse") = reverse);
    return T[perm]:
  end proc: # SigList

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilUnorderedList::static := proc(
    _self::LEM,
    $)::list;

    description "Return the unordered veiling list.";

    if type(_self:-m_UnveilTable, table) then
      return op(eval(_self:-m_UnveilTable));
    else
      return [];
    end if;
  end proc: # VeilUnorderedList

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SigUnorderedList::static := proc(
    _self::LEM,
    $)::list;

    description "Return the unordered signature list.";

    if type(_self:-m_SigTable, table) then
      return op(eval(_self:-m_SigTable));
    else
      return [];
    end if;
  end proc: # SigUnorderedList

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilTableSize::static := proc(
    _self::LEM,
    $)::nonnegint;

    description "Return the size of the internal veiling table.";

    return numelems(_self:-VeilUnorderedList(_self));
  end proc: # VeilTableSize

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SigTableSize::static := proc(
    _self::LEM,
    $)::nonnegint;

    description "Return the size of the internal signature table.";

    return numelems(_self:-SigUnorderedList(_self));
  end proc: # SigTableSize

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

  export SigTableImap::static := proc(
    _self::LEM,
    {reverse::boolean := false},
    $)::{[], list(anything = anything)}, {[], list(nonnegint)};

    description "Return the signature list and the permutation that sorts it.";

    local T, a, b, comparator;

    T := _self:-SigUnorderedList(_self);
    if reverse then
      comparator := (a, b) -> evalb(op(1, lhs(a)) > op(1, lhs(b)));
    else
      comparator := (a, b) -> evalb(op(1, lhs(a)) < op(1, lhs(b)));
    end if;
    return T, sort(T, comparator, output = 'permutation');
  end proc: # SigTableImap

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export TablesAppend::static := proc(
    _self::LEM,
    x::anything,
    {iter::nonnegint := 10},
    $)::indexed;

    description "Append the veiled expression <x> to the veiling table.";

    local k;

    k := 1;
    if type(_self:-m_UnveilTable, table) then
      k := numelems(_self:-m_UnveilTable) + 1;
      _self:-m_UnveilTable[_self:-m_VeilingLabel[k]] := x;
      if _self:-m_SigMode then
        _self:-m_SigTable[_self:-m_VeilingLabel[k]] := SIG(
          _self:-SubsSig(_self, x), _self:-m_SigValue,
          parse("iter")    = iter,
          parse("verbose") = _self:-m_WarningMode
        );
      else
        _self:-m_SigTable[_self:-m_VeilingLabel[k]] := 0;
      end if;
    else
      _self:-m_UnveilTable := table([_self:-m_VeilingLabel[1] = x]);
      if _self:-m_SigMode then
        _self:-m_SigTable   := table([_self:-m_VeilingLabel[1] = SIG(
          x, _self:-m_SigValue,
          parse("iter")    = iter,
          parse("verbose") = _self:-m_WarningMode
        )]);
      else
        _self:-m_SigTable   := table([_self:-m_VeilingLabel[1] = 0]);
      end if;
    end if;
    return _self:-m_VeilingLabel[k];
  end proc: # TablesAppend

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SubsSig::static := proc(
    _self::LEM,
    x::anything,
    $)::anything;

    description "Substitute the reversed signature values of the internal "
      "signature table in the expression <x>.";

    local out;
    try
      out := subs[eval](op(_self:-SigList(_self, parse("reverse") = true)), x);
    catch "division by zero":
      if _self:-m_WarningMode then
        WARNING("LEM:-SubsSig(...): division by zero, returning 'FAIL'.");
      end if;
      out := FAIL;
    catch:
      if _self:-m_WarningMode then
        WARNING(
          "LEM:-SubsSig(...): something went wrong in '%1', returning 'FAIL'.",
          lastexception
        );
      end if;
      out := FAIL;
    end try;
    return out;
  end proc: # SubsSig

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SubsVeil::static := proc(
    _self::LEM,
    x::anything,
    {
    reverse::boolean    := true,
    dependency::boolean := false
    }, $)::anything;

    description "Substitute the reversed veiling variables of the internal "
      "veiling table in the expression <x>, If dependency flag <dependency> "
      "is set to true, the veiling variables are substituted with their "
      "dependencies.";

    return subs[eval](op(_self:-VeilList(_self,
        parse("dependency") = dependency,
        parse("reverse")    = reverse
      )), x);
  end proc: # SubsVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ForgetVeil::static := proc(
    _self::LEM,
    $)

    description "Clear the internal veiling and signature tables.";

    _self:-m_UnveilTable := table([]);
    _self:-m_SigTable    := table([]);
    return NULL;
  end proc: # ForgetVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Signature::static := proc(
    _self::LEM,
    x::anything,
    p::prime := _self:-m_SigValue,
    {iter::nonnegint := 10},
    $)::anything;

    description "Compute the signature of the expression <x> modulo <p> (the "
      "default is the internal signature value) also by using the internal "
      "signature table of the veiled expressions.";

    return SIG(
      _self:-SubsSig(_self, x),
      _self:-m_SigValue,
      parse("iter")    = iter,
      parse("verbose") = _self:-m_WarningMode
    );
  end proc: # Signature

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # LEM

# That's all folks!
