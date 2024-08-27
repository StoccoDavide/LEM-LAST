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

  local m_VerboseMode       := false;
  local m_WarningMode       := true;
  local m_VeilingLabel      := parse(cat("V_", StringTools:-Random(5, 'alnum')));
  local m_UnveilTable       := table([]);
  local m_VeilingDependency := [];
  local m_ExprMaxCost       := 500;

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

    _self:-m_VerboseMode       := proto:-m_VerboseMode;
    _self:-m_WarningMode       := proto:-m_WarningMode;
    _self:-m_VeilingLabel      := proto:-m_VeilingLabel;
    _self:-m_UnveilTable       := copy(proto:-m_UnveilTable);
    _self:-m_VeilingDependency := proto:-m_VeilingDependency;
    _self:-m_ExprMaxCost       := proto:-m_ExprMaxCost;
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

    # Recognize zero if we can, so that we don't hide zeros
    c := Normalizer(x);

    # Check the veiling strategy
    if not force and not _self:-VeilingStrategy(_self, c) then
      return c;
    end if;

    # Remove the integer content and sign so that we don't hide them either
    i := icontent(c);

    # And we really mean sign, here, and not signum, because the interesting
    # case is when c may be a polynomial
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
    # Only if there is something complicated to hide we do actually hide it
    if (nops(_self:-m_VeilingDependency) > 0) then
      convert(_self:-m_VeilingDependency, set) intersect indets(s*c/i);
      deps := select[flatten](j -> evalb(j in %), _self:-m_VeilingDependency);
      if (nops(deps) > 0) then
        return s*i*_self:-VeilTableAppend(_self, s*c/i)(op(deps));
      end if;
    end if;
    return s*i*_self:-VeilTableAppend(_self, s*c/i);
  end proc; # Veil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ExpressionCost::static := proc(
    _self::LEM,
    x::anything,
    {
    dependency::boolean := true
    }, $)::integer;

    description "Compute the cost of the expression <x>. If dependency flag "
      "<dependency> is set to true, the veiling variables are considered with "
      "their dependencies (e.g., V(x,y)), otherwise they are considered as "
      "independent variables (e.g., V).";

    local tmp;

    if not (type(x, algebraic) or type(x, list)) then
      tmp := convert(x, list);
    else
      tmp := x;
    end if;

    if not dependency then
      _self:-VeilDependencyList(_self, parse("reverse") = true);
      tmp := subs(op(rhs~(%) =~ lhs~(%)), tmp);
    end if;

    return _self:-GraphComplexity(_self, tmp);
  end proc: # ExpressionCost

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GraphStatistics := proc(
    _self::LEM,
    expr::anything,
    $)::nonnegint, nonnegint, nonnegint;

    description "Given the expression <expr>, return the directed acyclic "
    "graph's edges, nodes, and leafs.";

    local i, edges, tmp_edges, nodes, tmp_nodes, leafs, tmp_leafs;

    edges := 0;
    nodes := 0;
    leafs := 0;
    if (nops(expr) <> 1) then
      nodes := nodes + 1;
      edges := edges + nops(expr);
      for i in op(expr) do
        tmp_edges, tmp_nodes, tmp_leafs := GraphStatistics(_self, i);
        edges := edges + tmp_edges;
        nodes := nodes + tmp_nodes;
        leafs := leafs + tmp_leafs;
      end do;
    else
      leafs := leafs + 1;
    end if;
    return edges, nodes, leafs;
  end proc: # GraphStatistics

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GraphComplexity := proc(
    _self::LEM,
    expr::anything,
    level::nonnegint := 0,
    $)::nonnegint;

    description "Given the expression <expr> and the starting subtree's level "
      "<level>, return the complexity of the expression.";

    local i, out;

    out := level + 1;
    if (nops(expr) <> 1) then
      for i in op(expr) do
        out := out + GraphComplexity(_self, i, level + 1);
      end do;
    end if;
    return out;
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilingStrategy::static := proc(
    _self::LEM,
    x::algebraic,
    $)::boolean;

    description "Evaluate the veiling strategy for the expression <x>.";

    return evalb(_self:-ExpressionCost(_self, x) > _self:-m_ExprMaxCost);
  end proc: # VeilingStrategy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetExprMaxCost::static := proc(
    _self::LEM,
    maxcost::nonnegint,
    $)

    description "Set the maximum expression cost <maxcost> for veiling.";

    _self:-m_ExprMaxCost := maxcost;
    return NULL;
  end proc: # SetExprMaxCost

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Unveil::static := proc(
    _self::LEM,
    x::anything,
    $)::anything;

    description "Unveil the expression <x> with the veiled variables.";

    return subs(op(_self:-VeilList(_self, parse("reverse") = true)), x);
  end proc: # Unveil

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

  export VeilDependencyList::static := proc(
    _self::LEM,
    {
    reverse::boolean := false
    }, $)::list(anything);

    description "Return a list of the veiling labels dependency substitution.";

    local veil, i, deps;

    # Extract variables and veils
    veil := Array(_self:-VeilList(
      _self, parse("reverse") = reverse, parse("dependency") = false
    ));

    # Iterate over veils to find dependencies
    for i from 1 to ArrayNumElems(veil) do
      convert(_self:-m_VeilingDependency, set) intersect indets(rhs(veil[i]));
      deps := select[flatten](j -> evalb(j in %), _self:-m_VeilingDependency);
      veil[i] := lhs(veil[i]) = `if`(
        nops(deps) > 0, lhs(veil[i])(op(deps)), lhs(veil[i])
      );
    end do;
    return convert(veil, list);
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
    if dependency then
      T := subs(op(
        _self:-VeilDependencyList(_self, parse("reverse") = reverse)
        ), lhs~(T)
      ) =~ rhs~(T);
    end if;
    return T[perm]:
  end proc: # VeilList

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

  export VeilTableSize::static := proc(
    _self::LEM,
    $)::nonnegint;

    description "Return the size of the internal veiling table.";

    return numelems(_self:-VeilUnorderedList(_self));
  end proc: # VeilTableSize

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VeilTableImap::static := proc(
    _self::LEM,
    {
    reverse::boolean := false
    }, $)::{[], list(anything = anything)}, {[], list(nonnegint)};

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

    description "Append the veiled expression <x> to the veiling table.";

    local k;

    k := 1;
    if type(_self:-m_UnveilTable, table) then
      k := numelems(_self:-m_UnveilTable) + 1;
      _self:-m_UnveilTable[_self:-m_VeilingLabel[k]] := x;
    else
      _self:-m_UnveilTable := table([_self:-m_VeilingLabel[1] = x]);
    end if;
    return _self:-m_VeilingLabel[k];
  end proc: # VeilTableAppend

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SubsVeil::static := proc(
    _self::LEM,
    x::anything,
    {
    reverse::boolean    := true,
    dependency::boolean := false
    }, $)::anything;

    description "Substitute the reversed veiling variables of the internal "
      "veiling table in the expression <x>. If dependency flag <dependency> "
      "is set to true, the veiling variables are substituted with their "
      "dependencies.";

    return subs[eval](op(_self:-VeilList(_self,
        parse("reverse")    = reverse,
        parse("dependency") = dependency
      )), x);
  end proc: # SubsVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ForgetVeil::static := proc(
    _self::LEM,
    $)

    description "Clear the internal veiling table.";

    _self:-m_UnveilTable := table([]);
    return NULL;
  end proc: # ForgetVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # LEM

# That's all folks!
