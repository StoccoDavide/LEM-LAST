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

  export Veil,
         UnVeil,
         UnVeilImap,
         VeilUnorderedList,
         VeilList,
         VeilTableSize,
         VeilTableImap,
         VeilTableAppend,
         VeilLabels,
         VeilSubs,
         VeilForget;

  local  ModuleLoad,
         ModuleUnload,
         UnVeilTables,
         UnVeilLabels,
         Auxiliary,
         InitLEM;

  option package,
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

    LEM:-InitLEM();

    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ModuleUnload := proc()

    description "Module 'LEM' module unload procedure.";

    LEM:-UnVeilTables := NULL;

    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  InitLEM := proc()

    description "Initialize 'LEM' module internal variables.";

    LEM:-UnVeilTables := table([]);

    return NULL;
  end proc: # InitLEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Veil := proc(
    x::{anything},
    $)::{anything};

    description "Veil an expression <x> and return a label to it.";

    local label, i, s, c;

	  label := `if`(procname::{indexed}, op(procname), '_V');

	  if (label <> eval(label, 2)) then
	    error "LEM::Veil(...): label %a is assigned a value already, please save"
            " its contents and unassign it.", label;
	  end if;

    # Recognize zero if we can, so that we don't hide zeros.
    c := Normalizer(x);

    # Remove the integer content and sign so that we don't hide them either.
    i := icontent(c);

    # And we really mean sign, here, and not signum, because the interesting
    # case is when c may be a polynomial.
    try
      s := sign(c); # sign is weak
      # Don't do anything if we can tell that the coefficient is just a number
      # or a simple multiple of a name (simple or indexed)
      if (s*i = c) or type(s*c/i, 'name') then
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

  UnVeil := proc(
    x::{anything},
    $)::{anything};

    description "Unveil the expression <x>.";

    local label, T, a, b;

	  label := `if`(procname::{indexed}, op(procname), '_V');
    T := LEM:-VeilUnorderedList(label);
    b := copy(x);
    while has(b, label) do
      a := b;
      b := eval(a, T);
    end do;
    return b;
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UnVeilImap := proc(
    x::{anything},
    $)::{anything};

    description "Unveil the expression <x> with internal permutation map.";

    local label, T, perm, a, b, k;

	  label := `if`(procname::{indexed}, op(procname), '_V');
    T, perm := LEM:-VeilTableImap(label, true);
    b := copy(x);
    for k from 1 to nops(perm) do
      a := b;
      b := subs[eval](T[perm[k]], a);
    end do;
    return b;
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilList := proc(
    label::{symbol, list(symbol)} := VeilLabels(),
    {reverse::{boolean} := false},
    $)::{list(anything)};

    description "Return a list of the veiling variables labelled as <label>. "
      "If <label> is not given, return a list of all veiling variables.";

    local T, perm;

    if type(label, list) then
      return map(x -> op(LEM:-VeilList(x, reverse)), label);
    else
      T, perm := LEM:-VeilTableImap(label, reverse);
      return T[perm]:
    end if;
  end proc: # VeilList

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilUnorderedList := proc(
    label::{symbol},
    $)::{anything};

    description "Return a list of the veiling variables labelled as <label>.";

    if type(LEM:-UnVeilTables[label], 'table') then
      return op(eval(LEM:-UnVeilTables[label]));
    else
      return [];
    end if;
  end proc: # VeilUnorderedList;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilTableSize := proc(
    label::{symbol},
    $)::{nonnegint};

    description "Return the size of the table for symbol <label>.";

    return numelems(LEM:-VeilUnorderedList(label));
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilTableImap := proc(
    label::{symbol},
    {reverse::{boolean} := false},
    $)::{table, list[nonnegint]};

    description "Return the table for symbol <label> and the permutation that "
      "sorts it.";

    local T, a, b, comparator;

    T := LEM:-VeilUnorderedList(label);
    if reverse then
      comparator := (a, b) -> evalb(op(1, lhs(a)) > op(1, lhs(b)));
    else
      comparator := (a, b) -> evalb(op(1, lhs(a)) < op(1, lhs(b)));
    end if;
    return T, sort(T, output = 'permutation', comparator);
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilTableAppend := proc(
    label::{symbol},
    x::{anything},
    $)::{anything};

    description "Append the veiled expression <x> to the veiling table with "
      "symbol <label>.";

    local k;

    k := 1;
    if type(LEM:-UnVeilTables[label], 'table') then
      k := numelems(op(eval(LEM:-UnVeilTables[label]))) + 1;
      LEM:-UnVeilTables[label][label[k]] := x;
    else
      LEM:-UnVeilTables[label] := table([label[1] = x]);
    end if;
    return label[k];
  end proc: # Auxiliary;

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

    return subs[eval](op(LEM:-VeilList(label, true)), x);
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
