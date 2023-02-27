# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                              _     _____ __  __                             #
#                             | |   | ____|  \/  |                            #
#                             | |   |  _| | |\/| |                            #
#                             | |___| |___| |  | |                            #
#                             |_____|_____|_|  |_|                            #
#                         Large Expressions Management                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Authors: Davide Stocco (University of Trento)
#          Matteo Larcher (University of Trento)
#          Enrico Bertolazzi (University of Trento)
#          Wenqin Zhou (University of Western Ontario) - former affiliation
#          David J. Jeffrey (University of Western Ontario)
#          Jacques Carette (McMaster University)
#          Robert M. Corless (University of Western Ontario)
#
# License: BSD 3-Clause License
#
# This is a module for the 'LEM' (Large Expression Management) package. It
# contains the functions to deal with large expressions management. The module
# can be used to veil and unveil large expressions to avoid expression swell.
#
# The module is very similar to the 'LargeExpression' module provided by Maple.
# The main difference is that the 'LEM' module can handle more than one vailing
# label at a time. Moreover, it has some built-in functions to display, list and
# substitute the veiled expressions.
#
# The following code is hopefully an improved version of the version provided in
# the following PhD thesis:
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
         VeilLabels,
         ListVeil,
         SubsVeil,
         ForgetVeil,
         LastUsed;

  local  ModuleLoad,
         ModuleUnload,
         auxiliary,
         InitLEM,
         UnVeilTable,
         lib_base_path;

  option package,
         load   = ModuleLoad,
         unload = ModuleUnload;

  description "Large Expressions Management module";

  LastUsed := table('sparse');

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ModuleLoad := proc()

    description "'LEM' module load procedure";

    local i;

    printf(cat(
      "'LEM' module version 1.0, ",
      "BSD 3-Clause License - Copyright (C) 2023, D. Stocco, M. Larcher, ",
      "E. Bertolazzi, W. Zhou, D. J. Jeffrey, J. Carette and R. M. Corless.\n"
      ));

    lib_base_path := null;
    for i in [libname] do
      if (StringTools[Search]("LEM", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if (lib_base_path = null) then
      error "Cannot find 'LEM' module" ;
    end if;

    InitLEM();

    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ModuleUnload := proc()

    description "Module 'LEM' module unload procedure";

    unprotect(LastUsed);
    LastUsed := NULL;
    UnVeilTable := NULL;
    unprotect(LEM);
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  InitLEM := proc()

    description "Initialize 'LEM' module internal variables";

    UnVeilTable := table('sparse' = table('sparse' = (0 = 0)));
  end proc: # InitLEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Veil := proc(
    x::{anything},
    $)::{anything};

    description "Veil an expression <x> and return a label to it.";

    local i, s, c, label;

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
      if (s*i = c) or type(s*c/i,'name') then
        return c
      end if;
    catch:
      s := 1;
    end try;
    # Only if there is something complicated to hide we do actually hide it and
    # return a label.
    return s * i * auxiliary(s*c/i, label);
  end proc; # Veil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UnVeil := proc(
    x::{anything},
    n::{nonnegint, infinity},
    $)::{anything};

    description "UnVeil the expression <x> up to <n> levels.";

    local
      a, b, level, label;

    label := `if`(procname::{indexed}, op(procname), '_V');
    level := `if`(nargs < 2, 1, min(LastUsed[label]+1, n));

    # Always do at least 1 unveiling
    a := copy(x);
    b := eval(a, op(eval(UnVeilTable[label]))[2]);
    from 2 to level while not Testzero(a - b) do
      a := b;
      b := eval(a, op(eval(UnVeilTable[label]))[2]);
    end do;
    return b;
  end proc;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auxiliary := proc(
    x::{anything},
    label::{symbol},
    $)::{anything};

    description "Auxiliary procedure to scope LastUsed etc and use option "
      "remember to detect duplicates. There is no nontrivial storage duplication "
      "because objects are hashed and stored uniquely. All this costs is a "
      "pointer.";

    option remember;

    unprotect(LastUsed);
    LastUsed[label] := LastUsed[label] + 1;
    protect(LastUsed);
    if LastUsed[label] = 1 then
      UnVeilTable[label] := table('sparse' = (0 = 0));
      UnVeilTable[label][label[LastUsed[label]]] := x;
    else
      UnVeilTable[label][label[LastUsed[label]]] := x;
    end if;
    return label[LastUsed[label]]
  end proc: # auxiliary;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VeilLabels := proc(
    $)::{list(symbol)};

    description "Return a list of the veiling labels.";

    return map(lhs, op(op(LastUsed))[2..-1]);
  end proc: # ListVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ListVeil := proc(
    label::{symbol, list(symbol)} := VeilLabels(),
    $)::{list(anything)};

    description "Return a list of the veiling variables labelled as <label>. "
      "If <label> is not given, return a list of all veiling variables.";

    if type(label, list) then
      return map(x -> op(ListVeil(x)), label);
    else
      return sort(op(eval(UnVeilTable[label]))[2]):
    end if;
  end proc: # ListVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SubsVeil := proc(
    x::{anything},
    label::{symbol, list(symbol)} := VeilLabels(),
    $)::{anything};

    description "Substitute the reversed veiling variables of the veiling label "
      "<label> in the expression <x>. If <label> is not given, substitute the "
      "reversed veiling variables of all veiling labels.";

    return subs[eval](op(ListTools[Reverse](ListVeil(label))), x);
  end proc: # SubsVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ForgetVeil := proc(
    label::{symbol, list(symbol)} := VeilLabels(),
    $)::{nothing};

    description "Clear all the veiling variables of the veiling label <label>. "
      "If <label> is not given, clear all the veiling variables.";

    if type(label, list) then
      map(x -> ForgetVeil(x), label);
    else
      unprotect(LastUsed);
      LastUsed[label] := 0;
      protect(LastUsed);
      UnVeilTable[label] := evaln(UnVeilTable[label]);
      forget(auxiliary);
    end if;
  end proc: # ForgetVeil

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # LEM

# That's all folks!