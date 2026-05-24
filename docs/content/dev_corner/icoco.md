@page Dev_ICoco ICoCo - Interface for Code Coupling

ICoCo stands for Interface for Code Coupling. It is a coupling standard to couple codes exhibiting a time loop.

The main source of information for this is found at:

- at the root of the [ICoCo GitHub website](https://github.com/cea-trust-platform/icoco-coupling/tree/main)
- and more specifically in the file [ICoCoProblem.hxx](https://github.com/cea-trust-platform/icoco-coupling/blob/main/include/ICoCoProblem.hxx)

In TRUST, the file ProblemTrio.h and ProblemTrio.cpp represent the TRUST implementation of this norm.

## Implementation of getInputFieldTemplate and setInputField from the ICoCo API.

* Creation of input field classes.
  When reading the dataset (JDD), the input fields register themselves
  in their associated `Pb_base`, which keeps a list of references.
  They can then be retrieved using the `findInputField` method.

  The syntax is always:

  ```
  Ch... {  nb_comp number_of_field_components
           nom field_name
           probleme associated_base_problem_name
           [sous_zone associated_subzone_name] }
  ```

  The following classes derive from `Champ_Input_Proto`:

  * `Ch_front_input_uniforme`: uniform over a boundary
  * `Ch_front_input`: variable over a boundary, values defined on faces
  * `Champ_Input_P0`: cell-wise P0 field in the volume

  All implement the `getTemplate` and `setValue` methods,
  analogous to `getInputFieldTemplate` and `setInputField`, except
  that there is no need to identify the field.

To avoid diamond inheritance, `Champ_Input_Proto` is not
an `Objet_U`. As a result, unless a new kind of reference is defined
outside the Trio style, it is not possible to have references to
`Champ_Input_Proto`.

Therefore, references to `Field_base` are used instead, and type checks
must be performed in order to cast to `Champ_input_*`...
It is not very elegant, but I do not see another solution.

* Inheritance link removed between `Problem` and `Probleme_U`
  This allows different function signatures
  (`Problem` uses `string` and `vector`).

* Creation of an `ICoCo` directory containing everything related to the API:
  `Problem` classes, `ICoCoTrioField`, exceptions, and driver examples.

* No `#include <string>` outside the `.cpp` files in the `ICoCo` directory.
  However, there are forward declarations:

  ```
  #include <bits/stringfwd.h>
  ```

  which are equivalent to:

  ```
  class string;
  ```

  considering that `string` is a typedef.

* `Champ_base` and `Champ_front_base` now share a common ancestor,
  `Field_base`, in order to allow mixed lists of them.

* `allocation()` is no longer moved up to the `Probleme_U` level,
  but remains at the `Pb_base` level.

* Addition to `Probleme_U` of a `futureTime` method similar to
  `presentTime`.
  Together they define the resolution time interval.

* User boundary conditions:
  Previously, the user BC was read, the actual BC was created,
  then copied and destroyed.
  In the end, this was incompatible with `Champ_Input`, because the field
  that had registered itself with the problem was destroyed...
  => addition of an `adopt` method in `Cond_lim` to adopt an existing
  `Cond_lim_base`, thus avoiding the copy-destruction sequence.
  `adopt` should be used with caution... (see comments)

* Test cases:
  `U_in_var_impl_ICoCo.data` should produce the same result as
  `U_in_var_impl.data`; it uses `Ch_front_input` and
  `Ch_front_input_uniforme`.

  `ChDonXYZ_ICoCo.data` should produce the same result as
  `ChDonXYZ.data`; it uses `Champ_input_P0`.

  They were placed in the directories of their non-ICoCo counterparts.

  They use the classes `Pilote_ICoCo1` and `Pilote_ICoCo_2`, which
  demonstrate examples of using the API to drive Trio.


