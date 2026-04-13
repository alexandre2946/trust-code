# User's documentation

Let's pretend you want to document keyword `the_kw` which inherits from `base_kw`. Imagine also that this keyword takes an attribute `tutu` which is a double.

The C++ source code implementing this keyword should be annotated with `// XD` tags as follows:


```
// XD  the_kw  base_kw  synonym_for_the_kw  with_brace  description
...
param.ajouter("tutu", &tutu);  // XD_ADD_P double A description for tutu attribute.

```

The `with_brace` flag can take one of the following values:
- 'INHERITS_BRACE' - the rule will be the same as the the parent class
- 'NO_BRACE' keyword does not expect braces when read in the dataset. Names of the attributes are not used for reading or writing.  Example : `Champ_Uniforme 3 0. 0. 0.`
- 'BRACE' keyword expects curly braces, the name of the attributes is explicitly provided. Example: `Lire sch { tinit 0. tmax 0. }`

It is important to place the attribute tags on the same line as the `param.ajouter` calls so that the name of the parameter used in the C++ code can be extracted.

The allowed types for attributes are:
- `int` a single integer value
- `double`a single floatting value
- `chaine` either a single word (i.e. a string without any space) **or** (this is weird!!) a full block encompassed in curly braces. Example: `toto` is a chaine, and `{ how do you { do } }` is also a `chaine`, but `titi toto` is not.
- `rien`: a void item - this type is usually used for flags, where just the presence/absence of the attribute suffices to indicate whether it is turned on or off


If you need to write really long lines, you can also use the `// XD_CONT` tag to continue your comment. This works with both the simple `XD` tag and the `XD_ADD_P` tag.
For example:

```
// XD  the_kw  base_kw  synonym_for_the_kw  with_brace A really long and verbose description, that you need to write on 
// XD_CONT several lines because you want to be expressive and this requires
// XD_CONT more space etc etc
...
param.ajouter("tutu", &tutu);  // XD_ADD_P double A also very long description for tutu attribute which also
                               // XD_CONT requires writing stuff on several lines 

```


Finally, note that you may use `// 2XD` or even `// 3XD ` tags if you want to nest keyword declarations, like in this example:

```
  Param param("lire_ceg");
  param.ajouter("frontiere",&la_surface_libre_nom_,Param::REQUIRED); // XD_ADD_P chaine  To specify the boundaries conditions representing the free surfaces
  param.ajouter("t_deb",&t_deb_,Param::REQUIRED); // XD_ADD_P double value of the CEG's initial calculation time
  param.ajouter("t_fin",&t_fin_); // XD_ADD_P double not_set time during which the CEG's calculation was stopped
  param.ajouter("dt_post",&dt_post_); // XD_ADD_P double periode refers to the printing period, this value is expressed in seconds
  param.ajouter("haspi",&haspi_,Param::REQUIRED); // XD_ADD_P double The suction height  required to calculate AREVA's criterion
  param.ajouter("debug",&debug_); // XD_ADD_P int not_set
  Param& param_cea_jaea =param.ajouter_param("CEA_JAEA"); // XD_ADD_P ceg_cea_jaea CEA_JAEA's criterion
  // 2XD ceg_cea_jaea objet_lecture nul -1 not_set
  param_cea_jaea.ajouter("normalise",&critere_cea_jaea_normalise_);   // 2XD_ADD_P int renormalize (1) or not (0) values alpha and gamma
  param_cea_jaea.ajouter("nb_mailles_mini",&nb_mailles_mini_);  // 2XD_ADD_P int Sets the minimum number of cells for the detection of a vortex.
  param_cea_jaea.ajouter("min_critere_Q_sur_max_critere_Q",&min_critere_Q_sur_max_critere_Q_);  // 2XD_ADD_P double Is an optional keyword used to correct the minimum values of Q's  criterion taken into account in the detection of a vortex

```

In this case the newly introduced keywords will be declared after the main one. Maximum 3 levels.
