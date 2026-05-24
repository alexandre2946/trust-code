@page DevTut_BaltikCD Build a BALTIK application — convection-diffusion tutorial

## Abstract

This tutorial describes the steps to follow when building a `BALTIK` application on top of the **TRUST** software.

To this end, a convection-diffusion equation is chosen as a working example of a mathematical problem to solve. In this equation, the convection velocity and the diffusion coefficient are user-provided inputs.

The development specifications are described and some recommendations are given for solving this convection-diffusion equation. The developments are then carried out within an environment provided by the **TRUST** software called a `BALTIK` application, whose setup is also described.

Finally, two verification and non-regression tests are performed using the tools provided by the **TRUST** software (`genererCourbes`, `make check_optim`).

---

## Introduction

This tutorial aims to describe the steps to follow when building a `BALTIK` application on top of the **TRUST** software.
To this end, the following convection-diffusion equation is chosen as a working example:

\f[
\frac{\partial c}{\partial t}+\nabla\cdot(\overrightarrow{V}c)+\nabla\cdot(D\nabla c)=0 \quad \textrm{on } \Omega=[a,b]\times[c,d]
\f]

with periodic boundary conditions, where:

- \f$c(\overrightarrow{X},t)\f$: quantity of interest (species concentration, temperature, ...),
- \f$\overrightarrow{V}(\overrightarrow{X},t)\f$: convection velocity of \f$c\f$ (user input),
- \f$D(\overrightarrow{X},t)\f$: diffusion coefficient of \f$c\f$ (user input),
- \f$a\f$, \f$b\f$, \f$c\f$ and \f$d\f$: bounds of the two-dimensional computational domain.

The spatial and temporal discretization of the convection-diffusion equation will be performed using the numerical schemes available in the **TRUST** software.

In the **Specification** section, we describe the development specifications and give some recommendations for solving this equation.

The developments specified in the **Specification** section are described in the **Developments** section. These developments will be carried out within a `BALTIK` application environment, whose setup is also described in that section.

The **Development Verification Process** section is dedicated to two verification and non-regression tests using the tools provided by the **TRUST** software (`genererCourbes`, `make check_optim`).

---

## Specification

### Concepts

- **Problem concept**: the role of a problem is to solve, over a domain, the equations that make up the problem. A time discretization scheme and a space discretization scheme are associated with it.

- **Equation concept**: the role of an equation is to compute one or more fields given the following choices:
  - a time discretization scheme,
  - a space discretization scheme,
  - boundary conditions,
  - source terms and operators.

  An equation is carried by a problem and holds a back-reference to the problem that owns it.

- **Operator concept**: operators are parts of an equation. Among the most commonly used operators are the convection and diffusion operators.

- **Medium concept**: description of the fluid or solid medium being modelled. The medium is defined within the problem and is accessible to the equations. The properties that characterise the medium include:
  - density,
  - diffusivity,
  - thermal conductivity,
  - heat capacity,
  - variation of density with temperature (thermal expansion).

In our case, the problem to be solved contains a single convection-diffusion equation composed of:

- an unsteady term \f$\frac{\partial c}{\partial t}\f$,
- a convection operator \f$\nabla\cdot(\overrightarrow{V}c)\f$,
- a diffusion operator \f$\nabla\cdot(D\nabla c)\f$,
- a source or sink term \f$f=0\f$.

### Input Data File — Requirements Analysis

When developing applications based on the **TRUST** software, it is recommended to start the specification by writing the input data file.
This recommendation stems from the fact that the keywords in the **TRUST** input data file correspond to C++ objects (classes). Starting from the data file thus serves, on one hand, to adopt a user perspective and, on the other hand, to identify the objects (problem, equation(s), ...) that need to be specified and developed.

The parts of the input data file specific to our problem are described below.

**Problem declaration**

This new problem will be named `Probleme_Convection_Diffusion`.


**Convection-diffusion equation declaration**

This new equation will be named `Convection_Diffusion`.
It is composed of two operators (convection and diffusion), initial conditions, and boundary conditions.


**User inputs**

For the convection-diffusion equation, the convection velocity \f$\overrightarrow{V}(\overrightarrow{X},t)\f$ and the diffusion coefficient \f$D(\overrightarrow{X},t)\f$ are provided by the user.
In **TRUST** development practice, these parameters must be associated with an object (equation, problem, ...) in the data file. They will be named `coefficient_diffusion` and `vitesse_convection`.


**Post-processing**

It is useful to post-process the quantity of interest \f$c\f$, the convection velocity \f$\overrightarrow{V}(\overrightarrow{X},t)\f$, and the diffusion coefficient \f$D(\overrightarrow{X},t)\f$.
The post-processing block for these variables takes the following form:


### TRUST Built-in Features

It is strongly recommended to make maximum use of the existing features in the **TRUST** software and to avoid duplicating code.
This minimises development time on one hand, and allows the application to benefit from the test coverage of the **TRUST** software on the other.

Features can be identified through the Doxygen documentation of the **TRUST** software.
We carry out this identification exercise below with respect to our requirements described in the **Input Data File — Requirements Analysis** section.

**Problem**

By examining the Doxygen graph of the `Probleme_base` class (figure below), no existing problem appears close enough to our needs.

![`Probleme_base` class](figures/classProbleme_base.png)

**Convection-diffusion equation**

By examining the Doxygen graph of the `Equation_base` class, we find that a convection-diffusion equation named `Convection_Diffusion_std` exists.
It is a base class for the transport equation of a scalar in laminar flow and holds a reference to the convecting velocity field.
This equation class is pure abstract and therefore cannot be instantiated. Its implementation is incomplete and serves as a base for other derived classes.

Based on the Doxygen graph of the `Convection_Diffusion_Concentration` class (figure below), the child class `Convection_Diffusion_Concentration`, derived from `Convection_Diffusion_std`, appears very close to the convection-diffusion equation. This class is a specific case of `Convection_Diffusion_std` for the transport of one or more constituents.

The specificities of the `Convection_Diffusion_Concentration` equation are:

- diffusion coefficient: this coefficient is a property of a constituent associated with the equation. It is provided by the user in the data file as follows:


- convection velocity: this velocity is obtained via the `vitesse_pour_transport` method, which retrieves it from the equation of the problem that owns it. In standard **TRUST** usage, this convection velocity comes from the Navier-Stokes hydraulics equation.

### Application-Specific Features

The analysis of the **TRUST** software features in the previous section identified a `Convection_Diffusion_Concentration` equation that is very close to the convection-diffusion equation, with a user-supplied diffusion coefficient provided through a constituent. We therefore choose to use this equation as a base to develop a new equation named `Convection_Diffusion`, and associate with it a new constituent `Constituant_Avec_Vitesse` derived from the existing one, extended with a convection velocity property.

The new constituent class `Constituant_Avec_Vitesse` will be declared in the data file as follows:


---

## Developments

### Development Environment: BALTIK Application

A `BALTIK` application development environment is available for applications built on top of either the **TRUST** software, its numerical kernel ("Kernel"), or another `BALTIK` application. This environment handles the compilation of such applications (application-specific sources) as well as modified and/or added **TRUST** sources.
In this framework, an application consists of:

- the numerical "Kernel", the **TRUST** software, or another `BALTIK` application,
- the application-specific sources.

In practice, the application-specific sources are gathered in a working directory called `src`, which is compiled following the standard **TRUST** compilation process. Sources from the "Kernel" or from **TRUST** that need to be modified from their standard versions can also be included.

#### Creating and Configuring the BALTIK Application

Creating a `BALTIK` application consists in creating a directory named after the application (e.g. `convection_diffusion`).

```bash
mkdir convection_diffusion
```

Setting up a `BALTIK` application involves creating a configuration file `project.cfg` inside the `convection_diffusion` directory.


@note In this configuration file, only the application name is mandatory. A number of optional parameters also exist, such as the author name and the name of the executable to generate.

To configure our `BALTIK` application, the following steps must be executed:

- initialise the source directory:

  ```bash
  cd convection_diffusion
  mkdir src
  ```

  `src` is the directory containing either modified **TRUST** or Kernel TRUST sources, or additional sources specific to our application.

- configure the `BALTIK` application:

  ```bash
  source $TRUST_ROOT/env_TRUST.sh
  cd convection_diffusion
  baltik_build_configure
  ```

  The `$TRUST_ROOT` environment variable is the path where the **TRUST** software is installed. The first instruction initialises the TRUST environment. The second generates a `configure` file and must be run from the `convection_diffusion` directory.

#### Compiling the BALTIK Application

The `Makefile` is generated by running the following instruction from the `convection_diffusion` directory:

```bash
cd convection_diffusion
./configure
```

@note When new sources are added to the `src` directory, this instruction must be re-run to update the `Makefile`.

Once the `Makefile` has been generated, two compilation modes are available:

- debug mode:

  ```bash
  cd convection_diffusion
  make debug
  ```

- optimised mode:

  ```bash
  cd convection_diffusion
  make optim
  ```

### Constituent Class `Constituant_Avec_Vitesse`

As specified in the **Specification** section, the user parameters for convection velocity and diffusion coefficient are assigned to the new constituent class `Constituant_Avec_Vitesse`, which inherits from the `Constituant` class of the **TRUST** software. The distinguishing feature of this new class is that it holds a convection velocity field.

Handling this additional attribute requires the following actions:

1. adding attributes:
   - `C`: will hold the content of `vitesse_convection` from the data file.
   - `vitesse_transport`: will hold the discretised convection velocity.

2. overriding the following methods (only the features specific to `Constituant_Avec_Vitesse` are described here):
   - `initialiser`: evaluate the convection velocity in space at the initial time.
   - `set_param`: assign the content of `vitesse_convection` from the data file to the attribute `C`.
   - `mettre_a_jour`: evaluate the convection velocity in space at the current time.
   - `discretiser`: discretise the convection velocity from the user-provided expression.

3. developing new methods:
   - `vit_convection_constituant`: returns the attribute holding the non-discretised convection velocity.
   - `vitesse_pour_transport`: returns the convection velocity discretised on the mesh faces.

4. developing the `printOn` and `readOn` methods: the `Constituant_Avec_Vitesse` class can be used in the data file and is therefore instantiable.
   To declare this type of class, the `Declare_instanciable` macro is used. This macro declares the signatures of the `printOn` and `readOn` methods, which are needed for reading objects from input streams (data file) and writing objects to output streams.

#### Header File

#### Source File

In the source file, the `Implemente_instanciable` macro serves, on one hand, to define the keyword used in the data file to represent the `Constituant_Avec_Vitesse` class, and on the other hand, to declare its parent class.
In our case, the class name `Constituant_Avec_Vitesse` is chosen as the keyword in the data file.


### Equation Class `Convection_Diffusion`

The distinguishing feature of this new equation class `Convection_Diffusion` compared to the base class `Convection_Diffusion_Concentration` is that `Convection_Diffusion` holds a reference to the new constituent class `Constituant_Avec_Vitesse`, and the convection velocity is provided by that class rather than by the resolution of a Navier-Stokes hydraulics equation.

Handling this new reference requires the following actions:

1. adding the reference `le_constituant`. To do so, the file `TRUST_Ref.h` must be included and our attribute added as a reference `REF(Constituant_Avec_Vitesse)`. This non-standard C++ technique is used in the **TRUST** software to prevent developers from manipulating pointers, which can be a source of memory leaks.

2. overriding the following methods (only the features specific to `Convection_Diffusion` are described here):
   - `associer_milieu_base`: associates the medium with the convection-diffusion equation, and more precisely associates the new constituent `Constituant_Avec_Vitesse` with the equation.
   - `associer_constituant`: associates the constituent with the `Convection_Diffusion` equation. This method is called by `associer_milieu_base`.
   - `lire_motcle_non_standard`: reads non-standard keywords from the data file, such as `convection`, in order to associate the convection velocity with the equation.

3. developing the `printOn` and `readOn` methods.

#### Header File

#### Source File

### Problem Class `Probleme_Convection_Diffusion`

In the **Specification** section, no existing problem was identified as a suitable base for the new `Probleme_Convection_Diffusion` problem. We therefore take the pure abstract class `Probleme_base` as the parent class.

To use this new problem class, the following steps are required:

1. adding the attribute `eq_conv_diff`, which represents the convection-diffusion equation.

2. completing the implementation of `Probleme_base` by developing the following methods:
   - `nombre_d_equations`: returns the number of equations that make up this problem (one in our case).
   - `equation`: returns the attribute holding the convection-diffusion equation (`eq_conv_diff`).

3. overriding the method `associer_milieu_base` to associate the medium with `eq_conv_diff`.

4. developing the `printOn` and `readOn` methods.

#### Header File

#### Source File

---

## Development Verification Process

In this section, we present:

- the setup of verification tests,
- the automatic generation of a verification report via the **`genererCourbes`** tool,
- the non-regression verification.

### Setting Up Verification Tests

Setting up the verification tests involves creating:

- a directory for each verification test:

  ```bash
  cd convection_diffusion
  mkdir -p share/Validation/Rapports_automatiques/Cas-tests/src/cas-test1
  mkdir -p share/Validation/Rapports_automatiques/Cas-tests/src/cas-test2
  ```

  `Cas-tests` is the directory containing the verification tests. `cas-test1` and `cas-test2` contain the two verification cases.

- a `PRM` file in the directory `convection_diffusion/share/Validation/Rapports_automatiques/Cas-tests/src`, used to automatically generate a verification report (cf. section **Verification Report**).

- an input data file for each verification case.

### Generating the Verification Report

Generation consists in running the command `make validation` from the `convection_diffusion` directory. In that same directory, this command will produce the verification report `validation.pdf`. This report is presented in the **Verification Report** section below.

### Non-Regression Verification

To use the previously set up tests as non-regression tests, it suffices to:

- create a directory to hold the reference results inside the `convection_diffusion` directory:

  ```bash
  cd convection_diffusion
  mkdir -p tests/Reference/Validation
  ```

- run the tests and check for non-regression by executing the command `make check_optim` from the `convection_diffusion` directory. The figure below shows that the sequential and parallel results exhibit no regression.

![Non-regression test results](figures/tnr.png)

---

## Verification Report

@note The verification report body is not included in this page. Run `make validation` from the `convection_diffusion` directory to generate `validation.pdf`.

---

## Conclusion

In this tutorial, we have described, through the resolution of a convection-diffusion equation, the following aspects:

- the development environment of the **TRUST** software, by building a new `BALTIK` application,

- the specification and development workflow. This workflow requires technical choices, such as associating the convection velocity with a constituent rather than with the equation or problem.
  These choices are neither unique nor necessarily better than alternatives. The main purpose of this tutorial is to present a workflow accompanied by its own technical choices, while providing some recommendations and methodologies to follow during the specification phase.

- the setup of non-regression tests and the automatic generation of a development verification report using tools provided by the **TRUST** software.
