/****************************************************************************
* Copyright (c) 2024, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/

#ifndef Champs_compris_interface_included
#define Champs_compris_interface_included

#include <TRUST_Ref.h>

class Champ_base;
class Motcle;
class Noms;
class Nom;

enum Option { NONE, DESCRIPTION };

/*! @brief class Champs_compris_interface This class contains an interface of methods intended to manage
 *
 *               the understood (post-processable) fields for a given class.
 *               The classes that inherit from this class are: Probleme_base, Equation_base, Milieu_base,
 *               Operateur_base, Source_base, Traitement_particulier_NS_base, Traitement_particulier_Solide_base,
 *               Modele_turbulence_hyd_base, Modele_turbulence_scal_base, Loi_Etat_base, Modele_Fonc_Bas_Reynolds_Base,
 *               Modele_Fonc_Bas_Reynolds_Thermique_Base and the interface is propagated in their derived classes
 *
 *         Interface methods:
 *         creer_champ()
 *         get_champ()
 *         get_noms_champs_postraitables()
 *
 */
class Champs_compris_interface
{
public :
  virtual inline ~Champs_compris_interface() {};
  virtual void creer_champ(const Motcle& motlu) =0;
  virtual const Champ_base& get_champ(const Motcle& nom) const=0;
  virtual void get_noms_champs_postraitables(Noms& nom, Option opt=NONE) const=0;

  virtual bool has_champ(const Motcle& nom, OBS_PTR(Champ_base)& ref_champ) const=0;
  virtual bool has_champ(const Motcle& nom) const=0;

  //To illustrate the description of the methods below, one can
  //refer to their implementation in Probleme_base, Equation_base and Navier_Stokes_std

  // the creer_champ() method
  /////////////////////////////////////////////////////////////////////////////////
  //This method allows creating a "calculated field". This field of type Champ_Fonc_base
  //(e.g.: vorticity) is estimated only for post-processing and is not used for
  //solving the problem.
  //The problem delegates to its medium and its equations the creation of the field.
  //A class that identifies the field as one of its attributes creates it
  //by launching its discretization. A class that carries a "calculated field" will
  //therefore have a creer_champ() method. This method will contain the call to the
  //method of a parent class that also has one or more "calculated fields".
  //The field update will only be done if a request via the get_champ() method
  //(see below) is made.
  //Currently only equations carry calculated fields because they
  //have a discretization object.
  //////////////////////////////////////////////////////////////////////////////

  // the get_champ() method
  ///////////////////////////////////////////////////////////////////////
  //This method allows launching a request to retrieve a REF to a
  //field from its identifier (name) or that of one of its components.
  //The general mechanism is as follows:
  //-The problem delegates the request to its medium then to its equations if the medium has not identified
  //the name.
  //-The class that makes a request (e.g.: Navier_Stokes) queries its parent class to check if
  //the identifier corresponds to a field understood by it or to one of its attributes.
  //The mechanism is recursive, so the parent class queries its own parent class...
  //If the upstream hierarchy of the class that launched the request (e.g.: Navier_Stokes) does not recognize
  //the identifier, it tests its own understood fields then delegates the request to its attributes
  //in case of a negative response.
  ////////////////////////////////////////////////////////////////////////

  // the get_noms_champs_postraitables() method
  ////////////////////////////////////////////////////////////////////////
  //This method allows informing the user about the post-processable fields for a given problem.
  //The class scanning principle is the same as for the get_champ() method (see above).
  //-If the DESCRIPTION option is activated:
  //Each class concerned by the problem writes the identifier of its fields (or their components)
  //in the output file (err).
  //-If the NONE option is activated (default value):
  //the list of names is added to nom.
  ////////////////////////////////////////////////////////////////////////

};

#endif
