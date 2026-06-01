/****************************************************************************
* Copyright (c) 2026, CEA
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


#ifndef Champ_Generique_Transformation_included
#define Champ_Generique_Transformation_included

#include <Champ_Gen_de_Champs_Gen.h>
#include <TRUST_Vector.h>
#include <Parser_U.h>

/*! @brief class Champ_Generique_Transformation
 *
 *  Field intended to post-process a "transformation" depending on generic fields and (or) on x,y,z and t
 *  The class carries:
 *    -vector of Parser (fxyz)
 *    -the expression(s) to express the transformation (les_fct)
 *    -the chosen transformation method (methode_)
 *
 */

//// Data file syntax to follow
//
// "field_name" Transformation {
//                 "transform_type" "type_info"
//                sources { ...{ ...source ref_Champ { Pb_champ "pb_name" "discrete_field_name1" } } ,
//                          ...{ ...source ref_Champ { Pb_champ "pb_name" "discrete_field_name2" } } ,
//                            ...
//                          }
//                  [localization "loc"]
//             }
//
// "field_name" set by the user will be the name of the generic field.
// "transform_type" can be:
//          "function"           allows to apply a formula using the values of the specified source fields and can depend on x,y,z and t.
//                              -> type_info: f(field_name1,field_name2...,x,y,z,t) analytical expression
//
//          "dot_product"        allows to calculate the dot product of two vectors.
//                              -> type_info: empty since the transformation expression is constructed automatically.
//                                 2 sources of vector nature must be specified.
//
//          "norm"               allows to calculate the norm of a vector.
//                              -> type_info: same as dot_product
//                                 1 source of vector nature must be specified.
//
//          "vector"             allows to construct a field of vector nature from the specified components.
//                              -> type_info: pb_name nb_compo f1(field_name1,field_name2...,x,y,z,t) ...fn(field_name1,field_name2...,x,y,z,t)
//                                 pb_name   : name of the problem.
//                                 nb_compo  : number of components of the vector field to create.
//                                 f1,...,fn : expression of the components of the vector field to create.
//                                If sources are specified, they must be of scalar nature.
//                                 Note: if no source is specified by the user, the current version assigns a
//                                     default source based on the unknown field of the first equation of the problem.
//                                     Eventually (after design review) this workaround will be eliminated.
//
//          "component"          allows to construct a scalar field by extracting a component of a vector field.
//                              -> type_info: num_component: number of the component to extract.
//                                1 source of vector nature must be specified.
//
// "generic_field_type" type of a generic field
// "discrete_field_name1" designates the name of the field constituting the first source
// "discrete_field_name2" designates the name of the field constituting the second source
// ...
// "loc" allows to specify a particular localization to evaluate the values of the storage space
// "unit" allows the user to give a unit to the obtained field
//  possible values: "elem", "som", "faces", "elem_dg" and "elem_som".
//  In the case where no localization is specified, the retained localization is that of the first source support.
//

class Champ_Generique_Transformation : public Champ_Gen_de_Champs_Gen
{

  Declare_instanciable(Champ_Generique_Transformation);

public:

  void set_param(Param& param) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void verifier_coherence_donnees();
  void verifier_localisation();
  const Noms get_property(const Motcle& query) const override;
  const Champ_base&  get_champ(OWN_PTR(Champ_base)& espace_stockage) const override;
  const Champ_base&   get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const override;

  Entity  get_localisation(const int index = -1) const override;
  const Motcle get_directive_pour_discr() const override;
  void completer(const Postraitement_base& post) override;
  void nommer_source() override;
  int preparer_macro();
  void creer_expression_macro();


protected:

  Nom methode_;                //Method indicating the type of transformation selected
  Noms les_fct;                //Contains the expression of the combination
  mutable VECT(Parser_U) fxyz; //Parser used to evaluate the value taken by the combination
  int nb_comp_ = 1;                //Number of components of the evaluated field
  Motcle localisation_;        //Localisation of the support for evaluating the expression
  Nom unite_;                 //unit of the obtained field (to be specified by the user)
  Nature_du_champ nature_ch = scalaire;   //Nature of the evaluated field
  bool fictive_source_ = false;

private:
  mutable OWN_PTR(Champ_Fonc_base) espace_stockage_;
};

#endif
