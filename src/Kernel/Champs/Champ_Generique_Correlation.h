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

#ifndef Champ_Generique_Correlation_included
#define Champ_Generique_Correlation_included

#include <Champ_Generique_Statistiques_base.h>
#include <TRUSTTabs_forward.h>
#include <Op_Correlation.h>

class Postraitement_base;

/*! @brief class Champ_Generique_Correlation
 *
 *  Field intended to post-process a correlation of a generic field
 *  The class carries a statistical operator "Op_Correlation"
 *
 */

//// Data file syntax to follow
//
// "field_name" correlation { t_deb "val_tdeb" t_fin "val_tfin"
//                source "generic_field_type" { ...source ref_Champ { Pb_champ "pb_name" "discrete_field_name1" } }
//                source "generic_field_type" { ...source ref_Champ { Pb_champ "pb_name" "discrete_field_name2" } }
//               }
// "field_name" set by the user will be the name of the generic field
// "val_tdeb" and "val_tfin" value of the start and end of statistics for this field
// "generic_field_type" type of a generic field
// Two sources are to be specified for this type of field

class Champ_Generique_Correlation : public Champ_Generique_Statistiques_base
{

  Declare_instanciable(Champ_Generique_Correlation);

public:

  const Noms get_property(const Motcle& query) const override;

  inline double temps() const override
  {
    return Op_Correlation_.integrale().le_champ_calcule().temps();
  };
  inline const Integrale_tps_Champ& integrale() const override
  {
    return Op_Correlation_.integrale();
  };

  inline const Operateur_Statistique_tps_base& Operateur_Statistique() const override;
  inline Operateur_Statistique_tps_base& Operateur_Statistique() override;
  void completer(const Postraitement_base& post) override;

  const Motcle get_directive_pour_discr() const override;
  const Champ_base&  get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const override;
  const Champ_base& get_champ(OWN_PTR(Champ_base)& espace_stockage) const override;
  void nommer_source() override;
  int get_info_type_post() const override;

protected:

  Op_Correlation Op_Correlation_;

private:
  mutable OWN_PTR(Champ_Fonc_base) espace_stockage_;
};

inline const Operateur_Statistique_tps_base& Champ_Generique_Correlation::Operateur_Statistique() const
{
  return Op_Correlation_;
}

inline Operateur_Statistique_tps_base& Champ_Generique_Correlation::Operateur_Statistique()
{
  return Op_Correlation_;
}

#endif
