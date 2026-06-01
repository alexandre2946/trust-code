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

#include <Flux_radiatif_base.h>
#include <Front_VF.h>
#include <Param.h>

Implemente_base(Flux_radiatif_base,"Flux_radiatif_base",Neumann_paroi);

Sortie& Flux_radiatif_base::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Flux_radiatif_base::readOn(Entree& is)
{
  Motcle motlu;
  Motcles les_motcles(2);
  {
    les_motcles[0] = "emissivite";
    les_motcles[1] = "A";
  }

  int ind = 0;
  while (ind < 2)
    {
      is >> motlu;
      int rang = les_motcles.search(motlu);

      switch(rang)
        {
        case 0:
          {
            is >> emissivite_;
            break;
          }
        case 1:
          {
            is >> A_;
            break;
          }
        default:
          {
            Cerr << "Error reading the boundary condition of type " << finl;
            Cerr << "Flux_radiatif_base " << finl;
            Cerr << "Expected " << les_motcles << "instead of " << motlu << finl;
            Process::exit();
          }
        }
      ind++;
    }

  le_champ_front.typer("Champ_front_fonc");
  champ_front().fixer_nb_comp(1);

  flux_radiatif_.typer("Champ_front_fonc");
  flux_radiatif().fixer_nb_comp(1);

  return is;
}

void Flux_radiatif_base::completer()
{
  Neumann_paroi::completer();
  // Type the champ_front flux_radiatif_ associated with the boundary condition
  const Front_VF& front_vf = ref_cast(Front_VF, le_champ_front->frontiere_dis());
  int nb_comp = 1;

  flux_radiatif().nommer(front_vf.le_nom());
  DoubleTab& tab_flux = flux_radiatif().valeurs();
  tab_flux.resize(front_vf.nb_faces(), nb_comp);

  champ_front().nommer(front_vf.le_nom());
  DoubleTab& tab = champ_front().valeurs();
  tab.resize(front_vf.nb_faces(), nb_comp);
  emissivite_->associer_fr_dis_base(front_vf);
}

/*! @brief Returns the imposed flux value at the radiative wall.
 *
 */
double Flux_radiatif_base::flux_impose(int i) const
{
  if (le_champ_front->valeurs().size() == 1)
    return le_champ_front->valeurs()(0, 0);
  else if (le_champ_front->valeurs().dimension(1) == 1)
    return le_champ_front->valeurs()(i, 0);
  else
    Cerr << "Flux_radiatif_base::flux_impose error" << finl;
  Process::exit();
  return 0.;
}

/*! @brief Returns the imposed flux value at the radiative wall.
 *
 */
double Flux_radiatif_base::flux_impose(int i, int j) const
{
  if (le_champ_front->valeurs().dimension(0) == 1)
    return le_champ_front->valeurs()(0, j);
  else
    return le_champ_front->valeurs()(i, j);
}
