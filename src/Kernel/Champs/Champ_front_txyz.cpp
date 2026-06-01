/****************************************************************************
* Copyright (c) 2025, CEA
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
//   New boundary field

#include <Champ_front_txyz.h>
#include <Domaine.h>
#include <Frontiere_dis_base.h>
#include <Domaine_VF.h>
#include <TRUSTTrav.h>

Implemente_instanciable(Champ_front_txyz,"Champ_front_fonc_txyz",Ch_front_var_instationnaire_indep);
// XD champ_front_fonc_txyz front_field_base champ_front_fonc_txyz NO_BRACE Boundary field which is not constant in
// XD_CONT space and in time.
// XD attr val listchaine val REQ Values of field components (mathematical expressions).

/*! @brief Print the field to an output stream.
 *
 * Prints the field size and the (constant) value on
 *     the boundary.
 *
 * @param (Sortie& os) an output stream
 * @return (Sortie&) the modified output stream
 */
Sortie& Champ_front_txyz::printOn(Sortie& os) const
{

  const DoubleTab& tab=valeurs();
  os << tab.size() << " ";
  for(int i=0; i<tab.size(); i++)
    os << tab(0,i);
  return os;
}


/*! @brief Read the field from an input stream.
 *
 * Format:
 *       Champ_front_txyz nb_compo vrel_1 ... [vrel_i]
 *
 * @param (Entree& is) an input stream
 * @return (Entree& is) the modified input stream
 */
Entree& Champ_front_txyz::readOn(Entree& is)
{
  //Cout << "Dans  Champ_front_txyz.readOn !!!!!!!!!" << finl;
  int dim;
  dim=lire_dimension(is,que_suis_je());
  fixer_nb_comp(dim);

  fxyz.dimensionner(dim);


  //        Cout << "dim = " << dim << finl;
  for (int i = 0; i<dim; i++)
    {
      Nom tmp;
      //Cout << "i = " << i << finl;
      is >> tmp;
      //Cout << "fonc = " << tmp << finl;
      Cerr << "Reading and interpretation of the function " << tmp << finl;
      fxyz[i].setNbVar(4);
      fxyz[i].setString(tmp);
      fxyz[i].addVar("t");
      fxyz[i].addVar("x");
      fxyz[i].addVar("y");
      fxyz[i].addVar("z");
      fxyz[i].parseString();
      Cerr << "Interpretation of function " << tmp << " Ok" << finl;
      //        Cout << "end = " << tmp << finl;
    }

  return is;
}



/*! @brief Returns the object upcast to Champ_front_base&
 *
 * @param (Champ_front_base& ch)
 * @return (Champ_front_base&) (*this) upcast to Champ_front_base&
 */
Champ_front_base& Champ_front_txyz::affecter_(const Champ_front_base& ch)
{
  return *this;
}


void Champ_front_txyz::mettre_a_jour(double temps)
{
  int dim=nb_comp();
  const Frontiere_dis_base& fr_dis=frontiere_dis();
  const Domaine_VF& zvf = ref_cast(Domaine_VF, fr_dis.domaine_dis());
  int nb_faces=ref_cast(Front_VF, fr_dis).nb_faces();
  int premiere_face = ref_cast(Front_VF, fr_dis).num_premiere_face();
  DoubleTab& tab=valeurs_au_temps(temps);
  for(int i=0; i<nb_faces; i++)
    {
      for(int k=0; k<dim; k++)
        {
          fxyz[k].setVar(0,temps);
          fxyz[k].setVar(1,zvf.xv(premiere_face + i, 0));
          fxyz[k].setVar(2,zvf.xv(premiere_face + i, 1));
          if (dimension >= 3)
            fxyz[k].setVar(3,zvf.xv(premiere_face + i, 2));
          tab(i,k)=fxyz[k].eval();
        }
    }
  tab.echange_espace_virtuel();
}

double Champ_front_txyz::valeur_au_temps_et_au_point(double temps,int som,double x,double y, double z, int k) const
{
  Parser_U& fxyzk=fxyz[k];
  fxyzk.setVar(0,temps);
  fxyzk.setVar(1,x);
  fxyzk.setVar(2,y);
  fxyzk.setVar(3,z);
  return fxyzk.eval();
}
