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

#include <Champ_front_debit_QC.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Fluide_Quasi_Compressible.h>
#include <Motcle.h>
#include <Domaine_VDF.h>
#include <Interprete.h>

Implemente_instanciable(Champ_front_debit_QC,"Champ_front_debit_QC_VDF",Ch_front_var_instationnaire_indep);
// XD Champ_front_debit_QC_VDF front_field_base Champ_front_debit_QC_VDF NO_BRACE This keyword is used to define a flow
// XD_CONT rate field for quasi-compressible fluids in VDF discretization. The flow rate is kept constant during a
// XD_CONT transient.
// XD attr dimension entier dim REQ Problem dimension
// XD attr liste bloc_lecture liste REQ List of the mass flow rate values [kg/s/m2] with the following syntaxe: { val1
// XD_CONT ... valdim }
// XD attr moyen chaine moyen OPT Option to use rho mean value
// XD attr pb_name chaine pb_name REQ Problem name



/*! @brief Prints to an output stream in the format: size
 *
 *     value(0) ... value(i)  ... value(size-1)
 *
 * @param os output stream
 * @return modified output stream
 */
Sortie& Champ_front_debit_QC::printOn(Sortie& os) const
{
  const DoubleTab& tab=valeurs();
  os << tab.size() << " ";
  for(int i=0; i<tab.size(); i++)
    os << tab(0,i);
  return os;
}

/*! @brief Reads from an input stream in the format: number_of_components
 *
 *     mean mean(0) ... mean(number_of_components-1)
 *     mean amplitude(0) ... amplitude(number_of_components-1)
 *
 * @param is input stream
 * @return modified input stream
 * @throws opening brace expected
 * @throws unknown keyword at this location
 * @throws closing brace expected
 */
Entree& Champ_front_debit_QC::readOn(Entree& is)
{
  Cerr<<"Champ_front_debit_QC_VDF usage : dim { val1 .. valdim } [ moyen ] nom_pb"<<finl;
  int dim;
  is >> dim;
  Motcle motlu;
  Motcles les_mots(2);
  les_mots[0]="{";
  les_mots[1]="}";
  is >> motlu;
  if (motlu != les_mots[0])
    {
      Cerr << "Error reading a Champ_front_debit_QC" << finl;
      Cerr << "Expected { instead of " << motlu << finl;
      exit();
    }
  fixer_nb_comp(dim);
  Debit.resize(dim);
  for(int i=0; i<dim; i++)
    is >> Debit(i);
  is >> motlu;
  if (motlu != les_mots[1])
    {
      Cerr << "Error reading a Champ_front_debit_QC" << finl;
      Cerr << "Expected } instead of " << finl;
      exit();
    }
  Nom nom_pb;
  is>>nom_pb;
  Motcle nom;
  nom=nom_pb;
  if (nom=="moyen")
    {
      ismoyen=1;
      is >> nom_pb;
    }
  Objet_U& ob1=Interprete::objet(nom_pb);
  const Probleme_base& pb=ref_cast(Probleme_base,ob1);
  fluide=ref_cast(Fluide_Quasi_Compressible,pb.equation(0).milieu());
  return is;
}


/*! @brief Not implemented!!
 *
 * @param ch source boundary field
 * @return reference to this boundary field
 */
Champ_front_base& Champ_front_debit_QC::affecter_(const Champ_front_base& ch)
{
  return *this;
}

/*! @brief Updates the time and re-draws random noise values.
 *
 * @param tps current time for the update
 */
void Champ_front_debit_QC::mettre_a_jour(double tps)
{


  const Frontiere& front=la_frontiere_dis->frontiere();
  int nb_faces=front.nb_faces();
  DoubleTab& tab=valeurs();
  int dim=nb_comp();
  const Domaine_VDF& le_dom_VDF = ref_cast(Domaine_VDF,domaine_dis());
  const IntTab& face_voisins=le_dom_VDF.face_voisins();
  const Front_VF& front_vf=ref_cast(Front_VF,la_frontiere_dis.valeur());
  int ndeb = front_vf.num_premiere_face();
  int nfin = ndeb + nb_faces;
  const DoubleTab& tab_rhonp1P0 =fluide->loi_etat()->rho_np1();
  if (ismoyen==0)

    for (int num_face=ndeb; num_face<nfin; num_face++)
      {
        int n0 = face_voisins(num_face, 0);

        if (n0 == -1)
          n0 = face_voisins(num_face, 1);
        for (int ori=0; ori<dim; ori++)
          tab(num_face-ndeb,ori)=Debit(ori)/tab_rhonp1P0(n0);
      }
  else
    {
      int num_face;
      double rho_moy=0,S=0,s;
      const DoubleVect& surface=le_dom_VDF.face_surfaces();
      for ( num_face=ndeb; num_face<nfin; num_face++)
        {
          int n0 = face_voisins(num_face, 0);

          if (n0 == -1)
            n0 = face_voisins(num_face, 1);
          s=surface(num_face);
          S+=s;
          rho_moy+=s*tab_rhonp1P0(n0);
        }
      // Optimization: combine 2 mp_sum into 1 collective call
      mp_sum_for_each(S, rho_moy);
      rho_moy/=S;
      for ( num_face=ndeb; num_face<nfin; num_face++)
        for (int ori=0; ori<dim; ori++)
          tab(num_face-ndeb,ori)=Debit(ori)/rho_moy;
    }

}
