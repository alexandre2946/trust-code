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

#include <Domaine_Cl_dis_base.h>
#include <Discretisation_base.h>
#include <Dirichlet_homogene.h>
#include <Format_Post_base.h>
#include <communications.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Schema_Comm.h>
#include <Domaine_VF.h>
#include <Champ_base.h>
#include <TRUSTVects.h>
#include <TRUSTTrav.h>
#include <Dirichlet.h>
#include <TRUSTList.h>
#include <Symetrie.h>
#include <strings.h>

Implemente_base_sans_constructeur(Champ_base,"Champ_base",Field_base);

Sortie& Champ_base::printOn(Sortie& os) const
{
  return os << le_nom() << finl;
}

Entree& Champ_base::readOn(Entree& is)
{
  return is >> nom_;
}

/*! @brief Default constructor of a Champ_base.
 *
 * Sets the field at time 0, specifies an empty unit,
 *     gives the name "anonymous" to the field and gives it a
 *     scalar nature.
 *
 */
Champ_base::Champ_base()
{
  changer_temps(0.);
  fixer_unite(".");
  nommer("anonyme");
  fixer_nature_du_champ(scalaire);
  nb_compo_ = 1; // By default, scalar field
}

/*! @brief Computes the "values" of the field at the point with coordinates "pos".
 *
 * In this base class, the implementation calls
 *   valeur_aux(const DoubleTab &, DoubleTab &)
 *
 * @param (DoubleVect&) the coordinates of the point where to evaluate the field
 * @param (DoubleVect& valeurs) On input: must have the right size (nb_comp), on output contains the components of the field.
 * @return (reference to "valeurs")
 */
DoubleVect& Champ_base::valeur_a(const DoubleVect& pos, DoubleVect& les_valeurs) const
{
  //assert(les_valeurs.size() == nb_comp());
  DoubleTrav values(1,les_valeurs.size());
  int taille=pos.size();
  DoubleTrav pos2(1,taille);
  for (int dir=0; dir<taille; dir++) pos2(0,dir)=pos(dir);
  valeur_aux(pos2,values);
  for (int comp=0; comp<les_valeurs.size(); comp++)
    les_valeurs(comp)=values(0,comp);
  return les_valeurs;
}

/*! @brief Causes an error! Must be overridden by derived classes
 *
 *     Not pure virtual for convenience of development!
 *     Returns the value of the field at the point specified
 *     by its coordinates, indicating that this point is
 *     located in a specified element.
 *
 * @param (DoubleVect&) the coordinates of the calculation point
 * @param (DoubleVect& les_valeurs) the value of the field at the specified point
 * @param (int) the element in which the calculation point is located
 * @return (DoubleVect&) the value of the field at the specified point
 */
DoubleVect& Champ_base::valeur_a_elem(const DoubleVect& ,
                                      DoubleVect& les_valeurs,
                                      int ) const
{
  Cerr << que_suis_je();
  Cerr << "::valeur_a_elem(...) is not coded " << finl ;
  exit();
  return les_valeurs;
}

/*! @brief Computes the point value of the component "compo" of the field at the point with coordinates pos.
 *
 *     In the base class, the implementation calls
 *     valeur_a(const DoubleVect &, DoubleVect &)
 *
 * @param (DoubleVect&) the coordinates of the calculation point
 * @param (int) the index of the component of the field to calculate
 * @return (double)
 */
double Champ_base::valeur_a_compo(const DoubleVect& pos, int compo) const
{
  DoubleVect values(nb_comp());
  valeur_a(pos, values);
  const double x = values[compo];
  return x;
}

/*! @brief Causes an error! Must be overridden by derived classes
 *
 *  Not pure virtual for convenience of development!
 *
 */
double Champ_base::valeur_a_elem_compo(const DoubleVect&, int ,int ) const
{
  Cerr << que_suis_je();
  Cerr << "::valeur_a_elem_compo(...) is not coded " << finl ;
  exit();
  return 0;
}

/*! @brief This method, generic but slow (calculation of gravity centers, filling les_poly, use of shape functions in the discretized field)
 * can be overridden by the discretized field for a much faster implementation
 */
DoubleTab& Champ_base::valeur_aux_centres_de_gravite(const Domaine& dom, DoubleTab& les_valeurs) const
{
#ifdef TRUST_USE_GPU
  Cerr << "Warning, try to implement a " << que_suis_je() << "::valeur_aux_centres_de_gravite() for a faster compute." << finl;
#endif
  int nb_elem = les_valeurs.dimension(0);
  DoubleTrav positions(nb_elem, dimension);
  if(sub_type(Champ_Inc_base, *this))
    {
      const Domaine_VF& zvf = ref_cast(Domaine_VF,ref_cast(Champ_Inc_base, *this).domaine_dis_base());
      // PL: ToDo Kokkos kernel host kept because difficult-to-find bug (case decroissance_ktau_jdd1 with TrioCFD):
      // stencil.append_line() allocates memory via a resize() on HOST memory already allocated on DEVICE!
      // Probably, DEVICE memory not properly deallocated in a PolyMAC_CDO mechanism not used in VEF...
      if (zvf.xp().isDataOnDevice())
        {
          // To avoid a resize by nb_elem_tot per call to xp()
          CDoubleTabView xp = zvf.xp().view_ro();
          DoubleTabView positions_v = positions.view_wo();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__),
          Kokkos::MDRangePolicy < Kokkos::Rank < 2 >> ({ 0, 0 },
          {nb_elem, dimension}), KOKKOS_LAMBDA(
            const int i,
            const int j)
          {
            positions_v(i, j) = xp(i, j);
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
      else
        {
          for (int i=0; i<nb_elem; i++)
            for (int k=0; k<dimension ; k++)
              positions(i,k) = zvf.xp(i,k);
        }
    }
  else
    {
      dom.calculer_centres_gravite(positions);
    }

  IntTrav les_polys(nb_elem);
  IntArrView les_polys_v = static_cast<ArrOfInt&>(les_polys).view_wo();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_elem, KOKKOS_LAMBDA(const int i)
  {
    les_polys_v(i) = i;
  });
  end_gpu_timer(__KERNEL_NAME__);

  return valeur_aux_elems(positions, les_polys, les_valeurs);
}

/*! @brief Causes an error! Must be overridden by derived classes
 *
 *     Not pure virtual for convenience of development!
 *     Returns the values of the field at the specified points
 *     by their coordinates.
 *
 * @param (DoubleTab&) the array of coordinates of the calculation points
 * @param (DoubleTab& les_valeurs) the array of field values at the specified points
 * @return (DoubleTab&) the array of field values at the specified points
 */
DoubleTab& Champ_base::valeur_aux(const DoubleTab& ,
                                  DoubleTab& les_valeurs) const
{
  Cerr << que_suis_je();
  Cerr << "::valeur_aux(const DoubleTab& ,DoubleTab& ) is not coded " << finl ;
  exit();
  return les_valeurs;
}

/*! @brief Same as valeur_aux(const DoubleTab &, DoubleTab &), but computes only the component compo of the field.
 *
 *   In the implementation of champ_base, we call
 *   valeur_aux(const DoubleTab &, DoubleTab &)
 *
 * @param (pos) the array of coordinates of the calculation points (we do not process the virtual space of the array)
 * @param (les_valeurs) destination array of values to calculate. The values array must have the right size on input, i.e. les_valeurs.size() == pos.dimension(0)
 * @param (compo) the index of the component of the field to calculate
 * @return (reference to the array les_valeurs)
 */
DoubleVect& Champ_base::valeur_aux_compo(const DoubleTab& pos ,
                                         DoubleVect& les_valeurs,
                                         int compo) const
{
  // Not optimal but works
  int nb_val=pos.dimension(0);
  // The array les_valeurs must have the right size on input:
//  assert(les_valeurs.size() == nb_val);
  DoubleTrav prov(nb_val,nb_comp());
  valeur_aux(pos,prov);
  for (int i=0; i<nb_val; i++)
    les_valeurs(i)=prov(i,compo);
  return les_valeurs;
}

/*! @brief Causes an error! Must be overridden by derived classes
 *
 *     Not pure virtual for convenience of development!
 *     Returns the values of the field at the specified points
 *     by their coordinates, indicating that the calculation
 *     points are located in the specified elements.
 *
 * @param (DoubleTab&) the array of coordinates of the calculation points
 * @param (IntVect&) the array of elements in which the calculation points are located
 * @param (DoubleTab& les_valeurs) the array of field values at the specified points
 * @return (DoubleTab&) the array of field values at the specified points
 */
DoubleTab& Champ_base::valeur_aux_elems(const DoubleTab&,
                                        const IntVect& ,
                                        DoubleTab& les_valeurs) const
{
  Cerr << que_suis_je();
  Cerr << "::valeur_aux_elems(...) is not coded " << finl ;
  exit();
  return les_valeurs;
}

/*! @brief Causes an error! Must be overridden by derived classes
 *
 *     Not pure virtual for convenience of development!
 *     Returns the values of a component of the field at the specified points
 *     by their coordinates, indicating that the calculation
 *     points are located in the specified elements.
 *
 * @param (DoubleTab&) the array of coordinates of the calculation points
 * @param (IntVect&) the array of elements in which the calculation points are located
 * @param (DoubleVect& les_valeurs) the array of values of the component of the field at the specified points
 * @param (int) the index of the component of the field to calculate
 * @return (DoubleVect&) the array of values of the component of the field at the specified points
 */
DoubleVect& Champ_base::valeur_aux_elems_compo(const DoubleTab&,
                                               const IntVect&,
                                               DoubleVect& les_valeurs,
                                               int ) const
{
  Cerr << que_suis_je();
  Cerr << "::valeur_aux_elems_compo(...) is not coded " << finl ;
  exit();
  return les_valeurs;
}

DoubleTab& Champ_base::valeur_aux_elems_smooth(const DoubleTab&,
                                               const IntVect& ,
                                               DoubleTab& les_valeurs)
{
  Cerr << que_suis_je();
  Cerr << "::valeur_aux_elems_smooth(...) is not coded " << finl ;
  Cerr << "The chsom option of probes does apply for the moment only for the field of type P1NC in VEF " << finl ;
  exit();
  return les_valeurs;
}

DoubleVect& Champ_base::valeur_aux_elems_compo_smooth(const DoubleTab&,
                                                      const IntVect&,
                                                      DoubleVect& les_valeurs,
                                                      int )
{
  Cerr << que_suis_je();
  Cerr << "::valeur_aux_elems_compo_smooth(...) is not coded " << finl ;
  Cerr << "The chsom option of probes does apply for the moment only for the field of type P1NC in VEF " << finl ;
  exit();
  return les_valeurs;
}

/*! @brief Time update.
 *
 * DOES NOTHING (in the base class)
 *
 * @param (double) update time
 */
DoubleVect& Champ_base::valeur_a_sommet(int sommet, const Domaine& dom, DoubleVect& val) const
{
  DoubleVect position(dimension);
  for(int i=0; i<dimension; i++)
    {
      position(i)=dom.coord(sommet,i);
    }
  return valeur_a(position, val);
}

/*! @brief Returns the compo-th coordinate of the values at the element le_poly at the vertex sommet
 *
 */
double Champ_base::valeur_a_sommet_compo(int sommet, int le_poly, int compo) const
{
  Cerr << que_suis_je();
  Cerr << "::valeur_a_sommet_compo(...) is not coded " << finl ;
  exit();
  return -1;
}

/*! @brief Returns the values at the vertices of the Domain dom
 *
 */
DoubleTab& Champ_base::valeur_aux_sommets(const Domaine& dom, DoubleTab& val) const
{
  const DoubleTab& positions=dom.coord_sommets();
  IntVect les_polys(positions.dimension(0));
  dom.chercher_elements(positions, les_polys);
  return valeur_aux_elems(positions, les_polys, val);
}

/*! @brief Returns the compo-th value at the vertices of dom.
 *
 */
DoubleVect& Champ_base::valeur_aux_sommets_compo(const Domaine& dom,
                                                 DoubleVect& val, int compo) const
{
  const DoubleTab& positions=dom.coord_sommets();
  IntVect les_polys(positions.dimension(0));
  dom.chercher_elements(positions, les_polys);
  return valeur_aux_elems_compo(positions, les_polys, val, compo);
}

/*! @brief Returns the field value at the faces
 *
 */
DoubleTab& Champ_base::valeur_aux_faces(DoubleTab& result) const
{
  return valeur_aux(ref_cast(Domaine_VF, domaine_dis_base()).xv(), result);
}

/*! @brief Returns the field value at the boundary faces
 *
 */
DoubleTab Champ_base::valeur_aux_bords() const
{
  const DoubleTab& xv_bord = ref_cast(Domaine_VF, domaine_dis_base()).xv_bord();
  DoubleTrav result(xv_bord.dimension_tot(0), valeurs().line_size());
  return valeur_aux(xv_bord, result);
}

/*! @brief Update of the base class Champ_base: does nothing!
 *
 */
void Champ_base::mettre_a_jour(double)
{

}
void Champ_base::abortTimeStep()
{
}

/*! @brief Assign a field to another.
 *
 * Returns the result of the assignment.
 *
 * @param (Champ_base& ch) right side of the assignment
 * @return (Champ_base&) the result of the assignment (*this)
 */
Champ_base& Champ_base::affecter(const Champ_base& ch)
{
  affecter_(ch);
  valeurs().echange_espace_virtuel();
  return *this;
}

// Factorize the error message for fields
// whose affecter method has no meaning.
void Champ_base::affecter_erreur()
{
  Nom message;
  message ="The method ";
  message += que_suis_je();
  message += "::affecter has no sense.\n";
  message += "TRUST has caused an error and will stop.\nUnexpected error during TRUST calculation.";
  Process::exit(message);
}

/*! @brief This method will set the units and the name of components, it is not actually const!!!
 *
 */
void Champ_base::corriger_unite_nom_compo()
{
  if(unite_.size() != nb_compo_)
    {
      Noms& unit=ref_cast_non_const(Noms,unite_);
      Nom nom=unite_[0];
      unit.dimensionner_force(nb_compo_);
      for(int compo=0; compo<nb_compo_; compo++)
        unit[compo] = nom;
    }
  int compo_nommees=(nom_compo(0) != Nom());
  if((!compo_nommees)&&(nb_compo_ != 1))
    {
      Noms ext(nb_compo_);
      if (nb_compo_<=3)
        {
          if(!axi)
            {
              ext[0]=le_nom()+"X";
              if (nb_compo_> 1) ext[1]=le_nom()+"Y";
              if (nb_compo_> 2) ext[2]=le_nom()+"Z";
            }
          else
            {
              ext[0]=le_nom()+"R";
              if (nb_compo_> 1) ext[1]=le_nom()+"teta";
              if (nb_compo_> 2) ext[2]=le_nom()+"Z";
            }
        }
      else for (int c=0; c<nb_compo_; c++)
          ext[c]=le_nom()+"_"+Nom(c);
      Champ_base& ch=ref_cast_non_const(Champ_base,*this);
      ch.fixer_noms_compo(ext);
    }
}

DoubleTab& Champ_base::eval_elem(DoubleTab& tab_valeurs) const
{
  Process::exit("This function is only for DG");

  return tab_valeurs;
}


void Champ_base::calculer_valeurs_elem_post(DoubleTab& les_valeurs,int nb_elem,Nom& nom_post,const Domaine& dom) const
{
  //nom_post=le_nom();
  Nom nom_dom=dom.le_nom();
  Nom nom_dom_inc= dom.le_nom();
  if(sub_type(Champ_Inc_base, *this) )
    {
      nom_dom_inc=ref_cast(Champ_Inc_base, *this).domaine().le_nom();
    }
  else if(sub_type(Champ_Fonc_base, *this) )
    {
      nom_dom_inc=ref_cast(Champ_Fonc_base, *this).domaine().le_nom();
    }

  bool isChamp_basis_function = is_basis_function();
  if (isChamp_basis_function)
    {
      int ndim = is_vectorial() ? Objet_U::dimension : 1;
      les_valeurs.resize(nb_elem, ndim);
    }
  else
    les_valeurs.resize(nb_elem, nb_compo_);

  if(nom_dom==nom_dom_inc)
    {
      valeur_aux_centres_de_gravite(dom, les_valeurs);
    }
  else
    {
      DoubleTrav centres_de_gravites(nb_elem, dimension);
      dom.calculer_centres_gravite(centres_de_gravites);
      valeur_aux(centres_de_gravites, les_valeurs);
    }


  if((axi) && (nb_compo_==dimension))
    {
      DoubleTrav centres_de_gravites(nb_elem, dimension);
      dom.calculer_centres_gravite(centres_de_gravites);
      double teta, vR, vT;
      for (int num_elem=0; num_elem<nb_elem; num_elem++)
        {
          teta = centres_de_gravites(num_elem,1);
          vR=les_valeurs(num_elem, 0);
          vT=les_valeurs(num_elem, 1);
          les_valeurs(num_elem, 0) =vR*cos(teta)-vT*sin(teta);
          les_valeurs(num_elem, 1) =vR*sin(teta)+vT*cos(teta);
        }
    }

  nom_post+= Nom("_elem_");
  nom_post+= nom_dom;
}

void Champ_base::calculer_valeurs_elem_compo_post(DoubleTab& les_valeurs,int ncomp,int nb_elem,Nom& nom_post,const Domaine& dom) const
{
  //nom_post=nom_compo(ncomp);
  Nom nom_dom=dom.le_nom();
  Nom nom_dom_inc= dom.le_nom();
  if(sub_type(Champ_Inc_base, *this) )
    {
      nom_dom_inc=ref_cast(Champ_Inc_base, *this).domaine().le_nom();
    }
  else if(sub_type(Champ_Fonc_base, *this) )
    {
      nom_dom_inc=ref_cast(Champ_Fonc_base, *this).domaine().le_nom();
    }
  ToDo_Kokkos("Critical; rewrite as Champ_base::calculer_valeurs_elem_post");
  DoubleTrav centres_de_gravites(nb_elem, dimension);
  les_valeurs.resize(nb_elem);
  if(nom_dom==nom_dom_inc)
    {
      if(sub_type(Champ_Inc_base, *this) )
        {
          const Domaine_VF& zvf = ref_cast(Domaine_VF,ref_cast(Champ_Inc_base, *this).equation().domaine_dis());
          // To avoid a resize by nb_elem_tot per call to xp()
          for (int i=0; i<nb_elem; i++)
            for (int k=0; k<dimension; k++)
              centres_de_gravites(i,k) = zvf.xp(i,k);
        }
      else
        dom.calculer_centres_gravite(centres_de_gravites);

      IntVect les_polys(nb_elem);
      {
        for(int elem=0; elem<nb_elem; elem++)
          {
            les_polys(elem)=elem;
          }
      }
      valeur_aux_elems_compo(centres_de_gravites, les_polys, les_valeurs, ncomp);
    }
  else
    {
      dom.calculer_centres_gravite(centres_de_gravites);
      valeur_aux_compo(centres_de_gravites, les_valeurs, ncomp);
    }
  nom_post+= Nom("_elem_");
  nom_post+= nom_dom;
}

// Adds the contribution of other processors to values and counter
inline void add_sommets_communs(const Domaine& dom, DoubleTab& les_valeurs, IntTab& compteur)
{

  //  if (Process::nproc()>9) return;
  char* theValue = getenv("TRUST_POST_SOM_NON_PARA");
  if (theValue != nullptr)
    return;
  int nb_compo_ = les_valeurs.line_size();

  // Iterate over joints
  const Joints& joints = dom.faces_joint();
  const int nb_joints = joints.size();
  for (int i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      const Joint& joint = joints[i_joint];
      const int PEvoisin = joint.PEvoisin();
      const Joint_Items& joint_item = joint.joint_item(JOINT_ITEM::SOMMET);
      const ArrOfInt& sommets_communs = joint_item.items_communs();

      // Temporary arrays
      IntTab envoie_sommets;
      DoubleTab envoie_valeurs;
      IntTab envoie_compteur;
      IntTab recoit_sommets;
      DoubleTab recoit_valeurs;
      IntTab recoit_compteur;
      int size=0;
      // Iterate over common vertices
      int size_sommets_communs = sommets_communs.size_array();
      for (int j=0; j<size_sommets_communs; j++)
        {
          int sommet = sommets_communs[j];
          // If this common vertex belongs to the list of Dirichlet vertices:
          if ( compteur(sommet)>0)
            {
              // The index of the remote neighboring vertex is:
              int sommet_voisin = joint_item.renum_items_communs()(j,0);
              // Verify that the local vertex is indeed a vertex:
              assert(joint_item.renum_items_communs()(j,1)==sommet);
              // Resize the send arrays and fill them
              size++;
              envoie_sommets.resize(size);
              envoie_sommets(size-1) = sommet_voisin;
              envoie_compteur.resize(size);
              envoie_compteur(size-1) = compteur(sommet);
              envoie_valeurs.resize(size, nb_compo_);
              for(int compo=0; compo<nb_compo_; compo++)
                envoie_valeurs(size-1,compo) = les_valeurs(sommet,compo);
            }
        }
      if (Process::me()<PEvoisin)
        {
          // Send arrays to PEvoisin
          envoyer(envoie_sommets,Process::me(),PEvoisin,Process::me()+1000);
          envoyer(envoie_valeurs,Process::me(),PEvoisin,Process::me()+2000);
          envoyer(envoie_compteur,Process::me(),PEvoisin,Process::me()+3000);
          // Receive arrays from PEvoisin
          recevoir(recoit_sommets,PEvoisin,Process::me(),PEvoisin+1000);
          recevoir(recoit_valeurs,PEvoisin,Process::me(),PEvoisin+2000);
          recevoir(recoit_compteur,PEvoisin,Process::me(),PEvoisin+3000);
        }
      else
        {
          // Receive arrays from PEvoisin
          recevoir(recoit_sommets,PEvoisin,Process::me(),PEvoisin+1000);
          recevoir(recoit_valeurs,PEvoisin,Process::me(),PEvoisin+2000);
          recevoir(recoit_compteur,PEvoisin,Process::me(),PEvoisin+3000);
          // Send arrays to PEvoisin
          envoyer(envoie_sommets,Process::me(),PEvoisin,Process::me()+1000);
          envoyer(envoie_valeurs,Process::me(),PEvoisin,Process::me()+2000);
          envoyer(envoie_compteur,Process::me(),PEvoisin,Process::me()+3000);

        }
      // Add the contribution of PEvoisin to the values and counter arrays
      // if the counters on both sides are non-zero
      int recoit_size = recoit_sommets.size_array();
      for (int i=0; i<recoit_size; i++)
        {
          int sommet = recoit_sommets(i);

          // Contribution received from a Dirichlet vertex
          if (recoit_compteur(i))
            {
              // If the local vertex is not Dirichlet, zero out values
              if (compteur(sommet)==0)
                for(int compo=0; compo<nb_compo_; compo++)
                  les_valeurs(sommet,compo)=0;

              compteur(sommet)+=recoit_compteur(i);
              for(int compo=0; compo<nb_compo_; compo++)
                les_valeurs(sommet,compo)+=recoit_valeurs(i,compo);
            }
        }
    }
}

void Champ_base::calculer_valeurs_som_post(DoubleTab& les_valeurs,int nb_som,Nom& nom_post,const Domaine& dom) const
{
  Nom nom_dom=dom.le_nom();
  Nom nom_dom_inc= dom.le_nom();
  if(sub_type(Champ_Inc_base, *this) )
    {
      nom_dom_inc=ref_cast(Champ_Inc_base, *this).domaine().le_nom();
    }
  else if(sub_type(Champ_Fonc_base, *this) )
    {
      nom_dom_inc=ref_cast(Champ_Fonc_base, *this).domaine().le_nom();
    }

  const DoubleTab& coord_sommets=dom.coord_sommets() ;

  les_valeurs.resize(nb_som, nb_compo_);

  if(nom_dom==nom_dom_inc)
    {
      valeur_aux_sommets(dom,les_valeurs);
    }
  else
    {
      valeur_aux(coord_sommets, les_valeurs);
    }
  int old_traitement_symetrie=0;
  {
    char* theValue = getenv("TRUST_POST_SOM_SYMETRIE_ERREUR");
    if (theValue != nullptr)
      {
        Cerr<<"results depend on order of bc or faces in bc ..."<<finl;
        old_traitement_symetrie=1;
      }
  }
  int impose_cl_diri=1;

  {
    char* theValue = getenv("TRUST_POST_SOM_NO_DIRICHLET");
    if (theValue != nullptr)
      {
        impose_cl_diri=0;
      }
  }

  // Taking boundary conditions into account:
  if (sub_type(Champ_Inc_base, *this)&&impose_cl_diri)
    {
      const Champ_Inc_base& chi=ref_cast(Champ_Inc_base, *this);
      if (!chi.mon_equation_non_nul())
        {
          Cerr<<"no equation associated to "<<que_suis_je()<<finl;
          impose_cl_diri=0;
        }
    }
  if (sub_type(Champ_Inc_base, *this)&&impose_cl_diri)
    {
      const Champ_Inc_base& chi=ref_cast(Champ_Inc_base, *this);

      const Equation_base& eqn=chi.equation();
      // GF we do not want to take into account the BCs in EF
      if (eqn.discretisation().que_suis_je()!="EF")
        if((eqn.inconnue().le_nom() == le_nom())
            && (sub_type(Champ_Inc_base, *this))
            && (dom==(ref_cast(Champ_Inc_base, *this).equation().probleme().domaine() ) ) )
          {
            IntTab compteur(dom.nb_som());
            compteur = 0;

            const Domaine_Cl_dis_base& zcl=eqn.domaine_Cl_dis();
            int nb_cond_lim=zcl.nb_cond_lim(),num_cl;
            for (num_cl=0; num_cl<nb_cond_lim; num_cl++)
              {
                const Cond_lim& la_cl=zcl.les_conditions_limites(num_cl);
                const Frontiere_dis_base& frontiere_dis=la_cl->frontiere_dis();
                const Frontiere& frontiere=frontiere_dis.frontiere();
                const Faces& faces=frontiere.faces();
                // modif bm: boundary faces now have virtual faces,
                //  we must not loop over virtual faces (old code: nb_faces_tot())
                int nb_faces=faces.nb_faces();
                int nb_som_faces=faces.nb_som_faces();
                if(sub_type(Dirichlet, la_cl.valeur()))
                  {
                    const Dirichlet& diri=ref_cast(Dirichlet, la_cl.valeur());
                    for(int num_face=0; num_face<nb_faces ; num_face++)
                      for(int num_som=0; num_som<nb_som_faces; num_som++)
                        {
                          int sommet=faces.sommet(num_face, num_som);

                          if (sommet < 0) continue;

                          for(int compo=0; compo<nb_compo_; compo++)
                            {
                              if (compteur(sommet) == 0)
                                les_valeurs(sommet, compo) = 0;

                              les_valeurs(sommet, compo) += diri.val_imp(num_face, compo);

                            }
                          compteur(sommet) += 1;
                        }
                  }
                else if(sub_type(Dirichlet_homogene, la_cl.valeur()))
                  {
                    for(int num_face=0; num_face<nb_faces ; num_face++)
                      for(int num_som=0; num_som<nb_som_faces; num_som++)
                        {
                          int sommet=faces.sommet(num_face, num_som);

                          if (sommet < 0) continue;

                          for(int compo=0; compo<nb_compo_; compo++)
                            {
                              if (compteur(sommet) == 0)
                                les_valeurs(sommet, compo) = 0;

                              les_valeurs(sommet, compo) += 0;

                            }
                          compteur(sommet) += 1;
                        }
                  }
                else if((sub_type(Symetrie, la_cl.valeur()))&&old_traitement_symetrie)
                  {
                    if(nb_compo_==dimension)
                      {
                        ArrOfDouble normale(dimension);
                        DoubleTab delta(dimension-1,dimension);
                        int nb_faces_tot=faces.nb_faces_tot();
                        for(int num_face=0; num_face<nb_faces_tot ; num_face++)
                          {
                            // compute the normal:
                            int sommet0=faces.sommet(num_face, 0);
                            for (int k=1; k<dimension; k++)
                              {
                                int sommet=faces.sommet(num_face, k);

                                if (sommet < 0) continue;

                                for(int compo=0; compo<dimension; compo++)
                                  delta(k-1,compo)=coord_sommets(sommet,compo)-coord_sommets(sommet0,compo);
                              }
                            if(dimension==2)
                              {
                                normale[0]=-delta(0,1);
                                normale[1]=delta(0,0);
                              }
                            else if(dimension==3)
                              {
                                normale[0]=delta(0,1)*delta(1,2) - delta(0,2)*delta(1,1);
                                normale[1]=delta(0,2)*delta(1,0) - delta(0,0)*delta(1,2);
                                normale[2]=delta(0,0)*delta(1,1) - delta(0,1)*delta(1,0);
                              }
                            else
                              {
                                Cerr << "We do not know treating the dimension : " << dimension << finl;
                                exit();
                              }
                            normale *= 1. / norme_array(normale);
                            for(int num_som=0; num_som<nb_som_faces; num_som++)
                              {
                                int sommet=faces.sommet(num_face, num_som);

                                if (sommet < 0) continue;

                                if (sommet<dom.nb_som()) // Real vertex
                                  {
                                    double psc=0;
                                    for(int k=0; k< dimension; k++)
                                      psc+=normale[k]*les_valeurs(sommet,k);
                                    for(int compo=0; compo<nb_compo_; compo++)
                                      les_valeurs(sommet, compo) -= psc*normale[compo];
                                  }
                              }
                          }
                      }
                  }
              }
            // Adds the contribution of other processors to values and counter
            add_sommets_communs(dom, les_valeurs, compteur);

            // Finish computing the average for cases that were modified
            // by a Dirichlet boundary condition
            int nb_som_l = dom.nb_som();
            //assert (nb_som_l==nb_som_l); // GF not sure about the assert, just to check
            for (int sommet = 0; sommet<nb_som_l; sommet++)
              for(int compo=0; compo<nb_compo_; compo++)
                if (compteur(sommet) != 0)
                  les_valeurs(sommet,compo) /= compteur(sommet);

          }
    }

  // Taking axisymmetry into account:
  if((axi) && (nb_compo_==dimension))
    {
      double teta, vR, vT;
      for(int sommet=0; sommet<nb_som; sommet++)
        {
          teta = coord_sommets(sommet,1);
          vR=les_valeurs(sommet, 0);
          vT=les_valeurs(sommet, 1);
          les_valeurs(sommet, 0) =vR*cos(teta)-vT*sin(teta);
          les_valeurs(sommet, 1) =vR*sin(teta)+vT*cos(teta);
        }
    }

  nom_post+= Nom("_som_");
  nom_post+= nom_dom;
}

void Champ_base::calculer_valeurs_som_compo_post(DoubleTab& les_valeurs,int ncomp,int nb_som,Nom& nom_post,const Domaine& dom,int appliquer_cl) const
{
  Nom nom_dom=dom.le_nom();
  Nom nom_dom_inc= dom.le_nom();
  if(sub_type(Champ_Inc_base, *this) )
    {
      nom_dom_inc=ref_cast(Champ_Inc_base, *this).domaine().le_nom();
    }
  else if(sub_type(Champ_Fonc_base, *this) )
    {
      nom_dom_inc=ref_cast(Champ_Fonc_base, *this).domaine().le_nom();
    }
  const DoubleTab& coord_sommets=dom.coord_sommets() ;
  les_valeurs.resize(nb_som);
  if(nom_dom==nom_dom_inc)
    {
      valeur_aux_sommets_compo(dom,les_valeurs, ncomp);
    }
  else
    {
      valeur_aux_compo(coord_sommets, les_valeurs, ncomp);
    }

  //int impose_cl_diri=1;

  {
    char* theValue = getenv("TRUST_POST_SOM_NO_DIRICHLET");
    if (theValue != nullptr)
      {
        appliquer_cl=0;
      }
  }

  if (appliquer_cl)
    {
      // List of boundary vertices in contact with a Dirichlet face:
      // Taking boundary conditions into account:
      if (sub_type(Champ_Inc_base, *this))
        {
          const Champ_Inc_base& chi=ref_cast(Champ_Inc_base, *this);
          const Equation_base& eqn=chi.equation();
          if((eqn.inconnue().le_nom() == le_nom())
              && (sub_type(Champ_Inc_base, *this))
              && (dom==(ref_cast(Champ_Inc_base, *this).equation().probleme().domaine() ) ) )
            {

              IntTab compteur(dom.nb_som());
              compteur = 0;
              int num_cl;
              const Domaine_Cl_dis_base& zcl=eqn.domaine_Cl_dis();
              int nb_cond_lim=zcl.nb_cond_lim();
              for (num_cl=0; num_cl<nb_cond_lim; num_cl++)
                {
                  const Cond_lim& la_cl=zcl.les_conditions_limites(num_cl);
                  const Frontiere_dis_base& frontiere_dis=la_cl->frontiere_dis();
                  const Frontiere& frontiere=frontiere_dis.frontiere();
                  const Faces& faces=frontiere.faces();
                  int nb_faces=faces.nb_faces_tot();
                  int nb_som_faces=faces.nb_som_faces();

                  if(sub_type(Dirichlet, la_cl.valeur()))
                    {
                      const Dirichlet& diri=ref_cast(Dirichlet, la_cl.valeur());
                      for(int num_face=0; num_face<nb_faces ; num_face++)
                        {
                          for(int num_som=0; num_som<nb_som_faces; num_som++)
                            {
                              int sommet=faces.sommet(num_face, num_som);
                              if (compteur(sommet) == 0)
                                les_valeurs(sommet) = 0;

                              les_valeurs(sommet) += diri.val_imp(num_face,ncomp);
                              compteur(sommet) += 1;
                            }
                        }
                    }
                  else if(sub_type(Dirichlet_homogene, la_cl.valeur()))
                    {
                      const Dirichlet_homogene& diri=ref_cast(Dirichlet_homogene, la_cl.valeur());
                      for(int num_face=0; num_face<nb_faces ; num_face++)
                        {
                          for(int num_som=0; num_som<nb_som_faces; num_som++)
                            {
                              int sommet=faces.sommet(num_face, num_som);

                              les_valeurs(sommet) = diri.val_imp(num_face,ncomp);
                              compteur(sommet) += 1;

                            }
                        }
                    }
                }
              // Adds the contribution of other processors to les_valeurs and counter
              add_sommets_communs(dom, les_valeurs, compteur);

              // Finish computing the average for cases that were modified
              // by a Dirichlet boundary condition
              int nb_som_l = dom.nb_som();
              assert (nb_som==nb_som_l); // GF not sure about the assert, just to check
              for (int sommet = 0; sommet<nb_som_l; sommet++)
                if (compteur(sommet) != 0)
                  les_valeurs(sommet) /= compteur(sommet);
            }
        }
    }
  nom_post+= Nom("_som_");
  nom_post+= nom_dom;
}

int Champ_base::completer_post_champ(const Domaine& dom,const int is_axi,const Nom& loc_post,
                                     const Nom& le_nom_champ_post,Format_Post_base& format) const
{

  const Nature_du_champ& la_nature = nature_du_champ();
  const int nb_compo = nb_comp();
  const Noms& noms_composante = noms_compo();

  format.completer_post(dom,is_axi,la_nature,nb_compo,noms_composante,loc_post,le_nom_champ_post);

  return 1;
}

// Sets the conditions for taking boundary conditions into account
// during field filtering/postprocessing operations
void Champ_base::completer(const Domaine_Cl_dis_base& zcl)
{
}

/*! @brief Sets the time at which the field is defined
 *
 * @param (double& t) the new time at which the field is defined
 * @return (double) the new time at which the field is defined
 */
double Champ_base::changer_temps(const double t)
{
  return temps_ = t ;
}

/*! @brief Returns the time of the field
 *
 * @return (double) the time of the field
 */
double Champ_base::temps() const
{
  return temps_ ;
}



int Champ_base::fixer_nb_valeurs_nodales(int n)
{
  Cerr << "Champ_base::fixer_nb_valeurs_nodales\n ";
  Cerr << que_suis_je() << " does not have nb_valeurs_nodales" << finl;
  exit();
  return 0;
}

// Not all fields have a domaine_dis_base. For those that
// do not, the call is invalid.
void Champ_base::associer_domaine_dis_base(const Domaine_dis_base& domaine_dis)
{
  Cerr << "Error in Champ_base::associer_domaine_dis_base\n";
  Cerr << " (field name : " << le_nom() << ")\n";
  Cerr << " The method " << que_suis_je();
  Cerr << "::associer_domaine_dis_base is not coded\n";
  Cerr << " or the field does not possess a domaine_dis_base." << finl;
  exit();
}

// Not all fields have a domaine_dis_base. For those that
// do not, the call is invalid.
const Domaine_dis_base& Champ_base::domaine_dis_base() const
{
  Cerr << "Error in Champ_base::domaine_dis_base\n";
  Cerr << " (field name : " << le_nom() << ")\n";
  Cerr << " The method " << que_suis_je();
  Cerr << "::domaine_dis_base is not coded\n";
  Cerr << " or the field does not possess a domaine_dis_base." << finl;
  exit();
  throw;
}


/*! @brief Computes the trace of a field on a boundary at time tps
 *
 *     WEC:
 *     The boundary passed as a parameter must be part of the domain
 *     on which the field is based.
 *     The result is computed on this boundary and the size of
 *     the DoubleTab x corresponds to the number of faces of the boundary.
 *     x can have a virtual space; the trace function calls
 *     echange_espace_virtuel.
 *
 *     Special case (unfortunately) of Champ_P0_VDF:
 *     If the boundary is a connector, the result is computed on the
 *     associated connector. In this case, the DoubleTab x must be
 *     dimensioned on the associated connector.
 *
 *
 * @param (Frontiere_dis_base&) discretized boundary on which we want to compute the trace of the field at time tps
 * @param (DoubleTab& x , double tps) the values of the field on the boundary at time tps
 * @return (DoubleTab&) the values of the field on the boundary at time tps
 */
DoubleTab& Champ_base::trace(const Frontiere_dis_base& , DoubleTab& x , double tps,int distant) const
{
  Cerr << que_suis_je() << "did not overloaded Champ_base::trace" << finl;
  exit();
  return x;
}
