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

#include <Reordonner_faces_periodiques.h>
#include <Domaine_Cl_dis_base.h>

#include <Domaine_VF.h>
#include <Periodique.h>
#include <Domaine.h>
#include <Scatter.h>

Implemente_instanciable(Periodique, "Periodique", Cond_lim_base);
// XD periodic condlim_base periodique INHERITS_BRACE 1). For Navier-Stokes equations, this keyword is used to indicate
// XD_CONT that the horizontal inlet velocity values are the same as the outlet velocity values, at every moment. As
// XD_CONT regards meshing, the inlet and outlet edges bear the same name.; 2). For scalar transport equation, this
// XD_CONT keyword is used to set a periodic condition on scalar. The two edges dealing with this periodic condition
// XD_CONT bear the same name.

Sortie& Periodique::printOn(Sortie& s) const { return s << que_suis_je() << finl; }

Entree& Periodique::readOn(Entree& s)
{
  le_champ_front.typer("Champ_front_uniforme");
  return s;
}

void Periodique::completer()
{
  Frontiere& frontiere = frontiere_dis().frontiere();
  Cerr << "Initialization for periodic on " << frontiere.le_nom() << finl;

  // Search for the periodicity direction:
  ArrOfDouble erreur;
  int ok = Reordonner_faces_periodiques::check_faces_periodiques(frontiere, direction_perio_, erreur, true /* verbose */);
  if (!ok)
    exit();

  distance_ = norme_array(direction_perio_);
  int i;
  const int dim = direction_perio_.size_array();
  direction_xyz_ = -2;
  for (i = 0; i < dim; i++)
    {
      if (std::fabs(direction_perio_[i]) > precision_geom)
        {
          if (direction_xyz_ == -2)
            direction_xyz_ = i;
          else
            // Second non-zero coordinate: direction vector is not aligned on an axis
            direction_xyz_ = -1;
        }
    }
  if (direction_xyz_ < 0)
    Cerr << "Periodic direction not aligned on an axis" << finl;
  else
    Cerr << "Periodic direction aligned on the axis " << direction_xyz_ << finl;

  const Domaine& domaine = frontiere.domaine();

  // Create an index array spanning all boundary faces
  IntTab tab_face_associee;
  const Domaine_VF& domainevf = ref_cast(Domaine_VF, domaine_Cl_dis().domaine_dis());
  domainevf.creer_tableau_faces_bord(tab_face_associee, RESIZE_OPTIONS::NOCOPY_NOINIT);
  tab_face_associee = -1;

  // Number of virtual faces on this boundary:
  const int nb_faces = frontiere.nb_faces();
  const int nb_faces_2_ = nb_faces / 2;
  const ArrOfInt& faces_virt = frontiere.get_faces_virt();
  const int nb_faces_virt = faces_virt.size_array();

  // Fill the real part of the face_associee array for the boundary of interest:
  // and the real part of the "associee" array
  const int i_premiere_face = frontiere.num_premiere_face();
  for (i = 0; i < nb_faces_2_; i++)
    {
      const int i1 = i_premiere_face + i;
      const int i2 = i_premiere_face + i + nb_faces_2_;
      tab_face_associee[i1] = i2;
      tab_face_associee[i2] = i1;
    }
  const MD_Vector& md_faces_front = tab_face_associee.get_md_vector();
  // Exchange virtual space with index translation:
  Scatter::construire_espace_virtuel_traduction(md_faces_front, md_faces_front, tab_face_associee, 1 /* fatal errors */);
  // Array giving for each virtual face of the domain -1 if it is not a
  // boundary face, otherwise its index in the boundaries.
  const ArrOfInt& ind_faces_virt_bord = domaine.ind_faces_virt_bord();
  // Create an array giving, for each virtual face of the boundaries (all boundaries)
  // the index of the face in the current periodic boundary (-1 otherwise)
  const int nb_faces_front_tot = tab_face_associee.size_totale();
  ArrOfInt index(nb_faces_front_tot);
  index = -2;
  const int nb_faces_domaine = domaine_Cl_dis().domaine_dis().face_sommets().dimension(0);
  for (i = 0; i < nb_faces_virt; i++)
    {
      const int face_domaine = frontiere.face_virt(i); // Index of a face in the domain
      const int face_front = ind_faces_virt_bord[face_domaine - nb_faces_domaine]; // Index in the boundaries
      index[face_front] = nb_faces + i;
    }

  // Fill the "face_front_associee_" array
  face_front_associee_.resize_array(nb_faces + nb_faces_virt);
  for (i = 0; i < nb_faces + nb_faces_virt; i++)
    {
      int resu = -1;
      if (i < nb_faces_2_)
        resu = i + nb_faces_2_;
      else if (i < nb_faces)
        resu = i - nb_faces_2_;
      else
        {
          // Virtual face:
          const int face_domaine = frontiere.face_virt(i - nb_faces); // Index of the face in the domain
          const int face_front = ind_faces_virt_bord[face_domaine - nb_faces_domaine]; // Index in the boundaries
          // Index of the associated face in the face_associee array:
          const int face_front_associee = tab_face_associee[face_front];
          if (face_front_associee >= 0)
            {
              // Index of the associated face in the boundary:
              const int face_asso = index[face_front_associee];
              assert(face_asso >= 0);
              resu = face_asso;
            }
          else
            {
              // The associated virtual face is not in the domain
              Cerr << "Error in Periodique::completer()" << finl;
              exit();
            }
        }
      face_front_associee_[i] = resu;
    }
}

int Periodique::direction_periodicite() const
{
  if (!est_periodique_selon_un_axe())
    {
      Cerr << "Error in Periodique::direction_periodicite():\n" << " An algorithm seems to assume that the periodic direction is aligned in X, Y or Z\n" << " and this is not the case !" << finl;
      exit();
    }
  return direction_xyz_;
}
