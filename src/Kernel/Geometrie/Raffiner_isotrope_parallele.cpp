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

#include <Raffiner_isotrope_parallele.h>
#include <Reordonner_faces_periodiques.h>
#include <LecFicDistribue.h>
#include <EcrFicCollecte.h>
#include <FichierHDFPar.h>
#include <LecFicDiffuse.h>
#include <TRUSTArrays.h>
#include <Schema_Comm.h>
#include <TRUSTTabs.h>
#include <SFichier.h>
#include <Domaine.h>
#include <Scatter.h>
#include <Param.h>
#include <Nom.h>

Implemente_instanciable( Raffiner_isotrope_parallele, "Raffiner_isotrope_parallele", Raffiner_Simplexes ) ;

Sortie& Raffiner_isotrope_parallele::printOn( Sortie& os ) const
{
  Raffiner_Simplexes::printOn( os );
  return os;
}

Entree& Raffiner_isotrope_parallele::readOn( Entree& is )
{
  Raffiner_Simplexes::readOn( is );
  return is;
}


Entree&  Raffiner_isotrope_parallele::interpreter(Entree& is)
{
  int form=0;
  bool format_hdf = false;
  Nom org,newd;
  Param param(que_suis_je());

  // XD Raffiner_isotrope_parallele interprete Raffiner_isotrope_parallele BRACE Refine parallel mesh in parallel
  param.ajouter("name_of_initial_domaines|name_of_initial_zones",&org,Param::REQUIRED); // XD_ADD_P chaine
  // XD_CONT name of initial Domaines
  param.ajouter("name_of_new_domaines|name_of_new_zones",&newd,Param::REQUIRED); // XD_ADD_P chaine
  // XD_CONT name of new Domaines
  param.ajouter("ascii",&form);  // XD_ADD_P flag
  // XD_CONT writing Domaines in ascii format
  param.ajouter_flag("single_hdf",&format_hdf); // XD_ADD_P rien
  // XD_CONT writing Domaines in hdf format
  param.lire_avec_accolades(is);

  // Force un fichier unique au dela d'un certain nombre de rangs MPI:
  if (Process::force_single_file(Process::nproc(), org+".Zones"))
    format_hdf = true;
  int binaire=!form;
  if (form && format_hdf)
    {
      Cerr << "Raffiner_isotrope_parallele::interpreter(): options 'ascii' and 'single_hdf' are mutually exclusive!" << finl;
      Process::exit(1);
    }
  Domaine dom_org;
  Noms& liste_bords_periodiques = dom_org.bords_perio();

  org+=".Zones";

  Nom copy(org);
  copy = copy.nom_me(Process::nproc(), "p", 1);
  //bool is_hdf = FichierHDF::is_hdf5(copy);
  LecFicDiffuse test;
  bool is_hdf = test.ouvrir(copy) && FichierHDF::is_hdf5(copy);
  bool has_perio = false;

  if (!is_hdf)
    {
      LecFicDistribue  fichier;
      fichier.set_bin(binaire);
      fichier.ouvrir(org);
      dom_org.readOn_has_perio(fichier, has_perio);
      dom_org.set_fichier_lu(org);
      if (!has_perio) // Old (pre TRUST 1.9.8) Domain format - periodic boundaries stored after:
        fichier >> liste_bords_periodiques;
    }
  else
    {
      FichierHDFPar fic_hdf;
      org = copy;
      fic_hdf.open(org, true);
      Entree_Brute data;
      fic_hdf.read_dataset("//zone", Process::me(), data);
      // Feed TRUST objects:
      dom_org.readOn_has_perio(data, has_perio);
      dom_org.set_fichier_lu(org);
      if (!has_perio) // Old (pre TRUST 1.9.8) Domain format - periodic boundaries stored after:
        data >> liste_bords_periodiques;
      fic_hdf.close();
    }

  Scatter::uninit_sequential_domain(dom_org);
  Domaine dom_new(dom_org);
  dom_new.typer(dom_org.type_elem()->que_suis_je());

  refine_domain(dom_org,dom_new);

  // After spliting the mesh and the boundaries, we reorder perdiodic faces:
  for (auto nom_bord : liste_bords_periodiques)
    {
      Cerr << "Reordering faces of the periodic boundary " << nom_bord << finl;
      ArrOfDouble direction_perio;
      Reordonner_faces_periodiques::chercher_direction_perio(direction_perio, dom_new, nom_bord);
      Bord& bord = dom_new.bord(nom_bord);
      IntTab& faces = bord.faces().les_sommets();
      Reordonner_faces_periodiques::reordonner_faces_periodiques(dom_new, faces, direction_perio, Objet_U::precision_geom);
    }

  if (nproc() > 1)
    {
      Scatter::uninit_sequential_domain(dom_new);
      int nb_sommet_avant_completion=dom_new.nb_som();
      Scatter::trier_les_joints(dom_new.faces_joint());

      // Rebuild the correspondance between vertices, knowing that the number of common items
      // has potentially changed.
      statistics().begin_count(STD_COUNTERS::parallel_meshing,statistics().get_last_opened_counter_level()+1);
      Scatter::construire_correspondance_sommets_par_coordonnees(dom_new, true /* allow resize of items_communs */);
#ifndef NDEBUG
      // In debug mode, we ensure exact matching of the updated items_communs by running the exchange a second time
      // without expecting a resize of items_communs:
      Scatter::construire_correspondance_sommets_par_coordonnees(dom_new, false);
#endif
      double maxtime = mp_max(statistics().get_time_since_last_open(STD_COUNTERS::parallel_meshing));
      statistics().end_count(STD_COUNTERS::parallel_meshing);
      Cerr << "Scatter::construire_correspondance_sommets_par_coordonnees fin, time:"
           << maxtime
           << finl;

      statistics().begin_count(STD_COUNTERS::parallel_meshing,statistics().get_last_opened_counter_level()+1);
      Scatter::construire_structures_paralleles(dom_new);
      maxtime = mp_max(statistics().get_time_since_last_open(STD_COUNTERS::parallel_meshing));
      statistics().end_count(STD_COUNTERS::parallel_meshing);
      Cerr << "Scatter::construire_structures_paralleles, time:" << maxtime << finl;

      int ecrit=1;
      if (ecrit)
        {
          int nb_elem_reel=dom_new.nb_elem();
          Scatter::uninit_sequential_domain(dom_new);
          dom_new.les_elems().resize( nb_elem_reel,dom_new.les_elems().dimension(1));
          dom_new.les_sommets().resize(nb_sommet_avant_completion,dimension);

          Scatter::uninit_sequential_domain(dom_new);
          newd+=".Zones";

          if( !format_hdf )
            {
              EcrFicCollecte os;
              os.set_bin(binaire);
              os.ouvrir(newd);
              if (!binaire)
                {
                  os.setf(ios::scientific);
                  os.precision(Objet_U::format_precision_geom);
                }
              os << dom_new;
            }
          else
            {
              Sortie_Brute os_hdf;
              os_hdf << dom_new;
              FichierHDFPar fic_hdf;
              newd = newd.nom_me(Process::nproc(), "p", 1);
              fic_hdf.create(newd);
              fic_hdf.create_and_fill_dataset_MW("/zone", os_hdf);
              fic_hdf.close();
            }
        }
    }
  else
    Scatter::init_sequential_domain(dom_new);

  return is;
}
