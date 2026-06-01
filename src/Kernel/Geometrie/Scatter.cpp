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


#include <Scatter.h>
#include <Domaine.h>
#include <LecFicDistribueBin.h>
#include <TRUSTTabs.h>
#include <Connectivite_som_elem.h>
#include <Schema_Comm.h>
#include <Faces_builder.h>
#include <Domaine_VF.h>
#include <Reordonner_faces_periodiques.h>
#include <communications.h>
#include <MD_Vector_tools.h>
#include <MD_Vector_std.h>
#include <MD_Vector_seq.h>
#include <unistd.h> // PGI
#include <Poly_geom_base.h>
#include <Entree_Brute.h>
#include <Comm_Group_MPI.h>
#include <FichierHDFPar.h>
#include <LecFicDiffuse.h>
#include <Format_Post_Lata.h>
#include <EFichierBin.h>
#include <Array_tools.h>
#include <Perf_counters.h>
#include <vector>
#include <numeric>

Implemente_instanciable(Scatter,"Scatter",Interprete);
// XD scatter interprete scatter NO_BRACE Class to read a partionned mesh from the files during a parallel calculation.
// XD_CONT The files are in binary format.
// XD attr file chaine file REQ Name of file.
// XD attr domaine ref_domaine domaine REQ Name of domain.

/*! @brief Simple call to: Interprete::printOn(Sortie&)
 *
 * @param os Output stream.
 * @return The modified output stream.
 */
Sortie& Scatter::printOn(Sortie& os) const
{
  return Interprete::printOn(os);
}


/*! @brief Simple call to: Interprete::readOn(Entree&)
 *
 * @param is Input stream.
 * @return The modified input stream.
 */
Entree& Scatter::readOn(Entree& is)
{
  return Interprete::readOn(is);
}

/*! @brief Returns the associated domain.
 *
 * @return The associated domain.
 */
Domaine& Scatter::domaine()
{
  return le_domaine.valeur();
}

namespace
{
// For debug:
void dump_lata(const Domaine& dom)
{
  Format_Post_Lata post;  // Lata V2
  Nom nom_fichier_lata("espaces_virtuels");

  const int nb_joints = dom.nb_joints();
  constexpr int IS_FIRST = 1;

  post.initialize_lata(nom_fichier_lata, Format_Post_Lata::BINAIRE, Format_Post_Lata::SINGLE_FILE);
  post.ecrire_entete(0.0, 0, IS_FIRST);
  post.ecrire_domaine(dom, IS_FIRST);
  post.ecrire_temps(0.0);

  Noms units, noms_compo;
  units.add("");
  noms_compo.add("I");
  DoubleTab data(dom.nb_elem());
  for(int ij = 0; ij < nb_joints; ij++)
    {
      const ArrOfInt& t1 = dom.joint(ij).joint_item(JOINT_ITEM::ELEMENT).items_distants();
      data = 0.;
      const int nt1 = t1.size_array();
      for (int i = 0; i < nt1; i++) data[t1[i]] += 1;

      post.ecrire_champ(dom,
                        units,
                        noms_compo,
                        1,           // ncomp,
                        0.0,         // time,
                        Nom("partition") + Nom(dom.joint(ij).PEvoisin()), // field_id,
                        dom.le_nom(), // domain_id
                        "ELEM",       // localisation,
                        "scalar",     // nature,
                        data          // values
                       );
    }
}
} // end anonymous namespace

/*! @brief Reads and completes a parallel domain according to the keywords read in the data set.
 *
 * Format:
 *    Scatter [debug] file_name domain_name
 *   Reads the vertices, elements and joint vertices and faces,
 *   builds the distant and virtual spaces according to
 *   the joint layer thickness.
 */
Entree& Scatter::interpreter(Entree& is)
{
  // Name of partition files: nomentree.xxxx
  Nom nomentree;
  is >> nomentree;
  if (Process::is_sequential())
    {
      Motcle n(nomentree);
      if (n != ";" && n != "unlock;")
        {
          Cerr << "Error ! You ran a sequential calculation and can't use Scatter keyword here. Run a parallel calculation or remove this keyword." << finl;
          exit();
        }
      Cerr << "Scatter: preparing domain structure\n"
           << " (this is workaround for bugged domain operators that don't do it)" << finl;
      Nom nomdomaine;
      is >> nomdomaine;
      Objet_U& obj = objet(nomdomaine);
      if(!sub_type(Domaine, obj))
        {
          Cerr << "obj : " << obj << " is not an object of type Domain !" << finl;
          exit();
        }
      Domaine& dom = ref_cast(Domaine, obj);
      if (n == ";")
        init_sequential_domain(dom);
      else
        uninit_sequential_domain(dom);
      return is;
    }
  // For debugging on linux in parallel
#ifdef linux
  static int gdb_non_lance=1;
  char* TRUST_GDB=getenv("TRUST_GDB");
  if (gdb_non_lance && ((Motcle)nomentree=="DEBUG" || TRUST_GDB!=nullptr))
    {
      gdb_non_lance=0;
      if ((Motcle)nomentree=="DEBUG") is >> nomentree;
      if (je_suis_maitre())
        {
          Cerr << "Enter \"return\" to this window after" << finl;
          Cerr << "typing \"cont\" in other gdb windows." << finl;
          Cerr << (int)system ("sh -c read ok") << finl;
        }
      else
        {
          Nom getpidn((int)getpid());
          Nom cmdfile=getpidn;
          Nom command0="echo attach ";
          command0+=getpidn;
          command0+=" > ";
          command0+=cmdfile;
          Cerr << (int)system(command0) << finl;
          command0=" ls -l /proc/";
          command0+=getpidn;
          command0+="/exe | awk '{print $NF}' > execname";
          Cerr << (int)system(command0) << finl;
          Nom command="[ -f /usr/X11R6/bin/xterm ] && x=\"/usr/X11R6/bin/xterm -exec gdb -x \";";
          command+="[ -f /usr/bin/konsole ] && x=\"/usr/bin/konsole -e gdb -x \";";
          command+="$x ";
          command+=cmdfile;
          command+=" `cat execname` ";
          command+=" &";
          Cerr<<"command: " <<command<<finl;
          Cerr << (int)system(command) << finl;
        }
    }
#endif
  barrier();

  if (Process::je_suis_maitre())
    Cerr << "Execution of the Scatter module." << finl;

  statistics().begin_count(STD_COUNTERS::interprete_scatter,statistics().get_last_opened_counter_level()+1);
  // Retrieve the domain:
  Nom nomdomaine;
  is >> nomdomaine;
  Objet_U& obj = objet(nomdomaine);
  if(!sub_type(Domaine, obj))
    {
      Cerr << "Error in Scatter: object of type '" << obj.que_suis_je() << "' when Domaine was expected!" << finl;
      exit();
    }
  Domaine& dom = ref_cast(Domaine, obj);
  le_domaine = dom;

  // Read the partition files:
  barrier();
  if (Process::je_suis_maitre())
    Cerr << "Reading the domain" << finl;

  lire_domaine(nomentree);

  barrier();
  Cerr << "Calculation of renum_items_communs for the nodes" << finl;
  calculer_renum_items_communs(dom.faces_joint(), JOINT_ITEM::SOMMET);

  // Not yet coded: verify that common vertices have identical coordinates
  // on all processors.
  // check_sommets_joints(dom);

  barrier();
  Cerr << "Construire_structures_paralleles" << finl;
  construire_structures_paralleles(dom);

  if (0)
    dump_lata(dom);

  barrier();
  Cerr << "End Distribue_domaines" << finl;

  Cerr << "\nQuality of partitioning --------------------------------------------" << finl;
  trustIdType total_nb_elem = Process::mp_sum(dom.nb_elem());
  Cerr << "\nTotal nb of elements = " << total_nb_elem << finl;
  Cerr << "Number of Domaines : " << Process::nproc() << finl;
  double min_element_domaine = mp_min(dom.nb_elem());
  double max_element_domaine = mp_max(dom.nb_elem());
  double mean_element_domaine = (double)(total_nb_elem / Process::nproc());
  Cerr << "Min number of elements on a Domaine = " << min_element_domaine << finl;
  Cerr << "Max number of elements on a Domaine = " << max_element_domaine << finl;
  Cerr << "Mean number of elements per Domaine = " << (int)(mean_element_domaine) << finl;
  double load_imbalance = max_element_domaine / mean_element_domaine;
  Cerr << "Load imbalance = " << load_imbalance << "\n" << finl;

  Elem_geom_base& elem=dom.type_elem().valeur();
  if (sub_type(Poly_geom_base,elem))
    ref_cast(Poly_geom_base,elem).compute_virtual_index();
  if(Process::me()==0)
    {
      double temps = statistics().get_time_since_last_open(STD_COUNTERS::interprete_scatter);
      Cerr << "Scatter time : " << temps << finl;
    }
  statistics().end_count(STD_COUNTERS::interprete_scatter);
  return is;
}

/*! @brief Merged domains receive joint information from their neighbours to ensure that their common items (vertices) appear in the same order
 *
 *  If it's not the case, the merged domain reorders its common items so that it matches the neighbour's order
 *  When 2 neighbouring domains have each been merged,
 *  only the processor with the lowest rank proceeds to reordering
 */
void Scatter::check_consistancy_remote_items(Domaine& dom, const ArrOfInt& mergedDomaines)
{
  const Joints& joints     = dom.faces_joint();
  const int nb_joints = joints.size();

  const DoubleTab& coords = dom.les_sommets();
  ArrOfInt liste_send;
  ArrOfInt liste_recv;



  const int moi = Process::me();
  const int myDomaineWasMerged = mergedDomaines[moi];

  for (int i_joint = 0; i_joint < nb_joints; i_joint++)
    {

      const int pe_voisin = joints[i_joint].PEvoisin();
      const int neighbourDomaineWasMerged = mergedDomaines[pe_voisin];
      if(myDomaineWasMerged && neighbourDomaineWasMerged)
        {
          if(pe_voisin < moi)
            liste_recv.append_array(pe_voisin);
          else
            liste_send.append_array(pe_voisin);
        }
      else if(myDomaineWasMerged && !neighbourDomaineWasMerged)
        liste_recv.append_array(pe_voisin);
      else if(!myDomaineWasMerged && neighbourDomaineWasMerged)
        liste_send.append_array(pe_voisin);
      else
        {
          //nothing to exchange
        }
    }

  DoubleTabs coord_items_locaux(nb_joints);
  DoubleTabs coord_items_distants(nb_joints);
  for (int i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      const Joint& joint     = joints[i_joint];
      const ArrOfInt& items_communs = joint.joint_item(JOINT_ITEM::SOMMET).items_communs();
      const int nb_items_communs = items_communs.size_array();

      DoubleTab&   coord   = coord_items_locaux[i_joint];
      coord.resize(nb_items_communs, dimension);
      for (int i = 0; i < nb_items_communs; i++)
        for (int j = 0; j < dimension; j++)
          coord(i,j) = coords(items_communs[i], j);
    }

  // Send local coordinates to the neighbouring processor
  {
    Schema_Comm schema_comm;
    schema_comm.set_send_recv_pe_list(liste_send, liste_recv);
    schema_comm.begin_comm();
    for (int i = 0; i < nb_joints; i++)
      {
        const int pe_voisin = joints[i].PEvoisin();
        const int neighbourDomaineWasMerged = mergedDomaines[pe_voisin];
        if( neighbourDomaineWasMerged && !(myDomaineWasMerged && pe_voisin<moi) )
          {
            Sortie& buffer = schema_comm.send_buffer(pe_voisin);
            buffer << coord_items_locaux[i];
          }
      }
    schema_comm.echange_taille_et_messages();

    if(myDomaineWasMerged)
      {
        for (int i = 0; i < nb_joints; i++)
          {
            const int pe_voisin = joints[i].PEvoisin();
            const int neighbourDomaineWasMerged = mergedDomaines[pe_voisin];
            if(!(neighbourDomaineWasMerged && pe_voisin>moi))
              {
                Entree& buffer = schema_comm.recv_buffer(pe_voisin);
                buffer >> coord_items_distants[i];
              }
          }
      }

    schema_comm.end_comm();
  }

  // check if the vertices in my joints appear in the same order as my neighbour joint
  if(myDomaineWasMerged)
    {
      for (int i_joint = 0; i_joint < nb_joints; i_joint++)
        {

          const int pe_voisin = joints[i_joint].PEvoisin();
          const int neighbourDomaineWasMerged = mergedDomaines[pe_voisin];
          if(neighbourDomaineWasMerged && pe_voisin>moi)
            continue;
          ArrOfInt& items_communs = dom.faces_joint()[i_joint].set_joint_item(JOINT_ITEM::SOMMET).set_items_communs();
          const ArrOfInt old_items_communs = joints[i_joint].joint_item(JOINT_ITEM::SOMMET).items_communs();
          const int     nb_items      = items_communs.size_array();
          const DoubleTab& coord_voisin = coord_items_distants[i_joint];
          const DoubleTab& my_coord = coord_items_locaux[i_joint];
          assert(my_coord.size_array() ==  coord_voisin.size_array());
          for(int i=0; i<nb_items; i++)
            {
              for(int j=0; j<nb_items; j++)
                {
                  int ok=1;
                  for (int dir=0; dir<Objet_U::dimension; dir++)
                    ok=ok&&(est_egal(coord_voisin(i,dir),my_coord(j,dir)));
                  if (ok)
                    {
                      items_communs[i] = old_items_communs[j];
                      break;
                    }
                }
            }
        }
    }
}


/*! @brief Does the exact same thing as the readOn of the class Domaine but without collective communication
 *
 *  Necessary when the processors don't have the same numbers of file to read
 */
void Scatter::read_domain_no_comm(Entree& fic, bool& read_perio)
{
  Domaine& dom = le_domaine.valeur();

  Cerr << "\treading vertices..." << finl;
  Domaine dom_tmp_for_vertices;
  dom_tmp_for_vertices.read_vertices(fic);

  Cerr << "\tDone !\n\treading elem infos (domaines)..." << finl;

  Nom accouverte="{";
  Motcle nom;
  fic >> nom;
  Domaine domaine_read;
  if(nom!=(const char*)"vide")
    {
      if (nom!=accouverte)
        Process::exit("Error: Scatter::read_domain_no_comm() -- One expected an opened bracket { to start.");
      domaine_read.read_former_domaine(fic, read_perio);
    }
  else
    Process::exit("Error: Scatter::read_domain_no_comm() -- Empty list ?! Should not happen?");
  Cerr << "Done!" << finl;

  //
  // Now merge the read domaine with the current domain
  //
  int nb_elems = dom.nb_elem();
  IntVect nums;
  Scatter::uninit_sequential_domain(dom);
  // Complete domain with new nodes and/or renumber nodes when we have doublons
  dom.ajouter(dom_tmp_for_vertices.les_sommets(), /*out*/ nums);
  if (nb_elems > 0)
    {
      domaine_read.renum(nums);
      domaine_read.renum_joint_common_items(nums, nb_elems);
    }

  // Merge domaine_read into current domain, w/o taking care of the joints.
  dom.merge_wo_vertices_with(domaine_read);

  if(nb_elems > 0)  // Current domain already had something, so joints will need update
    // Otherwise, joints were already read by "domaine_read.read_former_domaine(fic);" above and joints are OK.
    {
      //merging common vertices and remote items
      const int nb_joints = domaine_read.nb_joints();
      for (int i_joint = 0; i_joint < nb_joints; i_joint++)
        {
          const Joint& joint_to_add  = domaine_read.faces_joint()[i_joint];

          int my_joint_index = 0;
          while(joint_to_add.PEvoisin() != dom.faces_joint()[my_joint_index].PEvoisin())
            my_joint_index++;

          const ArrOfInt& sommets_to_add = joint_to_add.joint_item(JOINT_ITEM::SOMMET).items_communs();
          ArrOfInt& items_communs = dom.faces_joint()[my_joint_index].set_joint_item(JOINT_ITEM::SOMMET).set_items_communs();

          for(int index=0; index<sommets_to_add.size_array(); index++)
            items_communs.append_array(sommets_to_add[index]); // sommets_to_add is already renumbered with 'nums' - see call to renum_joint_common_items above
          array_trier_retirer_doublons(items_communs);

          const ArrOfInt& elements_to_add = joint_to_add.joint_item(JOINT_ITEM::ELEMENT).items_distants();
          ArrOfInt& items_distants = dom.faces_joint()[my_joint_index].set_joint_item(JOINT_ITEM::ELEMENT).set_items_distants();

          for(int index=0; index<elements_to_add.size_array(); index++)
            items_distants.append_array(elements_to_add[index]); // idem
        }
    }
}

/*! @brief Reads the domain from the file named "nomentree", of type LecFicDistribueBin or LecFicDistribue
 *
 *   Expected format: Domaine::ReadOn
 */
void Scatter::lire_domaine(Nom& nomentree)
{
  // Determine whether the file is in the new or old format
  if (Process::je_suis_maitre())
    Cerr << "Reading geometry from .Zones file(s) ..." << finl;
  barrier(); // Wait for the message to be displayed

  Domaine& dom = domaine();
  Noms& liste_bords_periodiques = dom.bords_perio();

  // Just in case - some dataset improperly build a Domain and then try to Scatter on it ...:
  dom.clear();

  Nom copy(nomentree);
  copy = copy.nom_me(Process::nproc(), "p", 1);

  LecFicDiffuse test;
  bool is_hdf = test.ouvrir(copy) && FichierHDF::is_hdf5(copy);
  if (test.ouvrir(nomentree) && FichierHDF::is_hdf5(nomentree))
    {
      Cerr << "Error: You probably made a single_hdf partitioning and using the wrong name of .Zones files in the scatter" << finl;
      Cerr << "You should remove '_p" << Process::nproc() << "' from the name of .Zones file (" << nomentree << ") in your datafile" << finl;
      Process::exit();
    }

  statistics().begin_count(STD_COUNTERS::read_scatter,statistics().get_last_opened_counter_level()+1);
  ArrOfInt mergedDomaines(Process::nproc());
  mergedDomaines = 0;
  bool domain_not_built = true;
  bool read_perio = false;
  if (is_hdf)
    {
      FichierHDFPar fic_hdf;

      nomentree = copy;
      fic_hdf.open(nomentree, true);

      std::string dname = "/zone_"  + std::to_string(Process::me());
      bool ok = fic_hdf.exists(dname.c_str());
      if(!ok)
        {
          mergedDomaines = 1;
          for(int i=0; i<Process::nproc(); i++)
            {
              Entree_Brute data_part;
              std::string tmp = dname + "_" + std::to_string(i);

              bool exists = fic_hdf.exists(tmp.c_str());
              if(exists)
                {
                  Nom dataset_name(dname);

                  fic_hdf.read_dataset(dataset_name, i, data_part);
                  read_domain_no_comm(data_part, read_perio);

                  // Record which file the domain was read from
                  dom.set_fichier_lu(nomentree);
                  if (!read_perio)  // are the periodic boundaries read from the Domain (new format) or after it?
                    data_part >> liste_bords_periodiques;
                  domain_not_built = false;
                }
              else
                break;

            }
        }
      else
        {
          Entree_Brute data;
          fic_hdf.read_dataset("/zone", Process::me(), data);

          // Feed TRUST objects:
          read_domain_no_comm(data, read_perio);
          dom.set_fichier_lu(nomentree);
          if (!read_perio)  // are the periodic boundaries read from the Domain (new format) or after it?
            data >> liste_bords_periodiques;
          domain_not_built = false;
        }

      fic_hdf.close();
    }
  else  // Not HDF
    {
      LecFicDistribueBin fichier_binaire;
      int isSingleDomaine = fichier_binaire.ouvrir(nomentree);
      if (!isSingleDomaine)
        {
          mergedDomaines = 1;
          Nom nomentree_with_suffix=nomentree.nom_me(Process::me());
          for(int i=0; i<Process::nproc(); i++)
            {
              EFichierBin fichier_binaire_part;
              std::string tmp = nomentree_with_suffix.getPrefix(".Zones").getString();
              tmp += "_";
              tmp += std::to_string(i);
              tmp += ".Zones";
              Nom nomentree_part(tmp);
              int ok = fichier_binaire_part.ouvrir(nomentree_part);
              if(ok)
                {
                  read_domain_no_comm(fichier_binaire_part, read_perio);

                  // Record which file the domain was read from
                  dom.set_fichier_lu(nomentree);
                  if (!read_perio)  // are the periodic boundaries read from the Domain (new format) or after it?
                    fichier_binaire_part >> liste_bords_periodiques;
                  fichier_binaire_part.close();
                  domain_not_built = false;
                }
              else
                break;
            }
        }
      else
        {
          read_domain_no_comm(fichier_binaire, read_perio);

          // Record which file the domain was read from
          dom.set_fichier_lu(nomentree);
          if (!read_perio)  // are the periodic boundaries read from the Domain (new format) or after it?
            fichier_binaire >> liste_bords_periodiques;
          fichier_binaire.close();
          domain_not_built = false;
        }
    }

  if(domain_not_built)
    {
      Cerr << "Error in Scatter::lire_domaine\n";
      Cerr << "The domain on the current process hasn't been built" << finl;
      Cerr << "The number of processes you mentionned is probaly higher than the number of domaines" << finl;
      Process::exit();
    }

  // Sanity check: number of processors = number of domains
  // (verify that there is no joint with a non-existent processor)
  // (the previous check is insufficient:
  // it only verifies that the number of processors does not exceed the number of domains)
  {

    const Joints& joints = dom.faces_joint();
    const int nb_joints = joints.size();
    int max_pe_voisin = 0;
    for (int i = 0; i < nb_joints; i++)
      {
        const int pe_voisin = joints[i].PEvoisin();
        if (pe_voisin >= max_pe_voisin)
          max_pe_voisin = pe_voisin;
      }

    max_pe_voisin = (int) mp_max(max_pe_voisin);
    double ok=1;
    if (max_pe_voisin >= nproc()) ok=0;
    if (!ok)
      {
        Cerr << "Error in Scatter::lire_domaine\n"
             << "The domain has been partitioned with at least " << max_pe_voisin << " "
             << "domaines whereas the number of processes asked is " << Process::nproc() << "." << finl;
        Cerr << "The number of domaines and number of processes must match." << finl;
        exit();
      }
  }

  // sort joints in increasing processor order
  Joints& joints = dom.faces_joint();
  trier_les_joints(joints);
  envoyer_all_to_all(mergedDomaines, mergedDomaines);
  check_consistancy_remote_items( dom, mergedDomaines );
  dom.check_domaine();

  // PL: not entirely exact: the displayed vertex count counts joint vertices multiple times...
  trustIdType nbsom = mp_sum(dom.les_sommets().dimension(0));
  Cerr << " Number of nodes: " << nbsom << finl;

  init_sequential_domain(dom);

  // merged domains need to reorder faces of periodic borders
  const int myDomaineWasMerged = mergedDomaines[Process::me()];
  if(myDomaineWasMerged)
    {
      for(auto& itr : liste_bords_periodiques)
        {
          Nom bp_nom = itr;
          Bord& bord = dom.bord(bp_nom);
          if(bord.nb_faces() == 0)
            continue;

          ArrOfDouble direction_perio(dimension);
          Reordonner_faces_periodiques::chercher_direction_perio(direction_perio, dom, bp_nom);
          IntTab& faces = bord.faces().les_sommets();
          double epsilon = precision_geom;
          Reordonner_faces_periodiques::reordonner_faces_periodiques(dom, faces, direction_perio, epsilon);
        }
    }
  statistics().end_count(STD_COUNTERS::read_scatter);
  barrier();
}

/*! @brief Builds the parallel structures of the domain (determination of distant elements as a function of joint thickness,
 *
 *    determination of distant vertices,
 *    creation of virtual vertices and elements)
 *
 */
void Scatter::construire_structures_paralleles(Domaine& dom)
{
  // First: remove the "sequential" structures associated with vertices and elements during reading:
  {
    MD_Vector md_nul;
    dom.les_sommets().set_md_vector(md_nul);
    dom.les_elems().set_md_vector(md_nul);
  }

  const Noms& liste_bords_periodiques = dom.bords_perio();

  // The call order is important:
  calculer_espace_distant_elements(dom);

  if (liste_bords_periodiques.size() > 0)
    corriger_espace_distant_elements_perio(dom);

  calculer_nb_items_virtuels(dom.faces_joint(), JOINT_ITEM::ELEMENT);

  // Determination of distant vertices from distant elements
  calculer_espace_distant_sommets(dom);

  // Creation of distant/virtual spaces and common items for vertex and element arrays:
  DoubleTab& sommets = dom.les_sommets();
  IntTab& elements = dom.les_elems();
  MD_Vector md_sommets, md_elements;
  construire_md_vector(dom, sommets.dimension(0), JOINT_ITEM::SOMMET, md_sommets);
  construire_md_vector(dom, elements.dimension(0), JOINT_ITEM::ELEMENT, md_elements);
  MD_Vector_tools::creer_tableau_distribue(md_sommets, sommets);
  sommets.echange_espace_virtuel();
  construire_espace_virtuel_traduction(md_elements /* type index */,
                                       md_sommets /* type valeur */,
                                       elements);
  // Reorder the joint faces (implicit correspondence with the neighboring pe)
  reordonner_faces_de_joint(dom);
}

/*! @brief Sort joints by increasing neighbor proc number
 */
void Scatter::trier_les_joints(Joints& joints)
{
  const int nb_joints = joints.size();
  ArrOfInt pe_voisins(nb_joints);
  for (int i = 0; i < nb_joints; i++)
    pe_voisins[i] = joints[i].PEvoisin();
  pe_voisins.ordonne_array();
  // Copy the joint list
  Joints anciens_joints(joints);
  for (int i = 0; i < nb_joints; i++)
    {
      // Process the neighbor processor pe_voisin:
      const int pe_voisin = pe_voisins[i];
      // Where is the joint with this processor in the old list?
      int i_old;
      for (i_old = 0; i_old < nb_joints; i_old++)
        if (anciens_joints[i_old].PEvoisin() == pe_voisin)
          break;
      assert(i_old < nb_joints);
      joints[i] = anciens_joints[i_old];
    }
}

// If a joint with "pe" exists, return its index,
// otherwise create a new joint and return its index.
static int ajouter_joint(Domaine& domaine, int pe)
{
  Joints& joints = domaine.faces_joint();
  const int i_joint = joints.size();

  {
    for (int i = 0; i < i_joint; i++)
      if (joints[i].PEvoisin() == pe)
        return i;
  }

  Joint& joint = joints.add(Joint());
  joint.nommer(Nom("Joint_")+Nom(pe));
  joint.associer_domaine(domaine);
  int ep = (i_joint > 0) ? joints[0].epaisseur() : 1;
  joint.affecte_epaisseur(ep);
  joint.affecte_PEvoisin(pe);

  // Initialise all arrays of additional joints.
  // Note BM: to be thorough, only
  // certain arrays should be initialised (those already initialised for
  // existing joints), but that is more complex to do...
  {
    for (int t = 0; t < 5; t++)
      {
        JOINT_ITEM type;
        switch(t)
          {
          case 0:
            type = JOINT_ITEM::SOMMET;
            break;
          case 1:
            type = JOINT_ITEM::ELEMENT;
            break;
          case 2:
            type = JOINT_ITEM::FACE;
            break;
          case 3:
            type = JOINT_ITEM::ARETE;
            break;
          case 4:
            type = JOINT_ITEM::FACE_FRONT;
            break;
          default:
            Cerr << "Error in Scatter.cpp : ajouter_joint" << finl;
            // To avoid the following warning on gcc 3.4:
            // Scatter.cpp:416: warning: 'type' might be used uninitialized in this function
            type = JOINT_ITEM::SOMMET;
            Process::exit();
          }
        Joint_Items& data = joint.set_joint_item(type);
        data.set_items_communs();
        data.set_items_distants();
        data.set_nb_items_virtuels(0);
        data.set_renum_items_communs().resize(0,2);
      }
  }

  return i_joint;
}


/*! @brief Determines the distant items from a list of items to send and lists of common items.
 *
 *   Example:
 *    calculer_espace_distant_sommets
 *    calculer_espace_distant_faces
 *   For vertices: the "items_to_send" are the vertices of distant elements.
 *    If processor A wants processor B to know vertex i,
 *    the processor that owns the vertex must send it to B.
 *    The "owning" processor is the smallest among the PEs
 *    sharing this vertex (common item) (required to perform
 *    echange_item_commun and echange_espace_virtuel in a single pass).
 *    Furthermore, if several processors request that the same vertex
 *    be sent to the same processor, it must only be inserted once in the
 *    distant space.
 *
 * @param (joints) the joints in which the distant space is to be computed
 * @param (nb_items_reels) the number of real items (vertices, faces, ...)
 * @param (items_to_send) a vector of "nproc()" arrays, for each processor, the list of items to send (e.g. all vertices of distant elements, or all faces)
 * @param (type_item) the items whose distant space is to be computed
 */
void Scatter::calculer_espace_distant(Domaine&                  domaine,
                                      const int           nb_items_reels,
                                      const ArrsOfInt& items_to_send,
                                      const JOINT_ITEM type_item)
{
  assert(items_to_send.size() == Process::nproc());

  Process::Journal() << "Scatter::calculer_espace_distant type_item="
                     << (int)type_item << finl;

  Joints& joints = domaine.faces_joint();

  // First, determine for all items the PE owner number:
  //  For each item of the domain:
  //   column 0: index of the item on the owning PE (index_on_pe_owner)
  //   column 1: number of the owning PE (the smallest PE sharing this item)
  IntTab num_global_items(nb_items_reels, 2);
  {
    int i;
    const int moi = Process::me();
    for (i = 0; i < nb_items_reels; i++)
      {
        num_global_items(i, 0) = i;
        num_global_items(i, 1) = moi;
      }
    const int nb_joints = joints.size();
    for (int i_joint = 0; i_joint < nb_joints; i_joint++)
      {
        const Joint&   joint              = joints[i_joint];
        const int   pe_voisin          = joint.PEvoisin();
        const IntTab& renum_items_communs= joint.joint_item(type_item).renum_items_communs();
        const int   nb_items_communs   = renum_items_communs.dimension(0);
        for (i = 0; i < nb_items_communs; i++)
          {
            const int num_item_distant = renum_items_communs(i, 0);
            const int num_item_local   = renum_items_communs(i, 1);
            const int pe_actuel = num_global_items(num_item_local, 1);
            if (pe_voisin < pe_actuel)
              {
                num_global_items(num_item_local, 0) = num_item_distant;
                num_global_items(num_item_local, 1) = pe_voisin;
              }
          }
      }
  }

  Schema_Comm schema_comm;
  const int nproc = Process::nproc();

  // First step: send to the owning processor of the items
  // the list of items to be sent and to which processor they
  // must be sent.
  // If processor A must send element E to processor B,
  // and that element uses item S belonging to processor C,
  // then send to C the message:
  // "put item S in the distant space of processor B"

  // Prepare a communication scheme between neighbors:
  //  Sending processor: the processor that owns the distant element,
  //  Receiving processor: the processor that owns an item of the element.
  // These processors are neighbors through existing joints.
  const int nb_joints = joints.size();
  ArrOfInt liste_voisins(nb_joints);
  //int i_joint;
  for (int i_joint = 0; i_joint < nb_joints; i_joint++)
    liste_voisins[i_joint] = joints[i_joint].PEvoisin();

  schema_comm.set_send_recv_pe_list(liste_voisins, liste_voisins, 1 /* me_to_me */);
  schema_comm.begin_comm();
  {
    // and for each item to be sent, send to the processor owning the item (pe_item_owner):
    //  - its local index on that processor (item_distant),
    //  - the number of the processor to which it must be sent (pe_destination)
    const int nb_procs = Process::nproc();
    for (int pe_destination = 0; pe_destination < nb_procs; pe_destination++)
      {
        const ArrOfInt& items    = items_to_send[pe_destination];
        const int     nb_items = items.size_array();
        for (int i_item = 0; i_item < nb_items; i_item++)
          {
            const int item          = items[i_item];
            const int item_distant  = num_global_items(item, 0);
            const int pe_item_owner = num_global_items(item, 1);
            // Send the index of the distant item and which joint it
            // must be placed in.
            // If pe_joint == pe_destination, the item is necessarily already
            // known by the other processor; no need to send it
            if (pe_item_owner != pe_destination)
              schema_comm.send_buffer(pe_item_owner) << item_distant << pe_destination;
          }
      }
  }

  // Exchange messages
  schema_comm.echange_taille_et_messages();

  // Receive the distant items. Read all buffers and
  // store items in "items_distants" by destination processor.
  // For each neighbor processor, the list of distant items to send:
  ArrsOfInt items_distants(nproc);

  // Loop over all processors (pe_source) that sent me messages:
  // Loop over neighbor processors plus myself:
  for (int i_source = 0; i_source < nb_joints + 1; i_source++)
    {
      const int pe_source =
        (i_source < nb_joints) ? liste_voisins[i_source] : Process::me();

      Entree& buffer = schema_comm.recv_buffer(pe_source);
      // Loop "while the buffer is not empty"
      while(1)
        {
          int item_distant; // Index of the distant item
          int pe_distant;   // Number of the pe to which the item must be sent
          buffer >> item_distant >> pe_distant;
          if (buffer.eof())
            break;
          assert(pe_distant != Process::me());
          ArrOfInt& array = items_distants[pe_distant];
          array.append_array(item_distant);
        }
    }
  schema_comm.end_comm();

  // Remove duplicates and items already known by the neighbor processor:
  {
    // List of joints corresponding to each pe
    ArrOfInt joint_of_pe(nproc);
    joint_of_pe = -1;
    for (int i_joint = 0; i_joint < nb_joints; i_joint++)
      {
        const int pe = joints[i_joint].PEvoisin();
        joint_of_pe[pe] = i_joint;
      }
    // List of items already known by the neighbor processor (common items)
    // sorted in increasing order
    ArrOfInt items_communs_tri;

    for (int pe = 0; pe < nproc; pe++)
      {
        ArrOfInt& items = items_distants[pe];
        // Remove duplicates:
        array_trier_retirer_doublons(items);
        // Remove items already known:
        const int i_joint = joint_of_pe[pe];
        if (i_joint >= 0)
          {
            items_communs_tri =
              joints[i_joint].joint_item(type_item).items_communs();
            items_communs_tri.ordonne_array();
            array_retirer_elements(items, items_communs_tri);
          }
        else
          {
            // No common item with this pe.
          }
      }
  }

  // Distant spaces may be created on processors with
  // which no joint exists yet. Add the new joints.
  {
    ArrOfInt nouveaux_voisins;

    int i;
    for (i = 0; i < nproc; i++)
      if (items_distants[i].size_array() > 0)
        nouveaux_voisins.append_array(i);

    // Add the new joints
    ajouter_joints(domaine, nouveaux_voisins);
    Process::Journal() << " News joints created : (ArrOfInt) "
                       << nouveaux_voisins << finl;
  }

  Joints& joints_non_const = domaine.faces_joint();
  const int nb_new_joints = joints_non_const.size();
  // Fill the distant items arrays
  for (int i_joint = 0; i_joint < nb_new_joints; i_joint++)
    {
      Joint& joint = joints_non_const[i_joint];
      const int pe = joint.PEvoisin();
      ArrOfInt& joint_items_distants = joint.set_joint_item(type_item).set_items_distants();
      joint_items_distants = items_distants[pe];
      Process::Journal() << " Joint with PE:" << pe
                         << " Number of remote items : "
                         << joint_items_distants.size_array() << finl;
    }
  // Fill the number of virtual items
  calculer_nb_items_virtuels(joints_non_const, type_item);
}
inline Nom endian()
{
  int x = 1;
  if(*(char *)&x == 1)
    return "little-endian";
  else
    return "big-endian";
}
/*! @brief Adds joints with all PEs in pe_voisins.
 *
 * To make the set of joints symmetric,
 *   a joint is also created on the destination processor:
 *    If A adds a joint with B, then B adds a joint with A.
 *   Joints are sorted in ascending order of PE number.
 *   WARNING: joints are therefore reordered!
 *   pe_voisins is updated with the list of joints actually created.
 *
 */
void Scatter::ajouter_joints(Domaine& domaine,
                             ArrOfInt& pe_voisins)
{
  Joints& joints = domaine.faces_joint();
  ArrOfInt liste_pe;


  // Make joints symmetric (if A->B then B->A):
  {
    // Put in liste_pe the "transpose" of the pe_voisins list:
    // list of processors that have my number in their "pe_voisins".
    reverse_send_recv_pe_list(pe_voisins, liste_pe);
    const int n = liste_pe.size_array();
    // Concatenate the two lists.
    for (int i = 0; i < n; i++)
      pe_voisins.append_array(liste_pe[i]);
    array_trier_retirer_doublons(pe_voisins);
    liste_pe.resize_array(0);
  }
  // Remove from pe_voisins the PEs for which a joint already exists
  {
    const int n = joints.size();
    liste_pe.resize_array(n);
    for (int i = 0; i < n; i++)
      liste_pe[i] = joints[i].PEvoisin();
    array_retirer_elements(pe_voisins, liste_pe);
  }
  // Add new joints and sort in ascending order
  // As of 2/11/2005, Liste::inserer does not allow inserting
  // at the beginning of the list. Unusable. Brute-force method:
  {
    const int n = pe_voisins.size_array();
    for (int i = 0; i < n; i++)
      ajouter_joint(domaine, pe_voisins[i]);
    trier_les_joints(joints);
  }
}

/*! @brief Generic method to compute the remote space of a geometric item type (vertex, face, edge) based on the remote space of elements:
 *
 *   The remote "type_item" items (for type_item = vertex, face or edge) are
 *   the "type_item" items attached to remote elements.
 *   Example: remote vertices are all vertices of all remote elements.
 *  @sa
 *   Scatter::calculer_espace_distant_sommets
 *   Scatter::calculer_espace_distant_faces
 *
 * @param (domaine) the domain
 * @param (type_item) the type of items whose remote space is to be computed
 * @param (connectivite_elem_item) the array giving for each domain element the indices of its items. Only the real part of the array is used (logically, the virtual part does not exist yet). (e.g. domaine().les_elems() for type_item==SOMMET or domaine_VF().face_sommets() for type_item==FACE)
 * @param (nb_items_reels) the number of real "type_item" items
 * @param (items_lies) if the array is non-empty, it must have size nb_items_reels. In that case, it enforces the property: "if item i is remote, then item items_lies[i] is remote too". This array is used to include associated virtual periodic vertices. (see calculer_espace_distant_sommets).
 */
static void calculer_espace_distant_item(Domaine& le_dom,
                                         const JOINT_ITEM type_item,
                                         const IntTab& connectivite_elem_item,
                                         const int nb_items_reels,
                                         const ArrOfInt& items_lies)
{
  if(Process::is_sequential())
    return;

  const Joints& joints                 = le_dom.faces_joint();
  const int   nb_joints              = joints.size();
  const int   nproc                  = Process::nproc();
  const int   nb_items_par_element   = connectivite_elem_item.dimension(1);
  // The type_item items to send to each processor:
  ArrsOfInt items_to_send(nproc);
  // A temporary array;
  ArrOfInt liste_items;


  // Are there linked items?
  const int flag_items_lies = (items_lies.size_array() > 0);
  assert(flag_items_lies == 0 || items_lies.size_array() == nb_items_reels);


  for (int i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      const Joint&     joint          = joints[i_joint];
      const int     pe_voisin      = joint.PEvoisin();
      const ArrOfInt& esp_dist_elems = joint.joint_item(JOINT_ITEM::ELEMENT).items_distants();
      const int     nb_elems_dist  = esp_dist_elems.size_array();
      liste_items.resize_array(0);
      // Put in liste_items all items of all elements
      // that are in esp_dist_elems:
      for (int i_elem = 0; i_elem < nb_elems_dist; i_elem++)
        {
          const int elem = esp_dist_elems[i_elem];
          for (int i_item = 0; i_item < nb_items_par_element; i_item++)
            {
              const int item = connectivite_elem_item(elem, i_item);
              if (item>-1)
                {
                  liste_items.append_array(item);
                  // If an item is linked to the current item, also send the linked item.
                  if (flag_items_lies)
                    {
                      const int item_lie = items_lies[item];
                      if (item_lie != item)
                        {
                          assert(item_lie >= 0 && item_lie < nb_items_reels);
                          assert(items_lies[item_lie] == item_lie); // chaining of links is forbidden
                          liste_items.append_array(item_lie);
                        }
                    }
                }
            }
        }
      array_trier_retirer_doublons(liste_items);
      // These items must be sent to the neighboring processor:
      items_to_send[pe_voisin] = liste_items;
    }
  // Compute remote spaces based on "items_to_send"
  Scatter::calculer_espace_distant(le_dom, nb_items_reels, items_to_send, type_item);
}

/*! @brief Based on the remote space of elements, computes the remote space of vertices.
 *
 * For each joint, the set of vertices of all joint elements is sent to the neighboring processor.
 *   It is the processor that owns the vertex
 *   (the smallest PE that holds it) that puts it in its remote space.
 *   Warning: new joints are created.
 *   The following arrays are filled:
 *    dom.faces_joint(i).joint_item(JOINT_ITEM::SOMMET).items_distants();
 *
 */
void Scatter::calculer_espace_distant_sommets(Domaine& dom)
{
  if (Process::je_suis_maitre())
    Cerr << "Scatter::calculer_espace_distant_sommets : start" << finl;

  const IntTab& connectivite_elem_som = dom.les_elems();
  const int   nb_sommets_reels      = dom.nb_som();

  ArrOfInt renum_som_perio(nb_sommets_reels);
  // Initialize the renum_som_perio array
  for (int i = 0; i < nb_sommets_reels; i++)
    renum_som_perio[i] = i;
  Reordonner_faces_periodiques::renum_som_perio(dom, renum_som_perio,
                                                0 /* do not compute for virtual vertices */);

  calculer_espace_distant_item(dom,
                               JOINT_ITEM::SOMMET,
                               connectivite_elem_som,
                               nb_sommets_reels,
                               renum_som_perio);
}

/*! @brief Same as Scatter::calculer_espace_distant_sommets for faces.
 *
 */
void Scatter::calculer_espace_distant_faces(Domaine& domaine,
                                            const int nb_faces_reelles,
                                            const IntTab& elem_faces)
{
  if (Process::je_suis_maitre())
    Cerr << "Scatter::calculer_espace_distant_faces : start" << finl;

  ArrOfInt tableau_vide;

  calculer_espace_distant_item(domaine,
                               JOINT_ITEM::FACE,
                               elem_faces,
                               nb_faces_reelles,
                               tableau_vide);
}

/*! @brief Same as Scatter::calculer_espace_distant_sommets for edges.
 *
 */
void Scatter::calculer_espace_distant_aretes(Domaine& domaine,
                                             const int nb_aretes_reelles,
                                             const IntTab& elem_aretes)
{
  if (Process::je_suis_maitre())
    Cerr << "Scatter::calculer_espace_distant_aretes : start" << finl;
  ArrOfInt tableau_vide;
  calculer_espace_distant_item(domaine,
                               JOINT_ITEM::ARETE,
                               elem_aretes,
                               nb_aretes_reelles,
                               tableau_vide);
}

/*! @brief Assumes that each joint[i].joint_item(type_item).items_communs() contains the local indices of common joint items in the same
 *   order on both processors (local and neighbor).
 *   Fills renum_items_communs:
 *    column 0 = content of the items_communs array on the neighboring PE
 *    column 1 = content of the items_communs array on the local PE
 *
 */
void Scatter::calculer_renum_items_communs(Joints& joints,
                                           const JOINT_ITEM type_item)
{
  // It suffices to send the _faces array to the neighbor in order
  // so it has the face indices on the other PE.

  const int nb_joints = joints.size();
  int       i_joint;
  Schema_Comm  schema_comm;
  ArrOfInt liste_voisins(nb_joints);
  for (i_joint = 0; i_joint < nb_joints; i_joint++)
    liste_voisins[i_joint] = joints[i_joint].PEvoisin();
  schema_comm.set_send_recv_pe_list(liste_voisins, liste_voisins);

  schema_comm.begin_comm();

  for (i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      const Joint& joint     = joints[i_joint];
      const int  pe_voisin = joint.PEvoisin();
      const ArrOfInt& items_communs =
        joint.joint_item(type_item).items_communs();
      schema_comm.send_buffer(pe_voisin) << items_communs;
    }

  schema_comm.echange_taille_et_messages();

  // The common items array received from the neighboring PE:
  ArrOfInt items_communs_voisin;


  for (i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      Joint&           joint         = joints[i_joint];
      const int     pe_voisin     = joint.PEvoisin();
      const ArrOfInt& items_communs = joint.joint_item(type_item).items_communs();
      const int     nb_items      = items_communs.size_array();
      schema_comm.recv_buffer(pe_voisin) >> items_communs_voisin;

      assert(nb_items == items_communs_voisin.size_array());

      IntTab& renum_items_communs = joint.set_joint_item(type_item).set_renum_items_communs();
      renum_items_communs.resize(nb_items, 2);
      // The index of the joint face on the other PE is in tmp(i,1)
      for (int i = 0; i < nb_items; i++)
        {
          renum_items_communs(i,0) = items_communs_voisin[i];
          renum_items_communs(i,1) = items_communs[i];
        }
    }

  schema_comm.end_comm();
}

/*! @brief Builds an MD_Vector_std from the joint information of the domain for the requested item type.
 *
 */
void Scatter::construire_md_vector(const Domaine& dom, int nb_items_reels, const JOINT_ITEM type_item, MD_Vector& md_vector)
{
  if(Process::is_sequential())
    {
      MD_Vector_seq mdseq(nb_items_reels);
      md_vector.copy(mdseq);
      return;
    }

  const Joints& joints  = dom.faces_joint();
  const int nb_joints = joints.size();

  ArrOfInt pe_voisins(nb_joints);
  ArrsOfInt items_to_send(nb_joints);
  ArrsOfInt items_to_recv(nb_joints);
  ArrsOfInt blocs_to_recv(nb_joints);

  // flag indicating whether the (common) item is received from a processor
  ArrOfBit flags(nb_items_reels);
  flags = 0;

  int nitems_tot = nb_items_reels;
  const int moi = me();

  for (int i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      const int pe = joints[i_joint].PEvoisin();
      pe_voisins[i_joint] = pe;
      const Joint_Items& joint = joints[i_joint].joint_item(type_item);
      {
        // Processing common items
        const ArrOfInt& items_communs = joint.items_communs();
        const int n = items_communs.size_array();

        // Joints must arrive in ascending PE number order,
        // otherwise the following algorithm does not work:
        assert((i_joint == 0) || (pe > joints[i_joint-1].PEvoisin()));
        if (pe > moi)
          {
            // I must send these items to the neighboring processor
            ArrOfInt& dest = items_to_send[i_joint];

            for (int i = 0; i < n; i++)
              {
                const int item = items_communs[i];
                if (!flags[item])
                  {
                    // item not received from a processor
                    dest.append_array(item);
                  }
              }
          }
        else
          {
            // I receive this item from another processor
            ArrOfInt& dest = items_to_recv[i_joint];

            for (int i = 0; i < n; i++)
              {
                const int item = items_communs[i];
                if (!flags.testsetbit(item))
                  {
                    // item not yet received from a processor
                    dest.append_array(item);
                  }
              }
          }
      }
      // Processing remote and virtual items
      {
        const int nitems_virt = joint.nb_items_virtuels();
        ArrOfInt& dest = blocs_to_recv[i_joint];
        if (nitems_virt > 0)
          {
            dest.resize_array(2, RESIZE_OPTIONS::NOCOPY_NOINIT);
            // Definition of the virtual items block for the neighboring processor
            dest[0] = nitems_tot;
            dest[1] = nitems_tot + nitems_virt;
            nitems_tot += nitems_virt;
          }
      }
      {
        const ArrOfInt& items_distants = joint.items_distants();
        const int n = items_distants.size_array();
        ArrOfInt& dest = items_to_send[i_joint];
        const int index = dest.size_array();
        dest.resize_array(index + n, RESIZE_OPTIONS::COPY_NOINIT); // copy the old values!
        dest.inject_array(items_distants, n, index /* dest index */, 0 /* src index */);
      }
    }

  MD_Vector_std md(nitems_tot, nb_items_reels, pe_voisins, items_to_send, items_to_recv, blocs_to_recv);
  md_vector.copy(md);

  // Verify that the md_vector is valid (send sizes match receive sizes)
  if (comm_check_enabled())
    {
      IntVect toto;
      MD_Vector_tools::creer_tableau_distribue(md_vector, toto);
      toto = 0;
      toto.echange_espace_virtuel();
    }
}

/*! @brief This class provides tools to build the virtual space of an array containing indices of geometric entities
 *
 *   (vertices, elements, faces). It handles in particular the
 *   renumbering of virtual elements.
 *
 */
class Traduction_Indice_Global_Local
{
public:
  Traduction_Indice_Global_Local() {};
  void initialiser(const MD_Vector& md_items);
  void reset();
  void traduire_indice_local_vers_global(const ArrOfInt& indices_locaux, ArrOfTID& indices_globaux, int n) const;
  int traduire_indice_global_vers_local(const ArrOfTID& indices_globaux, ArrOfInt& indices_locaux) const;
  int traduire_espace_virtuel(IntTab& tableau) const;

  int chercher_table_inverse(const trustIdType sommet_global) const;

private:
  // Metadata of the indices to be renumbered:
  MD_Vector md_items_;
  trustIdType premier_indice_global_ = -100;
  // Distributed array (with virtual spaces and common items)
  // containing a global index for all entities to be indexed (real and virtual).
  // (if type_table_==SOMMETS, table_[i] is the global index of vertex i)
  TIDVect table_;
  // Table for inverting the numbering, sorted in ascending order
  // of the global index:
  // * column 0: the global index of the entity
  // * column 1: the local index of the entity
  TIDTab table_inverse_;
};

/*! @brief Initializes the dictionary. Precondition:
 *
 *   The remote spaces of the entities used must have been computed.
 *
 */
void Traduction_Indice_Global_Local::initialiser(const MD_Vector& md_items)
{
  md_items_ = md_items;

  // Build "table": create a global number for real entities
  // (entity index + total number of entities on lower-rank processors)
  // then exchange the virtual space of this array, obtaining for each
  // real or virtual entity its global number.

  table_.reset();
  MD_Vector_tools::creer_tableau_distribue(md_items, table_);

  const int nb_entites = md_items->get_nb_items_tot();
  const trustIdType decal = Process::mppartial_sum(nb_entites);
  premier_indice_global_ = decal;

  for (int i = 0; i < nb_entites; i++)
    table_[i] = i + decal;
  table_.echange_espace_virtuel();

  // Build table_inverse containing non-trivial indices
  // (for which table_[i] != i + decal after the exchange)
  // sorted in ascending order of the global number.
  const int nb_entites_tot = table_.size_totale();
  table_inverse_.resize(0, 2);

  for (int i = 0; i < nb_entites_tot; i++)
    {
      if (table_[i] != i + decal)
        table_inverse_.append_line(table_[i], i);
    }
  // insure complains.. check if it is right
  if (table_inverse_.size_array()>0)
    {
      tri_lexicographique_tableau(table_inverse_);
    }
}

void Traduction_Indice_Global_Local::reset()
{
  md_items_.detach();
  table_.reset();
  table_inverse_.reset();
}

/*! @brief Searches for i such that table_inverse(i, 0) == sommet_global, and returns table_inverse(i, 1) (the local index of the vertex).
 *
 *   If the vertex is not found in the table, returns -1.
 *   table_inverse must be sorted in ascending order of column 0.
 *   table_inverse must not have a virtual space.
 *
 */
int Traduction_Indice_Global_Local::chercher_table_inverse(const trustIdType sommet_global) const
{
  // Algorithm: binary search:
  int imin = 0;
  int imax = table_inverse_.dimension(0) - 1;
  // If only one element in the table, the while loop is not entered
  //        (so initialize to table_inverse(0, 0))
  // Otherwise, if no element, valeur must not equal sommet_global,
  //        otherwise, any value will do as it will be overwritten in the while loop.
  trustIdType valeur;
  if (imax == 0)
    valeur = table_inverse_(0, 0);
  else
    valeur = sommet_global - 1;

  while (imax > imin)
    {
      const int milieu = (imin + imax) >> 1; // (min+max)/2
      valeur = table_inverse_(milieu, 0);
      const trustIdType compare = valeur - sommet_global;
      if (compare < 0)
        imin = milieu + 1;
      else if (compare > 0)
        imax = milieu - 1;
      else
        imin = imax = milieu;
    }
  int resu = -1;
  valeur = table_inverse_(imin, 0);
  if (valeur == sommet_global)
    resu = static_cast<int>(table_inverse_(imin, 1)); // 2nd col always an int
  return resu;
}

/*! @brief Transforms local indices into global indices using "table_" (see initialiser).
 *
 * Does:
 *   For debut <= i < debut+nb
 *    indices_globaux[i] = table_[indices_locaux[i]]
 *    if indices_locaux[i] < 0 then indices_globaux[i] = -1
 *
 */
void Traduction_Indice_Global_Local::traduire_indice_local_vers_global(const ArrOfInt& indices_locaux,
                                                                       ArrOfTID& indices_globaux, int nb_items_a_traiter) const
{
  for (int i = 0; i < nb_items_a_traiter; i++)
    {
      const int i_loc = indices_locaux[i];
      const trustIdType i_glob = (i_loc < 0) ? -1 : table_[i_loc];
      indices_globaux[i] = i_glob;
    }
}

/*! @brief For debut <= i < debut+nb, indices_locaux[i] = look up the local index of "indices_globaux[i]"
 *
 * @param (indices_globaux) the array of global indices to translate
 * @param (indices_locaux) on output, the local indices or -1 if the global index was not found. Return value: number of indices not found (global indices that do not correspond to any local index).
 */
int Traduction_Indice_Global_Local::traduire_indice_global_vers_local(const ArrOfTID& indices_globaux,
                                                                      ArrOfInt& indices_locaux) const
{
  assert(indices_globaux.size_array() == indices_locaux.size_array());
  int i;
  int nb_erreurs = 0;
  const int nb_indices = indices_globaux.size_array();
  const int size_table = table_.size_array();
  for (i = 0; i < nb_indices; i++)
    {
      const trustIdType i_glob = indices_globaux[i];
      int i_loc;
      if (i_glob < 0)
        {
          // Negative index, considered normal,
          // it is an "empty index" marker.
          i_loc = -1;
        }
      else
        {
          // Check whether the item is not renumbered
          i_loc = static_cast<int>(i_glob - premier_indice_global_); // the diff is local, hence small
          if (i_loc < 0 || i_loc >= size_table || table_[i_loc] != i_glob)
            {
              // no, need to invert the table:
              i_loc = chercher_table_inverse(i_glob);
            }
          if (i_loc < 0)
            nb_erreurs++;
        }
      indices_locaux[i] = i_loc;
    }
  return nb_erreurs;
}

/*! @brief Starting from an array whose virtual space structure is initialized (remote and virtual element descriptors, common items)
 *
 *   and containing indices compatible with the content of the tables
 *   (vertex or element indices depending on type_table_),
 *   fills the virtual elements of "tableau" based on remote elements
 *   and translates the indices to local indices.
 *   (example, see construire_espace_virtuel_elements and
 *    construire_espace_virtuel_faces).
 *  Return value: number of indices that could not be translated
 *   (e.g. the referenced vertex does not exist on the neighboring processor)
 *
 */
int Traduction_Indice_Global_Local::traduire_espace_virtuel(IntTab& tab) const
{
  // Create a copy of the tab in which we will store global indices
  // Can not use 'copy' since value types are different (int vs TID), so this a bit clumsy:
  // (TODO provide 'from_int_to_tid' in TRUSTTab.h)
  TIDTab ind_glob_tab;
  ArrOfInt sz(tab.nb_dim());
  for (int i=0; i < tab.nb_dim(); i++) sz[i] = tab.dimension_tot(i);
  ind_glob_tab.resize(sz, RESIZE_OPTIONS::NOCOPY_NOINIT);
  ind_glob_tab.set_md_vector(tab.get_md_vector());

  IntVect& tableau = tab; // tab seen as a Vect.
  TIDVect& indices_globaux = ind_glob_tab;

  const int nb_items_reels    = tableau.size_reelle();
  const int nb_items_tot      = tableau.size_totale();
  const int nb_items_virtuels = nb_items_tot - nb_items_reels;

  // Translate real items to global indices:
  traduire_indice_local_vers_global(tableau, indices_globaux, nb_items_reels);
  // Fill virtual slots
  indices_globaux.echange_espace_virtuel();

  // Translate back only the virtual items of "tableau" to local indices:
  ArrOfTID src;
  ArrOfInt dest;
  src.ref_array(indices_globaux, nb_items_reels /*debut*/, nb_items_virtuels /*taille*/);
  dest.ref_array(tableau, nb_items_reels /*debut*/, nb_items_virtuels /*taille*/);
  const int nb_erreurs = traduire_indice_global_vers_local(src, dest);
  return nb_erreurs;
}

/*! @brief Builds the items_communs + virtual space structure of an array containing indices of geometric items, indexed by another geometric item type.
 *
 *   Example: array indexed by md_indice, containing indices of md_valeur:
 *      type_indice  type_valeur   example array:
 *       element      vertex       domaine.les_elems()
 *       face         vertex       faces_sommets
 *       element      face         elem_faces
 *       face         element      faces_voisins
 *       element      element      ?
 *       element      edge         elem_aretes
 *   Nb_valeurs_max is the number of real items of type "type_valeur".
 *
 */
void Scatter::construire_espace_virtuel_traduction(const MD_Vector& md_indice,
                                                   const MD_Vector& md_valeur,
                                                   IntTab& tableau,
                                                   const int error_is_fatal)
{
  if(Process::is_sequential())
    {
      // MD_Vector is a MD_Vector_seq:
      assert( dynamic_cast<const MD_Vector_seq *>(&md_indice.valeur()) != nullptr);
      // The array should still get its (dummy sequential) MD_Vector, otherwise it will remain null.
      if (!(tableau.get_md_vector() == md_indice))
        tableau.set_md_vector(md_indice);
      return;
    }

  if (tableau.dimension_tot(0) != md_indice->get_nb_items_reels()
      && (tableau.dimension_tot(0) != md_indice->get_nb_items_tot()))
    {
      Cerr << "[PE " << Process::me()
           << "] Error in Scatter::construire_espace_virtuel_traduction\n"
           << " the array does not have the good dimension on input" << finl;
      exit();
    }
  // Build the global/local index dictionary
  // for the values of the array
  Traduction_Indice_Global_Local dictionnaire_indices;
  dictionnaire_indices.initialiser(md_valeur);

  // Build the virtual space structure of "tableau"
  if (!(tableau.get_md_vector() == md_indice))
    MD_Vector_tools::creer_tableau_distribue(md_indice, tableau, RESIZE_OPTIONS::COPY_NOINIT);

  // Fill the virtual values of "tableau"
  const int nb_erreurs = dictionnaire_indices.traduire_espace_virtuel(tableau);

  if (nb_erreurs > 0 && error_is_fatal)
    {
      Cerr << "[PE " << Process::me()
           << "] Error in Scatter::construire_espace_virtuel_traduction\n"
           << " some indices of values were not found in\n"
           << " the local area : it missing virtual items"
           << finl;
      exit();
    }
}


/*! @brief Reorders joint faces so that they appear in the same order on each pair of neighboring processors.
 *
 * In practice, for a pair
 *   pe1 < pe2, pe1 sends its joint faces to pe2 and pe2 translates them to local
 *   vertex indices. The joint faces of PE2 are therefore not used.
 *
 */
void Scatter::reordonner_faces_de_joint(Domaine& dom)
{
  // Build the global/local index dictionary
  // for the vertices of the domain:
  Traduction_Indice_Global_Local dictionnaire_indices;
  dictionnaire_indices.initialiser(dom.les_sommets().get_md_vector());

  Schema_Comm schema_comm;
  Joints&      joints    = dom.faces_joint();
  const int nb_joints = joints.size();
  const int moi       = Process::me();
  int       i_joint;

  // Fill the recipient lists:
  //  send to higher-rank neighbors and receive from lower-rank neighbors.
  ArrOfInt send_list;
  ArrOfInt recv_list;



  for (i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      const int pe_voisin = joints[i_joint].PEvoisin();
      if (pe_voisin > moi)
        send_list.append_array(pe_voisin);
      else
        recv_list.append_array(pe_voisin);
    }

  schema_comm.set_send_recv_pe_list(send_list, recv_list);

  schema_comm.begin_comm();
  // Send joint faces translated to global vertex indices
  TIDTab faces_num_global;

  for (i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      const Joint&   joint         = joints[i_joint];
      const int   pe_voisin     = joint.PEvoisin();
      if (pe_voisin > moi)
        {
          const IntTab& faces_sommets = joint.faces().les_sommets();
          if (faces_sommets.dimension(0) > 0)
            {
              faces_num_global.resize(faces_sommets.dimension(0), faces_sommets.dimension(1));
              dictionnaire_indices.traduire_indice_local_vers_global(faces_sommets,
                                                                     faces_num_global,
                                                                     faces_sommets.size_array());
            }
          else
            faces_num_global.resize(0);

          schema_comm.send_buffer(pe_voisin) << faces_num_global;
        }
    }
  schema_comm.echange_taille_et_messages();
  // Receive joint faces and translate to local indices
  for (i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      Joint&       joint     = joints[i_joint];
      const int pe_voisin = joint.PEvoisin();
      if (pe_voisin < moi)
        {
          IntTab& faces_sommets = joint.faces().les_sommets();
          schema_comm.recv_buffer(pe_voisin) >> faces_num_global;
          if (faces_sommets.dimension(0) != faces_num_global.dimension(0))
            {
              Cerr << "[PE " << moi
                   << "] Error in Scatter::reordonner_faces_de_joint:\n"
                   << " the number of joint faces is not identical to the PE "
                   << pe_voisin << finl;
              exit();
            }
          const int nb_erreurs =
            dictionnaire_indices.traduire_indice_global_vers_local(faces_num_global,
                                                                   faces_sommets);
          if (nb_erreurs > 0)
            {
              Cerr << "[PE " << moi
                   << "] Error in Scatter::reordonner_faces_de_joint:\n"
                   << " The faces of the joint with PE " << pe_voisin
                   << " use of unknown nodes" << finl;
              exit();
            }
        }
    }
  schema_comm.end_comm();
}

/*! @brief Utility method: returns a complete list of all joint vertices (face vertices + isolated vertices), sorted and
 *
 *  without duplicates.
 *
 */
static void calculer_liste_complete_sommets_joint(const Joint& joint, ArrOfInt& liste_sommets)
{
  liste_sommets = joint.joint_item(JOINT_ITEM::SOMMET).items_communs();
#if 0

  // Take all vertices of the joint faces:
  const IntTab& som_faces = joint.faces().les_sommets();
  liste_sommets = ref_cast(ArrOfInt,som_faces);
  // Add all isolated vertices:
  const ArrOfInt& som_isoles = joint.sommets();
  const int n = som_isoles.size_array();
  for (int i = 0; i < n; i++)
    liste_sommets.append_array(som_isoles[i]);
  // Remove duplicates from the list
  array_trier_retirer_doublons(liste_sommets);
#endif
}

inline int arete_de_sommets_Si_et_Sj(const int Si, const int Sj, const int arete, const IntTab& aretes_som)
{
  if ( (aretes_som(arete,0) == Si && aretes_som(arete,1) == Sj)
       || (aretes_som(arete,1) == Si && aretes_som(arete,0) == Sj) )
    return 1;
  else
    return 0;
}

/*! @brief Utility method: returns a complete list of all joint edges (face edges + isolated edges), sorted and
 *
 *  without duplicates.
 *
 */
static void calculer_liste_complete_aretes_joint(const Joint& joint, ArrOfInt& liste_aretes)
{
  // Build the list of common edges liste_aretes

  ///////////////////////////////////////////////////////
  // Search for joint edges on joint faces
  ///////////////////////////////////////////////////////
  int nb_faces_joint=joint.faces().nb_faces();
  int nb_som_faces=joint.faces().nb_som_faces();
  const IntTab& sommet=joint.faces().les_sommets();
  const Domaine& dom=joint.domaine();
  const DoubleTab& coord=dom.coord_sommets();
  const IntTab& aretes_som=joint.domaine().aretes_som();
  ArrOfInt aretes(1);
  int compteur=0;
  DoubleTab positions(1,Objet_U::dimension);
  ArrOfInt som_faces(nb_faces_joint*nb_som_faces);
  // Traverse vertices of each joint face pairwise
  for (int face=0; face<nb_faces_joint; face++)
    for (int i=0; i<nb_som_faces; i++)
      {
        int Si = sommet(face,i);
        som_faces[face*nb_som_faces+i]=Si;
        for (int j=i; j<nb_som_faces; j++)
          {
            int Sj = sommet(face,j);
            // Compute midpoint C between vertices Si and Sj
            for (int comp=0; comp<Objet_U::dimension; comp++)
              positions(0,comp)=0.5*(coord(Si,comp)+coord(Sj,comp));
            dom.chercher_aretes(positions,aretes);
            // If an edge is found whose center coincides with point C
            // and whose vertices are identical to Si and Sj, add the edge to the list
            if (aretes[0]>=0 && arete_de_sommets_Si_et_Sj(Si, Sj, aretes[0], aretes_som))
              {
                compteur++;
                liste_aretes.append_array(aretes[0]);
              }
          }
      }
  Process::Journal() << "common edges found on faces of joint with " << joint.PEvoisin() << " :" << compteur << finl;
  /////////////////////////////////////////////////////////////////////////
  // Search for isolated joint edges on isolated joint vertices
  /////////////////////////////////////////////////////////////////////////
  // joint.sommets() sometimes contains all vertices!
  // So we build a som_isoles array
  ArrOfInt som_isoles;
  // Put all vertices in som_isoles (isolated + from joint faces):
  calculer_liste_complete_sommets_joint(joint, som_isoles);
  // Sort som_faces and remove duplicates
  array_trier_retirer_doublons(som_faces);
  // Remove all vertices of som_isoles that are in som_faces
  array_retirer_elements(som_isoles, som_faces);
  // Remove joint face vertices
  const int n = som_isoles.size_array();
  Process::Journal() << "number of isolated nodes: " << n << finl;
  Process::Journal() << "number of nodes of faces of joint: " << 3*sommet.dimension(0) << finl;

  compteur=0;
  // Traverse isolated vertices pairwise
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      {
        // Compute midpoint C between vertices Si and Sj
        int Si = som_isoles[i];
        int Sj = som_isoles[j];
        for (int comp=0; comp<Objet_U::dimension; comp++)
          positions(0,comp)=0.5*(coord(Si,comp)+coord(Sj,comp));
        dom.chercher_aretes(positions,aretes);
        // If an edge is found whose center coincides with point C
        // and whose vertices are identical to Si and Sj, add the edge to the list
        if (aretes[0]>=0 && arete_de_sommets_Si_et_Sj(Si, Sj, aretes[0], aretes_som))
          {
            compteur++;
            liste_aretes.append_array(aretes[0]);
          }
      }
  Process::Journal() << "common edges found isolated on joint with " << joint.PEvoisin() << " :" << compteur << finl;
  // Remove duplicates from the list
  array_trier_retirer_doublons(liste_aretes);
}

static void calculer_liste_complete_items_joint(const Joint& joint, const JOINT_ITEM type_item, ArrOfInt& liste_items)
{
  switch(type_item)
    {
    case JOINT_ITEM::SOMMET:
      calculer_liste_complete_sommets_joint(joint, liste_items);
      break;
    case JOINT_ITEM::ARETE:
      calculer_liste_complete_aretes_joint(joint, liste_items);
      break;
    default:
      Cerr << "Error in Scatter::calculer_liste_complete_items_joint" << finl;
      Cerr << "Type of item not expected." << finl;
      Process::exit();
    }
}

/*! @brief Current periodic algorithms (P1B assembler, OpDivElem P1B) require that for each virtual periodic face, the opposite face is
 *
 *   also virtual. This is not guaranteed at the output of the
 *   calculer_elements_distants method. This method adds the missing elements
 *   to the remote spaces to ensure this condition:
 *   If a remote element for a given PE is adjacent to a periodic face,
 *   the element adjacent to the opposite face is added to the remote space.
 *
 */
void Scatter::corriger_espace_distant_elements_perio(Domaine& dom)
{
  if (Process::je_suis_maitre())
    Cerr << "Correction of remote spaces of the elements for the periodic faces" << finl;

  const Noms& liste_bords_periodiques = dom.bords_perio();

  const int nb_elem = dom.nb_elem();
  const IntTab& les_elems = dom.les_elems();

  // This array will contain, for a given periodic boundary:
  //  if element i is adjacent to a face of this boundary,
  //  element_oppose[i] is the number of the element adjacent to the
  //  opposite face on this boundary.
  // -1 otherwise.
  ArrOfInt element_oppose(nb_elem);

  Static_Int_Lists connectivite_som_elem;
  const int nb_sommets = dom.nb_som();
  construire_connectivite_som_elem(nb_sommets,
                                   les_elems,
                                   connectivite_som_elem,
                                   0 /* do not include virtual vertices */);

  const int nb_som_face = dom.type_elem()->nb_som_face();
  ArrOfInt une_face(nb_som_face);
  ArrOfInt elems_voisins;


  const int nb_joints = dom.nb_joints();

  // Markers for existing remote elements:
  ArrOfBit marqueurs_elements_distants(nb_elem);
  marqueurs_elements_distants = 0;

  // Adding an element to a remote space for a given boundary
  // may cause a problem on another boundary (the opposite element of the
  // newly added element may be missing for another periodicity direction).
  // We must therefore iterate until nothing changes.
  int nb_elements_ajoutes = 0;
  do
    {
      nb_elements_ajoutes = 0;
      for (auto& itr : liste_bords_periodiques)
        {
          const Nom& nom_bord = itr;
          const Bord& bord = dom.bord(nom_bord);
          const IntTab& faces_sommets = bord.les_sommets_des_faces();
          const int nb_faces = bord.nb_faces();

          // First step: identify the opposite elements for this periodic boundary.
          // Loop over the first half of the boundary.
          // WARNING: assumes that boundary faces are ordered: first all faces
          // on one side of the periodic domain, then, in the same order, the opposite faces.
          element_oppose = -1;

          for (int i_face = 0; i_face < nb_faces / 2; i_face++)
            {
              // For each face, find the element adjacent to this face and to the opposite face.
              // Loop over the face and the opposite face.
              int elem0 = -1; // The two opposite elements of this periodic face
              int elem1 = -1;
              for (int quel_cote = 0; quel_cote < 2; quel_cote++)
                {
                  const int face = i_face + quel_cote * nb_faces / 2;
                  int i;
                  for (i = 0; i < nb_som_face; i++)
                    une_face[i] = faces_sommets(face, i);
                  find_adjacent_elements(connectivite_som_elem, une_face, elems_voisins);
                  const int n = elems_voisins.size_array();
                  if (n != 1)
                    {
                      Cerr << "Error in Scatter::corriger_espace_distant_elements_perio: \n"
                           << " The face " << i_face << " of boundary " << nom_bord << " has "
                           << n << " neighbors." << finl;
                      Process::exit();
                    }
                  if (quel_cote == 0)
                    elem0 = elems_voisins[0];
                  else
                    elem1 = elems_voisins[0];
                }
              element_oppose[elem0] = elem1;
              element_oppose[elem1] = elem0;
            }
          // Second step: iterate over remote elements. If a remote element is
          // among the paired elements, add the other pair member to the remote elements.
          for (int i_joint = 0; i_joint < nb_joints; i_joint++)
            {
              ArrOfInt& elements_distants = dom.joint(i_joint).set_joint_item(JOINT_ITEM::ELEMENT).set_items_distants();
              int n = elements_distants.size_array();
              // Mark existing remote elements:
              int i;
              for (i = 0; i < n; i++)
                {
                  const int elem = elements_distants[i];
                  marqueurs_elements_distants.setbit(elem);
                }

              for (i = 0; i < n; i++)
                {
                  const int elem = elements_distants[i];
                  const int elem_oppose = element_oppose[elem];
                  if (elem_oppose >= 0 && (!marqueurs_elements_distants.testsetbit(elem_oppose)))
                    {
                      elements_distants.append_array(elem_oppose);
                      nb_elements_ajoutes++;
                    }
                }
              // Reset the markers array to zero:
              n = elements_distants.size_array();
              for (i = 0; i < n; i++)
                {
                  const int elem = elements_distants[i];
                  marqueurs_elements_distants.clearbit(elem);
                }
            }
        }

    }
  while (nb_elements_ajoutes > 0);
  // Final sort of remote elements in ascending order
  for (int i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      ArrOfInt& elements_distants = dom.joint(i_joint).set_joint_item(JOINT_ITEM::ELEMENT).set_items_distants();

      elements_distants.ordonne_array();
    }
}

/*! @brief Fills the "espace_distant()" array of elements in the joints.
 *
 * This is where joint elements are determined based on joint thickness.
 *  The espace_distant array contains the
 *  local indices of remote elements (to be sent to neighboring processors).
 *  For a joint of thickness 1, these are all elements neighboring a
 *  joint vertex (vertex on a joint face or isolated vertex).
 *  For a joint of thickness n>1, these are all elements neighboring a
 *  vertex of an element from the joint of thickness n-1.
 *  Neighborhood is defined on the global domain (all subdomains combined).
 *  History: first version by B.Mathieu on 16/01/2007.
 *    A method that determines remote elements at partition time also exists
 *    (DomaineCutter::construire_elements_distants_ssdom).
 *    The method below has been validated by comparison with the partitioner method.
 *    Outputs have been verified for thicknesses up to 5 on tetrahedral meshes.
 *  The difficulty of the algorithm is to obtain virtual elements of thickness > 1 that
 *  lie on subdomains not in direct contact with the local subdomain.
 *  This difficulty is resolved by the algorithm below.
 *
 */
void Scatter::calculer_espace_distant_elements(Domaine& dom)
{
  const int nbjoints    = dom.nb_joints();
  const int nb_som_elem  = dom.nb_som_elem();
  const IntTab& les_elems    = dom.les_elems();
  const int nproc        = Process::nproc();
  // PL: all processors must have the same epaisseur_joint for the following algorithm
  // that uses data exchange with schema_comm.begin_comm() schema_comm.end_comm()
  // otherwise there is a deadlock in debug mode due to the exit at line 1866 (see case Quasi_Comp_Coupl_Incomp)
  //const int epaisseur_joint = (nb_joints > 0) ? domaine.joint(0).epaisseur() : 1;
  const int epaisseur_joint = (int) mp_max((nbjoints > 0) ? dom.joint(0).epaisseur() : 1);

  if (Process::je_suis_maitre())
    Cerr << "Calculation of remote space of elements : thickness " << epaisseur_joint << finl;

  // The algorithm is based on the progressive construction of the "liste_sommets" array.
  // liste_sommets(pe) contains at a given moment the list of vertices owned by me()
  // for which "pe" wants to know the neighboring elements.
  // Each PE starts by requesting its direct neighboring elements, i.e. the elements
  // neighboring joint vertices. Knowing that it will subsequently want the elements
  // neighboring those found, we add to the lists the vertices of the remote elements
  // found in the previous iteration.
  ArrsOfInt liste_sommets(nproc);

  // For each processor, the list of local elements to send to it
  ArrsOfInt elements_distants(nproc);
  {
    // smart_resize because we will call append_array on these arrays
    for (int i = 0; i < nproc; i++)
      {


      }
  }

  Static_Int_Lists som_elem;
  {
    if (Process::je_suis_maitre())
      Cerr << "Nodes-elements connectivity ..." << finl;
    const int nb_sommets = dom.nb_som();
    construire_connectivite_som_elem(nb_sommets,
                                     les_elems,
                                     som_elem,
                                     0 /* do not include virtual vertices */);
  }

  ArrOfInt liste_pe_voisins(nbjoints);

  // Initialize liste_pe_voisins and
  // initialize the vertex list: each processor requests all elements
  // neighboring its joint vertices. The local processor knows that for thickness 1,
  // each neighboring processor via a joint wants to know all elements neighboring
  // the joint vertices. So we put in liste_sommets(pe) the joint vertices with that PE,
  // so as to send it the local elements neighboring those vertices.
  {
    for (int i_joint = 0; i_joint < nbjoints; i_joint++)
      {
        const Joint& joint = dom.joint(i_joint);
        const int pe = joint.PEvoisin();
        liste_pe_voisins[i_joint] = pe;
        const ArrOfInt& sommets_joint = joint.joint_item(JOINT_ITEM::SOMMET).items_communs();
        liste_sommets[pe] = sommets_joint;
      }
  }

  // We will need a communication scheme where each processor sends and receives
  // data to/from its direct neighbors (neighbors sharing a vertex)
  Schema_Comm schema_comm;
  schema_comm.set_send_recv_pe_list(liste_pe_voisins, liste_pe_voisins);

  // We will need fast access to the list of processors that share a
  // vertex and the index of that vertex on each processor.
  // Structure content:
  //  data_sommets_communs.get_list_size(sommet) = 2 * number of procs sharing the vertex
  //  data_sommets_communs(sommet, 2*i) = neighboring PE number
  //  data_sommets_communs(sommet, 2*i+1) = index of the vertex on that PE.
  Static_Int_Lists data_sommets_communs;
  // Filling: the structure is only used if joint thickness is > 1
  if (epaisseur_joint > 1)
    {
      ArrOfInt count(dom.nb_som());
      // Step 1: with how many processors is each vertex shared?
      for (int ijoint = 0; ijoint < nbjoints; ijoint++)
        {
          const Joint& joint = dom.joint(ijoint);
          const ArrOfInt& sommets_joint = joint.joint_item(JOINT_ITEM::SOMMET).items_communs();
          const int n = sommets_joint.size_array();
          for (int i = 0; i < n; i++)
            {
              const int som = sommets_joint[i];
              count[som] += 2; // 2 integers stored per common vertex
            }
        }
      data_sommets_communs.set_list_sizes(count);
      count = 0;
      // Step 2: fill the structure:
      for (int ijoint = 0; ijoint < nbjoints; ijoint++)
        {
          const Joint& joint = dom.joint(ijoint);
          const int pe = joint.PEvoisin();
          const IntTab& renum_sommets = joint.joint_item(JOINT_ITEM::SOMMET).renum_items_communs();
          const int n = renum_sommets.dimension(0);
          for (int i = 0; i < n; i++)
            {
              // Index of the shared vertex on the neighboring PE
              const int i_sommet_distant = renum_sommets(i, 0);
              // Index of the shared vertex on my local domain
              const int i_sommet_local = renum_sommets(i, 1);
              const int j = count[i_sommet_local]++;
              data_sommets_communs.set_value(i_sommet_local, j*2, pe);
              data_sommets_communs.set_value(i_sommet_local, j*2+1, i_sommet_distant);
            }
        }
    }

  // Loop over joint thickness:
  // At the start of the loop, liste_sommets is assumed to contain, for each processor
  // requesting virtual elements, the list of vertices of me() whose neighbors it wants.
  for (int epaisseur = 1; ; epaisseur++)
    {
      if (Process::je_suis_maitre())
        Cerr << " Calculation of the thickness " << epaisseur << finl;

      // For each vertex list, put in the remote elements of the same processor
      // the elements neighboring the vertices in the list.
      int pe;
      for (pe = 0; pe < nproc; pe++)
        {
          ArrOfInt& elems_dist = elements_distants[pe];
          const ArrOfInt& sommets = liste_sommets[pe];
          elems_dist.resize_array(0);
          const int nb_som_liste = sommets.size_array();
          for (int isom = 0; isom < nb_som_liste; isom++)
            {
              const int som = sommets[isom];
              if (som<0)
                continue;
              const int nb_elem_som = som_elem.get_list_size(som);
              for (int ielem = 0; ielem < nb_elem_som; ielem++)
                {
                  const int elem = som_elem(som, ielem);
                  elems_dist.append_array(elem);
                }
            }
          array_trier_retirer_doublons(elems_dist);
        }

      // The following updates liste_sommets for the next iteration.
      // No need to do this if we are at the last iteration:
      if (epaisseur == epaisseur_joint)
        break;

      // Put in the vertex lists the vertices of the remote elements found
      for (pe = 0; pe < nproc; pe++)
        {
          ArrOfInt& sommets = liste_sommets[pe];
          const ArrOfInt& elems_dist = elements_distants[pe];
          sommets.resize_array(0);
          const int nb_elems_dist = elems_dist.size_array();
          for (int ielem = 0; ielem < nb_elems_dist; ielem++)
            {
              const int elem = elems_dist[ielem];
              for (int isom = 0; isom < nb_som_elem; isom++)
                {
                  const int som = les_elems(elem, isom);
                  sommets.append_array(som);
                }
            }
          array_trier_retirer_doublons(sommets);
        }
      // Traverse the vertex lists. For each vertex, if it is a joint vertex,
      // send to the processors owning that vertex a request "processor i
      // wants all neighbors of this vertex".
      // Do not send the request to processor "i" for list "i": it already knows
      // its own elements!
      schema_comm.begin_comm();
      // First communication phase: stack data to send in buffers
      for (pe = 0; pe < nproc; pe++)
        {
          const ArrOfInt& sommets = liste_sommets[pe];
          const int nb_som_liste = sommets.size_array();
          for (int isom = 0; isom < nb_som_liste; isom++)
            {
              const int i_sommet_local = sommets[isom];
              if (i_sommet_local<0)
                continue;
              const int nb_pe_voisins = data_sommets_communs.get_list_size(i_sommet_local) / 2;
              for (int i = 0; i < nb_pe_voisins; i++)
                {
                  const int pe_voisin = data_sommets_communs(i_sommet_local, i*2);
                  // Index of the vertex on the neighboring processor.
                  const int i_sommet_distant = data_sommets_communs(i_sommet_local, i*2+1);
                  if (pe_voisin != pe)
                    {
                      // Send to processor "pe_voisin" the message: "processor PE needs
                      // the elements neighboring vertex i_sommet_distant"
                      schema_comm.send_buffer(pe_voisin) << pe << i_sommet_distant;
                    }
                }
            }
        }
      schema_comm.echange_taille_et_messages();
      for (int i_pevoisin = 0; i_pevoisin < nbjoints; i_pevoisin++)
        {
          const int pe_voisin = liste_pe_voisins[i_pevoisin];
          Entree& buffer = schema_comm.recv_buffer(pe_voisin);
          for (;;)
            {
              int pe2, sommet;
              // Retrieve the message "processor PE needs the elements neighboring SOMMET".
              buffer >> pe2 >> sommet;
              if (buffer.eof())
                break;
              liste_sommets[pe2].append_array(sommet);
            }
        }
      schema_comm.end_comm();
      // Remove duplicates from the vertex lists
      for (pe = 0; pe < nproc; pe++)
        {
          ArrOfInt& sommets = liste_sommets[pe];
          array_trier_retirer_doublons(sommets);
        }
    }

  // Create new joints if needed, and store remote elements in the joints
  {
    ArrOfInt voisins;


    for (int pe = 0; pe < nproc; pe++)
      if (elements_distants[pe].size_array() > 0)
        voisins.append_array(pe);
    ajouter_joints(dom, voisins);

#ifdef CHECK_ALGO_ESPACE_VIRTUEL
    // We do not use the virtual spaces computed above; we only compare them
    // to the virtual spaces computed at partition time by the sequential algorithm.
    // For now, the parallel algorithm appears to work correctly;
    // this test is disabled. (Benoit Mathieu)
    bool erreur = false;
    const int nbjoints = dom.nbjoints();
    for (int i = 0; i < nbjoints; i++)
      {
        Joint& joint = dom.joint(i);
        const int pe = joint.PEvoisin();
        if (!(joint.joint_item(JOINT_ITEM::ELEMENT).items_distants() == elements_distants[pe]))
          {
            Cerr << "Error in Scatter, PE " << Process::me() << finl;
            Process::Journal() << "Error scatter, remote elements pe " << pe << finl
                               << " Splitting algorithm: " << joint.joint_item(JOINT_ITEM::ELEMENT).items_distants()
                               << " Scatter algorithm  : " << elements_distants[pe] << finl;

            erreur = true;
          }
      }
    if (mp_or(erreur))
      Process::exit();
#else
    // Store the result
    const int nb_joints = dom.nb_joints();
    for (int i = 0; i < nb_joints; i++)
      {
        Joint& joint = dom.joint(i);
        const int pe = joint.PEvoisin();
        joint.set_joint_item(JOINT_ITEM::ELEMENT).set_items_distants() = elements_distants[pe];
      }
#endif
  }
}

static inline int fct_cmp_coordonnees(const double * s1, const double *s2, int dim, const double epsilon)
{
  assert(dim==2 || dim==3);
  if (s1[0] < s2[0] - epsilon)
    return -1;
  else if (s1[0] > s2[0] + epsilon)
    return 1;
  else if (s1[1] < s2[1] - epsilon)
    return -1;
  else if (s1[1] > s2[1] + epsilon)
    return 1;
  else if (dim<3)
    return 0;
  else if (s1[2] < s2[2] - epsilon)
    return -1;
  else if (s1[2] > s2[2] + epsilon)
    return 1;
  else
    return 0;
}

/*! @brief Builds the "correspondance" array such that for 0 <= i < sommets2.
 *
 * size_array(),
 *    If sommet2(i) exists in sommets1, then
 *       sommets2(i, ...) == sommets1(correspondance[i], ...)
 *    Otherwise
 *       correspondance[i] = -1
 *   Equality is checked to within epsilon in absolute value (i.e. abs(x1-x2) < epsilon)
 *   The algorithm is generally O(n1*log(n1) + n2*log(n1))
 *   (search based on quicksort).
 *   If the sort fails, an O(n1*n2) algorithm is used.
 *   Arrays sommets1 and sommets2 must be 2-dimensional.
 *   The correspondance array must have size sommets2.size_array().
 *  Return value: number of vertices of sommets2 not found in sommets1.
 *
 */
int Scatter::Chercher_Correspondance(const DoubleTab& sommets1, const DoubleTab& sommets2,
                                     ArrOfInt& correspondance, const double epsilon)
{
  const int nb_sommets1 = sommets1.dimension(0);
  const int nb_sommets2 = sommets2.dimension(0);
  // Precondition required for fct_cmp_index_coord
  assert(sommets1.nb_dim() == 2);
  assert(sommets2.nb_dim() == 2);
  assert(correspondance.size_array() == nb_sommets2);
  if (nb_sommets1 < 1)
    {
      correspondance = -1;
      return nb_sommets2;
    }

  // Sorted indirection array such that coordinates sommets1(index[i], .) are
  // sorted in lexicographic order.
  ArrOfInt index(nb_sommets1);
  {
    int i;
    for (i = 0; i < nb_sommets1; i++)
      index[i] = i;
  }

  // Sort the index array
  // Vertices are sorted in lexicographic coordinate order.
  // Since we test with a tolerance of epsilon, the sort may fail:
  //  If x=1, y=1.01, z=1.02 and epsilon=0.01, we have
  //   x==y (to within epsilon)
  //   y==z (to within epsilon)
  //  but x!=z
  // So the binary search may fail subsequently.
  tri_lexicographique_tableau_indirect(sommets1, index);

  // Build the correspondance array such that
  //   sommet1(correspondance[i], ...) == sommet2(i, ...)
  int nb_sommets_non_trouves = 0;
  int nb_echec_dichotomie = 0;
  {
    int i;
    int nb_dim = sommets1.dimension(1);
    for (i = 0; i < nb_sommets2; i++)
      {
        const double * s2 = & sommets2(i,0);
        int num_sommet = -1;

        // First search for the vertex in sommets1 using binary search (bsearch)
        int imin = 0;
        int imax = nb_sommets1 - 1;
        int resu_cmp = -1;
        int k = -1;
        while (imax > imin)
          {
            const int milieu = (imin + imax) >> 1; // (min+max)/2
            k = index[milieu];
            const double * s1 = & sommets1(k, 0);
            resu_cmp = fct_cmp_coordonnees(s1, s2, nb_dim, epsilon);
            switch(resu_cmp)
              {
              case -1:
                imin = milieu + 1;
                break; // s1<s2
              case 1 :
                imax = milieu - 1;
                break; // s1>s2
              default:
                imin = imax = milieu;
                break; // s1==s2
              }
          }
        if (resu_cmp != 0)
          {
            k = index[imin];
            const double * s1 = & sommets1(k, 0);
            resu_cmp = fct_cmp_coordonnees(s1, s2, nb_dim, epsilon);
          }
        if (resu_cmp == 0)
          {
            num_sommet = k;
          }
        else
          {
            nb_echec_dichotomie++;
            // If failure, the array may not be correctly ordered
            // => search for the vertex by scanning the entire array
            int j;
            for (j = 0; j < nb_sommets1; j++)
              {
                const double * s1 = & sommets1(j,0);
                resu_cmp = fct_cmp_coordonnees(s1, s2, nb_dim, epsilon);
                if (resu_cmp == 0)
                  break;
              }
            if (j < nb_sommets1)
              num_sommet = j;
            else
              nb_sommets_non_trouves++;
          }

        correspondance[i] = num_sommet;
      }
  }

  if (nb_echec_dichotomie > 0)
    Process::Journal() << "Chercher_Correspondance Dichotomy failure rate "
                       << nb_echec_dichotomie << " / " << nb_sommets2 << finl;

  return nb_sommets_non_trouves;
}

/*! @brief Generic method to build geometrical item correspondance between the local and the remote processor
 * around a joint.
 *
 * See also construire_correspondance_sommets_par_coordonnees() for the very specific usage of allow_resize.
 *
 */
void Scatter::construire_correspondance_items_par_coordonnees(Joints& joints, const JOINT_ITEM type_item,
                                                              const DoubleTab& coord_items, bool allow_resize)
{
  switch(type_item)
    {
    case JOINT_ITEM::SOMMET:
      break;
    case JOINT_ITEM::ARETE:
      break;
    default:
      Cerr << "Scatter::construire_correspondance_items_par_coordonnees unusable for item "
           << (int)type_item
           << finl;
      Process::exit();
    }
  const int dim = Objet_U::dimension;
  const int nb_joints = joints.size();

  // Indices of joint items in the domain on my processor
  ArrsOfInt  indices_items_locaux(nb_joints);
  // Indices of joint items in the domain on the neighboring processor
  ArrsOfInt  indices_items_distants(nb_joints);
  // Coordinates of the corresponding items (in the same order as indices_items_xxx)
  DoubleTabs coord_items_locaux(nb_joints);
  DoubleTabs coord_items_distants(nb_joints);

  // Fill indices_items_locaux
  // and coord_items_locaux
  for (int i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      const Joint& joint = joints[i_joint];
      ArrOfInt& items = indices_items_locaux[i_joint];
      // Note: **SORTED_LIST** indices_items_locaux are
      //  sorted in ascending order:
      calculer_liste_complete_items_joint(joint, type_item, items);

      const int n       = items.size_array();
      DoubleTab&   coord   = coord_items_locaux[i_joint];
      coord.resize(n, dim);
      for (int i = 0; i < n; i++)
        for (int j = 0; j < dim; j++)
          coord(i,j) = coord_items(items[i], j);
    }

  // Send local indices and coordinates to the neighboring processor
  {
    Schema_Comm schema_comm;
    ArrOfInt liste_pe_voisins(nb_joints);
    int i;
    for (i = 0; i < nb_joints; i++)
      liste_pe_voisins[i] = joints[i].PEvoisin();
    schema_comm.set_send_recv_pe_list(liste_pe_voisins, liste_pe_voisins);
    schema_comm.begin_comm();
    for (i = 0; i < nb_joints; i++)
      {
        const int pe = liste_pe_voisins[i];
        Sortie& buffer = schema_comm.send_buffer(pe);
        buffer << indices_items_locaux[i];
        buffer << coord_items_locaux[i];
      }
    schema_comm.echange_taille_et_messages();
    for (i = 0; i < nb_joints; i++)
      {
        const int pe = liste_pe_voisins[i];
        Entree& buffer = schema_comm.recv_buffer(pe);
        buffer >> indices_items_distants[i];
        buffer >> coord_items_distants[i];
      }
    schema_comm.end_comm();
  }

  // Loop over joints
  // This time, joints are modified (filling renum_virt_loc)
  const int moi = Process::me();
  for (int i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      Joint&           joint           = joints[i_joint];
      const int     PEvoisin        = joint.PEvoisin();
      const ArrOfInt& indices_locaux  = indices_items_locaux[i_joint];
      const DoubleTab& coord_locaux    = coord_items_locaux[i_joint];
      const DoubleTab& coord_distants  = coord_items_distants[i_joint];
      const int n = indices_locaux.size_array();

      // Search for correspondences between items
      ArrOfInt corresp(n);
      const double epsilon = Objet_U::precision_geom;
      Chercher_Correspondance(coord_distants, coord_locaux, corresp, epsilon);

      int nb_items_communs_trouves=0;
      for (int k = 0; k < n; k++)
        if (corresp[k]>=0) nb_items_communs_trouves++;

      ArrOfInt& items_communs = joint.set_joint_item(type_item).set_items_communs();
      if (allow_resize)
        items_communs.resize_array(nb_items_communs_trouves);
      else
        // If a resisze is not expected, everthing in items_communs should have been found:
        assert(items_communs.size_array() == nb_items_communs_trouves);
      items_communs = -1;

      // If a resize is expected, we need to shift the indices in 'corresp' to fit into
      // a (smaller) array of size 'nb_items_communs_trouves', avoiding the non-matching indices
      int n_dist = coord_distants.dimension(0);
      std::vector<bool> corres_ok(n_dist, false);
      std::vector<int> offset(n_dist, 0);
      if(allow_resize)
        {
          // corres_ok[j] is true iif the (remote) item 'j' was matched with a local one
          for(int k = 0; k < n; k++)
            if (corresp[k] >= 0)
              corres_ok[corresp[k]] = true;
          // offset[k] is the shift to be substracted to the remote item index once the invalid (=non
          // matched) remote items have been removed:
          int nb_holes = 0;  // nb of holes seen so far in corres_ok
          for(int k = 0; k < n_dist; k++)
            offset[k] = corres_ok[k] ? nb_holes : nb_holes++;
        }

      int i=0;
      for (int k = 0; k < n; k++)
        {
          const int i_local = indices_locaux[k];
          // The j-th remote item is identical to the k-th local item
          const int j = corresp[k];

          // Not found? Possible error
          if (j < 0)
            {
              if (!allow_resize)
                {
                  Cerr << "Error in Scatter::remplir_renum_virt_loc on PE " << moi << finl
                       << "The item of type " << (int)type_item << " number " << i_local << " with coordinates ";
                  for (int k2 = 0; k2 < dim; k2++)
                    Cerr << coord_locaux(i, k2) << " ";
                  Cerr << finl << "was not found in the joint with the PE " << PEvoisin << finl;
                  if (type_item==JOINT_ITEM::ARETE)
                    {
                      Cerr << "The searching algorithm of the isolated edges on a joint" << finl;
                      Cerr << "does not work yet in some cases. Two isolated nodes of a joint (example below" << finl;
                      Cerr << "joint between 0 and 2) can be those of an edge not belonging to this joint (below" << finl;
                      Cerr << "the edge belongs to the joint 0-1 but not 0-2):" << finl;
                      //Cerr << "  ________    " << finl;
                      //Cerr << "1\ 2/1\2 /1\  " << finl;
                      //Cerr << "__\/___\/___\ " << finl;
                      //Cerr << " 0/\ 0 /\ 0   " << finl;
                      //Cerr << " / 0\ / 0\    " << finl;
                      Cerr << "  ________      " << finl;
                      Cerr << "1\\ 2/1\\2 /1\\ " << finl;
                      Cerr << "__\\/___\\/___\\" << finl;
                      Cerr << " 0/\\ 0 /\\ 0   " << finl;
                      Cerr << " / 0\\ / 0\\    " << finl;
                      Cerr << finl;
                      Cerr << "One way to by-pass this problem is to split again your domain with" << finl;
                      Cerr << "different options of splitting or with another splitter to do not fall" << finl;
                      Cerr << "on the same configuration." << finl;
                    }
                  exit();
                }
            }
          else
            {
              if (moi < PEvoisin)
                {
                  // common items in the order of local joint items:
                  items_communs[i] = i_local;
                  // Verify that it is indeed in ascending order of the local index
                  // (see **SORTED_LIST**)
                  assert(i==0 || items_communs[i] > items_communs[i-1]);
                }
              else
                {
                  int j2 = j - offset[j];
                  // common items in the order of items on the neighbor:
                  assert(items_communs[j2] < 0);
                  items_communs[j2] = i_local;
                }
              i++;
            }
        }
      assert(i==nb_items_communs_trouves);
    }
  // Fill renum_items_communs:
  calculer_renum_items_communs(joints, type_item);
}

/*! @brief Builds the joint_item(JOINT_ITEM::SOMMET).items_communs arrays for all joints of the domain dom.
 *
 * @param dom The domain to process.
 * @param allow_resize may be set to True in some rare case (see Raffiner_isotrope_parallele)
 * when we know that the current size of 'items_communs' is wrong because part of the domain
 * was resized / changed.
 */
void Scatter::construire_correspondance_sommets_par_coordonnees(Domaine& dom, bool allow_resize)
{
  construire_correspondance_items_par_coordonnees(dom.faces_joint(), JOINT_ITEM::SOMMET, dom.coord_sommets(), allow_resize);
}

/*! @brief Builds the joint_item(JOINT_ITEM::ARETE).items_communs arrays for all joints of the domain.
 *
 * @param zvf The VF domain to process.
 */
void Scatter::construire_correspondance_aretes_par_coordonnees(Domaine_VF& zvf)
{
  construire_correspondance_items_par_coordonnees(zvf.domaine().faces_joint(), JOINT_ITEM::ARETE, zvf.xa());
}

/*! @brief For a geometric item "type_item", fills the nb_items_virtuels_ field of joints based on
 *  the number of remote items:
 *
 *   The number of virtual items on joint i of processor j is the
 *   number of remote items of joint j on processor i.
 *
 */
void Scatter::calculer_nb_items_virtuels(Joints& joints,
                                         const JOINT_ITEM type_item)
{
  Schema_Comm schema_comm;
  const int nb_joints = joints.size();
  ArrOfInt liste_voisins(nb_joints);
  int i_joint;
  for (i_joint = 0; i_joint < nb_joints; i_joint++)
    liste_voisins[i_joint] = joints[i_joint].PEvoisin();

  // Send the number of remote items to the neighboring PE
  schema_comm.set_send_recv_pe_list(liste_voisins, liste_voisins);
  schema_comm.begin_comm();
  for (i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      const Joint&        joint = joints[i_joint];
      const Joint_Items& items = joint.joint_item(type_item);
      const int n  = items.items_distants().size_array();
      const int pe = joint.PEvoisin();
      schema_comm.send_buffer(pe) << n;
    }
  // Exchange messages
  schema_comm.echange_taille_et_messages();
  // The neighboring PE receives this number of items and stores it.
  for (i_joint = 0; i_joint < nb_joints; i_joint++)
    {
      Joint&        joint = joints[i_joint];
      Joint_Items& items = joint.set_joint_item(type_item);
      const int pe = joint.PEvoisin();
      int n;
      schema_comm.recv_buffer(pe) >> n;
      items.set_nb_items_virtuels(n);
    }
  schema_comm.end_comm();
}

/*! @brief Create parallel descriptors for the vertex and element arrays of the domain (necessary because Scatter is
 *  never invoked in sequential).
 *
 *  In 64bit the corresponding number of items might be big. This is here the main justification for the need of the class
 *  MD_Vector_seq which unique useful argument is the total number of items (with type trustIdType).
 *  Alternative would have been to make all members of MD_Vector_std compatible with trustIdType ...
 */
template <typename _SIZE_>
void Scatter::init_sequential_domain(Domaine_32_64<_SIZE_>& dom)
{
  MD_Vector_seq mdseq_som(dom.les_sommets().dimension(0));
  MD_Vector md;
  md.copy(mdseq_som);
  dom.les_sommets().set_md_vector(md);
  MD_Vector_seq mdseq_elem(dom.les_elems().dimension(0));
  md.copy(mdseq_elem);
  dom.les_elems().set_md_vector(md);
}

/*! @brief Method used by interpreters that modify the domain (sequential), destroys the descriptors
 *  of vertices and elements to allow modification of these arrays.
 */
template <typename _SIZE_>
void Scatter::uninit_sequential_domain(Domaine_32_64<_SIZE_>& dom)
{
  MD_Vector md; // null descriptor
  dom.les_sommets().set_md_vector(md);
  dom.les_elems().set_md_vector(md);
}


// Explicit instanciation
template void Scatter::init_sequential_domain(Domaine_32_64<int>& dom);
template void Scatter::uninit_sequential_domain(Domaine_32_64<int>& dom);
#if INT_is_64_ == 2
template void Scatter::init_sequential_domain(Domaine_32_64<trustIdType>& dom);
template void Scatter::uninit_sequential_domain(Domaine_32_64<trustIdType>& dom);
#endif
