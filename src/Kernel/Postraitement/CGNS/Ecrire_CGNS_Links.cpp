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

#include <Comm_Group_MPI.h>
#include <communications.h>
#include <unordered_set>
#include <Ecrire_CGNS.h>
#include <Domaine.h>
#include <unistd.h>

#ifdef HAS_CGNS

/*
 * ***************** *
 * METHODS POUR LINK *
 * ***************** *
 */

void Ecrire_CGNS::init_proc_maitre_local_comm()
{
  assert(Process::is_parallel()
         && (Option_CGNS::LINKED_FILES_PER_COMM_GROUP || Option_CGNS::SINGLE_FILE_PER_COMM_GROUP)
         && PE_Groups::has_user_defined_group());

  const auto& grp = PE_Groups::get_user_defined_group();
  if (PE_Groups::enter_group(grp))
    {
      proc_maitre_local_comm_ = PE_Groups::groupe_TRUST().rank();
      envoyer_broadcast(proc_maitre_local_comm_, 0); // XXX should do this !
      PE_Groups::exit_group();
    }
}

void Ecrire_CGNS::cgns_open_grid_base_link_file()
{
  assert(Option_CGNS::USE_LINKS && !postraiter_domaine_);
  std::string fn;

  if (Process::is_parallel() && Option_CGNS::LINKED_FILES_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
    {
      fn = (Nom(baseFile_name_)).nom_me(proc_maitre_local_comm_).getString() + ".grid.cgns"; // file name
      unlink(fn.c_str());
      cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::PAR>(fn, fileId_, false);
      Cerr << "**** Multiple parallel CGNS files " << baseFile_name_ << "_XXXX.grid.cgns opened !" << finl;
    }
  else
    {
      fn = baseFile_name_ + ".grid.cgns"; // file name
      unlink(fn.c_str());

      if (Process::is_parallel())
        cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::PAR>(fn, fileId_);
      else
        cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_);
    }

  grid_file_opened_ = true;
}

void Ecrire_CGNS::cgns_close_grid_or_solution_link_file(const double t, const TYPE_LINK_CGNS type, bool is_cerr)
{
  assert((Option_CGNS::USE_LINKS && !postraiter_domaine_) || is_lagrangian_);

  std::string fn; // file name

  if (type == TYPE_LINK_CGNS::GRID)
    {
      if (Process::is_parallel() && Option_CGNS::LINKED_FILES_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
        fn =  baseFile_name_ + "_XXXX.grid.cgns";
      else
        fn = baseFile_name_ + ".grid.cgns";
    }
  else if (type == TYPE_LINK_CGNS::SOLUTION)
    {
      if (Process::is_parallel() && Option_CGNS::LINKED_FILES_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
        fn = baseFile_name_ + "_XXXX_" + ".solution." + cgns_helper_.convert_double_to_string(t) + ".cgns";
      else
        fn = baseFile_name_ + ".solution." + cgns_helper_.convert_double_to_string(t) + ".cgns";
    }
  else if (type == TYPE_LINK_CGNS::FINAL_LINK)
    fn = baseFile_name_ + ".cgns";
  else
    Process::exit("Error in Ecrire_CGNS::cgns_close_grid_or_solution_link_file !!! \n");

  if (Process::is_parallel() && (type != TYPE_LINK_CGNS::FINAL_LINK))
    {
      if ( Option_CGNS::LINKED_FILES_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
        {
          cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::PAR>(fn /* inutile */, fileId_, false);
          Cerr << "**** Multiple parallel CGNS files " << fn << " closed !" << finl;
        }
      else
        cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::PAR>(fn, fileId_, is_cerr);
    }
  else
    cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, is_cerr);

  if (type == TYPE_LINK_CGNS::GRID)
    grid_file_opened_ = false;

  if (type == TYPE_LINK_CGNS::SOLUTION)
    solution_file_opened_ = false;
}

void Ecrire_CGNS::cgns_fill_info_grid_link_file(const char* basename, const CGNS_TYPE& cgns_type_elem, const int icelldim, const int nb_som, const int nb_elem, const bool is_polyedre)
{
  cellDim_.push_back(icelldim);
  baseZone_name_.push_back(std::string(basename));
  sizeId_.push_back( { nb_som, nb_elem } );

  if (cgns_type_elem == CGNS_ENUMV(NGON_n)) // cas polyedre
    {
      if (is_polyedre)
        connectname_.push_back( { "NGON_n", "NFACE_n" });
      else
        connectname_.push_back({ "NGON_n" });
    }
  else
    connectname_.push_back({ "Elem" }); // autre cas
}

void Ecrire_CGNS::cgns_init_solution_link_file(const std::string& LOC, const Nom& nom_dom)
{
  assert (LOC == "ELEM" || LOC == "SOM");
  Cerr << "###  Building a new CGNS base with a linked zone to host the field located at : " << LOC << " !" << finl;
  doms_written_.push_back(nom_dom);
  baseId_.push_back(-123); // pour chaque dom, on a une baseId
  zoneId_.push_back(-123);
}

void Ecrire_CGNS::gather_local_sizeId_for_comm_group()
{
#ifdef MPI_
  if (!vec_proc_maitre_local_comm_.empty()) return; /* rien a faire */

  unique_vec_proc_maitre_local_comm_.clear();
  sizeId_som_local_comm_.clear();
  sizeId_elem_local_comm_.clear();

  vec_proc_maitre_local_comm_.assign(Process::nproc(), -123 /* default */);
  MPI_Allgather(&proc_maitre_local_comm_, 1, MPI_ENTIER, vec_proc_maitre_local_comm_.data(), 1, MPI_ENTIER, Comm_Group_MPI::get_trio_u_world());

  std::unordered_set<int> seen;

  for (int val : vec_proc_maitre_local_comm_)
    if (seen.insert(val).second)
      unique_vec_proc_maitre_local_comm_.push_back(val); // si val pas dedans

  std::vector<std::vector<cgsize_t>> sizeId_som_local_comm_tmp, sizeId_elem_local_comm_tmp;

  MPI_Datatype CGNS_MPI_SIZE;
  MPI_Type_match_size(MPI_TYPECLASS_INTEGER, sizeof(cgsize_t), &CGNS_MPI_SIZE);

  for (int i = 0; i < static_cast<int>(sizeId_.size()); i++)
    {
      assert(sizeId_[i].size() == 2);

      sizeId_som_local_comm_tmp.emplace_back(Process::nproc(), static_cast<cgsize_t>(-123));
      sizeId_elem_local_comm_tmp.emplace_back(Process::nproc(), static_cast<cgsize_t>(-123));

      MPI_Allgather(&sizeId_[i][0], 1, CGNS_MPI_SIZE, sizeId_som_local_comm_tmp[i].data(), 1, CGNS_MPI_SIZE, Comm_Group_MPI::get_trio_u_world());
      MPI_Allgather(&sizeId_[i][1], 1, CGNS_MPI_SIZE, sizeId_elem_local_comm_tmp[i].data(), 1, CGNS_MPI_SIZE, Comm_Group_MPI::get_trio_u_world());
    }

  const int nb_grps = static_cast<int>(unique_vec_proc_maitre_local_comm_.size());

  for (int i = 0; i < static_cast<int>(sizeId_.size()); i++)
    {
      sizeId_som_local_comm_.emplace_back(nb_grps, static_cast<cgsize_t>(-123));
      sizeId_elem_local_comm_.emplace_back(nb_grps, static_cast<cgsize_t>(-123));

      for (int j = 0; j < nb_grps; j++)
        {
          int proc_grp = unique_vec_proc_maitre_local_comm_[j];
          sizeId_som_local_comm_[i][j] = sizeId_som_local_comm_tmp[i][proc_grp];
          sizeId_elem_local_comm_[i][j] = sizeId_elem_local_comm_tmp[i][proc_grp];
        }
    }
#endif
}

/* Used to create a support in a single CGNS file
 * ie : called once at the beginning while filling fld_loc_map_ */
void Ecrire_CGNS::add_new_linked_base(const std::string& LOC, const Nom& nom_dom)
{
  assert (LOC == "ELEM" || LOC == "SOM");
  Cerr << "###  Building a new CGNS base with a linked zone to host the field located at : " << LOC << " !" << finl;

  doms_written_.push_back(nom_dom);
  baseId_.push_back(-123); // pour chaque dom, on a une baseId

  const Nom nom_dom_mod = TRUST_2_CGNS::modify_domaine_name_for_link(nom_dom, LOC);
  const int ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom_mod); // get index of orig dom

  if (cg_base_write(fileId_, nom_dom.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_.back()) != CG_OK)
    Cerr << "Error Ecrire_CGNS::add_new_linked_base : cg_base_write !" << finl, TRUST_CGNS_ERROR();

  if (Process::is_parallel() && (Option_CGNS::PARALLEL_OVER_ZONE || postraiter_domaine_))
    {
      zoneId_.clear(); // XXX commencons par ca
      TRUST_2_CGNS& TRUST2CGNS = T2CGNS_[ind_base];
      const std::vector<int>& global_nb_elem = TRUST2CGNS.get_global_nb_elem(),
                              &global_nb_som = TRUST2CGNS.get_global_nb_som(),
                               &proc_non_zero_elem = TRUST2CGNS.get_proc_non_zero_elem();

      const int nb_zones_to_write = TRUST2CGNS.nb_procs_writing();
      const bool all_write = TRUST2CGNS.all_procs_write(); // all procs will write !
      std::string zonename, zonename_link;

      for (int i = 0; i != nb_zones_to_write; i++)
        {
          const int indZ = all_write ? i : proc_non_zero_elem[i]; // procID
          const int ne_loc = global_nb_elem[indZ], ns_loc = global_nb_som[indZ]; /* nb_elem & nb_som local */
          assert (ne_loc > 0);

          zoneId_.push_back(-123);
          cgsize_t isize[3] = { ns_loc , ne_loc , 0 }; /* 0 => boundary vertex size (zero if elements not sorted) */

          zonename = nom_dom.nom_me(indZ).getString();
          zonename_link = nom_dom_mod.nom_me(indZ).getString();

          cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_, baseId_.back(), zonename, isize, zoneId_.back(), indZ + 1,
                                                         "" /* meme fichier */, nom_dom_mod.getString(), zonename_link, connectname_[ind_base],
                                                         "Ecrire_CGNS::add_new_linked_base");
        }

      zoneId_par_.push_back(zoneId_); // XXX : Dont touch
    }
  else
    {
      zoneId_.push_back(-123);
      cgsize_t isize[3] = { sizeId_[ind_base][0] , sizeId_[ind_base][1] , 0 };

      cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_, baseId_.back(), nom_dom.getString(), isize, zoneId_.back(), 1,
                                                     "" /* meme fichier */, baseZone_name_[ind_base], baseZone_name_[ind_base], connectname_[ind_base],
                                                     "Ecrire_CGNS::add_new_linked_base");
    }
}

void Ecrire_CGNS::cgns_open_solution_link_file(const double t)
{
  assert((Option_CGNS::USE_LINKS && !postraiter_domaine_) || is_lagrangian_);

  const bool enter_group_comm = Process::is_parallel() && Option_CGNS::LINKED_FILES_PER_COMM_GROUP && PE_Groups::has_user_defined_group();

  std::string fn;

  if (enter_group_comm)
    fn = (Nom(baseFile_name_)).nom_me(proc_maitre_local_comm_).getString() + ".solution." + cgns_helper_.convert_double_to_string(t) + ".cgns"; // file name
  else
    fn = baseFile_name_ + ".solution." + cgns_helper_.convert_double_to_string(t) + ".cgns"; // file name

  unlink(fn.c_str());

  if (Process::is_parallel())
    {
      cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::PAR>(fn, fileId_, enter_group_comm ? false : true);

      if (enter_group_comm)
        Cerr << "**** Multiple parallel CGNS files " << baseFile_name_ << "_XXXX.solution." + cgns_helper_.convert_double_to_string(t) + ".cgns opened !" << finl;
    }
  else
    cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);

  solution_file_opened_ = true;

  if (is_deformable_)
    return; /* Stop here if deformable */

  /* Otherwise, we have a grid file already written ... we link the zones in the opened solution files to it (coords + connectivity) ! */
  for (auto &itr : fld_loc_map_)
    {
      const std::string& LOC = itr.first;
      const Nom& nom_dom = itr.second;

      const int index_glob = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
      int ind_base = index_glob;

      if (LOC != "FACES")
        {
          const Nom nom_dom_mod = TRUST_2_CGNS::modify_domaine_name_for_link(nom_dom, LOC);
          ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom_mod);
        }

      if (cg_base_write(fileId_, nom_dom.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_open_solution_link_file : cg_base_write !" << finl, TRUST_CGNS_ERROR();

      /*
       * XXX this is done in //, not like final link file which is done only on proc 0
       * So no need to get sizes per local comm ... each proc available on its comm group take the good values
       */
      cgsize_t isize[3] = { sizeId_[ind_base][0] , sizeId_[ind_base][1] , 0 };

      std::string linkfile = baseFile_name_ + ".grid.cgns"; // file name

      if (enter_group_comm)
        linkfile = (Nom(baseFile_name_)).nom_me(proc_maitre_local_comm_).getString() + ".grid.cgns"; // file name

      TRUST_2_CGNS::remove_slash_linkfile(linkfile);

      cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_, baseId_[index_glob], nom_dom.getString(), isize, zoneId_[index_glob], 1,
                                                     linkfile, baseZone_name_[ind_base], baseZone_name_[ind_base], connectname_[ind_base],
                                                     "Ecrire_CGNS::cgns_open_solution_link_file");
    }
}

void Ecrire_CGNS::cgns_write_final_link_file_comm_group()
{
#ifdef MPI_
  if (vec_proc_maitre_local_comm_.empty())
    gather_local_sizeId_for_comm_group();

  if (!Process::me())
    {
      std::string fn = baseFile_name_ + ".cgns";

      unlink(fn.c_str());
      cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);

      const int nb_grps = static_cast<int>(unique_vec_proc_maitre_local_comm_.size());
      std::vector<int> zoneId_tmp(nb_grps, -123);
      std::vector<int> ind_doms_dumped;

      /* 1 : on iter juste sur le map fld_loc_map_; ie: pas domaine dis ... */
      for (auto& itr : fld_loc_map_)
        {
          const std::string& LOC = itr.first;
          const Nom& nom_dom = itr.second;

          const int index_glob = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
          ind_doms_dumped.push_back(index_glob);
          int ind_base = index_glob;

          if (LOC != "FACES")
            {
              const Nom nom_dom_mod = TRUST_2_CGNS::modify_domaine_name_for_link(nom_dom, LOC);
              ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom_mod);
            }

          if (cg_base_write(fileId_, nom_dom.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_base_write !" << finl, TRUST_CGNS_ERROR();

          for (int gid = 0; gid < nb_grps; gid++)
            {
              int proc_grp = unique_vec_proc_maitre_local_comm_[gid];
              std::string zone_name = Nom("Zone").nom_me(proc_grp).getString();

              std::string file_group_id = Nom(baseFile_name_).nom_me(proc_grp).getString();
              TRUST_2_CGNS::remove_slash_linkfile(file_group_id);

              cgsize_t isize[3];
              isize[0] = sizeId_som_local_comm_[ind_base][gid];
              isize[1] = sizeId_elem_local_comm_[ind_base][gid];
              isize[2] = 0;

              cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_, baseId_[index_glob], zone_name, isize, zoneId_tmp[gid], gid + 1,
                                                             file_group_id + ".grid.cgns" /* linkfile */, baseZone_name_[ind_base], baseZone_name_[ind_base], connectname_[ind_base],
                                                             "Ecrire_CGNS::cgns_write_final_link_file_comm_group");


              cgns_helper_.cgns_write_solution_classic_links(file_group_id, nom_dom.getString(), nom_dom.getString(), LOC, time_post_,
                                                             "Ecrire_CGNS::cgns_write_final_link_file");
            }

          cgns_helper_.cgns_write_iters<TYPE_ECRITURE_CGNS::SEQ>(true /* has_field */, nb_grps /* nb_zones_to_write */, fileId_, baseId_[index_glob], ind_base,
                                                                 zoneId_tmp, LOC, solname_som_, solname_elem_, solname_faces_, time_post_);

        }

      /* 2 : on iter sur les autres domaines; ie: domaine dis */
      for (int i = 0; i < static_cast<int>(doms_written_.size()); i++)
        {
          if (std::find(ind_doms_dumped.begin(), ind_doms_dumped.end(), i) == ind_doms_dumped.end()) // indice pas dans ind_doms_dumped
            {
              const Nom& nom_dom = doms_written_[i];
              const int ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
              assert(ind_base > -1);

              if (cg_base_write(fileId_, nom_dom.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[ind_base]) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_base_write !" << finl, TRUST_CGNS_ERROR();

              for (int gid = 0; gid < nb_grps; gid++)
                {
                  const int proc_grp = unique_vec_proc_maitre_local_comm_[gid];
                  std::string zone_name = Nom("Zone").nom_me(proc_grp).getString();

                  std::string file_group_id = Nom(baseFile_name_).nom_me(proc_grp).getString();
                  TRUST_2_CGNS::remove_slash_linkfile(file_group_id);

                  cgsize_t isize[3];
                  isize[0] = sizeId_som_local_comm_[ind_base][gid];
                  isize[1] = sizeId_elem_local_comm_[ind_base][gid];
                  isize[2] = 0;

                  const bool write_connectivity = (!(isize[0] == 1 && isize[1] == 1));

                  /* He we dont link to solutions since no fields ... just other domais dis ;) */
                  cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_, baseId_[ind_base], zone_name, isize, zoneId_tmp[gid], gid + 1,
                                                                 file_group_id + ".grid.cgns" /* linkfile */, baseZone_name_[ind_base], baseZone_name_[ind_base], connectname_[ind_base],
                                                                 "Ecrire_CGNS::cgns_write_final_link_file", write_connectivity);
                }

              cgns_helper_.cgns_write_iters<TYPE_ECRITURE_CGNS::SEQ>(false /* has_field */, nb_grps /* nb_zones_to_write */, fileId_, baseId_[ind_base], ind_base,
                                                                     zoneId_tmp, "rien", solname_som_, solname_elem_,solname_faces_, time_post_);
            }
        }
      cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);
    }
#endif
}

void Ecrire_CGNS::cgns_write_final_link_file()
{
  if (Process::is_parallel() && Option_CGNS::LINKED_FILES_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
    {
      cgns_write_final_link_file_comm_group();
      return;
    }

  /* Only master proc writes the link file ! */
  if (!Process::me())
    {
      // Fichier link maintenant
      std::string fn = baseFile_name_ + ".cgns"; // file name
      unlink(fn.c_str());
      cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);

      std::string base_link_file = baseFile_name_;
      TRUST_2_CGNS::remove_slash_linkfile(base_link_file);

      std::vector<int> ind_doms_dumped;

      /* 1 : on iter juste sur le map fld_loc_map_; ie: pas domaine dis ... */
      for (auto& itr : fld_loc_map_)
        {
          const std::string& LOC = itr.first;
          const Nom& nom_dom = itr.second;
          const int index_glob = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
          ind_doms_dumped.push_back(index_glob);
          assert(index_glob > -1);

          int ind_base = index_glob;
          if (LOC != "FACES")
            {
              const Nom nom_dom_mod = TRUST_2_CGNS::modify_domaine_name_for_link(nom_dom, LOC);
              ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom_mod);
            }

          if (cg_base_write(fileId_, nom_dom.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_open_solution_link_file : cg_base_write !" << finl, TRUST_CGNS_ERROR();

          cgsize_t isize[3] = { sizeId_[ind_base][0] , sizeId_[ind_base][1] , 0 };

          cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_, baseId_[index_glob], nom_dom.getString(), isize, zoneId_[index_glob], 1,
                                                         base_link_file + ".grid.cgns" /* linkfile */, baseZone_name_[ind_base], baseZone_name_[ind_base], connectname_[ind_base],
                                                         "Ecrire_CGNS::cgns_write_final_link_file");

          cgns_helper_.cgns_write_solution_classic_links(base_link_file, nom_dom.getString(), nom_dom.getString(), LOC, time_post_,
                                                         "Ecrire_CGNS::cgns_write_final_link_file");

          cgns_helper_.cgns_write_iters<TYPE_ECRITURE_CGNS::SEQ>(true /* has_field */, 1, fileId_, baseId_[index_glob], index_glob /* 1st Zone */,
                                                                 zoneId_, LOC, solname_som_, solname_elem_, solname_faces_, time_post_);

        }

      /* 2 : on iter sur les autres domaines; ie: domaine dis */
      for (int i = 0; i < static_cast<int>(doms_written_.size()); i++)
        {
          if (std::find(ind_doms_dumped.begin(), ind_doms_dumped.end(), i) == ind_doms_dumped.end()) // indice pas dans ind_doms_dumped
            {
              const Nom& nom_dom = doms_written_[i];
              const int ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
              assert(ind_base > -1);

              if (cg_base_write(fileId_, nom_dom.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[ind_base]) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_open_solution_link_file : cg_base_write !" << finl, TRUST_CGNS_ERROR();

              cgsize_t isize[3] = { sizeId_[ind_base][0] , sizeId_[ind_base][1] , 0 };

              /* He we dont link to solutions since no fields ... just other domais dis ;) */
              cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_, baseId_[ind_base], nom_dom.getString(), isize, zoneId_[ind_base], 1,
                                                             base_link_file + ".grid.cgns" /* linkfile */, baseZone_name_[ind_base], baseZone_name_[ind_base], connectname_[ind_base],
                                                             "Ecrire_CGNS::cgns_write_final_link_file");


              cgns_helper_.cgns_write_iters<TYPE_ECRITURE_CGNS::SEQ>(false /* has_field */, 1 /* nb_zones_to_write */, fileId_, baseId_[ind_base], ind_base,
                                                                     zoneId_, "rien", solname_som_, solname_elem_,solname_faces_, time_post_);
            }
        }

      cgns_close_grid_or_solution_link_file(-123. /* inutile*/, TYPE_LINK_CGNS::FINAL_LINK, true); // on ferme
    }
}

#endif /* HAS_CGNS */
