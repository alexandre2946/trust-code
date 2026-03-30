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
          const cgsize_t isize[3] = { ns_loc , ne_loc , 0 }; /* 0 => boundary vertex size (zero if elements not sorted) */

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
      const cgsize_t isize[3] = { sizeId_[ind_base][0] , sizeId_[ind_base][1] , 0 };

      cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_, baseId_.back(), nom_dom.getString(), isize, zoneId_.back(), 1,
                                                     "" /* meme fichier */, baseZone_name_[ind_base], baseZone_name_[ind_base], connectname_[ind_base],
                                                     "Ecrire_CGNS::add_new_linked_base");
    }
}

/* Used to open a new solution file, link all required supports to grid file (elem/som)
 * ie : called once at the beginning while filling fld_loc_map_, and then at each add_time */
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

      int index_glob = -123, ind_base = -123;
      TRUST_2_CGNS::init_glob_base_domain_idx(doms_written_, nom_dom, true /* has_fields */, LOC, index_glob, ind_base);

      if (cg_base_write(fileId_, nom_dom.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_open_solution_link_file : cg_base_write !" << finl, TRUST_CGNS_ERROR();

      /*
       * XXX this is done in //, not like final link file which is done only on proc 0
       * So no need to get sizes per local comm ... each proc available on its comm group take the good values
       */
      const cgsize_t isize[3] = { sizeId_[ind_base][0] , sizeId_[ind_base][1] , 0 };

      std::string linkfile = baseFile_name_ + ".grid.cgns"; // file name

      if (enter_group_comm)
        linkfile = (Nom(baseFile_name_)).nom_me(proc_maitre_local_comm_).getString() + ".grid.cgns"; // file name

      TRUST_2_CGNS::remove_slash_linkfile(linkfile);

      cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_, baseId_[index_glob], nom_dom.getString(), isize, zoneId_[index_glob], 1,
                                                     linkfile, baseZone_name_[ind_base], baseZone_name_[ind_base], connectname_[ind_base],
                                                     "Ecrire_CGNS::cgns_open_solution_link_file");
    }
}

void Ecrire_CGNS::cgns_write_final_link_file_for_single_file_comm_group()
{
#ifdef MPI_
  if (vec_proc_maitre_local_comm_.empty())
    gather_local_sizeId_for_comm_group();

  int fileId_link;
  std::string fn = baseFile_name_ + ".cgns";
  unlink(fn.c_str());
  cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::PAR>(fn, fileId_link, true);

  const int nb_grps = static_cast<int>(unique_vec_proc_maitre_local_comm_.size());
  std::vector<int> zoneId_tmp;//(nb_grps, -123);
  std::vector<int> baseId_tmp(baseId_);

  for (auto &itr : doms_written_)
    {
      bool has_field = false;
      std::string LOC = "rien";
      TRUST_2_CGNS::init_has_field_and_loc_iters(itr, fld_loc_map_, has_field, LOC);

      int index_glob = -123, ind_base = -123;
      TRUST_2_CGNS::init_glob_base_domain_idx(doms_written_, itr, has_field, LOC, index_glob, ind_base);

      if (cg_base_write(fileId_link, itr.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_tmp[index_glob]) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_base_write !" << finl, TRUST_CGNS_ERROR();

      int zone_goto_idx = 1;
      zoneId_tmp.clear();

      /* Loop on groups that have something to link to ;) */
      for (int gid = 0; gid < nb_grps; gid++)
        {
          int proc_grp = unique_vec_proc_maitre_local_comm_[gid];
          std::string zone_name = Nom("Zone").nom_me(proc_grp).getString();

          std::string file_group_id = Nom(baseFile_name_).nom_me(proc_grp).getString();
          TRUST_2_CGNS::remove_slash_linkfile(file_group_id);

          const cgsize_t isize[3] = { sizeId_som_local_comm_[ind_base][gid], sizeId_elem_local_comm_[ind_base][gid], 0 };

          if ((isize[0] == 0 && isize[1] == 0)) continue;

          zoneId_tmp.push_back(-123);

          if (is_deformable_)
            {
              cgns_helper_.cgns_write_zone_and_deformable_links(true /* write zone */, has_field, fileId_link, baseId_tmp[index_glob],
                                                                zone_name, isize, zoneId_tmp.back(), zone_goto_idx,
                                                                file_group_id, baseZone_name_[ind_base], baseZone_name_[ind_base],
                                                                connectname_[ind_base], itr, LOC, time_post_,
                                                                "Ecrire_CGNS::cgns_write_final_link_file_comm_group");
            }
          else
            {
              cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_link, baseId_tmp[index_glob], zone_name,
                                                             isize, zoneId_tmp.back(), zone_goto_idx,
                                                             file_group_id + ".cgns" /* linkfile */, baseZone_name_[ind_base],
                                                             baseZone_name_[ind_base], connectname_[ind_base],
                                                             "Ecrire_CGNS::cgns_write_final_link_file_for_single_file_comm_group");

              if(has_field)
                cgns_helper_.cgns_write_solution_classic_links(file_group_id, itr.getString(), itr.getString(), LOC, time_post_,
                                                               "Ecrire_CGNS::cgns_write_final_link_file_comm_group");
            }

          zone_goto_idx++;
        }

      /* Finally, write the iters ;) */
      if (is_deformable_)
        cgns_helper_.cgns_write_iters_deformable<TYPE_ECRITURE_CGNS::SEQ>(true, has_field, static_cast<int>(zoneId_tmp.size()) /* nb_zones_to_write */,
                                                                          fileId_link, baseId_tmp[index_glob], ind_base,
                                                                          zoneId_tmp, LOC, solname_som_, solname_elem_,
                                                                          solname_faces_, grid_name_, time_post_);
      else
        cgns_helper_.cgns_write_iters<TYPE_ECRITURE_CGNS::SEQ>(has_field, static_cast<int>(zoneId_tmp.size()) /* nb_zones_to_write */,
                                                               fileId_link, baseId_tmp[index_glob], ind_base,
                                                               zoneId_tmp, LOC, solname_som_, solname_elem_,
                                                               solname_faces_, time_post_);
    }
  cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::PAR>(fn, fileId_link, true);
#endif
}

void Ecrire_CGNS::cgns_write_final_link_file_comm_group()
{
#ifdef MPI_
  if (vec_proc_maitre_local_comm_.empty())
    gather_local_sizeId_for_comm_group();

  /* Fichier link : Only master proc writes the link file ! */
  if (!Process::me())
    {
      std::string fn = baseFile_name_ + ".cgns";
      unlink(fn.c_str());
      cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);

      const int nb_grps = static_cast<int>(unique_vec_proc_maitre_local_comm_.size());
      std::vector<int> zoneId_tmp;//(nb_grps, -123);

      for (auto &itr : doms_written_)
        {
          bool has_field = false;
          std::string LOC = "rien";
          TRUST_2_CGNS::init_has_field_and_loc_iters(itr, fld_loc_map_, has_field, LOC);

          int index_glob = -123, ind_base = -123;
          TRUST_2_CGNS::init_glob_base_domain_idx(doms_written_, itr, has_field, LOC, index_glob, ind_base);

          if (cg_base_write(fileId_, itr.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_base_write !" << finl, TRUST_CGNS_ERROR();

          int zone_goto_idx = 1;
          zoneId_tmp.clear();

          /* Loop on groups that have something to link to ;) */
          for (int gid = 0; gid < nb_grps; gid++)
            {
              int proc_grp = unique_vec_proc_maitre_local_comm_[gid];
              std::string zone_name = Nom("Zone").nom_me(proc_grp).getString();

              std::string file_group_id = Nom(baseFile_name_).nom_me(proc_grp).getString();
              TRUST_2_CGNS::remove_slash_linkfile(file_group_id);

              const cgsize_t isize[3] = { sizeId_som_local_comm_[ind_base][gid], sizeId_elem_local_comm_[ind_base][gid], 0 };

              if ((isize[0] == 0 && isize[1] == 0)) continue;

              zoneId_tmp.push_back(-123);

              if (is_deformable_)
                {
                  cgns_helper_.cgns_write_zone_and_deformable_links(true /* write zone */, has_field, fileId_, baseId_[index_glob],
                                                                    zone_name, isize, zoneId_tmp.back(), zone_goto_idx,
                                                                    file_group_id, baseZone_name_[ind_base], baseZone_name_[ind_base],
                                                                    connectname_[ind_base], itr, LOC, time_post_,
                                                                    "Ecrire_CGNS::cgns_write_final_link_file_comm_group");
                }
              else
                {
                  cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_, baseId_[index_glob], zone_name,
                                                                 isize, zoneId_tmp.back(), zone_goto_idx,
                                                                 file_group_id + ".grid.cgns" /* linkfile */, baseZone_name_[ind_base],
                                                                 baseZone_name_[ind_base], connectname_[ind_base],
                                                                 "Ecrire_CGNS::cgns_write_final_link_file_comm_group");

                  if(has_field)
                    cgns_helper_.cgns_write_solution_classic_links(file_group_id, itr.getString(), itr.getString(), LOC, time_post_,
                                                                   "Ecrire_CGNS::cgns_write_final_link_file_comm_group");
                }

              zone_goto_idx++;
            }

          /* Finally, write the iters ;) */
          if (is_deformable_)
            cgns_helper_.cgns_write_iters_deformable<TYPE_ECRITURE_CGNS::SEQ>(true, has_field, static_cast<int>(zoneId_tmp.size()) /* nb_zones_to_write */,
                                                                              fileId_, baseId_[index_glob], ind_base,
                                                                              zoneId_tmp, LOC, solname_som_, solname_elem_,
                                                                              solname_faces_, grid_name_, time_post_);
          else
            cgns_helper_.cgns_write_iters<TYPE_ECRITURE_CGNS::SEQ>(has_field, static_cast<int>(zoneId_tmp.size()) /* nb_zones_to_write */,
                                                                   fileId_, baseId_[index_glob], ind_base,
                                                                   zoneId_tmp, LOC, solname_som_, solname_elem_,
                                                                   solname_faces_, time_post_);
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

  /* Fichier link : Only master proc writes the link file ! */
  if (!Process::me())
    {
      std::string fn = baseFile_name_ + ".cgns"; // file name
      unlink(fn.c_str());
      cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);

      std::string base_link_file = baseFile_name_;
      TRUST_2_CGNS::remove_slash_linkfile(base_link_file);

      for (auto &itr : doms_written_)
        {
          bool has_field = false;
          std::string LOC = "rien";
          TRUST_2_CGNS::init_has_field_and_loc_iters(itr, fld_loc_map_, has_field, LOC);

          int index_glob = -123, ind_base = -123;
          TRUST_2_CGNS::init_glob_base_domain_idx(doms_written_, itr, has_field, LOC, index_glob, ind_base);

          if (cg_base_write(fileId_, itr.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_open_solution_link_file : cg_base_write !" << finl, TRUST_CGNS_ERROR();

          const cgsize_t isize[3] = { sizeId_[ind_base][0] , sizeId_[ind_base][1] , 0 };

          if (is_deformable_)
            {
              cgns_helper_.cgns_write_zone_and_deformable_links(true /* write zone */, has_field, fileId_, baseId_[index_glob],
                                                                itr.getString(), isize, zoneId_[index_glob], 1,
                                                                base_link_file, baseZone_name_[ind_base], baseZone_name_[ind_base],
                                                                connectname_[ind_base], itr, LOC, time_post_,
                                                                "Ecrire_CGNS::cgns_open_solution_link_file");

              cgns_helper_.cgns_write_iters_deformable<TYPE_ECRITURE_CGNS::SEQ>(true /* deformable */, has_field, 1 /* 1 zone per base */, fileId_,
                                                                                baseId_[index_glob], index_glob /* 1st Zone */,
                                                                                zoneId_, LOC, solname_som_, solname_elem_,
                                                                                solname_faces_, grid_name_, time_post_);
            }
          else
            {
              cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_, baseId_[index_glob],
                                                             itr.getString(), isize, zoneId_[index_glob], 1,
                                                             base_link_file + ".grid.cgns" /* linkfile */, baseZone_name_[ind_base],
                                                             baseZone_name_[ind_base], connectname_[ind_base],
                                                             "Ecrire_CGNS::cgns_write_final_link_file");

              if(has_field)
                cgns_helper_.cgns_write_solution_classic_links(base_link_file, itr.getString(), itr.getString(), LOC, time_post_,
                                                               "Ecrire_CGNS::cgns_write_final_link_file");

              cgns_helper_.cgns_write_iters<TYPE_ECRITURE_CGNS::SEQ>(has_field, 1, fileId_, baseId_[index_glob], index_glob /* 1st Zone */,
                                                                     zoneId_, LOC, solname_som_, solname_elem_, solname_faces_, time_post_);

            }
        }

      cgns_close_grid_or_solution_link_file(-123. /* inutile*/, TYPE_LINK_CGNS::FINAL_LINK, true); // on ferme
    }
}

void Ecrire_CGNS::link_multi_loc_support_pb_deformable()
{
  if (is_lagrangian_)
    {
      link_multi_loc_support_lagrangian();
      return;
    }

  const bool enter_group_comm = Process::is_parallel() && Option_CGNS::LINKED_FILES_PER_COMM_GROUP &&
                                PE_Groups::has_user_defined_group() && !postraiter_domaine_;

  // loop and write linked supports !
  if (Option_CGNS::USE_LINKS)
    {
      for (auto &itr : fld_loc_map_)
        {
          const std::string& LOC = itr.first;
          const Nom& nom_dom = itr.second;

          int index_glob = -123, ind_base = -123;
          TRUST_2_CGNS::init_glob_base_domain_idx(doms_written_, nom_dom, true /* has_field */, LOC, index_glob, ind_base);

          if (cg_base_write(fileId_, nom_dom.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::link_multi_loc_support_pb_deformable : cg_base_write !" << finl, TRUST_CGNS_ERROR();

          const cgsize_t isize[3] = { sizeId_[ind_base][0], sizeId_[ind_base][1], 0 };

          // XXX we use the helper but we dont write connectivity because it is a bit special : we link grid coords to same file while connectivity to 1st sol file ...
          cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_, baseId_[index_glob], nom_dom.getString(), isize, zoneId_[index_glob], 1,
                                                         "" /* this file */, baseZone_name_[ind_base], baseZone_name_[ind_base], connectname_[ind_base],
                                                         "Ecrire_CGNS::link_multi_loc_support_pb_deformable", false /* DONT WRITE CONN */);

          // Write conn : linked to 1st solution file since deformable ...
          std::string linkfile;
          if (!first_time_post_)
            linkfile = (enter_group_comm ? Nom(baseFile_name_).nom_me(proc_maitre_local_comm_).getString() : baseFile_name_) +
                       ".solution." + cgns_helper_.convert_double_to_string(time_post_[0]) + ".cgns";

          TRUST_2_CGNS::remove_slash_linkfile(linkfile);

          cgns_helper_.cgns_write_connectivity_deformable_links(fileId_, baseId_[index_glob], zoneId_[index_glob], linkfile, baseZone_name_[ind_base], baseZone_name_[ind_base],
                                                                connectname_[ind_base], "Ecrire_CGNS::link_multi_loc_support_pb_deformable");
        }
    }
  else
    {
      if (first_time_post_) return; /* Mais ouiiiiii car fait dans cgns_fill_field_loc_map !! */

      for (auto &itr : fld_loc_map_)
        {
          const std::string& LOC = itr.first;
          const Nom& nom_dom = itr.second;

          int index_glob = -123, ind_base = -123;
          TRUST_2_CGNS::init_glob_base_domain_idx(doms_written_, nom_dom, true /* has_field */, LOC, index_glob, ind_base);

          std::string linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/" + grid_name_loc_ + "/";

          if (cg_goto(fileId_, baseId_[index_glob], "Zone_t", 1, "end") != CG_OK)
            Cerr << "Error Ecrire_CGNS::link_multi_loc_support_pb_deformable : cg_goto !" << finl, TRUST_CGNS_ERROR();

          if (cg_link_write(grid_name_loc_.c_str(), "" /* this file */, linkpath.c_str()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::link_multi_loc_support_pb_deformable : cg_link_write !" << finl, TRUST_CGNS_ERROR();
        }
    }

  multi_loc_deformable_support_linked_ = true; // of course !
}

// Specifique FT !!
void Ecrire_CGNS::link_multi_loc_support_lagrangian()
{
  for (auto &itr : fld_loc_map_)
    {
      const std::string& LOC = itr.first;
      const Nom& nom_dom = itr.second;

      int index_glob = -123, ind_base = -123;
      TRUST_2_CGNS::init_glob_base_domain_idx(doms_written_, nom_dom, true /* has_field */, LOC, index_glob, ind_base);
      assert (ind_base == 0);

      if (cg_base_write(fileId_, nom_dom.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
        Cerr << "Error Ecrire_CGNS::link_multi_loc_support_lagrangian : cg_base_write !" << finl, TRUST_CGNS_ERROR();

      const int nb_current_post = static_cast<int> (time_post_.size());
      const cgsize_t isize[3] = { sizeId_[nb_current_post - 1][0], sizeId_[nb_current_post - 1][1], 0 };

      cgns_helper_.cgns_write_zone_and_classic_links(true /* write_zone */, fileId_, baseId_[index_glob], nom_dom.getString(), isize, zoneId_[index_glob], 1,
                                                     "" /* this file */, baseZone_name_[ind_base], baseZone_name_[ind_base], connectname_[ind_base],
                                                     "Ecrire_CGNS::link_multi_loc_support_lagrangian");
    }
  multi_loc_deformable_support_linked_ = true; // of course !
}

void Ecrire_CGNS::cgns_write_final_link_file_lagrangian()
{
  if (Process::me()) return; // seul le proc 0 ecrit le fichier link

  const int nsteps = static_cast<int>(time_post_.size());
  const cgsize_t nuse = static_cast<cgsize_t>(nsteps);

  if (nsteps == 0) return;

  std::string fn = baseFile_name_ + ".cgns"; // file name
  unlink(fn.c_str());
  cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);

  for (auto &itr : fld_loc_map_)
    {
      const std::string& LOC = itr.first;
      const Nom& nom_dom = itr.second;

      int index_glob = -123, ind_base = -123;
      TRUST_2_CGNS::init_glob_base_domain_idx(doms_written_, nom_dom, true /* has_field */, LOC, index_glob, ind_base);

      if (cg_base_write(fileId_, nom_dom.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_lagrangian : cg_base_write !" << finl, TRUST_CGNS_ERROR();

      if (cg_biter_write(fileId_, baseId_[index_glob], "TimeIterValues", nsteps) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_lagrangian : cg_biter_write !" << finl, TRUST_CGNS_ERROR();

      if (cg_goto(fileId_, baseId_[index_glob], "BaseIterativeData_t", 1, "end") != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_lagrangian : cg_goto BaseIterativeData_t !" << finl, TRUST_CGNS_ERROR();

      // TimeValues
      if (cg_array_write("TimeValues", CGNS_DOUBLE_TYPE, 1, &nuse, time_post_.data()) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_lagrangian : cg_array_write TimeValues !" << finl, TRUST_CGNS_ERROR();

      if (cg_simulation_type_write(fileId_, baseId_[index_glob], CGNS_ENUMV(TimeAccurate)) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_lagrangian : cg_simulation_type_write !" << finl, TRUST_CGNS_ERROR();

      // NumberOfZones : 1 zone active par step
      std::vector<int> number_of_zones(nsteps, 1);
      if (cg_array_write("NumberOfZones", CGNS_ENUMV(Integer), 1, &nuse, number_of_zones.data()) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_lagrangian : cg_array_write NumberOfZones !" << finl, TRUST_CGNS_ERROR();

      cgsize_t zpdims[3] = { CGNS_STR_SIZE, 1, nuse };
      std::string zone_ptrs;
      zone_ptrs.reserve(static_cast<size_t>(CGNS_STR_SIZE) * nsteps);

      bool first_zone = true;
      for (int i = 0; i < nsteps; i++)
        {
          std::string zname = nom_dom.getString();
          if (!first_zone)
            {
              zname += "_itr_";
              zname += std::to_string(i);
            }
          first_zone = false;

          zname.resize(CGNS_STR_SIZE, ' ');
          zone_ptrs += zname;
        }

      if (cg_array_write("ZonePointers", CGNS_ENUMV(Character), 3, zpdims, zone_ptrs.c_str()) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_lagrangian : cg_array_write ZonePointers !" << finl, TRUST_CGNS_ERROR();

      for (int i = 0; i < nsteps; i++)
        {
          cgsize_t isize[3] = { sizeId_[i][0] , sizeId_[i][1] , 0 };
          std::string zn = nom_dom.getString();

          if (i > 0)
            {
              zn += "_itr_";
              zn += std::to_string(i);
            }

          if (cg_zone_write(fileId_, baseId_[index_glob], zn.c_str(), isize, CGNS_ENUMV(Unstructured), &zoneId_[index_glob]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_lagrangian : cg_zone_write !" << finl, TRUST_CGNS_ERROR();

          std::string linkfile = baseFile_name_ + ".solution." + cgns_helper_.convert_double_to_string(time_post_[i]) + ".cgns";
          TRUST_2_CGNS::remove_slash_linkfile(linkfile);
          std::string linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/GridCoordinates/";

          if (cg_goto(fileId_, baseId_[index_glob], "Zone_t", zoneId_[index_glob], "end") != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_lagrangian : cg_goto Zone_t !" << finl, TRUST_CGNS_ERROR();

          if (cg_link_write("GridCoordinates", linkfile.c_str(), linkpath.c_str()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_lagrangian : cg_link_write GridCoordinates !" << finl, TRUST_CGNS_ERROR();

          for (auto& itr_conn : connectname_[ind_base])
            {
              linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/" + itr_conn + "/";
              if (cg_link_write(itr_conn.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_lagrangian : cg_link_write connectivity !" << finl, TRUST_CGNS_ERROR();
            }

          std::string solname = "FlowSolution_itr_" + std::to_string(i);
          linkpath = "/" + nom_dom.getString() + "/" + nom_dom.getString() + "/" + solname + "/";
          if (cg_link_write(solname.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_lagrangian : cg_link_write FlowSolution !" << finl, TRUST_CGNS_ERROR();
        }
    }

  cgns_close_grid_or_solution_link_file(-123., TYPE_LINK_CGNS::FINAL_LINK, true);
}

#endif /* HAS_CGNS */
