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

#include <Ecrire_CGNS.h>
#include <unordered_set>
#include <Domaine.h>
#include <unistd.h>

#ifdef HAS_CGNS

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
      for (double t : time_post_)
        {
          std::string zname = nom_dom.getString();
          if (!first_zone)
            zname += cgns_helper_.convert_double_to_string(t);
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
            zn += cgns_helper_.convert_double_to_string(time_post_[i]);

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

          std::string solname = "FlowSolution" + cgns_helper_.convert_double_to_string(time_post_[i]) + "_" + LOC;
          linkpath = "/" + nom_dom.getString() + "/" + nom_dom.getString() + "/" + solname + "/";
          if (cg_link_write(solname.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_lagrangian : cg_link_write FlowSolution !" << finl, TRUST_CGNS_ERROR();
        }
    }

  cgns_close_grid_or_solution_link_file(-123., TYPE_LINK_CGNS::FINAL_LINK, true);
}

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

          for (auto &itr_conn : connectname_[ind_base])
            {
              std::string linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/" + itr_conn + "/";

              if (cg_link_write(itr_conn.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::link_multi_loc_support_pb_deformable : cg_link_write !" << finl, TRUST_CGNS_ERROR();
            }
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

void Ecrire_CGNS::cgns_write_final_link_file_comm_group_pb_deformable()
{
  if (!Process::me())
    {
      std::string fn = baseFile_name_ + ".cgns";

      unlink(fn.c_str());
      cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);

      const int nb_grps = static_cast<int>(unique_vec_proc_maitre_local_comm_.size());
      std::vector<int> zoneId_tmp(nb_grps, -123);

      for (auto &itr : doms_written_)
        {
          bool has_field = false;
          std::string LOC = "rien";
          TRUST_2_CGNS::init_has_field_and_loc_iters(itr, fld_loc_map_, has_field, LOC);

          int index_glob = -123, ind_base = -123;
          TRUST_2_CGNS::init_glob_base_domain_idx(doms_written_, itr, true /* has_field */, LOC, index_glob, ind_base);

          if (cg_base_write(fileId_, itr.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_base_write !" << finl, TRUST_CGNS_ERROR();

          for (int gid = 0; gid < nb_grps; gid++)
            {
              int proc_grp = unique_vec_proc_maitre_local_comm_[gid];
              std::string zone_name = Nom("Zone").nom_me(proc_grp).getString();

              const cgsize_t isize[3] = { sizeId_som_local_comm_[ind_base][gid], sizeId_elem_local_comm_[ind_base][gid], 0 };

              std::string file_group_id = Nom(baseFile_name_).nom_me(proc_grp).getString();
              TRUST_2_CGNS::remove_slash_linkfile(file_group_id);

              cgns_helper_.cgns_write_zone_and_deformable_links(true /* write zone */, has_field, fileId_, baseId_[index_glob], zone_name, isize, zoneId_tmp[gid], gid + 1,
                                                                file_group_id, baseZone_name_[ind_base], baseZone_name_[ind_base], connectname_[ind_base], itr, LOC, time_post_,
                                                                "Ecrire_CGNS::cgns_write_final_link_file_comm_group_pb_deformable");
            }

          cgns_helper_.cgns_write_iters_deformable<TYPE_ECRITURE_CGNS::SEQ>(true, has_field, nb_grps /* nb_zones_to_write */, fileId_, baseId_[index_glob], ind_base,
                                                                            zoneId_tmp, LOC, solname_som_, solname_elem_, solname_faces_, grid_name_, time_post_);

        }
      cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);
    }
}

void Ecrire_CGNS::cgns_write_final_link_file_pb_deformable()
{
  if (Process::is_parallel() && Option_CGNS::LINKED_FILES_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
    {
#ifdef MPI_
      if (vec_proc_maitre_local_comm_.empty())
        gather_local_sizeId_for_comm_group();

      cgns_write_final_link_file_comm_group_pb_deformable();
#endif
      return;
    }

  if (!Process::me()) // seul le proc 0 ecrit le fichier link
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
          TRUST_2_CGNS::init_glob_base_domain_idx(doms_written_, itr, true /* has_field */, LOC, index_glob, ind_base);

          if (cg_base_write(fileId_, itr.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_base_write !" << finl, TRUST_CGNS_ERROR();

          const cgsize_t isize[3] = { sizeId_[ind_base][0] , sizeId_[ind_base][1] , 0 };

          cgns_helper_.cgns_write_zone_and_deformable_links(true /* write zone */, has_field, fileId_, baseId_[index_glob], itr.getString(), isize, zoneId_[index_glob], 1,
                                                            base_link_file, baseZone_name_[ind_base], baseZone_name_[ind_base], connectname_[ind_base], itr, LOC, time_post_,
                                                            "Ecrire_CGNS::cgns_write_final_link_file_pb_deformable");

          cgns_helper_.cgns_write_iters_deformable<TYPE_ECRITURE_CGNS::SEQ>(true /* deformable */, has_field, 1 /* 1 zone per base */, fileId_, baseId_[index_glob], index_glob /* 1st Zone */,
                                                                            zoneId_, LOC, solname_som_, solname_elem_, solname_faces_, grid_name_, time_post_);
        }
      cgns_close_grid_or_solution_link_file(-123., TYPE_LINK_CGNS::FINAL_LINK, true);
    }
}

/*
 * ******************** *
 * VERSION SEQUENTIELLE *
 * ******************** *
 */
void Ecrire_CGNS::cgns_write_domaine_deformable_seq(const Domaine * domaine,const Nom& nom_dom, const DoubleTab& les_som, const IntTab& les_elem, const Motcle& type_elem)
{
  const int ind = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
  TRUST_2_CGNS& TRUST2CGNS = T2CGNS_[ind];
  TRUST2CGNS.associer_domaine_TRUST(domaine, domaine_dis_.non_nul() ? &(domaine_dis_.valeur()) : nullptr, les_som, les_elem, postraiter_domaine_);

  CGNS_TYPE cgns_type_elem = TRUST2CGNS.convert_elem_type(type_elem);
  const bool is_polyedre = (type_elem == "POLYEDRE" || type_elem == "PRISME" || type_elem == "PRISME_HEXAG");
  const int icelldim = TRUST2CGNS.topo_dim_from_elem(cgns_type_elem, is_polyedre); // avant ca : icelldim = les_som.dimension(1)
  const int iphysdim = Objet_U::dimension, nb_som = les_som.dimension(0), nb_elem = les_elem.dimension(0);

  std::vector<double> xCoords, yCoords, zCoords;
  TRUST2CGNS.fill_coords(xCoords, yCoords, zCoords);

  const cgsize_t isize[3] = { (cgsize_t)nb_som, (cgsize_t)nb_elem, 0 }; /* 0 => boundary vertex size (zero if elements not sorted) */

  int coordsId = -1;
  char basename[CGNS_STR_SIZE];
  strcpy(basename, nom_dom.getChar()); // dom name

  if (!Option_CGNS::USE_LINKS && !is_lagrangian_)
    {
      if (nb_elem)
        {
          int G = -1;
          if (cg_grid_write(fileId_, baseId_[ind], zoneId_[ind], grid_name_loc_.c_str(), &G) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_seq : cg_grid_write !" << finl, TRUST_CGNS_ERROR();

          if (cg_goto(fileId_, baseId_[ind], "Zone_t",zoneId_[ind], "GridCoordinates_t",  G, "end") != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_seq : cg_goto Zone_t !" << finl, TRUST_CGNS_ERROR();

          const cgsize_t dims[1] = { (cgsize_t)nb_som };

          if (cg_array_write("CoordinateX", CGNS_ENUMV(RealDouble), 1, dims, xCoords.data()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_seq : cg_array_write CoordinateX !" << finl, TRUST_CGNS_ERROR();

          if (cg_array_write("CoordinateY", CGNS_ENUMV(RealDouble), 1, dims, yCoords.data()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_seq : cg_array_write CoordinateY !" << finl, TRUST_CGNS_ERROR();

          if (Objet_U::dimension > 2)
            if (cg_array_write("CoordinateZ", CGNS_ENUMV(RealDouble), 1, dims, zCoords.data()) != CG_OK)
              Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_seq : cg_array_write CoordinateZ !" << finl, TRUST_CGNS_ERROR();
        }
    }
  else
    {
      if (cg_base_write(fileId_, basename, icelldim, iphysdim, &baseId_[ind]) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_domaine_seq : cg_base_write !" << finl, TRUST_CGNS_ERROR();

      if (nb_elem)
        {
          /* Create zone & grid coords */
          cgns_helper_.cgns_write_zone_grid_coord<TYPE_ECRITURE_CGNS::SEQ>(icelldim, fileId_, baseId_[ind], basename /* Dom name */, isize,
                                                                           zoneId_[ind], xCoords, yCoords, zCoords, coordsId, coordsId, coordsId);

          if (is_lagrangian_)
            {
              sizeId_.push_back( { (cgsize_t)nb_som, (cgsize_t)nb_elem } ); // XXX required for links later !

              /* Set element connectivity : we rewrite since topology can change !! */
              int sectionId;
              cgsize_t start = 1, end;

              if (cgns_type_elem == CGNS_ENUMV(NGON_n)) // cas polyedre
                throw std::runtime_error("Ecrire_CGNS::cgns_write_domaine_deformable_seq => You should not be here !!! ");
              else
                {
                  std::vector<cgsize_t> elems;
                  const int nsom = TRUST2CGNS.convert_connectivity(cgns_type_elem, elems);

                  end = start + static_cast<cgsize_t>(elems.size()) / nsom - 1;

                  if (cg_section_write(fileId_, baseId_[ind], zoneId_[ind], "Elem", cgns_type_elem, start, end, 0, elems.data(), &sectionId) != CG_OK)
                    Cerr << "Error Ecrire_CGNS::cgns_write_domaine_seq : cg_section_write !" << finl, TRUST_CGNS_ERROR();
                }
            }
          else
            {
              /* Set element connectivity : by links */
              std::string linkfile = baseFile_name_ + ".solution." + cgns_helper_.convert_double_to_string(time_post_[0]) + ".cgns";
              TRUST_2_CGNS::remove_slash_linkfile(linkfile);

              if (cg_goto(fileId_, baseId_[ind], "Zone_t", zoneId_[ind], "end") != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_seq : cg_goto Zone_t !" << finl, TRUST_CGNS_ERROR();

              for (auto &itr_conn : connectname_[ind])
                {
                  const std::string linkpath = "/" + baseZone_name_[ind] + "/" + baseZone_name_[ind] + "/" + itr_conn + "/";

                  if (cg_link_write(itr_conn.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                    Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_seq : cg_link_write connectivity !" << finl, TRUST_CGNS_ERROR();
                }
            }
        }
    }
}

/*
 * ************************* *
 * VERSION PARALLELE IN ZONE *
 * ************************* *
 */
void Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone(const Domaine * domaine,const Nom& nom_dom, const DoubleTab& les_som, const IntTab& les_elem, const Motcle& type_elem)
{
#ifdef MPI_
  const int ind = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
  TRUST_2_CGNS& TRUST2CGNS = T2CGNS_[ind];
  TRUST2CGNS.associer_domaine_TRUST(domaine, domaine_dis_.non_nul() ? &(domaine_dis_.valeur()) : nullptr, les_som, les_elem, postraiter_domaine_);

  CGNS_TYPE cgns_type_elem = TRUST2CGNS.convert_elem_type(type_elem);
  const bool is_polyedre = (type_elem == "POLYEDRE" || type_elem == "PRISME" || type_elem == "PRISME_HEXAG");
  const int icelldim = TRUST2CGNS.topo_dim_from_elem(cgns_type_elem, is_polyedre); // avant ca : icelldim = les_som.dimension(1)
  const int nb_elem = les_elem.dimension(0), iphysdim = Objet_U::dimension;

  std::vector<double> xCoords, yCoords, zCoords;
  TRUST2CGNS.fill_coords(xCoords, yCoords, zCoords);

  char basename[CGNS_STR_SIZE];
  strcpy(basename, nom_dom.getChar()); // dom name

  TRUST2CGNS.fill_global_infos(); // XXX utile car info change en //

  if (cgns_type_elem == CGNS_ENUMV(NGON_n)) /*cas polygone/polyedre */
    TRUST2CGNS.fill_global_infos_poly(is_polyedre);

  const int ns_tot = TRUST2CGNS.get_ns_tot(), ne_tot = TRUST2CGNS.get_ne_tot();
  const bool enter_group_comm = Option_CGNS::LINKED_FILES_PER_COMM_GROUP && PE_Groups::has_user_defined_group() && !postraiter_domaine_;
  const int proc_me = enter_group_comm ? TRUST2CGNS.get_proc_me_local_comm() : Process::me();

  int coordsIdx = -123, coordsIdy = -123, coordsIdz = -123;

  if (!Option_CGNS::USE_LINKS && !is_lagrangian_)
    {
      int G = -1;
      if (cg_grid_write(fileId_, baseId_[ind], zoneId_[ind], grid_name_loc_.c_str(), &G) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cg_grid_write !" << finl, TRUST_CGNS_ERROR();

      if (cg_goto(fileId_, baseId_[ind], "Zone_t",zoneId_[ind], "GridCoordinates_t",  G, "end") != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cg_goto Zone_t !" << finl, TRUST_CGNS_ERROR();

      const cgsize_t dims[1] = { ns_tot };

      if (cgp_array_write("CoordinateX", CGNS_ENUMV(RealDouble), 1, dims, &coordsIdx) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cgp_array_write CoordinateX !" << finl, TRUST_CGNS_ERROR();

      if (cgp_array_write("CoordinateY", CGNS_ENUMV(RealDouble), 1, dims, &coordsIdy) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cgp_array_write CoordinateY !" << finl, TRUST_CGNS_ERROR();

      if (Objet_U::dimension > 2)
        if (cgp_array_write("CoordinateZ", CGNS_ENUMV(RealDouble), 1, dims, &coordsIdz) != CG_OK)
          Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cgp_array_write CoordinateZ !" << finl, TRUST_CGNS_ERROR();

      if (nb_elem > 0) // seulement si le proc a qlq chose a ecrire
        {
          const std::vector<int>& incr_max_som = TRUST2CGNS.get_global_incr_max_som(),
                                  &incr_min_som = TRUST2CGNS.get_global_incr_min_som();

          const cgsize_t min = incr_min_som[proc_me], max = incr_max_som[proc_me];
          assert (min < max);

          if (cgp_array_write_data(coordsIdx, &min, &max, xCoords.data()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cgp_array_write_data CoordinateX !" << finl, TRUST_CGNS_ERROR();

          if (cgp_array_write_data(coordsIdy, &min, &max, yCoords.data()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cgp_array_write_data CoordinateY !" << finl, TRUST_CGNS_ERROR();

          if (Objet_U::dimension > 2)
            if (cgp_array_write_data(coordsIdz, &min, &max, zCoords.data()) != CG_OK)
              Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cgp_array_write_data CoordinateZ !" << finl, TRUST_CGNS_ERROR();
        }
    }
  else
    {
      if (cg_base_write(fileId_, basename, icelldim, iphysdim, &baseId_[ind]) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cg_base_write !" << finl, TRUST_CGNS_ERROR();

      cgsize_t isize[3];
      isize[0] = (ns_tot == 0 && enter_group_comm) ? 1 : ns_tot; // si ns_tot = 0, on va juste creer une zone vide
      isize[1] = (ne_tot == 0 && enter_group_comm) ? 1 : ne_tot; // si ne_tot = 0, on va juste creer une zone vide
      isize[2] = 0; /* boundary vertex size (zero if elements not sorted) */

      cgns_helper_.cgns_write_zone_grid_coord<TYPE_ECRITURE_CGNS::PAR_IN>(icelldim, fileId_, baseId_[ind], basename /* Dom name */, isize,
                                                                          zoneId_[ind], xCoords, yCoords, zCoords, coordsIdx, coordsIdy, coordsIdz);

      int sectionId = -123;

      if (is_lagrangian_)
        {
          sizeId_.push_back( { isize[0], isize[1] } ); // XXX required for links later !

          if (ne_tot == 0 && ns_tot == 0) return; // XXX Elie Saikali : zone vide creer, rien a faire de plus ... (cas LINKED_FILES_PER_COMM_GROUP !!!)

          /* Construct the sections to host connectivity later */
          const cgsize_t start = 1, end = ne_tot;
          assert(start <= end);

          if (cgns_type_elem == CGNS_ENUMV(NGON_n)) // cas polyedre
            throw std::runtime_error("Ecrire_CGNS::cgns_write_domaine_deformable_seq => You should not be here !!! ");
          else
            {
              if (cgp_section_write(fileId_, baseId_[ind], zoneId_[ind], "Elem", cgns_type_elem, start, end, 0, &sectionId) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_in_zone : cgp_section_write !" << finl, TRUST_CGNS_ERROR();
            }
        }
      else
        {
          if (ne_tot == 0 && ns_tot == 0) return; // XXX Elie Saikali : zone vide creer, rien a faire de plus ... (cas LINKED_FILES_PER_COMM_GROUP !!!)
        }

      if (nb_elem > 0) // seulement si le proc a qlq chose a ecrire
        {
          const std::vector<int>& incr_max_som = TRUST2CGNS.get_global_incr_max_som(),
                                  &incr_min_som = TRUST2CGNS.get_global_incr_min_som();

          cgsize_t min = incr_min_som[proc_me], max = incr_max_som[proc_me];
          assert (min < max);

          /* Write grid coordinates */
          cgns_helper_.cgns_write_grid_coord_data<TYPE_ECRITURE_CGNS::PAR_IN>(icelldim, fileId_, baseId_[ind], zoneId_[ind],
                                                                              coordsIdx, coordsIdy, coordsIdz, min, max, xCoords, yCoords, zCoords);

          if (is_lagrangian_)
            {
              assert(cgns_type_elem != CGNS_ENUMV(NGON_n));

              std::vector<cgsize_t> elems;
              TRUST2CGNS.convert_connectivity(cgns_type_elem, elems);

              const std::vector<int>& incr_max_elem = TRUST2CGNS.get_global_incr_max_elem(),
                                      &incr_min_elem = TRUST2CGNS.get_global_incr_min_elem();

              min = incr_min_elem[proc_me], max = incr_max_elem[proc_me];
              assert (min <= max);

              if (cgp_elements_write_data(fileId_, baseId_[ind], zoneId_[ind], sectionId, min, max, elems.data()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_in_zone : cgp_elements_write_data !" << finl, TRUST_CGNS_ERROR();
            }
          else
            {
              /* Set element connectivity */
              std::string linkfile = (enter_group_comm ? Nom(baseFile_name_).nom_me(proc_maitre_local_comm_).getString() : baseFile_name_) +
                                     ".solution." + cgns_helper_.convert_double_to_string(time_post_[0]) + ".cgns";

              TRUST_2_CGNS::remove_slash_linkfile(linkfile);

              if (cg_goto(fileId_, baseId_[ind], "Zone_t", zoneId_[ind], "end") != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cg_goto Zone_t !" << finl, TRUST_CGNS_ERROR();

              for (auto &itr_conn : connectname_[ind])
                {
                  const std::string linkpath = "/" + baseZone_name_[ind] + "/" + baseZone_name_[ind] + "/" + itr_conn + "/";

                  if (cg_link_write(itr_conn.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                    Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cg_link_write connectivity !" << finl, TRUST_CGNS_ERROR();
                }
            }
        }
    }
#endif /*MPI_*/
}

#endif /* HAS_CGNS */
