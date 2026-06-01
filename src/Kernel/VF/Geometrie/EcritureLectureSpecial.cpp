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

#include <EcritureLectureSpecial.h>

#include <Champ_Fonc_base.h>
#include <Domaine_VF.h>
#include <Domaine.h>
#include <MD_Vector_tools.h>
#include <Octree_Double.h>
#include <MD_Vector_composite.h>
#include <MD_Vector_seq.h>
#include <TRUSTTab_parts.h>
#include <TRUST_2_PDI.h>

int EcritureLectureSpecial::mode_ecr=-1;
int EcritureLectureSpecial::mode_lec=0;
int EcritureLectureSpecial::Active=0;

#ifdef MPI_
Nom EcritureLectureSpecial::Input="LecFicDiffuseBin";    // "EFichierBin" was the past (<=1.6.9) and "LecFicPartageMPIIO" is may be the future
Nom EcritureLectureSpecial::Output="EcrFicPartageMPIIO"; // "EcrFicPartageBin" was the past (<=1.6.9) and "EcrFicPartageMPIIO" is after 1.7.0
#else
Nom EcritureLectureSpecial::Input="EFichierBin";
Nom EcritureLectureSpecial::Output="EcrFicPartageBin";
#endif

Implemente_instanciable(EcritureLectureSpecial,"EcritureLectureSpecial",Interprete);

const DoubleTab& get_ref_coordinates_items(const Domaine_VF& zvf, const MD_Vector& md)
{
  // I would have taken xv but it has no parallel structure !!!
  if (md == zvf.face_sommets().get_md_vector())
    return zvf.xv(); // Descriptor for faces
  else if (md == zvf.md_vector_faces_bord())
    return zvf.xv_bord();
  // I would have taken xp but it has no parallel structure !!!!!!$&
  else if (md == zvf.domaine().les_elems().get_md_vector())
    return zvf.xp(); // Descriptor for elements
  else if (md == zvf.xa().get_md_vector())
    return zvf.xa(); // Descriptor for edges
  else if (md == zvf.domaine().les_sommets().get_md_vector())
    return zvf.domaine().les_sommets();
  else
    {
      Cerr << "Error in get_ref_coordinates_items\n"
           << " descriptor not found in this domaine" << finl;
      Process::exit();
    }
  return zvf.xv(); // for compiler
}
Sortie& EcritureLectureSpecial::printOn(Sortie& os) const
{
  return os;
}

Entree& EcritureLectureSpecial::readOn(Entree& is)
{
  return is;
}
// XD ecriturelecturespecial interprete ecriturelecturespecial INHERITS_BRACE Class to write or not to write a .xyz file
// XD_CONT on the disk at the end of the calculation.
// XD attr type chaine type REQ If set to 0 (the default), no xyz file is created. If set to 1, the .xyz file is written
// XD_CONT at the end of the computation.
Entree& EcritureLectureSpecial::interpreter(Entree& is)
{
  Nom option;
  is >> option;
  if (option=="EFichierBin" || option=="LecFicDiffuseBin" || option=="LecFicPartageMPIIO")
    {
      Cerr << "ERROR: EcritureLectureSpecial: - option '" << option << "' is obsolete! You may only specify 0 or 1." <<finl;
      Process::exit();
    }
  else if (option=="EcrFicPartageMPIIO" || option=="EcrFicPartageBin" || option=="EcrFicPartage")
    {
      Cerr << "ERROR: EcritureLectureSpecial: - option '" << option << "' is obsolete! You may only specify 0 or 1." <<finl;
      Process::exit();
    }
  else if (option=="0")
    Active=0;
  else if (option=="1")
    Active=1;
  else
    {
      Cerr << "ERROR: \"EcritureLectureSpecial " << option << "\" is not recognized." << finl;
      Process::exit();
    }
  Cerr << "EcritureLectureSpecial::Active set to " << Active << finl;
  return is;
}

/*! @brief Indicates whether the special format was requested in active reading by xyz restart.
 *
 * ...
 *
 */
int EcritureLectureSpecial::is_lecture_special()
{
  return mode_lec;
}

/*! @brief Indicates whether the special format was requested in active writing by xyz save.
 *
 * ...
 *   if the write mode is special, i.e. if the save format is xyz:
 *     then special=1, a_faire=je_suis_maitre
 *     otherwise special=0, a_faire=1
 *
 *
 */
int EcritureLectureSpecial::is_ecriture_special(int& special,int& a_faire)
{
  // with PDI, no one is writing, the library handles the IO
  if(TRUST_2_PDI::is_PDI_checkpoint())
    {
      a_faire = 0;
      special = 0;
    }
  else
    {
      special=0;
      a_faire=1;
      assert(mode_ecr>=0); // mode_ecr is not set
      if (mode_ecr)
        {
          special=1;
          a_faire=Process::je_suis_maitre();
        }
    }
  return mode_ecr;
}

/*! @brief Simple call to EcritureLectureSpecial::ecriture_special (const Domaine_VF& zvf,Sortie& fich,int nbval,const DoubleTab& val)
 *
 *     after retrieving the val array
 *
 */
int EcritureLectureSpecial::ecriture_special(const Champ_base& ch, Sortie& fich)
{
  const Domaine_VF& zvf = ref_cast(Domaine_VF, ch.domaine_dis_base());
  const DoubleTab& val = ch.valeurs();
  return ecriture_special(zvf, fich, val);
}

int ecrit(Sortie& fich, const ArrOfBit& items_to_write, const DoubleTab& pos, const DoubleTab& val)
{
  const int nb_dim = val.nb_dim();
  const int nb_comp = (nb_dim == 2) ? val.dimension(1) : 1;
  const int nb_val = items_to_write.size_array();
  const int dim = pos.dimension(1);

  if (EcritureLectureSpecial::get_Output().finit_par("MPIIO"))
    {
      // No bufferisation needed for Parallel IO
      int jmax = (dim + nb_comp) * nb_val;
      ArrOfDouble tmp(jmax);
      int j=0;
      for (int p = 0; p < nb_val; p++)
        {
          // We only write the real non-shared items so as to have a .xyz file
          // of the same size regardless of the partitioning, and above all to be able to re-read it
          // regardless of the partitioning and the supports
          if (items_to_write[p])
            {
              for (int k = 0; k < dim; k++)
                tmp[j++] = pos(p, k);
              if (nb_dim == 1)
                tmp[j++] = val(p);
              else
                for (int k = 0; k < nb_comp; k++)
                  tmp[j++] = val(p, k);
            }
        }
      fich.put(tmp.addr(), j, dim + nb_comp /* nb colonnes en ascii */);
    }
  else
    {
      // Bufferisation needed for EcrFicPartage
      int jmax = (dim + nb_comp) * 128;
      ArrOfDouble tmp(jmax);
      int j = 0;
      for (int p = 0; p < nb_val; p++)
        {
          // We only write the real non-shared items so as to have a .xyz file
          // of the same size regardless of the partitioning, and above all to be able to re-read it
          // regardless of the partitioning and the supports
          if (items_to_write[p])
            {
              for (int k = 0; k < dim; k++)
                tmp[j++] = pos(p, k);
              if (nb_dim == 1)
                tmp[j++] = val(p);
              else
                for (int k = 0; k < nb_comp; k++)
                  tmp[j++] = val(p, k);
              if (j == jmax)
                {
                  fich.put(tmp.addr(), j, dim + nb_comp /* nb colonnes en ascii */);
                  // We flush regularly in sequential mode because on some very large meshes
                  // a stack overflow is possible...
                  if (Process::is_sequential()) fich.syncfile();
                  j = 0;
                }
            }
        }
      if (j)
        fich.put(tmp.addr(), j, dim + nb_comp /* nb colonnes en ascii */);
    }
  return 8 * (dim + nb_comp) * nb_val; // Bytes
}

/*! @brief "Inner" part of the write, called by the method below.
 *
 * Recursive method, if the array to write has an MD_Vector_composite descriptor
 *
 */
static int ecriture_special_part2(const Domaine_VF& zvf, Sortie& fich, const DoubleTab& val)
{
  const MD_Vector& md = val.get_md_vector();
  int bytes = 0;
  if (sub_type(MD_Vector_composite, md.valeur()))
    {
      // p1bubble fields and others: recursive call for the different sub-arrays:
      ConstDoubleTab_parts parts(val);
      int n = zvf.que_suis_je() == "Domaine_PolyMAC_MPFA" ? 1 : parts.size();//skip auxiliary variables of Champ_{P0,Face}_PolyMAC_MPFA
      for (int i = 0; i < n; i++)
        bytes += ecriture_special_part2(zvf, fich, parts[i]);
    }
  else if (sub_type(MD_Vector_std, md.valeur()) || sub_type(MD_Vector_seq, md.valeur()))
    {
      ArrOfBit items_to_write;
      md->get_sequential_items_flags(items_to_write);
      const DoubleTab& coords = get_ref_coordinates_items(zvf, md);
      bytes += ecrit(fich, items_to_write, coords, val);
      fich.syncfile();
    }
  else
    {
      Cerr << "EcritureLectureSpecial::ecriture_special_part: Error, unknown Metadata vector type : "
           << md->que_suis_je() << finl;
      Process::exit();
    }
  return bytes;
}

/*! @brief Encoding of the write of positions and values of val
 *
 */
int EcritureLectureSpecial::ecriture_special(const Domaine_VF& zvf, Sortie& fich, const DoubleTab& val)
{
  const MD_Vector& md = val.get_md_vector();
  if (!md)
    {
      Cerr << "EcritureLectureSpecial::ecriture_special: error, cannot save an array with no metadata" << finl;
      Process::exit();
    }
  const trustIdType nb_items_seq = md->nb_items_seq_tot();
  if (nb_items_seq == 0)
    return 0;

  const int nb_dim = val.nb_dim();
  const int nb_comp = (nb_dim == 2) ? val.dimension(1) : 1;
  const int dim = Objet_U::dimension;
  const trustIdType n = nb_items_seq * (nb_comp + dim);


  if (Process::je_suis_maitre())
    {
      fich << (int)1 << finl;
      fich	 << n << finl ;
      fich << (int)1 << finl;
      fich << n << finl ;
      fich << n <<finl;
    }

  int bytes = ecriture_special_part2(zvf, fich, val);

  if (Process::je_suis_maitre())
    {
      fich << n << finl;
      fich << (int)0 << finl;
      fich << (int)0 << finl;
      fich << (int)0 << finl;
      fich << (int)1 << finl;
      fich << (int)0 << finl;
      fich << n << finl;
      fich << finl;
      fich << (int)1 << finl;
      fich << (int)0 << finl;
      fich << (int)0 <<finl;
    }
  fich.syncfile();
  return bytes;
}

/*! @brief Simple call to EcritureLectureSpecial::lecture_special (const Domaine_VF& zvf,Entree& fich,int nbval, DoubleTab& val )
 *
 */
void EcritureLectureSpecial::lecture_special(Champ_base& ch, Entree& fich)
{
  const Domaine_VF& zvf=ref_cast(Domaine_VF,ch.domaine_dis_base());
  DoubleTab& val = ch.valeurs();
  lecture_special(zvf, fich, val);
}


/*! @brief Reciprocal of the ecrit(.
 *
 * ..) method, reads only the sequential items (i.e. not the shared items received from another processor).
 *    We verify at the end that we have read exactly the expected number of items; if any are missing
 *    it means epsilon is wrong (or the mesh has changed...).
 *  Return value: total number of sequential items read (across all procs)
 *
 */
static trustIdType lire_special(Entree& fich, const DoubleTab& coords, DoubleTab& val, const double epsilon)
{
  const int dim = coords.dimension(1);
  const int nb_dim = val.nb_dim();
  const int nb_comp = (nb_dim == 1) ? 1 : val.dimension(1);

  const MD_Vector& md_vect = val.get_md_vector();
  // Initially, 1 if the item is to be read, 0 if it is read by another processor.
  // Once the item is read, the flag is set to 2.
  ArrOfInt items_to_read;
  const int n_to_read = md_vect->get_sequential_items_flags(items_to_read);
  Octree_Double octree;
  // Build an octree with "thick" nodes (epsilon size)
  octree.build_nodes(coords, 0 /* do not include virtual elements */, epsilon);
  const ArrOfInt& floor_elements = octree.floor_elements();

  // The file contains this number of lines for this part of the array (total number of sequential items)
  const trustIdType ntot = Process::mp_sum(n_to_read);

  // We read from the file in blocks of buflines_max because there is
  //  a network broadcast at each comm:
  const int buflines_max = 2048; // not too large, so as to have several blocks in test cases
  DoubleTab buffer(buflines_max, dim + nb_comp);
  int bufptr = buflines_max;
  ArrOfInt items;


  double max_epsilon_needed = epsilon;
  // How many times did we find multiple candidates within epsilon?
  int error_too_many_matches = 0;
  // How many times did we encounter the same vertex to read multiple times?
  int error_duplicate_read = 0;
  // How many items have we read?
  int count_items_read = 0;

  // Loop over sequential items in the file:
  trustIdType pourcent=0;
  for (trustIdType i = 0; i < ntot; i++)
    {
      trustIdType tmp=(i*10)/(ntot-1);
      if (tmp>pourcent || i==0)
        {
          pourcent=tmp;
          Cerr<<"\r"<<pourcent*10<<"% of data has been found."<<flush;
        }
      if (bufptr == buflines_max)
        {
          bufptr = 0;
          trustIdType n = std::min(static_cast<trustIdType>(buflines_max), ntot - i);
          n *= (dim + nb_comp);
          assert(n <= static_cast<trustIdType>(buffer.size_array()));
          fich.get(buffer.addr(), static_cast<int>(n)); // n is small enough
        }
      const double x = buffer(bufptr, 0);
      const double y = buffer(bufptr, 1);
      const double z = (dim == 3) ? buffer(bufptr, 2) : 0.;
      // Search for items potentially corresponding to the point (x,y,z)
      int index = -1;
      int nb_items_proches = octree.search_elements(x, y, z, index);
      if (nb_items_proches > 0)
        {
          items.resize_array(nb_items_proches, RESIZE_OPTIONS::NOCOPY_NOINIT);
          // See doc of Octree_Double::search_elements: copy the indices of nearby items into items:
          for (int j = 0; j < nb_items_proches; j++)
            items[j] = floor_elements[index++];
          // Reduce the list to keep only items within epsilon
          const int item_le_plus_proche = octree.search_nodes_close_to(x, y, z, coords, items, epsilon);
          nb_items_proches = items.size_array();
          if (nb_items_proches == 1)
            {
              const int flag = items_to_read[item_le_plus_proche];
              if (flag == 1)
                {
                  // Ok, we need to read this value
                  items_to_read[item_le_plus_proche] = 2;
                  count_items_read++;
                  if (nb_dim == 1)
                    {
                      val(item_le_plus_proche) = buffer(bufptr, dim);
                    }
                  else
                    {
                      for (int j = 0; j < nb_comp; j++)
                        val(item_le_plus_proche, j) = buffer(bufptr, dim + j);
                    }
                }
              else if (flag == 0)
                {
                  // This item does not belong to me, do not read it
                }
              else
                {
                  // Error, this item has already been read!!! epsilon is too large (or a save error???)
                  error_duplicate_read++;
                }
            }
          else if (nb_items_proches == 0)
            {
              // ok, the vertex is on another processor (or epsilon too small??)
            }
          else
            {
              // Error: epsilon is too large, we have multiple candidates within epsilon
              // Compute the distance to the second closest to display an error message at the end:
              for (int ii = 0; ii < nb_items_proches; ii++)
                {
                  const int i_coord = items[ii];
                  if (i_coord == item_le_plus_proche)
                    continue; // that one is probably the correct match, we need an epsilon smaller than this distance...
                  double xx = 0;
                  for (int j = 0; j < dim; j++)
                    {
                      double yy = coords(i_coord, j) - buffer(bufptr, j);
                      xx += yy * yy;
                    }
                  // We propose to set epsilon at most equal to 1/10 of the distance to the second closest point:
                  xx = sqrt(xx) * 0.1;
                  if (max_epsilon_needed > xx)
                    max_epsilon_needed = xx;
                }
              error_too_many_matches++;
            }
        }
      bufptr++;
    }
  Cerr << finl;
  // Errors?
  int err = (count_items_read != n_to_read) || (error_too_many_matches > 0) || (error_duplicate_read > 0);
  err = static_cast<int>(Process::mp_sum(err));  // sum of 0 and 1, always 'int'
  if (err)
    {
      error_too_many_matches = static_cast<int>(Process::mp_sum(error_too_many_matches));
      error_duplicate_read = static_cast<int>(Process::mp_sum(error_duplicate_read));
      max_epsilon_needed = Process::mp_min(max_epsilon_needed);
      if (Process::je_suis_maitre())
        {
          if (error_too_many_matches)
            {
              Cerr << "Error in EcritureLectureSpecial: error_too_many_matches = " << error_too_many_matches
                   << ", epsilon is too large. Suggested value: " << max_epsilon_needed << finl;
              if (max_epsilon_needed==0)
                {
                  Cerr << "It could be because your mesh has two boundaries which match exactly." << finl;
                  Cerr << "It is possible to do calculation with this property but xyz restart process" << finl;
                  Cerr << "is impossible because it can't detect the differences between faces of the two boundaries..." << finl;
                  Cerr << "Try to do a classic restart with .sauv files." << finl;
                }
            }
          else if (error_duplicate_read)
            {
              Cerr << "Error in EcritureLectureSpecial: error_duplicate_read = " << error_duplicate_read
                   << ", probably epsilon too large. " << finl;
            }
          else
            {
              Cerr << "Error in EcritureLectureSpecial: Some items were not found: epsilon too small (or the mesh has changed?)" << finl;
            }
        }
      Process::barrier();
      Process::exit();
    }
  return ntot;
}

// Return value: total number of sequential items read (across all procs)
static trustIdType lecture_special_part2(const Domaine_VF& zvf, Entree& fich, DoubleTab& val)
{
  const MD_Vector& md = val.get_md_vector();

  trustIdType ntot = 0;
  if (sub_type(MD_Vector_composite, md.valeur()))
    {
      // p1bubble fields and others: recursive call for the different sub-arrays:
      DoubleTab_parts parts(val);
      const int n = parts.size();
      for (int i = 0; i < n; i++)
        ntot += lecture_special_part2(zvf, fich, parts[i]);
    }
  else if (sub_type(MD_Vector_std, md.valeur()) || sub_type(MD_Vector_seq, md.valeur()))
    {
      const DoubleTab& coords = get_ref_coordinates_items(zvf, md);
      const double epsilon = zvf.domaine().epsilon();
      ntot += lire_special(fich, coords, val, epsilon);
    }
  else
    {
      Cerr << "EcritureLectureSpecial::lecture_special_part2: Error, unknown Metadata vector type : "
           << md->que_suis_je() << finl;
      Process::exit();
    }
  return ntot;
}

/*! @brief Encoding of the re-reading of a field from a special positions/values file
 *
 */
void EcritureLectureSpecial::lecture_special(const Domaine_VF& zvf, Entree& fich, DoubleTab& val)
{

  const MD_Vector& md_vect = val.get_md_vector();
  if (!md_vect)
    {
      Cerr << "EcritureLectureSpecial::ecriture_special: error, cannot save an array with no metadata" << finl;
      Process::exit();
    }
  const trustIdType nb_items_seq = md_vect->nb_items_seq_tot();
  if (nb_items_seq == 0)
    return;

  trustIdType bidon;
  fich >> bidon >> bidon >> bidon >> bidon >> bidon;

  trustIdType ntot = lecture_special_part2(zvf, fich, val);
  if (ntot != nb_items_seq)
    {
      Cerr << "Internal error in EcritureLectureSpecial::lecture_special" << finl;
      exit();
    }

  fich >> bidon >> bidon >> bidon >> bidon >> bidon >> bidon >> bidon >> bidon >> bidon >> bidon;

  // Update the virtual parts of the val array
  val.echange_espace_virtuel();
}

/*! @brief Returns the write mode in use (so it can be modified).
 *
 *   This method is static.
 *
 */
Nom& EcritureLectureSpecial::get_Output()
{
  static Nom option=Output;

  // disable MPIIO in sequential mode
  if (Output=="EcrFicPartageMPIIO" && Process::is_sequential()) option="EcrFicPartageBin";

  // disable MPIIO if TRUST_DISABLE_MPIIO=1
  char* theValue = getenv("TRUST_DISABLE_MPIIO");
  if (theValue != nullptr)
    {
      if (option=="EcrFicPartageMPIIO" && strcmp(theValue,"1")==0) option="EcrFicPartageBin";
    }

  return option;
}
