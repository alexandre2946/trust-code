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

#include <EcrFicPartageMPIIO.h>
#include <Format_Post_Lata.h>
#include <EcrFicPartageBin.h>
#include <communications.h>
#include <Fichier_Lata.h>
#include <EFichier.h>
#include <sys/stat.h>
#include <Param.h>
#include <string> // Required with xlC for std::getline

Implemente_instanciable_sans_constructeur(Format_Post_Lata,"Format_Post_Lata",Format_Post_base);

#define _LATA_INT_TYPE_ trustIdType

/*! @brief Default constructor: format_ ASCII and options_para_ = SINGLE_FILE.
 *
 */
Format_Post_Lata::Format_Post_Lata()
{
  reset();
}

/*! @brief Resets the object to the state obtained by the default constructor.
 *
 */
void Format_Post_Lata::reset()
{
  lata_basename_ = "??";
  format_ = ASCII;
  options_para_ = SINGLE_FILE;
  status = RESET;
  restart_already_moved_ = false;
  tinit_ = -1.;
  temps_courant_ = -1.;
}

void Format_Post_Lata::resetTime(double t, const std::string dirname)
{
  temps_courant_ = -1; // not using t - this will come from outside when calling 'ecrire_temps'
}

Sortie& Format_Post_Lata::printOn(Sortie& os) const
{
  Process::exit("Format_Post_Lata::printOn : error");
  return os;
}

/*! @brief Reads post-processing parameters in "data set" format. The expected format is:
 *
 *   {
 *        nom_fichier nom                       required field
 *      [ format   ascii|binaire ]              default value: ascii
 *      [ parallel single_file|multiple_files ] default value: single_file
 *   }
 *
 */
Entree& Format_Post_Lata::readOn(Entree& is)
{
  assert(status == RESET);
  Format_Post_base::readOn(is);
  status = INITIALIZED;
  return is;
}

void Format_Post_Lata::set_param(Param& param) const
{
  Cerr << "Format_Post_Lata::set_param: Not implemented." << finl;
  Process::exit();
}

int Format_Post_Lata::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  Cerr << "Format_Post_Lata::lire_motcle_non_standard: Not implemented." << finl;
  Process::exit();
  return 0;
}

/*! @brief Returns the conventional extension for lata files: ".lata"
 *
 */
const char * Format_Post_Lata::extension_lata()
{
  static const char * ext = ".lata";
  return ext;
}

/*! @brief Returns the file name without its path: strips all characters before the last /
 *
 */
const char * Format_Post_Lata::remove_path(const char * filename)
{
  int i = (int)strlen(filename);
  while (i > 0 && filename[i-1] != '/')
    i--;
  return filename + i;
}

/*! @brief Opens the master file in ERASE mode and writes the lata file header (on the master processor only).
 *
 * void Format_Post_Lata::ecrire_entete_lata()
 *
 */
int Format_Post_Lata::ecrire_entete(const double temps_courant, const int reprise, const int est_le_premier_post)
{
  ecrire_entete_lata(lata_basename_, options_para_, format_, est_le_premier_post);
  return 1;
}

int Format_Post_Lata::completer_post(const Domaine&, const int , const Nature_du_champ&, const int , const Noms&, const Motcle&, const Nom&)
{
  return 1;
}

int Format_Post_Lata::preparer_post(const Nom& , const int , const int , const double )
{
  return 1;
}

namespace
{

template<typename TYP, typename LATA_TYP>
void fill_tmp_array(const TRUSTTab<TYP,int>& tab, int upper, int offset, bool decal_fort, trustIdType decalage_partiel, LATA_TYP *tmp);

template<typename TYP>
typename std::enable_if_t<std::is_same<TYP, int>::value || std::is_same<TYP, trustIdType>::value, void>
fill_tmp_array(const TRUSTTab<TYP,int>& tab, int upper, int offset, bool decal_fort, trustIdType decalage_partiel, _LATA_INT_TYPE_ *tmp)
{
  trustIdType decal_fort_val = decal_fort ? 1 : 0;
  const TYP *data = tab.addr();
  for (int i = 0; i < upper; i++)
    {
      // value to write (conversion to Fortran numbering if needed)
      _LATA_INT_TYPE_ x = data[i+offset];
      if (x > -1)
        x += decalage_partiel;
      else
        x += decal_fort_val;
      tmp[i] = x;
    }
}

template<>
void fill_tmp_array<double, float>(const TRUSTTab<double,int>& tab, int upper, int offset, bool , trustIdType , float *tmp)
{
  const double *data = tab.addr();
  for (int i = 0; i < upper; i++)
    tmp[i] = (float) data[i+offset];       // downcast to float
}

/** Generic method to write the block corresponding to an array of data in a LATA file
 */
template<typename TYP, typename LATA_TYP>
trustIdType write_T_tab(Fichier_Lata& fichier, bool decal_fort, trustIdType decalage_partiel, const TRUSTTab<TYP,int>& tab, int& nb_colonnes, const Format_Post_Lata::Options_Para& option)
{
  int nb_lignes = tab.dimension(0);
  int line_size = 1;
  const int nb_dim = tab.nb_dim();

  for (int i = 1; i < nb_dim; i++)
    line_size *= tab.dimension(i);

  trustIdType nb_lignes_tot = 0;

  const int tab_size = line_size * nb_lignes;
  trustIdType nb_octets = tab_size * (trustIdType) sizeof(LATA_TYP);
  switch(option)
    {
    case Format_Post_Lata::SINGLE_FILE_MPIIO:
    case Format_Post_Lata::SINGLE_FILE:
      nb_lignes_tot = Process::mp_sum(nb_lignes);
      // In parallel, all arrays must have the same number of columns (or be empty).
      nb_colonnes = Process::mp_max(line_size);
      nb_octets = nb_colonnes * nb_lignes_tot * (trustIdType) sizeof(LATA_TYP);
      assert(nb_lignes == 0 || line_size == nb_colonnes);
      break;
    case Format_Post_Lata::MULTIPLE_FILES:
      nb_lignes_tot = nb_lignes;
      nb_colonnes = line_size;
      break;
    default:
      Cerr << "Format_Post_Lata_write_tab: error nb_lignes_tot" << finl;
      Process::exit();
    }

  SFichier& sfichier = fichier.get_SFichier();

  // Start of Fortran block
  if (fichier.is_master())
    sfichier << nb_octets << finl;

  // Writing data.
  if (sub_type(EcrFicPartageMPIIO, sfichier))
    {
      // Convert everything to LATA_TYP
      LATA_TYP *tmp = new LATA_TYP[tab_size];
      fill_tmp_array(tab, tab_size, 0, decal_fort, decalage_partiel, tmp);
      sfichier.put(tmp, tab_size, line_size);
      delete[] tmp;
      // End of Fortran block
      if (fichier.is_master())
        sfichier << nb_octets << finl;
    }
  else
    {
      // Convert everything to _INT_TYPE_ in batches of N values
      // Buffer whose size is a multiple of line_size:
      const int N = 16384;
      int bufsize = (N / line_size + 1) * line_size;
      LATA_TYP *tmp = new LATA_TYP[bufsize];
      for (int i = 0; i < tab_size; i += bufsize)
        {
          int j_max = bufsize;
          if (j_max > tab_size - i)
            j_max = tab_size - i;

          // Convert the block to LATA_TYP:
          fill_tmp_array(tab, j_max, i, decal_fort, decalage_partiel, tmp);

          // Write with a newline at each row of the array
          sfichier.put(tmp, j_max, line_size);
        }
      delete[] tmp;
      fichier.syncfile();
      // End of Fortran block
      if (fichier.is_master())
        sfichier << nb_octets << finl;
      fichier.syncfile();
    }
  return nb_lignes_tot;
}

} // end anonymous namespace


/*! @brief fichier is a lata data file (not the master file). The array tab is written as-is (in binary or ASCII, on one or
 *
 *   several files in parallel).
 *   nb_colonnes is filled with the product of tab.dimension(i) for i>0.
 *
 */
trustIdType Format_Post_Lata::write_doubletab(Fichier_Lata& fichier, const DoubleTab& tab, int& nb_colonnes, const Options_Para& option)
{
  return ::write_T_tab<double, float>(fichier, 0, 0, tab, nb_colonnes, option);
}


/*! @brief Writing an integer array to the given file.
 *
 * The values written are the array values incremented by "decalage". This value is used to switch to Fortran numbering
 *   (add 1), or to switch to global numbering (add the number of elements on previous processors).
 *
 *  nb_colonnes is filled with the sum of dimension(i) for i>0.
 *  Return value: sum of the written dimension(0) values (depending on whether all processors write to the same file or not).
 */
trustIdType Format_Post_Lata::write_inttab(Fichier_Lata& fichier, bool decal_fort, trustIdType decalage_partiel, const IntTab& tab, int& nb_colonnes, const Options_Para& option)
{
  return ::write_T_tab<int, _LATA_INT_TYPE_>(fichier, decal_fort, decalage_partiel, tab, nb_colonnes, option);
}

/*! @brief Initializes the class with default parameters (ASCII format, SINGLE_FILE).
 *
 */
int Format_Post_Lata::initialize_by_default(const Nom& file_basename)
{
  assert(status == RESET);
  initialize_lata(file_basename, ASCII, SINGLE_FILE);
  return 1;
}

int Format_Post_Lata::initialize(const Nom& file_basename, const int format, const Nom& option_para)
{
  assert(status == RESET);
  // Change the LATA format (default is binary)
  format_ = BINAIRE;
  if (format == 0)
    format_ = ASCII;

  if (Motcle(option_para) == "MPI-IO")
    options_para_ = SINGLE_FILE_MPIIO;
  else if (Motcle(option_para) == "SIMPLE")
    options_para_ = SINGLE_FILE;
  else if (Motcle(option_para) == "MULTIPLE")
    options_para_ = MULTIPLE_FILES;
  else
    {
      Cerr << "The option " << option_para << " for lata format for the parallel is not correct." << finl;
      Process::exit();
    }

  initialize_lata(file_basename, format_, options_para_);

  return 1;
}

/*! @brief Initializes the class, opens the file and writes the header.
 *
 */
int Format_Post_Lata::initialize_lata(const Nom& file_basename, const Format format, const Options_Para options_para)
{
  assert(status == RESET);
  lata_basename_ = file_basename;
  format_ = format;
  options_para_ = options_para;
  status = INITIALIZED;
  return 1;
}

/*! @brief Modifying name of the post file, plus some clever management of previously saved data in case of restart.
 */
int Format_Post_Lata::modify_file_basename(const Nom file_basename, bool for_restart, const double tinit)
{
  Nom post_file;
  post_file = file_basename + extension_lata();
  // Check that the master file exists and has a correct header
  bool master_file_exists = false;
  if (Process::je_suis_maitre())
    {
      struct stat f;
      if (stat(post_file, &f) == 0)
        {
          EFichier tmp(post_file);
          Nom cle;
          tmp >> cle;
          if (cle.debute_par("LATA"))
            {
              if (tinit == -1)
                master_file_exists = true;
              else
                {
                  tmp >> cle;
                  while (!tmp.eof() && cle != "TEMPS")
                    tmp >> cle;
                  if (cle == "TEMPS")
                    {
                      double temps;
                      tmp >> temps;
                      // Check the time of the restart file
                      if (temps < tinit)
                        master_file_exists = true;
                    }
                }
            }
          tmp.close();
        }
    }

  // Saving what was there before:
  master_file_exists = (bool) Process::mp_max((int)master_file_exists);
  if (Process::je_suis_maitre() && master_file_exists && !restart_already_moved_ && for_restart)
    {
      Nom before_restart;
      before_restart = file_basename + ".before_restart" + extension_lata();
      rename(post_file, before_restart);
      Cerr << "File " << post_file << " is moved to " << before_restart << " with times<=tinit=" << tinit << finl;
      reconstruct(post_file, before_restart, tinit);
    }
  if (for_restart)
    restart_already_moved_ = true;

  Process::barrier();  // really necessary?

  lata_basename_ = file_basename;

  if (master_file_exists && for_restart)
    {
      lata_basename_ += ".after_restart";
      if (tinit == -1)
        finir_sans_reprise(post_file);
    }
  return 1;
}

// Copy the beginning of the before_restart file into post_file (up to TEMPS=tinit)
int Format_Post_Lata::reconstruct(const Nom post_file, const Nom before_restart, const double tinit)
{
  EFichier LataOld(before_restart);
  SFichier LataNew(post_file);
  LataNew.setf(ios::scientific);
  LataNew.precision(8);
  std::string line;
  Nom mot;
  double temps;
  tinit_ = tinit;
  while (!LataOld.eof())
    {
      LataOld >> mot; // read the first word of the line in .before_restart.lata
      if (mot != "FIN") // If the word is not FIN
        {
          if (mot == "TEMPS") // If the next word is TEMPS
            {
              LataOld >> temps;
              if (temps == tinit) // If the found time equals tinit then stop
                break;
              std::getline(LataOld.get_ifstream(), line); // read the rest of the line
              LataNew.get_ofstream() << mot << " " << temps << line << std::endl; // write it all to the lata
            }
          else
            {
              std::getline(LataOld.get_ifstream(), line); // read the rest of the line
              LataNew.get_ofstream() << mot << line << std::endl; // write it all to the lata
            }
        }
      else
        {
          break; // Otherwise stop
        }
    }
  LataNew.close();
  LataOld.close();
  Cerr << "File " << post_file << " is rebuilt but truncated to times<tinit=" << tinit << finl;
  return 1;
}

int Format_Post_Lata::finir_sans_reprise(const Nom file_basename)
{
  Nom post_file;
  post_file = lata_basename_ + extension_lata();
  if (Process::je_suis_maitre())
    {
      struct stat f;
      if (stat(post_file, &f) == 0)  // if file exists
        {
          EFichier Lata(file_basename);
          SFichier LataRep(post_file, ios::app);
          LataRep.setf(ios::scientific);
          LataRep.precision(8);
          std::string line;
          Nom mot;
          double temps;
          while (!Lata.eof())
            {
              Lata >> mot; // read the first word of the line in .lata
              if (mot == "TEMPS") // If the next word is not TEMPS
                {
                  Lata >> temps;
                  if (temps == tinit_) // If the found time equals tinit
                    {
                      std::getline(Lata.get_ifstream(), line);
                      LataRep.get_ofstream() << mot << " " << temps << line << std::endl; // can start writing
                      break;
                    }
                }
            }
          while (!Lata.eof())
            {
              std::getline(Lata.get_ifstream(), line); // read the rest of the .lata
              LataRep.get_ofstream() << line << std::endl; // write it all to after_restart.lata
            }
          Lata.close();
          LataRep.close();
          Cerr << "File " << post_file << " is built with times from tinit=" << tinit_ << finl;
        }
    }
  return 1;
}

/*! Write the "fileoffset=..." when needed in master lata file
 */
void Format_Post_Lata::ecrire_offset(SFichier& sfichier, long int offset_single_lata)
{
  if (un_seul_fichier_lata_)
    sfichier << " file_offset=" << offset_single_lata;
  else
    {
#ifdef INT_is_64_
      // 64b binary file start with the tag "INT64" which is 6 bytes long with the \0 terminating char.
      // In ASCII this is not needed.
      if(format_ == BINAIRE)
        sfichier << " file_offset=6";
#endif
    }
}

// E Saikali : adding this list which is useful for a coupled problem when writing to the same file
static Noms liste_single_lata_ecrit;

/*! @brief Low level routine to write a mesh into a LATA file.
 *
 * Also called directly in TrioCFD, by Postraitement_ft_lata for interface writing.
 */
void Format_Post_Lata::ecrire_domaine_low_level(const Nom& id_domaine, const DoubleTab& sommets, const IntTab& elements, const Motcle& type_element)
{
  const int dim = sommets.dimension(1);
  Motcle type_elem(type_element);

  // GF To ensure correct reading with the lata plugin
  if (type_element == "PRISME") type_elem = "PRISM6";

  trustIdType nb_som_tot, nb_elem_tot;  // Number of vertices/elements in the file

  // Build the geometry file name
  Nom basename_geom(lata_basename_), extension_geom(extension_lata());

  if (un_seul_fichier_lata_) extension_geom += "_single";
  else
    {
      extension_geom += ".";
      extension_geom += id_domaine;
      extension_geom += ".";
      char str_temps[100] = "0.0";
      if (temps_courant_ >= 0.)
        snprintf(str_temps, 100, "%.10f", temps_courant_);
      extension_geom += Nom(str_temps);
    }

  Nom nom_fichier_geom;
  trustIdType decalage_sommets = 1, decalage_elements = 1;

  {
    const bool not_in_list =  !liste_single_lata_ecrit.contient_(lata_basename_),
               should_erase = (!un_seul_fichier_lata_) ? true /* Always erase */ : (offset_elem_ < 0 && not_in_list);

    Fichier_Lata fichier_geom(basename_geom, extension_geom, should_erase ? Fichier_Lata::ERASE : Fichier_Lata::APPEND, format_, options_para_);

    // add to the list if not already in it and if un_seul_fichier_lata_ !!!
    if (not_in_list && un_seul_fichier_lata_) liste_single_lata_ecrit.add(lata_basename_); // BOOM !

    nom_fichier_geom = fichier_geom.get_filename();
    int nb_col;

    if (un_seul_fichier_lata_)
      if (fichier_geom.is_master())
        offset_som_ = fichier_geom.get_SFichier().get_ofstream().tellp();

    // Vertex coordinates
    if (axi)
      {
        DoubleTab sommets2(sommets);
        int ns = sommets2.dimension_tot(0);
        for (int s = 0; s < ns; s++)
          {
            double r = sommets(s, 0), theta = sommets(s, 1);
            sommets2(s, 0) = r * cos(theta);
            sommets2(s, 1) = r * sin(theta);
          }
        nb_som_tot = write_doubletab(fichier_geom, sommets2, nb_col, options_para_);
      }
    else
      nb_som_tot = write_doubletab(fichier_geom, sommets, nb_col, options_para_);

    assert(nb_som_tot == 0 || nb_col == dim);

    // Elements: vertex, element and other indices in the lata file start at 1:
    if (options_para_ == SINGLE_FILE || options_para_ == SINGLE_FILE_MPIIO)
      {
        // All processors write to a single file; indices must be renumbered to global numbering.
        // Processor 0 numbers its vertices from 1 to n0
        // Processor 1 numbers its vertices from n0+1 to n0+n1, etc...
        // Offset to add to indices for global numbering.
        const int nbsom = sommets.dimension(0);
        decalage_sommets += mppartial_sum(nbsom);
        const int nbelem = elements.dimension(0);
        decalage_elements += mppartial_sum(nbelem);
      }

    if (un_seul_fichier_lata_)
      {
        if (fichier_geom.is_master())
          offset_elem_ = fichier_geom.get_SFichier().get_ofstream().tellp();

        nb_elem_tot = write_inttab(fichier_geom, true, decalage_sommets, elements, nb_col, options_para_);
      }
    else
      {
        Fichier_Lata fichier_geom_elem(basename_geom, extension_geom + Nom(".elem"), Fichier_Lata::ERASE, format_, options_para_);
        nb_elem_tot = write_inttab(fichier_geom_elem, true, decalage_sommets, elements, nb_col, options_para_);
      }
  }

  {
    // Note: the file must be closed before calling ecrire_item_int, which will reopen it.
    Fichier_Lata_maitre fichier_lata(lata_basename_, extension_lata(), Fichier_Lata::APPEND, options_para_);
    SFichier& sfichier = fichier_lata.get_SFichier();

    if (fichier_lata.is_master())
      {
        sfichier << "GEOM " << id_domaine;
        sfichier << " type_elem=" << type_elem << finl;

        // SOMMETS support
        sfichier << "CHAMP SOMMETS " << remove_path(nom_fichier_geom);
        sfichier << " geometrie=" << id_domaine;
        sfichier << " size=" << nb_som_tot;
        sfichier << " composantes=" << dim;
        ecrire_offset(sfichier, offset_som_);
        sfichier << finl;

        // ELEMENTS support
        sfichier << "CHAMP ELEMENTS " << remove_path(nom_fichier_geom);
        if (!un_seul_fichier_lata_) sfichier << ".elem";
        sfichier << " geometrie=" << id_domaine;
        sfichier << " size=" << nb_elem_tot << " composantes=" << elements.dimension(1);

        ecrire_offset(sfichier, offset_elem_);


        switch(sizeof(_LATA_INT_TYPE_))
          {
          case 4:
            sfichier << " format=INT32" << finl;
            break;
          case 8:
            sfichier << " format=INT64" << finl;
            break;
          default:
            Cerr << "Error in Format_Post_Lata::ecrire_entete\n" << " sizeof(int) not supported" << finl;
            Process::exit();
          }

      }
    fichier_lata.syncfile();
  }

  // In parallel mode, additional files containing parallel data on vertices, elements and faces are written...
  if (Process::is_parallel())
    if (options_para_ == SINGLE_FILE || options_para_ == SINGLE_FILE_MPIIO)
      {
        TIDTab data(1,2);
        data(0, 0) = decalage_sommets;
        data(0, 1) = sommets.dimension(0);
        ecrire_item_tid("JOINTS_SOMMETS",
                        id_domaine,
                        "", /* id_domaine */
                        "", /* localisation */
                        "", /* reference */
                        data,
                        0); /* reference_size */
        data(0, 0) = decalage_elements;
        data(0, 1) = elements.dimension(0);
        ecrire_item_tid("JOINTS_ELEMENTS",
                        id_domaine,
                        "", /* id_domaine */
                        "", /* localisation */
                        "", /* reference */
                        data,
                        0); /* reference_size */
      }
}

/*! @brief See Format_Post_base::ecrire_domaine. Writing a domain within a time step is accepted, but
 *
 *   all id_domaines must be distinct.
 *   Writes the file "basename(_XXXXX).lata.nom_domaine", which contains the vertex list and element list.
 *   If the PE is master, opens the master file in APPEND mode and adds a reference to this file.
 */
int Format_Post_Lata::ecrire_domaine(const Domaine& domaine,const int est_le_premier_post)
{
  if (status == RESET)
    {
      Cerr << "Error in Format_Post_Lata::ecrire_domaine\n" << " status = RESET. Uninitialized object" << finl;
      Process::exit();
    }
  Motcle type_elem = domaine.type_elem()->que_suis_je();

  ecrire_domaine_low_level(domaine.le_nom(), domaine.les_sommets(), domaine.les_elems(), type_elem);

  // If there are domain boundary sub-domains, write them too
  const LIST(OBS_PTR(Domaine)) bords= domaine.domaines_frontieres();
  for (int i=0; i<bords.size(); i++)
    ecrire_domaine(bords[i].valeur(),est_le_premier_post);

  return 1; // ok all is well
}

/*! @brief Starts writing a new time step. For the LATA format specifically:
 *
 *   Opens the master file in APPEND mode and adds a line
 *    "TEMPS xxxxx" if this time has not yet been written.
 *
 */
int Format_Post_Lata::ecrire_temps(const double temps)
{
  ecrire_temps_lata(temps,temps_courant_,lata_basename_,status,options_para_);
  return 1;
}

/*! @brief voir Format_Post_base::ecrire_champ
 *
 */
int Format_Post_Lata::ecrire_champ(const Domaine& domaine, const Noms& unite_, const Noms& noms_compo, int ncomp, double temps, const Nom& id_du_champ, const Nom& id_du_domaine,
                                   const Nom& localisation, const Nom& nature, const DoubleTab& valeurs)
{
  Motcle id_du_champ_modifie(id_du_champ), iddomaine(id_du_domaine);

  //Using prefix with an uppercase argument
  if ((Motcle) localisation == "SOM")
    {
      id_du_champ_modifie.prefix(id_du_domaine);
      id_du_champ_modifie.prefix(iddomaine);
      id_du_champ_modifie.prefix("_SOM_");
    }
  else if ((Motcle) localisation == "ELEM")
    {
      id_du_champ_modifie.prefix(id_du_domaine);
      id_du_champ_modifie.prefix(iddomaine);
      id_du_champ_modifie.prefix("_ELEM_");
    }
  else if ((Motcle) localisation == "FACES")
    {
      id_du_champ_modifie.prefix(id_du_domaine);
      id_du_champ_modifie.prefix(iddomaine);
      id_du_champ_modifie.prefix("_FACES_");
    }
  Nom& id_champ = id_du_champ_modifie;

  // Build the file name
  Nom basename_champ(lata_basename_), extension_champ(extension_lata());

  if (un_seul_fichier_lata_) extension_champ += "_single";
  else
    {
      extension_champ += ".";
      extension_champ += id_champ;
      extension_champ += ".";
      extension_champ += localisation;
      extension_champ += ".";
      extension_champ += id_du_domaine;
      extension_champ += ".";
      char str_temps[100] = "0.0";
      if (temps >= 0.)
        snprintf(str_temps, 100, "%.10f", temps);
      extension_champ += str_temps;
    }

  Nom filename_champ;
  trustIdType size_tot;
  int nb_compo;
  {
    const bool not_in_list =  !liste_single_lata_ecrit.contient_(lata_basename_),
               should_erase = (!un_seul_fichier_lata_) ? true /* Always erase */ : (offset_elem_ < 0 && not_in_list);

    Fichier_Lata fichier_champ(basename_champ, extension_champ, should_erase ? Fichier_Lata::ERASE : Fichier_Lata::APPEND, format_, options_para_);

    // add to the list if not already in it and if un_seul_fichier_lata_ !!!
    if (not_in_list && un_seul_fichier_lata_) liste_single_lata_ecrit.add(lata_basename_); // BOOM !

    // XXX Elie Saikali : attention offset ici avant write_doubletab ! sinon decalage d'un champ !
    if (un_seul_fichier_lata_)
      if (fichier_champ.is_master())
        offset_elem_ = fichier_champ.get_SFichier().get_ofstream().tellp();

    filename_champ = fichier_champ.get_filename();
    size_tot = write_doubletab(fichier_champ, valeurs, nb_compo, options_para_);
  }

  // Opening the .lata file in append mode.
  // Adding the field reference
  Fichier_Lata_maitre fichier(lata_basename_, extension_lata(), Fichier_Lata::APPEND, options_para_);
  SFichier& sfichier = fichier.get_SFichier();
  if (fichier.is_master())
    {
      sfichier << "Champ " << id_champ << " ";
      sfichier << remove_path(filename_champ);
      sfichier << " geometrie=" << id_du_domaine;
      sfichier << " localisation=" << localisation;
      sfichier << " size=" << size_tot;
      sfichier << " nature=" << nature;
      sfichier << " noms_compo=" << noms_compo[0];
      for (int k = 1; k < noms_compo.size(); k++)
        sfichier << "," << noms_compo[k];

      sfichier << " composantes=" << nb_compo;

      ecrire_offset(sfichier, offset_elem_);

      sfichier << finl;
    }
  fichier.syncfile();

  return 1;
}

/*! @brief See Format_Post_base::ecrire_champ. WARNING: if "reference" is non-empty, 1 is added to all values to switch to Fortran numbering, and if furthermore a single lata file is written for all processors, an
 *
 *    offset is added to all values (renumbering of indices to switch to global numbering; see ecrire_domaine for an example)
 *
 */
template<typename TYP>
int Format_Post_Lata::ecrire_item_integral_T(const Nom& id_item, const Nom& id_du_domaine, const Nom& id_domaine, const Nom& localisation,
                                             const Nom& reference, const TRUSTVect<TYP, int>& val, const int reference_size)
{
  // Build the file name
  Nom basename_champ(lata_basename_), extension_champ(extension_lata());

  if (un_seul_fichier_lata_) extension_champ += "_single";
  else
    {
      extension_champ += ".";
      //extension_champ += id_champ;
      extension_champ += id_item;
      extension_champ += ".";
      extension_champ += localisation;
      extension_champ += ".";
      extension_champ += id_du_domaine;
      extension_champ += ".";
      char str_temps[100] = "0.0";
      if (temps_courant_ >= 0.)
        snprintf(str_temps, 100, "%.10f", temps_courant_);
      extension_champ += Nom(str_temps);
    }

  Nom filename_champ;
  trustIdType size_tot;
  int nb_compo = 0; // filled by ::write_T_tab() below, but my compiler doesn't see it since I templatized it.
  const TRUSTTab<TYP, int> valeurs = static_cast<const TRUSTTab<TYP, int>&>(val);
  {
    const bool not_in_list =  !liste_single_lata_ecrit.contient_(lata_basename_),
               should_erase = (!un_seul_fichier_lata_) ? true /* Always erase */ : (offset_elem_ < 0 && not_in_list);

    Fichier_Lata fichier_champ(basename_champ, extension_champ, should_erase ? Fichier_Lata::ERASE : Fichier_Lata::APPEND, format_, options_para_);

    // add to the list if not already in it and if un_seul_fichier_lata_ !!!
    if (not_in_list && un_seul_fichier_lata_) liste_single_lata_ecrit.add(lata_basename_); // BOOM !

    if (un_seul_fichier_lata_)
      if (fichier_champ.is_master())
        offset_elem_ = fichier_champ.get_SFichier().get_ofstream().tellp();

    filename_champ = fichier_champ.get_filename();
    // Assuming that if reference is non-empty, it is an index into another array, so Fortran numbering applies:
    bool decal = false;
    trustIdType decal_partiel = 0;
    if (reference != "")
      {
        decal = true;
        decal_partiel = 1;
        if (options_para_ == SINGLE_FILE || options_para_ == SINGLE_FILE_MPIIO)
          {
            // All processors write to a single file; indices must be renumbered to global numbering.
            // Offset to add to indices for global numbering.
            decal_partiel += mppartial_sum(reference_size);
          }
      }
    size_tot = ::write_T_tab<TYP, _LATA_INT_TYPE_>(fichier_champ, decal, decal_partiel, valeurs, nb_compo, options_para_);
  }

  {
    // Opening the .lata file in append mode. Adding the field reference.
    Fichier_Lata_maitre fichier(lata_basename_, extension_lata(), Fichier_Lata::APPEND, options_para_);
    SFichier& sfichier = fichier.get_SFichier();
    if (fichier.is_master())
      {
        sfichier << "Champ " << id_item << " ";
        sfichier << remove_path(filename_champ);
        sfichier << " geometrie=" << id_du_domaine;

        if (localisation != "") sfichier << " localisation=" << localisation;

        sfichier << " size=" << size_tot;
        sfichier << " composantes=" << nb_compo;
        if (reference != "") sfichier << " reference=" << reference;

        ecrire_offset(sfichier, offset_elem_);

        const int sz = (int) sizeof(_LATA_INT_TYPE_);
        switch(sz)
          {
          case 4:
            sfichier << " format=int32";
            break;
          case 8:
            sfichier << " format=int64";
            break;
          default:
            Cerr << "Error in Format_Post_Lata::ecrire_champ_lata\n" << " Integer type not supported: size=" << sz << finl;
            exit();
          }
        sfichier << finl;
      }
    fichier.syncfile();
  }

  // Trick for parallel face data:
  if ((id_item == "FACES" && Process::is_parallel()) && (options_para_ == SINGLE_FILE || options_para_ == SINGLE_FILE_MPIIO))
    {
      const int n = valeurs.dimension(0);
      TIDTab data(1,2);
      data(0, 0) = 1 + mppartial_sum(n);
      data(0, 1) = n;
      ecrire_item_tid("JOINTS_FACES",
                      id_du_domaine,
                      "", /* id_domaine */
                      "", /* localisation */
                      "", /* reference */
                      data,
                      0); /* reference_size */
    }
  return 1;
}

int Format_Post_Lata::ecrire_item_int(const Nom& id_item, const Nom& id_du_domaine, const Nom& id_domaine, const Nom& localisation,
                                      const Nom& reference, const IntVect& val, const int reference_size)
{
  return ecrire_item_integral_T(id_item, id_du_domaine, id_domaine, localisation, reference, val, reference_size);
}

int Format_Post_Lata::ecrire_item_tid(const Nom& id_item, const Nom& id_du_domaine, const Nom& id_domaine, const Nom& localisation,
                                      const Nom& reference, const TIDVect& val, const int reference_size)
{
  return ecrire_item_integral_T(id_item, id_du_domaine, id_domaine, localisation, reference, val, reference_size);
}

int Format_Post_Lata::ecrire_entete_lata(const Nom& base_name, const Options_Para& option, const Format& format, const int est_le_premier_post)
{
  if (est_le_premier_post)
    {
      // Determine the binary format:
      //  big endian => the 32-bit integer "1" is written as 0x00 0x00 0x00 0x01
      //  little endian =>                                   0x01 0x00 0x00 0x00
      const unsigned int one = 1;
      const int big_endian = (*((unsigned char*) &one) == 0) ? 1 : 0;

      // Erase the .lata file and write the header
      Fichier_Lata_maitre fichier(base_name, extension_lata(), Fichier_Lata::ERASE, option);

      SFichier& sfichier = fichier.get_SFichier();
      if (fichier.is_master())
        {
          sfichier << "LATA_V2.1 TRUST version " << TRUST_VERSION << finl;
          sfichier << Objet_U::nom_du_cas() << finl;
          sfichier << "Trio_U verbosity=0" << finl;

          sfichier << "Format ";
          switch(format)
            {
            case ASCII:
              sfichier << "ASCII,";
              break;
            case BINAIRE:
              sfichier << "BINAIRE,";
              if (big_endian)
                sfichier << "BIG_ENDIAN,";
              else
                sfichier << "LITTLE_ENDIAN,";
              break;
            default:
              Cerr << "Error in Format_Post_Lata::ecrire_entete\n" << " format not supported" << finl;
              exit();
            }
          switch(sizeof(_LATA_INT_TYPE_))
            {
            case 4:
              sfichier << "INT32,";
              break;
            case 8:
              sfichier << "INT64,";
              break;
            default:
              Cerr << "Error in Format_Post_Lata::ecrire_entete\n" << " sizeof(int) not supported" << finl;
              exit();
            }
          sfichier << "F_INDEXING,C_ORDERING,F_MARKERS_SINGLE,REAL32" << finl;
        }
      fichier.syncfile();

    }

  return 1;
}

int Format_Post_Lata::ecrire_temps_lata(const double temps, double& temps_format, const Nom& base_name, Status& stat, const Options_Para& option)
{
  assert(stat != RESET);
  // Write the time only if it has changed...
  if (stat != WRITING_TIME || temps_format != temps)
    {
      temps_format = temps;
      // Opening the .lata file in append mode
      Fichier_Lata_maitre fichier(base_name, extension_lata(), Fichier_Lata::APPEND, option);
      if (fichier.is_master())
        fichier.get_SFichier() << "TEMPS " << temps << finl;
      fichier.syncfile();
      stat = WRITING_TIME;
    }
  return 1;
}

int Format_Post_Lata::finir(const int est_le_dernier_post)
{
  if (est_le_dernier_post)
    {
      Fichier_Lata_maitre fichier(lata_basename_, extension_lata(), Fichier_Lata::APPEND, options_para_);
      SFichier& sfichier = fichier.get_SFichier();
      if (fichier.is_master()) sfichier << "FIN" << finl;
      fichier.syncfile();
    }
  return 1;
}
