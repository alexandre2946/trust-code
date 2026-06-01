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
#include <Format_Post_base.h>
#include <Domaine_VF.h>
#include <Param.h>

Implemente_base(Format_Post_base,"Format_Post_base",Objet_U);

/*! @brief erreur => exit
 *
 */
Sortie& Format_Post_base::printOn(Sortie& os) const
{
  Cerr << "Format_Post_base::printOn : error" << finl;
  exit();
  return os;
}

Entree& Format_Post_base::readOn(Entree& is)
{
  Cerr<<"Reading of data for a "<<que_suis_je()<<" post-processing format object"<<finl;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Format_Post_base::resetTime(double t, const std::string dirname)
{
  if (dirname.empty())
    {
      Cerr << "Format '" << que_suis_je() << " does not support resetTime()!!" << finl;
      Process::exit(-1);
      // but LATA does :-)
    }
}

int Format_Post_base::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  return -1;
}

/*! @brief Initializes the file with parameters appropriate to its format (e.g. ascii format, deletes the existing file, a
 *
 *   single file for all processors, etc.)
 *   Method to override in derived classes.
 *  Return value: 1 if the operation succeeded, 0 otherwise.
 *
 */
int Format_Post_base::initialize_by_default(const Nom& file_basename)
{
  Cerr << "Format_Post_base::initialize_by_default(" << file_basename
       << ")\n method not coded for " << que_suis_je() << finl;
  return 0;
}

int Format_Post_base::initialize(const Nom& file_basename, const int format, const Nom& option_para)
{
  Cerr << "Format_Post_base::initialize(" << file_basename
       << ")\n method not coded for " << que_suis_je() << finl;
  return 0;
}

/*! @brief Modification of the post processing file name.
 * For save/restart this might also move and rename files around to avoid overriding existing files.
 */
int Format_Post_base::modify_file_basename(const Nom file_basename, bool for_restart, const double tinit)
{
  return 0;
}

int Format_Post_base::ecrire_entete(const double temps_courant,const int reprise,const int est_le_premier_post)
{
  Cerr << "Format_Post_base::ecrire_entete method not coded for " << que_suis_je() << finl;
  return 0;
}

int Format_Post_base::finir(const int est_le_dernier_post)
{
  Cerr << "Format_Post_base::finir method not coded for " << que_suis_je() << finl;
  return 0;
}

int Format_Post_base::init_ecriture(double temps_courant,double temps_post,
                                    int est_le_premier_postraitement_pour_nom_fich_,const Domaine& domaine)
{
  return 1;
}

int Format_Post_base::finir_ecriture(double temps_courant)
{
  return 1;
}



int Format_Post_base::completer_post(const Domaine& dom,const int is_axi,
                                     const Nature_du_champ& nature,const int nb_compo,const Noms& noms_compo,
                                     const Motcle& loc_post,const Nom& le_nom_champ_post)
{

  Cerr << "Format_Post_base::preparer_post_champ(...)\n"
       << " method not coded for " << que_suis_je() << finl;
  return 0;

}

int Format_Post_base::preparer_post(const Nom& id_du_domaine,const int est_le_premier_post,
                                    const int reprise,
                                    const double t_init)
{

  Cerr << "Format_Post_base::preparer_post(...)\n"
       << " method not coded for " << que_suis_je() << finl;
  return 0;

}


/*! @brief Writing a mesh.
 *
 * The domain is written either at the beginning, before the first call to ecrire_temps, or later (dynamic mesh),
 *   but this is not necessarily supported by all post-processing formats.
 *
 *
 * @param (id_domaine) the name assigned to the domain in the lata file.
 * @param (type_elem) the type of geometric element (a type supported by the derived class; in general at least "TETRAEDRE", "HEXAEDRE", "TRIANGLE" and "RECTANGLE" are understood)
 * @param (dimension) the dimension of the domain (number of vertex coordinates). A 3D domain may contain triangle elements (post-processing of an interface or of the boundary of a volumetric domain).
 * @param (sommets) Vertex coordinates. If non-empty, dimension(1) must equal dimension.
 * @param (elements) Indices of vertices for each element. dimension(1) must match the element type (3 for a triangle, 4 for a rectangle or tetrahedron, etc.)
 */

int Format_Post_base::ecrire_domaine(const Domaine& domaine,const int est_le_premier_post)
{
  Cerr << "Format_Post_base::ecrire_domaine(...)\n"
       << " method not coded for " << que_suis_je() << finl;
  return 0;
}

void Format_Post_base::ecrire_domaine_dual(const Domaine& domaine,const int est_le_premier_post)
{
  Cerr << "ERROR: Format_Post_base::ecrire_domaine_dual(...) method not coded for " << que_suis_je() << finl;
  Process::exit();
}

int Format_Post_base::ecrire_domaine_dis(const Domaine& domaine,const OBS_PTR(Domaine_dis_base)& domaine_dis_base,const int est_le_premier_post)
{
  domaine_dis_ = domaine_dis_base;
  return ecrire_domaine(domaine, est_le_premier_post);
}

/*! @brief Starts writing a time step.
 *
 * The derived class must accept receiving multiple consecutive calls to this method with the same time.
 *
 */

int Format_Post_base::ecrire_temps(const double temps)
{
  Cerr << "Format_Post_base::ecrire_temps(const double temps)\n"
       << " method not coded for " << que_suis_je() << finl;
  return 0;
}

/*! @brief Writing a field to the post-processing file.
 *
 * @param (id_du_champ) field identifier (allows identifying a unique field when combined with a time step number)
 * @param (id_du_domaine) identifier of the domain on which the field is defined. This domain must have been written before by "ecrire_domaine".
 * @param (localisation) location of field values (SOMMETS, ELEMENTS, or any other id of an already written array) (not everything is necessarily supported by all post-processing formats)
 * @param (data) array of values to post-process. The number of rows must equal the number of rows of the "localisation" array (number of vertices, elements, faces, etc.). Return value: 1 if the operation succeeded, 0 otherwise (e.g., preconditions not met or feature not supported by the format).
 */

int Format_Post_base::ecrire_champ(const Domaine& domaine,const Noms& unite_,const Noms& noms_compo,
                                   int ncomp, double temps_,
                                   const Nom&   id_du_champ,
                                   const Nom&         id_du_domaine,
                                   const Nom&         localisation,
                                   const Nom&   nature,
                                   const DoubleTab&   data)
{
  Cerr << "Format_Post_base::ecrire_champ(...)\n"
       << " method not coded for " << que_suis_je() << finl;
  return 0;
}

/*! @brief Writing an integer array to the post-processing file.
 *
 * @sa ecrire_champ
 *
 * @param (reference) Name of another array already written. data[i] is an index into that array (e.g., the FACES_VOISINS array can be written with localisation=FACES, reference=ELEMENTS because the array is indexed by face number and contains element indices). Writing a domain automatically creates a SOMMETS array and an ELEMENTS array.
 * @param (reference_size) the (local) size of the referenced array (dimension(0) of the array). This dimension is used to renumber the contents of the data array to create a global numbering when all processors write to a single file.
 */
int Format_Post_base::ecrire_item_int(const Nom&   id_item,
                                      const Nom&   id_du_domaine,
                                      const Nom& id_domaine,
                                      const Nom&   localisation,
                                      const Nom&   reference,
                                      const IntVect& data,
                                      const int reference_size)
{
  Cerr << "Format_Post_base::ecrire_champ_int(...)\n"
       << " method not coded for " << que_suis_je() << finl;
  return 0;
}

