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

#include <Moyenne_volumique.h>
#include <communications.h>
#include <Equation_base.h>
#include <Postraitement.h>
#include <Octree_Double.h>
#include <Domaine_VF.h>
#include <TRUST_Ref.h>
#include <algorithm>
#include <Param.h>

Implemente_instanciable(Moyenne_volumique,"Moyenne_volumique",Interprete);
// XD moyenne_volumique interprete moyenne_volumique BRACE This keyword should be used after Resoudre keyword. It
// XD_CONT computes the convolution product of one or more fields with a given filtering function.

Sortie& Moyenne_volumique::printOn(Sortie& s ) const
{
  return s << que_suis_je() << finl;
}

/*! @brief Reading of the filtering function.
 *
 * @brief Expected format:
 * {
 *     type BOITE|CHAPEAU|QUADRA|GAUSSIENNE|PARSER
 *     demie-largeur L
 *     [ omega W ]
 *     [ expression FORMULE ]
 *   }
 *
 * @param is the input stream
 * @return the modified input stream
 */
Entree& Moyenne_volumique::readOn(Entree& is )
{
  expression_parser_ = "??";
  int type = -1;
  Param param(que_suis_je());
  param.ajouter("type", & type, Param::REQUIRED /* obligatoire */);
  param.dictionnaire("BOITE", BOITE);
  param.dictionnaire("CHAPEAU", CHAPEAU);
  param.dictionnaire("GAUSSIENNE", GAUSSIENNE);
  param.dictionnaire("QUADRA", QUADRA);
  param.dictionnaire("PARSER", PARSER);
  param.ajouter("demie-largeur", & box_size_, Param::REQUIRED /* obligatoire */);
  param.ajouter("omega", & l_);
  param.ajouter("expression", & expression_parser_);
  param.lire_avec_accolades_depuis(is);
  switch(type)
    {
    case BOITE:
      type_ = BOITE;
      l_ = box_size_;
      break;
    case CHAPEAU:
      type_ = CHAPEAU;
      l_ = box_size_;
      break;
    case GAUSSIENNE:
      type_ = GAUSSIENNE;
      if (l_ < 0.)
        {
          Cerr << "Error : OMEGA must be set to >= 0" << finl;
          barrier();
          exit();
        }
      break;
    case PARSER:
      type_ = PARSER;
      if (expression_parser_ == "??")
        {
          Cerr << "Error : EXPRESSION must be specified for the parser." << finl;
          barrier();
          exit();
        }
      {
        std::string s(expression_parser_);
        std::transform(s.begin(), s.end(), s.begin(), ::toupper);
        parser_.setString(s);
        parser_.setNbVar(Objet_U::dimension);
        parser_.addVar("x");
        parser_.addVar("y");
        if (Objet_U::dimension == 3)
          parser_.addVar("z");
        parser_.parseString();
      }
      break;
    case QUADRA:
      type_ = QUADRA;
      l_ = box_size_;
      Cerr << "l_ = " << l_ << "box_size_ = " << box_size_ << finl;
      break;
    default:
      exit(); // Internal error!
    }
  // Slight enlargement so that the octree correctly finds elements exactly on the boundary.
  box_size_ += precision_geom;
  return is;
}
inline double fonction_quadra(double x, double l_)
{
  assert(std::fabs(x) <= l_);
  double ax = 1. - std::fabs(x) / l_;
  ax = ax * ax;
  if (std::fabs(x) < (l_/3.))
    {
      double bx = -3. / l_ * std::fabs(x) + 1.;
      ax -= bx * bx / 3.; // signifie ax = ax - bx*bx/3
    }
  return ax * (27. / (16. * l_));
}
/*! @brief Evaluates the filter function at each coordinate in coords.
 *
 * @brief Method called from the Calcul_integrale_locale class.
 * @param coords the array of coordinates at which to evaluate the filter
 * @param result the output array of filter values
 */
void Moyenne_volumique::eval_filtre(const DoubleTab& coords, ArrOfDouble& result) const
{
  const int dim = Objet_U::dimension;
  assert(dim == coords.dimension(1));
  const int n = coords.dimension(0);
  assert(result.size_array() == n);
  switch(type_)
    {
    case PARSER:
      {
        for (int i = 0; i < n; i++)
          {
            for (int j = 0; j < dim; j++)
              parser_.setVar(j, coords(i,j));
            result[i] = parser_.eval();
          }
        break;
      }
    case BOITE:
      {
        double facteur = 0.;
        if (dim == 2)
          facteur = 1. / (l_ * l_ * 4.);
        else
          facteur = 1. / (l_ * l_ * l_ * 8.);

        for (int i = 0; i < n; i++)
          {
            const double x = coords(i,0);
            const double y = coords(i,1);
            const double z = (dim==3) ? coords(i,2) : 0.;
            double r = facteur;
            if (x > l_ || x < -l_
                || y > l_ || y < -l_
                || z > l_ || z < -l_)
              r = 0.;
            result[i] = r;
          }
        break;
      }
    case CHAPEAU:
      {
        const double L2D = l_*l_*l_*l_;
        const double L3D = L2D*l_*l_;
        double facteur;
        if (dim==3)
          facteur = 1. / L3D;
        else
          facteur = 1. / (L2D * l_);
        const int nbis = coords.dimension(0);
        for (int i = 0; i < nbis; i++)
          {
            const double x = coords(i, 0);
            const double y = coords(i, 1);
            const double z = (dim == 3) ? coords(i, 2) : 0.;
            double resu = 0.;
            const double ax = std::fabs(x);
            const double ay = std::fabs(y);
            const double az = std::fabs(z);
            if (ax <= l_ && ay <= l_ && az <= l_)
              resu = (l_-ax) * (l_-ay) * (l_-az) * facteur;
            result[i] = resu;
          }
        break;
      }
    case GAUSSIENNE:
      {
        double facteur1 = - 0.5 / (l_ * l_);
        double facteur2 = 1. / (l_ * sqrt(2 * M_PI));
        if (dim == 2)
          facteur2 = facteur2 * facteur2;
        else
          facteur2 = facteur2 * facteur2 * facteur2;
        for (int i = 0; i < n; i++)
          {
            const double x = coords(i, 0);
            const double y = coords(i, 1);
            const double z = (dim == 3) ? coords(i, 2) : 0.;
            const double k = (x*x + y*y + z*z) * facteur1;
            result[i] = exp(k) * facteur2;
          }
        break;
      }
    case QUADRA:
      {
        for (int i = 0; i < n; i++)
          {
            const double x = coords(i, 0);
            const double y = coords(i, 1);
            const double z = (dim == 3) ? coords(i, 2) : 0.;
            double resu = 0;
            if (std::fabs(x) < l_ && std::fabs(y) < l_ && std::fabs(z) < l_)
              {
                resu = fonction_quadra(x,l_) * fonction_quadra(y,l_);
                if (dim == 3)
                  resu *= fonction_quadra(z,l_);
              }
            result[i] = resu;
          }
        break;
      }
    default:
      {
        Cerr << "Error in Moyenne_volumique::eval() : filter function is not initialized." << finl;
        exit();
      }
    }
}

/*! @brief Searches for the field named "nom_champ" in the problem named "nom_pb" among the interpreter objects.
 *
 * @brief Method called by traiter_champs().
 * @param nom_pb the name of the problem
 * @param nom_champ the name of the field to retrieve
 * @param ref_champ reference to the field, set on output
 * @return 1 on success
 */
int Moyenne_volumique::get_champ(const Nom& nom_pb,
                                 const Nom& nom_champ,
                                 OBS_PTR(Champ_base) & ref_champ)
{
  Probleme_base& pb = ref_cast(Probleme_base, objet(nom_pb));
  // Le champ est-il defini dans les postraitements (statistiques) ?
  const int nb_post = pb.postraitements().size();
  Motcle mc_nom_champ(nom_champ);
  for (int i_post = 0; i_post < nb_post; i_post++)
    {
      if (sub_type(Postraitement_base, pb.postraitements()[i_post].valeur()))
        {
          Postraitement& post = ref_cast(Postraitement, pb.postraitements()[i_post].valeur());
          Operateurs_Statistique_tps& stats = post.les_statistiques();
          const int nstat = stats.size();
          for (int i_stat = 0; i_stat < nstat; i_stat++)
            {
              Motcle tmp(stats[i_stat]->le_nom() );

              if (tmp == mc_nom_champ)
                {
                  Operateur_Statistique_tps_base& stat = stats[i_stat].valeur();
                  ref_cast_non_const(DoubleTab, stat.integrale().le_champ_calcule().valeurs()) = stat.calculer_valeurs();
                  ref_champ = stat.integrale().le_champ_calcule();
                  return 1;
                }
            }
        }
    }

  ref_champ = pb.get_champ(nom_champ);
  return 1;
}

/*! @brief Helper function that performs the convolution calculations and writes the result to a lata file for all fields of a given type listed in noms_champs.
 *
 * @brief Method called by interpreter().
 *  type_champ=0 => process element fields
 *  type_champ=1 => process face fields
 * @param noms_champs list of field names to process
 * @param nom_pb name of the problem
 * @param nom_dom name of the destination domain
 * @param coords coordinates at which to evaluate the convolution
 * @param post the post-processing format object used for writing
 * @param temps current time
 * @param localisation field localisation (ELEM or SOM)
 */
void Moyenne_volumique::traiter_champs(const Motcles& noms_champs,
                                       const Nom& nom_pb, const Nom& nom_dom,
                                       const DoubleTab& coords,
                                       Format_Post_base& post,
                                       double temps,
                                       const Motcle& localisation)
{
  const Domaine& dom_post = ref_cast(Domaine, objet(nom_dom));
  const int nb_champs = noms_champs.size();
  if (nb_champs == 0)
    return;

  OBS_PTR(Champ_base) ref_champ;
  OBS_PTR(Domaine_VF) ref_domaine_vf;
  int i_champ;
  // ************************************
  // Compute the total number of components and ref_domaine_vf
  int nb_compo_tot = 0;
  for (i_champ = 0; i_champ < nb_champs; i_champ++)
    {
      get_champ(nom_pb, noms_champs[i_champ], ref_champ);
      const Champ_base& champ = ref_champ.valeur();
      const Domaine_VF& zvf = ref_cast(Domaine_VF, champ.domaine_dis_base());
      if (i_champ == 0)
        {
          ref_domaine_vf = zvf;
        }
      else
        {
          if (& (ref_domaine_vf.valeur()) != & zvf)
            {
              Cerr << "Error in Moyenne_volumique::traiter_champs all the fields must be discretized on the same Domaine." << finl;
              barrier();
              exit();
            }
        }
      const int nb_compo = champ.nb_comp();
      nb_compo_tot += nb_compo;
    }

  const Domaine_VF& domaine_source = ref_domaine_vf.valeur();

  // ************************************
  // Build a large array containing all the values to process plus the porosity
  DoubleTab valeurs_src;
  const int nb_lignes = domaine_source.nb_elem();
  valeurs_src.resize(nb_lignes, nb_compo_tot + 1);
  int count = 0;
  {
    DoubleTab tmp_val;
    IntVect liste_elems(nb_lignes);
    for (int i = 0; i < nb_lignes; i++)
      liste_elems[i] = i;
    const DoubleTab& xp = domaine_source.xp();
    for (i_champ = 0; i_champ < nb_champs; i_champ++)
      {
        get_champ(nom_pb, noms_champs[i_champ], ref_champ);
        const Champ_base& champ = ref_champ.valeur();
        const int nb_compo = champ.nb_comp();
        tmp_val.reset();
        tmp_val.resize(nb_lignes, nb_compo);
        champ.valeur_aux_elems(xp, liste_elems, tmp_val);

        for (int i = 0; i < nb_lignes; i++)
          for (int j = 0; j < nb_compo; j++)
            valeurs_src(i, count+j) = tmp_val(i, j);
        count += nb_compo;

        Cout << "Field name = " << ref_champ->le_nom()
             << " Field type = " << ref_champ->que_suis_je() << finl;
      }
  }
  // et une colonne de 1:
  for (int i = 0; i < nb_lignes; i++)
    valeurs_src(i, count) = 1.;

  // Tableau de resultats:
  const int nb_coords = coords.dimension(0);
  DoubleTab resu(nb_coords, nb_compo_tot + 1);

  // ************************************
  // Compute all convolution products

  calculer_convolution_champ_elem(domaine_source,
                                  valeurs_src,
                                  coords,
                                  resu);
  Noms nom_dir;
  nom_dir.add("_X");
  nom_dir.add("_Y");
  nom_dir.add("_Z");
  count = 0;
  for (i_champ = 0; i_champ < nb_champs; i_champ++)
    {
      get_champ(nom_pb, noms_champs[i_champ], ref_champ);
      const Champ_base& champ = ref_champ.valeur();
      const int nb_compo = champ.nb_comp();
      DoubleTab extrait(nb_coords, nb_compo);
      for (int i = 0; i < nb_coords; i++)
        for (int j = 0; j < nb_compo; j++)
          extrait(i, j) = resu(i, count + j);

      Cout << "Post writting " << champ.le_nom() << finl;
      Nom nature("scalar");
      if (champ.nature_du_champ()==vectoriel) nature="vector";
      post.ecrire_champ(dom_post,
                        champ.unites(),
                        champ.noms_compo(),
                        -1 /* ecrire toutes les composantes */, temps,
                        noms_champs[i_champ], nom_dom, localisation,nature, extrait);

      count += nb_compo;
    }
  // Derniere colonne (porosite)
  const int nb_compo = 1;
  DoubleTab extrait(nb_coords, nb_compo);
  for (int i = 0; i < nb_coords; i++)
    for (int j = 0; j < nb_compo; j++)
      extrait(i, j) = resu(i, count + j);

  Noms noms_compo;
  Noms unites;
  Nom nom_moyenne;
  unites.add("m3");
  noms_compo.add("porosite");
  nom_moyenne = "porosite";
  Cout << "Porosity post writing" << finl;

  post.ecrire_champ(dom_post, unites, noms_compo, -1 /* ecrire toutes les composantes */,
                    temps,
                    nom_moyenne, nom_dom, localisation, "scalar",extrait);
}

/*! @brief Reads the parameters from the data set.
 *
 * @brief Expected format: Moyenne_volumique {
 *     nom_pb NOM_DU_PROBLEME    (where to look for the source fields)
 *     nom_domaine DOMAINE_CIBLE (the convolution is evaluated at the elements of this domain)
 *     noms_champs N CHAMP1 CHAMP2 ... (names of the fields to filter in the problem)
 *     [ nom_fichier_post NOM_SANS_EXTENSION ] (either nom_fichier and format_post are given,
 *                                              or fichier_post is given)
 *     [ format_post lata|lml|med|... ] (default: lata)
 *     [ fichier_post Format_Post_XXX { ... } ] (read via readOn of Format_Post_XXX)
 *     fonction_filtre ...  (format: see Moyenne_volumique::readOn() )
 *     [ localisation ELEM|SOM ]
 *   }
 *
 * @param is the input stream
 * @return the modified input stream
 */
Entree& Moyenne_volumique::interpreter(Entree& is)
{
  Cerr << "Starting of Moyenne_volumique::interpreter" << finl;
  Nom nom_pb, nom_dom;
  Motcles noms_champs;
  Param param(que_suis_je() + Nom("::interpreter()"));
  const int id_elem = 0;
  const int id_som = 1;
  int localisation = id_elem; // default
  Motcle format_post("lata_v1");
  Nom nom_fichier_post;
  OWN_PTR(Format_Post_base) fichier_post;
  param.ajouter("nom_pb", & nom_pb, Param::REQUIRED); // XD_ADD_P ref_Pb_base
  // XD_CONT name of the problem where the source fields will be searched.
  param.ajouter("nom_domaine", & nom_dom, Param::REQUIRED); // XD_ADD_P ref_domaine
  // XD_CONT name of the destination domain (for example, it can be a coarser mesh, but for optimal performance in
  // XD_CONT parallel, the domain should be split with the same algorithm as the computation mesh, eg, same tranche
  // XD_CONT parameters for example)
  param.ajouter("noms_champs", & noms_champs, Param::REQUIRED); // XD_ADD_P listchaine
  // XD_CONT name of the source fields (these fields must be accessible from the postraitement) N source_field1
  // XD_CONT source_field2 ... source_fieldN
  param.ajouter("fichier_post", & fichier_post);
  param.ajouter("format_post", & format_post); // XD_ADD_P chaine
  // XD_CONT gives the fileformat for the result (by default : lata)
  param.ajouter("nom_fichier_post", & nom_fichier_post); // XD_ADD_P chaine
  // XD_CONT indicates the filename where the result is written
  // The Moyenne_volumique object is an interpreter, but it is also an object
  // whose only property is the filter function to use. Trick:
  // call the class readOn to read the filter function.
  param.ajouter("fonction_filtre", this, Param::REQUIRED); // XD_ADD_P bloc_lecture
  // XD_CONT to specify the given filter NL2 Fonction_filtre {NL2 type filter_typeNL2 demie-largeur lNL2 [ omega w ] NL2
  // XD_CONT [ expression string ]NL2 } NL2 NL2 type filter_type : This parameter specifies the filtering function.
  // XD_CONT Valid filter_type are:NL2 Boite is a box filter, $f(x,y,z)=(abs(x)<l)*(abs(y) <l)*(abs(z) <l) / (8 l^3)$NL2
  // XD_CONT Chapeau is a hat filter (product of hat filters in each direction) centered on the origin, the half-width
  // XD_CONT of the filter being l and its integral being 1.NL2 Quadra is a 2nd order filter.NL2 Gaussienne is a
  // XD_CONT normalized gaussian filter of standard deviation sigma in each direction (all field elements outside a
  // XD_CONT cubic box defined by clipping_half_width are ignored, hence, taking clipping_half_width=2.5*sigma yields an
  // XD_CONT integral of 0.99 for a uniform unity field).NL2 Parser allows a user defined function of the x,y,z
  // XD_CONT variables. All elements outside a cubic box defined by clipping_half_width are ignored. The parser is much
  // XD_CONT slower than the equivalent c++ coded function...NL2 NL2 demie-largeur l : This parameter specifies the half
  // XD_CONT width of the filterNL2 [ omega w ] : This parameter must be given for the gaussienne filter. It defines the
  // XD_CONT standard deviation of the gaussian filter.NL2 [ expression string] : This parameter must be given for the
  // XD_CONT parser filter type. This expression will be interpreted by the math parser with the predefined variables x,
  // XD_CONT y and z.
  param.ajouter("localisation", & localisation); // XD_ADD_P chaine(into=["elem","som"])
  // XD_CONT indicates where the convolution product should be computed: either on the elements or on the nodes of the
  // XD_CONT destination domain.
  param.dictionnaire("ELEM", id_elem);
  param.dictionnaire("SOM", id_som);
  param.lire_avec_accolades_depuis(is);

  // retrieve the domain
  const Domaine& dom = ref_cast(Domaine, objet(nom_dom));
  if (noms_champs.size() == 0)
    {
      Cerr << "Moyenne_volumique : no field to treat" << finl;
      return is;
    }
  Cerr << "Writing of the post-processing domain : " << nom_dom << finl;
  OBS_PTR(Champ_base) ref_champ;
  get_champ(nom_pb, noms_champs[0], ref_champ);
  const double temps = ref_champ->temps();
  if (!fichier_post)
    {
      if (nom_fichier_post == "??")
        {
          Cerr << "Error in Moyenne_volumique::interpreter:\n"
               << " missing NOM_FICHIER_POST or FICHIER_POST keyword" << finl;
          barrier();
          exit();
        }
      if (format_post == "lata_v2")
        format_post = "lata";
      // Trick to allow the non-regression test case to run in both sequential and parallel:
      //  (the output file name must match the case name)
      if (nom_fichier_post == "NOM_DU_CAS")
        {
          Cerr << "Post filename = NOM_DU_CAS => using " << nom_du_cas() << " instead" << finl;
          nom_fichier_post = nom_du_cas();
        }
      fichier_post.typer(Motcle("FORMAT_POST_") + format_post);
      fichier_post->initialize(nom_fichier_post, 1 /* binaire */, "SIMPLE");
    }
  else
    {
      if (nom_fichier_post != "??")
        {
          Cerr << "Error in Moyenne_volumique::interpreter:\n"
               << " you cannot give NOM_FICHIER_POST and FICHIER_POST. Choose one" << finl;
          barrier();
          exit();
        }
    }
  Format_Post_base& post = fichier_post.valeur();
  post.ecrire_entete(temps, 0 /*reprise*/, 1 /* premier post */);
  post.ecrire_domaine(dom, 1 /* premier_post */);
  post.ecrire_temps(temps);

  // Coordinates of the element centres of the destination domain
  DoubleTab coords;
  if (localisation == id_elem)
    {
      dom.calculer_centres_gravite(coords);
      // The array also contains virtual elements but without a virtual space. Beware.
      coords.resize(dom.nb_elem(), coords.dimension(1));
    }
  else
    {
      coords = dom.les_sommets();
    }

  Motcle loc("ELEM");
  if (localisation == id_som)
    loc = "SOM";
  traiter_champs(noms_champs,
                 nom_pb,
                 nom_dom,
                 coords,
                 post,
                 temps,
                 loc);
  int fin=1;
  post.finir(fin);
  return is;
}

/*! @brief Helper class used internally by calculer_convolution().
 *
 */
class Calcul_integrale_locale
{
public:
  Calcul_integrale_locale(const Domaine_VF& domaine_source,
                          const Moyenne_volumique& filter,
                          const DoubleTab& champ_source);
  void calculer(double x, double y, double z, ArrOfDouble& resu);

protected:
  const Domaine_VF& domaine_source_;
  Octree_Double octree_;
  const Moyenne_volumique& filter_;
  const DoubleTab& champ_source_;
  // Temporary arrays used in calculer():
  ArrOfInt liste_elems_;
  DoubleTab filter_coords_;
  ArrOfDouble filter_results_;
  int nb_items_reels_;
};

/*! @brief Constructor of the helper class.
 *
 * @brief See Moyenne_volumique::calculer_convolution().
 * @param domaine_source the source discretized domain
 * @param filter the volumetric average filter object
 * @param champ_source the source field values array
 */
Calcul_integrale_locale::Calcul_integrale_locale(const Domaine_VF& domaine_source,
                                                 const Moyenne_volumique& filter,
                                                 const DoubleTab& champ_source) :
  domaine_source_(domaine_source),
  filter_(filter),
  champ_source_(champ_source)
{



  // Build an octree containing the element centres.
  // The array is copied because it will be resized:
  DoubleTab coords = domaine_source.xp();
  nb_items_reels_ = domaine_source.nb_elem();
  // The xp array is dimensioned with dimension(0)=nb_elem_tot; resize it to nb_elem.
  coords.resize(nb_items_reels_, coords.dimension(1));
  if (champ_source.dimension(0) != nb_items_reels_)
    {
      Cerr << "Error in Calcul_integrale_locale::Calcul_integrale_locale() :\n"
           << " The source field is not discretized at the elements" << finl;
      Process::barrier();
      Process::exit();
    }
  octree_.build_nodes(coords, 0 /* no virtual items */);
}

/*! @brief Evaluates the convolution product "filter_ * champ_source_" at point x,y,z and stores the result in resu.
 *
 * @brief The elements to use are determined based on the filter size using an octree.
 *   The source field is assumed to be element-centred.
 *   Method called by Moyenne_volumique::calculer_convolution().
 * @param x x-coordinate of the evaluation point
 * @param y y-coordinate of the evaluation point
 * @param z z-coordinate of the evaluation point
 * @param resu output array receiving the convolution result
 */
void Calcul_integrale_locale::calculer(const double x, const double y, const double z,
                                       ArrOfDouble& resu)
{
  const double box_size = filter_.box_size();
  octree_.search_elements_box(x - box_size, y - box_size, z - box_size,
                              x + box_size, y + box_size, z + box_size,
                              liste_elems_);

  const DoubleTab& coord_items = domaine_source_.xp();

  const int nb_items = liste_elems_.size_array();
  const int dim = Objet_U::dimension;
  filter_coords_.resize(nb_items, dim);
  for (int i = 0; i < nb_items; i++)
    {
      const int item = liste_elems_[i];
      assert(item < nb_items_reels_);
      filter_coords_(i, 0) = coord_items(item, 0) - x;
      filter_coords_(i, 1) = coord_items(item, 1) - y;
      if (dim == 3)
        filter_coords_(i, 2) = coord_items(item, 2) - z;
    }
  filter_results_.resize_array(nb_items);
  filter_.eval_filtre(filter_coords_, filter_results_);

  const DoubleVect& volumes = domaine_source_.volumes();
  const int nb_comp = champ_source_.dimension(1);
  resu = 0.;
  for (int i = 0; i < nb_items; i++)
    {
      const int item = liste_elems_[i];
      const double valeur_filtre = filter_results_[i];
      const double volume = volumes(item);
      const double facteur = valeur_filtre * volume;
      for (int j = 0; j < nb_comp; j++)
        {
          // The integral is coarsely discretized as the product of the
          // values at the element centre times the element volume:
          const double valeur_champ = champ_source_(item, j);
          resu[j] += valeur_champ * facteur;
        }
    }
}

/*! @brief General method to compute a convolution from a field defined at elements or faces.
 *
 * @brief Method called by calculer_convolution_champ_elem() and calculer_convolution_champ_face().
 * @param domaine_source the source discretized domain
 * @param champ_source the source field values array
 * @param coords_to_compute coordinates at which to evaluate the convolution
 * @param resu output array of convolution results
 */
void Moyenne_volumique::calculer_convolution(const Domaine_VF& domaine_source,
                                             const DoubleTab& champ_source,
                                             const DoubleTab& coords_to_compute,
                                             DoubleTab& resu) const
{
  assert(resu.dimension(0) == coords_to_compute.dimension(0));
  const int dim = Objet_U::dimension;
  const int nbproc = Process::nproc();
  const int nb_coords_to_compute = coords_to_compute.dimension(0);
  const int nb_coords_max = mp_max(nb_coords_to_compute);

  int nb_comp;
  nb_comp = champ_source.line_size();
  assert(resu.line_size() == nb_comp);

  DoubleTab coords(nbproc, 3);
  ArrOfInt flag(nbproc);
  ArrOfDouble resu_partiel(nb_comp);
  DoubleTab resu_partiels(nbproc, nb_comp);

  Calcul_integrale_locale integrale_locale(domaine_source,
                                           *this,
                                           champ_source);

  // Loop over local coordinates x for which we want to compute the integral I(x).
  // If x is near the boundary, the filter function support covers neighbouring processors
  // whose contributions must be included. For each coordinate, all other processors are asked
  // to compute their contribution. We therefore loop over the maximum number of coordinates
  // to synchronize all processes:
  int i, j;
  for (int i_coord = 0; i_coord < nb_coords_max; i_coord++)
    {
      // Each processor sends the coordinate to compute to all other processors:
      if (i_coord < nb_coords_to_compute)
        {
          for (j = 0; j < dim; j++)
            {
              double x = coords_to_compute(i_coord, j);
              for (i = 0; i < nbproc; i++)
                coords(i, j) = x;
            }
          flag = 1;
        }
      else
        {
          flag = 0;
        }
      envoyer_all_to_all(coords, coords);
      envoyer_all_to_all(flag, flag);
      // coord(pe, i) holds the coordinates requested by each processor, and
      // flag[pe] indicates whether that processor has requested a coordinate or has exhausted
      // its coords_to_compute list.
      // We now compute the local processor's contribution to the integrals I(x) for these
      // coordinates. In general, the processor that owns the coordinate does most of the work
      // for it and almost nothing for other coordinates (if the coordinate is far from the
      // local domain, the octree quickly returns an empty list).
      // The workload is therefore approximately balanced across processors.
      for (i = 0; i < nbproc; i++)
        {
          if (flag[i])
            {
              integrale_locale.calculer(coords(i, 0), coords(i, 1), coords(i, 2), resu_partiel);
              for (j = 0; j < nb_comp; j++)
                resu_partiels(i, j) = resu_partiel[j];
            }
        }
      // Send back to each processor the local processor's contribution for the coordinate it requested:
      envoyer_all_to_all(resu_partiels, resu_partiels);
      // The result is the sum of contributions from all processors:
      if (i_coord < nb_coords_to_compute)
        {
          for (j = 0; j < nb_comp; j++)
            {
              double x = 0.;
              for (i = 0; i < nbproc; i++)
                x += resu_partiels(i, j);
              resu(i_coord, j) = x;
            }
        }
    }
}

/*! @brief Computes the convolution product between the filter function and the field "champ_source", which must be discretized at the elements of "domaine_source".
 *
 * @brief The resu array has the same number of columns as "champ_source" and the same number
 *   of rows as coords_to_compute.
 *   The filter function is assumed to have support contained in a cube of half-side box_size
 *   centred on the origin (contributions of elements outside this cube are ignored).
 * @param domaine_source the source discretized domain
 * @param champ_source the source field values array
 * @param coords_to_compute coordinates at which to evaluate the convolution
 * @param resu output array of convolution results
 */
void Moyenne_volumique::calculer_convolution_champ_elem(const Domaine_VF& domaine_source,
                                                        const DoubleTab& champ_source,
                                                        const DoubleTab& coords_to_compute,
                                                        DoubleTab& resu) const
{
  assert(champ_source.dimension(0) == domaine_source.domaine().nb_elem());
  calculer_convolution(domaine_source, champ_source,
                       coords_to_compute, resu);
}

/*! @brief Same as calculer_convolution_champ_elem but for a VDF face field.
 *
 * @brief The source field is assumed to be a vector field containing, for each face,
 *   the normal component of the field at that face.
 *   For each column of the champ_source array, "dimension" columns of resu are filled:
 *   the first using only faces with X-normal, the second with Y-normal faces, etc.
 * @param domaine_source the source discretized domain
 * @param champ_source the source face field values array
 * @param coords_to_compute coordinates at which to evaluate the convolution
 * @param resu output array of convolution results
 */
void Moyenne_volumique::calculer_convolution_champ_face(const Domaine_VF& domaine_source,
                                                        const DoubleTab& champ_source,
                                                        const DoubleTab& coords_to_compute,
                                                        DoubleTab& resu) const
{
  Cerr << " Moyenne_volumique::calculer_convolution_champ_face is not coded" << finl;
  exit();
}
