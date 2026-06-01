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

#include <Perte_Charge_Singuliere.h>
#include <Motcle.h>
#include <Domaine.h>
#include <Domaine_VF.h>
#include <Sous_Domaine.h>
#include <Octree.h>
#include <Source_base.h>
#include <Param.h>
#include <Equation_base.h>
#include <Milieu_base.h>
#include <Discretisation_base.h>
#include <EChaine.h>
#include <Pb_Multiphase.h>
#include <TRUSTTrav.h>

extern void convert_to(const char *s, double& ob);

/*! @brief Reads the specifications of a singular pressure drop from an input stream.
 *
 * @param is an input stream
 * @return the modified input stream
 * @throws unexpected keyword, expected "KX","KY", "KZ" or "K"
 */
// XD perte_charge_singuliere source_base perte_charge_singuliere BRACE Source term that is used to model a pressure
// XD_CONT loss over a surface area (transition through a grid, sudden enlargement) defined by the faces of elements
// XD_CONT located on the intersection of a subzone named subzone_name and a X,Y, or Z plane located at X,Y or Z =
// XD_CONT location.
// XD attr dir chaine(into=["kx","ky","kz","K"]) dir REQ KX, KY or KZ designate directional pressure loss coefficients
// XD_CONT for respectively X, Y or Z direction. Or in the case where you chose a target flow rate with regul. Use K for
// XD_CONT isotropic pressure loss coefficient
// XD attr coeff floattant coeff OPT Value (float) of friction coefficient (KX, KY, KZ).
// XD attr regul bloc_lecture regul OPT option to have adjustable K with flowrate target NL2 { K0 valeur_initiale_de_k
// XD_CONT deb debit_cible eps intervalle_variation_mutiplicatif}.
// XD attr surface bloc_lecture surface REQ Three syntaxes are possible for the surface definition block: NL2 For VDF
// XD_CONT and VEF: { X|Y|Z = location subzone_name } NL2 Only for VEF: { Surface surface_name }. NL2 For PolyMAC_CDO {
// XD_CONT Surface surface_name Orientation champ_uniforme }

Entree& Perte_Charge_Singuliere::lire_regul(Entree& is)
{
  regul_ = 1;
  Nom obsolete(""), alpha_str(""), deb_str;
  Param param("regul");
  param.ajouter("K0", &K_, Param::REQUIRED);
  param.ajouter("deb", &deb_str, Param::REQUIRED);
  param.ajouter("eps", &obsolete);
  param.ajouter("alpha", &alpha_str);
  param.lire_avec_accolades(is);
  if (obsolete != "")
    {
      Cerr << "Syntax error!" << finl;
      Cerr << "The flow rate regulation formula has changed from:" << finl;
      Cerr << "K_ *= min(max(|debit|/deb_cible)^2,(1-eps)^dt),(1+eps)^dt);" << finl;
      Cerr << "to:" << finl;
      Cerr << "K_ += dt*alpha*(debit-debit_cible);" << finl;
      Cerr << "The keyword eps (value between 0 and 1) is no longer valid." << finl;
      Cerr << "The amplitude of the flow rate regulation must be specified by the keyword alpha. If you had a non-zero value for eps, e.g. 0.5, replace with alpha 10" << finl;
      Cerr << "and eps 0 must be replaced by alpha 0" << finl;
      Process::exit();
    }
  if (alpha_str=="") Process::exit("alpha must be specified for the regulation.");
  deb_cible_.setNbVar(1), alpha_.setNbVar(1);
  deb_cible_.setString(deb_str), alpha_.setString(alpha_str);
  deb_cible_.addVar("t"), alpha_.addVar("t");
  deb_cible_.parseString(), alpha_.parseString();
  return is;
}

Entree& Perte_Charge_Singuliere::lire_donnees(Entree& is)
{
  Motcle motlu;
  Motcle acc_ouverte("{");
  Motcles les_motcles(4);
  les_motcles[0] = "KX";
  les_motcles[1] = "KY";
  les_motcles[2] = "KZ";
  les_motcles[3] = "K";
  regul_ = 0;
  is >> motlu;
  if (motlu != acc_ouverte)
    {
      Cerr << "We expected the keyword" << acc_ouverte << " instead of " << motlu << finl;
      Process::exit();
    }
  is >> motlu;
  if (motlu != "dir")
    {
      Cerr << "We expected the keyword dir instead of " << motlu << finl;
      Process::exit();
    }
  else
    {
      is >> motlu;
      int rang = les_motcles.search(motlu);
      if (rang == -1)
        {
          Cerr << "Error while reading the data of Perte_Charge_Singuliere" << finl;
          Cerr << "We expected one of the keywords " << les_motcles << finl;
          Cerr << "instead of " << motlu << finl;
          Process::exit();
        }
      direction_perte_charge_ = rang == 3 ? -1 : rang;

    }
  is >> motlu;
  if (motlu == "regul")
    lire_regul(is);
  else if (motlu == "coeff")
    {
      is >> motlu;
      convert_to(motlu.getChar(), K_);
    }
  else
    {
      Cerr << "We expected the keyword coeff or regul instead of " << motlu << finl;
      Process::exit();
    }

  Cerr << " direction_perte_charge_ " << direction_perte_charge_ << finl;
  Cerr << " pressure drop K_ " << K_ << finl;

  return is;
}

void Perte_Charge_Singuliere::lire_surfaces(Entree& is, const Domaine& le_domaine,
                                            const Domaine_dis_base& domaine_dis, IntVect& les_faces, IntVect& sgn, int lire_derniere_accolade)
{
  const Domaine_VF& zvf = ref_cast(Domaine_VF,domaine_dis);
  const IntTab& elem_faces = zvf.elem_faces();
  const DoubleTab& xv = zvf.xv();
  int nfe = zvf.domaine().nb_faces_elem();
  IntTab face_tab; //1 for faces in the surface
  zvf.creer_tableau_faces(face_tab);
  OWN_PTR(Champ_Don_base) orientation;

  int algo=-1;
  Motcle method;

  Motcle motlu;
  Motcle acc_ouverte("{");
  Motcle acc_fermee("}");

  Nom nom_ss_domaine, nom_surface;
  double position;

  /* Read plan */
  is >> motlu;
  if (motlu != "surface")
    {
      Cerr << "We expected the keyword surface instead of " << motlu << finl;
      Process::exit();
    }
  is >> motlu;
  if (motlu != acc_ouverte)
    {
      Cerr << "We expected the keyword" << acc_ouverte << " instead of " << motlu << finl;
      Process::exit();
    }
  is >> method;
  if (method=="X" || method=="Y" || method=="Z")
    {
      /* Subdomaine algorithm */
      algo=0;
      Nom direction=method, egal;
      if (direction_perte_charge() >= 0 && (((direction=="X") && (direction_perte_charge()!=0))
                                            || ((direction=="Y") && (direction_perte_charge()!=1))
                                            || ((direction=="Z") && (direction_perte_charge()!=2))))
        {
          Nom dir,sect;
          if (direction_perte_charge() == 0)
            {
              dir="KX";
              sect="X";
            }
          else if (direction_perte_charge() == 1)
            {
              dir="KY";
              sect="Y";
            }
          else if (direction_perte_charge() == 2)
            {
              dir="KZ";
              sect="Z";
            }
          Cerr << "The section specified " << direction << " is not coherent with the coefficient pressure loss direction indicated by " << dir << finl;
          Cerr << "Here the section to specify is " << sect << "." << finl;
          Process::exit();

        }
      is >> egal >> position;
      is >> nom_ss_domaine;
      Cerr << " position " << direction << " " << position << finl;
      identifiant_ = nom_ss_domaine + "_" + direction + "=" + Nom(position);
    }
  else if (method=="Surface")
    {
      /* Surface algorithm */
      algo=1;
      /*
            if (zvf.dimension!=3)
              {
                Cerr << "Error in Perte_Charge_Singuliere::lire_surfaces" << finl;
                Cerr << "The keyword " << method << " specified for plan definition is only possible in 3D !" << finl;
                Cerr << "You must used the method of intersection between subdomaine and location of plane." << finl;
                Process::exit();
              }
      */
      if (domaine_dis.que_suis_je().debute_par("Domaine_VDF"))
        {
          Cerr << "Error in Perte_Charge_Singuliere::lire_surfaces" << finl;
          Cerr << "The keyword " << method << " specified for plan definition is only possible in VEF discretization !" << finl;
          Cerr << "You must used the method of intersection between subdomaine and location of plane." << finl;
          Process::exit();
        }
      is >> nom_surface;
      Cerr << " surface " << nom_surface << finl;
      identifiant_ = nom_surface;
    }
  else if (method=="Face_group")
    {
      /* Surface algorithm */
      algo=2;
      is >> nom_surface;
      Cerr << " surface " << nom_surface << finl;
      identifiant_ = nom_surface;
    }
  else
    {
      Cerr << "Error in Perte_Charge_Singuliere::lire_surfaces" << finl;
      Cerr << "The keyword " << method << " specified for plan definition is not coherent !" << finl;
      Process::exit();
    }
  is >> motlu;
  if (motlu == "orientation") is >> orientation, is >> motlu; //read a field to orient the surface

  if (motlu != acc_fermee)
    {
      Cerr << "We expected the keyword" << acc_fermee << " instead of " << motlu << finl;
      Process::exit();
    }
  if (lire_derniere_accolade)
    {
      is >> motlu;
      if (motlu != acc_fermee)
        {
          Cerr << "We expected the keyword" << acc_fermee << " instead of " << motlu << finl;
          Process::exit();
        }
    }
  /* Found faces */
  int compteur=0;

  if (algo==0)
    {
      /* Subdomaine algorithm */
      const Sous_Domaine& ssz = le_domaine.ss_domaine(nom_ss_domaine);
      int coord = method == "X" ? 0 : method == "Y" ? 1 : 2;
      for (int poly=0; poly<ssz.nb_elem_tot(); poly++)
        {
          for (int j=0; j<nfe; j++)
            {
              int numfa = elem_faces(ssz(poly),j);
              // numfa might be negative in case of polyhedron where number of faces for an elem varies (padding)
              if (numfa >= 0 && est_egal(xv(numfa,coord),position) )
                {
                  bool trouve=0;
                  for (int i=0; i<compteur; i++)
                    if (les_faces[i]==numfa)
                      {
                        trouve=1;
                        break;
                      }
                  if (trouve==0) les_faces[compteur++] = numfa, face_tab(numfa) = 1;
                }
            }
        }
    }

  else if (algo==1)
    {
      /* Surface algorithm */
      const Domaine& le_domaine2D = ref_cast(Domaine,le_domaine.interprete().objet(nom_surface));
      const DoubleTab& coord_sommets_2D=le_domaine2D.coord_sommets();

      DoubleTab xv2D;
      le_domaine2D.calculer_centres_gravite(xv2D);
      const int nb_elem_2D = le_domaine2D.nb_elem_tot();
      int nse2D = le_domaine2D.nb_som_elem();

      const DoubleTab& coord_sommets=le_domaine.coord_sommets();
      const IntTab& face_sommets = zvf.face_sommets();
      const OctreeRoot& octree_vol = le_domaine.construit_octree();

      // Loop on faces on the surface domain
      Cerr << " Surface " << nom_surface << " with " << nb_elem_2D << " faces" << finl;
      for (int ind_face=0; ind_face<nb_elem_2D; ind_face++)
        {
          // Fill elem_list_vol by elements surrounding the center of gravity
          ArrOfInt elem_list_vol;
          double zv2d=0;
          if (Objet_U::dimension>2)
            zv2d=xv2D(ind_face,2);
          octree_vol.rang_elems_sommet(elem_list_vol,xv2D(ind_face,0),xv2D(ind_face,1),zv2d);
          int nb_elem_vol = elem_list_vol.size_array();
          // Loop on elements on the volume domain
          for (int ind_elem=0; ind_elem<nb_elem_vol; ind_elem++)
            {
              int elem = elem_list_vol[ind_elem];
              // Loop on faces of each element on the volume domain
              for (int j=0, numfa; j < nfe && (numfa = elem_faces(elem,j)) >= 0; j++)
                {
                  // Test if the same face between surface domain and volum domain
                  int coincide=0;
                  // Loop on vertex of each face of each element on the volume domain

                  int nse_reel = 0;
                  for (int k = 0, numso; k < nse2D && (numso = face_sommets(numfa, k)) >= 0; k++)
                    {
                      nse_reel++;
                      double xcoord_vol=coord_sommets(numso,0);
                      double ycoord_vol=coord_sommets(numso,1);
                      //double zcoord_vol=coord_sommets(numso,2);
                      double  zcoord_vol=0;
                      if (Objet_U::dimension>2)
                        zcoord_vol=coord_sommets(numso,2);

                      for (int i = 0, numso2D; i < nse2D && (numso2D = le_domaine2D.sommet_elem(ind_face, i)) >= 0; i++)
                        {
                          double xcoord_2D=coord_sommets_2D(numso2D,0);
                          double ycoord_2D=coord_sommets_2D(numso2D,1);
                          double zcoord_2D=0;
                          if (Objet_U::dimension>2)
                            zcoord_2D=coord_sommets_2D(numso2D,2);
                          if ( est_egal(xcoord_vol,xcoord_2D)
                               && est_egal(ycoord_vol,ycoord_2D)
                               && est_egal(zcoord_vol,zcoord_2D))
                            {
                              coincide=coincide+1;
                              break;
                            }
                        }
                    }
                  if (coincide == nse_reel)
                    {
                      bool trouve=0;
                      for (int i=0; i<compteur ; i++)
                        if (les_faces[i]==numfa)
                          {
                            trouve=1;
                            break;
                          }
                      if (trouve==0) les_faces[compteur++] = numfa, face_tab(numfa) = 1;
                    }
                }
            }
        }
    }
  else if (algo==2)
    {
      const Groupe_Faces grp_faces = le_domaine.groupe_faces(nom_surface);
      int nb_faces = grp_faces.nb_faces();
      Cerr << "Face group " << nom_surface << " with " << nb_faces << " faces" << finl;
      const ArrOfInt indice_faces = grp_faces.get_indices_faces();
      for (int k=0; k < nb_faces; k++)
        les_faces[k] = indice_faces[k], face_tab(indice_faces[k]) = 1;
      compteur = nb_faces;
    }

  trustIdType faces_found=mp_somme_vect(face_tab);
  if (faces_found==0)
    {
      Cerr << "Error in Perte_Charge_Singuliere::lire_surfaces" << finl;
      Cerr << "No faces has been found !" << finl;
      Process::exit();
    }
  else
    Cerr << " " << faces_found << " faces have been found for the section." << finl;

  les_faces.resize(compteur);
  if (orientation)
    {
      sgn.resize(compteur);
      DoubleTrav xvf(compteur, Objet_U::dimension), ori(compteur, Objet_U::dimension);
      for (int i = 0; i < compteur; i++)
        for (int j = 0; j < Objet_U::dimension; j++) xvf(i, j) = xv(les_faces(i), j);
      orientation->valeur_aux(xvf, ori);
      for (int i = 0; i < compteur; i++)
        {
          double scal = 0;
          for (int j = 0; j < Objet_U::dimension; j++) scal += ori(i, j) * zvf.face_normales(les_faces(i), j);
          sgn(i) = (scal > 0 ? 1 : -1);
        }
    }
}

double Perte_Charge_Singuliere::calculate_Q(const Equation_base& eqn, const IntVect& num_faces, const IntVect& sgn) const
{
  if (sub_type(Pb_Multiphase, eqn.probleme()) && !ref_cast(Pb_Multiphase, eqn.probleme()).equation_masse().has_champ_conserve()) return 0;
  const Domaine_VF& zvf = ref_cast(Domaine_VF, eqn.domaine_dis());
  const DoubleTab& vit = eqn.inconnue().valeurs(),
                   &fac = sub_type(Pb_Multiphase, eqn.probleme()) ? ref_cast(Pb_Multiphase, eqn.probleme()).equation_masse().champ_conserve().passe()
                          : eqn.probleme().get_champ("masse_volumique").valeurs(); // get_champ for flica5 because masse_volumique is not in the milieu...
  const DoubleVect& pf = eqn.milieu().porosite_face(), &fs = zvf.face_surfaces();
  const IntTab& f_e = zvf.face_voisins();
  int cF = fac.dimension_tot(0) == 1, i, n, N = fac.line_size(), d, D = Objet_U::dimension;

  DoubleTrav deb_vect;
  zvf.creer_tableau_faces(deb_vect); //to correctly sum flow rates in parallel
  for (i = 0; i < num_faces.size(); i++)
    {
      int f = num_faces(i), e = f_e(f, f_e(f, 0) < 0); //todo : evaluer fac du bon cote
      if (vit.line_size() > N)
        for (d = 0; d < D; d++)
          for (n = 0; n < N; n++) deb_vect(f) += zvf.face_normales(f, d) * vit(f, N * d + n) * fac(!cF * e, n); //vecteur complet : v . nf
      else for (n = 0; n < N; n++) deb_vect(f) += fs(f) * vit(f, n) * fac(!cF * e, n); //normale aux faces seule
      deb_vect(f) *= pf(f) * (sgn.size() ? sgn(i) : deb_vect(f) > 0 ? 1 : -1) ; //multiply by porosity + orientation (if not defined, count everything as positive)
    }
  return mp_somme_vect(deb_vect);
}

void Perte_Charge_Singuliere::update_K(const Equation_base& eqn, double deb, DoubleVect& bilan)
{
  if (!regul_) return;
  double t = eqn.probleme().schema_temps().temps_courant(), dt = eqn.probleme().schema_temps().pas_de_temps();
  deb_cible_.setVar(0, t), alpha_.setVar(0, t);
  double deb_cible = deb_cible_.eval();
  if (std::abs(deb_cible) > 1e-10)
    {
      const double alpha = alpha_.eval(), error = (deb - deb_cible) / deb_cible;
      K_ += dt * alpha * error;
    }

  //for the monitoring file: only on the master, because Source_base::imprimer() sums over all processes
  if (!Process::me()) bilan(0) = K_, bilan(1) = deb, bilan(2) = deb_cible;
}
