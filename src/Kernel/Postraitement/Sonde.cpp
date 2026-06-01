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

#include <Champ_Generique_Interpolation.h>
#include <Champ_Generique_refChamp.h>
#include <Domaine_Cl_dis_base.h>
#include <Entree_complete.h>

#include <Process.h>
#include <Champ_Inc_base.h>
#include <Postraitement.h>
#include <Domaine_VF.h>
#include <sys/stat.h>
#include <Sonde.h>

Implemente_instanciable_sans_constructeur_ni_destructeur(Sonde,"Sonde",Objet_U);
// XD sonde objet_lecture nul NO_BRACE Keyword is used to define the probes. Observations: the probe coordinates should
// XD_CONT be given in Cartesian coordinates (X, Y, Z), including axisymmetric.
// XD attr nom_sonde chaine nom_sonde REQ Name of the file in which the values taken over time will be saved. The
// XD_CONT complete file name is nom_sonde.son.
// XD attr special chaine(into=["grav","som","nodes","chsom","gravcl"]) special OPT Option to change the positions of
// XD_CONT the probes. Several options are available: NL2 grav : each probe is moved to the nearest cell center of the
// XD_CONT mesh; NL2 som : each probe is moved to the nearest vertex of the mesh NL2 nodes : each probe is moved to the
// XD_CONT nearest face center of the mesh; NL2 chsom : only available for P1NC sampled field. The values of the probes
// XD_CONT are calculated according to P1-Conform corresponding field. NL2 gravcl : Extend to the domain face boundary a
// XD_CONT cell-located segment probe in order to have the boundary condition for the field. For this type the extreme
// XD_CONT probe point has to be on the face center of gravity.
// XD attr nom_inco chaine nom_inco REQ Name of the sampled field.
// XD attr mperiode chaine(into=["periode"]) mperiode REQ Keyword to set the sampled field measurement frequency.
// XD attr prd floattant prd REQ Period value. Every prd seconds, the field value calculated at the previous time step
// XD_CONT is written to the nom_sonde.son file.
// XD attr type sonde_base type REQ Type of probe.

// XD sonde_base objet_lecture sonde_base NO_BRACE Basic probe. Probes refer to sensors that allow a value or several
// XD_CONT points of the domain to be monitored over time. The probes may be a set of points defined one by one (keyword
// XD_CONT Points) or a set of points evenly distributed over a straight segment (keyword Segment) or arranged according
// XD_CONT to a layout (keyword Plan) or according to a parallelepiped (keyword Volume). The fields allow all the values
// XD_CONT of a physical value on the domain to be known at several moments in time.

// XD segmentfacesx sonde_base segmentfacesx NO_BRACE Segment probe where points are moved to the nearest x faces
// XD attr nbr entier nbr REQ Number of probe points of the segment, evenly distributed.
// XD attr point_deb un_point point_deb REQ First outer probe segment point.
// XD attr point_fin un_point point_fin REQ Second outer probe segment point.

// XD segmentfacesy sonde_base segmentfacesy NO_BRACE Segment probe where points are moved to the nearest y faces
// XD attr nbr entier nbr REQ Number of probe points of the segment, evenly distributed.
// XD attr point_deb un_point point_deb REQ First outer probe segment point.
// XD attr point_fin un_point point_fin REQ Second outer probe segment point.

// XD segmentfacesz sonde_base segmentfacesz NO_BRACE Segment probe where points are moved to the nearest z faces
// XD attr nbr entier nbr REQ Number of probe points of the segment, evenly distributed.
// XD attr point_deb un_point point_deb REQ First outer probe segment point.
// XD attr point_fin un_point point_fin REQ Second outer probe segment point.

// XD radius sonde_base radius NO_BRACE not_set
// XD attr nbr entier nbr REQ Number of probe points of the segment, evenly distributed.
// XD attr point_deb un_point point_deb REQ First outer probe segment point.
// XD attr radius floattant radius REQ not_set
// XD attr teta1 floattant teta1 REQ not_set
// XD attr teta2 floattant teta2 REQ not_set

// XD un_point objet_lecture nul NO_BRACE A point.
// XD attr pos listf pos REQ Point coordinates.

// XD listpoints listobj nul NO_BRACE un_point NO_COMMA Points.
// XD points sonde_base points NO_BRACE Keyword to define the number of probe points. The file is arranged in columns.
// XD attr points listpoints points REQ Probe points.

// XD numero_elem_sur_maitre sonde_base numero_elem_sur_maitre NO_BRACE Keyword to define a probe at the special
// XD_CONT element. Useful for min/max sonde.
// XD attr numero entier numero REQ element number

// XD segmentpoints points segmentpoints NO_BRACE This keyword is used to define a probe segment from specifics points.
// XD_CONT The nom_champ field is sampled at ns specifics points.

// XD position_like sonde_base position_like NO_BRACE Keyword to define a probe at the same position of another probe
// XD_CONT named autre_sonde.
// XD attr autre_sonde chaine autre_sonde REQ Name of the other probe.

// XD plan sonde_base plan NO_BRACE Keyword to set the number of probe layout points. The file format is type .lml
// XD attr nbr entier nbr REQ Number of probes in the first direction.
// XD attr nbr2 entier nbr2 REQ Number of probes in the second direction.
// XD attr point_deb un_point point_deb REQ First point defining the angle. This angle should be positive.
// XD attr point_fin un_point point_fin REQ Second point defining the angle. This angle should be positive.
// XD attr point_fin_2 un_point point_fin_2 REQ Third point defining the angle. This angle should be positive.

// XD volume sonde_base volume NO_BRACE Keyword to define the probe volume in a parallelepiped passing through 4 points
// XD_CONT and the number of probes in each direction.
// XD attr nbr entier nbr REQ Number of probes in the first direction.
// XD attr nbr2 entier nbr2 REQ Number of probes in the second direction.
// XD attr nbr3 entier nbr3 REQ Number of probes in the third direction.
// XD attr point_deb un_point point_deb REQ Point of origin.
// XD attr point_fin un_point point_fin REQ Point defining the first direction (from point of origin).
// XD attr point_fin_2 un_point point_fin_2 REQ Point defining the second direction (from point of origin).
// XD attr point_fin_3 un_point point_fin_3 REQ Point defining the third direction (from point of origin).

// XD circle sonde_base circle NO_BRACE Keyword to define several probes located on a circle.
// XD attr nbr entier nbr REQ Number of probes between teta1 and teta2 (angles given in degrees).
// XD attr point_deb un_point point_deb REQ Center of the circle.
// XD attr direction entier(into=[0,1,2]) direction OPT Axis normal to the circle plane (0:x axis, 1:y axis, 2:z axis).
// XD attr radius floattant radius REQ Radius of the circle.
// XD attr theta1 floattant theta1 REQ First angle.
// XD attr theta2 floattant theta2 REQ Second angle.

// XD circle_3 sonde_base circle_3 NO_BRACE Keyword to define several probes located on a circle (in 3-D space).
// XD attr nbr entier nbr REQ Number of probes between teta1 and teta2 (angles given in degrees).
// XD attr point_deb un_point point_deb REQ Center of the circle.
// XD attr direction entier(into=[0,1,2]) direction REQ Axis normal to the circle plane (0:x axis, 1:y axis, 2:z axis).
// XD attr radius floattant radius REQ Radius of the circle.
// XD attr theta1 floattant theta1 REQ First angle.
// XD attr theta2 floattant theta2 REQ Second angle.

static int fichier_sondes_cree=0;
static SFichier fichier_sondes;

/*! @brief Constructor for a probe, given its name.
 *
 * @param nom the name of the probe to construct
 */
Sonde::Sonde(const Nom& nom)  :
  nom_(nom),
  dim(-1),
  ncomp(-1),
  numero_elem_(-1),
  periode(1.e10),   // default period initialization
  nb_bip(0.),
  orientation_faces_(-1)
{}

/*! @brief Default constructor for a probe.
 *
 */
Sonde::Sonde() :
  Sonde(Nom())    // thank you C++11
{}

/*! @brief Prints the type of the object to an output stream.
 *
 * @param (Sortie& s) an output stream
 * @return (Sortie&) the modified output stream
 */
Sortie& Sonde::printOn(Sortie& s ) const
{
  return s << que_suis_je();
}

void Sonde::completer()
{
  // Recherche du champ sonde
  // Remplissage de la reference au champ
  Cerr << "Reading of the probe: " << nom_ << " on the field: " << nom_champ_lu_ << finl;
  //We should access the domain through the generic field
  //But the reference has not been made yet
  Probleme_base& Pb = mon_post->probleme();
  const Motcle nom_domaine = mon_post->domaine()->le_nom();
  Motcle nom_champ_ref;
  Motcle nom_macro=nom_champ_lu_;
  Motcle nom_macro_test, nom1("_not_def"), nom2("_not_def");
  const Motcles& noms_champs_postraitables = mon_post->les_sondes().get_noms_champs_postraitables();
  if (!(noms_champs_postraitables.contient_(nom_champ_lu_)||mon_post->comprend_champ_post(nom_champ_lu_)))
    {
      if (nom_champ_lu_.debute_par("MOYENNE_"))
        nom_macro = nom_macro.suffix("MOYENNE_");
      else if (nom_champ_lu_.debute_par("ECART_TYPE_"))
        nom_macro = nom_macro.suffix("ECART_TYPE_");

      if (nom_champ_lu_.debute_par("CORRELATION_"))
        {
          nom_macro_test = nom_champ_lu_;
          nom_macro_test = nom_macro_test.suffix("CORRELATION_");
          for (auto& itr : noms_champs_postraitables)
            {
              nom_macro_test = nom_macro_test.suffix(itr);
              if (nom_macro_test.debute_par("_"))
                {
                  nom1 = itr;
                  nom2 = nom_macro_test.suffix("_");
                }
            }
        }
    }

  if (noms_champs_postraitables.contient_(nom_macro))
    {
      if (nom_macro==Motcle(nom_champ_lu_))
        {
          OBS_PTR(Champ_base) champ_ref = Pb.get_champ(nom_champ_lu_);
          const Nom& nom_champ_base = champ_ref->le_nom();
          const Noms& compos_base = champ_ref->noms_compo();
          int ncomposante = Champ_Generique_base::composante(nom_champ_lu_,nom_champ_base,compos_base,champ_ref->get_synonyms());

          if (ncomposante==-1)
            nom_champ_ref = nom_champ_base +"_natif_"+ nom_domaine;
          else
            {
              Nom nume(ncomposante);
              nom_champ_ref = nom_champ_base+"_natif_"+ nom_domaine+nume;
            }
        }
      else
        nom_champ_ref = nom_champ_lu_+"_natif_"+nom_domaine;
    }
  else
    {
      //Case of correlations
      if ((noms_champs_postraitables.contient_(nom1)) && (noms_champs_postraitables.contient_(nom2)))
        {
          nom_champ_ref = "Correlation_";
          nom_champ_ref += nom1+"_natif_"+nom_domaine+"_"+nom2+"_natif_"+nom_domaine;
        }
      else
        //Case of a general generic field
        nom_champ_ref = nom_champ_lu_;
    }

  mon_champ = mon_post->get_champ_post(nom_champ_ref);

  // Filling the ncomp attribute (its default value is -1)

  const Noms nom_champ = mon_champ->get_property("nom");
  const Noms noms_comp = mon_champ->get_property("composantes");
  ncomp = Champ_Generique_base::composante(nom_champ_ref,nom_champ[0],noms_comp,mon_champ->get_property("synonyms"));

  initialiser();
}
/*! @brief Reads the specifications of a probe from an input stream.
 *
 * Format:
 *     Sondes
 *      {
 *        [nom_sonde nom_champ Periode dts Points n x1 y1 [z1] ... xn yn [zn]
 *        [nom_sonde nom_champ Periode dts Segment ns x1 y1 [z1] x2 y2 [z2]
 *        ...
 *      }
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the modified input stream
 * @throws probe data not defined
 * @throws format error, unknown keyword
 * @throws probe data not correctly defined
 */
Entree& Sonde::readOn(Entree& is)
{
  assert(mon_post);
  Motcle motlu;
  //Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  int nbre_points;

  is >> motlu;
  // Reading a keyword that is not the field
  if (motlu=="nodes")
    {
      nodes = true;
      is >> motlu;
    }
  else if (motlu=="chsom")
    {
      chsom = true;
      is >> motlu;
    }
  else if (motlu=="grav")
    {
      grav = true;
      is >> motlu;
    }
  else if (motlu=="gravcl")
    {
      gravcl = true;
      is >> motlu;
    }
  else if (motlu=="som")
    {
      som = true;
      is >> motlu;
    }
  // Assign the field name
  nom_champ_lu_ = motlu;

  validate_position();
  create_champ_generique(is, motlu);

  // Reading the probe characteristics
  IntVect fait(2);

  Motcles les_motcles(16);
  {
    les_motcles[0] = "periode";
    les_motcles[1] = "point";
    les_motcles[2] = "points";
    les_motcles[3] = "segment";
    les_motcles[4] = "plan";
    les_motcles[5] = "volume";
    les_motcles[6] = "segmentxdx";
    les_motcles[7] = "planxdxdy";
    les_motcles[8] = "circle";
    les_motcles[9] = "position_like";
    les_motcles[10] = "numero_elem_sur_maitre";
    les_motcles[11] = "segmentpoints";
    les_motcles[12] = "segmentfacesx";
    les_motcles[13] = "segmentfacesy";
    les_motcles[14] = "segmentfacesz";
    les_motcles[15] = "radius";
  }

  while ((fait(0) != 1) || (fait(1) != 1))
    {
      is >> motlu;
      if (motlu == accolade_fermee)
        {
          Cerr << "Error while reading the probe " << nom_ <<finl;
          Cerr << "The data of the probe were not defined" << finl;
          exit();
        }
      int rang=les_motcles.search(motlu);
      if (rang == -1)
        {
          Cerr << "Error while reading the probe " << nom_ <<finl;
          Cerr << motlu << " is not understood; the keywords understood are : " << finl;
          Cerr << les_motcles;
          exit();
        }

      if(rang != 0)
        validate_type(motlu); // restrict valid localisation for IJK notably

      switch(rang)
        {
        case 0:
          {
            is >> periode;
            break;
          }
        case 1:
        case 2:
        case 11:
          {
            type_ = les_motcles[rang];
            rang = 1;
            dim = 0;
            is >> nbre_points;
            les_positions_sondes_initiales_.resize(nbre_points,dimension);

            for (int i=0; i<nbre_points; i++)
              for (int j=0; j<dimension; j++)
                is >> les_positions_sondes_initiales_(i,j);

            break;
          }
        case 10:
          {
            type_ = les_motcles[rang];
            rang=1;
            dim=0;
            gravcl = false;
            is >> numero_elem_;
            les_positions_sondes_initiales_.resize(1,dimension);
            break;
          }
        case 3:
        case 6:
        case 12:
        case 13:
        case 14:
          {
            type_ = les_motcles[rang];
            int rang2=rang;
            rang = 1;
            dim = 1;
            if (rang2==12 || rang2==13 || rang2==14)
              {
                nodes = true;
                grav = false;
                gravcl = false;
                som = false;
              }
            const Motcle mot = type_.getSuffix("SEGMENTFACES");
            Motcles nom_dir(3);

            nom_dir[0] = "X";
            nom_dir[1] = "Y";
            nom_dir[2] = "Z";

            orientation_faces_ = nom_dir.search(mot);

            DoubleVect origine(dimension);
            DoubleVect extremite(dimension);
            DoubleVect dx(dimension);
            is >> nbre_points;
            les_positions_sondes_initiales_.resize(nbre_points,dimension);

            for (int i=0; i<dimension; i++)
              {
                is >> origine(i);
                type_+=" ";
                type_+=(Nom)origine(i);
              }
            for (int i=0; i<dimension; i++)
              {
                is >> extremite(i);
                type_+=" ";
                type_+=(Nom)extremite(i);
              }
            for (int i=0; i<dimension; i++)
              if (rang2==6)
                dx(i)=(extremite(i))/(nbre_points-1);
              else
                dx(i)=(extremite(i)-origine(i))/(nbre_points-1);
            for (int i=0; i<nbre_points; i++)
              for (int j=0; j<dimension; j++)
                les_positions_sondes_initiales_(i,j)=origine(j)+i*dx(j);
            break;
          }
        case 4:
        case 7:
          {
            type_ = les_motcles[rang];
            int rang2=rang;
            rang = 1;
            dim = 2;
            gravcl = false;
            DoubleVect origine(dimension);
            DoubleVect extremite1(dimension);
            DoubleVect extremite2(dimension);
            DoubleVect dx1(dimension);
            DoubleVect dx2(dimension);
            is >> nbre_points1;
            is >> nbre_points2;
            nbre_points=nbre_points1*nbre_points2;
            les_positions_sondes_initiales_.resize(nbre_points,dimension);

            for (int i=0; i<dimension; i++)
              is >> origine(i);
            for (int i=0; i<dimension; i++)
              is >> extremite1(i);
            for (int i=0; i<dimension; i++)
              is >> extremite2(i);
            if (rang2==7)
              {
                for (int i=0; i<dimension; i++)
                  dx1(i)=(extremite1(i))/(nbre_points1-1);
                for (int i=0; i<dimension; i++)
                  dx2(i)=(extremite2(i))/(nbre_points2-1);
              }
            else
              {
                for (int i=0; i<dimension; i++)
                  dx1(i)=(extremite1(i)-origine(i))/(nbre_points1-1);
                for (int i=0; i<dimension; i++)
                  dx2(i)=(extremite2(i)-origine(i))/(nbre_points2-1);
              }
            for (int i=0; i<nbre_points1; i++)
              for (int j=0; j<nbre_points2; j++)
                for (int k=0; k<dimension; k++)
                  les_positions_sondes_initiales_(i*nbre_points2+j,k)=origine(k)+i*dx1(k)+j*dx2(k);
            break;
          }

        case 5:
          {
            type_ = les_motcles[rang];
            rang = 1;
            dim = 3;
            gravcl = false;
            ArrOfDouble origine(dimension);
            ArrOfDouble extremite1(dimension);
            ArrOfDouble extremite2(dimension);
            ArrOfDouble extremite3(dimension);
            ArrOfDouble dx1(dimension);
            ArrOfDouble dx2(dimension);
            ArrOfDouble dx3(dimension);
            is >> nbre_points1;
            is >> nbre_points2;
            is >> nbre_points3;
            nbre_points=nbre_points1*nbre_points2*nbre_points3;
            les_positions_sondes_initiales_.resize(nbre_points,dimension);

            for (int i=0; i<dimension; i++)
              is >> origine[i];
            for (int i=0; i<dimension; i++)
              is >> extremite1[i];
            for (int i=0; i<dimension; i++)
              is >> extremite2[i];
            for (int i=0; i<dimension; i++)
              is >> extremite3[i];
            for (int i=0; i<dimension; i++)
              dx1[i]=(extremite1[i]-origine[i])/(nbre_points1-1);
            for (int i=0; i<dimension; i++)
              dx2[i]=(extremite2[i]-origine[i])/(nbre_points2-1);
            for (int i=0; i<dimension; i++)
              dx3[i]=(extremite3[i]-origine[i])/(nbre_points3-1);
            for (int i=0; i<nbre_points1; i++)
              for (int j=0; j<nbre_points2; j++)
                for (int m=0; m<nbre_points3; m++)
                  for (int k=0; k<dimension; k++)
                    les_positions_sondes_initiales_(i+j*nbre_points1+m*nbre_points1*nbre_points2,k)=origine[k]+i*dx1[k]+j*dx2[k]+m*dx3[k];
            break;
          }
        case 8:
          {
            type_ = les_motcles[rang];
            // circle nbre_points x0 y0 [z0 dir] radius teta1 teta2
            rang = 1;
            dim = 1;
            gravcl = false;
            int dir;
            double radius, teta1, teta2;
            DoubleVect origine(dimension);
            is >> nbre_points;
            les_positions_sondes_initiales_.resize(nbre_points,dimension);
            for (int i=0; i<dimension; i++)
              is >> origine(i);
            if (dimension==3) is >> dir;
            is >> radius >> teta1 >> teta2;
            // Adding information
            for (int i=0; i<dimension; i++)
              {
                type_+=" ";
                type_+=(Nom)origine(i);
              }
            type_+=" ";
            type_+=(Nom)radius;
            type_+=" ";
            type_+=(Nom)teta1;
            type_+=" ";
            type_+=(Nom)teta2;
            // We calculate the positions
            for (int i=0; i<nbre_points; i++)
              {
                double angle=teta1+(teta2-teta1)*i/(nbre_points-1);
                angle*=M_PI/180;
                if (dimension==2)
                  {
                    les_positions_sondes_initiales_(i,0)=origine(0)+radius*cos(angle);
                    les_positions_sondes_initiales_(i,1)=origine(1)+radius*sin(angle);
                  }
                else if (dimension==3)
                  {
                    if (dir==0)
                      {
                        les_positions_sondes_initiales_(i,0)=origine(0);
                        les_positions_sondes_initiales_(i,1)=origine(1)+radius*cos(angle);
                        les_positions_sondes_initiales_(i,2)=origine(2)+radius*sin(angle);
                      }
                    else if (dir==1)
                      {
                        les_positions_sondes_initiales_(i,0)=origine(0)+radius*cos(angle);
                        les_positions_sondes_initiales_(i,1)=origine(1);
                        les_positions_sondes_initiales_(i,2)=origine(2)+radius*sin(angle);
                      }
                    else if (dir==2)
                      {
                        les_positions_sondes_initiales_(i,0)=origine(0)+radius*cos(angle);
                        les_positions_sondes_initiales_(i,1)=origine(1)+radius*sin(angle);
                        les_positions_sondes_initiales_(i,2)=origine(2);
                      }
                  }
              }
            break;
          }
        case 15:
          {
            assert(dimension==3);
            type_ = les_motcles[rang];
            // radius nbre_points x0 y0 z0 theta radius1 radius2
            rang = 1;
            dim = 1;
            gravcl = false;
            double theta, radius1, radius2;
            DoubleVect origine(dimension);
            is >> nbre_points;
            les_positions_sondes_initiales_.resize(nbre_points,dimension);
            for (int i=0; i<dimension; i++)
              is >> origine(i);
            is >> theta >> radius1 >> radius2;
            // Adding information
            for (int i=0; i<dimension; i++)
              {
                type_+=" ";
                type_+=(Nom)origine(i);
              }
            type_+=" ";
            type_+=(Nom)theta;
            type_+=" ";
            type_+=(Nom)radius1;
            type_+=" ";
            type_+=(Nom)radius2;
            // We calculate the positions
            theta*=M_PI/180;
            for (int i=0; i<nbre_points; i++)
              {
                double radius = radius1+(radius2-radius1)*i/(nbre_points-1);
                les_positions_sondes_initiales_(i,0)=origine(0)+radius*cos(theta);
                les_positions_sondes_initiales_(i,1)=origine(1)+radius*sin(theta);
                les_positions_sondes_initiales_(i,2)=origine(2);
              }
            break;
          }
        case 9:
          {
            Motcle autre_sonde;
            is >> autre_sonde;
            // on cherche la sonde correspondante
            int m=-1;
            const Sondes& les_sondes=mon_post->les_sondes();
            for (int i=0; i<les_sondes.size(); i++)
              if (les_sondes(i)->get_nom()==autre_sonde)
                {
                  m=i;
                  break;
                }
            //else Cerr<<les_sondes(i).get_nom()<<finl;
            if (m==-1)
              {
                Cerr<<" The probe name "<<autre_sonde<< " was not found"<<finl;
                exit();
              }
            // on recupere  les_positions_
            const Sonde& la_sonde_ref=les_sondes(m);
            type_ = la_sonde_ref.get_type();
            les_positions_sondes_initiales_=la_sonde_ref.les_positions_sondes_initiales();
            dim =  la_sonde_ref.get_dim();
            rang=1;
            break;
          }
        default:
          {
            Cerr << motlu <<"is not yet understood!" << finl;
            exit();
          }
        }
      fait(rang) = 1;
    }
  if ( (fait[0] == 0) || (fait[1] == 0) || (dim==-1))
    {
      Cerr << "Error while reading the probe " << nom_ << finl;
      Cerr << "The data of the probe have not been properly defined" << finl;
      exit();
    }

  // Build the file name associated with the probe
  nom_fichier_=nom_du_cas();
  nom_fichier_+= "_";
  nom_fichier_+= nom_;
  if (dim==2)
    nom_fichier_+= ".plan";
  else
    nom_fichier_+= ".son";

  if (Process::je_suis_maitre())
    {
      // Adding the probe file name to the file listing probes
      // This file will be used by Run_sonde
      // SFichier fichier_sondes;
      Nom nom_fich = nom_du_cas();
      nom_fich += ".sons";
      if (!fichier_sondes.is_open())
        {
          if (!fichier_sondes_cree)
            fichier_sondes.ouvrir(nom_fich);
          else
            fichier_sondes.ouvrir(nom_fich, ios::app);
        }
      fichier_sondes << nom_fichier_ << finl;
      fichier_sondes_cree = 1;
    }
  return is;
}

/** Creation of Champ_Generique_refChamp objects needed for initializing the REF to Champ_Generique_base.
 * If the requested field is a Champ_base known to the problem, the corresponding Champ_Generique_refChamp is created.
 */
void Sonde::create_champ_generique(Entree& is, const Motcle& motlu)
{
  Probleme_base& Pb = mon_post->probleme();
  const Motcles& noms_champs_postraitables = mon_post->les_sondes().get_noms_champs_postraitables();
  if (noms_champs_postraitables.contient_(nom_champ_lu_))
    {
      Pb.creer_champ(nom_champ_lu_);
      //We will create a Champ_Generique_refChamp whose name is based on
      //the name of the field it references, not a component of that field
      OBS_PTR(Champ_base) champ_ref = Pb.get_champ(nom_champ_lu_);
      const Nom& le_nom_champ = champ_ref->le_nom();
      const Motcle nom_domaine = mon_post->domaine()->le_nom();
      Motcle identifiant;

      identifiant = Motcle(le_nom_champ)+"_natif_"+nom_domaine;
      if (!mon_post->comprend_champ_post(identifiant))
        {
          mon_post->creer_champ_post(le_nom_champ,"natif",is);
        }
    }
  else
    {
      Nom expression;
      int is_champ_predefini = Pb.expression_predefini(motlu,expression);
      if ((is_champ_predefini) && (!mon_post->comprend_champ_post(nom_champ_lu_)))
        {
          OWN_PTR(Champ_Generique_base) champ;
          Entree_complete s_complete(expression,is);
          s_complete>>champ;
          mon_post->complete_champ(champ,nom_champ_lu_);
        }
    }
}

/*! @brief Associates the post-processing object with the probe.
 *
 * @param (Postraitement& le_post) the post-processing object to associate
 */
void Sonde::associer_post(const Postraitement& le_post)
{
  mon_post=le_post;
  if (mon_post->noms_fichiers_sondes().contient(nom_)==0)
    mon_post->noms_fichiers_sondes().add(nom_);
  else
    {
      Cerr << "The filename " << nom_ << " is used several times for probes." << finl;
      Cerr << "A probe must have a unique filename to avoid writing conflicts." << finl;
      exit();
    }
}

const Domaine& Sonde::get_domaine_geom() const
{
  return mon_champ->get_ref_domain();
}

const Noms Sonde::get_noms_champ() const
{
  return mon_champ->get_property("nom");
}

int Sonde::get_nb_compo_champ() const
{
  const Noms noms_comp = mon_champ->get_property("composantes");
  Motcle directive = mon_champ->get_directive_pour_discr();
  int nb_comp = (directive == "champ_Fonc_Quad_DG" || directive == "champ_elem_DG") ? 1 : noms_comp.size();
  return nb_comp;
}

double Sonde::get_temps_champ() const
{
  // TODO ABN : when is it different from Pb::temps_courant ?? Inc fields maybe ?
  return mon_champ->get_time();
}

void Sonde::fix_probe_position_grav()
{
  int nbre_points_tot = les_positions_sondes_.dimension(0);
  DoubleTab coords_bords(2,dimension);
  for (int idim=0; idim<dimension; idim++)
    {
      coords_bords(0,idim) = les_positions_sondes_(0,idim);
      coords_bords(1,idim) = les_positions_sondes_(nbre_points_tot-1,idim);
    }

  Cerr<<"The location of the probe named "<<nom_<<" are modified (to centers of gravity). Check the .log files to see the new location."<<finl;
  const Domaine_VF& domaineVF = ref_cast(Domaine_VF,mon_champ->get_ref_domaine_dis_base());
  const DoubleTab& xp = domaineVF.xp();
  for (int i=0; i<nbre_points_tot; i++)
    {
      if(elem_[i]!=-1)
        {
          Journal()<<"The point " << i << " of the probe "<<nom_<<" is moved:";
          for (int dir=0; dir<dimension; dir++)
            {
              Journal() << " x(" << dir << "): " << les_positions_sondes_(i,dir) << " -> " << xp(elem_[i],dir);
              les_positions_sondes_(i,dir)=xp(elem_[i],dir);
            }
          Journal() << finl;
        }
    }
  if (gravcl)
    ajouter_bords(coords_bords);
}

/*! @brief Initialises the probe.
 *
 * Sizes the value arrays, and checks that the specified points lie within the computation domain.
 *
 * @param domaine_geom the computation domain to be probed
 * @throws probe point outside the computation domain
 */
void Sonde::initialiser()
{
  // Compute exact positions from the initial probe positions:
  les_positions_sondes_ = les_positions_sondes_initiales_;
  // Dimension the elem_ array:
  int nbre_points_tot = les_positions_sondes_.dimension(0);
  if(elem_.size() != nbre_points_tot)
    elem_.resize(nbre_points_tot);

  const Domaine& domaine_geom = get_domaine_geom();

  if (numero_elem_==-1)
    {
      // Location of probes is given by coordinates in the les_positions_ array:
      int nb_coord = les_positions_sondes_.dimension(1);
      if (nb_coord != Objet_U::dimension)
        {
          Cerr << "You can't specify the probe named " << nom_ << " with "<< nb_coord << " coordinates on the domain named " <<domaine_geom.le_nom()<<finl;
          Cerr << "which has spatial dimension " << Objet_U::dimension << finl;
          Cerr << "Change the probe coordinates or use numero_elem_sur_maitre keyword (see documentation)" << finl;
          Cerr << "to specify a cell containing the probe and not its coordinates." << finl;
          // [ABN] : we should exit, otherwise we just don't see the warning:
          Process::exit();
        }
      // Fill the elem_ array (which list real cells containing all the probes):
      domaine_geom.chercher_elements(les_positions_sondes_,elem_,1);
    }
  else
    {
      // Location of the probe is given by an element number:
      elem_[0]= numero_elem_;
      const IntTab& les_elems=mon_champ->get_ref_domaine_dis_base().domaine().les_elems();
      if (numero_elem_<les_elems.dimension_tot(0))
        {
          const DoubleTab& coord=mon_champ->get_ref_domaine_dis_base().domaine().les_sommets();
          int nb_som=les_elems.dimension(1);
          // Fill les_positions_ with the cog of the cell numero_elem_
          for (int s=0; s<nb_som; s++)
            {
              int soml=les_elems(numero_elem_,s);
              if (soml>-1)
                for (int dir=0; dir<dimension; dir++)
                  les_positions_sondes_(0,dir)+=coord(soml,dir)/nb_som;
            }
        }
      else
        {
          if (je_suis_maitre())
            {
              Cerr<<" On the probe named " << nom_ << " , the element number "<<numero_elem_<<" does not exist on the master processor, we put the position to zero"<<finl;
              les_positions_sondes_=0;
            }
          elem_[0]=-1;
        }
    }
  // Check if some probes are outside the domain:
  ArrOfDouble tmp(nbre_points_tot);
  for (int i=0; i<nbre_points_tot; i++)
    tmp[i] = elem_[i];
  mp_max_for_each_item(tmp);
  for (int i=0; i<nbre_points_tot; i++)
    if (tmp[i]==-1)
      Cerr << "WARNING: The point number " << i+1 << " of the probe named " << nom_ << " is outside the computational domain " << domaine_geom.le_nom() << finl;

  const Noms noms_champ = get_noms_champ();

  if (grav || gravcl)
    fix_probe_position_grav();
  else if (nodes)
    {
      const Domaine_VF& domaineVF = ref_cast(Domaine_VF,mon_champ->get_ref_domaine_dis_base());
      const DoubleTab& xv = domaineVF.xv();
      const IntTab& elem_faces = domaineVF.elem_faces();
      if (mp_max(elem_faces.size_array())==0)
        {
          Cerr << "Error: the domain " << domaineVF.domaine().le_nom() << " is not discretized." << finl;
          exit();
        }
      if (sub_type(Champ_Generique_Interpolation,mon_champ.valeur()))
        {
          Motcle dom_interp=mon_champ->get_ref_domain().le_nom();
          Cerr << finl;
          Cerr << "Error in your probe : " << nom_ << finl;
          Cerr << "You can not project to nodes, the field " << noms_champ[0] << finl;
          Cerr << "which is interpolated on the domain " << dom_interp << finl;
          exit();
        }
      Cerr<<"The location of the probe named "<<nom_<<" are modified (to faces). Check the .log files to see the new location."<<finl;
      const int nfaces_par_element = domaine_geom.nb_faces_elem() ;
      for (int i=0; i<nbre_points_tot; i++)
        {
          double dist_min=DMAXFLOAT;
          int face_min=-1;
          if(elem_[i]!=-1)
            {
              Journal()<<"The point " << i << " of the probe "<<nom_<<" is moved:";
              for(int fac=0; fac<nfaces_par_element; fac++)
                {
                  int face=elem_faces(elem_[i],fac);
                  // segmentfaces probe type: search only among faces oriented along
                  // orientation_faces_
                  if (face >= 0 &&
                      (orientation_faces_ == -1 || domaineVF.orientation_si_definie(face)==orientation_faces_))
                    {
                      double dist=0.;

                      for (int dir=0; dir<dimension; dir++)
                        dist+=(xv(face,dir)-les_positions_sondes_(i,dir))*(xv(face,dir)-les_positions_sondes_(i,dir));

                      if(dist<=dist_min)
                        {
                          dist_min=dist;
                          face_min=face;
                        }
                    }
                }

              for (int dir=0; dir<dimension; dir++)
                {
                  Journal() << " x(" << dir << "): " << les_positions_sondes_(i,dir) << " -> " << xv(face_min,dir);
                  les_positions_sondes_(i,dir)=xv(face_min,dir);
                }
              Journal() << " (" << (face_min < domaineVF.premiere_face_int() ? "boundary face " : "internal face ") << face_min<<")" << finl;
            }
        }
    }
  else if (som)
    {
      if (sub_type(Champ_Generique_Interpolation,mon_champ.valeur()))
        {
          Motcle dom_interp=mon_champ->get_ref_domain().le_nom();
          Cerr << finl;
          Cerr << "Error in your probe : " << nom_ << finl;
          Cerr << "You can not project to vertexes, the field " << noms_champ[0] << finl;
          Cerr << "which is interpolated on the domain " << dom_interp << finl;
          exit();
        }
      Cerr<<"The location of the probe named "<<nom_<<" are modified (to vertexes). Check the .log files to see the new location."<<finl;
      const IntTab& sommet_elem = domaine_geom.les_elems();
      const int sommets_par_element = domaine_geom.les_elems().dimension(1);
      const DoubleTab& coord = domaine_geom.les_sommets();
      for (int i=0; i<nbre_points_tot; i++)
        {
          double dist_min=DMAXFLOAT;
          int sommet_min=-1;
          if(elem_[i]!=-1)
            {
              Journal()<<"The point " << i << " of the probe "<<nom_<<" is moved:";
              for(int isom=0, sommet; isom<sommets_par_element && (sommet = sommet_elem(elem_[i],isom)) >= 0; isom++)
                {
                  double dist=0.;
                  for (int dir=0; dir<dimension; dir++)
                    dist+=(coord(sommet,dir)-les_positions_sondes_(i,dir))*(coord(sommet,dir)-les_positions_sondes_(i,dir));

                  if(dist<=dist_min)
                    {
                      dist_min=dist;
                      sommet_min=sommet;
                    }
                }

              for (int dir=0; dir<dimension; dir++)
                {
                  Journal() << " x(" << dir << "): " << les_positions_sondes_(i,dir) << " -> " << coord(sommet_min,dir);
                  les_positions_sondes_(i,dir)=coord(sommet_min,dir);
                }
              Journal() << finl;
            }
        }
    }

  // See if some positions should be modified (to elem with "grav" for instance)
  // Should be done before supprime_doublons, because it may create duplicated values.
  fix_probe_position();

  bool supprime_doublons = true; // New in 1.8.4 (uniqueness of probe points)
  if (mon_post->DeprecatedKeepDuplicatedProbes) supprime_doublons = false; // option to keep as P-E needs to preserve duplicated probes in 1.9.3
  ArrOfInt doublon(elem_.size_array());
  if (supprime_doublons)
    {
      int doublons=0;
      int size = elem_.size();
      for (int i=0; i<size; i++)
        {
          for (int j=i+1; j<size; j++)
            {
              bool same = true;
              if (elem_[i]!=elem_[j])
                same = false;
              else
                for (int dir = 0; dir < dimension; dir++)
                  if (les_positions_sondes_(i, dir) != les_positions_sondes_(j, dir))
                    {
                      same = false;
                      break;
                    }
              if (same)
                {
                  doublon[i] = 1;
                  doublons++;
                  break;
                }
            }
        }
      if (doublons)
        Cerr << "We remove " << doublons << " duplicated points from the probe " << nom_ << finl;
    }

  // each processor has checked if it holds the point
  // the master builds a table (prop) that determines which processor will provide the value to the master
  // the master also builds the ArrsOfInt participant
  // giving for a proc the various elements it owns
  IntVect prop(elem_);
  if (je_suis_maitre())
    {
      ArrOfInt elems2(elem_);
      prop=0; // By default the master is the owner (especially for probes outside the domain)
      ArrOfInt elem_recu, doublon_recu;
      DoubleTab positions_recu;
      int nbproc=Process::nproc();
      for(int p=1; p<nbproc; p++)
        {
          recevoir(doublon_recu,p,2003+p);
          recevoir(elem_recu,p,0,2002+p);
          recevoir(positions_recu,p,0,2001+p);
          for (int el=0; el<nbre_points_tot; el++)
            if (elems2[el]==-1 && elem_recu[el]!=-1)
              {
                elems2[el]=elem_recu[el];
                doublon[el]=doublon_recu[el];
                prop[el]=p;
                for (int j=0; j<dimension; j++)
                  les_positions_sondes_(el,j)=positions_recu(el,j);
              }
        }
      // OK the prop table has been filled;
      // the master sizes participant;
      participant.dimensionner(nbproc);
      for(int p=0; p<nbproc; p++)
        {
          int size=0;
          for (int el=0; el<nbre_points_tot; el++)
            if (prop[el]==p&&!doublon[el]) size++;
          participant[p].resize_array(size);
          participant[p]=-1;
          int pos=0;
          int pt=0;
          for (int el=0; el<nbre_points_tot; el++)
            if (!doublon[el])
              {
                if (prop[el] == p)
                  {
                    participant[p][pos] = pt;
                    pos++;
                  }
                pt++;
              }
          assert((size==0)||(min_array(participant[p])>-1));
        }
    }
  else
    {
      envoyer(doublon, Process::me(),0,2003+Process::me());
      envoyer(elem_,Process::me(),0,2002+Process::me());
      envoyer(les_positions_sondes_,Process::me(),0,2001+Process::me());
    }
  envoyer_broadcast(prop, 0);
  //
  // resize arrays to local size from prop, removing duplicates
  int nbre_points=0;
  nbre_points_tot=0;
  les_positions_=les_positions_sondes_;
  for (int el=0; el<prop.size(); el++)
    {
      if (!doublon[el])
        {
          if (prop(el) == me())
            {
              elem_(nbre_points) = elem_(el);
              for (int dir = 0; dir < dimension; dir++)
                les_positions_(nbre_points, dir) = les_positions_(el, dir);
              nbre_points++;
            }
          for (int dir = 0; dir < dimension; dir++)
            les_positions_sondes_(nbre_points_tot, dir) = les_positions_sondes_(el, dir);
          nbre_points_tot++;
        }
    }
#ifndef NDEBUG
  int test = static_cast<int>(mp_sum(nbre_points));
  //cerr << "Remove " << Process::me() << " " << nbre_points << " " << nbre_points_tot << " = " << test << endl;
  if (je_suis_maitre()) assert(test==nbre_points_tot);
#endif
  elem_.resize(nbre_points);
  les_positions_.resize(nbre_points, dimension);
  if (je_suis_maitre()) // Array only on the master:
    les_positions_sondes_.resize(nbre_points_tot, dimension);
  else
    les_positions_sondes_.resize(0,0);
  // Sizing valeurs_sur_maitre
  if (je_suis_maitre()&&(nproc()>1))
    {
      if (ncomp == -1)
        {
          int nb_comp = get_nb_compo_champ();
          valeurs_sur_maitre.resize(nbre_points_tot,nb_comp);
        }
      else
        valeurs_sur_maitre.resize(nbre_points_tot, 1);
    }
  // Sizing valeurs_locales
  if (ncomp == -1)
    {
      int nb_comp = get_nb_compo_champ();
      valeurs_locales.resize(nbre_points,nb_comp);
    }
  else
    valeurs_locales.resize(nbre_points, 1);
}


/*! @brief Opens the file associated with the probe.
 *
 * (*.son)
 *
 */
void Sonde::ouvrir_fichier()
{
  if (je_suis_maitre() && !le_fichier_.is_open())
    {
      struct stat f;
      int reprise = !stat(nom_fichier_, &f) && mon_post->probleme().reprise_effectuee();
      if (reprise)
        {
          // Restarting a computation, appending to file:
          le_fichier_.ouvrir(nom_fichier_, ios::app);
          le_fichier_.setf(ios::scientific);
          le_fichier_.precision(8);
        }
      else
        {
          // Starting computation or file does not exist:
          le_fichier_.ouvrir(nom_fichier_);
          le_fichier_.setf(ios::scientific);
          le_fichier_.precision(8);
          // Write header:
          SFichier& s = le_fichier_;
          // Write the header of the probe files:
          if (dim == 0 || dim == 1)
            {
              const DoubleTab& p = les_positions_sondes();
              int nbre_points = les_positions_sondes_.dimension(0);
              s << "# " << nom_fichier_ << finl;
              s << "# Temps";
              for (int i = 0; i < nbre_points; i++)
                {
                  s << " x= " << p(i, 0) << " y= " << p(i, 1);
                  if (dimension == 3) fichier() << " z= " << p(i, 2);
                }
              s << finl;
              if (mon_champ)
                {
                  const Noms unites = mon_champ->get_property("unites");
                  s << "# Champ " << nom_champ_lu_ << " [" << unites[ncomp == -1 ? 0 : ncomp] << "]" << finl;
                }
              else
                s << "# Champ " << nom_champ_lu_ << " [??]" << finl;
              s << "# Type " << get_type() << finl;
            }
          // Write the header of the plan files:
          else
            {
              s << "TRUST   Version1  01/09/96" << finl;
              s << nom_du_cas() << finl;
              s << "TRUST" << finl;
              s << "GRILLE";
              const DoubleTab& p = les_positions_sondes();
              // Nouveau on ajoute des informations pour Run_sonde
              s << " " << get_type() << " " << nom_champ_lu_ << " " << nbre_points1 << " " << nbre_points2;
              for (int j = 0; j < dimension; j++) s << " " << p(0, j);
              for (int j = 0; j < dimension; j++) s << " " << p((nbre_points1 - 1) * nbre_points2, j);
              for (int j = 0; j < dimension; j++) s << " " << p(nbre_points2 - 1, j);
              s << finl;
              Nom nom_grille("Grille");
              Nom nom_topologie("Topologie");
              nom_grille += "_";
              nom_grille += nom_;
              nom_topologie += "_";
              nom_topologie += nom_;
              int k, kn;
              int nbre_points = les_positions_sondes_.dimension(0);
              if (dim == 2)
                {
                  double xn = 0., yn = 0., zn = 0., norme = 0.;
                  int p1 = 1, p_nbre_points2 = nbre_points2;
                  while (p(0, 0) == p(p1, 0) && p(0, 1) == p(p1, 1) && p(0, 2) == p(p1, 2))
                    {
                      p1 += 1;
                      assert(p1 < nbre_points);
                    }
                  while (p(0, 0) == p(p_nbre_points2, 0) && p(0, 1) == p(p_nbre_points2, 1) &&
                         p(0, 2) == p(p_nbre_points2, 2))
                    {
                      p_nbre_points2 += 1;
                      assert(p_nbre_points2 < nbre_points);
                    }
                  while (p(p1, 0) == p(p_nbre_points2, 0) && p(p1, 1) == p(p_nbre_points2, 1) &&
                         p(p1, 2) == p(p_nbre_points2, 2))
                    {
                      p_nbre_points2 += 1;
                      assert(p_nbre_points2 < nbre_points);
                    }

                  if (dimension == 3)
                    {
                      xn = (p(p1, 2) - p(0, 2)) * (p(p_nbre_points2, 1) - p(0, 1))
                           - (p(p1, 1) - p(0, 1)) * (p(p_nbre_points2, 2) - p(0, 2));
                      yn = (p(p1, 0) - p(0, 0)) * (p(p_nbre_points2, 2) - p(0, 2))
                           - (p(p1, 2) - p(0, 2)) * (p(p_nbre_points2, 0) - p(0, 0));
                    }
                  else if (dimension == 2)
                    {
                      xn = 0.;
                      yn = 0.;
                    }
                  zn = (p(p1, 1) - p(0, 1)) * (p(p_nbre_points2, 0) - p(0, 0))
                       - (p(p1, 0) - p(0, 0)) * (p(p_nbre_points2, 1) - p(0, 1));
                  norme = std::fabs(xn) + std::fabs(yn) + std::fabs(zn);
                  xn /= norme;
                  yn /= norme;
                  zn /= norme;
                  s << nom_grille << " 3 " << 2 * nbre_points << finl;
                  int i;
                  for (i = 0; i < nbre_points; i++)
                    {
                      s << p(i, 0) << " " << p(i, 1);
                      if (dimension == 3) fichier() << " " << p(i, 2) << finl;
                      else if (dimension == 2) s << " 0." << finl;
                    }
                  for (i = 0; i < nbre_points; i++)
                    {
                      s << p(i, 0) + xn << " " << p(i, 1) + yn;
                      if (dimension == 3) s << " " << p(i, 2) + zn << finl;
                      else if (dimension == 2) s << " " << zn << finl;
                    }
                  s << "TOPOLOGIE" << finl;
                  s << nom_topologie << " " << nom_grille << finl;
                  s << "MAILLE" << finl;
                  s << (nbre_points1 - 1) * (nbre_points2 - 1) << finl;
                  for (int j = 0; j < nbre_points2 - 1; j++)
                    for (i = 0; i < nbre_points1 - 1; i++)
                      {
                        k = j * nbre_points1 + i + 1;
                        kn = k + nbre_points;
                        s << "VOXEL8 " << k << " " << k + 1 << " ";
                        s << k + nbre_points1 << " " << k + nbre_points1 + 1;
                        s << " " << kn << " " << kn + 1 << " " << kn + nbre_points1;
                        s << " " << kn + nbre_points1 + 1 << finl;
                      }
                  s << "FACE" << finl;
                  s << "0" << finl;
                }
              else if (dim == 3)
                {
                  const DoubleTab& pbis = les_positions_sondes();
                  s << nom_grille << " 3 " << nbre_points << finl;
                  int i;
                  for (i = 0; i < nbre_points; i++)
                    {
                      s << pbis(i, 0) << " " << pbis(i, 1);
                      s << " " << pbis(i, 2) << finl;
                    }
                  s << "TOPOLOGIE" << finl;
                  s << nom_topologie << " " << nom_grille << finl;
                  s << "MAILLE" << finl;
                  s << (nbre_points1 - 1) * (nbre_points2 - 1) * (nbre_points3 - 1) << finl;
                  for (int m = 0; m < nbre_points3 - 1; m++)
                    for (int j = 0; j < nbre_points2 - 1; j++)
                      for (i = 0; i < nbre_points1 - 1; i++)
                        {
                          k = m * nbre_points2 * nbre_points1 + j * nbre_points1 + i + 1;
                          kn = k + nbre_points2 * nbre_points1;
                          s << "VOXEL8 " << k << " " << k + 1 << " ";
                          s << k + nbre_points1 << " " << k + nbre_points1 + 1;
                          s << " " << kn << " " << kn + 1 << " " << kn + nbre_points1;
                          s << " " << kn + nbre_points1 + 1 << finl;
                        }
                  s << "FACE" << finl;
                  s << "0" << finl;
                }
            }
        }
    }
}

void Sonde::update_source(double un_temps)
{
  // Cache mechanism for the source field behind the post-processed field (mon_champ)
  // Implemented at the Sondes level
  OBS_PTR(Champ_base) ma_source = mon_post->les_sondes().get_from_cache(mon_champ, nom_champ_lu_);
  ma_source->mettre_a_jour(un_temps);
}

/*! @brief Updates the probe in time and performs the post-processing.
 *
 * @param un_temps the update time
 * @param tinit the initial time of the probe
 */
void Sonde::mettre_a_jour(double un_temps, double tinit)
{
  double temps_courant = get_temps_champ();
  double dt=mon_post->probleme().schema_temps().pas_de_temps();
  double nb;

  // The *(1+Objet_U::precision_geom) factor avoids rounding errors depending on the machine
  // 21/01/25: replaced Objet_U::precision_geom with 1e-15 (from M3D patch, bttrust #242895)
  // because precision_geom is a user-settable parameter via the data file
  if (periode<=dt)
    nb=nb_bip+1.;
  else
    modf(temps_courant*(1+1e-15)/periode, &nb);

  // Time to write the probes
  if (nb>nb_bip)
    {
      update_source(un_temps);

      // Si le maillage est deformable il faut reconstruire les sondes
      if (mon_post->les_sondes().get_update_positions())
        {
          if (mon_post->probleme().domaine().deformable())
            initialiser();
        }
      nb_bip=nb;
      postraiter();
    }
}


void Sonde::fill_local_values()
{
  // Cache mechanism for the source field behind the post-processed field (mon_champ)
  // Implemented at the Sondes level
  OBS_PTR(Champ_base) ma_source = mon_post->les_sondes().get_from_cache(mon_champ, nom_champ_lu_);
  if (chsom)
    {
      Champ_base& ma_source_mod =ref_cast_non_const(Champ_base,ma_source.valeur());
      if (ncomp == -1)
        ma_source_mod.valeur_aux_elems_smooth(les_positions(),elem_, valeurs_locales);
      else
        ma_source_mod.valeur_aux_elems_compo_smooth(les_positions(),elem_,valeurs_locales, ncomp);
    }
  else
    {
      if (ncomp == -1)
        ma_source->valeur_aux_elems(les_positions(),elem_, valeurs_locales);
      else
        ma_source->valeur_aux_elems_compo(les_positions(),elem_,valeurs_locales, ncomp);
    }
}

/*! @brief Performs post-processing.
 *
 * Computes field values at the requested positions and writes them to the associated file.
 *
 */
void Sonde::postraiter()
{
  fill_local_values();

  if (gravcl)
    mettre_a_jour_bords();
  const int N = valeurs_locales.line_size();

  // the master reconstructs the values array from the various contributors
  if(je_suis_maitre())
    {
      ouvrir_fichier();
      int nbproc = Process::nproc();
      DoubleTab valeurs_pe;
      DoubleTab& valeurs=(nbproc==1?valeurs_locales:valeurs_sur_maitre);
      if (nbproc==1) { /* Do nothing */} //valeurs=valeurs_locales;
      else
        {
          int nb_val = valeurs_locales.dimension(0);
          for (int i = 0; i < nb_val; i++)
            for (int n = 0; n < N; n++)
              valeurs(participant[0][i], n) = valeurs_locales(i, n);
        }
      for(int p=1; p<nbproc; p++)
        {
          // le message n'est envoye que si le proc participe
          if (participant[p].size_array()!=0)
            {
              recevoir(valeurs_pe,p,0,2002+p);
              int nb_val1=valeurs_pe.dimension(0);
              assert(N == valeurs_pe.line_size());
              for(int i=0; i<nb_val1; i++)
                for(int n = 0; n < N; n++)
                  {
                    valeurs(participant[p][i], n)=valeurs_pe(i, n);
                    //  val_max = std::max(std::fabs(valeurs(i,k)),std::fabs(valeurs_pe(i,k)));
                    //                     if(val_max==(std::fabs(valeurs_pe(i,k))))
                    //                       valeurs(i,k)=valeurs_pe(i,k);
                  }
            }
        }

      double temps_courant = mon_post->probleme().schema_temps().temps_courant();

      if (dim==0 || dim==1)
        {
          fichier() << temps_courant;
          for(int i=0; i<valeurs.dimension(0); i++)
            for(int k=0; k<N; k++)
              fichier() << " " << valeurs(i,k);
          fichier() << finl;
        }
      // For plan-type probes, output in lml format:
      // num_sommet comp1 [comp2] [comp3]
      // and in the third direction:
      else if (dim==2 || dim==3)
        {
          Nom nom_post;
          int nbre_points = valeurs.dimension(0);
          const Noms noms_comp = mon_champ->get_property("composantes");
          int nb_comp = noms_comp.size();
          const Noms unites = mon_champ->get_property("unites");
          if (ncomp==-1)
            {
              const Noms noms_champ = mon_champ->get_property("nom");
              nom_post = noms_champ[0];
            }
          else
            nom_post = noms_comp[ncomp];

          Nom nom_topologie("Topologie");
          nom_topologie += "_";
          nom_topologie += nom_;

          fichier() << "TEMPS " << temps_courant << "\n";
          fichier() << "CHAMPPOINT " << nom_post << " " << nom_topologie
                    << " " << temps_courant << "\n";
          fichier() << nom_post << " " << nb_comp << " " << unites[0] << "\n";

          int nbp=nbre_points;
          if (dim==2) nbp*=2;
          if (nb_comp>1)
            fichier() << "type1 " << nbp << "\n";
          else fichier() << "type0 " << nbp << "\n";
          int i;
          for(i=0; i<nbre_points; i++)
            {
              fichier() << i+1;
              for(int j=0; j<N; j++)
                fichier() << " " << valeurs(i,j);
              // Pour ne pas flusher :
              fichier() << "\n";
            }
          // For 2D, add an extra direction
          if (dim==2)
            {
              for(i=0; i<nbre_points; i++)
                {
                  fichier() << nbre_points+i+1;
                  for(int j=0; j<N; j++)
                    fichier() << " " << valeurs(i,j);
                  // Pour ne pas flusher :
                  fichier() << "\n";
                }
            }
        }
      fichier().flush();
    }
  else
    {
      // the processor sends a message only if it contributes
      if (valeurs_locales.dimension(0)!=0)
        envoyer(valeurs_locales,Process::me(),0,2002+Process::me());
    }
}

void Sonde::init_bords()
{
  if (!gravcl)
    return;

  const int nbf = faces_bords_.size_array();
  rang_cl_.resize(nbf);

  if (sub_type(Champ_Generique_refChamp,mon_champ.valeur()))
    {
      Probleme_base& Pb = mon_post->probleme();
      OBS_PTR(Champ_base) chref = Pb.get_champ(nom_champ_lu_);
      const Champ_Inc_base& ch_inc = ref_cast(Champ_Inc_base,chref.valeur());
      const Domaine_Cl_dis_base& zcl = ch_inc.domaine_Cl_dis();

      for (int iface=0; iface<nbf; iface++)
        {
          int face = faces_bords_[iface];
          for (int icl=0; icl<zcl.nb_cond_lim(); icl++)
            {
              const Cond_lim& cl = zcl.les_conditions_limites(icl);
              const Front_VF& le_bord = ref_cast(Front_VF,cl->frontiere_dis());
              const int ndeb = le_bord.num_premiere_face();
              const int nfin = ndeb + le_bord.nb_faces();
              for (int i=ndeb; i<nfin; i++)
                {
                  if (i==face)
                    rang_cl_[iface] = icl;
                }
            }
        }
    }
}

void Sonde::ajouter_bords(const DoubleTab& coords_bords)
{
  double eps = mon_champ->get_ref_domain().epsilon();
  const Domaine_VF& domaineVF = ref_cast(Domaine_VF,mon_champ->get_ref_domaine_dis_base());

  const DoubleTab& xv = domaineVF.xv();
  const IntTab& e_f = domaineVF.elem_faces(), &f_e = domaineVF.face_voisins();
  int nbre_points = les_positions_sondes_.dimension(0), face;



  IntTab tmp(2);
  tmp[0] = 0;
  tmp[1] = nbre_points-1;

  for (int i = 0; i<coords_bords.dimension(0); i++)
    {
      int e = tmp[i];
      int elem = elem_[e];
      if (elem !=-1)
        {
          for (int j = 0; j < e_f.dimension(1) && (face = e_f(elem, j)) >= 0; j++)
            {
              if(f_e(face, 1) == -1 || f_e(face, 0) == -1)
                {
                  double dist=0.;
                  for (int idim=0; idim<dimension; idim++)
                    dist+=((xv(face,idim)-coords_bords(i,idim))*(xv(face,idim)-coords_bords(i,idim)));

                  if (sqrt(dist) < eps)
                    {
                      faces_bords_.append_array(face);
                      for (int idim=0; idim<dimension; idim++)
                        les_positions_sondes_(e,idim) = coords_bords(i,idim);
                    }
                }
            }
        }
    }
// [ABN] skiping the assert below. On funny shaped domains (e.g. U-shape) a probe might
// re-enter the domain, and hence hit more than two faces overall:
// assert(mp_sum(faces_bords_.size_array())<=2);
}

void Sonde::mettre_a_jour_bords()
{
  if (sub_type(Champ_Generique_refChamp,mon_champ.valeur()))
    {
      Probleme_base& Pb = mon_post->probleme();
      OBS_PTR(Champ_base) chref = Pb.get_champ(nom_champ_lu_);
      const Champ_Inc_base& ch_inc = ref_cast(Champ_Inc_base,chref.valeur());
      const Domaine_Cl_dis_base& zcl = ch_inc.domaine_Cl_dis();
      const DoubleTab& vals_ch = ch_inc.valeurs();
      const Domaine_VF& domaineVF = ref_cast(Domaine_VF,mon_champ->get_ref_domaine_dis_base());
      const IntTab& face_voisins = domaineVF.face_voisins();

      int nbval = valeurs_locales.dimension(0);

      for (int iface=0; iface<faces_bords_.size_array(); iface++)
        {
          int face = faces_bords_[iface];
          int icl = rang_cl_[iface];
          const Cond_lim& cl = zcl.les_conditions_limites(icl);
          const Front_VF& le_bord = ref_cast(Front_VF,cl->frontiere_dis());
          const int ndeb = le_bord.num_premiere_face();
          DoubleVect valcl;
          cl->champ_front().valeurs_face(face-ndeb,valcl);

          if (valcl.size() == 0)
            {
              int elem0 = face_voisins(face,0);
              int elem1 = face_voisins(face,1);

              if (elem0 != -1)
                valeurs_locales(nbval-1) = vals_ch[elem0];
              else
                valeurs_locales(0) = vals_ch[elem1];
            }
          else
            {
              if (face_voisins(face,0) != -1)
                valeurs_locales(nbval-1) = valcl[0];
              else
                valeurs_locales(0) = valcl[0];
            }
        }
    }
}
