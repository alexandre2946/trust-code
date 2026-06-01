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

#include <Perte_Charge_Reguliere.h>
#include <Motcle.h>

// XD perte_charge_reguliere source_base perte_charge_reguliere NO_BRACE Source term modelling the presence of a bundle
// XD_CONT of tubes in a flow.
// XD attr spec spec_pdcr_base spec REQ Description of longitudinale or transversale type.
// XD attr zone_name chaine name_of_zone REQ Name of the sub-area occupied by the tube bundle. A Sous_Zone (Sub-area)
// XD_CONT type object called zone_name should have been previously created.

// XD spec_pdcr_base objet_lecture spec_pdcr_base NO_BRACE Class to read the source term modelling the presence of a
// XD_CONT bundle of tubes in a flow. Cf=A Re-B.

// XD longitudinale spec_pdcr_base longitudinale NO_BRACE Class to define the pressure loss in the direction of the tube
// XD_CONT bundle.
// XD attr dir chaine(into=["x","y","z"]) dir REQ Direction.
// XD attr dd floattant dd REQ Tube bundle hydraulic diameter value. This value is expressed in m.
// XD attr ch_a chaine(into=["a","cf"]) ch_a REQ Keyword to be used to set law coefficient values for the coefficient of
// XD_CONT regular pressure losses.
// XD attr a floattant a REQ Value of a law coefficient for regular pressure losses.
// XD attr ch_b chaine(into=["b"]) ch_b OPT Keyword to be used to set law coefficient values for regular pressure
// XD_CONT losses.
// XD attr b floattant b OPT Value of a law coefficient for regular pressure losses.

// XD transversale spec_pdcr_base transversale NO_BRACE Class to define the pressure loss in the direction perpendicular
// XD_CONT to the tube bundle.
// XD attr dir chaine(into=["x","y","z"]) dir REQ Direction.
// XD attr dd floattant dd REQ Value of the tube bundle step.
// XD attr chaine_d chaine(into=["d"]) chaine_d REQ Keyword to be used to set the value of the tube external diameter.
// XD attr d floattant d REQ Value of the tube external diameter.
// XD attr ch_a chaine(into=["a","cf"]) ch_a REQ Keyword to be used to set law coefficient values for the coefficient of
// XD_CONT regular pressure losses.
// XD attr a floattant a REQ Value of a law coefficient for regular pressure losses.
// XD attr ch_b chaine(into=["b"]) ch_b OPT Keyword to be used to set law coefficient values for regular pressure
// XD_CONT losses.
// XD attr b floattant b OPT Value of a law coefficient for regular pressure losses.


//Add a constructor for initialization
Perte_Charge_Reguliere::Perte_Charge_Reguliere() : Perte_Charge( )
{

  D_ = - 1e3;
  d_ = - 1e3;
  a_ = - 1e3;
  b_ = - 1e3;
  Cf_ =  -1e3;

  couronne_tube = -1e3; ;
  Cf_utilisateur = -1e3 ;
  xo = -1e3;
  yo = -1e3;
  zo = -1e3;
  dir[ 0 ] = -1e3 ;
  dir[ 1 ] = -1e3 ;
  dir[ 2 ] = -1e3 ;
  h0 = -1e3;
  h = -1e3 ;

}

/*! @brief Reads the specifications of a regular pressure drop from an input stream.
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the modified input stream
 * @throws unexpected keyword, expected "X","Y" or "Z"
 * @throws unexpected keyword, expected "d"
 * @throws unexpected keyword, expected "longitudinale" or "transversale"
 */
Entree& Perte_Charge_Reguliere::lire_donnees(Entree& is)
{
  int dimension = Objet_U::dimension ;

  for (int i =0; i<dimension; i++ ) dir[i] = 0 ;

  couronne_tube = 0 ;
  Cf_utilisateur = 0 ;

  Motcle motlu;
  Motcles les_motcles(6);
  les_motcles[0] = "X";
  les_motcles[1] = "Y";
  les_motcles[2] = "Z";
  les_motcles[3] = "A";
  les_motcles[4] = "B";
  les_motcles[5] = "CF";

  is >> motlu;
  if (motlu == "longitudinale")
    {
      is >> motlu;
      int rang = les_motcles.search(motlu);
      if (rang == -1)
        {
          Cerr << "Error reading Perte_Charge_Reguliere data" << finl;
          Cerr << "Expected one of the keywords" << les_motcles << finl;
          Cerr << "instead of " << motlu << finl;
          Process::exit();
        }
      else
        {
          direction_perte_charge_ = rang;
          is >> D_;
          Cerr << "D_ " << D_ << finl;
        }
      d_ = D_;
    }
  else if (motlu == "transversale")
    {
      is >> motlu;
      int rang = les_motcles.search(motlu);
      if (rang == -1)
        {
          Cerr << "Error reading Perte_Charge_Reguliere data" << finl;
          Cerr << "Expected one of the keywords" << les_motcles << finl;
          Cerr << "instead of " << motlu << finl;
          Process::exit();
        }
      else
        {
          direction_perte_charge_ = rang;
          is >> D_;
        }
      is >> motlu;
      if (motlu != "d")
        {
          Cerr << "Error reading Perte_Charge_Reguliere data" << finl;
          Cerr << "Expected keyword d instead of " << motlu << finl;
          Process::exit();
        }
      else
        {
          is >> d_;
        }
    }
  else if (motlu == "couronne")
    {
      couronne_tube = 1 ;

      if(dimension!=2)
        {
          // A ring (couronne) is in 2D only!
          Cerr << "A ring (couronne) is in 2D only: in 3D use Tube "<< finl;
          Process::exit();
        }
      is >> motlu;
      if(motlu!=Motcle("Origine"))
        {
          Cerr << "Expected keyword ORIGINE" << finl;
          Process::exit();
        }

      is >> xo >> yo ;
      is >> motlu;
      if (motlu != "d")
        {
          Cerr << "Error reading Perte_Charge_Reguliere data" << finl;
          Cerr << "Expected keyword d instead of " << motlu << finl;
          Process::exit();
        }
      else
        {
          is >> d_;
          D_ = d_;
        }
      direction_perte_charge_ = 4;
      dir[0] = 1 ;
      dir[1] = 1 ;
    }
  else if (motlu == "tube")
    {
      couronne_tube = 1 ;
      if(dimension!=3)
        {
          // A tube is in 3D only!
          Cerr << "A tube is in 3D only: in 2D use couronne "<< finl;
          Process::exit();
        }
      is >> motlu;
      if(motlu!=Motcle("Origine"))
        {
          Cerr << "Expected keyword ORIGINE" << finl;
          Process::exit();
        }
      is >> xo >> yo >> zo;
      is >> motlu;
      if(motlu!=Motcle("DIR"))
        {
          Cerr << "Expected the tube direction keyword: DIR" << finl;
          Process::exit();
        }
      int idir;
      is >> idir;
      switch(idir)
        {
        case 0:
          dir[0]=0;
          dir[1]=dir[2]=1;
          h0=xo;
          break;
        case 1:
          dir[1]=0;
          dir[0]=dir[2]=1;
          h0=yo;
          break;
        case 2:
          dir[2]=0;
          dir[1]=dir[0]=1;
          h0=zo;
          break;
        default:
          h0=-1;
          Cerr << "DIR is 0 for tube // OX ; 1 for a tube // OY and 2 for a tube // OZ " << finl;
          Process::exit();
        }
      is >> motlu;
      if (motlu != "d")
        {
          Cerr << "Error reading Perte_Charge_Reguliere data" << finl;
          Cerr << "Expected keyword d instead of " << motlu << finl;
          Process::exit();
        }
      else
        {
          is >> d_;
          D_ = d_;
        }
      direction_perte_charge_ = 4;
    }
  else
    {
      Cerr << "Error reading Perte_Charge_Reguliere data" << finl;
      Cerr << "Expected keyword longitudinale or transversale or couronne (2D) or tube (3D)" << finl;
      Cerr << "instead of " << motlu << finl;
      Process::exit();
    }

  // Mandatory reading of a_ and b_:
  int err=0;
  is >> motlu;
  if (motlu == "A")
    is >> a_;
  else
    err=1;
  if (err)
    {
      if (motlu == "CF")
        {
          is >> Cf_;
          Cf_utilisateur = 1;
          Cerr << " Cf_ " << Cf_ << finl;
        }
    }
  else
    {
      is >> motlu;
      if (motlu == "B")
        is >> b_;
      else
        err=1;

    }
  if (err)
    {
      if (!Cf_utilisateur)
        {
          Cerr << "Error reading Perte_Charge_Reguliere data" << finl;
          Cerr << "You must now specify the values of the pressure drop law:" << finl;
          Cerr << "Example: Perte_Charge_Reguliere longitudinale X val A val B val" << finl;
          Cerr << "With A and B, for a pressure drop law A.Re**(-B)" << finl;
          Cerr << "For Blasius law: A=0.316, B=0.25" << finl;
          Process::exit();
        }
    }
  else
    Cerr << a_ <<  " " << b_ << finl;

  return is;
}
