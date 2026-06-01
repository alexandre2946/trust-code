/****************************************************************************
* Copyright (c) 2025, CEA
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

#include <SolveurSys_base.h>
#include <EcrFicCollecte.h>
#include <Matrice_Base.h>
#include <MD_Vector_tools.h>
#include <SChaine.h>
#include <Motcle.h>

Implemente_base_sans_constructeur(SolveurSys_base,"SolveurSys_base",Objet_U);
// XD solveur_sys_base class_generic solveur_sys_base INHERITS_BRACE Basic class to solve the linear system.

SolveurSys_base::SolveurSys_base()
{
}

Sortie& SolveurSys_base::printOn(Sortie& s ) const
{
  return s  << que_suis_je();
}

Entree& SolveurSys_base::readOn(Entree& is )
{
  return is;
}

int SolveurSys_base::resoudre_systeme(const Matrice_Base& A,
                                      const DoubleVect& b,
                                      DoubleVect& x,
                                      int niter_max)
{
  save_matrice_secmem_conditionnel(A, b, x);
  return resoudre_systeme(A,b,x);
}


void SolveurSys_base::save_matrice_secmem_conditionnel(const Matrice_Base& la_matrice, const DoubleVect& secmem, const DoubleVect& solution, int binaire2 )
{
  int binaire=binaire2;
  if (save_matrice_==1)
    {
      // simulate writing a Matrice rather than a Matrice_Base
      // to facilitate re-reading
      {
        EcrFicCollecte sortie;
        sortie.set_bin(binaire);
        sortie.ouvrir("Matrice.sa");
        sortie<<la_matrice.que_suis_je()<<finl;
        sortie<<la_matrice;
      }
      // save the RHS vector (Secmem) with virtual spaces !!!!!
      {
        EcrFicCollecte sortie;
        sortie.set_bin(binaire);
        sortie.ouvrir("Secmem.sa");
        MD_Vector_tools::dump_vector_with_md(secmem, sortie);
      }
      // save the solution with virtual spaces !!!!!
      {
        EcrFicCollecte sortie;
        sortie.set_bin(binaire);
        sortie.ouvrir("Solution.sa");
        MD_Vector_tools::dump_vector_with_md(solution, sortie);
      }
      Cout <<"Saving of the matrix, secmem and solution ended."<<finl;
    }
}


// Read solver parameters
// Ex: solveur type { ... }
// chaine_lue_ = "type { ... }"
void SolveurSys_base::lecture(Entree& is)
{
  // Read the keyword string
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;
  Nom nomlu;
  is >> motlu;
  SChaine prov;
  prov<<motlu;
  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "Error while reading parameters of Petsc: " << finl;
      Cerr << "We expected " << accolade_ouverte << " instead of " << motlu << finl;
      Process::exit();
    }
  prov<<" { ";
  int nb_acco=1;
  while (nb_acco!=0)
    {
      is >> nomlu;
      if ((Motcle)nomlu=="READ_MATRIX") set_read_matrix(true);
      prov<<nomlu<<" ";
      if (nomlu==accolade_ouverte)
        nb_acco++;
      else if (nomlu==accolade_fermee)
        nb_acco--;
    }
  chaine_lue_=prov.get_str();
}

