/****************************************************************************
* Copyright (c) 2024, CEA
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
#include <Ensemble_Lagrange_base.h>
#include <Sous_Domaine.h>
#include <Domaine.h>
#include <communications.h>
Implemente_base_sans_constructeur(Ensemble_Lagrange_base,"Ensemble_Lagrange_base",Objet_U);

Ensemble_Lagrange_base::Ensemble_Lagrange_base()
{
  nb_marqs_sz.resize(0);
  nb_marqs_par_dir.resize(0,0);
}

Sortie& Ensemble_Lagrange_base::printOn(Sortie& os) const
{
  return os;
}

Entree& Ensemble_Lagrange_base::readOn(Entree& is)
{
  return is;
}

void Ensemble_Lagrange_base::associer_domaine(const Domaine& un_domaine)
{
  mon_dom_=un_domaine;
}

void Ensemble_Lagrange_base::remplir_sommets_tmp(DoubleTab& soms_tmp)
{
  //First case: the vertices were read from a file
  //Copy sommets_lu_ into a temporary array sommets_tmp

  if (sommets_lu_.dimension(0)!=0)
    {
      const DoubleTab& som_lu = sommets_lu() ;
      int dim0=som_lu.dimension(0);
      int dim1=som_lu.dimension(1);
      soms_tmp.resize(dim0,dim1);
      soms_tmp = som_lu;
    }
  else if (nb_marqs_par_sz().size()!=0)
    //Second case: the vertices are created in sub-domains
    //Fill the temporary array
    generer_marqueurs_sz(soms_tmp);
  else
    soms_tmp.resize(0,0);
}

//Generation of particle coordinates over one or several sub-domains
//For each sub-domain:
// - Step 1: determine the min and max for x, y and z, taking all processors into account
// - Step 2: the master processor computes the particle coordinates
//                -either with a uniform distribution
//                -or with a random distribution
// - Step 3: the master processor sends the computed coordinate values to the other processors

void Ensemble_Lagrange_base::generer_marqueurs_sz(DoubleTab& soms_tmp)
{
  const Domaine& madomaine = mon_dom_.valeur();
  const Domaine& dom = madomaine;
  int nb_sz,nb_marq_sz, old_size, dim;
  soms_tmp.resize(0);
  dim=Objet_U::dimension;
  nb_sz=nb_marqs_sz.size();

  for (int i=0; i<nb_sz; i++)
    {
      const Nom& nom = nom_sz[i];
      const Sous_Domaine& sz = dom.ss_domaine(nom);
      int nb_elem = sz.nb_elem_tot();

      int nb_som = madomaine.nb_som_elem();
      int le_poly,poly,le_som;

      double xmoy,ymoy,zmoy,deltax,deltay,deltaz;
      double xmin,xmax,ymin,ymax,zmin,zmax;
      xmin = ymin = zmin = 1.e8;
      xmax = ymax = zmax= -1.e8;

      if (dim==2)
        {
          zmin = zmax = 0.;
        }

      //Step 1: determine the min and max taking all processors into account
      if (nb_elem!=0)
        {
          poly = sz(0);
          le_som = 0;

          xmin=dom.coord(madomaine.sommet_elem(poly,le_som),0);
          ymin=dom.coord(madomaine.sommet_elem(poly,le_som),1);
          if (dim==3)
            zmin=dom.coord(madomaine.sommet_elem(poly,le_som),2);
          xmax=xmin;
          ymax=ymin;
          if (dim==3)
            zmax=zmin;

          for (le_poly=0; le_poly<nb_elem; le_poly++)
            for(le_som=0; le_som<nb_som; le_som++)
              {
                poly = sz(le_poly);
                xmin=std::min(xmin,dom.coord(madomaine.sommet_elem(poly,le_som),0));
                ymin=std::min(ymin,dom.coord(madomaine.sommet_elem(poly,le_som),1));
                if (dim==3)
                  zmin=std::min(zmin,dom.coord(madomaine.sommet_elem(poly,le_som),2));
                xmax=std::max(xmax,dom.coord(madomaine.sommet_elem(poly,le_som),0));
                ymax=std::max(ymax,dom.coord(madomaine.sommet_elem(poly,le_som),1));
                if (dim==3)
                  zmax=std::max(zmax,dom.coord(madomaine.sommet_elem(poly,le_som),2));
              }
        } //End if (nb_elem!=0)


      xmin = Process::mp_min(xmin);
      xmax = Process::mp_max(xmax);
      ymin = Process::mp_min(ymin);
      ymax = Process::mp_max(ymax);
      if (dim==3)
        {
          zmin = Process::mp_min(zmin);
          zmax = Process::mp_max(zmax);
        }

      //Step 2: the master processor will compute the particle coordinates
      if (je_suis_maitre())
        {

          xmoy = (xmax+xmin)/2.;
          ymoy = (ymax+ymin)/2.;
          zmoy = (zmax+zmin)/2.;

          deltax = (xmax-xmin)/2.;
          deltay = (ymax-ymin)/2.;
          deltaz = (zmax-zmin)/2.;

          nb_marq_sz = nb_marqs_sz(i);
          old_size = soms_tmp.dimension(0);
          soms_tmp.resize(old_size+nb_marq_sz,dim);
          int marq;

          //Case of a random distribution
          if (nb_marqs_par_dir(i,0)==0)
            {
              //double a,b,c;
              srand48(1);
              for (marq=0; marq<nb_marq_sz; marq++)
                {
                  soms_tmp(old_size+marq,0) = xmoy + deltax*(-1.+2.*drand48());
                  soms_tmp(old_size+marq,1) = ymoy + deltay*(-1.+2.*drand48());
                  if (dim==3)
                    soms_tmp(old_size+marq,2) = zmoy + deltaz*(-1.+2.*drand48());
                }
            }
          //Case of a uniform distribution
          else
            {
              double dx, dy, dz;
              int nx_max = nb_marqs_par_dir(i,0);
              int ny_max = nb_marqs_par_dir(i,1);
              dx = 2.*deltax/double(nx_max+1);
              dy = 2.*deltay/double(ny_max+1);

              if (dim==2)
                {
                  for (int ny=1; ny<ny_max+1; ny++)
                    for (int nx=1; nx<nx_max+1; nx++)
                      {
                        marq = nx + (ny-1)*nx_max-1;
                        soms_tmp(old_size+marq,0) = xmin + nx*dx;
                        soms_tmp(old_size+marq,1) = ymin + ny*dy;
                      }
                }
              else if (dim==3)
                {
                  int nz_max = nb_marqs_par_dir(i,2);
                  dz = 2.*deltaz/double(nz_max+1);
                  for (int nz=1; nz<nz_max+1; nz++)
                    for (int ny=1; ny<ny_max+1; ny++)
                      for (int nx=1; nx<nx_max+1; nx++)
                        {
                          marq = nx + (ny-1)*nx_max + (nz-1)*nx_max*ny_max-1;
                          soms_tmp(old_size+marq,0) = xmin + nx*dx;
                          soms_tmp(old_size+marq,1) = ymin + ny*dy;
                          soms_tmp(old_size+marq,2) = zmin + nz*dz;
                        }
                }
            }
        }//End of if je_suis_maitre()

      //Step 3: the master processor sends the computed coordinate values (soms_tmp)
      //             to the other processors
      if (je_suis_maitre())
        for(int p=1; p<Process::nproc(); p++)
          envoyer(soms_tmp,0, p, 0);
      for(int p=1; p<Process::nproc(); p++)
        recevoir(soms_tmp,0,p,0);
    }
}
