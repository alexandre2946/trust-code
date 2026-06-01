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

#include <Transport_Interfaces_base.h>
#include <MD_Vector_composite.h>
#include <Matrice_Morse_Sym.h>
#include <Navier_Stokes_std.h>
#include <Frontiere_dis_base.h>
#include <Schema_Temps_base.h>
#include <MD_Vector_tools.h>
#include <Assembleur_base.h>
#include <TRUSTTab_parts.h>
#include <communications.h>
#include <Champ_Uniforme.h>
#include <MD_Vector_std.h>
#include <solv_iteratif.h>
#include <Probleme_base.h>
#include <Discret_Thyd.h>
#include <Fluide_base.h>
#include <Domaine_VF.h>
#include <TRUSTTrav.h>
#include <SFichier.h>
#include <Domaine.h>
#include <Param.h>
#include <Avanc.h>
#include <Debog.h>

#include <TRUST_2_PDI.h>

Implemente_instanciable_sans_constructeur(Navier_Stokes_std,"Navier_Stokes_standard",Equation_base);
// XD navier_stokes_standard eqn_base navier_stokes_standard INHERITS_BRACE Navier-Stokes equations.

Navier_Stokes_std::Navier_Stokes_std():methode_calcul_pression_initiale_(0),div_u_nul_et_non_dsurdt_divu_(0),postraitement_gradient_P_(0)
{
  projection_initiale = 1;
  dt_projection=1.e12;
  max_div_U=1.e12;
  seuil_projection = 0.;
  seuil_uzawa = 1.e-12;
  seuil_divU = 100.;
  cumulative_=0;
  raison_seuil_divU=-100;
}

Sortie& Navier_Stokes_std::printOn(Sortie& is) const
{
  return Equation_base::printOn(is);
}

/*! @brief Calls Equation_base::readOn(Entree& is). On exit, verifies that the following have been read:
 *
 *         - the diffusive term,
 *         - the convective term,
 *         - the pressure solver
 *
 * @param (Entree& is) input stream
 * @return (Entree&) the modified input stream
 * @throws diffusive term not specified in data file; specify a negligeable type if it should be neglected
 * @throws convective term not specified in data file; specify a negligeable type if it should be neglected
 * @throws pressure solver not defined in data file
 */
Entree& Navier_Stokes_std::readOn(Entree& is)
{
  Equation_base::readOn(is);
  if (est_egal(seuil_projection,0.) && (sub_type(solv_iteratif,solveur_pression_.valeur())))
    {
      solv_iteratif& solv_iter = ref_cast(solv_iteratif,solveur_pression_.valeur());
      seuil_projection = solv_iter.get_seuil();
    }

  terme_convectif.set_fichier("Convection_qdm");
  terme_convectif.set_description("Momentum flow rate=Integral(rho*u*u*ndS) [N] if SI units used");
  terme_diffusif.set_fichier("Contrainte_visqueuse");
  terme_diffusif.set_description("Friction drag exerted by the fluid=Integral(-mu*(grad(u) +grad(u)^T)*ndS) [N] if SI units used");
  divergence.set_fichier("Debit");
  divergence.set_description((Nom)"Volumetric flow rate=Integral(u*ndS) [m"+(Nom)(dimension+bidim_axi)+".s-1] if SI units used");
  gradient.set_fichier("Force_pression");
  gradient.set_description("Pressure drag exerted by the fluid=Integral(P*ndS) [N] if SI units used");

  return is;
}

void Navier_Stokes_std::set_param(Param& param) const
{
  Equation_base::set_param(param);
  param.ajouter_non_std("diffusion",(this));
  param.ajouter_non_std("convection",(this));
  param.ajouter_condition("is_read_diffusion","The diffusion operator must be read, select negligeable type if you want to neglect it.");
  param.ajouter_condition("is_read_convection","The convection operator must be read, select negligeable type if you want to neglect it.");
  param.ajouter_non_std("solveur_pression",(this),Param::REQUIRED); // XD attr solveur_pression solveur_sys_base solveur_pression OPT Linear pressure system resolution method.
  param.ajouter_non_std("dt_projection",(this));  // XD attr dt_projection deuxmots dt_projection OPT nb value : This
  // XD_CONT keyword checks every nb time-steps the equality of velocity divergence to zero. value is the criteria
  // XD_CONT convergency for the solver used.
  param.ajouter_non_std("Traitement_particulier",(this)); // XD attr traitement_particulier traitement_particulier traitement_particulier OPT Keyword to post-process particular values.
  param.ajouter_non_std("Erreur_max_DivU",(this));
  param.ajouter("uzawa",&seuil_uzawa);
  //param.ajouter_non_std("vitesse_transportante",(this));
  param.ajouter_non_std("seuil_divU",(this));  // XD attr seuil_divU floatfloat seuil_divU OPT value factor : this
  // XD_CONT keyword is intended to minimise the number of iterations during the pressure system resolution. The
  // XD_CONT convergence criteria during this step (\'seuil\' in solveur_pression) is dynamically adapted according to
  // XD_CONT the mass conservation. At tn , the linear system Ax=B is considered as solved if the residual
  // XD_CONT ||Ax-B||<seuil(tn). For tn+1, the threshold value seuil(tn+1) will be evualated as: NL2 If (
  // XD_CONT |max(DivU)*dt|<value ) NL2 Seuil(tn+1)= Seuil(tn)*factor NL2 Else NL2 Seuil(tn+1)= Seuil(tn)*factor NL2
  // XD_CONT Endif NL2 The first parameter (value) is the mass evolution the user is ready to accept per timestep, and
  // XD_CONT the second one (factor) is the factor of evolution for \'seuil\' (for example 1.1, so 10% per timestep).
  // XD_CONT Investigations has to be lead to know more about the effects of these two last parameters on the behaviour
  // XD_CONT of the simulations.
  param.ajouter_non_std("solveur_bar",(this)); // XD attr solveur_bar solveur_sys_base solveur_bar OPT This keyword is
  // XD_CONT used to define when filtering operation is called (typically for EF convective scheme, standard diffusion
  // XD_CONT operator and Source_Qdm_lambdaup ). A file (solveur.bar) is then created and used for inversion procedure.
  // XD_CONT Syntax is the same then for pressure solver (GCP is required for multi-processor calculations and, in a
  // XD_CONT general way, for big meshes).
  param.ajouter("projection_initiale",&projection_initiale); // XD attr projection_initiale entier projection_initiale OPT Keyword to suppress, if boolean equals 0, the initial projection which checks DivU=0. By default, boolean equals 1.
  param.ajouter_flag("postraiter_gradient_pression_sans_masse", &postraiter_gradient_pression_sans_masse_); // XD attr postraiter_gradient_pression_sans_masse rien postraiter_gradient_pression_sans_masse OPT Avoid mass matrix multiplication for the gradient postprocessing
  param.ajouter_non_std("methode_calcul_pression_initiale",(this));  // XD attr methode_calcul_pression_initiale chaine(into=["avec_les_cl","avec_sources","avec_sources_et_operateurs","sans_rien"]) methode_calcul_pression_initiale OPT Keyword to select an option for the pressure calculation before the fist time step. Options are : avec_les_cl (default option lapP=0 is solved with Neuman boundary conditions on pressure if any), avec_sources (lapP=f is solved with Neuman boundaries conditions and f integrating the source terms of the Navier-Stokes equations) and avec_sources_et_operateurs (lapP=f is solved as with the previous option avec_sources but f integrating also some operators of the Navier-Stokes equations). The two last options are useful and sometime necessary when source terms are implicited when using an implicit time scheme to solve the Navier-Stokes equations.
}

int Navier_Stokes_std::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="diffusion")
    {
      Cerr << "Reading and typing of the diffusion operator : " << finl;
      terme_diffusif.associer_diffusivite(diffusivite_pour_transport());
      is >> terme_diffusif;
      // the field for dt_stab is the same as that of the operator
      terme_diffusif.associer_diffusivite_pour_pas_de_temps(diffusivite_pour_pas_de_temps());
      return 1;
    }
  else if (mot=="convection")
    {
      Cerr << "Reading and typing of the convection operator : " << finl;
      const Champ_base& vitesse_transportante = vitesse_pour_transport();
      terme_convectif.associer_vitesse(vitesse_transportante);
      is >> terme_convectif;
      return 1;
    }
  else if (mot=="solveur_pression")
    {
      Cerr << "Reading and typing of pressure solver : " << finl;
      is >> solveur_pression_;
      Cerr<<"Pressure solver type : "<<solveur_pression_->que_suis_je()<< finl;
      solveur_pression_.nommer("solveur_pression");
      return 1;
    }
  else if (mot=="dt_projection")
    {
      Cerr << "Reading projection time step " << finl;
      is >> dt_projection;
      is >> seuil_projection;
      return 1;
    }
  else if (mot=="Traitement_particulier")
    {
      Cerr << "Reading and typing of Traitement_particulier considered : " << finl;
      Nom type="Traitement_particulier_NS_";
      Motcle motbidon;
      Motcle accouverte = "{" , accfermee = "}" ;
      is >> motbidon;
      if (motbidon == accouverte)
        {
          Motcle le_cas;
          is >> le_cas;
          if (le_cas == accfermee)
            le_cas ="";
          else
            {
              type+= le_cas;
              type+= "_";
            }
          Nom discr=discretisation().que_suis_je();
          if (discr == "VEFPreP1B")
            discr = "VEF";
          type+=discr;
          Cerr << type << finl;
          le_traitement_particulier.typer(type);
          le_traitement_particulier->associer_eqn(*this);
          le_traitement_particulier->lire(is);
        }
      else
        {
          Cerr << "Error while reading Traitement_particulier for Navier_Stokes_standard equation";
          Cerr << "A { was expected." << finl;
          exit();
        }
      return 1;
    }
  else if (mot=="Erreur_max_DivU")
    {
      Cerr << "Reading DivU maximum" << finl;
      is >> max_div_U;
      is >> seuil_projection;
      return 1;
    }
  else if (mot=="seuil_divU")
    {
      Cerr << "Reading the threshold value for the velocity divergence " << finl;
      is >> seuil_divU;
      is >> raison_seuil_divU;
      return 1;
    }
  else if (mot=="solveur_bar")
    {
      Motcle motlu;
      Motcle accouverte = "{" , accfermee = "}" ;
      Nom type_solv("");
      int acc=0;
      int ok=0;
      while (acc!=0 || !ok)
        {
          is >> motlu;
          if (motlu==accouverte)
            {
              ok=1;
              acc++;
            };
          if (motlu==accfermee) acc--;
          type_solv+=" ";
          type_solv+=motlu;
        }
      if (je_suis_maitre())
        {
          SFichier s("solveur.bar");
          s<<type_solv<<finl;
          s.close();
        }
      return 1;
    }
  else if (mot=="methode_calcul_pression_initiale")
    {
      Motcle methode;
      is >> methode;
      Motcles compris(4);
      compris[0]="avec_les_cl";
      compris[1]="avec_sources";
      compris[2]="avec_sources_et_operateurs";
      compris[3]="sans_rien";
      methode_calcul_pression_initiale_=compris.rang(methode);
      if (methode_calcul_pression_initiale_<0)
        {
          Cerr<<methode<<" is not understood."<<finl;
          Cerr<<" Allowed keywords are :"<<compris<<finl;
          exit();
        }
      return 1;
    }
  else
    return Equation_base::lire_motcle_non_standard(mot,is);
}

const Champ_Don_base& Navier_Stokes_std::diffusivite_pour_transport() const
{
  return fluide().viscosite_cinematique();
}

const Champ_base& Navier_Stokes_std::diffusivite_pour_pas_de_temps() const
{
  return terme_diffusif.diffusivite();
}

const Champ_base& Navier_Stokes_std::vitesse_pour_transport() const
{
  return la_vitesse;
}


/*! @brief Associates with the problem: calls Equation_base::associer_pb_base(const Probleme_base&)
 *
 *      and associates with the divergence and gradient operators.
 *
 * @param (Probleme_base& pb) the problem to associate with
 */
void Navier_Stokes_std::associer_pb_base(const Probleme_base& pb)
{
  Equation_base::associer_pb_base(pb);
  divergence.associer_eqn(*this);
  gradient.associer_eqn(*this);
}

/*! @brief Completes the base equation, associates the pressure with the equation,
 *
 *     completes the divergence, gradient and pressure solver.
 *     Adds 2 source terms: one representing the centrifugal force
 *     in the axisymmetric case, the other involved in the resolution
 *     in 2D axisymmetric.
 *     Associates an interface transport equation with the set
 *     of tracked points if the fluid is marked.
 *
 */
void Navier_Stokes_std::completer()
{
  if (axi == 1)
    {
      Source t;
      Source& so=les_sources.add(t);
      Cerr << "Centrifuge force term creation for Axi case."<< finl;
      Nom type_so = "Force_Centrifuge_";
      Nom disc = discretisation().que_suis_je();
      Nom champ = inconnue().que_suis_je();
      champ.suffix("Champ_");
      type_so+=disc;
      type_so+="_";
      type_so+=champ;
      type_so += "_Axi";
      so.typer_direct(type_so);
      so->associer_eqn(*this);
    }

  Equation_base::completer();

  // Do not split the residual by components for the velocity:
  initialise_residu(1);

  la_pression->associer_eqn(*this);
  la_pression->completer(le_dom_Cl_dis.valeur());
  // [ABN] make sure the pressure knows the domaine_Cl_dis to be able to use specific postreatment like 'gravcl'
  la_pression->associer_domaine_cl_dis(le_dom_Cl_dis);

  divergence_U->associer_eqn(*this);
  if (gradient_P) gradient_P->associer_eqn(*this);
  la_pression_en_pa->associer_eqn(*this);
  la_pression_en_pa->completer(le_dom_Cl_dis.valeur());
  la_pression_en_pa->associer_domaine_cl_dis(le_dom_Cl_dis);
  divergence.completer();
  gradient.completer();
  assembleur_pression_->associer_domaine_cl_dis_base(domaine_Cl_dis());
  assembleur_pression_->completer(*this);

  if (distance_paroi_globale)// Initialize the wall distance at the start of the computation if needed; it will no longer be updated since the mesh is fixed; done late because boundary conditions must be read first.
    {
      Domaine_dis_base& domaine = domaine_dis();
      domaine.init_dist_paroi_globale(domaine_Cl_dis().les_conditions_limites());
      Cerr << "Initializing distance_paroi_globale ... " << finl;
      const DoubleTab& dist_calc = domaine.y_elem();
      for (int e = 0 ; e < domaine.nb_elem() ; e++) distance_paroi_globale->valeurs()(e, 0) = dist_calc(e);
      distance_paroi_globale->valeurs().echange_espace_virtuel();
    }

}

int Navier_Stokes_std::verif_Cl() const
{
  return Equation_base::verif_Cl();
}

/*! @brief Discretizes the equation.
 *
 */
void Navier_Stokes_std::discretiser()
{
  Cerr << "Hydraulic equation discretization (Navier_Stokes_std::discretiser)" << finl;
  const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());

  discretiser_vitesse();
  la_vitesse->add_synonymous(Nom("velocity"));
  champs_compris_.ajoute_champ(la_vitesse);

  dis.pression(schema_temps(), domaine_dis(), la_pression);
  la_pression->add_synonymous(Nom("P_star"));
  champs_compris_.ajoute_champ(la_pression);

  dis.pression_en_pa(schema_temps(), domaine_dis(), la_pression_en_pa);
  la_pression_en_pa->add_synonymous(Nom("Pressure"));
  champs_compris_.ajoute_champ(la_pression_en_pa);


  dis.divergence_U(schema_temps(), domaine_dis(), divergence_U);
  discretiser_grad_p();
  divergence.typer();
  divergence.l_op_base().associer_eqn(*this);
  gradient.typer();
  gradient.l_op_base().associer_eqn(*this);


  champs_compris_.ajoute_champ(divergence_U);


  // Call the virtual method for pressure assembler discretization:
  discretiser_assembleur_pression();

  Equation_base::discretiser();
}

void Navier_Stokes_std::discretiser_vitesse()
{
  const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
  dis.vitesse(schema_temps(), domaine_dis(), la_vitesse);
}

void Navier_Stokes_std::discretiser_grad_p()
{
  const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
  dis.gradient_P(schema_temps(), domaine_dis(), gradient_P);
  champs_compris_.ajoute_champ(gradient_P);
}

/*! @brief Types the pressure assembler.
 *
 * The assembler name is built as:
 *   "Assembleur_P_xxx" where "xxx" is the name of the discretization.
 *   This method is virtual and overridden in front-tracking.
 *   It is called by Navier_Stokes_std::discretiser()
 *
 */

void Navier_Stokes_std::discretiser_assembleur_pression()
{
  Nom type = "Assembleur_P_";
  type += discretisation().que_suis_je();
  Cerr << "Navier_Stokes_std::discretiser_assembleur_pression : type="<< type << finl;
  assembleur_pression_.typer(type);
  assembleur_pression_->associer_domaine_dis_base(domaine_dis());
}

void Navier_Stokes_std::reassembler_pression_si_necessaire()
{
  if (!probleme().domaine().mesh_update_required()) return;

  if (!probleme().is_dilatable())
    assembleur_pression_->assembler(matrice_pression_);
  else
    assembleur_pression_->assembler_QC(fluide().masse_volumique().valeurs(), matrice_pression_);

  solveur_pression_->reinit();
}

/*! @brief Returns the number of operators in the equation: 2 for the standard Navier-Stokes equations.
 *
 * @return (int) the number of operators in the equation
 */
int Navier_Stokes_std::nombre_d_operateurs() const
{
  return 2;
}

int Navier_Stokes_std::nombre_d_operateurs_tot() const
{
  return 4;
}

/*! @brief Returns the i-th operator of the equation: - terme_diffusif if i = 0
 *
 *       - terme_convectif if i = 1
 *      exits if i > 1
 *      (const version)
 *
 * @param (int i) the index of the operator to return
 * @return (Operateur&) the operator at index i
 */
const Operateur& Navier_Stokes_std::operateur(int i) const
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    case 1:
      return terme_convectif;
    default :
      Cerr << "Error for Navier_Stokes_std::operateur(int i)" << finl;
      Cerr << "Navier_Stokes_std has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // For compilers!
  return terme_diffusif;
}

/*! @brief Returns the i-th operator of the equation: - terme_diffusif if i = 0
 *
 *       - terme_convectif if i = 1
 *      exits if i > 1
 *
 * @param (int i) the index of the operator to return
 * @return (Operateur&) the operator at index i
 */
Operateur& Navier_Stokes_std::operateur(int i)
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    case 1:
      return terme_convectif;
    default :
      Cerr << "Error for Navier_Stokes_std::operateur(int i)" << finl;
      Cerr << "Navier_Stokes_std has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // For compilers!
  return terme_diffusif;
}

const Operateur& Navier_Stokes_std::operateur_fonctionnel(int i) const
{
  switch(i)
    {
    case 0:
      return gradient;
    case 1:
      return divergence;
    default :
      Cerr << "Error for Navier_Stokes_std::operateur_fonctionnel(int i)" << finl;
      Cerr << "Navier_Stokes_std has " << nombre_d_operateurs() <<" functional operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // For compilers!
  return Navier_Stokes_std::operateur_fonctionnel(i);
}

Operateur& Navier_Stokes_std::operateur_fonctionnel(int i)
{
  switch(i)
    {
    case 0:
      return gradient;
    case 1:
      return divergence;
    default :
      Cerr << "Error for Navier_Stokes_std::operateur_fonctionnel(int i)" << finl;
      Cerr << "Navier_Stokes_std has " << nombre_d_operateurs() <<" functional operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // For compilers!
  return Navier_Stokes_std::operateur_fonctionnel(i);
}


/*! @brief Returns the divergence operator associated with the equation.
 *
 * @return (Operateur_Div&) the divergence operator.
 */
Operateur_Div& Navier_Stokes_std::operateur_divergence()
{
  return divergence;
}

/*! @brief Returns the divergence operator associated with the equation.
 *
 *     (const version)
 *
 * @return (Operateur_Div&) the divergence operator
 */
const Operateur_Div& Navier_Stokes_std::operateur_divergence() const
{
  return divergence;
}

/*! @brief Returns the gradient operator associated with the equation.
 *
 * @return (Operateur_Grad&) the gradient operator
 */
Operateur_Grad& Navier_Stokes_std::operateur_gradient()
{
  return gradient;
}

Operateur_Diff& Navier_Stokes_std::operateur_diff()
{
  return terme_diffusif;
}

const Operateur_Diff& Navier_Stokes_std::operateur_diff() const
{
  return terme_diffusif;
}

/*! @brief Returns the gradient operator associated with the equation.
 *
 *     (const version)
 *
 * @return (Operateur_Grad&) the gradient operator
 */
const Operateur_Grad& Navier_Stokes_std::operateur_gradient() const
{
  return gradient;
}


/*! @brief Returns the velocity (unknown field of the equation) (const version)
 *
 * @return (Champ_Inc_base&) the unknown field representing the velocity
 */
const Champ_Inc_base& Navier_Stokes_std::inconnue() const
{
  return la_vitesse.valeur();
}

/*! @brief Returns the velocity (unknown field of the equation)
 *
 * @return (Champ_Inc_base&) the unknown field representing the velocity
 */
Champ_Inc_base& Navier_Stokes_std::inconnue()
{
  return la_vitesse.valeur();
}

/*! @brief Returns the pressure solver (const version)
 *
 * @return (SolveurSys&) the pressure solver
 */
SolveurSys& Navier_Stokes_std::solveur_pression()
{
  return solveur_pression_;
}

/*! @brief Returns the incompressible fluid (physical medium of the equation) associated with the equation.
 *
 *     (const version)
 *
 * @return (Fluide_base&) the incompressible fluid associated with the equation
 */
const Fluide_base& Navier_Stokes_std::fluide() const
{
  return le_fluide.valeur();
}


/*! @brief Returns the incompressible fluid (physical medium of the equation) associated with the equation.
 *
 * @return (Fluide_base&) the incompressible fluid associated with the equation
 */
Fluide_base& Navier_Stokes_std::fluide()
{
  return le_fluide.valeur();
}

Entree& Navier_Stokes_std::lire_cond_init(Entree& is)
{
  Cerr << "Reading of initial conditions\n";
  Nom nom;
  Motcle motlu;
  is >> nom;
  motlu = nom;
  if(motlu!=Motcle("{"))
    {
      Cerr << "We expected a { while reading " << que_suis_je() << finl;
      Cerr << "and not : " << nom << finl;
      exit();
    }
  Motcles compris(3);
  compris[0]="}";
  compris[1]="vitesse";
  compris[2]="pression";
  int ind = -1;
  while (ind!=0)
    {
      is >> nom;
      motlu = nom;
      ind = compris.rang(motlu);
      if (ind==1)
        {
          OWN_PTR(Champ_Don_base) ch_init;
          is >> ch_init;
          verifie_ch_init_nb_comp(inconnue(),ch_init->nb_comp());
          inconnue().affecter(ch_init.valeur());
        }
      else if (ind==2)
        {
          OWN_PTR(Champ_Don_base) ch_init;
          is >> ch_init;
          verifie_ch_init_nb_comp(pression(),ch_init->nb_comp());
          pression().affecter(ch_init.valeur());
        }
      else if (ind==-1)
        {
          Cerr << nom << " is not understood. Keywords are:" << finl;
          Cerr << compris << finl;
          exit();
        }
    }
  return is;
}

/*! @brief Add a specific term for Navier Stokes (-gradP(n)) if necessary
 *
 */
DoubleTab& Navier_Stokes_std::corriger_derivee_expl(DoubleTab& derivee)
{
  if (assembleur_pression_->get_resoudre_increment_pression())
    {
      // PL: To avoid computing this gradient, the following must hold:
      // A) postraitement_gradient_P_==0 car sinon grad contient alors M-1BtP
      // B) pressure boundary conditions are stationary (hard to detect: Orlansky, P(t), imposed pressure gradient...)
      // Moreover, this causes discrepancies with the iterative CN scheme
      const DoubleTab& tab_pression = la_pression->valeurs();
      DoubleTab& gradP = gradient_P->valeurs();
      gradient.calculer(tab_pression, gradP);
      derivee -= gradP;
    }
  return derivee;
}

/*! @brief Resolution de la pression, inconnue implicitee de Navier Stokes
 *
 */
DoubleTab& Navier_Stokes_std::corriger_derivee_impl(DoubleTab& derivee)
{
  // We want to solve:
  // dU/dt + M-1 Bt Cp = M-1(F - BtP)
  // B dU/dt = 0
  // with F explicit terms: sum(operators)+sources
  // In:  derivee = M-1(F - BtP(n))
  // Out: derivee = M-1(F - BtP(n+1)), P(n+1)=P(n)+Cp

  DoubleTab& tab_pression=la_pression->valeurs();
  DoubleTab& gradP=gradient_P->valeurs();
  DoubleTrav secmemP(tab_pression);

  const bool is_ALE = probleme().domaine().deformable();

  if (div_u_nul_et_non_dsurdt_divu_ && is_ALE) Process::exit();

  const double dt = schema_temps().pas_de_temps();
  if (div_u_nul_et_non_dsurdt_divu_)
    {
      // we want div u = 0 rather than d/dt(div u) = 0 to avoid accumulating errors
      // this only works with explicit Euler-type schemes
      DoubleTab derivee2(derivee);
      derivee2*=dt;
      derivee2+=la_vitesse->passe();
      derivee2/=dt;
      divergence.calculer(derivee2, secmemP); // Div(M-1(F - BtP))
    }
  else if (is_ALE)
    {
      // ALE explicit volumetric correction:
      // In ALE formulation, the conservative form of the transient term is:
      // (V^{n+1} U^{n+1} - V^n U^n ) / dt
      // For an explicit Euler scheme, this leads to the additional term:
      // (V^n / V^{n+1}) * (U^n / dt) which is added to the explicit RHS (M-1(F - BtP)).
      // This ensures discrete mass conservation (Geometric Conservation Law) when the mesh is deformable.
      // This correction is required only for explicit schemes because
      // implicit schemes naturally incorporate the volume variation inside the mass matrix formulation.
      DoubleTab derivee2(derivee);
      // Conservative ALE form: deriveeALE= (Volume_n+1/Volume_n)*Un/delta t + derivee
      probleme().domaine().ajouter_correctif_volumique(la_vitesse->valeurs(), derivee, dt, derivee2);
      divergence.calculer(derivee2, secmemP);
    }
  else
    divergence.calculer(derivee, secmemP); // Div(M-1(F - BtP))

  secmemP *= -1; // because div = -B
  // Correction of the right-hand side according to the boundary conditions:
  assembleur_pression_->modifier_secmem(secmemP);

  // Set print of the linear system solve according to dt_impr:
  solveur_pression_->fixer_schema_temps_limpr(schema_temps().limpr());

  if (assembleur_pression_->get_resoudre_increment_pression())
    {
      // Solve B M-1 Bt Cp = M-1(F - BtP)
      DoubleTrav Cp(tab_pression);
      solveur_pression_.resoudre_systeme(matrice_pression_.valeur(), secmemP, Cp);

      // P(n+1) = P(n) + Cp
      tab_pression += Cp;
      assembleur_pression_->modifier_solution(tab_pression);

      // M-1 Bt P(n+1)
      solveur_masse->appliquer(gradP);
      derivee += gradP; // M-1 F
    }
  else
    {
      // Solve B M-1 Bt P(n+1) = B M-1 F
      solveur_pression_.resoudre_systeme(matrice_pression_.valeur(), secmemP, tab_pression);
      assembleur_pression_->modifier_solution(tab_pression);
      // It is not done anymore cause:
      // Iterative solvers are less accurate
      // Time converges in O(sqrt(dt)) and not O(dt)
      // See: http://www.sciencedirect.com/science/article/pii/S0021999108004518
    }

  // (BM) gradient operator requires updated virtual space in source vector
  // Calculate Bt P(n+1)
  tab_pression.echange_espace_virtuel();
  gradient.calculer(tab_pression, gradP);

  // gradP = Bt P(n+1) is kept and
  // M-1Bt P(n+1) is calculated:
  DoubleTrav Mmoins1gradP(gradP);
  Mmoins1gradP = gradP;
  solveur_masse->appliquer(Mmoins1gradP);

  // dU/dt = M-1(F-Bt P(n+1))
  derivee -= Mmoins1gradP;

  return derivee;
}

/*! @brief Computes the solution U of the equations: | M(U-V)/dt + BtP = 0
 *
 *           |-BU=0
 *     The pressure problem is solved: -BM-1BtP = -BV/dt
 *     where -BV represents the divergence of V.
 *     The velocity problem is solved by applying the mass solver
 *     to the gradient of P:  U=V - dt*M-1BtP
 *
 * @throws time step too small
 */
void Navier_Stokes_std::projeter()
{
  if (probleme().is_dilatable() && probleme().reprise_effectuee())
    Cerr << "WARNING: Quasi compressible model --> no projection (except the first time step)." << finl;
  else
    {
      Cerr << "Projection of initial and boundaries conditions " << finl;
      DoubleTab& tab_vitesse = la_vitesse->valeurs();
      tab_vitesse.echange_espace_virtuel();
      la_pression->valeurs().echange_espace_virtuel();

      double normal_seuil = 0.;

      //  M u + eBt l = M v
      //  B u = 0
      //  => e B(M-1)Bt l = Bv
      //
      modify_initial_variable();

      DoubleTrav secmem(la_pression->valeurs());
      divergence.calculer(tab_vitesse, secmem);
      // The time step is now computed before the projection.
      // Previously, dt=dt_min at the start of the computation.
      double dt = std::max(le_schema_en_temps->pas_temps_min(),calculer_pas_de_temps());
      dt = std::min(dt, le_schema_en_temps->pas_temps_max());

      secmem*=(-1./dt);
      secmem.echange_espace_virtuel();

      double bilan=mp_norme_vect(secmem); // TODO DG surcharger cette fonction pour avoir \sum \int_T || \sum_{nfunc_p} secmem*fbase_p ||
      Cout << "-------------  Projection  -----------------" << finl;
      Cout << "--------------------------------------------" << finl;
      Cout << "Bilan de masse avant projection : " << bilan << finl;

      if( sub_type(solv_iteratif,solveur_pression_.valeur()) )
        {
          solv_iteratif& solv_iter=ref_cast(solv_iteratif,solveur_pression_.valeur());
          normal_seuil=solv_iter.get_seuil();
          solv_iter.set_seuil(seuil_projection);
        }

      // Correction of the right-hand side according to the boundary conditions:
      // There is no point initializing lagrange with the pressure;
      // it actually penalizes the computation in p1B with CL p<>0.
      // Use a DoubleTrav instead of a DoubleTab to initialize lagrange=0.
      DoubleTrav lagrange(la_pression->valeurs());
      solveur_pression_.resoudre_systeme(matrice_pression_.valeur(),secmem,lagrange);
      assembleur_pression_->modifier_solution(lagrange);
      lagrange.echange_espace_virtuel();

      // M-1 Bt l
      DoubleTrav gradP(gradient_P->valeurs());
      gradient->multvect(lagrange, gradP);
      modify_initial_gradP(gradP);
      gradP.echange_espace_virtuel();
      solveur_masse->appliquer(gradP);
      gradP.echange_espace_virtuel();

      if (tab_vitesse.dimension_tot(0) == gradP.dimension_tot(0))
        tab_vitesse.ajoute(-dt,gradP);
      else
        {
          DoubleTab_parts partv(tab_vitesse);
          partv[0].ajoute(-dt,gradP);
        }
      tab_vitesse.echange_espace_virtuel();
      solveur_masse->corriger_solution(tab_vitesse, tab_vitesse);

      Debog::verifier("Navier_Stokes_std::projeter, vitesse", tab_vitesse);

      // Verif ...
      divergence.calculer(tab_vitesse, secmem);
      secmem.echange_espace_virtuel();

      bilan=mp_norme_vect(secmem);
      Cout << "Bilan de masse apres projection : " << bilan << finl;
      Cout << "-------------  Projection  OK---------------" << finl;
      Cout << "--------------------------------------------" << finl;

      if( sub_type(solv_iteratif,solveur_pression_.valeur()) )
        {
          solv_iteratif& solv_iter=ref_cast(solv_iteratif,solveur_pression_.valeur());
          solv_iter.set_seuil(normal_seuil);
        }
    }
  projection_initiale = 0;
}

int Navier_Stokes_std::projection_a_faire()
{
  // No projection if the equation is not solved
  // if (equation_non_resolue()) return 0;

  double temps = le_schema_en_temps->temps_courant()+le_schema_en_temps->pas_de_temps();
  // See Schema_Temps_base::limpr for information on modf
  double nb_proj_int;
  modf(temps/dt_projection, &nb_proj_int);
  static double nb_proj = nb_proj_int;
  if (projection_initiale)
    return 1;
  else if (inf_ou_egal((nb_proj+1.)*dt_projection,temps))
    {
      nb_proj=nb_proj+1.;
      return 1;
    }
  return 0;
}

/*! @brief cf Equation_base::preparer_calcul() Assembly of the pressure solver and
 *
 *      initialization of the pressure.
 *
 *      Assembly of the pressure system.
 *
 * @return (int) always returns 1
 */
int Navier_Stokes_std::preparer_calcul()
{
  const double temps = schema_temps().temps_courant();
  sources().mettre_a_jour(temps);
  Equation_base::preparer_calcul();
  bool is_dilatable = probleme().is_dilatable();
  if (!is_dilatable)
    assembleur_pression_->assembler(matrice_pression_);
  else
    {
      Cerr << "Assembling for quasi-compressible" << finl;
      assembleur_pression_->assembler_QC(fluide().masse_volumique().valeurs(), matrice_pression_);
    }

  // GF en cas de reprise on conserve la valeur de la pression
  // avant elle ne servait qu' a initialiser le lagrange pour la projection
  // C'est important pour le Simpler/Piso de bien repartir de la pression
  // sauvegardee...
  //la_pression->valeurs()=0.;
  Debog::verifier("Navier_Stokes_std::preparer_calcul, la_pression av projeter", la_pression->valeurs());
  if (projection_a_faire())
    projeter();

  // Au cas ou une cl de pression depend de u que l'on vient de modifier
  le_dom_Cl_dis->mettre_a_jour(temps);
  Debog::verifier("Navier_Stokes_std::preparer_calcul, la_pression ap projeter", la_pression->valeurs());

  // Initialisation du champ de pression (resolution de Laplacien(P)=0 avec les conditions limites en pression)
  // Permet de demarrer la resolution avec une bonne approximation de la pression (important pour le Piso ou P!=0)
  if (!probleme().reprise_effectuee() && methode_calcul_pression_initiale_ != 3)
    {
      Cout << "Estimation du champ de pression au demarrage:" << finl;
      DoubleTrav secmem(la_pression->valeurs());
      DoubleTrav vpoint(gradient_P->valeurs());
      gradient.calculer(la_pression->valeurs(), gradient_P->valeurs());
      vpoint -= gradient_P->valeurs();

      if (methode_calcul_pression_initiale_ >= 2)
        for (int op = 0; op < nombre_d_operateurs(); op++)
          operateur(op).ajouter(vpoint);
      if (methode_calcul_pression_initiale_ >= 1)
        {
          int mod = 0;
          if (le_schema_en_temps->pas_de_temps() == 0)
            {
              double dt = std::max(le_schema_en_temps->pas_temps_min(), calculer_pas_de_temps());
              dt = std::min(dt, le_schema_en_temps->pas_temps_max());
              le_schema_en_temps->set_dt() = (dt);
              mod = 1;
            }
          sources().ajouter(vpoint);
          if (mod)
            le_schema_en_temps->set_dt() = 0;
        }

      solveur_masse->appliquer(vpoint);
      vpoint.echange_espace_virtuel();
      divergence.calculer(vpoint, secmem);
      secmem *= -1;
      secmem.echange_espace_virtuel();

      assembleur_pression_->modifier_secmem_pour_incr_p(la_pression->valeurs(), 1, secmem);
      DoubleTrav inc_pre(la_pression->valeurs());
      solveur_pression_.resoudre_systeme(matrice_pression_.valeur(), secmem, inc_pre);
      Cerr << "Pressure increment computed successfully" << finl;

      // We want the virtual space to be up-to-date, hence all_items
      operator_add(la_pression->valeurs(), inc_pre, VECT_ALL_ITEMS);
    }
  // Update pressure
  la_pression->changer_temps(temps);
  calculer_la_pression_en_pa();
  // Compute pressure forces:
  gradient->calculer_flux_bords();

  // Compute gradient_P (ToDo: make consistent with ::mettre_a_jour()):
  gradient.calculer(la_pression->valeurs(), gradient_P->valeurs());
  gradient_P->changer_temps(temps);

  // Compute divergence_U
  divergence.calculer(la_vitesse->valeurs(), divergence_U->valeurs());
  divergence_U->changer_temps(temps);

  if (le_traitement_particulier)
    le_traitement_particulier->preparer_calcul_particulier();

  Debog::verifier("Navier_Stokes_std::preparer_calcul, vitesse", inconnue());
  Debog::verifier("Navier_Stokes_std::preparer_calcul, pression", la_pression);

  return 1;
}

/*! @brief Performs a time update of the equation.
 *
 * Calls Equation_base::mettre_a_jour(double)
 *      and updates the pressure.
 *      Integrates tracked points if the fluid is marked.
 *      Updates the corresponding post-processable field.
 *
 * @param (double temps) the update time
 */
void Navier_Stokes_std::mettre_a_jour(double temps)
{
  // Update the parent class (advance the time wheel).
  Equation_base::mettre_a_jour(temps);

  // Update pressure
  la_pression->mettre_a_jour(temps);
  calculer_la_pression_en_pa();
  // Compute pressure forces:
  gradient->calculer_flux_bords();

  // Update the divergence of the velocity div(U)
//  statistics().end_count(STD_COUNTERS::update_variables,0,0);
  divergence.calculer(la_vitesse->valeurs(),divergence_U->valeurs());
  //statistics().begin_count(STD_COUNTERS::update_variables,statistics().get_last_opened_counter_level()+1);
  divergence_U->mettre_a_jour(temps);

  // For post-processing, we want M-1BtP and not BtP
  if (postraitement_gradient_P_)
    {
      gradient.calculer(la_pression->valeurs(), gradient_P->valeurs());
      if (!postraiter_gradient_pression_sans_masse_)
        solveur_masse->appliquer(gradient_P->valeurs());
      gradient_P->mettre_a_jour(temps);
    }

  // PQ: 04/03: procedure for dynamic determination of the pressure convergence threshold
  if(sub_type(solv_iteratif,solveur_pression_.valeur()) && seuil_divU < 1.)
    {
      // Dynamic threshold computation for the iterative pressure solver
      solv_iteratif& solv_iter=ref_cast(solv_iteratif,solveur_pression_.valeur());
      double seuil_dyn=solv_iter.get_seuil();

      if(LocalFlowRateRelativeError()<seuil_divU)
        seuil_dyn*=raison_seuil_divU;
      else
        seuil_dyn/=raison_seuil_divU;
      double seuil_dyn_max = 1.e-10;
      seuil_dyn=std::max(seuil_dyn,seuil_dyn_max);
      solv_iter.set_seuil(seuil_dyn);
    }
  // end of the dynamic pressure convergence threshold procedure

  if (projection_a_faire())
    projeter();

  if (le_traitement_particulier)
    le_traitement_particulier->post_traitement_particulier();
  Debog::verifier("Navier_Stokes_std::mettre_a_jour : pression", la_pression->valeurs());
  Debog::verifier("Navier_Stokes_std::mettre_a_jour : vitesse", la_vitesse->valeurs());

  if (la_vorticite) la_vorticite->mettre_a_jour(temps);
  if (critere_Q) critere_Q->mettre_a_jour(temps);
  if (Reynolds_maille) Reynolds_maille->mettre_a_jour(temps);
  if (Taux_cisaillement) Taux_cisaillement->mettre_a_jour(temps);
  if (grad_u) grad_u->mettre_a_jour(temps);
}

double Navier_Stokes_std::LocalFlowRateRelativeError() const
{
  // Estimation of a flow rate relative error
  DoubleTrav array(divergence_U->valeurs()); // array(i)=sum(u.ndS)
  divergence.volumique(array); // array(i)=sum(u.ndS)/vol(i)
  return mp_max_abs_vect(array) * schema_temps().pas_de_temps(); // =max|sum(u.ndS)/(vol(i)/dt)|
}

void Navier_Stokes_std::abortTimeStep()
{
  // Restore the pressure from the beginning of the time step.
  // Useful if the time step is retried because the pressure diverged (otherwise very poor accuracy)
  // and when using Piso (which assumes correct pressure at the start of the time step).
  pression().valeurs()=P_n;
  //pression().valeurs()=0;
  Equation_base::abortTimeStep();
}

/* @brief Override. Reset pression too !
 */
void Navier_Stokes_std::resetTime(double time)
{
  pression().resetTime(time);
  Equation_base::resetTime(time);
}

bool Navier_Stokes_std::initTimeStep(double dt)
{
  P_n=pression().valeurs();

  // Check that dt_max is correctly set for a zero
  // velocity field with active implicit diffusion <=> dt_conv=INF
  const Schema_Temps_base& sch_tps = le_schema_en_temps.valeur();
  bool ddt = Equation_base::initTimeStep(dt);

  for (int i=1; i<=sch_tps.nb_valeurs_futures(); i++)
    if (i <= pression().nb_valeurs_temporelles())
      {
        double tps=sch_tps.temps_futur(i);
        // Update time in pressure fields
        pression().changer_temps_futur(tps,i);
        pression_pa().changer_temps_futur(tps,i);
        pression().futur(i)=pression().valeurs();
        pression_pa().futur(i)=pression_pa().valeurs();
      }

  return ddt;
}

/*! @brief Computes "la_pression_en_pa" from "la_pression".
 *
 * If the field milieu().masse_volumique() is uniform, it is assumed that
 *   la_pression is P* = P/rho, and the value is multiplied by rho. Otherwise,
 *   la_pression is already in Pa.
 *   This method is overridden in front-tracking.
 *
 */
void Navier_Stokes_std::calculer_la_pression_en_pa()
{
  DoubleTab& Pa=la_pression_en_pa->valeurs();
  DoubleTab& tab_pression=la_pression->valeurs();
  const Champ_base& rho=milieu().masse_volumique();
  if (Pa.get_md_vector() == tab_pression.get_md_vector())
    Pa = tab_pression; //Pa and tab_pression share the same support
  else
    {
      ConstDoubleTab_parts ppart(tab_pression);
      assert(Pa.get_md_vector() == ppart[0].get_md_vector());
      Pa = ppart[0]; //tab_pression has one extra piece
    }
  // Multiply by rho if uniform, otherwise already in Pa...
  if (sub_type(Champ_Uniforme,rho))
    Pa *= rho.valeurs()(0,0);
  la_pression_en_pa->mettre_a_jour(pression().temps());
}

/*! @brief for PDI IO: retrieve name, type and dimensions of the fields to save/restore
 *
 */
std::vector<YAML_data> Navier_Stokes_std::data_a_sauvegarder() const
{
  std::vector<YAML_data> data = Equation_base::data_a_sauvegarder();
  std::vector<YAML_data> pression = la_pression->data_a_sauvegarder();
  data.insert(data.end(), pression.begin(), pression.end());
  return data;
}

/*! @brief Calls Equation_base::sauvegarder(Sortie&) and saves the pressure to an output stream.
 *
 * @param (Sortie& os) output stream to save to
 * @return (int) always returns 1
 */
int Navier_Stokes_std::sauvegarder(Sortie& os) const
{
  int bytes=0;
  bytes += Equation_base::sauvegarder(os);
  bytes += la_pression->sauvegarder(os);
  // The sauver() method ensuring save for the particular treatment
  // is now called here instead of being called in specific problems.
  sauver();

  return bytes;
}

/*! @brief Performs a restart from an input stream.
 *
 * Calls Equation_base::reprendre()
 *      and restores the pressure.
 *
 * @param (Entree& is) input stream
 * @return (int) always returns 1
 * @throws restart failed, pressure identifier not found
 */
int Navier_Stokes_std::reprendre(Entree& is)
{
  Equation_base::reprendre(is);
  if(!TRUST_2_PDI::is_PDI_restart())
    {
      double temps = schema_temps().temps_courant();
      Nom ident_pression(la_pression->le_nom());
      ident_pression += la_pression->que_suis_je();
      ident_pression += probleme().domaine().le_nom();
      ident_pression += Nom(temps,probleme().reprise_format_temps());
      if (probleme().discretisation().is_poly_family())
        {
          Nom field_tag_syno = create_polymacfamily_syno(ident_pression);
          avancer_fichier_with_syno(is,ident_pression,field_tag_syno);
        }
      // end of the backward compatibility
      else
        avancer_fichier(is,ident_pression);
    }
  la_pression->reprendre(is);

  if (le_traitement_particulier)
    le_traitement_particulier->reprendre_stat();

  return 1;
}

/*! @brief Associates a physical medium with the equation by dynamically constructing (casting) an object of type Fluide_base
 *
 *     from the Milieu_base object passed as parameter.
 *
 * @param (Milieu_base& un_milieu) the medium to associate with the equation
 */
void Navier_Stokes_std::associer_milieu_base(const Milieu_base& un_milieu)
{
  if (sub_type(Fluide_base, un_milieu))
    {
      const Fluide_base& un_fluide = ref_cast(Fluide_base,un_milieu);
      associer_fluide(un_fluide);
    }
  else
    {
      Cerr << "Error of fluid type for the method Navier_Stokes_std::associer_milieu_base" << finl;
      exit();
    }
}

/*! @brief Returns the physical medium of the equation (Fluide_base upcast to Milieu_base)
 *
 * @return (Milieu_base&) the Fluide_base of the equation upcast to Milieu_base
 */
const Milieu_base& Navier_Stokes_std::milieu() const
{
  if (!le_fluide)
    {
      Cerr << "You forgot to associate a fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide.valeur();
}

/*! @brief Returns the physical medium of the equation (Fluide_base upcast to Milieu_base)
 *
 *     (const version)
 *
 * @return (Milieu_base&) the Fluide_base of the equation upcast to Milieu_base
 */
Milieu_base& Navier_Stokes_std::milieu()
{
  if (!le_fluide)
    {
      Cerr << "You forgot to associate a fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide.valeur();
}

void Navier_Stokes_std::creer_champ(const Motcle& motlu)
{
  Equation_base::creer_champ(motlu);

  if (motlu == "vorticite")
    {
      if (!la_vorticite)
        {
          const Discret_Thyd& dis=ref_cast(Discret_Thyd,discretisation());
          dis.creer_champ_vorticite(schema_temps(),la_vitesse,la_vorticite);
          champs_compris_.ajoute_champ(la_vorticite);
        }
    }
  else if (motlu == "critere_Q")
    {
      if (!critere_Q)
        {
          const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
          dis.critere_Q(domaine_dis(),domaine_Cl_dis(),la_vitesse,critere_Q);
          champs_compris_.ajoute_champ(critere_Q);
        }
    }
  else if (motlu == "y_plus")
    {
      if (!y_plus)
        {
          const Discret_Thyd& dis=ref_cast(Discret_Thyd,discretisation());
          dis.y_plus(domaine_dis(),domaine_Cl_dis(),la_vitesse,y_plus);
          champs_compris_.ajoute_champ(y_plus);
        }
    }
  else if (motlu == "distance_paroi_globale")
    {
      if (!distance_paroi_globale)
        {
          const Discret_Thyd& dis=ref_cast(Discret_Thyd,discretisation());
          dis.distance_paroi_globale(schema_temps(), domaine_dis(), distance_paroi_globale);
          champs_compris_.ajoute_champ(distance_paroi_globale);
        }
    }
  else if (motlu == "reynolds_maille")
    {
      if (!Reynolds_maille)
        {
          const Discret_Thyd& dis=ref_cast(Discret_Thyd,discretisation());
          dis.reynolds_maille(domaine_dis(),fluide(),la_vitesse,Reynolds_maille);
          champs_compris_.ajoute_champ(Reynolds_maille);
        }
    }
  else if (motlu == "courant_maille")
    {
      if (!Courant_maille)
        {
          const Discret_Thyd& dis=ref_cast(Discret_Thyd,discretisation());
          dis.courant_maille(domaine_dis(),schema_temps(),la_vitesse,Courant_maille);
          champs_compris_.ajoute_champ(Courant_maille);
        }
    }
  else if (motlu == "taux_cisaillement")
    {
      if (!Taux_cisaillement)
        {
          const Discret_Thyd& dis=ref_cast(Discret_Thyd,discretisation());
          dis.taux_cisaillement(domaine_dis(),domaine_Cl_dis(),la_vitesse,Taux_cisaillement);
          champs_compris_.ajoute_champ(Taux_cisaillement);
        }
    }
  else if (motlu == "pression_hydrostatique")
    {
      if (!pression_hydrostatique_)
        {
          const Discret_Thyd& dis=ref_cast(Discret_Thyd,discretisation());
          dis.discretiser_champ("Champ_sommets",domaine_dis(),"pression_hydrostatique","Pa",1,0.,pression_hydrostatique_);
          champs_compris_.ajoute_champ(pression_hydrostatique_);
        }
    }

  else if (motlu == "gradient_vitesse")
    {
      if (!grad_u)
        {
          const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
          dis.grad_u(domaine_dis(),domaine_Cl_dis(),la_vitesse,grad_u);
          champs_compris_.ajoute_champ(grad_u);
        }
    }

  if (le_traitement_particulier)
    le_traitement_particulier->creer_champ(motlu);

  if (Taux_cisaillement)
    if (!grad_u) creer_champ("gradient_vitesse");
}

void Navier_Stokes_std::calculer_pression_hydrostatique(Champ_base& pression_hydro) const
{
  DoubleTab& val= pression_hydro.valeurs();
  const DoubleTab& coords = domaine_dis().domaine().les_sommets();
  if (!milieu().a_gravite())
    {
      Cerr<<"postprocessing of presion_hydrostatique needs gravity"<<finl;
      exit();
    }
  const Champ_base& rho = milieu().masse_volumique();
  if (!sub_type(Champ_Uniforme,rho))
    {
      Cerr<<"postprocessing of presion_hydrostatique availabe only for incompressible flow"<<finl;
      exit();
    }
  const DoubleTab& gravite = milieu().gravite().valeurs();

  val=rho.valeurs()(0,0);
  const int nb_som=val.dimension(0);

  for (int som=0; som<nb_som; som++)
    {
      double gz=0;
      for (int dir=0; dir<dimension; dir++)
        gz+=coords(som,dir)*gravite(0,dir);
      val[som]*=gz;
    }
  val.echange_espace_virtuel();
}

bool Navier_Stokes_std::has_champ(const Motcle& nom, OBS_PTR(Champ_base)& ref_champ) const
{
  if (nom == "gradient_pression")
    {
      ref_champ = Navier_Stokes_std::get_champ(nom);
      return true;
    }

  if (nom == "vorticite" && la_vorticite)
    {
      ref_champ = Navier_Stokes_std::get_champ(nom);
      return true;
    }

  if (nom == "critere_Q" && critere_Q)
    {
      ref_champ = Navier_Stokes_std::get_champ(nom);
      return true;
    }

  if (nom == "y_plus" && y_plus)
    {
      ref_champ = Navier_Stokes_std::get_champ(nom);
      return true;
    }

  if (nom == "reynolds_maille" && Reynolds_maille)
    {
      ref_champ = Navier_Stokes_std::get_champ(nom);
      return true;
    }

  if (nom == "courant_maille" && Courant_maille)
    {
      ref_champ = Navier_Stokes_std::get_champ(nom);
      return true;
    }

  if (nom == "taux_cisaillement" && Taux_cisaillement)
    {
      ref_champ = Navier_Stokes_std::get_champ(nom);
      return true;
    }

  if (nom == "gradient_vitesse" && grad_u)
    {
      ref_champ = Navier_Stokes_std::get_champ(nom);
      return true;
    }

  if (nom == "pression_hydrostatique" && pression_hydrostatique_)
    {
      ref_champ = Navier_Stokes_std::get_champ(nom);
      return true;
    }

  if (Equation_base::has_champ(nom, ref_champ))
    return true;

  if (le_traitement_particulier)
    if (le_traitement_particulier->has_champ(nom, ref_champ))
      return true;

  return false; /* nothing found */
}

bool Navier_Stokes_std::has_champ(const Motcle& nom) const
{
  if (nom == "gradient_pression")
    return true;

  if (nom == "vorticite" && la_vorticite)
    return true;

  if (nom == "critere_Q" && critere_Q)
    return true;

  if (nom == "y_plus" && y_plus)
    return true;

  if (nom == "reynolds_maille" && Reynolds_maille)
    return true;

  if (nom == "courant_maille" && Courant_maille)
    return true;

  if (nom == "taux_cisaillement" && Taux_cisaillement)
    return true;

  if (nom == "gradient_vitesse" && grad_u)
    return true;

  if (nom == "pression_hydrostatique" && pression_hydrostatique_)
    return true;

  if (Equation_base::has_champ(nom))
    return true;

  if (le_traitement_particulier)
    if (le_traitement_particulier->has_champ(nom))
      return true;

  return false; /* nothing found */
}

const Champ_base& Navier_Stokes_std::get_champ(const Motcle& nom) const
{
  double temps_init = schema_temps().temps_init();
  if (nom == "gradient_pression")
    postraitement_gradient_P_ = 1;

  if (nom == "vorticite")
    {
      if (!la_vorticite)
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      Champ_Fonc_base& ch = ref_cast_non_const(Champ_Fonc_base, la_vorticite.valeur());
      if ((ch.temps() == temps_init) && (la_vitesse->mon_equation_non_nul()))
        ch.mettre_a_jour(la_vitesse->temps());
      return champs_compris_.get_champ(nom);
    }

  if (nom == "critere_Q")
    {
      if (!critere_Q)
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      Champ_Fonc_base& ch = ref_cast_non_const(Champ_Fonc_base, critere_Q.valeur());
      if ((ch.temps() == temps_init) && (la_vitesse->mon_equation_non_nul()))
        ch.mettre_a_jour(la_vitesse->temps());
      return champs_compris_.get_champ(nom);
    }

  if (nom == "y_plus")
    {
      if (!y_plus)
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      Champ_Fonc_base& ch = ref_cast_non_const(Champ_Fonc_base, y_plus.valeur());
      if (((ch.temps() != la_vitesse->temps()) || (ch.temps() == temps_init)) && (la_vitesse->mon_equation_non_nul()))
        ch.mettre_a_jour(la_vitesse->temps());
      return champs_compris_.get_champ(nom);
    }

  if (nom == "reynolds_maille")
    {
      if (!Reynolds_maille)
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      Champ_Fonc_base& ch = ref_cast_non_const(Champ_Fonc_base, Reynolds_maille.valeur());
      if ((ch.temps() == temps_init) && (la_vitesse->mon_equation_non_nul()))
        ch.mettre_a_jour(la_vitesse->temps());
      return champs_compris_.get_champ(nom);
    }

  if (nom == "courant_maille")
    {
      if (!Courant_maille)
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      Champ_Fonc_base& ch = ref_cast_non_const(Champ_Fonc_base, Courant_maille.valeur());
      if (((ch.temps() != la_vitesse->temps()) || (ch.temps() == temps_init)) && (la_vitesse->mon_equation_non_nul()))
        ch.mettre_a_jour(la_vitesse->temps());
      return champs_compris_.get_champ(nom);
    }

  if (nom == "taux_cisaillement")
    {
      if (!Taux_cisaillement)
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      Champ_Fonc_base& ch = ref_cast_non_const(Champ_Fonc_base, Taux_cisaillement.valeur());
      if ((ch.temps() == temps_init) && (la_vitesse->mon_equation_non_nul()))
        ch.mettre_a_jour(la_vitesse->temps());
      return champs_compris_.get_champ(nom);
    }

  if (nom == "gradient_vitesse")
    {
      if (!grad_u)
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      Champ_Fonc_base& ch = ref_cast_non_const(Champ_Fonc_base, grad_u.valeur());
      if ((ch.temps() == temps_init) && (la_vitesse->mon_equation_non_nul()))
        ch.mettre_a_jour(la_vitesse->temps());
      return champs_compris_.get_champ(nom);
    }

  if (nom == "pression_hydrostatique")
    {
      if (!pression_hydrostatique_)
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      Champ_Fonc_base& ch = ref_cast_non_const(Champ_Fonc_base, pression_hydrostatique_.valeur());
      if (((ch.temps() != la_vitesse->temps()) || (ch.temps() == temps_init)) && (la_vitesse->mon_equation_non_nul()))
        {
          calculer_pression_hydrostatique(ch);
          ch.mettre_a_jour(la_vitesse->temps());
        }
      return champs_compris_.get_champ(nom);
    }

  OBS_PTR(Champ_base) ref_champ;

  if (Equation_base::has_champ(nom, ref_champ))
    return ref_champ;

  if (le_traitement_particulier)
    if (le_traitement_particulier->has_champ(nom, ref_champ))
      return ref_champ;

  throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));
}

void Navier_Stokes_std::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  Equation_base::get_noms_champs_postraitables(nom, opt);

  if (le_traitement_particulier)
    le_traitement_particulier->get_noms_champs_postraitables(nom, opt);

  Noms noms_compris = champs_compris_.liste_noms_compris();
  noms_compris.add("vorticite");
  noms_compris.add("critere_Q");
  noms_compris.add("y_plus");
  noms_compris.add("reynolds_maille");
  noms_compris.add("courant_maille");
  noms_compris.add("taux_cisaillement");
  noms_compris.add("pression_hydrostatique");
  noms_compris.add("gradient_vitesse");

  if (opt == DESCRIPTION)
    Cerr << " Navier_Stokes_std : " << noms_compris << finl;
  else
    nom.add(noms_compris);
}

/*! @brief Prints some information to an output stream: - maximum of div U
 *
 *        - convective term
 *        - diffusive term
 *        - divergence
 *        - gradient
 *
 * @param (Sortie& os) output stream
 * @return (int) always returns 1
 */
int Navier_Stokes_std::impr(Sortie& os) const
{
  // Display volumetric balance if not in quasi-compressible mode or Front Tracking
  if (!probleme().is_dilatable() && probleme().que_suis_je()!="Probleme_FT_Disc_gen")
    {
      double LocalFlowRateError=mp_max_abs_vect(divergence_U->valeurs());
      os << finl;
      os << "Cell balance flow rate control for the problem " << probleme().le_nom() << " : " << finl;
      os << "Absolute value : " << LocalFlowRateError << " m"<<dimension+bidim_axi<<"/s" << finl;
      os << "Relative value : " << LocalFlowRateRelativeError() << finl; // max|sum(u.ndS)i/(vol(i)/dt)|=max|div(U)i/dt|
      // Calculation as OpenFOAM: http://foam.sourceforge.net/docs/cpp/a04190_source.html
      // It is relative errors (normalized by the volume/dt)
      double dt = schema_temps().pas_de_temps();
      double local = LocalFlowRateError / ( probleme().domaine().volume_total() / dt );
      double global = mp_somme_vect(divergence_U->valeurs()) / ( probleme().domaine().volume_total() / dt );
      cumulative_ += global;
      os << "time step continuity errors : sum local = " << local << ", global = " << global << ", cumulative = " << cumulative_ << finl;
      // New in 1.6.1: stop if mass balance is bad and threshold < 1.e20
      if (local>0.01 && sub_type(solv_iteratif,solveur_pression_.valeur()))
        {
          if (ref_cast(solv_iteratif,solveur_pression_.valeur()).get_seuil()<1e10)
            {
              Cerr << "The mass balance is too bad (relative value>1%)." << finl;
              Cerr << "Please check and lower the convergence value of the pressure solver." << finl;
              exit();
            }
        }
#ifndef TRUST_USE_GPU
      // Since 1.6.6, warning to use PETSc Cholesky instead of an iterative method for pressure solver
      int nw=100;
      if (solveur_pression_->solveur_direct()==0 && le_schema_en_temps->nb_pas_dt()<nw && Process::nproc()<256 && la_pression->valeurs().size_array()<40000)
        {
          Cerr << finl << "********************** Advice (printed only on the first " << nw << " time steps) *********************" << finl;
          Cerr << "You should use PETSc Cholesky solver instead of an iterative method for the pressure solver." << finl;
          Cerr << "For the caracteristics of your problem, it will be faster and give a better mass flow balance." << finl;
          Cerr << "**********************************************************************************************" << finl << finl;
        }
#endif
    }

  if ((seuil_divU < 1.) && (sub_type(solv_iteratif,solveur_pression_.valeur())))
    {
      const solv_iteratif& solv_iter=ref_cast(solv_iteratif,solveur_pression_.valeur());
      os << " convergence threshold of the iterative solver  : " << solv_iter.get_seuil() << finl;
    }
  Equation_base::impr(os);
  divergence.impr(os);
  gradient.impr(os);
  return 1;
}


/*! @brief Returns the name of the application domain: "Hydraulique".
 *
 * @return the name representing the application domain
 */
const Motcle& Navier_Stokes_std::domaine_application() const
{
  static Motcle domaine = "Hydraulique";
  return domaine;
}

static void construire_matrice_implicite(Operateur_base& op,
                                         const DoubleTab& valeurs_inconnue,
                                         const Solveur_Masse_base& solv_masse,
                                         const double dt)
{
  Matrice& mat = op.set_matrice();
  if(!mat)
    mat.typer("Matrice_Morse");

  if(op.get_decal_temps()==1)
    {
      Matrice_Morse& matrice = ref_cast(Matrice_Morse, mat.valeur());
      op.dimensionner(matrice);
      op.contribuer_a_avec(valeurs_inconnue, matrice);
      solv_masse.ajouter_masse(dt, matrice);
      matrice *= dt;

      // If the solver is Cholesky or GCP, a matrix of type
      // Matrice_Morse_Sym is expected. Convert the matrix type:
      const Nom& type_solveur = op.get_solveur()->que_suis_je();
      if(type_solveur == "Solv_Cholesky" || type_solveur == "Solv_GCP")
        {
          Matrice_Morse_Sym new_mat(matrice);
          new_mat.set_est_definie(1);
          // mat is destroyed and then rebuilt:
          mat = new_mat;
          // Reinitialize the solver (recompute preconditioners, factorization, etc...)
          //ref_cast_non_const(SolveurSys_base,op.get_solveur().valeur()).reinit();
          op.set_solveur()->reinit();
        }
    }
}

/* In PolyMAC_HFV, the gradient contributes to the N-S equation matrix */
void Navier_Stokes_std::dimensionner_matrice_sans_mem(Matrice_Morse& matrice)
{
  Equation_base::dimensionner_matrice_sans_mem(matrice);
  if (gradient->has_interface_blocs())
    gradient->dimensionner_blocs({{ "vitesse", &matrice }});
}

int Navier_Stokes_std::has_interface_blocs() const
{
  return Equation_base::has_interface_blocs() && gradient->has_interface_blocs();
}

/* le gradient passe en dernier */
void Navier_Stokes_std::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  Equation_base::dimensionner_blocs(matrices, semi_impl);
  gradient->dimensionner_blocs(matrices, semi_impl);
}

void Navier_Stokes_std::assembler_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  Equation_base::assembler_blocs(matrices, secmem, semi_impl);
  gradient->ajouter_blocs(matrices, secmem, semi_impl);
}

DoubleTab& Navier_Stokes_std::derivee_en_temps_inco(DoubleTab& derivee)
{
  reassembler_pression_si_necessaire();
  // Computation of the time derivative:
  if(!implicite_)
    {
      // Explicit computation, use the standard time derivative:
      return Equation_base::derivee_en_temps_inco(derivee);
    }
  else
    {
      // Implicit computation of one or more operators (rarely used)
      // Dataset syntax: operateur { implicite solveur cholesky|gcp ... }
      derivee = 0;
      for(int i=0; i<nombre_d_operateurs(); i++)
        operateur(i).ajouter(derivee);

      les_sources.ajouter(derivee);
      derivee.echange_espace_virtuel();

      const double dt=schema_temps().pas_de_temps();
      static double dt_old=dt;

      for(int i=0; i<nombre_d_operateurs(); i++)
        {
          Operateur_base& op=operateur(i).l_op_base();
          // If matrix not build or matrix time dependant:
          if(!op.get_matrice() || !sys_invariant_)
            construire_matrice_implicite(op, inconnue().valeurs(), solv_masse(), dt);

          if(op.get_decal_temps()==1)
            {
              if(sys_invariant_ && dt!=dt_old)
                {
                  // The matrix does not change but the time step changes.
                  // The matrix is written A =
                  // Simplified update of the matrix
                  Matrice_Morse& matrice=ref_cast(Matrice_Morse, op.set_matrice().valeur());
                  matrice/=dt_old;
                  solv_masse().ajouter_masse(-dt_old, op.set_matrice().valeur());
                  solv_masse().ajouter_masse(dt, op.set_matrice().valeur());
                  matrice*=dt;
                  ref_cast_non_const(SolveurSys_base,op.get_solveur().valeur()).reinit();
                }
              Matrice_Morse& matrice=ref_cast(Matrice_Morse, op.set_matrice().valeur());
              if(implicite_==1)
                {
                  // A single implicit operator.
                  DoubleTrav secmem(derivee);
                  secmem=derivee;
                  DoubleTrav incre_pre(la_pression->valeurs());
                  gradient.calculer(la_pression->valeurs(),gradient_P->valeurs());
                  secmem-=gradient_P->valeurs();
                  uzawa(secmem, matrice,op.set_solveur(),derivee, incre_pre);

                  la_pression->valeurs()+=incre_pre;
                  gradient.calculer(la_pression->valeurs(),gradient_P->valeurs());
                }
              else
                {
                  // plusieurs operateurs implicites ...
                  Cerr << "To be developped ... " << finl;
                  exit();
                }
            }
        }
      dt_old=dt;
      return derivee;
    }
}

void Navier_Stokes_std::uzawa(const DoubleTab& secmem, const Matrice_Base& A, SolveurSys& solveur, DoubleTab& U, DoubleTab& P)
{
  // A U + Bt P = secmem
  // B U        = G
  // Start from the current pressure and
  // secmem = inertia + conv + sources + diff BCs
  // Start from P0 and U0 satisfying BCs and BU0=G

  // AU + Bt Cp = secmem
  // BU         = 0

  // Write a CG on B(A-1)Bt Cp = B(A-1)(secmem)

  DoubleTrav Cu(U);
  DoubleTrav grad(U);
  DoubleTrav grad0(U);
  DoubleTrav resu(P);
  DoubleTrav residu(P);
  DoubleTrav Cp(P);
  double dold,dnew,alfa;
  double seuil=seuil_uzawa;

  //Cu = A(-1) secmem
  Cerr << "Begining Uzawa, secmem norm value : " << mp_norme_vect(secmem) << finl;
  P=0.;
  gradient->multvect(P, grad0);
  solveur.nommer("uzawa_solver");
  solveur.resoudre_systeme(A, secmem, U);
  solv_masse().corriger_solution(U,Cu); // for Dirichlet boundary conditions!

  // residu=BCu
  divergence->multvect(U, resu);

  residu.copy(resu);
  residu*=-1.;

  // Cp = -residu;
  Cp = residu;
  Cp*=-1;
  Cp.echange_espace_virtuel();

  // Square of the norm
  dold = mp_norme_vect(residu);
  dold = dold * dold;
  dnew = dold;

  double s=0;
  int niter=0;
  int nmax=Cp.size();
  Cerr << "Uzawa, initial residue : " << dnew << finl;
  //     seuil=std::max(seuil, dnew*1.e-12);
  while ( ( dnew > seuil ) && (niter++ < nmax) )
    {
      gradient->multvect(Cp, grad);
      grad-=grad0;
      grad*=-1;
      solveur.resoudre_systeme(A, grad, Cu);
      solv_masse().corriger_solution(Cu,U); // for Dirichlet boundary conditions!
      divergence->multvect(Cu, resu);
      resu*=-1;

      s = mp_prodscal(resu, Cp);
      alfa = dold/(s);
      P.ajoute(alfa,Cp);
      residu.ajoute(alfa,resu);
      U.ajoute(alfa,Cu);
      dnew = mp_norme_vect(residu);
      dnew = dnew * dnew;
      assert(dnew >= 0);
      Cp *= (dnew/dold);
      Cp -= residu;
      dold = dnew;
      //       Cerr << "Uzawa, Residu apres "
      //            << niter << " iterations = " << dnew << finl;
    }

  if(dnew > seuil)
    {
      Cerr << "######## Uzawa, No convergence after : " << niter << " iterations\n";
      Cerr << "######## Uzawa,  Residue : "<< dnew << "\n";
      exit();
    }

  else if ((je_suis_maitre()))
    {
      Cerr << finl << "Uzawa, convergence reached after " << niter << " iterations" << finl;
    }
  {
    // On verifie :

    DoubleTrav R(resu);
    divergence->multvect(U, R);

    gradient->multvect(P, grad);
    grad-=grad0;
    DoubleTrav F(secmem);
    DoubleTrav UU(U);
    F=secmem;
    F-=grad;

    solveur.resoudre_systeme(A, F, UU);
    UU-=U;

    // Optimization: combine 2 mp_norme_vect into 1 collective call
    double R_carre = local_carre_norme_vect(R);
    double UU_carre = local_carre_norme_vect(UU);
    Process::mp_sum_for_each(R_carre, UU_carre);
    Cerr << "Ending Uzawa : mass residue : " << sqrt(R_carre) <<finl;
    Cerr << "Ending Uzawa : Qdm residue : " << sqrt(UU_carre)<<finl;
  }
}

void Navier_Stokes_std::sauver() const
{
  if (le_traitement_particulier)
    le_traitement_particulier->sauver_stat();
}

const Champ_Inc_base& Navier_Stokes_std::rho_la_vitesse() const
{
  Cerr<<" Navier_Stokes_std::rho_la_vitesse() must be overloaded "<<finl;
  assert(0);
  exit();
  throw;
}

void Navier_Stokes_std::update_y_plus(const DoubleTab& tab)
{
  if (!y_plus) Process::exit(que_suis_je() + " : y_plus must be initialised so it can be updated") ;
  DoubleTab& tab_y_p = y_plus->valeurs();
  if (tab.nb_dim()==2)
    for (int i = 0 ; i < tab_y_p.dimension_tot(0) ; i++)
      for (int n = 0 ; n < tab_y_p.dimension_tot(1) ; n++) tab_y_p(i,n) = tab(i,n);
  if (tab.nb_dim()==3)
    for (int i = 0 ; i < tab_y_p.dimension_tot(0) ; i++)
      for (int n = 0 ; n < tab_y_p.dimension_tot(1) ; n++) tab_y_p(i,n) = tab(i,0,n);
}

void Navier_Stokes_std::setPressureTimeN()
{
  if(probleme().domaine().getCouplingMethod()) //implicit IFS coupling
    // Equivalent of setting present = past for velocity done in for eg. Schema_Euler_Implicite::test_stationnaire.
    // Ensures sub-iterations start with the correct pressure for implicit coupling
    pression().valeurs()= P_n;
}

bool Navier_Stokes_std::getCouplingInfoForFiltering() const
{
  return probleme().domaine().getCouplingMethod();
}

void Navier_Stokes_std::updateFluidForce(DoubleTab& velocity)
{
  // in case of implicit coupling with a structural code: update the fluxes (used for computing the fluid force) during implicit sub-iterations

  if(probleme().domaine().getCouplingMethod()) //implicit IFS coupling case
    {
      Cout<<" Implicit coupling: Navier_Stokes_std_ALE::updateFluidForce "<<finl;
      //update diffusion operator
      DoubleTab field_value = velocity;
      field_value = 0.;
      operateur_diff().ajouter(velocity, field_value);


      //update gradient operator
      pression().mettre_a_jour(schema_temps().temps_courant());
      calculer_la_pression_en_pa();
      gradient->calculer_flux_bords();
    }
}
