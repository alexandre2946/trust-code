
#include <Pb_Euler.h>
#include <Discretisation_base.h>
#include <Milieu_composite_Euler.h>
#include <Loi_Fermeture_base.h>
#include <Interprete_bloc.h>
#include <Domaine.h>
#include <EChaine.h>
#include <Debog.h>
#include <SETS.h>

//#include <Fluide_Incompressible.h>


Implemente_instanciable(Pb_Euler, "Pb_Euler", Pb_Fluide_base);

Sortie& Pb_Euler::printOn(Sortie& os) const
{

  return Pb_Fluide_base::printOn(os);
}

Entree& Pb_Euler::readOn(Entree& is)
{
  if (!discretisation().is_coloc())
    {
      Cerr << "Error: Problem of type " << que_suis_je() << " is only available for Coloc discretizations\n";
      Process::exit();
    }
  return Pb_Fluide_base::readOn(is);
}

void Pb_Euler::typer_lire_milieu(Entree& is)
{
  le_milieu_.resize(1);
  is >> le_milieu_[0];
  if (!sub_type(Milieu_composite_Euler, le_milieu_[0].valeur()))
    {
      Cerr << "Error: Fluid of type " << le_milieu_[0]->le_type() << " is not compatible with " << que_suis_je() << " problem which accepts only Milieu_composite_Euler medium" << finl;
      Cerr << "Check your datafile!" << finl;
      Process::exit();
    }
  noms_phases_ = ref_cast(Milieu_composite_Euler,le_milieu_[0].valeur()).noms_phases();
  associer_milieu_base(le_milieu_[0].valeur());
  Probleme_base::discretiser_equations();
  for (int i = 0; i < nombre_d_equations(); i++)
    equation(i).associer_milieu_equation();
  equation(0).milieu().discretiser((*this), la_discretisation_.valeur());
}

Entree& Pb_Euler::lire_equations(Entree& is, Motcle& mot)
{
  bool already_read {true};
  is >> mot;
  if (mot == "correlations" || mot == "models")
    lire_correlations(is), already_read = false;

  // typer_lire_correlation_hem(); // enforce an interfacial flux correlation with constant coefficient if HEM

  Cerr << "Reading of the equations" << finl;
  for (int i = 0; i < nombre_d_equations(); i++, already_read = false)
    {
      if (!already_read)
        is >> mot;
      is >> getset_equation_by_name(mot);
    }

  //correction des donnee lu depuis .data
  eq_masse_.init_alpha_rho();
  eq_qdm_.init_alpha_rho_u();
  eq_energie_.init_energie_tot();

  read_optional_equations(is, mot);
  return is;
}

const Equation_base& Pb_Euler::equation(int i) const
{
  if (i == 0)
    return eq_qdm_;
  else if (i == 1)
    return eq_masse_;
  else if (i == 2)
    return eq_energie_;
  else if (i == 3)
    return eq_fraction_;
  else
    {
      Cerr << "Pb_Euler::equation() : Wrong equation number" << i << "!" << finl;
      Process::exit();
    }
  return eq_qdm_;
}

Equation_base& Pb_Euler::equation(int i)
{
  if (i == 0)
    return eq_qdm_;
  else if (i == 1)
    return eq_masse_;
  else if (i == 2)
    return eq_energie_;
  else if (i == 3)
    return eq_fraction_;
  else
    {
      Cerr << "Pb_Euler::equation() : Wrong equation number" << i << "!" << finl;
      Process::exit();
    }
  return eq_qdm_; //pour renvoyer quelque chose
}


void Pb_Euler::associer_milieu_base(const Milieu_base& mil)
{
  /* controler le type de milieu ici */
  equation_qdm().associer_milieu_base(mil);
  equation_energie().associer_milieu_base(mil);
  equation_masse().associer_milieu_base(mil);
  equation_fraction().associer_milieu_base(mil);
}


//int Pb_Euler::verifier()
//{
//  // const Domaine_Cl_dis_base& domaine_Cl_hydr = equation_qdm().domaine_Cl_dis();
//  // const Domaine_Cl_dis_base& domaine_Cl_th = equation_energie().domaine_Cl_dis();
//  return 1; // tester_compatibilite_hydr_thermique(domaine_Cl_hydr, domaine_Cl_th);
//}


void Pb_Euler::preparer_calcul()
{
  Pb_Fluide_base::preparer_calcul();
  const double temps = schema_temps().temps_courant();
  mettre_a_jour(temps);
}

void Pb_Euler::mettre_a_jour(double temps)
{
  Probleme_base::mettre_a_jour(temps);
  equation_qdm().mettre_a_jour_p_c();
}





