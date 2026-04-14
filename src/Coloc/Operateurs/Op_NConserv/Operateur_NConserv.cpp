
#include <Discretisation_base.h>
#include <Operateur_NConserv.h>
#include <stat_counters.h>
#include <Champ_base.h>

Implemente_instanciable(Operateur_NConserv,"Operateur_NConserv",OWN_PTR(Operateur_NConserv_base));

Sortie& Operateur_NConserv::printOn(Sortie& os) const
{
  return Operateur::ecrire(os);
}

Entree& Operateur_NConserv::readOn(Entree& is)
{
  Operateur::lire(is);
  return is;
}

void Operateur_NConserv::typer()
{
  if (Motcle(typ)==Motcle("negligeable"))
    {
      OWN_PTR(Operateur_NConserv_base)::typer("Op_NConserv_negligeable");
    }
  else
    {
      Equation_base& eqn=mon_equation.valeur();
      Nom nom_type=eqn.discretisation().get_name_of_type_for(que_suis_je(),typ,eqn);
      OWN_PTR(Operateur_NConserv_base)::typer(nom_type);
    }
  Cerr << valeur().que_suis_je() << finl;
}


DoubleTab& Operateur_NConserv::ajouter(const DoubleTab& donnee,
                                       DoubleTab& resu) const
{
  //statistiques().begin_count(convection_counter_);
  DoubleTab& tmp = valeur().ajouter(donnee, resu);
  //statistiques().end_count(convection_counter_);
  return tmp;
}


DoubleTab& Operateur_NConserv::calculer(const DoubleTab& donnee,
                                        DoubleTab& resu) const
{
  //statistiques().begin_count(convection_counter_);
  DoubleTab& tmp = valeur().calculer(donnee, resu);
  //statistiques().end_count(convection_counter_);
  return tmp;
}

