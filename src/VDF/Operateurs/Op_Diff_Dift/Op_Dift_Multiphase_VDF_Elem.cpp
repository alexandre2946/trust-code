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

#include <Op_Dift_Multiphase_VDF_Elem.h>
#include <Op_Dift_Multiphase_VDF_Face.h>
#include <Pb_Multiphase.h>

Implemente_instanciable_sans_constructeur(Op_Dift_Multiphase_VDF_Elem,"Op_Diff_VDFTURBULENTE_P0_VDF|Op_Diff_VDFTURBULENT_P0_VDF",Op_Dift_VDF_Elem_base);

Sortie& Op_Dift_Multiphase_VDF_Elem::printOn(Sortie& s ) const { return s << que_suis_je() ; }
Entree& Op_Dift_Multiphase_VDF_Elem::readOn(Entree& is)
{
  //lecture de la correlation de diffusivite turbulente
  Correlation_base::typer_lire_correlation(corr_, equation().probleme(), "transport_turbulent", is);
  associer_corr_impl<Type_Operateur::Op_DIFT_MULTIPHASE_ELEM, Eval_Dift_Multiphase_VDF_Elem>(corr_);
  associer_proto(equation().probleme(), champs_compris_);
  ajout_champs_proto_elem();
  return is;
}

void Op_Dift_Multiphase_VDF_Elem::get_noms_champs_postraitables(Noms& nom,Option opt) const
{
  Op_Dift_VDF_Elem_base::get_noms_champs_postraitables(nom,opt);
  get_noms_champs_postraitables_proto(que_suis_je(), nom, opt);
}

void Op_Dift_Multiphase_VDF_Elem::creer_champ(const Motcle& motlu)
{
  Op_Dift_VDF_Elem_base::creer_champ(motlu);
  creer_champ_proto_elem(motlu);
}

void Op_Dift_Multiphase_VDF_Elem::completer()
{
  assert(corr_.non_nul());
  completer_Op_Dift_VDF_base();
  associer_pb<Eval_Dift_Multiphase_VDF_Elem>(equation().probleme());
  completer_proto_elem(*this);
  set_nut_impl<Type_Operateur::Op_DIFT_MULTIPHASE_ELEM, Eval_Dift_Multiphase_VDF_Elem>(nu_ou_lambda_turb_);
}

void Op_Dift_Multiphase_VDF_Elem::mettre_a_jour(double temps)
{
  Op_Dift_VDF_Elem_base::mettre_a_jour(temps);

  const Operateur_base& op_qdm = equation().probleme().equation(0).operateur(0).l_op_base();
  if (!sub_type(Op_Dift_Multiphase_VDF_Face, op_qdm))
    {
      Cerr << "Error in " << que_suis_je() << ": no turbulent momentum diffusion found!" << finl;
      Process::exit();
    }

  const Correlation_base& corr_visc_qdm = ref_cast(Op_Dift_Multiphase_VDF_Face, op_qdm).correlation();
  if (!sub_type(Viscosite_turbulente_base, corr_visc_qdm))
    {
      Cerr << "Error in " << que_suis_je() << ": no turbulent viscosity correlation found!" << finl;
      Process::exit();
    }

  // on calcule d_t_
  nu_ou_lambda_turb_ = 0.; // XXX : pour n'avoir pas la partie laminaire
  call_compute_diff_turb(ref_cast(Convection_Diffusion_std, equation()), ref_cast(Viscosite_turbulente_base, corr_visc_qdm));
  set_nut_impl<Type_Operateur::Op_DIFT_MULTIPHASE_ELEM, Eval_Dift_Multiphase_VDF_Elem>(nu_ou_lambda_turb_);
  mettre_a_jour_proto_elem(temps);
}

double Op_Dift_Multiphase_VDF_Elem::calculer_dt_stab() const
{
  double dt_stab, coef = -1.e10;
  const Domaine_VDF& domaine_VDF = iter_->domaine();
  const IntTab& elem_faces = domaine_VDF.elem_faces();
  const DoubleTab& lambda = alpha_() /* comme mu */, &diffu = diffusivite_pour_pas_de_temps().valeurs() /* comme nu */;
  const DoubleTab* alp = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()).equation_masse().inconnue().passe() : nullptr;
  const DoubleTab& rho = equation().milieu().masse_volumique().passe();
  const int cL = (lambda.dimension(0) == 1), cD = (diffu.dimension(0) == 1), cR = (rho.dimension(0) == 1), dim = Objet_U::dimension;

  double mu_turbulent, mu_physique, nu_physique, alfa;

  ArrOfInt numfa(2 * dim);
  for (int elem = 0; elem < domaine_VDF.nb_elem(); elem++)
    {
      double diflo = 0.;
      double deltax = 0.;
      for (int i = 0; i < 2 * dim; i++)
        numfa[i] = elem_faces(elem, i);

      for (int d = 0; d < dim; d++)
        {
          const double hd = domaine_VDF.dist_face(numfa[d], numfa[dim + d], d);
          deltax += 1. / (hd * hd);
        }

      // nu_ou_lambda_turb = alpha * nut * sigma
      for (int ncomp = 0; ncomp < nu_ou_lambda_turb_.line_size(); ncomp++)
        {
          //if (elem==0) cout << "ncomp "<< ncomp << " nu_ou_lambda_turb = alpha * nut * sigma "<< nu_ou_lambda_turb_(elem, ncomp) << endl;
          alfa = (alp ? (*alp)(elem, ncomp) : 1.0);
          mu_turbulent = rho(!cR * elem, ncomp) * nu_ou_lambda_turb_(elem, ncomp);
          //if (elem==0) cout << "ncomp "<< ncomp << " mu_turbulent avec alpha "<< mu_turbulent << endl;
          if (alfa != 0.0 ) mu_turbulent = mu_turbulent/alfa;
          //if (elem==0) cout << "ncomp "<< ncomp << " mu_turbulent sans alpha "<< mu_turbulent << endl;
          mu_physique = rho(!cR * elem, ncomp) * lambda(!cL * elem, ncomp);
          nu_physique = diffu(!cD * elem, ncomp);
          //if (elem==0) cout << "ncomp "<< ncomp << " mu_physique "<< mu_physique << endl;
          //if (elem==0) cout << "ncomp "<< ncomp << " nu_physique "<< nu_physique << endl;

          // le pas de temps de stab est alpha(nu+nu_t), on calcule a(mu+mu_t)*(nu/mu)=a(mu+mu_t)/rho=a(nu+nu_t) (avantage par rapport a la division par rho ca marche aussi pour alpha et lambda et en VEF
          diflo = deltax * (mu_physique + mu_turbulent) * (nu_physique / mu_physique);
          //if (elem==0) cout << "ncomp "<< ncomp << " dt "<< 0.5 / (diflo + DMINFLOAT) << endl;
          coef = std::max(coef, diflo);
        }
    }

  coef = Process::mp_max(coef);
  dt_stab = 0.5 / (coef + DMINFLOAT);
  //cout << "dt_stab "<< dt_stab << endl;

  return dt_stab;
}

