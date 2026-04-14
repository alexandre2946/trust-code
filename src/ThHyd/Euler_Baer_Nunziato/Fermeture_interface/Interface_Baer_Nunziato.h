#ifndef Interface_Baer_Nunziato_included
#define Interface_Baer_Nunziato_included

#include <Interface_base.h>

class Interface_Baer_Nunziato : public Interface_base
{
  Declare_instanciable(Interface_Baer_Nunziato);
public:
//  virtual const DoubleTab vitesse() const;
//  virtual const DoubleTab pression() const;
//  virtual DoubleTab vitesse();
//  virtual DoubleTab pression();
  void sigma_(const SpanD T, const SpanD P, SpanD res, int ncomp = 1, int ind = 0) const override {Process::exit();};
  void sigma_h_(const SpanD H, const SpanD P, SpanD res, int ncomp = 1, int ind = 0) const override {Process::exit();};
  // void discretiser_sigma(const Nom& sig_nom, double temps) override {};
  void mettre_a_jour(double ) override {};
  void set_param(Param& param) override;
  int id_phase_vitesse_inter() const { return id_vitesse_interface_; }
  int id_phase_pression_inter() const { return id_pression_interface_; }
private:
  int id_vitesse_interface_, id_pression_interface_;
};

#endif /* Interface_Baer_nunziato_included */
