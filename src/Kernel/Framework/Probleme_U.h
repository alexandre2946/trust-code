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

#ifndef Probleme_U_included
#define Probleme_U_included

#include <ScalarRegister.h>
#include <TRUST_Ref.h>

class Champ_Generique_base;
class Field_base;
class Noms;

namespace ICoCo
{
class TrioField;
}

/*! @brief Probleme_U
 *
 *      This class adds TRUST-specific features to the general
 *      interface defined in the Problem class.
 *
 *      Derived from Probleme_U:
 *      * Probleme_base for all individual problems
 *      * Couplage_U for couplings of several Probleme_U objects
 *
 *      All these classes must satisfy the Problem API
 *      extended by the Probleme_U API.
 *
 *
 */
class Probleme_U : public Objet_U
{
  Declare_base(Probleme_U);
public :

  // Implementation of the Problem API

  virtual void initialize();
  virtual void terminate();
  virtual double presentTime() const;
  virtual double computeTimeStep(bool& stop) const;
  virtual bool initTimeStep(double dt);
  virtual bool solveTimeStep();
  virtual void validateTimeStep();
  virtual bool isStationary() const;
  virtual std::string newCompute() { return ""; };
  virtual void setStationary(bool);
  virtual void abortTimeStep();
  virtual void resetTime(double time);
  virtual bool iterateTimeStep(bool& converged);
  virtual void getInputFieldsNames(Noms& noms) const;
  virtual void getInputFieldTemplate(const Nom& name, ICoCo::TrioField& afield) const;
  virtual void setInputField(const Nom& name, const ICoCo::TrioField& afield);
  virtual void getOutputFieldsNames(Noms& noms) const;
  virtual void getOutputField(const Nom& nameField, ICoCo::TrioField& afield) const;
  virtual void setInputIntValue(const Nom& name, const int& val);
  virtual int getOutputIntValue(const Nom& name) const;
  virtual bool checkOutputIntEntry(const Nom& name) const;

  virtual void getOutputPointValues(const Nom& name,
                                    const std::vector<double>& x,
                                    const std::vector<double>& y,
                                    const std::vector<double>& z,
                                    std::vector<double>& vals, int compo)
  {
    Cerr << "Probleme_U::" << __func__ << " must be overloaded in " << que_suis_je() << finl;
    Process::exit();
  }

  double getOutputPointValues(const Nom& name, const double x, const double y, const double z, int compo);

  virtual void setInputDoubleValue(const Nom& name, const double val);

  virtual void setInputStringValue(const std::string& name, const std::string& val) { str_params_[name] = val; }
  virtual std::string getOutputStringValue(const std::string& name);

  virtual void post_initialize() { }

  // Complements to the Problem API

  virtual int postraiter(int force=1);
  virtual int limpr() const;
  virtual int lsauv() const;
  virtual void sauver() const;
  virtual bool updateGivenFields();

  virtual double futureTime() const;
  virtual bool run();
  virtual bool runUntil(double time);

  virtual OBS_PTR(Field_base) findInputField(const Nom& name) const;
  virtual OBS_PTR(Champ_Generique_base) findOutputField(const Nom& name) const;

  inline void nommer(const Nom& name) override { nom_ = name; }
  inline const Nom& le_nom() const override { return nom_; }

protected :
  Nom nom_;
  ScalarRegister reg_;
  std::map<std::string, std::string> str_params_;

};

#endif /* Probleme_U_included */
