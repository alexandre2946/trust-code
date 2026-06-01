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

#include <TRUSTTab.h>
#include <Roue.h>

// NOTE: OF 19/6/98
// No Doxygen comments for Roue_Ptr — it is not a class.

/*! @brief
 *
 */
Roue_ptr::Roue_ptr() : ptr(new(Roue))
{
}


/*! @brief Copy constructor.
 *
 * Copies the entire wheel. Useful for copy constructors of fields.
 * (This constructor is called when doing OWN_PTR(Champ_base) = another field;
 *  the previous version copied the reference, leading to a crash at destruction.)
 *
 * @param x The Roue_ptr to copy.
 */
Roue_ptr::Roue_ptr(const Roue_ptr& x)
{
  ptr=new Roue(x.valeur());
}

/*! @brief Copy assignment operator.
 *
 * Copies the entire wheel.
 *
 * @param x The Roue_ptr to copy.
 * @return Reference to the modified Roue_ptr.
 */
Roue_ptr& Roue_ptr::operator=(const Roue_ptr& x)
{
  Cerr << "We pass through Roue_ptr=Roue_ptr" << finl;
  assert(0);
  ptr=new Roue(x.valeur());
  return *this;
}

/*! @brief
 *
 */
Roue_ptr::Roue_ptr(Roue& x):ptr(&x)
{
}

/*! @brief
 *
 */
Roue_ptr::~Roue_ptr()
{
  if(ptr)
    {
      delete ptr;
      ptr=0;
    }
}

/*! @brief
 *
 */
Roue& Roue_ptr::operator[](int i)
{
  return ptr->futur(i);
}


/*! @brief
 *
 */
const Roue& Roue_ptr::operator[](int i) const
{
  return ptr->futur(i);
}


/*! @brief Default constructor. Builds a wheel with a single slot at time t=0.
 *
 */
Roue::Roue() : temps_(0), nb_cases_(1), valeurs_(), passe_(this), futur_(this)
{
}


/*! @brief Copy constructor. Copies the time, values, and values of future wheels.
 *
 * @param roue The Roue to copy.
 */
Roue::Roue(const Roue& roue) : temps_(roue.temps_), nb_cases_(1), valeurs_(roue.valeurs_), passe_(this), futur_(this)
{
  fixer_nb_cases(roue.nb_cases_);
  for(int i=1; i<nb_cases_; i++)
    {
      futur(i).valeurs_ = roue.futur(i).valeurs_;
      // abort();
    }
}


/*! @brief Destructor.
 *
 */
Roue::~Roue()
{
  supprimer_case(nb_cases_);
  assert(nb_cases_==1);
  assert(futur_.ptr==this);
  assert(passe_.ptr==this);
  assert(valeurs_.size()==0);
  futur_.ptr=passe_.ptr=0;
}

/*! @brief Returns the Roue corresponding to the i-th future slot.
 *
 * @param i Index of the future slot.
 * @return Reference to the i-th future Roue.
 */
const Roue& Roue::futur(int i) const
{
  assert(i>=0);
  assert(i<=nb_cases_);
  if(i==0)
    return *this;
  else
    return futur_->futur(--i);
}


/*! @brief Returns the Roue corresponding to the i-th future slot.
 *
 * @param i Index of the future slot.
 * @return Reference to the i-th future Roue.
 */
Roue& Roue::futur(int i)
{
  assert(i>=0);
  assert(i<=nb_cases_);
  if(i==0)
    return *this;
  else
    return futur_->futur(--i);
}


/*! @brief Returns the Roue corresponding to the i-th past slot.
 *
 * @param i Index of the past slot.
 * @return Reference to the i-th past Roue.
 */
const Roue& Roue::passe(int i) const
{
  assert(i>=0);
  assert(i<=nb_cases_);
  if(i==0)
    return *this;
  else
    return passe_->passe(--i);
}


/*! @brief Returns the Roue corresponding to the i-th past slot.
 *
 * @param i Index of the past slot.
 * @return Reference to the i-th past Roue.
 */
Roue& Roue::passe(int i)
{
  assert(i>=0);
  assert(i<=nb_cases_);
  if(i==0)
    return*this;
  else
    return passe_->passe(--i);
}

//    Roue::Roue(const Roue& roue) : passe_(roue.passe_),
//    futur_(roue.passe_), temps_(roue.temps_)
//    {
//       valeurs_.ref(roue.valeurs_);
//    }


/*! @brief Resizes (1D) the values of the Roue and its future wheels.
 *
 * @param nb_val Number of values.
 */
void Roue::dimensionner(int nb_val)
{
  if (valeurs_.size() != nb_val )
    {
      valeurs_.resize(nb_val);
      for(int j=1; j<nb_cases_; j++)
        futur(j).valeurs_.resize(nb_val);
    }
}


/*! @brief Resizes (2D) the values of the Roue and its future wheels.
 *
 * @param nb_val First dimension of the values.
 * @param nb_comp Second dimension (number of components) of the values.
 */
void Roue::dimensionner(int nb_val, int nb_comp)
{
  if ((valeurs_.size() != nb_val*nb_comp )||(valeurs_.nb_dim()!=2) ||  (nb_comp!=valeurs_.dimension(1)))
    {
      valeurs_.resize(nb_val, nb_comp);
      for(int j=1; j<nb_cases_; j++)
        futur(j).valeurs_.resize(nb_val, nb_comp);
    }
}


/*! @brief Changes the number of slots in the Roue.
 *
 * @param nb_case The new number of slots.
 * @return The new number of slots.
 */
int Roue::fixer_nb_cases(int nb_case)
{
  if(nb_case==nb_cases_)
    return nb_cases_;
  int nb_case_a_ajouter=nb_case-nb_cases_;
  if(nb_case_a_ajouter>0)
    ajouter_case(nb_case_a_ajouter);
  else
    supprimer_case(-nb_case_a_ajouter);
  return nb_cases_;
}


/*! @brief Adds n slots to the Roue.
 *
 * @param n Number of slots to add.
 */
void Roue::ajouter_case(int n)
{
  for (int n_to_add = 0; n_to_add < n; n_to_add++)
    {
      Roue * old_futur = futur_.ptr;
      Roue * new_futur = new Roue;
      this     ->futur_.ptr = new_futur;
      old_futur->passe_.ptr = new_futur;
      new_futur->futur_.ptr = old_futur;
      new_futur->passe_.ptr = this;
      // Update nb_cases_ for all slots:
      {
        const int new_nb_cases = nb_cases_ + 1;
        Roue * ptr_roue = this;
        for(int i = 0; i < new_nb_cases; i++, ptr_roue = ptr_roue->futur_.ptr)
          ptr_roue->nb_cases_ = new_nb_cases;
      }
    }
}

/*! @brief Removes n slots from the Roue.
 *
 * @param n Number of slots to remove.
 */
void Roue::supprimer_case(int n)
{
  for (int n_to_kill = 0; n_to_kill < n; n_to_kill++)
    {
      if (nb_cases_ == 1)
        {
          assert(futur_.ptr == this);
          assert(passe_.ptr == this);
          valeurs_.reset();
          temps_=0;
        }
      else
        {
          // Save the future-of-future:
          Roue * new_futur = futur_->futur_.ptr;
          // Destroy the future slot:
          {
            Roue * to_kill = futur_.ptr;
            to_kill->nb_cases_ = 1;
            to_kill->futur_.ptr = to_kill->passe_.ptr = to_kill;
            delete to_kill;
          }
          // Close the linked list:
          new_futur->passe_.ptr = this;
          futur_.ptr = new_futur;

          const int new_nb_cases = nb_cases_ - 1;
          if(new_nb_cases == 1)
            passe_.ptr=this;

          // Update nb_cases_ for all slots:
          {
            Roue * ptr_roue = this;
            for(int i = 0; i < new_nb_cases; i++, ptr_roue = ptr_roue->futur_.ptr)
              ptr_roue->nb_cases_ = new_nb_cases;
          }
        }
    }
}


/*! @brief Assignment operator for a Roue. Copies the time, values, number of slots, and future slot values.
 *
 * @param roue The Roue to copy.
 * @return Reference to the modified Roue.
 */
Roue& Roue::operator=(const Roue& roue)
{
  temps_ = roue.temps_;
  valeurs_ = roue.valeurs_;
  fixer_nb_cases(roue.nb_cases_);
  for(int i=1; i<nb_cases_; i++)
    {
      futur(i).valeurs_ = roue.futur(i).valeurs_;
      // abort();
    }
  return *this;
}

