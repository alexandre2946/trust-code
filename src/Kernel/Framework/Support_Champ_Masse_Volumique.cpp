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

#include <Support_Champ_Masse_Volumique.h>
#include <Nom.h>

/*! @brief Constructor of the class.
 *
 * By default, an already coded derived class calls the constructor without argument. support_ok is set to zero
 *   and an error is produced if Associer_champ_masse_volumique is called.
 *   To indicate that the density field is supported by
 *   the derived class, call Declare_support_masse_volumique
 *
 */
Support_Champ_Masse_Volumique::Support_Champ_Masse_Volumique() :
  support_ok_(0)
{
}

/*! @brief Virtual destructor (to avoid warnings)
 *
 */
Support_Champ_Masse_Volumique::~Support_Champ_Masse_Volumique()
{
}

/*! @brief The constructor of a derived class that uses the density field must call this function with the value 1.
 *
 * If a client class
 *   (Navier Stokes for example) tries to associate the density field to a
 *   class that has not set ok=1, it stops: function not implemented.
 *
 */
void Support_Champ_Masse_Volumique::declare_support_masse_volumique(int ok)
{
  support_ok_ = ok;
}

/*! @brief Method to be called during problem preparation to ask the object to take into account the density field passed as parameter.
 *
 */
void Support_Champ_Masse_Volumique::associer_champ_masse_volumique(const Champ_base& ch)
{
  if (! support_ok_)
    {
      Cerr << "Error in Support_Champ_Masse_Volumique::Associer_champ_masse_volumique" << finl;
      Cerr << " Support from a density field was not coded for this object" << finl;
      assert(0);
      Process::exit();
    }
  ref_champ_rho_ = ch;
}

/*! @brief Cancels the reference to the density field.
 *
 */
void Support_Champ_Masse_Volumique::dissocier_champ_masse_volumique()
{
  OBS_PTR(Champ_base) ref_nulle;
  ref_champ_rho_ = ref_nulle;
}

/*! @brief Returns 1 if the density field has been associated, 0 otherwise.
 *
 */
int Support_Champ_Masse_Volumique::has_champ_masse_volumique() const
{
  int ref_non_nulle = bool(ref_champ_rho_);
  return ref_non_nulle;
}

/*! @brief Returns the density field
 *
 */
const Champ_base& Support_Champ_Masse_Volumique::get_champ_masse_volumique() const
{
  assert(ref_champ_rho_);
  const Champ_base& ch = ref_champ_rho_.valeur();
  return ch;
}

