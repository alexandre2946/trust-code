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

#include <Solv_tools.h>
#include <Objet_U.h>

void Solv_tools::construit_renum(const DoubleVect& b)
{
  // Initialize the items_to_keep_ array if not already done
  nb_items_to_keep_ = b.get_md_vector()->get_sequential_items_flags(items_to_keep_, b.line_size());

  // Compute important value:
  secmem_sz_ = b.size_totale();
  nb_rows_ = nb_items_to_keep_;
  nb_rows_tot_ = Process::mp_sum(nb_rows_);
  decalage_local_global_ = Process::mppartial_sum(nb_rows_);
  //Journal()<<"nb_rows_=" << nb_rows_ << " nb_rows_tot_=" << nb_rows_tot_ << " decalage_local_global_=" << decalage_local_global_ << finl;

  /**********************/
  /* Build renum_ array */
  /**********************/
  //if (MatricePetsc_==nullptr)
  {
    const MD_Vector& md = b.get_md_vector();
    renum_.reset();
    renum_.resize(0, b.line_size());
    MD_Vector_tools::creer_tableau_distribue(md, renum_, RESIZE_OPTIONS::NOCOPY_NOINIT);
  }
  int cpt=0;
  int size=items_to_keep_.size_array();
  renum_ = INT_MAX; //to crash if the MD_Vector is inconsistent
  ArrOfTID& renum_array = renum_;  // array viewed as linear
  for(int i=0; i<size; i++)
    if(items_to_keep_[i])
      {
        renum_array[i]=cpt+decalage_local_global_;
        cpt++;
      }

  renum_.echange_espace_virtuel();
  // Build index_
  index_.resize(size);
  int index = 0;
  // ToDo OpenMP: factorize with ix since index_=ix-decalage_local_global_
  for (int i=0; i<size; i++)
    {
      if (items_to_keep_[i])
        {
          index_[i] = index;
          index++;
        }
      else
        index_[i] = -1;
    }
  // Build ix
  size=b.size_array();
  auto colonne_globale=decalage_local_global_;
  ix.resize(size);
  for (int i=0; i<size; i++)
    if (items_to_keep_[i])
      {
        ix[i] = colonne_globale;
        colonne_globale++;
      }
    else
      ix[i] = -1;
}

