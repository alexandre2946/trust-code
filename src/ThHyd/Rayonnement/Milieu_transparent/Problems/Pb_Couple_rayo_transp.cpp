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

#include <Pb_Couple_rayo_transp.h>
#include <Pb_Fluide_base.h>

Implemente_instanciable(Pb_Couple_rayo_transp, "Pb_Couple_rayo_transp", Probleme_Couple);

Entree& Pb_Couple_rayo_transp::readOn(Entree& is) { return is; }

Sortie& Pb_Couple_rayo_transp::printOn(Sortie& os) const { return Probleme_Couple::printOn(os); }

void Pb_Couple_rayo_transp::initialize()
{
  int nb_ray_mod_associated = 0;
  for (int l = 0; l < nb_problemes(); l++)
    {
      Probleme_base& pb = ref_cast(Probleme_base, probleme(l));

      if (pb.has_mod_rayo_transp())
        {
          nb_ray_mod_associated++;

          if (nb_ray_mod_associated > 1)
            Process::exit("We can not treat at present several associated radiation models ... \n");

          assert (sub_type(Pb_Fluide_base, pb));
          ref_cast(Pb_Fluide_base, pb).assoscier_rayo_model_CL();
        }
    }

  Probleme_Couple::initialize();
}
