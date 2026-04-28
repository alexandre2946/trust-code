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

#ifndef Interpolation_IBM_power_law_tbl_proto_included
#define Interpolation_IBM_power_law_tbl_proto_included

/*! @brief : class Interpolation_IBM_power_law_tbl_proto
 *
 *  Pure C++ class to allow multiple inheritance in Interpolation_IBM_power_law_tbl
 *
 *
 */

class Interpolation_IBM_power_law_tbl_proto
{

public :
// Constantes
  inline double get_A_pwl(int a)
  {
    if (a == 1)
      return A_pwl_WJSP_;
    else
      return  A_pwl_;
  };

  inline double get_C_pwl_WJSP() { return  C_pwl_WJSP_; }
  inline double get_D_pwl_WJSP() { return  D_pwl_WJSP_; }


// Puissances
  inline double get_B_pwl() { return  B_pwl_; }
  inline double get_p_pwl_WJSP() { return  p_pwl_WJSP_; }

// Limites
  inline double get_y_c_p_pwl() { return y_c_p_pwl_; }
  inline double get_y_c1_p_pwl_WJSP() { return y_c1_p_pwl_WJSP_; }
  inline double get_y_c2_p_pwl_WJSP() { return y_c2_p_pwl_WJSP_; }

protected :
  double A_pwl_ = 8.3;
  double A_pwl_WJSP_ = 8.438565669851982;
  double C_pwl_WJSP_ = 20.197802756630782;
  double D_pwl_WJSP_ = -34.3779384724374;
  double B_pwl_ = 1./7.;
  double p_pwl_WJSP_ = 0.49570985985551774;
  double y_c_p_pwl_ = 11.81;
  double y_c1_p_pwl_WJSP_ = 6.549332667516647;
  double y_c2_p_pwl_WJSP_ = 54.75271424217823;
  friend class Source_PDF_base;
};

#endif /* Interpolation_IBM_power_law_tbl_proto_included */
