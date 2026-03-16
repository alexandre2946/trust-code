/****************************************************************************
* Copyright (c) 2022, CEA
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

#ifndef TRUST_Post_Loader_included
#define TRUST_Post_Loader_included

#include <medcoupling++.h>
#include <LataFilter.h>
#include <vector>

namespace MEDCoupling
{
  class MEDCouplingMesh;
  class DataArray;
  class MEDCouplingFieldDouble;
}

class TRUST_Post_Loader
{
public:
  TRUST_Post_Loader(const char*, bool print = true);
  ~TRUST_Post_Loader() { /* Do nothing */ }

  int getNumberOfTimeSteps();

  const char* getPostExtension();

  inline std::vector<double> getTimes()
  {
    std::vector<double> a;
    getTimes_(a);
    return a;
  }

  MEDCoupling::MEDCouplingMesh* getMesh(const char *varname, int timestate, int block = -1);
  MEDCoupling::MEDCouplingFieldDouble* getFieldDouble(const char *varname, int timestate, int block = -1);

  std::vector<std::string> getMeshNames();
  std::vector<std::string> getFieldNames();
  std::vector<std::string> getFieldNamesOnMesh(const std::string&);

private:

  void populateDatabaseMetaData_(int);
  void getTimes_(std::vector<double> &times);
  void register_fieldname_(const char *visit_name, const Field_UName&, int component);
  void register_meshname_(const char *visit_name, const char *latafilter_name);
  void get_field_info_from_visitname_(const char *varname, Field_UName&, int &component) const;
  MEDCoupling::DataArray* getVectorVar_(int, int, const char*);

  LataDB lata_db_; // Source database
  LataFilter filter_; // Data processor and cache
  Noms field_username_;
  Field_UNames field_uname_;

  Noms mesh_username_, mesh_latafilter_name_;
  Nom filename_;
  // For each name, which component is it in the source field:
  LataVector<int> field_component_;
  bool print_ = false;
};

#endif /* TRUST_Post_Loader_included */
