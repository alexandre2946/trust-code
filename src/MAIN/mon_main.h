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

#ifndef mon_main_included
#define mon_main_included


#include <Process.h>
#include <Interprete_bloc.h>
#include <Comm_Group.h>
#include <list>

/*! @brief Class created and run by main() and during a TRUST execution through Python.
 *
 * Usage:
 *    - create a mon_main instance
 *    - initialise the parallel layer
 *    - call dowork(nom_du_cas) (read and interpret the dataset)
 *    At this point Python objects can be manipulated
 *    (see Interprete_bloc::objet_global(nom))
 *    - destroy the mon_main instance
 *    The main interpreter keeps its objects until the mon_main instance is destroyed.
 *
 */
class mon_main
{
public:
  mon_main(int verbose_level=9, bool journal_master=false, Nom log_directory = "",  bool apply_verification=true, bool disable_stop=false);
  ~mon_main();

  void init_parallel(const int argc, char **argv, bool with_mpi, bool check_enabled=false, bool with_petsc=true);
  void finalize();
  void dowork(const Nom& nom_du_cas);

private:
  int verbose_level_;
  bool journal_master_;

  Nom log_directory_;
  bool trio_began_mpi_;
  bool disable_stop_;
  OWN_PTR(Comm_Group) groupe_trio_;
  OWN_PTR(Comm_Group) node_group_;
  OWN_PTR(Comm_Group) node_master_;
  Interprete_bloc interprete_principal_;
};
extern bool error_handlers;

#endif
