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

#include <LecFicDiffuse_JDD.h>
#include <Execute_parallel.h>
#include <Interprete_bloc.h>
#include <TRUST_Error.h>
#include <TRUSTArray.h>
#include <PE_Groups.h>
#include <Journal.h>
#include <Param.h>

Implemente_instanciable(Execute_parallel,"Execute_parallel",Interprete);
// XD execute_parallel interprete execute_parallel BRACE This keyword allows to run several computations in parallel on
// XD_CONT processors allocated to TRUST. The set of processors is split in N subsets and each subset will read and
// XD_CONT execute a different data file. Error messages usualy written to stderr and stdout are redirected to .log
// XD_CONT files (journaling must be activated).

Entree& Execute_parallel::readOn(Entree& is)
{
  Cerr << "Execute_parallel::readOn() not coded" << finl;
  exit();
  return is;
}

Sortie& Execute_parallel::printOn(Sortie& os) const
{
  Cerr << "Execute_parallel::printOn() not coded" << finl;
  exit();
  return os;
}

/*! @brief Creates a partition of the nproc processors for the computation in order to interpret N different data sets.
 *
 * The syntax of the data set is as follows:
 *   Execute_parallel {
 *     liste_cas N cas1 cas2 cas3 ...
 *     [ nb_procs N nproc1 nproc2 nproc3 ... ]
 *   }
 *   "cas1" is the name of the case (the file cas1.data is read from disk)
 *   nproc1 is the number of processors to use for that case
 *   By default, 1 processor is used for each case
 *   Cerr and Cout outputs are redirected to the journal of the master
 *   processor of each case.
 *
 */
Entree& Execute_parallel::interpreter(Entree& is)
{
  Cerr << "Execute_parallel::interpreter to run several cases:" << finl;

  Noms liste_cas;
  ArrOfInt nb_procs;

  bool disable_journal = false;

  Param param(que_suis_je());
  param.ajouter("liste_cas", &liste_cas, Param::REQUIRED); // XD_ADD_P listchaine
  // XD_CONT N datafile1 ... datafileN. datafileX the name of a TRUST data file without the .data extension.
  param.ajouter("nb_procs", &nb_procs); // XD_ADD_P listentier
  // XD_CONT nb_procs is the number of processors needed to run each data file. If not given, TRUST assumes that
  // XD_CONT computations are sequential.
  param.ajouter_flag("disable_journal", &disable_journal);
  param.lire_avec_accolades_depuis(is);
  // If nb_procs was not given, assume it is 1
  const int n_calculs = liste_cas.size();
  if (nb_procs.size_array() == 0)
    {
      Cerr << "Nb_procs not given, we assume that calculations are sequential."
           << finl;
      nb_procs.resize_array(n_calculs);
      nb_procs = 1;
    }
  // Verify array sizes
  if (nb_procs.size_array() != n_calculs)
    {
      Cerr << "Error : nb_procs array must have " << n_calculs
           << " values." << finl;
      barrier();
      exit();
    }
  if (n_calculs == 0)
    return is;

  // Verify the content:
  if (min_array(nb_procs) < 1)
    {
      Cerr << "Error : processor numbers must be >= 1" << finl;
      barrier();
      exit();
    }
  int count = 0;
  for (int i = 0; i < n_calculs; i++)
    count += nb_procs[i];
  if (count > nproc())
    {
      Cerr << "Error : computations require " << count << " processors." << finl;
      Cerr << "but only " << nproc() << " processors has been asked." << finl;
      barrier();
      exit();
    }

  // Create the N processor groups
  // (groups are destroyed when the VECT is destroyed)
  VECT(OWN_PTR(Comm_Group)) groupes(n_calculs);
  count = 0;
  Nom log_courant("");
  for (int i = 0; i < n_calculs; i++)
    {
      const int n = nb_procs[i];
      Nom log("");
      Nom log1(Objet_U::nom_du_cas());
      log1+="_";
      char s[20];
      snprintf(s, 20, "%05d", (int)count);
      log1+=s;
      log1+=".log";
      if (n==1)
        {
          log+=log1;
          log+=" file";
        }
      else
        {
          Nom log2(Objet_U::nom_du_cas());
          log2+="_";
          char s2[20];
          snprintf(s2, 20, "%05d", (int)(count+n-1));
          log2+=s2;
          log2+=".log";
          log="log files from ";
          log+=log1;
          log+=" to ";
          log+=log2;
        }
      Cerr << "Error and standard outputs are redirected into " << log << " for case " << liste_cas[i] << finl;
      // Store in log_courant as it is reused later
      if (Process::me()>=count && Process::me()<=count+n-1)
        log_courant=log;

      ArrOfInt tab(n);
      for (int j = 0; j < n; j++)
        tab[j] = count++;
      PE_Groups::create_group(tab, groupes[i]);
    }
  // Throw an exception if a computation stops, so
  // change the default Process behaviour (MPI_Abort)
  Process::exception_sur_exit=1;
  Cerr << n_calculs << " cases are running..." << finl;
  // Each processor enters its group and interprets the data set
  Nom ancien_nom_du_cas(nom_du_cas());
  const int old_journal_level = get_journal_level();
  for (int i = 0; i < n_calculs; i++)
    {
      if (PE_Groups::enter_group(groupes[i].valeur()))
        {
          set_Cerr_to_journal(1);
          if (disable_journal)
            change_journal_level(0);
          Nom nom_fichier(liste_cas[i]);
          get_set_nom_du_cas() = nom_fichier;
          Journal(1) << "Execute_parallel: Entering subgroup " << i
                     << " to run case " << nom_fichier << finl;

          nom_fichier += ".data";
          {
            // Open the file (the LecFicDiffuse object is created inside
            // the braces so it is destroyed before exiting the group)
            LecFicDiffuse_JDD data_file(nom_fichier);
            data_file.set_check_types(1);
            // Create a new interpreter. At the end of reading the case,
            // the objects will be destroyed.
            Interprete_bloc interp;
            // Use exceptions for computations that stop
            int ok=1;
            try
              {
                interp.interpreter_bloc(data_file,
                                        Interprete_bloc::FIN /* we expect FIN at the end of the file */,
                                        0 /* verifie_sans_interpreter=0 */);
              }
            catch (TRUST_Error& err)
              {
                if (err.get_pe()!=Process::me())
                  {
                    Cerr << err.get_pe() << " <> " << Process::me() << " in Execute_parallel::interpreter." << finl;
                    Process::exit();
                  }
                ok=0;
              }
            set_Cerr_to_journal(0);
            if (ok)
              Cerr << "Case " << liste_cas[i] << " has finished. See " << log_courant << finl;
            else
              Cerr << "!!! Case " << liste_cas[i] << " has failed. See " << log_courant << " !!!" << finl;
          }
          Journal(1) << "Execute_parallel: Exiting subgroup " << i << finl;
          PE_Groups::exit_group();
        }
    }
  // return to the standard behaviour
  Process::exception_sur_exit=0;
  change_journal_level(old_journal_level);
  // Wait for all processors to finish executing their computation.
  barrier();
  get_set_nom_du_cas() = ancien_nom_du_cas;
  Cerr << finl << "End of Execute_parallel::interpreter" << finl;
  return is;
}
