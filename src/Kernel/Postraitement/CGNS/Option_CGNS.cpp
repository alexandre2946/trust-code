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

#include <Option_CGNS.h>
#include <Param.h>

Implemente_instanciable(Option_CGNS, "Option_CGNS", Interprete);
// XD Option_CGNS interprete Option_CGNS BRACE Class for CGNS options.

bool Option_CGNS::PARALLEL_OVER_ZONE = false; /* NOT BY DEFAULT */
bool Option_CGNS::USE_LINKS = false; /* NOT BY DEFAULT */
bool Option_CGNS::LINKED_FILES_PER_COMM_GROUP = false; /* NOT BY DEFAULT */
bool Option_CGNS::SINGLE_FILE_PER_COMM_GROUP = false; /* NOT BY DEFAULT */
bool Option_CGNS::KEEP_FILES_BEFORE_RESET_TIME = false; /* NOT BY DEFAULT */
int Option_CGNS::CLOSE_EVERY_N = -1; /* -1 BY DEFAULT => never opened/closed */
int Option_CGNS::FLUSH_EVERY_N = 1; /* 1 BY DEFAULT => flush each dt post */

Sortie& Option_CGNS::printOn(Sortie& os) const { return Interprete::printOn(os); }
Entree& Option_CGNS::readOn(Entree& is) { return Interprete::readOn(is); }

Entree& Option_CGNS::interpreter(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter_flag("PARALLEL_OVER_ZONE", &PARALLEL_OVER_ZONE); // XD_ADD_P rien If used, data will be written in separate zones (ie: one zone per processor). This is not so performant but easier to read later ...
  param.ajouter_flag("USE_LINKS", &USE_LINKS); // XD_ADD_P rien If used, data will be written in separate files; one file for mesh, and then one file for solution time. Links will be used.
  param.ajouter_flag("LINKED_FILES_PER_COMM_GROUP", &LINKED_FILES_PER_COMM_GROUP); // XD_ADD_P rien If used, data will be written (at each comm group) in separate files; one file for mesh, and then one file for solution time. Links will be used.
  param.ajouter_flag("SINGLE_FILE_PER_COMM_GROUP", &SINGLE_FILE_PER_COMM_GROUP); // XD_ADD_P rien If used, data will be written (at each comm group) in a single file.
  param.ajouter_flag("KEEP_FILES_BEFORE_RESET_TIME", &KEEP_FILES_BEFORE_RESET_TIME); // XD_ADD_P rien If used with resetTime, CGNS files will be kept. Otherwise the files are overwritten.
  param.ajouter("CLOSE_EVERY_N", &CLOSE_EVERY_N); // XD_ADD_P entier Used to fix the opening/closing frequency when writing in a single CGNS file (choice by defaut).
  param.ajouter("FLUSH_EVERY_N", &FLUSH_EVERY_N); // XD_ADD_P entier Used to fix the flush-to-disc frequency when writing in a single CGNS file (choice by defaut).
  param.lire_avec_accolades_depuis(is);

  const bool single_file = (!USE_LINKS && !LINKED_FILES_PER_COMM_GROUP);

  if ((PARALLEL_OVER_ZONE || SINGLE_FILE_PER_COMM_GROUP) && !single_file)
    {
      Cerr << finl << "Error in Option_CGNS :" << finl;
      Cerr << "   You can not activate the option 'PARALLEL_OVER_ZONE' and/or 'SINGLE_FILE_PER_COMM_GROUP' with 'USE_LINKS' and/or 'LINKED_FILES_PER_COMM_GROUP' !!!" << finl;
      Cerr << "   The default CGNS post behavior is a single file for all times/fields." << finl;
      Cerr << "   PARALLEL_OVER_ZONE and/or SINGLE_FILE_PER_COMM_GROUP options remain in this context; ie: single file." << finl;
      Cerr << "   USE_LINKS and/or LINKED_FILES_PER_COMM_GROUP options require several linked files." << finl;
      Cerr << "   See the doc for more details or contact the TRUST support ... " << finl << finl;
      Process::exit();
    }

  Cerr << finl << "*********************************" << finl;
  Cerr << "********** Option_CGNS **********" << finl;
  Cerr << "*********************************" << finl << finl;

  if (single_file)
    {
      if (FLUSH_EVERY_N > 0)
        Cerr << " Data will be flushed into a single CGNS file each " << FLUSH_EVERY_N << " time post ..." << finl;

      if (CLOSE_EVERY_N > 0)
        Cerr << " A single CGNS file will be closed and opened each " << CLOSE_EVERY_N << " time post ..." << finl;

      if (PARALLEL_OVER_ZONE)
        Cerr << " PARALLEL_OVER_ZONE => CGNS data will be written in separate zones ..." << finl;

      if (SINGLE_FILE_PER_COMM_GROUP)
        Cerr << " SINGLE_FILE_PER_COMM_GROUP => A single CGNS file will be written in each COMM group ... ..." << finl;
    }
  else
    {
      if (LINKED_FILES_PER_COMM_GROUP)
        {
          USE_LINKS = true;
          Cerr << " LINKED_FILES_PER_COMM_GROUP => Several linked CGNS files will be written in each COMM group ..." << finl;
          Cerr << " This will activate also the option USE_LINKS ..." << finl;
        }

      if (USE_LINKS)
        Cerr << " USE_LINKS => CGNS data will be written in separate files (mesh, solution ...)" << finl;
    }

  Cerr << finl << "*********************************" << finl << finl;

  return is;
}
