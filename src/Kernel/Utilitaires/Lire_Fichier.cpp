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

#include <Lire_Fichier.h>
#include <Interprete_bloc.h>
#include <EFichierBin.h>
#include <LecFicDiffuse_JDD.h>
#include <LecFicDiffuse.h>
#include <Read_unsupported_ASCII_file_from_ICEM.h>

Implemente_instanciable(Lire_Fichier,"Lire_Fichier|Read_file",Interprete);
// XD read_file interprete lire_fichier INHERITS_BRACE Keyword to read the object name_obj contained in the file
// XD_CONT filename. NL2 This is notably used when the calculation domain has already been meshed and the mesh contains
// XD_CONT the file filename, simply write read_file dom filename (where dom is the name of the meshed domain). NL2 If
// XD_CONT the filename is ;, is to execute a data set given in the file of name name_obj (a space must be entered
// XD_CONT between the semi-colon and the file name).
// XD attr name_obj chaine name_obj REQ Name of the object to be read.
// XD attr filename chaine filename REQ Name of the file.

/*! @brief Not implemented: calls exit().
 *
 */
Sortie& Lire_Fichier::printOn(Sortie& os) const
{
  Cerr << "Error in Lire_Fichier::printOn()" << finl;
  exit();
  return os;
}

/*! @brief Not implemented: calls exit().
 *
 */
Entree& Lire_Fichier::readOn(Entree& is)
{
  Cerr << "Error in Lire_Fichier::readOn()" << finl;
  exit();
  return is;
}

/*! @brief Two possible syntaxes in the data file: Lire_fichier NOM_OBJET NOM_FICHIER
 *
 *     (reads the file content into the object using the object's readOn method)
 *    Lire_fichier NOM_FICHIER ;
 *     (interprets the file in a local interpreter: objects declared
 *      in the file are destroyed at the end of the file reading)
 *
 */
Entree& Lire_Fichier::interpreter(Entree& is)
{
  Nom nom1, nom2;
  is >> nom1 >> nom2;
  if (nom2 != ";")
    {
      if (is_a_binary_file(nom2))
        {
          check_ICEM_binary_file(nom2,nom1);
          Cerr << "Lire_Fichier: reading binary file " << nom2 << " to object " << nom1 << finl;
          Objet_U& ob1=objet(nom1);
          EFichierBin fic(nom2);
          if(!fic.get_ifstream())
            {
              Cerr << "Unable to open the file " << nom2 << finl;
              Cerr << "Enter a different file name please ... ";
              exit();
            }
          fic >> ob1;
        }
      else
        {
          check_ICEM_ascii_file(nom2,*this);
          Cerr << "Lire_Fichier: reading ASCII file " << nom2 << " to object " << nom1 << finl;
          LecFicDiffuse fic(nom2);
          // Cannot activate check_types because some geom files have a strange format
          // where everything is concatenated like: -0.12000000E+01-0.25000000E+00-0.25000000E+00
          // (exemple croix.geom)
          fic.set_check_types(0);
          Objet_U& ob1 = objet(nom1);
          fic >> ob1;
        }
    }
  else
    {
      Cerr << "Lire_Fichier: interpreting file " << nom1 << finl;
      // Not counting lines inside this file
      LecFicDiffuse_JDD data_file(nom1);
      data_file.track_lines(false);
      data_file.set_check_types(1);
      // Create a new interpreter. At the end of reading,
      // the objects will be destroyed.
      Interprete_bloc interp;
      interp.interpreter_bloc(data_file,
                              Interprete_bloc::BLOC_EOF /* end of block at end of file */,
                              0 /* check_without_interpreting=0 */);
      Cerr << "Lire_Fichier: end of file " << nom1 << finl;
    }
  Cerr << "Lire_Fichier: end of file " << nom1 << finl;
  return is;
}
