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

#include <Read_unsupported_ASCII_file_from_ICEM.h>
#include <Lire_Fichier.h>
#include <SFichierBin.h>
#include <EFichierBin.h>
#include <Interprete.h>
#include <TRUSTTab.h>
#include <Objet_U.h>
#include <Motcle.h>


// PLEASE DO NOT REMOVE Read_unsupported_ASCII_file_from_ICEM KEYWORD FROM TRUST
// IT IS STILL USEFUL (in v1.9.7) TO READ OLD ICEM FILES (ND, 07/11/2025)

Implemente_instanciable(Read_unsupported_ASCII_file_from_ICEM,"Read_unsupported_ASCII_file_from_ICEM",Lire_Fichier);
// XD read_unsupported_ascii_file_from_icem read_file read_unsupported_ascii_file_from_icem INHERITS_BRACE not_set

Sortie& Read_unsupported_ASCII_file_from_ICEM::printOn(Sortie& os) const
{
  return os;
}

Entree& Read_unsupported_ASCII_file_from_ICEM::readOn(Entree& is)
{
  return is;
}

// Used function
inline char read_write_string_from(EFichierBin& s,SFichierBin& o, const char* first_caracter)
{
  // The rule is:
  // Read up to first_caracter if it is defined
  // Pad with the separator 00
  // Return the last character (before the space)
  char octet;
  char zero=0;
  char last_caracter;
  // If first_caracter is defined, read up to it
  if (*first_caracter!=0)
    {
#ifndef NDEBUG
      Cerr << " -> Look for the byte:";
      fprintf(stderr," %02x <=> ",*first_caracter);
      Cerr << (Nom)*first_caracter << finl;
      Cerr << " -> Read the bytes :";
#endif
      //while (s.good() && s.get_istream().read(&octet,1)!=0 && octet!=*first_caracter)
      while (s.good() && s.get_istream().read(&octet,1) && octet!=*first_caracter)
        {
#ifndef NDEBUG
          fprintf(stderr," %02x",octet);
#endif
        };
    }
  else
    {
      // If first_caracter is not defined, simply skip
      // spaces (32, 20 in hex) and separators (00) and delete (127, 7F in hex)
#ifndef NDEBUG
      Cerr << " -> Read the bytes :";
#endif
      //while (s.good() && s.get_istream().read(&octet,1)!=0 && (octet==32 || octet==0 || octet==127) )
      while (s.good() && s.get_istream().read(&octet,1) && (octet==32 || octet==0 || octet==127) )
        {
#ifndef NDEBUG
          fprintf(stderr," %02x",octet);
#endif
        };
    }
#ifndef NDEBUG
  fprintf(stderr," %02x",octet);
#endif
  o.get_ostream().write(&octet,1);
  last_caracter = octet;
  Nom chaine(octet);
  // Read up to the space:
  //while (s.good() && s.get_istream().read(&octet,1)!=0 && octet!=32)
  while (s.good() && s.get_istream().read(&octet,1) && octet!=32)
    {
#ifndef NDEBUG
      fprintf(stderr," %02x",octet);
#endif
      o.get_ostream().write(&octet,1);
      last_caracter = octet;
      chaine+=octet;
    }
#ifndef NDEBUG
  fprintf(stderr," %02x",octet);
#endif
  o.get_ostream().write(&zero,1);
#ifndef NDEBUG
  Cerr << finl;
  Cerr << " -> Read the string: " << chaine << finl;
#endif
  return last_caracter;
}
// Function that checks if we are reading an ICEM binary file and if so
// creates a TRUST-compatible file
void check_ICEM_binary_file(Nom& filename, const Nom& nom_objet_lu)
{
  // Check that the read object is a domain
  Objet_U& objet_lu = Interprete::objet(nom_objet_lu);
  // [ABN] hmmmm ... hopefully the only place where we have this sort of things....:
  if (!(objet_lu.que_suis_je()=="Domaine" || objet_lu.que_suis_je()=="Domaine")) return;

  EFichierBin tmp(filename);
  // An ASCII binary file is recognizable by the fact that the domain name
  // is followed by a space, so the byte is 32 (= 20 in hex = space)
  // In a standard binary file the domain name is followed by 0
  char octet=(char)-1;
  while(octet!=32 && octet!=0)
    tmp.get_istream().read(&octet,1);
  tmp.close();
  if (octet==32) // This is an ICEM binary file
    {
      char zero=0;
      Cerr << "==============================================" << finl;
      Cerr << filename << " is an ICEM binary file." << finl;
      Cerr << "To save space, you can now delete this file and use" << finl;
      // Removal of certain spaces from the binary file
      Nom new_filename(filename);
      new_filename+=".cleaned";
      Cerr << "instead, the newly created " << new_filename << " file." << finl;
      Cerr << "==============================================" << finl;
      EFichierBin s(filename);
      SFichierBin o(new_filename);
      s.get_istream().read(&octet,1);
      while(octet!=32)
        {
          o.get_ofstream().write(&octet,1);
          s.get_istream().read(&octet,1);
        }
      // Write 00 in place of 32
      o.get_ostream().write(&zero,1);
      // Then go up to integer 2
      s.get_istream().read(&octet,1);
      while(octet!=2)
        s.get_istream().read(&octet,1);
      // Go back and start reading the DoubleTab
      s.get_istream().unget();
      // Read nodes
      DoubleTab sommets;
      s >> sommets;
      o << sommets;
      Cerr << "End of the read of the nodes." << finl;
      // "{"
      if (read_write_string_from(s,o,"{")!=123)
        {
          Cerr << "Error in is_a_ICEM_binary_file : { is waited." << finl;
          Process::exit();
        }
      // Domain name
      read_write_string_from(s,o,"");
      // Read element type
      read_write_string_from(s,o,"T"); // TETRAEDRE
      // Read elements
      IntTab elems;
      s >> elems;
      o << elems;
      Cerr << "End of the read of the cells." << finl;
      // Loop over boundaries
      //  {  : 7b in hex, 123 in decimal
      //  ,  : 2C in hex, 44  in decimal
      //  }  : 7D in hex, 125 in decimal
      const char* sep="{";
      while(read_write_string_from(s,o,sep)!=125)
        {
          // Boundary name
          read_write_string_from(s,o,"");
          // Element type
          read_write_string_from(s,o,"T"); // TRIANGLE_3D
          // Boundary elements
          s >> elems;
          o << elems;
          // Neighbor faces
          IntTab faces_voisins;
          s >> faces_voisins;
          o << faces_voisins;
          Cerr << "End of the read of a boundary." << finl;
          sep=""; // , or } cannot be known in advance...
        }
      // Read the 3 "vide" entries
      for (int i=0; i<3; i++)
        {
          if (read_write_string_from(s,o,"v")!=101)
            {
              Cerr << "Error in is_a_ICEM_binary_file : 'vide' is waited." << finl;
              Process::exit();
            }
        }
      // Read }
      if (read_write_string_from(s,o,"}")!=125)
        {
          Cerr << "Error in is_a_ICEM_binary_file : } is waited." << finl;
          Process::exit();
        }
      // Read the last "vide" entry
      read_write_string_from(s,o,"v");
      o.close();
      s.close();
      // Point to the new file
      filename=new_filename;
      Cerr << "=============================================" << finl;
      Cerr << "You can use now: " << finl;
      Cerr << "Lire_fichier " << nom_objet_lu << " "<< filename  << finl;
      Cerr << "=============================================" << finl;
    }
}

void check_ICEM_ascii_file(const Nom& filename, const Lire_Fichier& keyword)
{
  Nom ICEM(filename);
  if (ICEM.finit_par(".asc") && !sub_type(Read_unsupported_ASCII_file_from_ICEM,keyword))
    {
      Cerr << "Since the 1.6.7 version, the read of an ASCII file from ICEM mesh tool is not" << finl;
      Cerr << "supported anymore cause insufficient number of digits to define the node coordinates." <<finl;
      Cerr << "Please, use the TRUST binary export from ICEM." << finl;
      Process::exit();
    }
  else if (sub_type(Read_unsupported_ASCII_file_from_ICEM,keyword))
    {
      Cerr << "====================================================" << finl;
      Cerr << "Warning: You are using an obsolete feature of TRUST" << finl;
      Cerr << "with Unsupported_ASCII_file_read_from_ICEM keyword." << finl;
      Cerr << "Use binary format for ICEM file!" << finl;
      Cerr << "====================================================" << finl;
    }
}
