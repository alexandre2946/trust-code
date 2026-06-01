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

#include <Interprete_bloc.h>
#include <Type_Verifie.h>

// `export <Type> <name>` is a TRUST-specific prefix on forward
// declarations — it tells Interprete_bloc::lire to add the resulting
// object to the parent Interprete_bloc's scope instead of the local
// one, so the name survives Read_file's sub-block destruction. There
// is no dedicated C++ class for it; the parsing is inlined below.
// trustify recognises the `export` keyword directly in
// Dataset_Parser.ReadFromTokens — it treats `export <Type> <name>`
// as a regular forward declaration and stashes the `export` token on
// the resulting Declaration_Parser for round-trip emission.

Implemente_instanciable_sans_constructeur_ni_destructeur(Interprete_bloc,"Interprete_bloc",Liste_bloc);

// See Interprete_bloc::interprete_courant()
static OBS_PTR(Interprete_bloc) interprete_courant_;

/*! @brief Returns the Interprete_bloc currently being read from the data set.
 *
 * The current interpreter changes
 *   when an object of type Interprete_bloc is created or destroyed
 *   (for example when entering or leaving a { } block).
 *
 */
Interprete_bloc& Interprete_bloc::interprete_courant()
{
  return interprete_courant_.valeur();
}

Interprete_bloc::Interprete_bloc()
{
  // If a current interpreter exists, it becomes the parent
  // of the interpreter being constructed:
  if (interprete_courant_) pere_ = interprete_courant_;

  interprete_courant_ = *this;
}

Interprete_bloc::~Interprete_bloc()
{
  // Restore the value of the previous current interpreter
  interprete_courant_ = pere_;
}

Sortie& Interprete_bloc::printOn(Sortie& os) const
{
  Process::exit();
  return os;
}

Entree& Interprete_bloc::readOn(Entree& is)
{
  Process::exit();
  return is;
}

/*! @brief Interprets a block of instructions read from input is.
 *
 * If the block begins with {, the opening brace is assumed to have already been read.
 *   The block ends with }, end of file, or the keyword FIN depending on bloc_type
 *   (see mon_main.cpp for example).
 *   If verifie_sans_interpreter!=0, no objects are created and no interpreter
 *   is started; only the consistency of braces is checked (same number of { as }).
 *
 */
Entree& Interprete_bloc::interpreter_bloc(Entree& is, Bloc_Type bloc_type, int verifie_sans_interpreter)
{
  // The output level for journal messages
  const int jlevel = 3;
  Journal(jlevel) << "Interprete_bloc::interpreter_bloc bloc_type=" << (int) bloc_type << finl;
  Motcle motlu;
  is >> motlu;
  while (1)
    {
      // Are we at the end of the block?
      if (is.eof())
        {
          // Some .data files do not end with FIN; reaching eof without finding "FIN" is also accepted.
          if (bloc_type == BLOC_EOF || bloc_type == FIN)
            {
              Journal(jlevel) << "Interprete_bloc: end of file => end of bloc" << finl;
              break;
            }
          else
            {
              Cerr << "Error in Interprete_bloc: unexpected end of file\n" << "check for missig \"}\" or missing FIN keyword" << finl;
              Process::exit();
            }
        }
      // Other possible endings:
      if (motlu == "}" || motlu == "FIN|END")
        {
          if ((motlu == "}" && bloc_type == ACCOLADE) || (motlu == "FIN|END" && bloc_type == FIN))
            {
              Journal(jlevel) << "Interprete_bloc: reading " << motlu << " => end of bloc" << finl;
              break;
            }
          else
            {
              if (motlu == "}")
                Cerr << "Error in Interprete_bloc: extra \"}\" in data file" << finl;
              else
                Cerr << "Error in Interprete_bloc: unexpected FIN, check for missing \"}\"" << finl;
              Process::exit();
            }
        }

      // Is it a comment?
      if (motlu == "#")
        {
          Nom commentaire("# ");
          Nom nom_lu;
          int countleft = 8; // Keep the first 8 words of the comment
          do
            {
              is >> nom_lu;
              if (is.eof())
                break;
              if ((countleft--) > 0)
                {
                  commentaire += nom_lu;
                  commentaire += " ";
                }
            }
          while (nom_lu != "#");
          if (is.eof())
            {
              Cerr << "Error in Interprete_bloc: end of file while reading a comment :\n" << commentaire << "..." << finl;
            }
          // Next keyword to interpret:
          motlu = nom_lu;
        }
      else if (motlu == "{")
        {
          Journal(jlevel) << "Interprete_bloc: reading { => creating new Interprete_bloc" << finl;
          // Is it the beginning of a block?
          // New interpreter:
          Interprete_bloc inter;
          inter.interpreter_bloc(is, ACCOLADE, verifie_sans_interpreter);
        }
      else
        {
          if (verifie_sans_interpreter)
            {
              // If we only want to verify the number of braces,
              // it is enough to ignore the keyword... recursively open
              // new interpreters when { is encountered and close them
              // when } is found. If a } is missing at the end, we get the message
              // "check for missing }" and if there is an extra one we get
              // the message "extra }"
              Journal(jlevel) << "Interprete_bloc: just checking {} => ignore keyword " << motlu << finl;
              verifie(motlu);
            }
          else
            {
              // In the following block, if there is a read error on is, it is an error:
              is.set_error_action(Entree::ERROR_EXIT);
              int export_object = 0;
              if (motlu == "export")
                {
                  Journal(jlevel) << "Exporting object, reading object type..." << finl;
                  export_object = 1;
                  is >> motlu;
                }
              // The keyword must be an instantiable object type
              Journal(jlevel) << "Interprete_bloc: reading " << motlu << " => trying to instanciate object" << finl;
              DerObjU objet;
              objet.typer(motlu);

              if (sub_type(Interprete, objet.valeur()))
                {
                  // If it is an interpreter, call the interpreter method of the object:
                  Journal(jlevel) << " Calling object.interpreter()" << finl;
                  Interprete& inter = ref_cast(Interprete, objet.valeur());
                  inter.interpreter(is);
                }
              else
                {
                  // Not an interpreter; read the object name and store it
                  Journal(jlevel) << " Reading object name" << finl;
                  Nom nom_objet;
                  is >> nom_objet;
                  Journal(jlevel) << " Storing object " << nom_objet << " of type " << motlu << finl;
                  if (!export_object || !pere_)
                    ajouter(nom_objet, objet);
                  else
                    pere_->ajouter(nom_objet, objet);
                }
              is.set_error_action(Entree::ERROR_CONTINUE);
            }
        }
      is >> motlu;
    }
  return is;
}

/*! @brief Returns the Objet_U corresponding to nom contained in this interprete_bloc. If the object does not exist, exit() (no search in the parent).
 *
 */
Objet_U& Interprete_bloc::objet_local(const Nom& nom)
{
  const int i = les_noms_.search(nom);
  if (i < 0)
    {
      Cerr << "Internal error in Interprete::objet_local() : object '" << nom << "' does not exist" << finl;
      Process::exit();
    }
  return operator[](i);
}

/*! @brief Returns a flag indicating whether an object with this name is registered in this interpreter (does not check the parent).
 *
 */
int Interprete_bloc::objet_local_existant(const Nom& nom)
{
  const int i = les_noms_.search(nom);
  return (i >= 0);
}

/*! @brief Adds object ob to the interpreter's object list and names it with nom.
 *
 * If the object already exists, exit().
 *
 */
Objet_U& Interprete_bloc::ajouter(const Nom& nom, DerObjU& ob)
{
  Journal(3) << "Interprete::ajouter(" << nom << ") of type " << ob->que_suis_je() << finl;
  if (les_noms_.search(nom) >= 0)
    {
      Cerr << "Error in Interprete::ajouter: object " << nom << " already exists." << finl;
      Process::exit();
    }
  les_noms_.add(nom);
  Objet_U& obu = add_deplace(ob);
  obu.nommer(nom);
  return obu;
}

/*! @brief Searches for the requested object in the current Interprete_bloc (Interprete_bloc::interprete_courant()) and in all
 *
 *   its successive parents. If the object does not exist, exit().
 *
 */
Objet_U& Interprete_bloc::objet_global(const Nom& nom)
{
  OBS_PTR(Interprete_bloc) ptr(interprete_courant());
  while (ptr)
    {
      Interprete_bloc& interp = ptr.valeur();
      if (interp.objet_local_existant(nom))
        {
          Objet_U& objet = interp.objet_local(nom);
          return objet;
        }
      ptr = interp.pere_;
    }
  // Object not found!
  Cerr << "Error in Interprete::objet: object " << nom << " does not exist." << finl;
  Process::exit();
  return objet_global(nom); // Pour le compilo
}

/*! @brief Returns a flag indicating whether an object with this name exists in interprete_courant() or one of its parents.
 *
 */
int Interprete_bloc::objet_global_existant(const Nom& nom)
{
  OBS_PTR(Interprete_bloc) ptr(interprete_courant());
  while (ptr)
    {
      Interprete_bloc& interp = ptr.valeur();
      if (interp.objet_local_existant(nom))
        return 1;
      ptr = interp.pere_;
    }
  // Object not found!
  return 0;
}
