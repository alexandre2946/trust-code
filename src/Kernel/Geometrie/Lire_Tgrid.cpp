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

#include <NettoieNoeuds.h>
#include <Lire_Tgrid.h>
#include <Frontiere.h>
#include <TRUST_Ref.h>
#include <EFichier.h>
#include <Domaine.h>
#include <ctype.h>
#include <math.h>

// Method to be moved to a base class above all lire_... classes,
// or to be placed in the EFichier constructor with optional deletion
// of the decompressed file in the EFichier destructor.
inline void decompression(Nom& nom_fichier)
{
  Nom tmp(nom_fichier);
  if (tmp.prefix(".gz")!=nom_fichier)
    {
      Cerr << "Unzipping file " << nom_fichier << " ..." << finl;
      Nom cmd="gunzip -c ";
      cmd+=nom_fichier+" > "+tmp;
      Cerr << (int)system(cmd) << finl;
      nom_fichier=tmp;
    }
}
Implemente_instanciable(Lire_Tgrid,"Lire_Tgrid",Interprete_geometrique_base);
// XD read_tgrid interprete lire_tgrid INHERITS_BRACE Keyword to reaf Tgrid/Gambit mesh files. 2D (triangles or
// XD_CONT quadrangles) and 3D (tetra or hexa elements) meshes, may be read by TRUST.
// XD attr dom ref_domaine dom REQ Name of domaine.
// XD attr filename chaine filename REQ Name of file containing the mesh.

Sortie& Lire_Tgrid::printOn(Sortie& os) const { return Interprete::printOn(os); }

Entree& Lire_Tgrid::readOn(Entree& is) { return Interprete::readOn(is); }

int chartoint(char c)
{
  if( (0>(c-'0')) || ((c-'0')>9) )
    return -1;
  else
    return (c-'0');
}

// Converts a hex-format string (e.g. "00aad22") to an integer value
int htoi(const char * szChaine)
{
  int lResult = 0;
  int iLength = 0;
  // Null pointer: return -1
  if (szChaine == nullptr)
    return -1;
  // Compute the string length
  iLength = (int)strlen(szChaine);

  // Convert the string to uppercase in a new string (allocated by strdup)
  char * szHexaString = strdup(szChaine);
  // For each character, compute its integer value
  for (int i = iLength - 1; i >= 0; i--)
    {
      char cCharacter = szHexaString[i];
      int iValue = 0;
      // It's a digit, convert it to integer
      if (isdigit(cCharacter))
        {
          //iValue = atoi(&cCharacter);
          // atoi seems badly buggy!
          // Using a correct chartoint implementation:
          iValue = chartoint(cCharacter);
          if (iValue>9)
            {
              Cerr << "iValue is worth " << iValue << " ! " << finl;
              Process::exit();
            }
        }
      // It's a letter, assign it an integer value
      else if (isalpha(cCharacter))
        {
          switch(cCharacter)
            {
            case 'A' :
            case 'a' :
              iValue = 10;
              break;
            case 'B' :
            case 'b' :
              iValue = 11;
              break;
            case 'C' :
            case 'c' :
              iValue = 12;
              break;
            case 'D' :
            case 'd' :
              iValue = 13;
              break;
            case 'E' :
            case 'e' :
              iValue = 14;
              break;
            case 'F' :
            case 'f' :
              iValue = 15;
              break;
            default  :
              return -3;              // Invalid hex character.
            }
        }
      // Neither a letter nor a digit.
      else
        return -4;   // Invalid character.
      //lResult += iValue * pow(16, iLength - i - 1);
      for (int puissance=0; puissance<iLength - i - 1; puissance++)
        iValue *= 16;
      lResult += iValue;
    }
  // Free the allocated string
  free (szHexaString);
  //Cerr << lResult << finl;
  return lResult;
}

inline void va_a_la_parenthese_fermante(EFichier& fic)
{
  int parenthese_ouverte=1;
  Nom lu;
  while ((parenthese_ouverte!=0)&&(fic.good()))
    {
      fic >> lu;
      const char* chaine = lu.getChar();
      Process::Journal()<<"|"<<chaine<<"|"<<finl;
      size_t iLength = strlen(chaine);
      for (size_t i=0; i<iLength; i++)
        {
          char c = chaine[i];
          switch(c)
            {
            case 40 :
              parenthese_ouverte++;
              break;
            case 41 :
              parenthese_ouverte--;
              break;
            }
        }
    }

  if (parenthese_ouverte!=0)
    {
      Cerr<< "Error in file"<<finl;
      Process::exit();
    }

}
inline void va_a_la_parenthese_ouvrante(EFichier& fic)
{
  int ok=0;
  Nom lu;
  while (ok!=1)
    {
      fic >> lu;
      const char* chaine = lu.getChar();
      size_t iLength = strlen(chaine);
      for (size_t i=0; i<iLength; i++)
        {
          char c = chaine[i];
          switch(c)
            {
            case 40 :
              ok=1;
              break;
            }
        }
    }
}

/*! @brief Main function of the Lire_Tgrid interpreter. Reads a Tgrid mesh file.
 *
 * With 2 arguments nom1 and nom2, reads the object from file nom2 into object nom1.
 * With a single argument nom1, interprets the file named nom1.
 *
 * @param is An input stream.
 * @return The modified input stream.
 */
Entree& Lire_Tgrid::interpreter_(Entree& is)
{
  Cerr << "Reading a mesh which comes from Tgrid" << finl;
  associer_domaine(is);
  Domaine& dom=domaine();
  DoubleTab& coord_sommets=dom.les_sommets();
  // Variable declarations
  int dim = -1;
  int nb_som = 0;
  int nb_elem = 0;
  int nb_face = 0;
  int nb_som_elem = 0;
  int type_elements = 0;
  int compteur = 0;
  // Work array
  ArrOfInt nb_som_lu_elem;
  Motcle motlu;
  // File management
  Nom nom_fichier;
  is >> nom_fichier;
  decompression(nom_fichier);
  Cerr << "Reading of the file " << nom_fichier << " ..." << finl;
  EFichier lecture(nom_fichier);
  // Check immediately that the file exists, because if it does not,
  // execution hangs in eof();
  if (!lecture.good())
    {
      Cerr << "Problem to open the file " << nom_fichier << finl;
      Cerr << "There is maybe an error in the filename." << finl;
      exit();
    }

  // First pass through the .msh file
  // to find the number of elements, as it is sometimes placed at the end of the file!
  while (!lecture.eof())
    {
      lecture >> motlu;
      if (motlu=="(12")
        {
          lecture >> motlu;
          if (motlu=="(0")
            {
              lecture >> motlu;        // Index of the first element
              lecture >> motlu;        // Number of elements
              nb_elem=htoi(motlu);
              Cerr << "The total number of elements to read is " << nb_elem << finl;
              lecture >> motlu;        // Domain type (0=dead; 1=active; 32=inactive) or "0))"
              if (motlu!="0))") lecture >> motlu; // Skip the type if it exists
            }
          else if (motlu=="(id") lecture >> motlu;        // skip reading the description of the (12 tag
          else
            {
              lecture >> motlu;        // Index of the first element
              lecture >> motlu;        // Number of elements
              lecture >> motlu;        // Domain type (1:fluid or 0x11:solid)
              lecture >> motlu;        // Element type
              if (motlu=="1))")
                {
                  // Reading triangles
                  type_elements=1;
                  Cerr << "2D elements of type Triangle" << finl;
                }
              else if (motlu=="3))")
                {
                  // Reading quadrangles
                  type_elements=3;
                  Cerr << "2D elements of type Quadrangle" << finl;
                }
              else if (motlu=="2))")
                {
                  // Reading tetrahedra
                  type_elements=2;
                  Cerr << "3D elements of type Tetrahedron" << finl;
                }
              else if (motlu=="4))")
                {
                  // Reading hexahedra
                  type_elements=4;
                  Cerr << "3D elements of type Hexahedron" << finl;
                }
              else
                {
                  // Unknown element type
                  Cerr << "Elements unknown !" << finl;
                  Cerr << "It should probably crashed !!!!!" << finl;
                }
              lecture.close();
            }
        }
    }
  // Second pass through the .msh file
  EFichier fic(nom_fichier);
  while (!fic.eof())
    {
      fic >> motlu;
      if (motlu=="(0")
        {
          Cerr << "Reading a comment:" << finl;
          Cerr << motlu;
          va_a_la_parenthese_fermante(fic);
          Cerr << finl << finl;
        }
      else if (motlu=="(1")
        {
          Cerr << "Reading a header:" << finl;
          Cerr << motlu;
          va_a_la_parenthese_fermante(fic);
          Cerr << finl << finl;
        }
      else if (motlu=="(2")
        {
          Cerr << "Reading of the dimension of the case:" << finl;
          fic >> motlu;
          dim=atoi(motlu.prefix(")"));
          if (dim==3) Cerr << "Dimension 3." << finl;
          else if (dim==2) Cerr << "Dimension 2." << finl;
          else
            {
              Cerr << "Dimension " << dim << " of the mesh not provided." << finl;
              exit();
            }
          Cerr << finl;
        }
      else if (motlu=="(10")
        {
          fic >> motlu;
          if (motlu=="(0")
            {
              fic >> motlu;        // Index of the first vertex
              fic >> motlu;        // Number of vertices
              nb_som=htoi(motlu);
              Cerr << "The total number of nodes to read is " << nb_som << finl;
              // Resize the vertex array
              coord_sommets.resize(nb_som,dim);
              // Format-dependent
              fic >> motlu;
              Nom tmp=motlu;
              if (tmp==motlu.prefix("))"))
                fic >> motlu;
            }
          else
            {
              int idomaine=htoi(motlu.suffix("("));
              fic >> motlu;        // Start index
              int ideb=htoi(motlu);
              fic >> motlu;        // End index
              int ifin=htoi(motlu);
              Cerr << ifin-ideb+1 << " nodes are read in the area " << idomaine << finl;
              // Format-dependent: advance to the opening parenthesis
              va_a_la_parenthese_ouvrante(fic);
              /*
                fic >> motlu;        // Type (0: virtual, 1:any, 2:boundary)
                if (motlu.prefix(")"))
                fic >> motlu;        // Dimension
                assert(htoi(motlu.prefix(")"))==dim);
                fic >> motlu;        // ( */

              for (int i=ideb-1; i<ifin; i++)
                for (int j=0; j<dim; j++)
                  fic >> coord_sommets(i,j);
              fic >> motlu;        // ))
            }
          Cerr << finl;
        }
      else if (motlu=="(12")
        {
          fic >> motlu;
          if (motlu=="(0")
            {
              // Information already obtained in the first pass
              fic >> motlu;        // Index of the first element
              fic >> motlu;        // Number of elements
              fic >> motlu;        // Domain type (0=dead; 1=active; 32=inactive) or "0))"
              if (motlu!="0))") fic >> motlu; // Skip the type if it exists
            }
          else
            {
              int idomaine=htoi(motlu.suffix("("));
              fic >> motlu;        // Debut
              //int ideb=htoi(motlu);
              fic >> motlu;        // Fin
              //int ifin=htoi(motlu);
              fic >> motlu;        // Type (1:fluid, Ox11:solid)
              int type=htoi(motlu);
              Cerr << "The type of area " << idomaine << " is " << type << " (1:fluid, 17:solid)" << finl;
              fic >> motlu;        // Type cell (0:mixed,1:tri,2:tetra,3:quad,4:hexa,5:pyramid,6:wedge)
              if (motlu=="2))")
                {
                  // Reading tetrahedra
                  dom.type_elem().typer("Tetraedre");
                }
              else if (motlu=="4))")
                {
                  // Reading hexahedra
                  dom.type_elem().typer("Hexaedre_VEF");
                }
              else if (motlu=="1))")
                {
                  // Reading triangles
                  dom.type_elem().typer("Triangle");
                }
              else if (motlu=="3))")
                {
                  // Reading quadrangles
                  dom.type_elem().typer("Quadrangle");
                }
              else
                {
                  Cerr << "Reading the elements is not provided in this interpreter." << finl;
                  Cerr << "Indeed, we read faces to reconstruct the elements." << finl;
                  Cerr << "Contact TRUST support." << finl;
                  exit();
                }
              dom.type_elem()->associer_domaine(dom);
            }
          Cerr << finl;
        }
      else if (motlu.debute_par("(13"))
        {
          motlu.suffix("(13");
          if (motlu=="")
            fic >> motlu;
          if (motlu=="(0")
            {
              fic >> motlu;        // Index of the first face
              fic >> motlu;        // Number of faces
              nb_face=htoi(motlu);
              Cerr << "The total number of faces to read is " << nb_face << finl;
              fic >> motlu;        // Face type or "0))"
              if (motlu != "0))") fic >> motlu;        // Skip the type if it exists
            }
          else
            {
              int idomaine=htoi(motlu.suffix("("));
              fic >> motlu;        // Start index
              int ideb=htoi(motlu);
              fic >> motlu;        // End index
              int ifin=htoi(motlu);
              fic >> motlu;        // Type (2: interior, >2: boundary condition): this does not seem accurate (see page C-8)!
              int type=htoi(motlu);
              fic >> motlu;        // Face type (0:mixed, 2:linear, 3:triangular, 4:quadrilateral)
              // Format-dependent:
              Nom tmp=motlu;
              int nb_som_face,mixte=0;
              if (tmp!=motlu.prefix(")"))
                fic >> tmp;        // (
              else
                motlu.prefix(")(");
              nb_som_face=htoi(motlu);
              if (nb_som_face==0)
                {
                  mixte=1;                // On va lire des elements mixtes en esperant qu'ils sont de meme type
                  fic >> motlu;        // On lit la premiere ligne
                  nb_som_face=htoi(motlu);
                }
              if (nb_som_face==3)
                nb_som_elem=4;         // On va lire des triangles donc les elements sont des tetras (4 sommets)
              else if (nb_som_face==4)
                nb_som_elem=8;        // On va lire des quadrangles donc les elements sont des hexas (8 sommets)
              else if (nb_som_face==2 && type_elements==1)
                nb_som_elem=3;        // On va lire des segments dont les elements sont des triangles (3 sommets)
              else if (nb_som_face==2 && type_elements==3)
                nb_som_elem=4;        // On va lire des segments dont les elements sont des quadrangles (4 sommets)
              else
                {
                  if (nb_som_face==0)
                    Cerr << "It seems that we try to read faces with different types..." << finl;
                  //                        if (nb_som_face==2)                                                        // N'est plus necessaire, vu qu'on sait lire
                  //                           Cerr << "Il semble que l'on essaie de lire un segment..." << finl;        // maintenant des segments
                  Cerr << "The case of faces that are not triangles or quadrangles" << finl;
                  Cerr << "or a mixture of several types of faces" << finl;
                  Cerr << "is not yet provided. Your mesh is not only composed" << finl;
                  Cerr << "of tetrahedra, hexahedra, triangles or quadrangles." << finl;
                  exit();
                }
              // Tableaux de travail
              IntTab& les_elems=dom.les_elems();
              if (les_elems.size()==0)
                {
                  // Si les_elems n'est pas dimensionne (1ere lecture de faces)
                  les_elems.resize(nb_elem,nb_som_elem);
                  // On initialise a -1 pour voir les sommets non definis
                  les_elems=-1;
                  nb_som_lu_elem.resize_array(nb_elem);
                  nb_som_lu_elem=0;
                }
              ArrOfInt elem(2),som(nb_som_face);
              int nb_face_lu=ifin-ideb+1;
              OBS_PTR(Frontiere) nouveau_bord;
              // Read the face vertices and the 2 elements adjacent to the face
              for (int i=0; i<nb_face_lu; i++)
                {
                  if (mixte && i>0)
                    {
                      fic >> motlu;
                      if (htoi(motlu)!=nb_som_face)
                        {
                          Cerr << "We read an element to " << htoi(motlu) << " faces." << finl;
                          Cerr << "So it was planned to read elements to " << nb_som_face << " faces."<< finl;
                          exit();
                        }
                    }
                  for (int j=0; j<nb_som_face; j++)
                    {
                      fic >> motlu;
                      som[j]=htoi(motlu)-1;
                    }
                  // For hexahedra due to TRUST numbering, swap vertices 2 and 3
                  if (nb_som_face==4)
                    {
                      int tmp2=som[2];
                      som[2]=som[3];
                      som[3]=tmp2;
                    }
                  fic >> motlu;        // first element neighboring the face (zero if boundary)
                  elem[0]=htoi(motlu)-1;
                  fic >> motlu;        // second element neighboring the face (zero if boundary)
                  // Start addition by Cyril MALOD : 15-06-2006
                  Nom tmp2=motlu;
                  if (tmp2!=motlu.prefix("))"))
                    {
                      tmp2=motlu.prefix("))");
                      elem[1]=htoi(tmp2)-1;
                      compteur=1; // This counter prevents reading the next "motlu"
                    }
                  else
                    {
                      if (tmp2!=motlu.prefix(")"))
                        {
                          tmp2=motlu.prefix(")");
                          elem[1]=htoi(motlu)-1;
                          compteur=-1; // This counter triggers reading the next "motlu"
                        }
                      else
                        {
                          elem[1]=htoi(motlu)-1;
                          compteur=0; // This counter triggers reading the next "motlu"
                        }
                    }

                  //                        elem(1)=htoi(motlu)-1;
                  // End addition by Cyril MALOD : 15-06-2006
                  // First pass: verify that we are reading a boundary
                  if (i==0)
                    {
                      if (elem[0]<0 || elem[1]<0)
                        {
                          // Confirmed as a boundary: allocate the necessary structures
                          type=3;
                          Cerr << nb_face_lu << " faces are read from the boundary number " << idomaine << finl;
                          nouveau_bord=dom.faces_bord().add(Bord());
                          nouveau_bord->nommer((Nom)idomaine);
                          if (nb_som_face==3)
                            nouveau_bord->faces().typer(Type_Face::triangle_3D);
                          else if (nb_som_face==4)
                            nouveau_bord->faces().typer(Type_Face::quadrangle_3D);
                          else if (nb_som_face==2)
                            nouveau_bord->faces().typer(Type_Face::segment_2D);
                          else
                            {
                              Cerr << "Type of boundary face not provided for nb_som_face=" << nb_som_face << finl;
                              exit();
                            }
                          nouveau_bord->faces().dimensionner(nb_face_lu);
                        }
                      else
                        {
                          type=2;
                          Cerr << nb_face_lu << " internal faces are read in the area " << idomaine << finl;
                        }
                    }
                  // Internal faces are read but not stored; boundary faces are stored
                  if (type!=2)
                    for (int j=0; j<nb_som_face; j++)
                      nouveau_bord->faces().sommet(i,j)=som[j];

                  // Build the les_elems array from the faces read
                  for (int i2=0; i2<2; i2++)
                    {
                      if (elem[i2]>=0) // Skip elements with index -1 (neighbors of boundary faces)
                        {
                          // First fill of elem(i)
                          if (nb_som_lu_elem[elem[i2]]==0)
                            {
                              for (int j=0; j<nb_som_face; j++)
                                les_elems(elem[i2],nb_som_lu_elem[elem[i2]]++)=som[j];
                            }
                          else
                            {
                              int face_opposee=1;
                              for (int j=0; j<nb_som_face; j++)
                                {
                                  // Add the vertex if it is not already in les_elems
                                  int k=0,trouve=0;
                                  while (k<nb_som_lu_elem[elem[i2]] && trouve==0)
                                    {
                                      if (les_elems(elem[i2],k)==som[j])
                                        {
                                          face_opposee=0;
                                          trouve=1;
                                        }
                                      else
                                        k++;
                                    }
                                  // Only complete for tetrahedra or triangles
                                  if ((trouve==0 && nb_som_face==3) || (trouve==0 && nb_som_face==2 && type_elements==1))
                                    les_elems(elem[i2],nb_som_lu_elem[elem[i2]]++)=som[j];

                                  assert(elem[i2]<nb_elem);
                                  if (nb_som_lu_elem[elem[i2]]>nb_som_elem)
                                    {
                                      Cerr << "Problem on reading the element " << elem[i2] << finl;
                                      Cerr << "There is more than " << nb_som_elem << " nodes !" << finl;
                                      Cerr << "Check that the read file contains only tetrahedra or triangles" << finl;
                                      Cerr << "or contact TRUST support." << finl;
                                      exit();
                                    }
                                }
                              // Fill the opposite face of the hexahedron or quadrangle
                              if ((nb_som_face==4 && face_opposee==1) || (nb_som_face==2 && face_opposee==1 && type_elements==3))
                                {
                                  for (int j=0; j<nb_som_face; j++)
                                    les_elems(elem[i2],nb_som_lu_elem[elem[i2]]++)=som[j];
                                }
                            }
                        }
                    }
                }
              if (compteur==0 || compteur==-1)
                {
                  fic >> motlu;        // ) or ))
                  if (motlu==")" && compteur==0)
                    fic >> motlu;        // )
                }
            }
          Cerr << finl;
        }
      else if ((motlu=="(45") || (motlu=="(39"))
        {
          Cerr << "Reading of a name:" << finl;
          fic >> motlu;        // Domain number
          // Note: the domain number is in decimal!
          //int idomaine=htoi(motlu.suffix("("));
          int idomaine=atoi(motlu.suffix("("));
          fic >> motlu;        // Domain type
          Nom Nomdomaine;
          fic >> Nomdomaine;        // Domain name)())
          Nom nom_domaine=Nomdomaine;
          nom_domaine.prefix(")())");
          if (nom_domaine==Nomdomaine)
            if (1)
              {
                // line break?
                Nom app;
                fic >> app;
                Nomdomaine+=app;
                nom_domaine=Nomdomaine;
                nom_domaine.prefix(")())");

              }
          Cerr << "The area " << idomaine << " is called " << nom_domaine << finl;
          // Iterate through boundaries to rename them
          Bords& les_bords=dom.faces_bord();
          les_bords.associer_domaine(dom);
          if (les_bords.est_vide())
            {
              Cerr << "Reading a name before reading the boundaries..." << finl;
              Cerr << "Case not provided, contact TRUST support." << finl;
              exit();
            }

          for (auto& itr : les_bords)
            if (itr.le_nom()==(Nom)idomaine) itr.nommer(nom_domaine);

          Cerr << finl;
        }
      else
        {
          if (motlu.debute_par("("))
            {
              Cerr << "Reading a tag:" << finl;
              Cerr << motlu;
              va_a_la_parenthese_fermante(fic);
              Cerr << finl << finl;
            }
          else if ( (motlu==" ") || (motlu=="") )
            {
              Cerr << "End of file ?" << finl;
            }
          else if (motlu!="??")
            {
              Cerr << "Tag " << motlu << " unrecognized." << finl;
              exit();
            }
        }
    }
  // Verify that the les_elems array has been fully filled
  IntTab& les_elems=dom.les_elems();
  for (int i=0; i<nb_elem; i++)
    for (int j=0; j<nb_som_elem; j++)
      if (les_elems(i,j)==-1)
        {
          Cerr << "The array of connectivity elements-nodes is wrong filled in Lire_Tgrid::interpreter." << finl;
          Cerr << "Contact TRUST support." << finl;
          exit();
        }

  // Reorder the domain (useful especially for hexahedra)
  dom.type_elem()->reordonner();

  // Clean the domain to remove unused nodes
  // A method should be added to Domaine::nettoie
  // Note: the following lines are not compatible with TRUST < v1.4.6

  if (Process::is_sequential() && (NettoieNoeuds::NettoiePasNoeuds==0) )
    NettoieNoeuds::nettoie(dom);

  return is;
}
