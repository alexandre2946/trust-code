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

#include <Domaine.h>
#include <TRUSTList.h>
#include <TRUSTTabs.h>
#include <Sous_Domaine.h>
#include <Scatter.h>
#include <Poly_geom_base.h>
#include <Octree.h>
#include <Octree_Double.h>
#include <Periodique.h>
#include <Reordonner_faces_periodiques.h>
#include <Frontiere_dis_base.h>
#include <Frontiere.h>
#include <Conds_lim.h>
#include <NettoieNoeuds.h>
#include <Polyedre.h>
#include <TRUST_2_MED.h>
#include <Comm_Group_MPI.h>
#include <Option_Interpolation.h>
#include <Array_tools.h>
#include <Schema_Comm.h>
#include <Interprete_bloc.h>
#include <Extraire_surface.h>
#include <Domaine_VF.h>
#include <MD_Vector_std.h>
#include <MD_Vector_seq.h>
#include <Reorder_Mesh.h>
#include <Perf_counters.h>

Implemente_instanciable_sans_constructeur_32_64( Domaine_32_64, "Domaine", Domaine_base );
// XD domaine domaine_base domaine INHERITS_BRACE Keyword to create a domain.
// XD domaine_64 domaine_base domaine_64 INHERITS_BRACE Keyword to create a big (64b) domain.

// Anonymous namespace for all local methods to this translation unit
namespace
{

static double cached_memory = 0;

template<class LIST_FRONTIERE>
void check_frontiere(const LIST_FRONTIERE& list, const char *msg)
{
  int n = list.size();
  if (!is_parallel_object(n))
    {
      Cerr << " Fatal error: processors don't have the same number of boundaries " << msg << finl;
      Process::barrier();
      Process::exit();
    }
  for (int i = 0; i < n; i++)
    {
      const Nom& nom = list[i].le_nom();
      Cerr << "  Boundary " << msg << " : " << nom << finl;
      if (!is_parallel_object(nom))
        {
          Cerr << " Fatal error: processors don't have the same number of boundaries " << msg << finl;
          Process::barrier();
          Process::exit();
        }
    }
}

// If there is a proc with a type different from vide_OD, set the face type on all
// processors to that type:
template <class _SIZE_>
void corriger_type(Faces_32_64<_SIZE_>& faces, const OWN_PTR(Elem_geom_base_32_64<_SIZE_>)& type_elem)
{
  Type_Face typ = faces.type_face();
  const int pe = (faces.type_face() == Type_Face::vide_0D) ? Process::nproc() - 1 : Process::me();
  const int min_pe = Process::mp_min(pe);
  // Processor min_pe broadcasts its type to all others
  int typ_commun_i = static_cast<int>(typ);
  envoyer_broadcast(typ_commun_i, min_pe);
  Type_Face typ_commun = static_cast<Type_Face>(typ_commun_i);

  if (typ_commun != typ)
    {
      if (typ != Type_Face::vide_0D)
        {
          Cerr << "Error in Domaine.cpp corriger_type: invalid boundary face type" << finl;
          Process::exit();
        }
      faces.typer(typ_commun);
      int n = type_elem->nb_som_face();
      faces.les_sommets().resize(0, n);
    }
}

}  // end anonymous namespace

/*! @brief Reset the Domaine completely except for its name.
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::clear()
{
  sommets_.reset();
  renum_som_perio_.reset();
  elem_ = OWN_PTR(Elem_geom_base_32_64<_SZ_>)();
  mes_elems_.reset();
  aretes_som_.reset();
  elem_aretes_.reset();
  mes_faces_bord_.vide();
  mes_faces_raccord_.vide();
  mes_bords_int_.vide();
  mes_groupes_faces_.vide();
  mes_faces_joint_.vide();
  ind_faces_virt_bord_.reset();
  cg_moments_.reset();
  elem_virt_pe_num_.reset();
  domaines_frontieres_.vide();
  les_ss_domaines_.vide();

  moments_a_imprimer_ = 0;
  bords_a_imprimer_.vide();
  bords_a_imprimer_sum_.vide();

  bords_perio_.clear();

  epsilon_ = Objet_U::precision_geom;
  fichier_lu_ = Nom();

#ifdef MEDCOUPLING_
  mc_mesh_.nullify();
  rmps.clear();
#endif

  volume_total_ = -1;
}

/*! @brief Writes the Domain to an output stream.
 *
 * Writes the name, element type, elements,
 *     boundaries, periodic boundaries, joints,
 *     connections and internal boundaries.
 *
 * @param (Sortie& s) an output stream
 * @return (Sortie&) the modified output stream
 */
template<typename _SZ_>
Sortie& Domaine_32_64<_SZ_>::printOn(Sortie& s) const
{
  Cerr << "Writing of " << nb_som() << " nodes." << finl;
#ifdef SORT_POUR_DEBOG
  s.setf(ios::scientific);
  s.precision(20);
#endif
  s << nom_ << finl;
  s << sommets_;

  // Now write what was formerly the "Zon-e-s" (before TRUST 1.9.2):
  // Write them in the form of a list with a single element, for backward compat (Domains used to contain a list of Zon-e-s)
  s << "{" << finl;
  Cerr << "Writing of " << nb_elem() << " elements." << finl;
  s << "DUMMY_ZONE" << finl; // really just to keep a name here for backward compat
  s << elem_ << finl;
  s << mes_elems_;
  s << mes_faces_bord_;
  s << mes_faces_joint_;
  s << mes_faces_raccord_;
  s << mes_bords_int_;
  if (nb_groupes_faces() !=0)
    {
      s << finl << "groupes_faces" << finl;
      s << mes_groupes_faces_;
    }

  // New in TRUST 1.9.8 - list of periodic boundaries are directly stored in Domain class
  s << finl << "bords_perio" << finl;
  s << bords_perio_;

  s << "}" << finl;

  return s;
}

/*! @brief See readOn_has_perio()
 */
template<typename _SZ_>
Entree& Domaine_32_64<_SZ_>::readOn(Entree& s)
{
  bool dnu;
  return readOn_has_perio(s, dnu);
}

/*! @brief Reads the objects constituting a Domain from an input stream.
 *
 * Once the objects are read they are associated to the domain.
 *
 * @param (Entree& s) an input stream
 * @param (bool& has_perio) set to True if periodic boundaries were read, false otherwise.
 * @return (Entree&) the modified input stream
 */
template<typename _SZ_>
Entree& Domaine_32_64<_SZ_>::readOn_has_perio(Entree& s, bool& has_perio)
{
#ifdef SORT_POUR_DEBOG
  s.setf(ios::scientific);
  s.precision(20);
#endif
  has_perio = false;
  // BM addition: reset the structure (this has the effect of unlocking the parallel structure)
  sommets_.reset();
  renum_som_perio_.reset();
  // do not reset the name (already read)
  // for deformable I don't know...

  Nom tmp;
  s >> tmp;
  // If the domain is not yet named, use the name that was read
  if (nom_=="??") nom_=tmp;
  Cerr << "Reading domain " << le_nom() << finl;
  s >> sommets_;
  // PL : not quite exact the number of nodes displayed, joint nodes are counted multiple times...
  trustIdType nbsom = mp_sum(sommets_.dimension(0));
  Cerr << " Number of nodes: " << nbsom << finl;

  // Reading element description (what was fomerly the "domaine" part) - this used to be a list so check for '{ }'
  Nom acc;
  s >> acc;
  assert (acc == "{");
  read_former_domaine(s, has_perio);
  check_domaine();

  if (Process::is_sequential() && (NettoieNoeuds_32_64<_SZ_>::NettoiePasNoeuds==0) )
    {
      NettoieNoeuds_32_64<_SZ_>::nettoie(*this);
      nbsom = mp_sum(sommets_.dimension(0));
      Cerr << " Number of nodes after node-cleanup: " << nbsom << finl;
    }

  // Initialize the "sequential" descriptors (warning: this blocks the resize of the vertex and element arrays!)
  Scatter::init_sequential_domain(*this);
  check_domaine();
  return s;
}


/*! @brief read what was (before TRUST 1.9.2) the "domaine" part from the input stream
 * i.e. (roughly) the element description.
 *
 * @param read_perio set to True if periodic boundaries were read in the domain (from TRUST 1.9.8)
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::read_former_domaine(Entree& s, bool& read_perio)
{
  Nom dnu, acc;
  read_perio = false;
  Cerr << " Reading part of domain " << le_nom() << finl;
  s >> dnu; // Name of the Domaine, now unused ...
  s >> elem_;
  mes_elems_.reset();
  s >> mes_elems_;
  mes_faces_bord_.vide();
  s >> mes_faces_bord_;
  mes_faces_joint_.vide();
  s >> mes_faces_joint_;
  mes_faces_raccord_.vide();
  s >> mes_faces_raccord_;
  mes_bords_int_.vide();
  s >> mes_bords_int_;
  mes_groupes_faces_.vide();
  s >> acc;
  if (acc == "groupes_faces")
    {
      s >> mes_groupes_faces_;
      s >> acc;
    }
  // Tries to read list of periodic boundaries:
  bords_perio_.clear();
  if (acc == "bords_perio")
    {
      s >> bords_perio_;
      s >> acc;
      read_perio = true;
    }
  if (acc != "}")
    Process::exit( "misformatted domain file : One expected a closing bracket } to end. ");
}

template<class LIST_FRONTIERE>
void check_frontiere(const LIST_FRONTIERE& list, const char *msg)
{
  int n = list.size();
  if (!is_parallel_object(n))
    {
      Cerr << " Fatal error: processors don't have the same number of boundaries " << msg << finl;
      Process::exit();
    }
  for (int i = 0; i < n; i++)
    {
      const Nom& nom = list[i].le_nom();
      Cerr << "  Boundary " << msg << " : " << nom << finl;
      if (!is_parallel_object(nom))
        {
          Cerr << " Fatal error: processors don't have the same number of boundaries " << msg << finl;
          Process::exit();
        }
    }
}

template<class LIST_FRONTIERE>
void check_frontiere_own_ptr(const LIST_FRONTIERE& list, const char *msg)
{
  int n = list.size();
  if (!is_parallel_object(n))
    {
      Cerr << " Fatal error: processors don't have the same number of boundaries " << msg << finl;
      Process::exit();
    }
  for (int i = 0; i < n; i++)
    {
      const Nom& nom = list[i]->le_nom();
      Cerr << "  Boundary " << msg << " : " << nom << finl;
      if (!is_parallel_object(nom))
        {
          Cerr << " Fatal error: processors don't have the same number of boundaries " << msg << finl;
          Process::exit();
        }
    }
}


/*! @brief associate the read objects to the domaine and check that the reading objects are coherent
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::check_domaine()
{
  // replace Type_Face::vide_0D with the correct type for processors that have no boundary faces:
  {
    int i;
    int n = nb_front_Cl();
    for (i = 0; i < n; i++)
      ::corriger_type(frontiere(i).faces(), type_elem());
  }

  if (mes_faces_bord_.size() == 0 && mes_faces_raccord_.size() == 0 && Process::is_sequential())
    Cerr << "Warning, the reread domaine " << nom_ << " has no defined boundaries (none boundary or connector)." << finl;

  mes_faces_bord_.associer_domaine(*this);
  mes_faces_joint_.associer_domaine(*this);
  mes_faces_raccord_.associer_domaine(*this);
  mes_bords_int_.associer_domaine(*this);
  mes_groupes_faces_.associer_domaine(*this);
  elem_->associer_domaine(*this);
  fixer_premieres_faces_frontiere();

  const trustIdType nb_elem = mp_sum(mes_elems_.dimension(0));
  Cerr << "  Number of elements: " << nb_elem << finl;

  // Sanity checks:
  // All processors must have the same number of boundaries and the same names
  ::check_frontiere(mes_faces_bord_, "(Bord)");
  ::check_frontiere_own_ptr(mes_faces_raccord_, "(Raccord)");
  ::check_frontiere(mes_bords_int_, "(Bord_Interne)");
  ::check_frontiere(mes_groupes_faces_, "(Groupe_Faces)");
}

/*! @brief Searches the indices of elements containing the vertices specified by the "sommets" parameter.
 *
 *     Uses:
 *      ArrOfInt_t& Domaine_32_64<_SZ_>::chercher_elements(const DoubleTab&,ArrOfInt_t&) const
 *
 * @param (IntTab& sommets) the array of vertex indices whose containing elements are searched
 * @param (ArrOfInt_t& elem_) the array containing the indices of elements containing the specified vertices
 * @return (ArrOfInt_t&) the array of vertex indices whose containing elements are searched
 */
template<typename _SZ_>
typename Domaine_32_64<_SZ_>::SmallArrOfTID_t& Domaine_32_64<_SZ_>::indice_elements(const IntTab& sommets, SmallArrOfTID_t& elem, int reel) const
{
  int i, j, k;
  const DoubleTab_t& les_coord = sommets_;
  int sz_sommets = sommets.dimension(0);
  DoubleTab xg(sz_sommets, Objet_U::dimension);
  for (i = 0; i < sz_sommets; i++)
    for (j = 0; j < nb_som_elem(); j++)
      for (k = 0; k < Objet_U::dimension; k++)
        xg(i, k) += les_coord(sommets(i, j), k);

  xg /= (double)nb_som_elem();
  return chercher_elements(xg, elem, reel);
}

/*! @brief Searches for the elements containing the points whose coordinates are specified.
 *
 * @param (DoubleTab& positions) the coordinates of the points whose containing element is sought
 * @param (ArrOfInt_t& elements) the array of indices of the elements containing the specified points
 * @return (ArrOfInt_t&) the array of indices of the elements containing the specified points
 */
template<typename _SZ_>
typename Domaine_32_64<_SZ_>::SmallArrOfTID_t& Domaine_32_64<_SZ_>::chercher_elements(const DoubleTab& positions, SmallArrOfTID_t& elements, int reel) const
{
  bool set_cache = false;
  // PL: We should call chercher_elements(x,y,z,elem) if positions.dimension(0)=1 ...
  if (!deformable() && positions.dimension(0) > 1)
    {
      set_cache = true;
      if (!deriv_octree_ || !deriv_octree_->construit())
        {
          // Flush the cache
          cached_elements_.reset();
          cached_positions_.reset();
        }
      else
        {
          // Search in the cache:
          for (int i = 0; i < cached_positions_.size(); i++)
            if (sameDoubleTab(positions, cached_positions_[i]))
              {
                int size = cached_positions_[i].dimension(0);
                if (elements.size_array() != size)
                  elements.resize_tab(size);
                elements = cached_elements_[i];
                // elements.ref_array(cached_elements_[i]); // No - triggers an assert (ex Sondes.data) and also in parallel, elements is modified in probes....
                return elements;
              }
        }
    }
  const OctreeRoot_t& octree = construit_octree(reel);
  int sz = positions.dimension(0);
  const int dim = positions.dimension_int(1);
  // resize_tab is virtual; if it is a Vect or a Tab it calls the resize
  // method of the derived class:
  elements.resize_tab(sz, RESIZE_OPTIONS::NOCOPY_NOINIT);
  double y = 0, z = 0;
  for (int i = 0; i < sz; i++)
    {
      double x = positions(i, 0);
      if (dim > 1)
        y = positions(i, 1);
      if (dim > 2)
        z = positions(i, 2);
      elements[i] = octree.rang_elem(x, y, z);
    }
  if (set_cache)
    {
      // Safety measure: observed on an FT calculation (cache growing indefinitely, variable number of particles...)
      // if (cached_memory>1e8) // 100Mo/proc
      if (cached_positions_.size()>100) // Change heuristic cause 100Mo on GPU is tiny !
        {
          // Flush the cache
          Cerr << "Warning, cache flushed in Domaine_32_64<_SZ_>::chercher_elements() cause too much lines used !" << finl;
          cached_elements_.reset();
          cached_positions_.reset();
          cached_memory = 0;
        }
      else
        {
          // Store in cache
          cached_positions_.add(positions);
          cached_elements_.add(elements);
          // Send cached arrays to device:
          int last = cached_positions_.size();
          mapToDevice(cached_positions_[last-1]);
          mapToDevice(cached_elements_[last-1]);
          cached_memory += (double)(positions.size_array() * sizeof(double));
          cached_memory += (double)(elements.size_array() * sizeof(int));
          if (cached_memory > 1e7)   // 10Mo
            {
              Cerr << 2 * cached_positions_.size() << " arrays cached in memory for Domaine_32_64<_SZ_>::chercher_elements(...): ";
              if (cached_memory < 1e6)
                Cerr << int(cached_memory / 1024) << " KBytes" << finl;
              else
                Cerr << int(cached_memory / 1024 / 1024) << " MBytes" << finl;
            }
        }
    }
  return elements;
}

/*! @brief Searches for the elements containing the points whose coordinates are specified.
 *
 * @param (DoubleVect_t<_SZ_>& positions) the coordinates of the point whose containing element is sought
 * @param (ArrOfInt_t& elements) the array of indices of the elements containing the specified points
 * @return (ArrOfInt_t&) the array of indices of the elements containing the specified points
 */
template<typename _SZ_>
typename Domaine_32_64<_SZ_>::SmallArrOfTID_t& Domaine_32_64<_SZ_>::chercher_elements(const DoubleVect& positions, SmallArrOfTID_t& elements, int reel) const
{
  int n = positions.size();
  if (n != dimension)
    {
      Cerr << "Domaine_32_64::chercher_elements(const DoubleVect& positions, ArrOfInt& elements, int reel) const -> Coding is made to copy a doublevect(dimesnion) in a DoubleTab(1,dimension)" << finl;
      Cerr << "But, it comes with a DoubleVect of size " << n << " instead of " << dimension << finl;
      assert(0);
      Process::exit();
    }
  DoubleTab positions2(1, n);
  for (int ii = 0; ii < n; ii++)
    positions2(0, ii) = positions(ii);
  return chercher_elements(positions2, elements, reel);
}


/*! @brief Returns -1 if face is not an internal boundary face, or the index of the duplicated face otherwise.
 *
 * @param (int face) the index of the internal boundary face to search for
 * @return (int) -1 if the specified face is not an internal boundary face, or the index of the duplicated face otherwise
 * @throws TRUST error (face not found)
 */
template<typename _SZ_>
typename Domaine_32_64<_SZ_>::int_t Domaine_32_64<_SZ_>::face_bords_interne_conjuguee(int_t face) const
{
  if ((face) >= nb_faces_frontiere())
    return -1;
  int_t compteur = nb_faces_frontiere() - nb_faces_bords_int();
  if ((face) < compteur)
    return -1;

  for (const auto& itr : mes_bords_int_)
    {
      const Faces_32_64<_SZ_>& les_faces = itr.faces();
      int_t nbf = les_faces.nb_faces();
      if (face < nbf + compteur)
        {
          nbf /= 2;
          if ((face - compteur) < nbf)
            return face + nbf;
          else
            return face - nbf;
        }
      compteur += (2 * nbf);
    }

  Cerr << "TRUST error in Domaine_32_64<_SZ_>::face_bords_interne_conjuguee " << finl;
  Process::exit();
  return -1;
}

/*! @brief Merges boundaries with the same name for: boundaries, periodic boundaries, internal boundaries and face groups.
 */
template<typename _SZ_>
int Domaine_32_64<_SZ_>::comprimer()
{
  {
    // Boundaries
    auto& list = mes_faces_bord_.get_stl_list();

    // first loop over list elements
    for (auto it = list.begin(); it != list.end(); ++it)
      {
        Frontiere_t& front = *it;
        front.associer_domaine(*this); // In case the boundary's domain is not the correct one
        Journal() << "Domaine_32_64<_SZ_>::comprimer() bord : " << front.le_nom() << finl;

        // second loop over list elements, starting from an incremented position
        for (auto it2 = std::next(it); it2 != list.end(); )
          {
            Frontiere_t& front2 = *it2;
            if (front.le_nom() == front2.le_nom())
              {
                Journal() << "Merging boundary: " << front.le_nom() << finl;
                front.add(front2);
                it2 = list.erase(it2);
              }
            else
              ++it2;

            Journal() << front.le_nom() << " is associated with: " << front.domaine().le_nom() << finl;
          }
      }
  }

  {
    // Internal boundaries:
    auto& list = mes_bords_int_.get_stl_list();
    for (auto it = list.begin(); it != list.end(); ++it)
      {
        Frontiere_t& front = *it;
        for (auto it2 = std::next(it); it2 != list.end(); )
          {
            Frontiere_t& front2 = *it2;
            if (front.le_nom() == front2.le_nom())
              {
                front.add(front2);
                it2 = list.erase(it2);
              }
            else
              ++it2;
          }
      }
  }

  {
    // Face groups:
    auto& list = mes_groupes_faces_.get_stl_list();
    for (auto it = list.begin(); it != list.end(); ++it)
      {
        Frontiere_t& front = *it;
        for (auto it2 = std::next(it); it2 != list.end(); )
          {
            Frontiere_t& front2 = *it2;
            if (front.le_nom() == front2.le_nom())
              {
                front.add(front2);
                it2 = list.erase(it2);
              }
            else
              ++it2;
          }
      }
  }

  {
    // Connections
    auto& list = mes_faces_raccord_.get_stl_list();
    for (auto it = list.begin(); it != list.end(); ++it)
      {
        Frontiere_t& front = (*it).valeur();
        Journal() << "Raccord : " << front.le_nom() << finl;
        for (auto it2 = std::next(it); it2 != list.end(); )
          {
            Frontiere_t& front2 = (*it2).valeur();
            if (front.le_nom() == front2.le_nom())
              {
                front.add(front2);
                it2 = list.erase(it2);
              }
            else
              ++it2;
          }
      }
  }
  return 1;
}

/*! @brief Returns the index of the element containing the point whose coordinates are specified.
 *
 * @param (double x) X coordinate
 * @param (double y) Y coordinate
 * @param (double z) Z coordinate
 * @return (int) the index of the element containing the point whose coordinates are specified.
 */
template<typename _SZ_>
typename Domaine_32_64<_SZ_>::int_t Domaine_32_64<_SZ_>::chercher_elements(double x, double y, double z, int reel) const
{

  const OctreeRoot_t& octree = construit_octree(reel);
  return octree.rang_elem(x, y, z);
}

/*! @brief
 *
 * @param (DoubleTab& pos)
 * @param (ArrOfInt_t& som)
 * @return (ArrOfInt_t&)
 */
template<typename _SZ_>
typename Domaine_32_64<_SZ_>::SmallArrOfTID_t& Domaine_32_64<_SZ_>::chercher_sommets(const DoubleTab& pos, SmallArrOfTID_t& som, int reel) const
{
  const OctreeRoot_t& octree = construit_octree(reel);
  octree.rang_sommet(pos, som);
  return som;
}

/*! @brief
 *
 * @param (DoubleTab& pos)
 * @param (IntTab& aretes_som) the definition of edges by their vertices
 * @return (ArrOfInt_t& aretes) list of edges found
 */
template<typename _SZ_>
typename Domaine_32_64<_SZ_>::SmallArrOfTID_t& Domaine_32_64<_SZ_>::chercher_aretes(const DoubleTab& pos, SmallArrOfTID_t& aretes, int reel) const
{
  const OctreeRoot_t& octree = construit_octree(reel);
  octree.rang_arete(pos, aretes);
  return aretes;
}

/*! @brief
 *
 * @param (double x) X coordinate
 * @param (double y) Y coordinate
 * @param (double z) Z coordinate
 */
template<typename _SZ_>
typename Domaine_32_64<_SZ_>::int_t Domaine_32_64<_SZ_>::chercher_sommets(double x, double y, double z, int reel) const
{
  const OctreeRoot_t& octree = construit_octree(reel);
  return octree.rang_sommet(x, y, z);
}

/*! Builds the elem_virt_pe_num_ array from the mes_elems array
* (using the distant and virtual spaces of mes_elems).
* Non-optimal memory algorithm: mes_elems is duplicated whereas only
* a two-column array is needed.
* See Domaine.h: elem_virt_pe_num_
*/
template<typename _SZ_>
void Domaine_32_64<_SZ_>::construire_elem_virt_pe_num()
{
  this->construire_elem_virt_pe_num(elem_virt_pe_num_);
}

template<typename _SZ_>
void Domaine_32_64<_SZ_>::construire_elem_virt_pe_num(IntTab_t& elem_virt_pe_num_cpy) const
{
  IntTab_t tableau_echange(mes_elems_);
  assert(tableau_echange.dimension(1) >= 2);
  const int_t n = nb_elem();
  const int_t n_virt = nb_elem_tot() - n;
  const int moi = me();
  for (int_t i = 0; i < n; i++)
    {
      tableau_echange(i, 0) = moi;
      tableau_echange(i, 1) = i;
    }
  tableau_echange.echange_espace_virtuel();

  elem_virt_pe_num_cpy.resize(n_virt, 2);
  for (int_t i = 0; i < n_virt; i++)
    {
      elem_virt_pe_num_cpy(i, 0) = tableau_echange(n + i, 0);
      elem_virt_pe_num_cpy(i, 1) = tableau_echange(n + i, 1);
    }
}


/*! @brief Computes the center of gravity of the domain
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::calculer_mon_centre_de_gravite(ArrOfDouble& c)
{
  c = 0;
  // Volumes computed cause stored in Domaine_VF and so not available in Domaine...
  DoubleVect_t volumes;
  DoubleVect_t inverse_volumes;
  calculer_volumes(volumes, inverse_volumes);
  DoubleTab_t xp;
  calculer_centres_gravite(xp);
  double volume = 0;
  for (int_t i = 0; i < nb_elem(); i++)
    for (int j = 0; j < dimension; j++)
      {
        c[j] += xp(i, j) * volumes(i);
        volume += volumes(i);
      }
  // Case of an empty Domain:
  if (volume > 0)
    c /= volume;
  cg_moments_ = c;
  volume_total_ = mp_somme_vect(volumes);
}

/*! @brief Computes the volumes of the domain elements.
 *
 * @param (DoubleVect& volumes) the array containing the volumes of the domain elements
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::calculer_volumes(DoubleVect_t& volumes, DoubleVect_t& inverse_volumes) const
{
  if (!volumes.get_md_vector())
    creer_tableau_elements(volumes, RESIZE_OPTIONS::NOCOPY_NOINIT);
  elem_->calculer_volumes(volumes); // Sizes and computes the DoubleVect volumes
  // Check and fill inverse_volumes
  if (!inverse_volumes.get_md_vector())
    creer_tableau_elements(inverse_volumes, RESIZE_OPTIONS::NOCOPY_NOINIT);
  int_t size = volumes.size_totale();
  for (int_t i = 0; i < size; i++)
    {
      double v = volumes(i);
      if (v <= 0.)
        {
          Cerr << "Volume[" << i << "]=" << v << finl;
          Cerr << "Several volumes of the mesh are not positive." << finl;
          Cerr << "Something is wrong in the mesh..." << finl;
          Process::exit();
        }
      inverse_volumes(i) = 1. / v;
    }
}

/*! @brief Computes the centers of gravity of the domain edges.
 *
 * @param (DoubleTab& xa) the array containing the centers of gravity of the domain edges
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::calculer_centres_gravite_aretes(DoubleTab_t& xa) const
{
  const DoubleTab_t& coord = sommets_;
  // Computes the centers of gravity of real edges only
  xa.resize(nb_aretes(), dimension);
  for (int_t i = 0; i < nb_aretes(); i++)
    for (int j = 0; j < dimension; j++)
      xa(i, j) = 0.5 * (coord(aretes_som_(i, 0), j) + coord(aretes_som_(i, 1), j));
}

template<typename _SZ_>
void Domaine_32_64<_SZ_>::rang_elems_sommet(SmallArrOfTID_t& elems, double x, double y, double z) const
{
  const OctreeRoot_t& octree = construit_octree();
  octree.rang_elems_sommet(elems, x, y, z);
}

template<typename _SZ_>
void Domaine_32_64<_SZ_>::invalide_octree()
{
  if (deriv_octree_)
    deriv_octree_.detach();
}

template<typename _SZ_>
const typename Domaine_32_64<_SZ_>::OctreeRoot_t& Domaine_32_64<_SZ_>::construit_octree() const
{
  if (!deriv_octree_)
    {
      if constexpr (sizeof(_SZ_) == 4)
        deriv_octree_.typer("OctreeRoot");
      else if constexpr (sizeof(_SZ_) == 8)
        deriv_octree_.typer("OctreeRoot_64");
    }
  OctreeRoot_t& octree = deriv_octree_.valeur();
  if (!octree.construit())
    {
      octree.associer_Domaine(*this);
      octree.construire();
    }
  return octree;
}

/*! @brief Build the octree if not already done
 */
template<typename _SZ_>
const typename Domaine_32_64<_SZ_>::OctreeRoot_t& Domaine_32_64<_SZ_>::construit_octree(int& reel) const
{
  if (!deriv_octree_)
    {
      if constexpr (sizeof(_SZ_) == 4)
        deriv_octree_.typer("OctreeRoot");
      else if constexpr (sizeof(_SZ_) == 8)
        deriv_octree_.typer("OctreeRoot_64");
    }
  OctreeRoot_t& octree = deriv_octree_.valeur();
  if (!octree.construit() || (reel != octree.reel()))
    {
      octree.associer_Domaine(*this);
      octree.construire(reel);
    }
  return octree;
}

/*! @brief Creates a parallel array of values at elements.
 *
 * See MD_Vector_tools::creer_tableau_distribue()
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::creer_tableau_elements(Array_base& x, RESIZE_OPTIONS opt) const
{
  const MD_Vector& md = md_vector_elements();
  MD_Vector_tools::creer_tableau_distribue(md, x, opt);
}

/*! @brief Returns the parallel descriptor of element arrays of the domain
 */
template<typename _SZ_>
const MD_Vector& Domaine_32_64<_SZ_>::md_vector_elements() const
{
  const MD_Vector& md = mes_elems_.get_md_vector();
  if (!md)
    {
      Cerr << "Internal error in Domaine_32_64<_SZ_>::md_vector_elements(): descriptor for elements not initialized\n"
           << " You might use a buggy Domain constructor that does not build descriptors,\n"
           << " Use the following syntax to finish the domain construction\n"
           << "  Scatter ; " << le_nom() << finl;
      Process::exit();
    }
  // For now the descriptor is taken from the mes_elems array, but we could
  // store a copy in the domain if that would be useful...
  return md;
}

template<typename _SZ_>
double Domaine_32_64<_SZ_>::volume_total() const
{
  assert(volume_total_ >= 0.); // Not computed yet ???
  return volume_total_;
}

template<typename _SZ_>
DoubleTab Domaine_32_64<_SZ_>::getBoundingBox() const
{
  DoubleTab BB(dimension, 2);
  int_t nbsom=sommets_.dimension(0);
  for (int j=0; j<dimension; j++)
    {
      double min_=0.5*DMAXFLOAT;
      double max_=-0.5*DMAXFLOAT;
      for (int_t i=0; i<nbsom; i++)
        {
          double c = sommets_(i,j);
          min_ = (c < min_ ? c : min_);
          max_ = (c > max_ ? c : max_);
        }
      BB(j,0) = min_;
      BB(j,1) = max_;
    }
  return BB;
}

/*! @brief Adds nodes (or vertices) to the domain (without checking for duplicates)
 *
 * @param (DoubleTab& soms) the array containing the coordinates of the nodes to add to the domain
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::ajouter(const DoubleTab_t& soms)
{
  int_t oldsz=sommets_.dimension(0);
  int_t ajoutsz=soms.dimension(0);
  int dim = soms.dimension_int(1);
  sommets_.resize(oldsz+ajoutsz,dim);
  for(int_t i=0; i<ajoutsz; i++)
    for(int k=0; k<dim; k++)
      sommets_(oldsz+i,k)=soms(i,k) ;
}

/*! @brief Adds nodes to the domain with elimination of duplicate nodes. On return, nums contains the new indices of the nodes from soms
 *
 *     after elimination of duplicates.
 *
 * @param (DoubleTab& soms) the array containing the coordinates of the nodes to add to the domain
 * @param (IntVect& nums) the array of new indices after adding new nodes and eliminating duplicates.
 * @throws duplicate nodes were found
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::ajouter(const DoubleTab_t& soms, IntVect_t& nums)
{
  int_t oldsz = sommets_.dimension(0);
  int_t ajoutsz = soms.dimension(0);
  int dim = soms.dimension_int(1);
  nums.resize(ajoutsz);
  nums=-1;
  if(oldsz!=0)
    {
      assert(dim==sommets_.dimension(1));
      Octree_Double_32_64<_SZ_> octree;
      octree.build_nodes(les_sommets(), 0 /* do not include virtual vertices */);

      int compteur=0;
      ArrOfDouble tab_coord(dim);
      ArrOfInt_t liste_sommets;
      for(int_t i=0; i< ajoutsz; i++)
        {
          for (int j = 0; j < dim; j++)
            tab_coord[j] = soms(i,j);
          octree.search_elements_box(tab_coord, epsilon_, liste_sommets);
          octree.search_nodes_close_to(tab_coord, les_sommets(), liste_sommets, epsilon_);
          const int_t nb_sommets_proches = liste_sommets.size_array();
          if (nb_sommets_proches == 0)
            {
              // No vertex of the first domain is close to vertex i.
              // Keep i.
            }
          else if (nb_sommets_proches == 1)
            {
              // One vertex coincides with vertex i within epsilon_.
              // Do not keep the vertex
              nums(i) = liste_sommets[0];
              compteur++;
            }
          else
            {
              // Several vertices of the initial domain are within radius epsilon.
              // epsilon is too large.
              Cerr << "Error : several nodes of the domain 1 are within radius epsilon="
                   << epsilon_ << " of point " << tab_coord << ". We must reduce epsilon. " << finl;
              Process::exit();
            }
        }
      Cerr << compteur << " double nodes were found \n";
      sommets_.resize(oldsz+ajoutsz-compteur,dim);
      compteur=0;
      for(int_t i =0; i<ajoutsz; i++)
        if(nums(i)==-1)
          {
            nums(i)=oldsz+compteur;
            compteur++;
            for(int k=0; k<dim; k++)
              sommets_(nums(i),k)=soms(i,k) ;
          }
    }
  else
    {
      sommets_=soms;
      // if som has a descriptor, delete it:
      sommets_.set_md_vector(MD_Vector());
      for(int_t i=0; i<ajoutsz; i++)
        nums(i)=i;
    }
}

/*! @brief Creates an array with one "row" per mesh vertex.
 *
 * See MD_Vector_tools::creer_tableau_distribue()
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::creer_tableau_sommets(Array_base& v, RESIZE_OPTIONS opt) const
{
  const MD_Vector& md = md_vector_sommets();
  MD_Vector_tools::creer_tableau_distribue(md, v, opt);
}


/*! @brief only read vertices from the stream s
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::read_vertices(Entree& s)
{
  // BM addition: reset the structure (this has the effect of unlocking the parallel structure)
  sommets_.reset();
  renum_som_perio_.reset();

  Nom tmp;
  s >> tmp;
  // If the domain is not yet named, use the name that was read
  if (nom_=="??") nom_=tmp;
  Cerr << "Reading vertices for domain " << le_nom() << finl;
  s >> sommets_;
}


/*! @brief Writes the boundary names to an output stream.
 *
 * Writes the names of: boundaries, periodic boundaries, connections and face groups.
 *
 * @param (Sortie& os) an output stream
 */
template <typename _SZ_>
void Domaine_32_64<_SZ_>::ecrire_noms_bords(Sortie& os) const
{
  // Boundaries
  for (const auto &itr : mes_faces_bord_)
    os << itr.le_nom() << finl;

  // Connections:
  for (const auto &itr : mes_faces_raccord_)
    os << itr->le_nom() << finl;

  // Internal boundaries:
  for (const auto &itr : mes_bords_int_)
    os << itr.le_nom() << finl;

  // Face groups:
  for (const auto &itr : mes_groupes_faces_)
    os << itr.le_nom() << finl;
}

template <typename _SZ_>
int Domaine_32_64<_SZ_>::rang_frontiere(const Nom& un_nom) const
{
  int i = 0;
  for (const auto &itr : mes_faces_bord_)
    {
      if (itr.le_nom() == un_nom)
        return i;
      ++i;
    }

  for (const auto &itr : mes_faces_raccord_)
    {
      if (itr->le_nom() == un_nom)
        return i;
      ++i;
    }

  for (const auto &itr : mes_bords_int_)
    {
      if (itr.le_nom() == un_nom)
        return i;
      ++i;
    }

  for (const auto &itr : mes_groupes_faces_)
    {
      if (itr.le_nom() == un_nom)
        return i;
      ++i;
    }
  Cerr << "Domaine_32_64<_SZ_>::rang_frontiere(): We have not found a boundary with name " << un_nom << finl;
  Process::exit();
  return -1;
}

template <typename _SZ_>
const typename Domaine_32_64<_SZ_>::Frontiere_t& Domaine_32_64<_SZ_>::frontiere(const Nom& un_nom) const
{
  int i = rang_frontiere(un_nom);
  return frontiere(i);
}

template <typename _SZ_>
typename Domaine_32_64<_SZ_>::Frontiere_t& Domaine_32_64<_SZ_>::frontiere(const Nom& un_nom)
{
  int i = rang_frontiere(un_nom);
  return frontiere(i);
}

template <typename _SZ_>
void Domaine_32_64<_SZ_>::fixer_premieres_faces_frontiere()
{
  Journal() << "Domaine_32_64<_SZ_>::fixer_premieres_faces_frontiere()" << finl;
  int_t compteur = 0;
  for (auto &itr : mes_faces_bord_)
    {
      itr.fixer_num_premiere_face(compteur);
      compteur += itr.nb_faces();
      Journal() << "Boundary " << itr.le_nom() << " starts at face: " << itr.num_premiere_face() << finl;
    }
  for (auto &itr : mes_faces_raccord_)
    {
      itr->fixer_num_premiere_face(compteur);
      compteur += itr->nb_faces();
      Journal() << "Connection " << itr->le_nom() << " starts at face: " << itr->num_premiere_face() << finl;
    }
  if (std::is_same<_SZ_, int>::value)
    for (auto &itr : mes_faces_joint_)
      {
        itr.fixer_num_premiere_face((int)compteur);
        compteur += itr.nb_faces();
        Journal() << "Joint " << itr.le_nom() << " starts at face: " << itr.num_premiere_face() << finl;
      }
  for (auto &itr : mes_groupes_faces_)
    itr.fixer_num_premiere_face(-1);
}


template<typename _SZ_>
template<typename _BORD_TYP_>
void Domaine_32_64<_SZ_>::correct_type_single_border_type(std::list<_BORD_TYP_>& list)
{
  // first loop over list elements
  for (auto it = list.begin(); it != list.end(); ++it)
    {
      Frontiere_t& front = *it;
      if (front.faces().type_face() == Type_Face::vide_0D)
        {
          // second loop over list elements, starting from an incremented position
          for (auto it2 = std::next(it); it2 != list.end();)
            {
              Frontiere_t& front2 = *it2;
              if (front.le_nom() == front2.le_nom())
                {
                  front.faces().typer(front2.faces().type_face());
                  break;
                }
              else
                ++it2;
            }
        }
    }
}

/*! @brief Correcting type of borders if they were empty before merge (ie equal to vide_0D)
 *
 * Difference with corriger_type is that we don't want to delete faces inside borders afterwards.
 * Int version handles joints, not the 64b one.
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::correct_type_of_borders_after_merge()
{
  correct_type_single_border_type(mes_faces_bord_.get_stl_list());
  correct_type_single_border_type(mes_bords_int_.get_stl_list());
  correct_type_single_border_type(mes_faces_raccord_.get_stl_list());
  correct_type_single_border_type(mes_groupes_faces_.get_stl_list());
}


template<>
void Domaine_32_64<int>::correct_type_of_borders_after_merge()
{
  correct_type_single_border_type(mes_faces_bord_.get_stl_list());
  correct_type_single_border_type(mes_bords_int_.get_stl_list());
  correct_type_single_border_type(mes_faces_raccord_.get_stl_list());
  correct_type_single_border_type(mes_groupes_faces_.get_stl_list());
  // The joints - only for 32b Domaine:
  correct_type_single_border_type(mes_faces_joint_.get_stl_list());
}


template<typename _SZ_>
void Domaine_32_64<_SZ_>::imprimer() const
{
  Cerr << "==============================================" << finl;
  Cerr << "The extreme coordinates of the domain " << le_nom() << " are:" << finl;
  // There is no min/max search method in DoubleTab so it is coded here:
  DoubleTab BB = getBoundingBox();
  ArrOfDouble bb_min(dimension), bb_max(dimension);
  for (int j=0; j<dimension; j++)
    {
      bb_min[j] = BB(j,0);
      bb_max[j] = BB(j,1);
    }
  Process::mp_min_for_each_item(bb_min);
  Process::mp_max_for_each_item(bb_max);
  for (int j=0; j<dimension; j++)
    {
      if (j==0) Cerr << "x ";
      if (j==1) Cerr << "y ";
      if (j==2) Cerr << "z ";
      Cerr << "is between " << bb_min[j] << " and " << bb_max[j] << finl;
    }
  Cerr << "==============================================" << finl;
  // We recompute volumes (cause stored in Domaine_VF and so not available from Domaine...):
  DoubleVect_t volumes;
  DoubleVect_t inverse_volumes;
  calculer_volumes(volumes,inverse_volumes);
  Cerr << "==============================================" << finl;
  Cerr << "The volume cells of the domain " << le_nom() << " are:" << finl;
  const int_t i_vmax = imax_array(volumes);
  const int_t i_vmin = imin_array(volumes);
  const double vmin_local = (i_vmin < 0) ? 1e40 : volumes[i_vmin];
  const double vmax_local = (i_vmax < 0) ? -1e40 : volumes[i_vmax];
  const double volmin = mp_min(vmin_local);
  const double volmax = mp_max(vmax_local);
  double volume_total = mp_somme_vect(volumes);
  const int_t nbe = nb_elem();
  double volmoy = volume_total / Process::mp_sum_as_double(nbe);
  Cerr << "sum(volume cells)= "  << volume_total << finl;
  Cerr << "mean(volume cells)= " << volmoy << finl;
  Cerr << "min(volume cells)= "  << volmin << finl;
  Cerr << "max(volume cells)= "  << volmax << finl;
  if (volmin*1000<volmoy)
    {
      Cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
      Cerr << "Warning, a cell volume is more than 1000 times smaller than the average cell volume. Check your mesh." << finl;
      Cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
    }
  Cerr << "==============================================" << finl;
}

/*! @brief Merge another Domaine into this, without considering vertices which are handled separately
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::merge_wo_vertices_with(Domaine_32_64<_SZ_>& dom2)
{
  Cerr << "   Merging elem info for domain "<< nom_ << " with " << dom2.nom_ << finl;

  // Prepare type if first merge:
  if (!elem_)
    elem_ = dom2.elem_;

  // Prepare correct initial elem size if first merge
  if (nb_elem() == 0)
    les_elems().resize(0, dom2.les_elems().dimension_int(1));

  int_t sz1 = les_elems().dimension(0);
  int_t sz2 = dom2.les_elems().dimension(0);
  int nb_ccord = les_elems().dimension_int(1);
  IntTab_t& elems1 = les_elems();
  IntTab_t& elems2 = dom2.les_elems();
  elems1.resize(sz1+sz2, nb_ccord);
  for(int_t i=0; i<sz2; i++)
    for(int j=0; j<nb_ccord; j++)
      elems1(sz1+i,j)=elems2(i,j);

  dom2.faces_bord().associer_domaine(*this);
  faces_bord().add(dom2.faces_bord());

  // Take care of the joints only in 32 bits instance when this is called by Scatter
  dom2.faces_joint().associer_domaine(*this);
  faces_joint().add(dom2.faces_joint());

  dom2.faces_raccord().associer_domaine(*this);
  faces_raccord().add(dom2.faces_raccord());

  dom2.bords_int().associer_domaine(*this);
  bords_int().add(dom2.bords_int());

  dom2.groupes_faces().associer_domaine(*this);
  groupes_faces().add(dom2.groupes_faces());

  // Add periodic boundary names if not already there:
  auto& bp = bords_perio();
  const auto& bp2 = dom2.bords_perio();
  for (const auto &b : bp2)
    {
      if (bp.rang(b) < 0)
        bp.add(b);
    }

  correct_type_of_borders_after_merge();
  comprimer();
  comprimer_joints();
  invalide_octree();
}

/*! @brief Associates a Sous_Domaine to the Domain.
 *
 * The interface accepts any Objet_U but only handles (dynamically) the
 *     association of an object derived from Sous_Domaine.
 *
 * @param (Objet_U& ob) the object to associate
 * @return (int) 1 if the association succeeded, 0 otherwise (the object was not derived from Sous_Domaine)
 */
template<typename _SZ_>
int Domaine_32_64<_SZ_>::associer_(Objet_U& ob)
{
  if( sub_type(Sous_Domaine_t, ob))
    {
      add(ref_cast(Sous_Domaine_t, ob));
      ob.associer_(*this);
      return 1;
    }
  return 0;
}

/*! @brief Initialize the renumerotation array for periodicity
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::init_renum_perio()
{
  const int_t nb_s = sommets_.dimension(0);
  IntTab_t renum(nb_s);
  for (int_t i = 0; i < nb_s; i++)
    renum[i] = i;
  set_renum_som_perio(renum);
}


/*! @brief Build the MEDCoupling mesh corresponding to the TRUST mesh.
 */
template<typename _SZ_>
void Domaine_32_64<_SZ_>::build_mc_mesh(bool virt) const
{
#ifdef MEDCOUPLING_
  Cerr << "Domaine: Creating a MEDCouplingUMesh object for the domain '" << le_nom() << "'" << finl;

  using MEDCoupling::DataArrayInt;
  using MEDCoupling::DataArrayDouble;

  // Initialize mesh
  Nom type_ele = elem_->que_suis_je();
  int mesh_dim;
  INTERP_KERNEL::NormalizedCellType cell_type = type_geo_trio_to_type_medcoupling(type_ele, mesh_dim);
  MCAuto<MEDCouplingUMesh>& mc_mesh = virt ? mc_mesh_virt_ : mc_mesh_;
  mc_mesh = MEDCouplingUMesh::New(nom_.getChar(), mesh_dim);

  //
  // Nodes
  //
  int_t nnodes = sommets_.dimension(0);
  MCAuto<DataArrayDouble> coord(DataArrayDouble::New());
  if (nnodes==0)
    coord->alloc(0, Objet_U::dimension);
  else
    // Avoid deep copy of vertices:
    coord->useArray(sommets_.addr(), false, MEDCoupling::DeallocType::CPP_DEALLOC, nnodes, Objet_U::dimension);

  coord->setInfoOnComponent(0, "x");
  coord->setInfoOnComponent(1, "y");
  if (Objet_U::dimension == 3) coord->setInfoOnComponent(2, "z");
  mc_mesh->setCoords(coord);

  //
  // Connectivity
  //
  int_t ncells = virt ? mes_elems_.dimension_tot(0) : mes_elems_.dimension(0);
  int nverts = (int)mes_elems_.dimension(1);

  // TRUST -> MED connectivity
  IntTab_t les_elems2(mes_elems_);
  conn_trust_to_med(les_elems2, type_ele, true);

  mc_mesh->allocateCells(ncells);
  if (cell_type == INTERP_KERNEL::NORM_POLYHED)
    {
      // Polyedron is special, see page 10:
      // http://trac.lecad.si/vaje/chrome/site/doc8.3.0/extra/Normalisation_pour_le_couplage_de_codes.pdf
      const Polyedre_32_64<_SZ_>& poly = ref_cast(Polyedre_32_64<_SZ_>, elem_.valeur());
      ArrOfInt_t nodes_glob;
      poly.remplir_Nodes_glob(nodes_glob, les_elems2);
      const ArrOfInt_t& facesIndex = poly.getFacesIndex();
      const ArrOfInt_t& polyhedronIndex = poly.getPolyhedronIndex();
      assert(ncells <= polyhedronIndex.size_array() - 1);

      for (int_t i = 0; i < ncells; i++)
        {
          int size = 0;
          for (int_t face = polyhedronIndex[i]; face < polyhedronIndex[i + 1]; face++)
            size += (int)(facesIndex[face + 1] - facesIndex[face] + 1);
          size--; // No -1 at the end of the cell
          ArrOfTID cell_def(size);  // ArrOfTID whatever the template parameter, since TID == mcIdType.
          size = 0;
          for (int_t face = polyhedronIndex[i]; face < polyhedronIndex[i + 1]; face++)
            {
              for (int_t node = facesIndex[face]; node < facesIndex[face + 1]; node++)
                cell_def[size++] = nodes_glob[node];
              if (size < cell_def.size_array())
                // Add -1 to mark the end of a face:
                cell_def[size++] = -1;
            }
          mc_mesh->insertNextCell(cell_type, cell_def.size_array(), cell_def.addr());
        }
    }
  else
    {
      // Other cells:
      if (std::is_same<_SZ_, trustIdType>::value) // 64b version of the Domaine, or TRUST compiled in 32b
        {
          // We can directly point into les_elems2, types are compatible
          for (int_t i = 0; i < ncells; i++)
            {
              int nvertices = nverts;
              // Polygons don't have a constant number of vertices - need to discard -1 values:
              for (int j = nverts-1; j >= 0 && les_elems2(i, j) < 0; j--) nvertices--;
              // Brutal pointer cast below, just so that the compiler does not complain when instanciating for _SZ_ = int:
              mc_mesh->insertNextCell(cell_type, nvertices, (trustIdType *)(les_elems2.addr() + i * nverts));
            }
        }
      else
        {
          // Need to upcast from int to mcIdType:
          for (int_t i = 0; i < ncells; i++)
            {
              ArrOfTID cell_def(nverts);
              int j = 0;
              for (; j<nverts && les_elems2(i, j) >= 0; j++)
                cell_def[j] = (trustIdType)les_elems2(i, j);
              mc_mesh->insertNextCell(cell_type, j, cell_def.addr());  // j is the final numb of vertices
            }
        }
    }
  *(virt ? &mc_mesh_virt_ready_ : &mc_mesh_ready_) = true;

#endif // MEDCOUPLING_
}


template<typename _SZ_>
void Domaine_32_64<_SZ_>::prepare_rmp_with(const Domaine_32_64& other_domain, bool virt) const
{
#ifdef MEDCOUPLING_
  using namespace MEDCoupling;

  // Retrieve mesh upfront to possibly build them if they were not already:
  get_mc_mesh();
  const MEDCouplingUMesh* oth_msh = other_domain.get_mc_mesh(virt);

  Cerr << "Building remapper between " << le_nom() << " (" << (int)mc_mesh_->getSpaceDimension() << "D) mesh with " << (int)mc_mesh_->getNumberOfCells()
       << " cells and " << other_domain.le_nom() << " (" << (int)oth_msh->getSpaceDimension() << "D) mesh with "
       << (int)oth_msh->getNumberOfCells() << " cells" << finl;
  rmps[&other_domain].prepare(oth_msh, mc_mesh_, "P0P0");
  Cerr << "remapper prepared with " << rmps.at(&other_domain).getNumberOfColsOfMatrix() << " columns in matrix, with max value = " << rmps.at(&other_domain).getMaxValueInCrudeMatrix() << finl;
#else
  Process::exit("Domaine_32_64<_SZ_>::prepare_rmp_with should not be called since it requires a TRUST version compiled with MEDCoupling !");
#endif
}

template <typename _SIZE_>
void Domaine_32_64<_SIZE_>::prepare_dec_with(const Domaine_32_64& other_domain, MEDCouplingFieldDouble *dist, MEDCouplingFieldDouble *loc) const
{
#if defined(MEDCOUPLING_) && defined(MPI_)
  using namespace MEDCoupling;

  Perf_counters::time_point t0 = statistics().start_clock();
  Cerr << "Building DEC of nature" << MEDCouplingNatureOfField::GetRepr(dist->getNature())
       << "from " << other_domain.le_nom() << " (" << Process::mp_sum(dist->getMesh()->getNumberOfCells())
       << " cells) to " << le_nom() << " (" << Process::mp_sum(loc->getMesh()->getNumberOfCells()) << " cells) : ";
  std::set<int> pcs;
  for (int i=0; i<Process::nproc(); i++) pcs.insert(i);
  /* a bit technical */
  decs.emplace(std::piecewise_construct,
               std::forward_as_tuple(&other_domain, dist->getNature()),
               std::forward_as_tuple(pcs, ref_cast(Comm_Group_MPI,PE_Groups::current_group()).get_trio_u_world()));
  OverlapDEC& dec = decs.at({ &other_domain, dist->getNature()});
  dec.setWorkSharingAlgo(Option_Interpolation::SHARING_ALGO);
  dec.attachSourceLocalField(dist);
  dec.attachTargetLocalField(loc);
  dec.synchronize();

  Cerr << statistics().compute_time(t0) << " s" << finl;
#else
  Process::exit("Domaine::prepare_dec_with() should not be called since it requires a TRUST version compiled with MEDCoupling and MPI!");
#endif
}

#ifdef MEDCOUPLING_

template <typename _SIZE_>
MEDCoupling::MEDCouplingRemapper* Domaine_32_64<_SIZE_>::get_remapper(const Domaine_32_64& other_domain, bool virt) const
{
  if (!rmps.count(&other_domain))
    prepare_rmp_with(other_domain, virt);
  return &rmps.at(&other_domain);
}

#ifdef MPI_
template <typename _SIZE_>
MEDCoupling::OverlapDEC* Domaine_32_64<_SIZE_>::get_dec(const Domaine_32_64& other_domain, MEDCouplingFieldDouble *dist, MEDCouplingFieldDouble *loc) const
{
  if (!decs.count({ &other_domain, dist->getNature() } ))
    prepare_dec_with(other_domain, dist, loc);
  return &decs.at({ &other_domain, dist->getNature() });
}
#endif

#endif


/*! @brief Fills the Domaine from a list of Domaine objects by aggregating them.
 *
 * See Mailler for example
 */
template <typename _SIZE_>
void Domaine_32_64<_SIZE_>::fill_from_list(std::list<Domaine_32_64*>& lst)
{
  Cerr << "Filling domain from list of domains in progress... " << finl;
  if (Process::is_parallel())
    Process::exit("Error in Domaine_32_64<_SIZE_>::fill_from_list() : compression prohibited in parallel mode");
  if (lst.size() == 0)
    Process::exit("Error in Domaine_32_64<_SIZE_>::fill_from_list() : compression prohibited in parallel mode");

  for(auto& elem: lst)
    elem->comprimer();

#ifndef NDEBUG
  Domaine_32_64& fst_dom = *lst.front();
  Nom typ_elem = fst_dom.type_elem()->que_suis_je();
#endif
  for(auto& it: lst)
    {
      Domaine_32_64& dom2 = *it;
      Cerr << "   Concatenating Domains "<< nom_ << " and " << dom2.nom_ << finl;;
      // Check single geometrical type:
      assert(typ_elem == dom2.type_elem()->que_suis_je());
      // Handle nodes:
      IntVect_t les_nums;
      // Copy sommets to this
      ajouter(dom2.sommets_, les_nums);  // les_nums: out parameter
      // Renumber current Domaine to prepare addition of elements
      dom2.renum(les_nums);
      // Merge elem info:
      merge_wo_vertices_with(dom2);
    }

  Cerr << "Filling from list - End!" << finl;
}

/*! @brief Renumbers the nodes and elements present in the common items of joints.
*
* Node number k becomes node number Les_Nums[k], and element number e becomes element number e+elem_offset.
*
* @param (IntVect& Les_Nums) the vector containing the new numbering: New_node_number_i = Les_Nums[Old_node_number_i]
*/
template <typename _SIZE_>
void Domaine_32_64<_SIZE_>::renum_joint_common_items(const IntVect_t& Les_Nums, const int_t elem_offset)
{
  for (int i_joint = 0; i_joint < nb_joints(); i_joint++)
    {
      ArrOfInt_t& sommets_communs = mes_faces_joint_[i_joint].set_joint_item(JOINT_ITEM::SOMMET).set_items_communs();
      for (int_t index = 0; index < sommets_communs.size_array(); index++)
        sommets_communs[index] = Les_Nums[sommets_communs[index]];

      ArrOfInt_t& elements_distants = mes_faces_joint_[i_joint].set_joint_item(JOINT_ITEM::ELEMENT).set_items_distants();
      elements_distants += elem_offset;
    }
}

/*! @brief Merges joints with the same name
 *
 */
template <typename _SIZE_>
int Domaine_32_64<_SIZE_>::comprimer_joints()
{
  auto& list = mes_faces_joint_.get_stl_list();
  for (auto it = list.begin(); it != list.end(); ++it)
    {
      Frontiere_t& front = *it;
      for (auto it2 = std::next(it); it2 != list.end();)
        {
          Frontiere_t& front2 = *it2;
          if (front.le_nom() == front2.le_nom())
            {
              front.add(front2);
              it2 = list.erase(it2);
            }
          else
            ++it2;
        }
    }
  return 1;
}


/////////////////////////////////////////////////
//// Methods only used in the 32 bits version
/////////////////////////////////////////////////

namespace  // Anonymous namespace - only 32 bits stuff here
{

/*! @brief This method performs a virtual space exchange of an edge array without going through the edge descriptor.
 *
 * The elem_aretes array and the virtual space exchange of elements are used.
 *
 */
void echanger_tableau_aretes(const IntTab& elem_aretes, int nb_aretes_reelles, ArrOfInt& tab_aretes)
{
  const int moi = Process::me();

  const int nb_elem = elem_aretes.dimension(0);
  const int nb_elem_tot = elem_aretes.dimension_tot(0);
  const int nb_aretes_elem = elem_aretes.dimension(1);
  int i;

  // **********************
  // I) Exchange to update common items
  //  Slightly complex algorithm to update common items: for each real edge,
  //  the value of tab_aretes must equal the initial value of tab_arete given by
  //  the processor with the smallest rank among those sharing the edge (i.e.
  //  processors that have an element adjacent to this edge).

  // Array to identify the owner processor of a real edge
  ArrOfInt pe_arete(nb_aretes_reelles);
  pe_arete = moi;
  // Array giving, for each element, the owner processor
  IntVect pe_elem(nb_elem_tot);
  pe_elem = moi; // initialized with "me"
  {
    pe_elem.set_md_vector(elem_aretes.get_md_vector());
    pe_elem.echange_espace_virtuel();
    // Store in pe_arete the number of the smallest-rank owner processor among
    // the processors owning elements adjacent to this edge.
    // No need to iterate over real elements, we would find pe_elem[i]==moi...
    // If the edge is on a processor with a lower rank, assign it that rank
    for (i = nb_elem; i < nb_elem_tot; i++)
      for (int pe = pe_elem[i], j = 0, a; j < nb_aretes_elem && (a = elem_aretes(i, j)) >= 0; j++)
        if (a < nb_aretes_reelles && pe_arete[a] > pe)
          pe_arete[a] = pe;
  }
  // Assuming the virtual element space contains at least one layer of virtual elements
  //   (all neighbors of real elements through vertices), the real edges are exchanged
  //   (virtual edges not yet).
  // In this case, pe_arete is now correctly filled for real edges.

  IntTab tmp;
  tmp.copy(elem_aretes, RESIZE_OPTIONS::NOCOPY_NOINIT); // copy structure only

  // Copy tab_aretes into the tmp structure (we can exchange tmp, not tab_aretes)
  for (i = 0; i < nb_elem; i++)
    for (int j = 0, a; j < nb_aretes_elem && (a = elem_aretes(i, j)) >= 0; j++)
      tmp(i, j) = tab_aretes[a];

  // 2) Exchange the array
  tmp.echange_espace_virtuel();

  // 3) Copy back into the real part of tab_aretes the values taken from tmp:
  //    for an edge shared by several procs, the proc with the smallest rank
  //    provides the value.
  // No need to iterate over real elements, the value would not change
  for (i = nb_elem; i < nb_elem_tot; i++)
    for (int pe = pe_elem[i], j = 0, a; j < nb_aretes_elem && (a = elem_aretes(i, j)) >= 0; j++)
      if (a < nb_aretes_reelles && pe_arete[a] == pe)
        tab_aretes[a] = tmp(i, j);

  // tab_aretes now contains correct values for all real edges
  //  (common items are up to date). We do one more exchange via tmp to
  //  update the virtual items:

  // ******************
  // II) Exchange to update the virtual space of edges

  // Copy tab_aretes into the tmp structure again
  for (i = 0; i < nb_elem; i++)
    for (int j = 0, a; j < nb_aretes_elem && (a = elem_aretes(i, j)) >= 0; j++)
      tmp(i, j) = tab_aretes[a];

  // Exchange the array
  tmp.echange_espace_virtuel();
  // Copy tmp back into tab_aretes
  for (i = nb_elem; i < nb_elem_tot; i++)
    for (int j = 0, a; j < nb_aretes_elem && (a = elem_aretes(i, j)) >= 0; j++)
      tab_aretes[a] = tmp(i, j);
}

} // end anonymous namespace

/*! Selects a unique item (vertex, face ...) from a list (item_possible)
* in order to ensure the parallelism of certain algorithms.
* The selection is made by testing the distance between the coordinates (coord_possible)
* locating these items with respect to the coordinates (coord_ref) of a reference point.
* The retained item is the one with the minimum distance to the reference point.
* If several items remain at the same distance from the reference point,
* the test is repeated by translating the reference point.
*/
template <>
int Domaine_32_64<int>::identifie_item_unique(IntList& item_possible, DoubleTab& coord_possible, const DoubleVect& coord_ref)
{
  int it_selection = -1;
  DoubleTab decentre_face(4, Objet_U::dimension);
  decentre_face = 0.;
  for (int t = 1; t < 4; t++)
    for (int dir = 0; dir < Objet_U::dimension; dir++)
      if (dir == (t - 1))
        decentre_face(t, dir) = 1.;
  // decentre_face(0,0:dim)={0,0,0}
  // decentre_face(1,0:dim)={1,0,0}
  // decentre_face(2,0:dim)={0,1,0}
  // decentre_face(3,0:dim)={0,0,1}

  //At the first pass (t=0) no translation is performed
  DoubleVect dist;
  assert(item_possible.size() != 0);
  int t = 0;
  while ((item_possible.size() != 1) && (t < 4))
    {
      double distmin = DMAXFLOAT;
      int size_initiale = item_possible.size();
      dist.resize(size_initiale);
      dist = 0.;

      for (int ind_it = 0; ind_it < size_initiale; ind_it++)
        {
          for (int dir = 0; dir < Objet_U::dimension; dir++)
            dist[ind_it] += (coord_possible(ind_it, dir) - (coord_ref(dir) + decentre_face(t, dir))) * (coord_possible(ind_it, dir) - (coord_ref(dir) + decentre_face(t, dir)));
          if (dist[ind_it] <= distmin)
            distmin = dist[ind_it];
        }

      int ind_it = 0;
      int nb_it_suppr = 0;
      while (ind_it < size_initiale)
        {
          if (!est_egal(dist[ind_it], distmin))
            {
              int ind_it_suppr = ind_it - nb_it_suppr;
              int it_suppr = item_possible[ind_it_suppr];
              item_possible.suppr(it_suppr);

              int size_actuelle = item_possible.size();
              for (int ind = ind_it_suppr; ind < size_actuelle; ind++)
                for (int dir = 0; dir < dimension; dir++)
                  coord_possible(ind, dir) = coord_possible(ind + 1, dir);
              coord_possible.resize(size_actuelle, dimension, RESIZE_OPTIONS::COPY_NOINIT);
              nb_it_suppr++;
            }
          ind_it++;
        }
      t++;
    }
  if (item_possible.size() == 1)
    it_selection = item_possible[0];
  else
    {
      Cerr << "Domaine::identifie_item_unique()" << finl;
      Cerr << "An item has not been found among the list." << finl;
      Cerr << "Please contact TRUST support." << finl;
      Process::exit();
    }
  return it_selection;
}

template <typename _SIZE_>
int Domaine_32_64<_SIZE_>::identifie_item_unique(IntList& item_possible, DoubleTab& coord_possible, const DoubleVect& coord_ref)
{
  assert(false);
  throw;
}

/*! @brief Method called by Domaine_VF::discretiser().
 *
 * Builds the descriptor for boundary faces.
 *   Fills ind_faces_virt_bord and the get_faces_virt() arrays of the boundaries
 *   from the parallel descriptor of faces.
 *   Note B.M.: having placed faces in Domaine_VF, edges in Domaine,
 *    some face boundary properties in Domaine_VF and others in Domaine
 *    makes the initialization follow somewhat convoluted paths... this should be cleaned up.
 *
 */
template <>
void Domaine_32_64<int>::init_faces_virt_bord(const MD_Vector& md_vect_faces, MD_Vector& md_vect_faces_front)
{
  if (Process::is_sequential()) // Much simpler in this case:
    {
      ind_faces_virt_bord_.resize_array(0);
      MD_Vector_seq mdseq(nb_faces_frontiere());
      md_vect_faces_front.copy(mdseq);

      // Build the MD_Vector_seq of each boundary:
      const int nb_frontieres = nb_front_Cl() + nb_groupes_faces();
      for (int i_frontiere = 0; i_frontiere < nb_frontieres; i_frontiere++)
        {
          Frontiere& front = frontiere(i_frontiere);
          IntTab& faces_sommets_frontiere = front.les_sommets_des_faces();
          // Some problems have multiple Domaine_VF objects attached to the same Domaine (radiation)
          // If we already went through here, don't redo the work:
          if (faces_sommets_frontiere.get_md_vector())
            continue;
          const int nb_faces_front = front.nb_faces();
          // Build a descriptor containing the subset of faces of this boundary
          MD_Vector md_frontiere;
          MD_Vector_seq mdseq_front(nb_faces_front);
          md_frontiere.copy(mdseq_front);
          faces_sommets_frontiere.set_md_vector(md_frontiere);
        }

      return;
    }

  // ***************************************
  // 1) Build array structures for all boundary faces
  //   (faces from 0 to nb_faces_frontiere())
  const int nb_faces_fr = nb_faces_frontiere();
  //  Mark boundary faces (-1=>not a boundary face, 0=>boundary face)
  IntVect vect_renum;
  MD_Vector_tools::creer_tableau_distribue(md_vect_faces, vect_renum, RESIZE_OPTIONS::NOCOPY_NOINIT);
  vect_renum = -1;
  for (int i = 0; i < nb_faces_fr; i++)
    vect_renum[i] = 0;
  vect_renum.echange_espace_virtuel();

  // Create the descriptor for boundary faces (by extracting a subset of the face descriptor).
  // The default numbering in ascending order is used:
  MD_Vector_tools::creer_md_vect_renum_auto(vect_renum, md_vect_faces_front);

  //  Fill the ind_faces_virt_bord array. It is just the virtual part of the renum array.
  //  (the real part is trivial: it is a contiguous numbering from 0 to nb_faces_frontiere())
  const int nb_faces = vect_renum.size();
  const int nb_faces_tot = vect_renum.size_totale();
  const int nb_faces_virt = nb_faces_tot - nb_faces;
  ind_faces_virt_bord_.resize_array(nb_faces_virt, RESIZE_OPTIONS::NOCOPY_NOINIT);
  for (int i = 0; i < nb_faces_virt; i++)
    ind_faces_virt_bord_[i] = vect_renum[nb_faces + i];

  // **************************************
  // 2) Build array structures for each boundary

  // Fill the arrays
  //   frontiere(i).get_faces_virt() for 0 <= i < nb_front_Cl()
  // This array contains the indices in Domaine_VF of the virtual faces
  // that are on boundary i.
  // Compute the virtual space of faces for each boundary

  // Number of boundaries:
  const int nb_frontieres = nb_front_Cl();
  int i_frontiere;
  // Fill the get_faces_virt() arrays:
  // and build the MD_Vector of each boundary (associated with the face array)
  for (i_frontiere = 0; i_frontiere < nb_frontieres; i_frontiere++)
    {
      Frontiere& front = frontiere(i_frontiere);
      IntTab& faces_sommets_frontiere = front.les_sommets_des_faces();
      // Some problems have multiple Domaine_VF objects attached to the same Domaine (radiation)
      // If we already went through here, don't redo the work:
      if (faces_sommets_frontiere.get_md_vector())
        continue;
      //the faces_sommets_frontiere arrays must have the same width on all procs before exchange
      int nb_som_faces = Process::mp_max(faces_sommets_frontiere.dimension(1));
      if (faces_sommets_frontiere.dimension(1) < nb_som_faces)
        {
          IntTab fsf_old;
          fsf_old = faces_sommets_frontiere;
          faces_sommets_frontiere.resize(fsf_old.dimension_tot(0), nb_som_faces);
          faces_sommets_frontiere = -1;
          for (int i = 0, j; i < fsf_old.dimension_tot(0); i++)
            for (j = 0; j < fsf_old.dimension(1); j++)
              faces_sommets_frontiere(i, j) = fsf_old(i, j);
        }

      vect_renum = -1;
      const int i_premiere_face = front.num_premiere_face();
      const int nb_faces_front = front.nb_faces();
      // Mark the faces of this boundary
      for (int i = i_premiere_face; i < i_premiere_face + nb_faces_front; i++)
        vect_renum[i] = 0;
      vect_renum.echange_espace_virtuel();
      // Build a descriptor containing the subset of faces of this boundary
      MD_Vector md_frontiere;
      MD_Vector_tools::creer_md_vect_renum_auto(vect_renum, md_frontiere);

      // Create the virtual space of boundary faces
      // (this is where the md_frontiere descriptor is associated with the face array)
      const MD_Vector& md_sommets = les_sommets().get_md_vector();
      Scatter::construire_espace_virtuel_traduction(md_frontiere, /* array indexed by boundary face numbers */
                                                    md_sommets, /* containing vertex indices of the domain */
                                                    faces_sommets_frontiere, /* array to process */
                                                    1 /* fatal error: if a vertex is missing, it is an error */);

      // Retrieve from renum the renumbered index of each face:
      //  extract the indices of the virtual faces of this boundary
      ArrOfInt& tab = front.get_faces_virt();
      assert(faces_sommets_frontiere.dimension(0) == nb_faces_front);
      const int nb_faces_tot_frontiere = faces_sommets_frontiere.dimension_tot(0);
      const int nb_faces_virt_frontiere = nb_faces_tot_frontiere - nb_faces_front;
      tab.resize_array(nb_faces_virt_frontiere);
      const int ndebut = nb_faces; // number of faces in the Domain!
      const int nfin = nb_faces_tot; // idem!
      for (int i = ndebut; i < nfin; i++)
        {
          const int j = vect_renum[i];
          if (j >= 0)
            {
              assert(j >= nb_faces_front && j < nb_faces_tot_frontiere);
              // Face i is virtual and on this boundary
              tab[j - nb_faces_front] = i;
            }
        }
    }
}

template <typename _SIZE_>
void Domaine_32_64<_SIZE_>::init_faces_virt_bord(const MD_Vector& md_vect_faces, MD_Vector& md_vect_faces_front)
{
  assert(false);
  throw;
}

/*! Version of creer_aretes compatible with polyhedra
  */
template <>
void Domaine_32_64<int>::creer_aretes()
{
  const IntTab& elem_som = les_elems();
  // Number of real elements:
  const int nbelem = elem_som.dimension(0);
  // Virtual elements are already built:
  const int nbelem_tot = elem_som.dimension_tot(0);

  aretes_som_.resize(0, 2);
  bool is_poly = sub_type(Poly_geom_base, type_elem().valeur());

  std::vector<std::vector<int> > v_e_a(nbelem_tot);  //list of edges for each element
  int nb_aretes_reelles = 0, i;
  int j;
  {
    // A linked list to retrieve, for each vertex, the list of edges
    // attached to that vertex. The array has the same size as Aretes_som.dimension(0).
    // chaine_aretes_sommets[i] contains the index of the next edge attached to
    // the same vertex, or -1 if it is the last one.
    ArrOfInt chaine_aretes_sommets;
    // Index of the first edge attached to each vertex in chaine_aretes_sommets
    ArrOfInt premiere_arete_som(nb_som_tot());
    premiere_arete_som = -1;

    std::map<std::array<double, 3>, std::array<int, 2> > aretes_loc; //edges of the current element: aretes_loc[{xa, ya, za}] = { s1, s2}
    //using a map ensures that edges are in the same order on all procs!
    for (int i_elem = 0; i_elem < nbelem_tot; aretes_loc.clear(), i_elem++)
      {
        /* 1. retrieve the edges of the element by iterating over its faces */
        const Elem_geom_base& elem_g = ref_cast(Elem_geom_base, type_elem().valeur());
        IntTab f_e_r;
        if (is_poly)
          {
            const Poly_geom_base& poly_g = ref_cast(Poly_geom_base, type_elem().valeur());
            poly_g.get_tab_faces_sommets_locaux(f_e_r, i_elem);
          }
        else
          elem_g.get_tab_faces_sommets_locaux(f_e_r);

        for (i = 0; i < f_e_r.dimension(0) && f_e_r(i, 0) >= 0; i++)
          for (j = 0; j < f_e_r.dimension(1) && f_e_r(i, j) >= 0; j++)
            {
              int s1 = elem_som(i_elem, f_e_r(i, j)), s2 = elem_som(i_elem, f_e_r(i, j + 1 < f_e_r.dimension(1) && f_e_r(i, j + 1) >= 0 ? j + 1 : 0));
              std::array<double, 3> key;
              for (int l = 0; l < 3; l++)
                key[l] = (sommets_(s1, l) + sommets_(s2, l)) / 2;
              aretes_loc[key] = {{ std::min(s1, s2), std::max(s1, s2) }};
            }

        for (auto &&kv : aretes_loc)
          {
            //have we already seen this edge ?
            int k = premiere_arete_som[kv.second[0]];
            while (k >= 0 && (aretes_som_(k, 0) != kv.second[0] || aretes_som_(k, 1) != kv.second[1]))
              k = chaine_aretes_sommets[k];
            if (k < 0) //edge not yet found -> update premiere_arete_som and chaine_arete_sommets
              {
                // The edge does not exist yet
                k = chaine_aretes_sommets.size_array();
                assert(k == aretes_som_.dimension(0));
                aretes_som_.append_line(kv.second[0], kv.second[1]);
                // Insert the edge at the head of the linked list
                int old_head = premiere_arete_som[kv.second[0]];
                // Index of the new edge
                int new_head = chaine_aretes_sommets.size_array();
                chaine_aretes_sommets.append_array(old_head);
                premiere_arete_som[kv.second[0]] = new_head;
              }
            v_e_a[i_elem].push_back(k); //add the edge to the element's edge list
          }
        if (i_elem == nbelem - 1)
          {
            // We have just finished the real edges
            nb_aretes_reelles = aretes_som_.dimension(0);
          }
      }
  }
  /* fill the elem_aretes array using v_e_a */
  int nb_aretes_elem = 0;
  for (i = 0; i < nbelem_tot; i++)
    nb_aretes_elem = std::max(nb_aretes_elem, (int) v_e_a[i].size());
  nb_aretes_elem = mp_max(nb_aretes_elem);
  elem_aretes_.resize(0, nb_aretes_elem);
  creer_tableau_elements(elem_aretes_, RESIZE_OPTIONS::NOCOPY_NOINIT);
  for (i = 0, elem_aretes_ = -1; i < nbelem_tot; i++)
    for (j = 0; j < (int) v_e_a[i].size(); j++)
      elem_aretes_(i, j) = v_e_a[i][j];

  // Adjust the size of the Aretes_som array
  const int n_aretes_tot = aretes_som_.dimension(0); // note: nb_aretes_tot is a method!
  aretes_som_.append_line(-1, -1); // because the following resize only does something if the size changes
  aretes_som_.resize(n_aretes_tot, 2);

  Journal() << "Domaine " << le_nom() << " nb_aretes=" << nb_aretes_reelles << " nb_aretes_tot=" << n_aretes_tot << finl;

  // Build the parallel descriptor
  {
    // For each edge, index of the processor owning the edge
    const int moi = Process::me();
    ArrOfInt pe_aretes(n_aretes_tot);
    pe_aretes = moi;
    echanger_tableau_aretes(elem_aretes_, nb_aretes_reelles, pe_aretes);

    // For each edge, index of the edge on the owning processor
    ArrOfInt indice_aretes_owner;
    indice_aretes_owner.resize_array(n_aretes_tot, RESIZE_OPTIONS::NOCOPY_NOINIT);
    for (i = 0; i < nb_aretes_reelles; i++)
      indice_aretes_owner[i] = i;
    echanger_tableau_aretes(elem_aretes_, nb_aretes_reelles, indice_aretes_owner);

    // Build pe_voisins
    ArrOfInt pe_voisins;
    for (i = 0; i < n_aretes_tot; i++)
      if (pe_aretes[i] != moi)
        pe_voisins.append_array(pe_aretes[i]);

    ArrOfInt liste_pe;
    reverse_send_recv_pe_list(pe_voisins, liste_pe);

    // Concatenate the two lists.
    for (i = 0; i < liste_pe.size_array(); i++)
      pe_voisins.append_array(liste_pe[i]);
    array_trier_retirer_doublons(pe_voisins);

    int nb_voisins = pe_voisins.size_array();
    ArrOfInt indices_pe(nproc());
    indices_pe = -1;
    for (i = 0; i < nb_voisins; i++)
      indices_pe[pe_voisins[i]] = i;

    ArrsOfInt aretes_communes_to_recv(nb_voisins);
    ArrsOfInt blocs_aretes_virt(nb_voisins);
    ArrsOfInt aretes_to_send(nb_voisins);
    // Iterate over edges: look for edges to receive from another processor.
    // Real edges (common items)
    for (i = 0; i < nb_aretes_reelles; i++)
      {
        const int pe = pe_aretes[i];
        if (pe != moi)
          {
            const int indice_pe = indices_pe[pe];
            if (indice_pe < 0)
              {
                Cerr << "Error: indice_pe=" << indice_pe << " shouldn't be negative in Domaine_32_64<_SZ_>::creer_aretes." << finl;
                Cerr << "It is a TRUST bug on this mesh with the Pa discretization, contact support." << finl;
                Cerr << "You could also try another partitioned mesh to get around this issue." << finl;
                Process::exit();
              }
            // I receive this edge from another proc
            const int indice_distant = indice_aretes_owner[i];
            aretes_to_send[indice_pe].append_array(indice_distant); // index on the neighboring pe
            aretes_communes_to_recv[indice_pe].append_array(i); // local index of the edge
          }
      }
// Virtual edges
    for (i = nb_aretes_reelles; i < n_aretes_tot; i++)
      {
        const int pe = pe_aretes[i];
        assert(pe < nproc() && pe != moi);
        const int indice_pe = indices_pe[pe];
        if (indice_pe < 0)
          {
            Cerr << "Error: indice_pe=" << indice_pe << " shouldn't be negative in Domaine_32_64<_SZ_>::creer_aretes." << finl;
            Cerr << "It is a TRUST bug on this mesh with the Pa discretization, contact support." << finl;
            Cerr << "You could also try another partitioned mesh to get around this issue." << finl;
            Process::exit();
          }
        const int indice_distant = indice_aretes_owner[i];
        aretes_to_send[indice_pe].append_array(indice_distant); // index on the neighboring pe
        MD_Vector_base::append_item_to_blocs(blocs_aretes_virt[indice_pe], i);
      }
    {
      Schema_Comm schema;
      schema.set_send_recv_pe_list(pe_voisins, pe_voisins);
      schema.begin_comm();
      // Push the aretes_to_send array and the number of edges shared with this pe:
      for (i = 0; i < nb_voisins; i++)
        schema.send_buffer(pe_voisins[i]) << aretes_to_send[i];
      schema.echange_taille_et_messages();
      // Receive
      for (i = 0; i < nb_voisins; i++)
        schema.recv_buffer(pe_voisins[i]) >> aretes_to_send[i];
      schema.end_comm();
    }

    MD_Vector md;
    // Build the descriptor object
    if (Process::is_parallel())
      {
        MD_Vector_std md_aretes(n_aretes_tot, nb_aretes_reelles, pe_voisins, aretes_to_send, aretes_communes_to_recv, blocs_aretes_virt);
        md.copy(md_aretes);
      }
    else
      {
        MD_Vector_seq md_aretes(n_aretes_tot);
        md.copy(md_aretes);
      }
    Cerr << "Total number of edges = " << md->get_nb_items_tot() << finl;

    // Attach the descriptor to the array
    aretes_som_.set_md_vector(md);
  }
}

template <typename _SIZE_>
void Domaine_32_64<_SIZE_>::creer_aretes()
{
  assert(false);
  throw;
}

/*! Creation of boundary domains (called during discretisation).
 * Currently a static list of Domains where we need to know
 * the first element for each domain.
 */
template <>
void Domaine_32_64<int>::creer_mes_domaines_frontieres(const Domaine_VF& domaine_vf)
{
  const Nom expr_elements("1");
  const Nom expr_faces("1");
  int nb_frontieres = nb_front_Cl();
  domaines_frontieres_.vide();

  for (int i=0; i<nb_frontieres; i++)
    {
      // Name of the boundary
      Noms nom_frontiere(1);
      nom_frontiere[0]=frontiere(i).le_nom();
      // Name of the surface domain to be built
      Nom nom_domaine_surfacique=le_nom();
      nom_domaine_surfacique+="_boundaries_";
      nom_domaine_surfacique+=frontiere(i).le_nom();
      // Creation
      Cerr << "Creating a surface domain named " << nom_domaine_surfacique << " for the boundary " << nom_frontiere[0] << " of the domain " << le_nom() << finl;

      Interprete_bloc& interp = Interprete_bloc::interprete_courant();
      if (interp.objet_global_existant(nom_domaine_surfacique))
        {
          Cerr << "Domain " << nom_domaine_surfacique
               << " already exists, writing to this object." << finl;

          Domaine& dom_new = ref_cast(Domaine, interprete().objet(nom_domaine_surfacique));
          Scatter::uninit_sequential_domain(dom_new);
        }
      else
        {
          DerObjU ob;
          ob.typer("Domaine");
          interp.ajouter(nom_domaine_surfacique, ob);
        }
      Domaine& dom_new = ref_cast(Domaine, interprete().objet(nom_domaine_surfacique));

      Extraire_surface::extraire_surface(dom_new,*this,nom_domaine_surfacique,domaine_vf,expr_elements,expr_faces,0,nom_frontiere);
      OBS_PTR(Domaine)& ref_dom_new=domaines_frontieres_.add(OBS_PTR(Domaine)());
      ref_dom_new=dom_new;
    }
}

template <typename _SIZE_>
void Domaine_32_64<_SIZE_>::creer_mes_domaines_frontieres(const Domaine_VF& domaine_vf)
{
  assert(false);
  throw;
}


/*! @brief Renumbering of nodes: node number k becomes node number Les_Nums[k]
 *
 * @param (IntVect& Les_Nums) vector containing the new numbering: New_node_number_i = Les_Nums[Old_node_number_i]
 */
template <typename _SIZE_>
void Domaine_32_64<_SIZE_>::renum(const IntVect_t& Les_Nums)
{
  int_t dim0 = mes_elems_.dimension(0);
  int dim1 = mes_elems_.dimension_int(1);

  for (int_t i = 0; i < dim0; i++)
    for (int j = 0; j < dim1; j++)
      mes_elems_(i, j) = Les_Nums[mes_elems_(i, j)];

  for (int i = 0; i < nb_bords(); i++)
    mes_faces_bord_(i).renum(Les_Nums);
  for (int i = 0; i < nb_joints(); i++)
    mes_faces_joint_(i).renum(Les_Nums);
  for (int i = 0; i < nb_raccords(); i++)
    mes_faces_raccord_(i)->renum(Les_Nums);
  for (int i = 0; i < nb_frontieres_internes(); i++)
    mes_bords_int_(i).renum(Les_Nums);
  for (int i = 0; i < nb_groupes_faces(); i++)
    mes_groupes_faces_(i).renum(Les_Nums);
}

template<>
void Domaine_32_64<int>::construire_renum_som_perio(const Conds_lim& les_cl, const Domaine_dis_base& domaine_dis)
{
  // Sanity check - make sure that the periodic BC are put on a periodic boundary
  // (the opposite is allowed, eventhough it is probably stupid: a non periodic BC on a periodic boundary)
  const int nb_bords = les_cl.size();
  const Noms& bords_per = this->bords_perio();
  for (int n_bord = 0; n_bord < nb_bords; n_bord++)
    {
      if (sub_type(Periodique, les_cl[n_bord].valeur()))
        {
          const Nom& nom_b =les_cl[n_bord]->frontiere_dis().frontiere().le_nom();
          if(bords_per.rang(nom_b) < 0)
            {
              Cerr << "ERROR: you have put a periodic boundary condition on a boundary ('" << nom_b << "') which is not periodic." << finl;
              Cerr << "Use the keyword: declarer_bord_perio { domaine " << le_nom() << " bord " << nom_b << " }" << finl;
              Cerr << "after the loading of the domain to declare this boundary as being periodic." << finl;
              Process::exit();
            }
        }
    }

  Reordonner_faces_periodiques::renum_som_perio(*this, renum_som_perio_,
                                                1 /* Compute values for virtual vertices */);
}


template<typename _SZ_>
void Domaine_32_64<_SZ_>::construire_renum_som_perio(const Conds_lim& les_cl, const Domaine_dis_base& domaine_dis)
{
  assert(false);
  throw;
}



/////////////////////////////////////////////////
//// Template instanciations
/////////////////////////////////////////////////

template class Domaine_32_64<int>;
#if INT_is_64_ == 2
template class Domaine_32_64<trustIdType>;
#endif
