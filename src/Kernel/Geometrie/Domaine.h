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

#ifndef Domaine_32_64_included
#define Domaine_32_64_included

#include <Domaine_base.h>

#include <TRUST_Deriv.h>
#include <Octree.h>
#include <Elem_geom.h>

#include <medcoupling++.h>

#ifdef MEDCOUPLING_
#include <MEDCouplingFieldTemplate.hxx>
#include <MEDCouplingUMesh.hxx>
#include <MEDCouplingRemapper.hxx>
using MEDCoupling::MEDCouplingRemapper;
using MEDCoupling::MEDCouplingUMesh;
using MEDCoupling::MCAuto;
using MEDCoupling::MEDCouplingFieldDouble;
#ifdef MPI_
#include <OverlapDEC.hxx>
using MEDCoupling::OverlapDEC;
#endif
#endif

// Forward decl:
class Domaine_dis_base;
class Conds_lim;
class Reorder_Mesh;

template <typename _SIZE_> class OctreeRoot_32_64;
template <typename _SIZE_> class Sous_Domaine_32_64;

/*! @brief class Domaine_32_64 A Domain is a mesh composed of a set of geometric elements of the same type.
 *
 * The different types of elements are objects of classes derived from Elem_geom_base.
 * A domain is made up of nodes, elements, boundaries, periodic boundaries,
 * joints, connections and internal boundaries.
 *
 * This class is templatized on the 32/64 bit configuration.
 * All the methods/members not sensitive to this are in Domaine_base.
 *
 * @sa Domaine, Domaine_64, Sous_Domaine, Frontiere, Elem_geom, Elem_geom_base, Bord, Bord_perio, Joint, Raccord, Bords_Interne
 */
template<typename _SIZE_>
class Domaine_32_64 : public Domaine_base
{

  Declare_instanciable_sans_constructeur_32_64( Domaine_32_64 );

public:
  //
  // Types
  //
  using int_t = _SIZE_;
  using ArrOfInt_t = ArrOfInt_T<_SIZE_>;
  using IntVect_t = IntVect_T<_SIZE_>;
  using IntTab_t = IntTab_T<_SIZE_>;
  using SmallArrOfTID_t = SmallArrOfTID_T<_SIZE_>;
  using ArrOfDouble_t= ArrOfDouble_T<_SIZE_>;
  using DoubleVect_t = DoubleVect_T<_SIZE_>;
  using DoubleTab_t = DoubleTab_T<_SIZE_>;

  using ArrsOfInt_t = ArrsOfInt_T<_SIZE_>;
  using DoubleTabs_t = DoubleTabs_T<_SIZE_>;

  using OctreeRoot_t = OctreeRoot_32_64<_SIZE_>;
  using Sous_Domaine_t = Sous_Domaine_32_64<_SIZE_>;
  using Bord_t = Bord_32_64<_SIZE_>;
  using Bords_t = Bords_32_64<_SIZE_>;
  using Bord_Interne_t = Bord_Interne_32_64<_SIZE_>;
  using Bords_Internes_t = Bords_Internes_32_64<_SIZE_>;
  using Groupe_Faces_t = Groupe_Faces_32_64<_SIZE_>;
  using Groupes_Faces_t = Groupes_Faces_32_64<_SIZE_>;
  using Frontiere_t = Frontiere_32_64<_SIZE_>;
  using Raccord_t = OWN_PTR(Raccord_base_32_64<_SIZE_>);
  using Raccords_t = Raccords_32_64<_SIZE_>;
  using Joint_t = Joint_32_64<_SIZE_>;
  using Joints_t = Joints_32_64<_SIZE_>;

  //
  // General
  //
  Domaine_32_64() { clear(); }
  Entree& readOn_has_perio(Entree& s, bool& has_perio);
  inline void typer(const Nom&);
  inline const OWN_PTR(Elem_geom_base_32_64<_SIZE_>)& type_elem() const { return elem_; }
  inline OWN_PTR(Elem_geom_base_32_64<_SIZE_>)& type_elem() {  return elem_; }
  inline void reordonner() { elem_->reordonner(); }


  //
  // Nodes
  //
  inline double coord(int_t i, int j) const { return sommets_(i,j); }
  inline double& coord(int_t i, int j) { return sommets_(i,j); }
  inline const DoubleTab_t& coord_sommets() const { return sommets_; }
  inline DoubleTab_t& les_sommets()  { return sommets_; }
  inline const DoubleTab_t& les_sommets() const  { return sommets_; }
  inline void saveSommetsCoordinates() { sommets_n=sommets_; }
  inline void resetSommetsCoordinates() {sommets_=sommets_n; }
  DoubleTab getBoundingBox() const;
  void ajouter(const DoubleTab_t& soms);
  void ajouter(const DoubleTab_t& soms, IntVect_t& nums);
  /// Returns the number of vertices of the domain.
  int_t nb_som() const { return sommets_.dimension(0); }
  /// Returns the total number of vertices of the domain i.e. the number of real and virtual vertices on the current processor.
  int_t nb_som_tot() const { return sommets_.dimension_tot(0); }
  void read_vertices(Entree& s);

  //
  // Elements
  //
  inline IntTab_t& les_elems()   { return mes_elems_; }
  inline const IntTab_t& les_elems() const  {  return mes_elems_; }
  inline int_t nb_elem() const  { return mes_elems_.dimension(0); }
  inline int_t nb_elem_tot() const { return mes_elems_.dimension_tot(0); }
  inline int nb_som_elem() const;
  inline int nb_faces_elem(int=0) const;
  /// @brief Returns the (global) number of the j-th vertex of the i-th element
  inline int_t sommet_elem(int_t i, int j) const  {  return mes_elems_(i,j); }

  //
  // Aretes
  //

  /// Returns the number of real edges.
  inline int_t nb_aretes() const { return aretes_som_.dimension(0); }
  /// returns the total number of edges (real+virtual).
  inline int_t nb_aretes_tot() const { return aretes_som_.dimension_tot(0); }

  //
  // Correspondances
  //

  /// returns the number of the j-th vertex of the i-th edge.
  inline int_t arete_sommets(int_t i, int j) const { return aretes_som_(i, j); }
  /// returns the number of the j-th edge of the i-th element.
  inline int_t elem_aretes(int_t i, int j) const {  return elem_aretes_(i, j); }
  /// returns the connectivity array edges/vertices.
  inline const IntTab_t& aretes_som() const {  return aretes_som_; }
  /// returns the connectivity array elements/edges.
  inline const IntTab_t& elem_aretes() const {   return elem_aretes_; }
  inline IntTab_t& set_aretes_som()   {  return aretes_som_; }
  inline IntTab_t& set_elem_aretes()  {  return elem_aretes_; }

  //
  // Faces
  //
  inline int_t nb_faces_bord() const { return mes_faces_bord_.nb_faces(); }
  inline int_t nb_faces_bord(int num_bord) const  { return mes_faces_bord_(num_bord).nb_faces(); }
  inline int_t nb_faces_bord(Type_Face type) const { return mes_faces_bord_.nb_faces(type); }
  inline int_t nb_faces_joint() const { return mes_faces_joint_.nb_faces(); }
  inline int_t nb_faces_joint(int num_joint) const { return mes_faces_joint_(num_joint).nb_faces(); }
  inline int_t nb_faces_joint(Type_Face type) const { return mes_faces_joint_.nb_faces(type); }
  inline int_t nb_faces_raccord() const { return mes_faces_raccord_.nb_faces(); }
  inline int_t nb_faces_raccord(int num_rac) const { return mes_faces_raccord_(num_rac)->nb_faces(); }
  inline int_t nb_faces_raccord(Type_Face type) const { return mes_faces_raccord_.nb_faces(type); }
  inline int_t nb_faces_bords_int() const  { return mes_bords_int_.nb_faces(); }
  inline int_t nb_faces_bords_int(int num_bord) const { return mes_bords_int_(num_bord).nb_faces(); }
  inline int_t nb_faces_bords_int(Type_Face type) const { return mes_bords_int_.nb_faces(type); }
  inline int_t nb_faces_groupes_faces() const { return mes_groupes_faces_.nb_faces(); }
  inline int_t nb_faces_groupes_faces(int num_g) const { return mes_groupes_faces_(num_g).nb_faces(); }
  inline int_t nb_faces_groupes_faces(Type_Face type) const { return mes_groupes_faces_.nb_faces(type); }

  inline int_t nb_faces_frontiere() const;
  inline int_t nb_faces_frontiere(Type_Face type) const;

  inline int_t nb_faces_specifiques() const;
  inline int_t nb_faces_specifiques(Type_Face type) const;

  int_t face_bords_interne_conjuguee(int_t face) const;

  ///
  /// Bords
  ///
  inline int nb_bords() const {  return mes_faces_bord_.nb_bords(); }
  inline Bord_t& bord(int i) {  return mes_faces_bord_(i); }
  inline const Bord_t& bord(int i) const {   return mes_faces_bord_(i); }
  inline Bord_t& bord(const Nom& nom) {  return mes_faces_bord_(nom); }
  inline const Bord_t& bord(const Nom& nom) const {   return mes_faces_bord_(nom); }

  inline Bords_t& faces_bord() {  return mes_faces_bord_; }
  inline const Bords_t& faces_bord() const { return mes_faces_bord_; }

  void fixer_premieres_faces_frontiere();
  void correct_type_of_borders_after_merge();
  void ecrire_noms_bords(Sortie& ) const;

  ///
  /// Bords internes
  ///
  inline Bord_Interne_t& bords_interne(int i) {   return mes_bords_int_(i); }
  inline const Bord_Interne_t& bords_interne(int i) const {  return mes_bords_int_(i); }
  inline Bord_Interne_t& bords_interne(const Nom& nom) {   return mes_bords_int_(nom); }
  inline const Bord_Interne_t& bords_interne(const Nom& nom) const {   return mes_bords_int_(nom); }

  inline Bords_Internes_t& bords_int() { return mes_bords_int_; }
  inline const Bords_Internes_t& bords_int() const { return mes_bords_int_; }


  ///
  /// Groupes_Faces
  ///
  inline int nb_groupes_faces() const {  return mes_groupes_faces_.nb_groupes_faces(); }
  inline Groupe_Faces_t& groupe_faces(int i) {   return mes_groupes_faces_(i); }
  inline const Groupe_Faces_t& groupe_faces(int i) const {   return mes_groupes_faces_(i); }
  inline const Groupe_Faces_t& groupe_faces(const Nom& nom) const {   return mes_groupes_faces_(nom); }
  inline Groupes_Faces_t& groupes_faces() {  return mes_groupes_faces_; }
  inline const Groupes_Faces_t& groupes_faces() const { return mes_groupes_faces_; }

  ///
  /// Frontieres
  ///
  inline const Frontiere_t& frontiere(int i) const;
  inline Frontiere_t& frontiere(int i);
  int rang_frontiere(const Nom& ) const;
  const Frontiere_t& frontiere(const Nom&) const;
  Frontiere_t& frontiere(const Nom&);
  inline int nb_frontieres_internes() const {  return mes_bords_int_.nb_bords_internes(); }
  inline int nb_front_Cl() const  {   return nb_bords() +nb_raccords() + nb_frontieres_internes(); }

  ///
  /// Domaines frontieres
  ///
  inline const LIST(OBS_PTR(Domaine_32_64))& domaines_frontieres() const { return domaines_frontieres_; }
  inline Domaine_32_64& domaine_frontiere(int i) { return domaines_frontieres_(i).valeur(); }

  ///
  /// Raccords
  ///
  inline int nb_raccords() const {  return mes_faces_raccord_.nb_raccords(); }
  inline Raccord_t& raccord(int i) {  return mes_faces_raccord_(i); }
  inline const Raccord_t& raccord(int i) const { return mes_faces_raccord_(i); }
  inline Raccord_t& raccord(const Nom& nom) {  return mes_faces_raccord_(nom); }
  inline const Raccord_t& raccord(const Nom& nom) const {   return mes_faces_raccord_(nom); }

  inline Raccords_t& faces_raccord() { return mes_faces_raccord_; }
  inline const Raccords_t& faces_raccord() const { return mes_faces_raccord_; }

  ///
  /// Joints
  ///
  inline int nb_joints() const { return mes_faces_joint_.nb_joints(); }

  inline Joint_t& joint(int i) {   return mes_faces_joint_(i); }
  inline const Joint_t& joint(int i) const { return mes_faces_joint_(i); }
  inline Joint_t& joint(const Nom& nom) {  return mes_faces_joint_(nom); }
  inline const Joint_t& joint(const Nom& nom) const {   return mes_faces_joint_(nom); }
  inline Joints_t& faces_joint() { return mes_faces_joint_; }
  inline const Joints_t& faces_joint() const {  return mes_faces_joint_; }

  inline Joint_t& joint_of_pe(int);
  inline const Joint_t& joint_of_pe(int) const;
  int comprimer_joints();

  void renum_joint_common_items(const IntVect_t& nums, const int_t elem_offset);


  ///
  /// Periodicity
  ///
  inline const Noms& bords_perio() const { return bords_perio_; }
  inline Noms& bords_perio() { return bords_perio_; }
  void init_renum_perio();
  inline int_t get_renum_som_perio(int_t i) const { return renum_som_perio_[i]; }
  void construire_renum_som_perio(const Conds_lim&, const Domaine_dis_base&);
  inline void set_renum_som_perio(IntTab_t& renum)  {    renum_som_perio_=renum;   };
  const ArrOfInt_t& get_renum_som_perio() const { return renum_som_perio_; }

  //
  // Sous Domaines
  //
  inline int nb_ss_domaines() const { return les_ss_domaines_.size(); }
  inline const Sous_Domaine_t& ss_domaine(int i) const { return les_ss_domaines_[i].valeur(); }
  inline Sous_Domaine_t& ss_domaine(int i) { return les_ss_domaines_[i].valeur(); }
  inline const Sous_Domaine_t& ss_domaine(const Nom& nom) const { return les_ss_domaines_(nom).valeur(); }
  inline Sous_Domaine_t& ss_domaine(const Nom& nom)  { return les_ss_domaines_(nom).valeur(); }
  void add(const Sous_Domaine_t& sd)  { les_ss_domaines_.add(sd); }
  int associer_(Objet_U&) override;  // Associate sous_domaine
  inline const LIST(OBS_PTR(Sous_Domaine_t))& ss_domaines() const { return les_ss_domaines_; }


  ///
  /// Geometric computations
  ///
  inline void calculer_centres_gravite(DoubleTab_t& xp) const;
  void calculer_centres_gravite_aretes(DoubleTab_t& xa) const;
  virtual void calculer_volumes(DoubleVect_t& volumes, DoubleVect_t& inv_volumes) const;
  void calculer_mon_centre_de_gravite(ArrOfDouble& c);
  double volume_total() const;
  inline const ArrOfDouble& cg_moments() const  {  return cg_moments_;  }
  inline ArrOfDouble& cg_moments() {  return cg_moments_;  }
  inline void exporter_mon_centre_de_gravite(ArrOfDouble c)  {  cg_moments_ = c; }

  ///
  /// Lookup methods and mapping arrays
  ///
  SmallArrOfTID_t& chercher_elements(const DoubleTab& pos, SmallArrOfTID_t& elem, int reel=0) const;
  SmallArrOfTID_t& chercher_elements(const DoubleVect& pos, SmallArrOfTID_t& elem, int reel=0) const;
  int_t chercher_elements(double x, double y=0, double z=0,int reel=0) const;
  SmallArrOfTID_t& indice_elements(const IntTab& som, SmallArrOfTID_t& elem, int reel=0) const;

  SmallArrOfTID_t& chercher_sommets(const DoubleTab& pos, SmallArrOfTID_t& som, int reel=0) const;
  int_t chercher_sommets(double x, double y=0, double z=0,int reel=0) const;

  SmallArrOfTID_t& chercher_aretes(const DoubleTab& pos, SmallArrOfTID_t& arr, int reel=0) const;
  void rang_elems_sommet(SmallArrOfTID_t& elems, double x, double y=0, double z=0) const;

  const OctreeRoot_t& construit_octree() const;
  const OctreeRoot_t& construit_octree(int&) const;
  void invalide_octree();

  ///
  /// Various
  ///
  virtual void clear();
  void renum(const IntVect_t& nums);
  void check_domaine();
  void imprimer() const;
  int comprimer();
  void read_former_domaine(Entree& s, bool& read_perio);   // used in Scatter
  void merge_wo_vertices_with(Domaine_32_64& z);
  void fill_from_list(std::list<Domaine_32_64*>& lst);

  ///
  /// MEDCoupling:
  ///
#ifdef MEDCOUPLING_
  inline const MEDCouplingUMesh* get_mc_mesh(bool virt = false) const;
  inline void set_mc_mesh(MCAuto<MEDCouplingUMesh> m) const  { mc_mesh_ = m;    }
  // remapper with other domains
  MEDCouplingRemapper* get_remapper(const Domaine_32_64& other_dom, bool virt=false) const;
  // DEC with other domains
#ifdef MPI_
  OverlapDEC* get_dec(const Domaine_32_64& other_dom, MEDCouplingFieldDouble *dist, MEDCouplingFieldDouble *loc) const;
#endif
#endif
  void build_mc_mesh(bool virt = false) const;
  bool is_mc_mesh_ready() const { return mc_mesh_ready_; }
  // Used in Interprete_geometrique_base for example - method is const, because mc_mesh_ready_ is mutable:
  void set_mc_mesh_ready(bool flag) const { mc_mesh_ready_ = flag; }

  ///
  /// Parallelism and virtual items management
  ///
  inline const ArrOfInt_t& ind_faces_virt_bord() const  { return ind_faces_virt_bord_; }
  void construire_elem_virt_pe_num();
  void construire_elem_virt_pe_num(IntTab_t& elem_virt_pe_num_cpy) const;
  const IntTab_t& elem_virt_pe_num() const  { return elem_virt_pe_num_; }
  virtual void creer_tableau_elements(Array_base&, RESIZE_OPTIONS opt=RESIZE_OPTIONS::COPY_INIT) const;
  virtual const MD_Vector& md_vector_elements() const;
  virtual void creer_tableau_sommets(Array_base&, RESIZE_OPTIONS opt=RESIZE_OPTIONS::COPY_INIT) const;
  virtual const MD_Vector& md_vector_sommets() const { return sommets_.get_md_vector(); }


  ///
  /// Methods only used in 32 bits (i.e. after Scatter)
  ///
  static int identifie_item_unique(IntList& item_possible, DoubleTab& coord_possible, const DoubleVect& coord_ref);
  void init_faces_virt_bord(const MD_Vector& md_vect_faces, MD_Vector& md_vect_faces_bord);
  void creer_aretes();
  void creer_mes_domaines_frontieres(const Domaine_VF& domaine_vf);


protected:
  // Geometric element type of this domain
  OWN_PTR(Elem_geom_base_32_64<_SIZE_>) elem_;

  // Array of vertices
  DoubleTab_t sommets_;
  DoubleTab_t sommets_n; // vertex coordinates at beginning of step, needed for implicit iterations
  // Renumbering array for periodicity
  ArrOfInt_t renum_som_perio_;
  // Description of elements (for multi-element, the array may contain -1 values !!!)
  IntTab_t mes_elems_;
  // Definition of element edges (for each edge, indices of the two vertices)
  //  (this array is not always filled, depending on the discretisation)
  IntTab_t aretes_som_;
  // For each element, indices of its edges in Aretes_som (see Elem_geom_base::get_tab_aretes_sommets_locaux())
  IntTab_t elem_aretes_;
  // For the virtual faces of Domaine_VF, indices of the same face in the boundary face array
  // (see Domaine_32_64<_SIZE_>::init_faces_virt_bord())
  ArrOfInt_t ind_faces_virt_bord_; // contains the indices of the virtual boundary faces // BigArrOfTID
  // For each virtual element i with nb_elem<=i<nb_elem_tot:
  // elem_virt_pe_num_(i-nb_elem,0) = number of the PE owning the element
  // elem_virt_pe_num_(i-nb_elem,1) = local number of this element on the owning PE
  IntTab_t elem_virt_pe_num_;

  // The octree is mutable: it is built on-the-fly when used in const methods
  mutable OWN_PTR(OctreeRoot_t) deriv_octree_;
  ArrOfDouble cg_moments_;  // max dim 3

  // List of references to sub domains
  LIST(OBS_PTR(Sous_Domaine_t)) les_ss_domaines_;

  // Bords, raccords and Bords_Internes form the "faces_frontiere" on which
  //  boundary conditions are defined.
  Bords_t mes_faces_bord_;
  Raccords_t mes_faces_raccord_;
  Bords_Internes_t mes_bords_int_;
  // Groupes_Faces represent the groups of faces read from input files
  Groupes_Faces_t mes_groupes_faces_;
  // Joint faces are the faces shared with other processors (local domain boundaries
  //  that connect to a neighboring processor)
  Joints_t mes_faces_joint_;

  LIST(OBS_PTR(Domaine_32_64)) domaines_frontieres_;

  Noms bords_perio_;  ///< List of periodic boundaries - this is filled by Interprete 'Declarer_bord_perio'

#ifdef MEDCOUPLING_
  ///! MEDCoupling version of the domain:
  mutable MCAuto<MEDCouplingUMesh> mc_mesh_, mc_mesh_virt_;
  // One remapper per distant domain...
  mutable std::map<const Domaine_32_64*, MEDCoupling::MEDCouplingRemapper> rmps;
#ifdef MPI_
  // ... but one DEC per (distant domain, field nature)
  mutable std::map<std::pair<const Domaine_32_64*, MEDCoupling::NatureOfField>, OverlapDEC> decs;
#endif
  mutable bool mc_mesh_ready_ = false, mc_mesh_virt_ready_ = false;
#endif

private:
  void prepare_rmp_with(const Domaine_32_64& other_dom, bool virt) const;
  void prepare_dec_with(const Domaine_32_64& other_dom, MEDCouplingFieldDouble *dist, MEDCouplingFieldDouble *loc) const;

  template<typename _BORD_TYP_>
  void correct_type_single_border_type(std::list<_BORD_TYP_>& list);

  // Cached infos to accelerate chercher_elements():
  mutable DoubleTabs cached_positions_;
  mutable TRUST_Vector<SmallArrOfTID_t> cached_elements_;
};

/*! @brief Sets the element type of the domain using the name passed as parameter.
 *
 * Associates the element type to the domain.
 * @param (Nom& typ) the name of the geometric element type of the domain.
 */
template<typename _SIZE_>
inline void Domaine_32_64<_SIZE_>::typer(const Nom& typ)
{
  const std::string suff = !std::is_same<_SIZE_, int>::value ? "_64" : "";
  Nom typ_32_64 = typ + suff;
  elem_.typer(typ_32_64);  // e.g. 'Point' or 'Point_64'
  elem_->associer_domaine(*this);
}


/*! @brief Returns the number of vertices of the geometric elements that make up the domain.
 *
 * Since all elements of the domain are of the same type they all have the same number of vertices
 * which is the number of vertices of the type of the geometric elements of the domain.
 *
 * @return (int) the number of vertices per element
 */
template<typename _SIZE_>
inline int Domaine_32_64<_SIZE_>::nb_som_elem() const {   return elem_->nb_som(); }

/*! @brief Returns the number of faces of type i of the geometric elements that make up the domain.
 *
 *     Ex: objects of the Prism class have 2 types of faces: triangle or quadrangle.
 *
 * @param (int i) the type of face
 * @return (int) the number of faces of type i of the geometric elements that make up the domain
 */
template<typename _SIZE_>
inline int Domaine_32_64<_SIZE_>::nb_faces_elem(int i) const {  return elem_->nb_faces(i); }

/// Returns the number of boundary faces of the domain (sum of boundaries, connections, and internal boundaries)
template<typename _SIZE_>
inline typename Domaine_32_64<_SIZE_>::int_t Domaine_32_64<_SIZE_>::nb_faces_frontiere() const { return nb_faces_bord() + nb_faces_raccord() + nb_faces_bords_int(); }


/*! @brief Returns the number of special faces of the domain.
 *
 * It is the sum of the numbers of boundaries, connections, internal boundaries and specified face groups
 */
template<typename _SIZE_>
inline typename Domaine_32_64<_SIZE_>::int_t Domaine_32_64<_SIZE_>::nb_faces_specifiques() const { return nb_faces_bord() + nb_faces_raccord() + nb_faces_bords_int() + nb_faces_groupes_faces(); }

/*! @brief Calculates the centers of gravity of the domain elements.
 *
 * @param (DoubleTab_t& xp) the array containing the centers of gravity of the domain elements
 */
template<typename _SIZE_>
inline void Domaine_32_64<_SIZE_>::calculer_centres_gravite(DoubleTab_t& xp) const {  elem_->calculer_centres_gravite(xp); }

/*! @brief Returns the number of boundary faces of the domain of the specified type.
 *
 *     It is the sum of the numbers of boundaries, connections and internal boundaries of the specified type.
 *
 * @param (Type_Face type) a type of face (some geometric elements have multiple types of faces)
 * @return (int) the number of boundary faces of the domain of the specified type
 */
template<typename _SIZE_>
inline typename Domaine_32_64<_SIZE_>::int_t Domaine_32_64<_SIZE_>::nb_faces_frontiere(Type_Face type) const
{
  return
    nb_faces_bord(type) +
    nb_faces_bords_int(type) +
    nb_faces_raccord(type);
}

/*! @brief Returns the number of specific faces of the domain of the specified type.
 *
 *     It is the sum of the numbers of boundaries, connections of internal boundaries and groups of faces of the specified type.
 *
 * @param (Type_Face type) a face type (some geometric elements have several types of faces)
 * @return (int) the number of boundary faces of the domain of the specified type
 */
template<typename _SIZE_>
inline typename Domaine_32_64<_SIZE_>::int_t Domaine_32_64<_SIZE_>::nb_faces_specifiques(Type_Face type) const
{
  return
    nb_faces_bord(type) +
    nb_faces_bords_int(type) +
    nb_faces_raccord(type) +
    nb_faces_groupes_faces(type);
}

template<typename _SIZE_>
inline const typename Domaine_32_64<_SIZE_>::Frontiere_t& Domaine_32_64<_SIZE_>::frontiere(int i) const
{
  int fin=nb_bords();
  if(i<fin)
    return mes_faces_bord_(i);
  i-=fin;
  fin=nb_raccords();
  if(i<fin)
    return mes_faces_raccord_(i).valeur();
  i-=fin;
  fin=nb_frontieres_internes();
  if(i<fin)
    return mes_bords_int_(i);
  i-=fin;
  fin=nb_groupes_faces();
  if(i<fin)
    return mes_groupes_faces_(i);
  assert(0);
  Process::exit();
  return frontiere(i);
}


template<typename _SIZE_>
inline typename Domaine_32_64<_SIZE_>::Frontiere_t& Domaine_32_64<_SIZE_>::frontiere(int i)
{
  int fin=nb_bords();
  if(i<fin)
    return mes_faces_bord_(i);
  i-=fin;
  fin=nb_raccords();
  if(i<fin)
    return mes_faces_raccord_(i).valeur();
  i-=fin;
  fin=nb_frontieres_internes();
  if(i<fin)
    return mes_bords_int_(i);
  i-=fin;
  fin=nb_groupes_faces();
  if(i<fin)
    return mes_groupes_faces_(i);
  assert(0);
  Process::exit();
  return frontiere(i);
}

template<typename _SIZE_>
inline const typename Domaine_32_64<_SIZE_>::Joint_t& Domaine_32_64<_SIZE_>::joint_of_pe(int pe) const
{
  int i;
  for(i=0; i<nb_joints(); i++)
    if(mes_faces_joint_(i).PEvoisin()==pe)
      break;
  return mes_faces_joint_(i);
}

template<typename _SIZE_>
inline typename Domaine_32_64<_SIZE_>::Joint_t& Domaine_32_64<_SIZE_>::joint_of_pe(int pe)
{
  int i;
  for(i=0; i<nb_joints(); i++)
    if(mes_faces_joint_(i).PEvoisin()==pe)
      break;
  return mes_faces_joint_(i);
}


#ifdef MEDCOUPLING_
template<typename _SIZE_>
inline const MEDCouplingUMesh* Domaine_32_64<_SIZE_>::get_mc_mesh(bool virt) const
{
  if (virt ? !mc_mesh_virt_ready_ : !mc_mesh_ready_)
    build_mc_mesh(virt);
  return virt ? mc_mesh_virt_ : mc_mesh_;
}

#endif  // MEDCOUPLING_


using Domaine = Domaine_32_64<int>;
using Domaine_64 = Domaine_32_64<trustIdType>;

#endif
