/****************************************************************************
 * Copyright (c) 2025, CEA
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

#include <Reorder_Mesh.h>
#include <TRUSTTab.h>
#include <Param.h>
#include <Scatter.h>

#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <numeric>
#include <fstream>

Implemente_instanciable(Reorder_Mesh, "Reorder_Mesh", Objet_U);
// XD reorder_mesh objet_u reorder_mesh 1 Reordering option to be used in a discretisation : the geometrical entities (nodes, elems, faces) of a mesh can be reordered to follow a Z-curve (Hilbert or Morton) improving data locality in memory.

namespace // Anonymous namespace
{
// Number of bits to partition a single dimension (x or y) in 2D:
constexpr unsigned int NB_BITS_2D = 32;
// Number of bits to partition a single dimension (x,y or z) in 3D:
// -> NB: nb of bits for the quantization of each coordinate - 21*3=63 < 64 -> will fit in uint64_t
constexpr unsigned int NB_BITS_3D = 21;

// Maximal bitwise values for the binary quantization space:
constexpr uint32_t MAX_VAL_2D = std::numeric_limits<uint32_t>::max(); // binary: 11111...11 (32 times)
constexpr uint32_t MAX_VAL_3D = (1 << NB_BITS_3D) - 1; // binary: 00...0011...11 (eleven 0 and 21 ones)


using cube_pos_t = std::array<int8_t, 3>; // 3 signed integers representing the orientation of the unit x,y,z trihedron
using ar8_t = std::array<uint8_t, 8>;
const ar8_t unit_hilbert = { 0,7,3,4,1,6,2,5 };  // initial indexing of the unit Hilbert 3D curve


/**! Computation of the Morton code for a given point.
 *
 * The Morton code is 'just' the clever bitwise assembly of all the bits from the 3 coordinates.
 * In 2D for example, if x = 0  1  0  1  1
 *                  and  y = 1  0  1  1  0
 * the resutling code is c = 10 01 10 11 01
 * where bits from x and y have been interleaved. Thus two points close in 2D space will also have close (numerical) Morton codes! So clever.
 */
uint64_t mortonCode(uint32_t x, uint32_t y, uint32_t z)
{
  uint64_t result = 0;
  uint64_t xx=x, yy=y, zz=z;
  if (Objet_U::dimension == 2)
    for (unsigned i = 0; i < ::NB_BITS_2D; i++)
      {
        result |= (xx & (1ULL << i)) << i;
        result |= (yy & (1ULL << i)) << (i+1);
      }
  else // dim 3
    for (unsigned i = 0; i < ::NB_BITS_3D; i++)
      {
        result |= (xx & (1ULL << i)) << (2*i);
        result |= (yy & (1ULL << i)) << (2*i+1);
        result |= (zz & (1ULL << i)) << (2*i+2);
      }
  return result;
}

/**! Computation of the Hilbert code. Idea is similar as Morton, but this time points are organized along a Hilbert curve.
 *
 * See wikipedia for example for a nice picture.
 * We proceed by swaping/flipping the bits of each 'slot' (2 bits in 2D, 3 bits in 3D) corresponding to each depth level
 * of the Morton code. See hilbertCode_3D below which explains the spirit in more details.
 */
uint64_t hilbertCode_2D(uint32_t x, uint32_t y)
{
  assert(Objet_U::dimension == 2);

  uint64_t result = 0;
  uint64_t mort = mortonCode(x,y,0.0);

  const uint8_t ud[4] = {0, 3, 1, 2};  // up - down
  const uint8_t rl[4] = {3, 2, 0, 1};  // right - left
  unsigned int orient = 0; // 0: up, 1: right, 2: down, 3:left
  bool flip = false;
  for(int i = NB_BITS_2D-1; i >= 0; i--)
    {
      result <<= 2;
      uint64_t two_bits = (mort >> (2*i)) & 0b11;
      uint8_t new_two_bits = 0b00;
      switch(orient)
        {
        case 0:  // up
          new_two_bits = ud[two_bits];
          break;
        case 1:  // right
          new_two_bits = rl[two_bits];
          break;
        case 2:  // down
          new_two_bits = ud[3-two_bits];
          break;
        case 3:  // left
          new_two_bits = rl[3-two_bits];
          break;
        default:
          throw;
        }
      result |= flip ? (0b11 - new_two_bits) : new_two_bits;
      // Rotate:
      if (new_two_bits == 0b00) { orient = (orient+1) % 4; flip = !flip; }
      if (new_two_bits == 0b11) { orient = (orient+3) % 4; flip = !flip; }
    }
  return result;
}

/**! Computation of the Hilbert code. Idea is similar as Morton, but this time points are organized along a Hilbert curve.
 *
 * See wikipedia for example for a nice picture.
 *
 * We proceed by first computing the Morton code, and then swaping/flipping the bits of each 'slot' (2 bits in 2D, 3 bits in 3D)
 * corresponding to each depth level of the Morton code.
 *
 * Detailed explanation for 3D:
 * Going from one level to the next finer level implies rotating the unit Hilbert curve that we will apply.
 * This rotation corresponds to the rotation of a cube in 3D. In the algorithm below, the notation U, R, F corresponds
 * to the usual notations given when solving the Rubik's cube:
 *    U is for Up (the top face of the cube),
 *    F if for Front (the face pointing towards the viewer)
 *    R is for Right
 *  A prime like in U' or in F' means the rotation is done counter-clockwise.
 *
 * The cube orientation at any time is given by a triplet of numbers between -3 and 3 (0 excluded)
 * At the begining the cube is oriented like this
 *
 *      z
 *      |
 *      |   y
 *      | /
 *      |/________ x
 *
 * which corresponds to the triplet {1,2,3}
 *
 * For example a cube orientation of {3, 2, -1} means that what was the orignal X axis is now the Y axis, what was Y is now Z, and what was Z is now
 * the opposite of X:
 *        y
 *       /
 *      /________ z
 *      |
 *      |
 *      |
 *      x
 *
 * Each time we go to the next finer level, we rotate the cube according to the position we are at in the current Hilbert curve.
 * Sometimes we also need to go through the unit Hilbert curve numbering in reverse.
 *
 * The needed rotation were written based on the nice images found there:
 *     https://pypi.org/project/numpy-hilbert-curve/
 */
uint64_t hilbertCode_3D(uint32_t x, uint32_t y, uint32_t z)
{
  uint64_t result = 0;
  uint64_t mort = mortonCode(x,y,z);

  // Compute p1 . p2 --- mathematical definition : (p1.p2)(x) = p1(p2(x)) = res(x)
  auto compose = [](const cube_pos_t& p1, const cube_pos_t& p2) -> cube_pos_t
  {
    cube_pos_t res;
    for (int i = 0; i < 3; i++)
      {
        uint8_t idx = std::abs(p1[i])-1;
        int8_t sig = p1[i] > 0 ? 1 : -1;
        res[i] = p2[idx]*sig;
      }
    return res;
  };

  auto apply_perm = [](const cube_pos_t& perm, const uint64_t& cod) -> uint8_t
  {
    std::array<uint8_t, 3> cod2;  // x,y,z
    cod2[0] = (uint8_t)(cod & 0b001);
    cod2[1] = (uint8_t)((cod & 0b010) >> 1);
    cod2[2] = (uint8_t)((cod & 0b100) >> 2);

    std::array<uint8_t, 3> res;

    for(int i=0; i < 3; i++)
      {
        int8_t v = perm[i];
        uint8_t idx = std::abs(v)-1;
        int8_t mask = (v < 0) ? 0b1 : 0b0;
        res[idx] = (uint8_t)(cod2[i] ^ mask); // XOR
      }
    return (uint8_t)(res[0] | (res[1] << 1) | (res[2] << 2));
  };

  cube_pos_t curr_perm = {1,2,3};   // 1 represents X axis, 2 the Y axis, 3 the Z axis.

  bool flip = false;  // do we need to scan the unit Hilbert curve in reverse?
  for(int i = NB_BITS_3D-1; i >= 0; i--)  // for each level
    {
      result <<= 3;
      uint64_t three_bits = (mort >> (3*i)) & 0b111;
      // Apply current permutation to pos
      uint64_t new_three_bits = apply_perm(curr_perm, three_bits);

      // Apply Hilbert ordering:
      new_three_bits = unit_hilbert[new_three_bits];

      uint64_t orig3b = new_three_bits; // without the flip!!!

      if(flip)
        new_three_bits = new_three_bits ^ 0b111;

      result |= new_three_bits;
      if (i==0) break;

      // Prepare cube position for next level at depth i+1.
      // The permutation to apply depends on where we are on the current Hilbert curve at level i:
      switch(orig3b)
        {
        case 0: // U R -> {2,-1,3} {1,3,-2} -> {2,3,1}
          curr_perm = compose(curr_perm, {2,3,1});
          break;
        case 7: // U' R -> {-2,1,3} {1,3,-2} -> {-2,3,-1}
          curr_perm = compose(curr_perm, {-2,3,-1});
          break;
        case 3: // F' flip
          curr_perm = compose(curr_perm, {-3,2,1});
          flip = !flip;
          break;
        case 4: // F flip
          curr_perm = compose(curr_perm, {3,2,-1});
          flip = !flip;
          break;
        case 1: // U flip
          curr_perm = compose(curr_perm, {2,-1,3});
          flip = !flip;
          break;
        case 6: // U' flip
          curr_perm = compose(curr_perm, {-2,1,3});
          flip = !flip;
          break;
        case 2: // nothing to do!
          break;
        case 5: // nothing to do!
          break;
        default:
          throw;
        }
    }
  return result;
}

} // end anonymous namespace


Sortie& Reorder_Mesh::printOn(Sortie& os) const { return os; }

// Default is Morton, with all geometrical entities handled.
Entree& Reorder_Mesh::readOn(Entree& is)
{
  int algo;

  Param p(que_suis_je());
  algo = 0;
  p.ajouter("algo", &algo);     // XD_ADD_P dico Z-Curve algorithm to use for geometrical entity renumbering.
  p.dictionnaire("none", 0);    // XD_ADD_DICO No reordering performed (the default).
  p.dictionnaire("morton", 1);  // XD_ADD_DICO Morton scheme for reordering.
  p.dictionnaire("hilbert", 2); // XD_ADD_DICO Hilbert scheme for reordering.

  p.ajouter_flag("dump", &dump_);  // XD_ADD_P flag if set, will dump text files giving the numbering of the various geometrical entities before and after renumbering. Mainly used for debugging.
  // Values of those flags will be set to false by default:
  p.ajouter_flag("no_nodes", &no_nodes_);  // XD_ADD_P flag Whether to avoid node reordering.
  p.ajouter_flag("no_elems", &no_elems_);  // XD_ADD_P flag Whether to avoid element reordering.
  p.ajouter_flag("no_faces", &no_faces_);  // XD_ADD_P flag Whether to avoid face reordering.

  p.lire_avec_accolades(is);

  algo_ = static_cast<Reorder_Algo>(algo);

  return is;
}

/**! Performs the reordering of the nodes and elements of a domain according to the chosen algorithm
 */
template<typename _SIZE_>
void Reorder_Mesh::reorder_domain(Domaine_32_64<_SIZE_>& dom) const
{
  using int_t = _SIZE_;
  using ArrOfInt_t = ArrOfInt_T<_SIZE_>;
  using DoubleTab_t = DoubleTab_T<_SIZE_>;

  // No re-ordering requested or nothing to do:
  if (algo() == Reorder_Algo::None) return;
  if (skip_nodes() && skip_elems()) return;

  Cerr << "****************************************************************" << finl;

  std::string algon = algo() == Reorder_Algo::Morton ? "Morton" : "Hilbert";

  // Renumbering utilities:
  auto renum_tab_indices = [] (auto& tab, const auto& renum)
  {
    assert(tab.nb_dim() == 2);
    auto new_tab(tab);
    for (int_t i = 0; i < tab.dimension(0); i++)
      for (int j = 0; j < tab.dimension(1); j++)
        new_tab(renum(i), j) = tab(i, j);
    tab = new_tab;
  };
  auto renum_tab_values = [] (auto& tab, const auto& renum)
  {
    assert(tab.nb_dim()==2);
    int_t sz_renum = renum.size_array();
    for (int_t i=0; i<tab.dimension_tot(0); i++)  // with virtuals
      for (int j=0; j<tab.dimension(1); j++)
        {
          auto val = tab(i,j);
          if (val < 0 || val >= sz_renum) continue;  // skip uninitialized or virtual values
          tab(i,j) = renum(tab(i,j));
        }
  };
  auto renum_vect_values = [] (auto& vect, const auto& renum)
  {
    for (int_t i=0; i<vect.size_array(); i++)
      vect(i) = renum(vect(i));
  };
  auto compute_how_many = [] (ArrOfInt_t renum) -> int_t
  {
    int_t cnt = 0;
    for(int_t i = 0; i < renum.size_array(); i++)
      if (renum[i] != i) cnt++;
    return cnt;
  };

  // Nodes and Cells renumbering
  ArrOfInt_t renum_nodes, renum_elems;
  int_t nnodes=0, nelems=0;

  if (!skip_nodes())
    {
      Cerr << "[Reordering] mesh *nodes* using " << algon << " scheme ..." << finl;
      compute_renumbering(dom.les_sommets(), renum_nodes);
      nnodes = compute_how_many(renum_nodes);

      dump_to_file(dom.les_sommets(), "reordering_som_before.txt");
      renum_tab_indices(dom.les_sommets(), renum_nodes);
      dump_to_file(dom.les_sommets(), "reordering_som_after.txt");

      // renum_vect_values(renum_som_perio_, renum_nodes);
      renum_tab_values(dom.les_elems(), renum_nodes);

      for (int i=0; i<dom.faces_bord().size(); i++)
        renum_tab_values(dom.bord(i).les_sommets_des_faces(), renum_nodes);
      for (int i=0; i<dom.faces_raccord().size(); i++)
        renum_tab_values(dom.raccord(i)->les_sommets_des_faces(), renum_nodes);
      for (int i=0; i<dom.bords_int().size(); i++)
        renum_tab_values(dom.bords_interne(i).les_sommets_des_faces(), renum_nodes);
      for (int i=0; i<dom.groupes_faces().size(); i++)
        renum_tab_values(dom.groupe_faces(i).les_sommets_des_faces(), renum_nodes);
    }

  if (!skip_elems())
    {
      Cerr << "[Reordering] mesh *elements* using " << algon << " scheme ..." << finl;
      DoubleTab_t xp;

      dom.calculer_centres_gravite(xp);
      // grrrrr .... xp also contains virtuals, but without a proper // structure, hence dimension(0) == dimension_tot(0).
      // And we just want to reorder real elements, not virtual ones (they must remain at the end)
      // We must trim :
      xp.resize(dom.les_elems().dimension(0), xp.dimension_int(1));

      compute_renumbering(xp, renum_elems);
      nelems = compute_how_many(renum_elems);

      dump_to_file(xp, "reordering_elem_before.txt");

      renum_tab_indices(dom.les_elems(), renum_elems);
      renum_tab_indices(xp, renum_elems);

      dump_to_file(xp, "reordering_elem_after.txt");
    }

  if (dom.nb_ss_domaines())
    Process::exit("Reorder not impl for sub-domains");

  for (int i=0; i<dom.domaines_frontieres().size(); i++)
    reorder_domain(dom.domaine_frontiere(i));

  if (Process::nproc()>1)
    {
      // If this piece of code is reached, we are necessarily after a Scatter, and hence with a Domaine_32 object:
      if constexpr (std::is_same<_SIZE_, trustIdType>::value)
        Process::exit("Should never happen!");
      else
        {
          Domaine_32_64<int>& this32 = static_cast<Domaine_32_64<int>&>(dom);
          Cerr << "[Reordering] Updating joints and parallel structures ..." << finl;
          // Local bits of the joints needs renumbering so that they will be the items to be sent when we
          // update the parallel structures:
          if (!skip_nodes())
            {
              for (Joint & j: dom.faces_joint())
                {
                  ArrOfInt& ic = j.set_joint_item(JOINT_ITEM::SOMMET).set_items_communs();
                  renum_vect_values(ic, renum_nodes);
                  // Sort them by increasing order (requirement of Scatter::construire_correspondance_sommets_par_coordonnees ?)
                  ic.ordonne_array();

                  IntTab& soms = j.faces().les_sommets();
                  renum_tab_values(soms, renum_nodes);
                  if (!skip_elems())
                    {
                      IntTab& fv = j.faces().voisins();
                      renum_tab_values(fv, renum_elems);
                    }
                }

              // real only, not virtual
              int nb_som = this32.les_sommets().dimension(0);
              int nb_elem = this32.les_elems().dimension(0);

              // Then update joints by exchanging with neighbor procs to recompute correspondances
              // Logic here is the same as what is done in Raffiner_isotrope_parallele:


              Scatter::uninit_sequential_domain(this32);

              // Remove virtual parts from sommets and elems, they will be recomputed by construire_structures_paralleles()
              this32.les_sommets().resize(nb_som, this32.les_sommets().dimension(1));  // cut virtual
              this32.les_elems().resize(nb_elem, this32.les_elems().dimension(1));  // cut virtual


              Scatter::trier_les_joints(this32.faces_joint());
              Scatter::construire_correspondance_sommets_par_coordonnees(this32, false /* does not allow resize of items_communs */);
              Noms liste_bords_periodiques; // TODO!
              Scatter::construire_structures_paralleles(this32, liste_bords_periodiques);
            }
        }
    }
  else
    Scatter::init_sequential_domain(dom);

  Cerr << "[Reordering] " << nnodes << " nodes and " << nelems << " cells were permuted." << finl;
  Cerr << "****************************************************************" << finl;
}


/**! Compute the renumbering of the given array so that the points follow either a Morton curve or a Hilbert curve.
 *
 * @param renum an array 'renum' so that renum[i] gives the new index of the point orginally numbered 'i' (old-to-new).
 */
template <typename _SIZE_>
void Reorder_Mesh::compute_renumbering(const DoubleTab_T<_SIZE_>& points, ArrOfInt_T<_SIZE_>& renum) const
{
  using int_t = _SIZE_;

  const int dim = Objet_U::dimension;
  assert(points.dimension_int(1) == dim);

  int_t nb_pts = points.dimension(0);  // without virtuals! They must remain unchanged at the end of the arrays.

  // Find bounding box
  std::array<double, 3> minV = { points(0,0), points(0,1), dim == 3 ? points(0,2) : 0.0 };
  std::array<double, 3> maxV = minV;

  for(int_t i = 0; i < nb_pts; i++)
    for (int j = 0; j < dim; j++)
      {
        if (points(i,j) < minV[j]) minV[j] = points(i,j);
        if (points(i,j) > maxV[j]) maxV[j] = points(i,j);
      }

  // Quantization function - this function turns a double value within a range [minVal, maxVal]
  // into an integer such that the whole range of possible integers will be covered between min and maxVal
  auto quantize = [&](double value, double minVal, double maxVal) -> uint32_t
  {
    const uint32_t MAX_VAL = dim == 2 ? MAX_VAL_2D : MAX_VAL_3D;
    double normalized = (value - minVal) / (maxVal - minVal);
    normalized = std::clamp(normalized, 0.0, 1.0);
    return static_cast<uint32_t>(std::round(normalized * MAX_VAL));
  };

  // Prepare renumbering vector
  // new-2-old format: new2old[i] gives the index in the initial mesh of the point now located at index 'i'.
  ArrOfInt_T<_SIZE_> new2old(nb_pts);

  renum.resize(nb_pts);
  std::iota(new2old.begin(), new2old.end(), 0);

  // Compute code vector
  std::vector<uint64_t> codes(nb_pts);
  if (dim == 2)
    for (int_t i = 0; i < nb_pts; i++)
      {
        uint32_t  ax = quantize(points(i,0), minV[0], maxV[0]),
                  ay = quantize(points(i,1), minV[1], maxV[1]);
        codes[i] = (algo_ == Reorder_Algo::Morton) ? mortonCode(ax, ay, 0.0) : hilbertCode_2D(ax, ay);
      }
  else
    for (int_t i = 0; i < nb_pts; i++)
      {
        uint32_t  ax = quantize(points(i,0), minV[0], maxV[0]),
                  ay = quantize(points(i,1), minV[1], maxV[1]),
                  az = quantize(points(i,2), minV[2], maxV[2]);
        codes[i] = (algo_ == Reorder_Algo::Morton) ? mortonCode(ax, ay, az) : hilbertCode_3D(ax, ay, az);
      }

  // Sort indices based on code comparison
  std::sort(new2old.begin(), new2old.end(), [&](int a, int b)
  {
    return codes[a] < codes[b];
  });

  // Translate as old-2-new in renum:
  for (int_t i = 0; i < nb_pts; i++) renum[new2old[i]] = i;
}

template<typename _SIZE_>
void Reorder_Mesh::dump_to_file(const DoubleTab_T<_SIZE_>& points, const std::string& filename) const
{
  if(!dump_) return;

  std::ofstream outFile(filename);
  if (!outFile.is_open())
    throw std::runtime_error("Failed to open file: " + filename);

  for (_SIZE_ i = 0; i < points.dimension(0); i++)
    {
      for (int j = 0; j < points.dimension_int(1); j++)
        outFile << points(i, j) << " ";
      outFile << "\n";
    }

  outFile.close();
}


// Instanciate
template void Reorder_Mesh::reorder_domain(Domaine_32_64<int>& dom) const;
template void Reorder_Mesh::dump_to_file(const DoubleTab_T<int>& points, const std::string& filename) const;
template void Reorder_Mesh::compute_renumbering(const DoubleTab_T<int>& points, ArrOfInt_T<int>& renum) const;

#ifdef INT_is_64_
template void Reorder_Mesh::reorder_domain(Domaine_32_64<trustIdType>& dom) const;
template void Reorder_Mesh::dump_to_file(const DoubleTab_T<trustIdType>& points, const std::string& filename) const;
template void Reorder_Mesh::compute_renumbering(const DoubleTab_T<trustIdType>& points, ArrOfInt_T<trustIdType>& renum) const;
#endif

