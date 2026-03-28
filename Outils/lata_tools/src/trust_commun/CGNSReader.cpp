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

#include <CGNSReader.h>
#include <LataFilter.h>
#include <LataDB.h>

#ifndef WITH_CGNSLOADER

void cgns_reader(const char*, const char*, LataDB&)
{
  Journal() << "CGNS PLUGIN not compiled!" << endl;
  throw;
}
void cgns_to_lata(const char *, const char *, bool , bool , bool , bool )
{
  Journal() << "CGNS PLUGIN not compiled!" << endl;
  throw;
}

#else

#include <cgns++.h>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>

#ifdef HAS_CGNS
  #pragma GCC diagnostic push
  #if __GNUC__ > 5 || __clang_major__ > 10
    #pragma GCC diagnostic ignored "-Wsuggest-override"
  #endif

  #ifdef MPI_
    #include <pcgnslib.h>
  #else
    #include <cgnslib.h>
  #endif

  #pragma GCC diagnostic pop

using std::cerr;
using std::endl;

namespace
{
  static inline float double_to_float_cgns(double x)
  {
    if (!(x < 1.e38 && x > -1.e38))
      {
        cerr << "cgns_reader: Error converting double value " << x << " to float" << endl;
        throw LataDBError(LataDBError::READ_ERROR);
      }
    return (float) x;
  }

  static void cgns_check(int ierr, const char *where)
  {
    if (ierr != CG_OK)
      {
        cerr << "CGNS error in " << where << " : " << cg_get_error() << endl;
        throw LataDBError(LataDBError::READ_ERROR);
      }
  }

  static Nom build_filename_in_master_file(const char *data_filename)
  {
    if (!data_filename)
      return "DATA_NOT_WRITTEN";
    return data_filename;
  }

  static const char* zone_type_to_string(ZoneType_t t)
  {
    switch(t)
      {
      case Structured:
        return "Structured";
      case Unstructured:
        return "Unstructured";
      default:
        return "UnknownZoneType";
      }
  }

  static const char* grid_location_to_string(GridLocation_t loc)
  {
    switch(loc)
      {
      case Vertex:
        return "Vertex";
      case CellCenter:
        return "CellCenter";
      default:
        return "UnknownGridLocation";
      }
  }

  static const char* data_type_to_string(DataType_t t)
  {
    switch(t)
      {
      case Integer:
        return "Integer";
      case LongInteger:
        return "LongInteger";
      case RealSingle:
        return "RealSingle";
      case RealDouble:
        return "RealDouble";
      case Character:
        return "Character";
      default:
        return "UnknownDataType";
      }
  }

  static const char* elem_type_to_string_dbg(ElementType_t t)
  {
    switch(t)
      {
      case NODE:
        return "POINT";
      case BAR_2:
        return "SEGMENT";
      case TRI_3:
        return "TRIANGLE";
      case QUAD_4:
        return "RECTANGLE";
      case TETRA_4:
        return "TETRAEDRE";
      case HEXA_8:
        return "HEXAEDRE";
      case NGON_n:
        return "POLYGONE";
      case NFACE_n:
        return "POLYEDRE";
      default:
        return "UnknownElementType";
      }
  }

  static int elem_dimension_from_type(ElementType_t t)
  {
    switch(t)
      {
      case NODE:
//        return 0;
      case BAR_2:
        return 1;
      case TRI_3:
      case QUAD_4:
      case NGON_n:
        return 2;
      case TETRA_4:
      case PYRA_5:
      case PENTA_6:
      case HEXA_8:
      case NFACE_n:
        return 3;
      default:
        return -1;
      }
  }

  static const char* map_cgns_elemtype_to_lata(ElementType_t t, int &nb_comp)
  {
    switch(t)
      {
      case NODE:
        nb_comp = 1;
        return "POINT";
      case BAR_2:
        nb_comp = 2;
        return "SEGMENT";
      case TRI_3:
        nb_comp = 3;
        return "TRIANGLE";
      case QUAD_4:
        nb_comp = 4;
        return "RECTANGLE";
      case TETRA_4:
        nb_comp = 4;
        return "TETRAEDRE";
      case HEXA_8:
        nb_comp = 8;
        return "HEXAEDRE";
      case PENTA_6:
        nb_comp = 6;
        return "PRISM6";
      case NGON_n:
        nb_comp = -1;
        return "POLYGONE";
      case NFACE_n:
        nb_comp = -1;
        return "POLYEDRE";
      default:
        nb_comp = 0;
        return nullptr;
      }
  }

  static void dump_zone_sections(int fn, int ibase, int izone, const char *zonename)
  {
    int nsections = 0;
    cgns_check(cg_nsections(fn, ibase, izone, &nsections), "cg_nsections(debug)");

    Journal(2) << "cgns_reader: zone=" << zonename << " nsections=" << nsections << endl;

    for (int isec = 1; isec <= nsections; isec++)
      {
        char secname[33];
        ElementType_t elem_type;
        cgsize_t start = 0, end = 0;
        int nbndry = 0;
        int parent_flag = 0;
        cgns_check(cg_section_read(fn, ibase, izone, isec, secname, &elem_type, &start, &end, &nbndry, &parent_flag), "cg_section_read(debug)");

        Journal(2) << "cgns_reader:   section[" << isec << "]" << " name=" << secname << " type=" << elem_type_to_string_dbg(elem_type) << " elem_dim=" << elem_dimension_from_type(elem_type)
            << " start=" << start << " end=" << end << " nbndry=" << nbndry << " parent_flag=" << parent_flag << endl;
      }
  }

  static void dump_zone_solutions(int fn, int ibase, int izone, const char *zonename)
  {
    int nsols = 0;
    cgns_check(cg_nsols(fn, ibase, izone, &nsols), "cg_nsols(debug)");

    Journal(2) << "cgns_reader: zone=" << zonename << " nsols=" << nsols << endl;

    for (int isol = 1; isol <= nsols; isol++)
      {
        char solname[33];
        GridLocation_t location = Vertex;
        cgns_check(cg_sol_info(fn, ibase, izone, isol, solname, &location), "cg_sol_info(debug)");

        Journal(2) << "cgns_reader:   sol[" << isol << "]" << " name=" << solname << " location=" << grid_location_to_string(location) << endl;

        int nfields = 0;
        cgns_check(cg_nfields(fn, ibase, izone, isol, &nfields), "cg_nfields(debug)");

        Journal(2) << "cgns_reader:   sol[" << isol << "] nfields=" << nfields << endl;

        for (int ifield = 1; ifield <= nfields; ifield++)
          {
            DataType_t dtype;
            char field_name[33];
            cgns_check(cg_field_info(fn, ibase, izone, isol, ifield, &dtype, field_name), "cg_field_info(debug)");

            Journal(2) << "cgns_reader:     field[" << ifield << "]" << " name=" << field_name << " dtype=" << data_type_to_string(dtype) << endl;
          }
      }
  }

  static int get_max_number_of_solutions(int fn)
  {
    int nbases = 0;
    cgns_check(cg_nbases(fn, &nbases), "cg_nbases(get_max_number_of_solutions)");

    int max_nsols = 0;

    for (int ibase = 1; ibase <= nbases; ibase++)
      {
        int nzones = 0;
        cgns_check(cg_nzones(fn, ibase, &nzones), "cg_nzones(get_max_number_of_solutions)");

        for (int izone = 1; izone <= nzones; izone++)
          {
            int nsols = 0;
            cgns_check(cg_nsols(fn, ibase, izone, &nsols), "cg_nsols(get_max_number_of_solutions)");
            max_nsols = std::max(max_nsols, nsols);
          }
      }

    return max_nsols;
  }

  static std::vector<double> read_cgns_time_values(int fn, int ibase)
  {
    std::vector<double> times;

    int nsteps = 0;
    char bitername[33] = "";
    if (cg_biter_read(fn, ibase, bitername, &nsteps) != CG_OK || nsteps <= 0)
      return times; // pas de BaseIterativeData_t

    // XXX : Aller sous BaseIterativeData_t puis lire TimeValues
    cgns_check(cg_goto(fn, ibase, "BaseIterativeData_t", 1, "end"), "cg_goto(BaseIterativeData_t)");

    int narrays = 0;
    cgns_check(cg_narrays(&narrays), "cg_narrays");

    for (int ia = 1; ia <= narrays; ia++)
      {
        char array_name[33];
        DataType_t dtype;
        int dim;
        cgsize_t dims[12] = { 0 };

        cgns_check(cg_array_info(ia, array_name, &dtype, &dim, dims), "cg_array_info");

        if (std::string(array_name) == "TimeValues")
          {
            times.resize((size_t) nsteps);

            if (dtype == RealDouble)
              cgns_check(cg_array_read_as(ia, RealDouble, times.data()), "cg_array_read_as(TimeValues)");
            else if (dtype == RealSingle)
              {
                std::vector<float> tmp((size_t) nsteps);
                cgns_check(cg_array_read_as(ia, RealSingle, tmp.data()), "cg_array_read_as(TimeValues)");
                for (int i = 0; i < nsteps; i++)
                  times[i] = tmp[i];
              }
            else
              {
                cerr << "cgns_reader: TimeValues has unsupported type" << endl;
                throw LataDBError(LataDBError::READ_ERROR);
              }

            return times;
          }
      }

    return times; // NumberOfSteps existe mais pas TimeValues
  }

  struct ZonePartInfo
  {
    int izone = -1;
    int rank = -1;
    std::string zonename;

    bool has_main_section = false;

    trustIdType nb_nodes = 0;
    trustIdType nb_cells = 0;
    trustIdType nb_cells_effective = 0;
  };

  struct ZoneConnectivityData
  {
    Nom lata_elem_type;
    BigTIDTab elements; // elem -> som
    BigTIDTab faces; // face -> som
    BigTIDTab elem_faces; // elem -> faces
  };

  // XXX just declare here ;)
  static bool try_choose_main_section(int fn, int ibase, int izone, int cell_dim,
                                      const Nom &geom_name, int &best_sec_out);

  static bool is_all_digits(const std::string &s)
  {
    if (s.empty())
      return false;
    for (char c : s)
      if (c < '0' || c > '9')
        return false;
    return true;
  }

  static bool split_parallel_zone_name(const std::string &basename, const std::string &zonename, int &rank)
  {
    {
      const std::string prefix = basename + "_";
      if (zonename.rfind(prefix, 0) == 0)
        {
          const std::string suffix = zonename.substr(prefix.size());
          if (is_all_digits(suffix))
            {
              rank = std::atoi(suffix.c_str());
              return true;
            }
        }
    }

    {
      const std::string generic_prefix = "Zone_";
      if (zonename.rfind(generic_prefix, 0) == 0)
        {
          const std::string suffix = zonename.substr(generic_prefix.size());
          if (is_all_digits(suffix))
            {
              rank = std::atoi(suffix.c_str());
              return true;
            }
        }
    }

    return false;
  }

  static std::string strip_parallel_suffix_if_any(const std::string& name)
  {
    const size_t p = name.rfind('_');
    if (p == std::string::npos || p + 1 >= name.size())
      return name;

    for (size_t i = p + 1; i < name.size(); i++)
      if (name[i] < '0' || name[i] > '9')
        return name;

    return name.substr(0, p);
  }

  // XXX just declare
  static trustIdType read_zone_connectivity(int fn, int ibase, int izone, int cell_dim, const Nom &geom_name, trustIdType nb_nodes, ZoneConnectivityData &zc);

  static std::vector<ZonePartInfo> collect_parallel_zone_parts(int fn, int ibase, const std::string &basename, int cell_dim)
  {
    int nzones = 0;
    cgns_check(cg_nzones(fn, ibase, &nzones), "cg_nzones(collect_parallel_zone_parts)");

    std::vector<ZonePartInfo> parts;
    parts.reserve((size_t) nzones);

    for (int izone = 1; izone <= nzones; izone++)
      {
        char zonename_c[33];
        cgsize_t size[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        cgns_check(cg_zone_read(fn, ibase, izone, zonename_c, size), "cg_zone_read(collect_parallel_zone_parts)");

        ZoneType_t ztype;
        cgns_check(cg_zone_type(fn, ibase, izone, &ztype), "cg_zone_type(collect_parallel_zone_parts)");
        if (ztype != Unstructured)
          continue;

        const std::string zonename(zonename_c);
        int rank = -1;
        if (!split_parallel_zone_name(basename, zonename, rank))
          continue;

        ZonePartInfo p;
        p.izone = izone;
        p.rank = rank;
        p.zonename = zonename;

        int isec = -1;
        p.has_main_section = try_choose_main_section(fn, ibase, izone, cell_dim, zonename.c_str(), isec);

        if (p.has_main_section)
          {
            p.nb_nodes = (trustIdType) size[0];
            p.nb_cells = (trustIdType) size[1];

            ZoneConnectivityData zc_tmp;
            p.nb_cells_effective = read_zone_connectivity(fn, ibase, izone, cell_dim, zonename.c_str(), p.nb_nodes, zc_tmp);
          }
        else
          {
            p.nb_nodes = 0;
            p.nb_cells = 0;
            p.nb_cells_effective = 0;
          }

        parts.push_back(p);
      }

    std::sort(parts.begin(), parts.end(), [](const ZonePartInfo &a, const ZonePartInfo &b)
      {
        return a.rank < b.rank;
      });

    for (size_t i = 1; i < parts.size(); i++)
      {
        if (parts[i].rank == parts[i - 1].rank)
          {
            cerr << "cgns_reader: duplicate parallel rank " << parts[i].rank << " in base " << basename << endl;
            throw LataDBError(LataDBError::READ_ERROR);
          }
      }

    return parts;
  }

  static void read_zone_coordinates(int fn, int ibase, int izone, int phys_dim, trustIdType nb_nodes, BigFloatTab &nodes)
  {
    nodes.resize(nb_nodes, phys_dim);

    std::vector<double> cx((size_t) nb_nodes), cy, cz;
    if (phys_dim >= 2)
      cy.resize((size_t) nb_nodes);
    if (phys_dim >= 3)
      cz.resize((size_t) nb_nodes);

    cgsize_t rmin[3] = { 1, 1, 1 };
    cgsize_t rmax[3] = { (cgsize_t) nb_nodes, 1, 1 };

    cgns_check(cg_coord_read(fn, ibase, izone, "CoordinateX", RealDouble, rmin, rmax, cx.data()), "cg_coord_read(CoordinateX)");

    if (phys_dim >= 2)
      cgns_check(cg_coord_read(fn, ibase, izone, "CoordinateY", RealDouble, rmin, rmax, cy.data()), "cg_coord_read(CoordinateY)");

    if (phys_dim >= 3)
      cgns_check(cg_coord_read(fn, ibase, izone, "CoordinateZ", RealDouble, rmin, rmax, cz.data()), "cg_coord_read(CoordinateZ)");

    for (trustIdType i = 0; i < nb_nodes; i++)
      {
        nodes(i, 0) = double_to_float_cgns(cx[(size_t) i]);
        if (phys_dim >= 2)
          nodes(i, 1) = double_to_float_cgns(cy[(size_t) i]);
        if (phys_dim >= 3)
          nodes(i, 2) = double_to_float_cgns(cz[(size_t) i]);
      }
  }

  static void add_sommets_field(const Nom &filename_in_master_file, const Nom &geom_name, int tstep, Size_t &file_offset, trustIdType nb_nodes, int phys_dim, const BigFloatTab &nodes,
                                const char *data_filename, LataDB &lata_db)
  {
    LataDBField sommets;
    sommets.name_ = "SOMMETS";
    sommets.timestep_ = tstep;
    sommets.filename_ = filename_in_master_file;
    sommets.geometry_ = geom_name;
    sommets.uname_ = Field_UName(geom_name, "SOMMETS", "");
    sommets.nb_comp_ = phys_dim;
    sommets.size_ = nb_nodes;
    sommets.nature_ = LataDBField::VECTOR;
    sommets.datatype_ = lata_db.default_type_float();
    sommets.datatype_.file_offset_ = file_offset++;

    lata_db.add_field(sommets);
    if (data_filename)
      lata_db.write_data(tstep, sommets.uname_, nodes);
  }

  static void add_elements_field(const Nom &filename_in_master_file, const Nom &geom_name, int tstep, Size_t &file_offset,
                                 const Nom &lata_elem_type, trustIdType nb_elem, int nb_comp, const BigTIDTab &elems,
                                 const char *data_filename, LataDB &lata_db)
  {
    lata_db.set_elemtype(tstep, geom_name, lata_elem_type);

    LataDBField elements;
    elements.name_ = "ELEMENTS";
    elements.timestep_ = tstep;
    elements.filename_ = filename_in_master_file;
    elements.geometry_ = geom_name;
    elements.uname_ = Field_UName(geom_name, "ELEMENTS", "");
    elements.nb_comp_ = nb_comp;
    elements.size_ = nb_elem;
    elements.nature_ = LataDBField::SCALAR;
    elements.datatype_ = lata_db.default_type_int_;
    elements.datatype_.array_index_ = LataDBDataType::C_INDEXING;
    elements.datatype_.file_offset_ = file_offset++;

    lata_db.add_field(elements);
    if (data_filename)
      lata_db.write_data(tstep, elements.uname_, elems);
  }

  static void add_faces_field(const Nom &filename_in_master_file, const Nom &geom_name, int tstep, Size_t &file_offset,
                              trustIdType nb_faces, int nb_comp, const BigTIDTab &faces,
                              const char *data_filename, LataDB &lata_db)
  {
    LataDBField field;
    field.name_ = "FACES";
    field.timestep_ = tstep;
    field.filename_ = filename_in_master_file;
    field.geometry_ = geom_name;
    field.localisation_ = "FACES";
    field.uname_ = Field_UName(geom_name, "FACES", "FACES");
    field.nb_comp_ = nb_comp;
    field.size_ = nb_faces;
    field.nature_ = LataDBField::SCALAR;
    field.datatype_ = lata_db.default_type_int_;
    field.datatype_.array_index_ = LataDBDataType::C_INDEXING;
    field.datatype_.file_offset_ = file_offset++;

    lata_db.add_field(field);
    if (data_filename)
      lata_db.write_data(tstep, field.uname_, faces);
  }

  static void add_elem_faces_field(const Nom &filename_in_master_file, const Nom &geom_name, int tstep, Size_t &file_offset,
                                   trustIdType nb_elem, int nb_comp, const BigTIDTab &elem_faces,
                                   const char *data_filename, LataDB &lata_db)
  {
    LataDBField field;
    field.name_ = "ELEM_FACES";
    field.timestep_ = tstep;
    field.filename_ = filename_in_master_file;
    field.geometry_ = geom_name;
    field.localisation_ = "ELEM";
    field.uname_ = Field_UName(geom_name, "ELEM_FACES", "ELEM");
    field.nb_comp_ = nb_comp;
    field.size_ = nb_elem;
    field.nature_ = LataDBField::SCALAR;
    field.datatype_ = lata_db.default_type_int_;
    field.datatype_.array_index_ = LataDBDataType::C_INDEXING;
    field.datatype_.file_offset_ = file_offset++;

    lata_db.add_field(field);
    if (data_filename)
      lata_db.write_data(tstep, field.uname_, elem_faces);
  }

  static void add_scalar_field_to_lata(const Nom &filename_in_master_file, const Nom &geom_name, const Nom &field_name,
                                       const Nom &lata_loc, int tstep_field, Size_t &file_offset,
                                       trustIdType field_size, const BigFloatTab &tab,
                                       const char *data_filename, LataDB &lata_db)
  {
    LataDBField field;
    field.name_ = field_name;
    field.timestep_ = tstep_field;
    field.filename_ = filename_in_master_file;
    field.geometry_ = geom_name;
    field.localisation_ = lata_loc;
    field.uname_ = Field_UName(geom_name, field.name_, lata_loc);
    field.nb_comp_ = 1;
    field.size_ = field_size;
    field.nature_ = LataDBField::SCALAR;
    field.datatype_ = lata_db.default_type_float();
    field.datatype_.file_offset_ = file_offset++;

    lata_db.add_field(field);
    if (data_filename)
      lata_db.write_data(tstep_field, field.uname_, tab);
  }

  static bool try_choose_main_section(int fn, int ibase, int izone, int cell_dim, const Nom &geom_name, int &best_sec_out)
  {
    int nsections = 0;
    cgns_check(cg_nsections(fn, ibase, izone, &nsections), "cg_nsections");

    if (nsections <= 0)
      {
        Journal(2) << "cgns_reader: zone " << geom_name << " has no Elements_t section -> treated as empty" << endl;
        best_sec_out = -1;
        return false;
      }

    int best_sec = -1;
    cgsize_t best_count = -1;

    for (int isec = 1; isec <= nsections; isec++)
      {
        char secname[33];
        ElementType_t elem_type;
        cgsize_t start = 0, end = 0;
        int nbndry = 0;
        int parent_flag = 0;

        cgns_check(cg_section_read(fn, ibase, izone, isec, secname, &elem_type, &start, &end, &nbndry, &parent_flag), "cg_section_read(try_choose_main_section)");

        int nb_comp = 0;
        const char *lata_elem_type = map_cgns_elemtype_to_lata(elem_type, nb_comp);
        const int edim = elem_dimension_from_type(elem_type);
        const cgsize_t count = end - start + 1;

        if (lata_elem_type && edim == cell_dim)
          {
            if (count > best_count)
              {
                best_count = count;
                best_sec = isec;
              }
          }
      }

    if (best_sec < 0)
      {
        Journal(2) << "cgns_reader: zone " << geom_name << " has no supported main section -> treated as empty" << endl;
        best_sec_out = -1;
        return false;
      }

    Journal(2) << "cgns_reader: zone=" << geom_name << " selected main section=" << best_sec << endl;

    best_sec_out = best_sec;
    return true;
  }

  static void read_scalar_field(int fn, int ibase, int izone, int isol, const char *field_name, trustIdType nitems, BigFloatTab &tab)
  {
    tab.resize(nitems, 1);

    std::vector<double> values((size_t) nitems);
    cgsize_t rmin[3] = { 1, 1, 1 };
    cgsize_t rmax[3] = { (cgsize_t) nitems, 1, 1 };

    cgns_check(cg_field_read(fn, ibase, izone, isol, field_name, RealDouble, rmin, rmax, values.data()), "cg_field_read");

    for (trustIdType i = 0; i < nitems; i++)
      tab(i, 0) = double_to_float_cgns(values[(size_t) i]);
  }

  static void read_zone_solution_fields_and_add(int fn, int ibase, int izone, int isol, const Nom &filename_in_master_file, const Nom &geom_name, int tstep_field, Size_t &file_offset,
                                                trustIdType nb_nodes, trustIdType nb_elem, const char *data_filename, LataDB &lata_db)
  {
    char solname[33];
    GridLocation_t location = Vertex;
    cgns_check(cg_sol_info(fn, ibase, izone, isol, solname, &location), "cg_sol_info");

    Nom lata_loc;
    trustIdType field_size = -1;

    if (location == Vertex)
      {
        lata_loc = "SOM";
        field_size = nb_nodes;
      }
    else if (location == CellCenter)
      {
        lata_loc = "ELEM";
        field_size = nb_elem;
      }
    else
      {
        cerr << "cgns_reader: skipping solution " << solname << " in zone " << geom_name << " unsupported GridLocation=" << grid_location_to_string(location) << endl;
        return;
      }

    int nfields = 0;
    cgns_check(cg_nfields(fn, ibase, izone, isol, &nfields), "cg_nfields");

    Journal(2) << "cgns_reader: reading zone=" << geom_name << " solution=" << solname << " isol=" << isol << " -> timestep=" << tstep_field << " location=" << lata_loc << " nfields=" << nfields
        << endl;

    for (int ifield = 1; ifield <= nfields; ifield++)
      {
        DataType_t dtype;
        char field_name[33];
        cgns_check(cg_field_info(fn, ibase, izone, isol, ifield, &dtype, field_name), "cg_field_info");

        BigFloatTab tab;
        read_scalar_field(fn, ibase, izone, isol, field_name, field_size, tab);

        Journal(2) << "cgns_reader:   field name=" << field_name << " dtype=" << data_type_to_string(dtype) << " size=" << field_size;
        if (field_size > 0)
          Journal(2) << " first=" << tab(0, 0) << " last=" << tab(field_size - 1, 0);
        Journal(2) << endl;

        LataDBField field;
        field.name_ = field_name;
        field.timestep_ = tstep_field;
        field.filename_ = filename_in_master_file;
        field.geometry_ = geom_name;
        field.localisation_ = lata_loc;
        field.uname_ = Field_UName(geom_name, field.name_, lata_loc);
        field.nb_comp_ = 1;
        field.size_ = field_size;
        field.nature_ = LataDBField::SCALAR;
        field.datatype_ = lata_db.default_type_float();
        field.datatype_.file_offset_ = file_offset++;

        lata_db.add_field(field);
        if (data_filename)
          lata_db.write_data(tstep_field, field.uname_, tab);
      }
  }

  static void read_merged_zone_solution_fields_and_add(int fn, int ibase,
                                                       const std::vector<ZonePartInfo> &parts,
                                                       const Nom &filename_in_master_file,
                                                       const Nom &geom_name,
                                                       int tstep_field,
                                                       Size_t &file_offset,
                                                       const std::vector<trustIdType> &node_offsets,
                                                       const std::vector<trustIdType> &elem_offsets,
                                                       trustIdType total_nb_nodes,
                                                       trustIdType total_nb_elem,
                                                       const char *data_filename,
                                                       LataDB &lata_db)
  {
    if (parts.empty())
      return;

    // XXX Premiere zone non vide ==> reference
    int iref = -1;
    for (int ip = 0; ip < (int)parts.size(); ip++)
      if (parts[ip].has_main_section)
        {
          iref = ip;
          break;
        }

    if (iref < 0)
      return;

    int nsols0 = 0;
    cgns_check(cg_nsols(fn, ibase, parts[iref].izone, &nsols0), "cg_nsols(merged)");

    if (tstep_field < 1 || tstep_field > nsols0)
      return;

    char solname0[33];
    GridLocation_t location0 = Vertex;
    cgns_check(cg_sol_info(fn, ibase, parts[iref].izone, tstep_field, solname0, &location0), "cg_sol_info(merged)");

    Nom lata_loc;
    trustIdType merged_size = -1;
    bool on_nodes = false;

    if (location0 == Vertex)
      {
        lata_loc = "SOM";
        merged_size = total_nb_nodes;
        on_nodes = true;
      }
    else if (location0 == CellCenter)
      {
        lata_loc = "ELEM";
        merged_size = total_nb_elem;
        on_nodes = false;
      }
    else
      {
        cerr << "cgns_reader: skipping merged solution " << solname0
             << " in geom " << geom_name
             << " unsupported GridLocation=" << grid_location_to_string(location0) << endl;
        return;
      }

    if (merged_size <= 0)
      return;

    int nfields0 = 0;
    cgns_check(cg_nfields(fn, ibase, parts[iref].izone, tstep_field, &nfields0), "cg_nfields(merged)");

    for (int ifield = 1; ifield <= nfields0; ifield++)
      {
        DataType_t dtype0;
        char field_name0[33];
        cgns_check(cg_field_info(fn, ibase, parts[iref].izone, tstep_field, ifield, &dtype0, field_name0), "cg_field_info(merged)");

        BigFloatTab merged_tab;
        merged_tab.resize(merged_size, 1);

        for (size_t ip = 0; ip < parts.size(); ip++)
          {
            if (!parts[ip].has_main_section)
              continue;

            const int izone = parts[ip].izone;

            int nsols = 0;
            cgns_check(cg_nsols(fn, ibase, izone, &nsols), "cg_nsols(merged.part)");
            if (tstep_field > nsols)
              {
                cerr << "cgns_reader: zone " << parts[ip].zonename
                     << " has fewer solutions than expected for merged geometry " << geom_name << endl;
                throw LataDBError(LataDBError::READ_ERROR);
              }

            char solname[33];
            GridLocation_t location = Vertex;
            cgns_check(cg_sol_info(fn, ibase, izone, tstep_field, solname, &location), "cg_sol_info(merged.part)");
            if (location != location0)
              {
                cerr << "cgns_reader: inconsistent GridLocation for merged geometry " << geom_name
                     << " field " << field_name0 << endl;
                throw LataDBError(LataDBError::READ_ERROR);
              }

            int nfields = 0;
            cgns_check(cg_nfields(fn, ibase, izone, tstep_field, &nfields), "cg_nfields(merged.part)");
            if (ifield > nfields)
              {
                cerr << "cgns_reader: inconsistent number of fields for merged geometry " << geom_name << endl;
                throw LataDBError(LataDBError::READ_ERROR);
              }

            DataType_t dtype;
            char field_name[33];
            cgns_check(cg_field_info(fn, ibase, izone, tstep_field, ifield, &dtype, field_name), "cg_field_info(merged.part)");

            if (std::string(field_name) != std::string(field_name0))
              {
                cerr << "cgns_reader: inconsistent field ordering while merging geometry " << geom_name
                     << " expected=" << field_name0 << " got=" << field_name << endl;
                throw LataDBError(LataDBError::READ_ERROR);
              }

            const trustIdType local_size = on_nodes ? parts[ip].nb_nodes : parts[ip].nb_cells_effective;
            const trustIdType offset = on_nodes ? node_offsets[ip] : elem_offsets[ip];

            if (local_size <= 0)
              continue;

            BigFloatTab local_tab;
            read_scalar_field(fn, ibase, izone, tstep_field, field_name0, local_size, local_tab);

            for (trustIdType i = 0; i < local_size; i++)
              merged_tab(offset + i, 0) = local_tab(i, 0);
          }

        add_scalar_field_to_lata(filename_in_master_file, geom_name, field_name0, lata_loc,
                                 tstep_field, file_offset, merged_size, merged_tab,
                                 data_filename, lata_db);
      }
  }

  static bool is_generic_parallel_zone_name(const std::string& zonename)
  {
    const std::string generic_prefix = "Zone_";
    if (zonename.rfind(generic_prefix, 0) != 0)
      return false;

    const std::string suffix = zonename.substr(generic_prefix.size());
    return is_all_digits(suffix);
  }

  static std::string build_geometry_name_from_base_and_zone(const std::string& basename,
                                                            const std::string& zonename)
  {
    // Cas master/base parallele avec Zone_0000, Zone_0002, ... => le vrai nom logique est le nom de base
    if (is_generic_parallel_zone_name(zonename))
      return strip_parallel_suffix_if_any(basename);

    // Cas classique: dom_boundaries_xxx_0000 => dom_boundaries_xxx
    return strip_parallel_suffix_if_any(zonename);
  }

  // *******************
  // *** Pour polys ...
  // *******************
  static bool zone_has_section_type(int fn, int ibase, int izone, ElementType_t wanted_type, int *section_id = nullptr)
  {
    int nsections = 0;
    cgns_check(cg_nsections(fn, ibase, izone, &nsections), "cg_nsections(zone_has_section_type)");

    for (int isec = 1; isec <= nsections; isec++)
      {
        char secname[33];
        ElementType_t elem_type;
        cgsize_t start = 0, end = 0;
        int nbndry = 0;
        int parent_flag = 0;

        cgns_check(cg_section_read(fn, ibase, izone, isec, secname, &elem_type, &start, &end, &nbndry, &parent_flag), "cg_section_read(zone_has_section_type)");

        if (elem_type == wanted_type)
          {
            if (section_id)
              *section_id = isec;
            return true;
          }
      }
    return false;
  }

  static void read_cgns_poly_raw(int fn, int ibase, int izone, int isec, std::vector<cgsize_t> &conn, std::vector<cgsize_t> &offsets,
                                 trustIdType &nb_elem, ElementType_t expected_type)
  {
    char secname[33];
    ElementType_t elem_type;
    cgsize_t start = 0, end = 0;
    int nbndry = 0;
    int parent_flag = 0;

    cgns_check(cg_section_read(fn, ibase, izone, isec, secname, &elem_type, &start, &end, &nbndry, &parent_flag), "cg_section_read(read_cgns_poly_raw)");

    if (elem_type != expected_type)
      {
        cerr << "cgns_reader: wrong poly section type in section " << secname << endl;
        throw LataDBError(LataDBError::READ_ERROR);
      }

    nb_elem = (trustIdType) (end - start + 1);

    cgsize_t total_conn_size = 0;
    cgns_check(cg_ElementDataSize(fn, ibase, izone, isec, &total_conn_size), "cg_ElementDataSize(read_cgns_poly_raw)");

    conn.resize((size_t) total_conn_size);
    offsets.resize((size_t) nb_elem + 1);

    cgns_check(cg_poly_elements_read(fn, ibase, izone, isec, conn.data(), offsets.data(), nullptr), "cg_poly_elements_read");
  }

  static BigTIDTab build_padded_table_from_cgns_poly(const std::vector<cgsize_t>& offsets,
                                                     const std::vector<cgsize_t>& conn,
                                                     bool take_abs = false)
  {
    const trustIdType nb_elem = (trustIdType)offsets.size() - 1;

    trustIdType max_deg = 0;
    for (trustIdType i = 0; i < nb_elem; i++)
      {
        const trustIdType n = (trustIdType)(offsets[(size_t)i + 1] - offsets[(size_t)i]);
        if (n > max_deg) max_deg = n;
      }

    BigTIDTab tab;
    tab.resize(nb_elem, (int)max_deg);
    tab = -1;

    for (trustIdType i = 0; i < nb_elem; i++)
      {
        const trustIdType begin = (trustIdType)offsets[(size_t)i];
        const trustIdType end   = (trustIdType)offsets[(size_t)i + 1];

        for (trustIdType j = begin; j < end; j++)
          {
            trustIdType v = (trustIdType)conn[(size_t)j];
            if (take_abs && v < 0) v = -v;
            tab(i, (int)(j - begin)) = v - 1; // CGNS -> C indexing
          }
      }

    return tab;
  }

  static BigTIDTab build_polyhedron_elements_from_elem_faces_and_faces(const BigTIDTab& elem_faces,
                                                                       const BigTIDTab& faces)
  {
    const trustIdType nb_elem = elem_faces.dimension(0);

    std::vector<std::vector<trustIdType> > elems((size_t) nb_elem);
    trustIdType max_deg = 0;

    for (trustIdType e = 0; e < nb_elem; e++)
      {
        std::vector<trustIdType> nodes;
        std::unordered_set<trustIdType> seen;

        const int nb_faces_elem = (int) elem_faces.dimension(1);
        for (int jf = 0; jf < nb_faces_elem; jf++)
          {
            const trustIdType f = elem_faces(e, jf);
            if (f < 0)
              break;

            if (f >= faces.dimension(0))
              {
                cerr << "cgns_reader: invalid face index in elem_faces: " << f << endl;
                throw LataDBError(LataDBError::READ_ERROR);
              }

            const int nb_nodes_face = (int) faces.dimension(1);
            for (int js = 0; js < nb_nodes_face; js++)
              {
                const trustIdType s = faces(f, js);
                if (s < 0)
                  break;

                if (seen.insert(s).second)
                  nodes.push_back(s);
              }
          }

        elems[(size_t) e] = nodes;
        if ((trustIdType) nodes.size() > max_deg)
          max_deg = (trustIdType) nodes.size();
      }

    BigTIDTab elements;
    elements.resize(nb_elem, (int) max_deg);
    elements = -1;

    for (trustIdType e = 0; e < nb_elem; e++)
      for (int j = 0; j < (int) elems[(size_t) e].size(); j++)
        elements(e, j) = elems[(size_t) e][(size_t) j];

    return elements;
  }

  static bool read_zone_poly_connectivity(int fn, int ibase, int izone, const Nom &geom_name, ZoneConnectivityData &zc)
  {
    int isec_ngon = -1;
    int isec_nface = -1;

    const bool has_ngon = zone_has_section_type(fn, ibase, izone, NGON_n, &isec_ngon);
    const bool has_nface = zone_has_section_type(fn, ibase, izone, NFACE_n, &isec_nface);

    if (!has_ngon)
      return false;

    std::vector<cgsize_t> ngon_conn, ngon_offsets;
    trustIdType nb_faces = 0;
    read_cgns_poly_raw(fn, ibase, izone, isec_ngon, ngon_conn, ngon_offsets, nb_faces, NGON_n);

    BigTIDTab faces = build_padded_table_from_cgns_poly(ngon_offsets, ngon_conn, false);

    if (!has_nface)
      {
        zc.lata_elem_type = "POLYGONE";
        zc.elements = faces;
        return true;
      }

    std::vector<cgsize_t> nface_conn, nface_offsets;
    trustIdType nb_elem = 0;
    read_cgns_poly_raw(fn, ibase, izone, isec_nface, nface_conn, nface_offsets, nb_elem, NFACE_n);

    zc.lata_elem_type = "POLYEDRE";
    zc.faces = faces;
    zc.elem_faces = build_padded_table_from_cgns_poly(nface_offsets, nface_conn, true);
    zc.elements = build_polyhedron_elements_from_elem_faces_and_faces(zc.elem_faces, zc.faces);

    return true;
  }

  static trustIdType read_zone_fixed_elements(int fn, int ibase, int izone, int cell_dim, const Nom &geom_name,
                                              trustIdType nb_nodes, ZoneConnectivityData &zc)
  {
    int isec = -1;
    if (!try_choose_main_section(fn, ibase, izone, cell_dim, geom_name, isec))
      {
        zc.lata_elem_type = "";
        zc.elements.resize(0, 0);
        return 0;
      }

    char secname[33];
    ElementType_t elem_type;
    cgsize_t start = 0, end = 0;
    int nbndry = 0;
    int parent_flag = 0;

    cgns_check(cg_section_read(fn, ibase, izone, isec, secname, &elem_type, &start, &end, &nbndry, &parent_flag), "cg_section_read(main)");

    if (elem_type == NGON_n || elem_type == NFACE_n)
      return 0; // traité ailleurs

    int nb_comp = 0;
    const char *lata_elem_type = map_cgns_elemtype_to_lata(elem_type, nb_comp);
    if (!lata_elem_type)
      {
        cerr << "cgns_reader: unsupported CGNS element type in section " << secname << endl;
        throw LataDBError(LataDBError::READ_ERROR);
      }

    zc.lata_elem_type = lata_elem_type;

    const trustIdType nb_elem = (trustIdType) (end - start + 1);
    std::vector<cgsize_t> connectivity((size_t) nb_elem * (size_t) nb_comp);

    cgns_check(cg_elements_read(fn, ibase, izone, isec, connectivity.data(), nullptr), "cg_elements_read");

    zc.elements.resize(nb_elem, nb_comp);

    for (trustIdType i = 0; i < nb_elem; i++)
      {
        const cgsize_t *c = &connectivity[(size_t) i * nb_comp];

        if (elem_type == QUAD_4)
          {
            zc.elements(i, 0) = c[0] - 1;
            zc.elements(i, 1) = c[1] - 1;
            zc.elements(i, 2) = c[3] - 1;
            zc.elements(i, 3) = c[2] - 1;
          }
        else if (elem_type == HEXA_8)
          {
            zc.elements(i, 0) = c[0] - 1;
            zc.elements(i, 1) = c[1] - 1;
            zc.elements(i, 2) = c[3] - 1;
            zc.elements(i, 3) = c[2] - 1;
            zc.elements(i, 4) = c[4] - 1;
            zc.elements(i, 5) = c[5] - 1;
            zc.elements(i, 6) = c[7] - 1;
            zc.elements(i, 7) = c[6] - 1;
          }
        else
          {
            for (int j = 0; j < nb_comp; j++)
              zc.elements(i, j) = c[j] - 1;
          }

        for (int j = 0; j < nb_comp; j++)
          {
            if (zc.elements(i, j) < 0 || zc.elements(i, j) >= nb_nodes)
              {
                cerr << "cgns_reader: bad node index in connectivity for zone " << geom_name << " elem(" << i << "," << j << ")=" << zc.elements(i, j) << endl;
                throw LataDBError(LataDBError::READ_ERROR);
              }
          }
      }

    return nb_elem;
  }

  static trustIdType read_zone_connectivity(int fn, int ibase, int izone, int cell_dim, const Nom &geom_name, trustIdType nb_nodes, ZoneConnectivityData &zc)
  {
    zc.lata_elem_type = "";
    zc.elements.resize(0, 0);
    zc.faces.resize(0, 0);
    zc.elem_faces.resize(0, 0);

    if (read_zone_poly_connectivity(fn, ibase, izone, geom_name, zc))
      return zc.elements.dimension(0);

    return read_zone_fixed_elements(fn, ibase, izone, cell_dim, geom_name, nb_nodes, zc);
  }

} // namespace

#endif /* HAS_CGNS */

void cgns_reader(const char *cgnsfilename, const char *data_filename, LataDB &lata_db)
{
#ifndef HAS_CGNS
  cerr << "cgns_reader: code compiled without CGNS support" << endl;
  throw LataDBError(LataDBError::READ_ERROR);
#else
  const Nom filename_in_master_file = build_filename_in_master_file(data_filename);

  Journal() << "TRUST cgns_reader: reading file " << cgnsfilename << endl;

  int fn = -1;
  cgns_check(cg_open(cgnsfilename, CG_MODE_READ, &fn), "cg_open");

  lata_db.header_ = "CGNS";
  lata_db.case_ = cgnsfilename;
  lata_db.software_id_ = "Trio_U verbosity=0";

  int nbases = 0;
  cgns_check(cg_nbases(fn, &nbases), "cg_nbases");

  Journal(2) << "cgns_reader: nb_bases=" << nbases << endl;

  // Dump complet de la structure
  for (int ibase = 1; ibase <= nbases; ibase++)
    {
      char basename[33];
      int cell_dim = 0, phys_dim = 0;
      cgns_check(cg_base_read(fn, ibase, basename, &cell_dim, &phys_dim), "cg_base_read(debug)");

      Journal(2) << "cgns_reader: base[" << ibase << "]" << " name=" << basename << " cell_dim=" << cell_dim << " phys_dim=" << phys_dim << endl;

      int nzones = 0;
      cgns_check(cg_nzones(fn, ibase, &nzones), "cg_nzones(debug)");

      Journal(2) << "cgns_reader: base[" << ibase << "] nb_zones=" << nzones << endl;

      for (int izone = 1; izone <= nzones; izone++)
        {
          char zonename[33];
          cgsize_t size[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
          cgns_check(cg_zone_read(fn, ibase, izone, zonename, size), "cg_zone_read(debug)");

          ZoneType_t ztype;
          cgns_check(cg_zone_type(fn, ibase, izone, &ztype), "cg_zone_type(debug)");

          Journal(2) << "cgns_reader: zone[" << izone << "]" << " name=" << zonename << " type=" << zone_type_to_string(ztype) << " sizes=(" << size[0] << "," << size[1] << "," << size[2] << ","
              << size[3] << "," << size[4] << "," << size[5] << "," << size[6] << "," << size[7] << "," << size[8] << ")" << endl;

          dump_zone_sections(fn, ibase, izone, zonename);
          dump_zone_solutions(fn, ibase, izone, zonename);
        }
    }

  Journal(2) << "cgns_reader: timestep[0] = global definitions" << endl;
  lata_db.add_timestep(-1.); // timestep 0 = defs globales

  std::vector<double> time_values;
  if (nbases >= 1)
    time_values = read_cgns_time_values(fn, 1);

  if (!time_values.empty())
    {
      Journal(2) << "cgns_reader: nb_physical_timesteps=" << (int)time_values.size() << " (read from CGNS TimeValues)" << endl;
      for (int it = 0; it < (int)time_values.size(); it++)
        {
          lata_db.add_timestep(time_values[(size_t)it]);
          Journal(2) << "cgns_reader: timestep[" << (it + 1) << "] = " << time_values[(size_t)it] << " (from CGNS TimeValues)" << endl;
        }
    }
  else
    {
      int nb_phys_steps = get_max_number_of_solutions(fn);
      if (nb_phys_steps <= 0)
        nb_phys_steps = 1;

      Journal(2) << "cgns_reader: nb_physical_timesteps=" << nb_phys_steps << " (fallback from FlowSolution count)" << endl;
      for (int it = 0; it < nb_phys_steps; it++)
        {
          const double time_value = (double) it;
          lata_db.add_timestep(time_value);
          Journal(2) << "cgns_reader: timestep[" << (it + 1) << "] = " << time_value << " (synthetic fallback)" << endl;
        }
    }

  Journal(2) << "cgns_reader: lata_db.nb_timesteps()=" << lata_db.nb_timesteps() << endl;

  const int tstep_geom = 0;
  Size_t file_offset = 0;

  // Remplissage LataDB
  for (int ibase = 1; ibase <= nbases; ibase++)
    {
      char basename_c[33];
      int cell_dim = 0, phys_dim = 0;
      cgns_check(cg_base_read(fn, ibase, basename_c, &cell_dim, &phys_dim), "cg_base_read");

      const std::string basename_str(basename_c);
      const std::string basename_str_mod = strip_parallel_suffix_if_any(basename_str);
      const Nom geom_name(basename_str_mod.c_str());

      std::vector<ZonePartInfo> parts = collect_parallel_zone_parts(fn, ibase, basename_str, cell_dim);

      int nzones = 0;
      cgns_check(cg_nzones(fn, ibase, &nzones), "cg_nzones");

      // XXX cas comm_group => nb_bases/zones pas unique dans les fichiers ....
      // Donc on merge des qu'on a au moins une zone parallele reconnue.
      // XXX march eaussi pour les bases qui ne contiennent qu'une seule zone generique Zone_xxxx !!!
      const bool merge_parallel_over_zone = !parts.empty();

      if (!merge_parallel_over_zone)
        {
          for (int izone = 1; izone <= nzones; izone++)
            {
              char zonename[33];
              cgsize_t size[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
              cgns_check(cg_zone_read(fn, ibase, izone, zonename, size), "cg_zone_read");

              ZoneType_t ztype;
              cgns_check(cg_zone_type(fn, ibase, izone, &ztype), "cg_zone_type");

              if (ztype != Unstructured)
                {
                  cerr << "cgns_reader: skipping structured zone " << zonename << " for now" << endl;
                  continue;
                }

              const trustIdType nb_nodes = (trustIdType) size[0];
              const trustIdType nb_cells = (trustIdType) size[1];

              Journal(2) << "cgns_reader: filling zone=" << zonename
                         << " nb_nodes=" << nb_nodes
                         << " nb_cells=" << nb_cells
                         << " phys_dim=" << phys_dim << endl;

              const std::string zone_name_str(zonename);
              const std::string geom_name_str = build_geometry_name_from_base_and_zone(basename_str, zone_name_str);

              LataDBGeometry geom;
              geom.name_ = geom_name_str.c_str();
              geom.timestep_ = tstep_geom;
              lata_db.add_geometry(geom);

              BigFloatTab nodes;
              read_zone_coordinates(fn, ibase, izone, phys_dim, nb_nodes, nodes);
              add_sommets_field(filename_in_master_file, geom.name_, tstep_geom, file_offset,
                                nb_nodes, phys_dim, nodes, data_filename, lata_db);

              ZoneConnectivityData zc;
              const trustIdType nb_elem = read_zone_connectivity(fn, ibase, izone, cell_dim, geom.name_,
                                                                 nb_nodes, zc);

              if (nb_elem <= 0 || zc.lata_elem_type == "" || zc.elements.dimension(1) <= 0)
                {
                  Journal() << "cgns_reader: zone " << geom.name_
                            << " has no supported element section, skipping geometry-dependent fields" << endl;
                  continue;
                }

              add_elements_field(filename_in_master_file, geom.name_, tstep_geom, file_offset,
                                 zc.lata_elem_type, nb_elem, (int)zc.elements.dimension(1), zc.elements,
                                 data_filename, lata_db);

              if (zc.faces.dimension(0) > 0 && zc.faces.dimension(1) > 0)
                add_faces_field(filename_in_master_file, geom.name_, tstep_geom, file_offset,
                                zc.faces.dimension(0), (int)zc.faces.dimension(1), zc.faces,
                                data_filename, lata_db);

              if (zc.elem_faces.dimension(0) > 0 && zc.elem_faces.dimension(1) > 0)
                add_elem_faces_field(filename_in_master_file, geom.name_, tstep_geom, file_offset,
                                     zc.elem_faces.dimension(0), (int)zc.elem_faces.dimension(1), zc.elem_faces,
                                     data_filename, lata_db);

              int nsols = 0;
              cgns_check(cg_nsols(fn, ibase, izone, &nsols), "cg_nsols");

              for (int isol = 1; isol <= nsols; isol++)
                {
                  const int tstep_field = isol;
                  read_zone_solution_fields_and_add(fn, ibase, izone, isol,
                                                    filename_in_master_file, geom.name_,
                                                    tstep_field, file_offset,
                                                    nb_nodes, nb_elem, data_filename, lata_db);
                }
            }

          continue;
        }

      Journal(2) << "cgns_reader: base " << basename_c
                << " detected as parallel-over-zone, merging " << (int)parts.size()
                << " zones into one geometry" << endl;

      for (size_t ip = 0; ip < parts.size(); ip++)
        Journal(2) << "cgns_reader:   part[" << ip << "] zone=" << parts[ip].zonename
                   << " rank=" << parts[ip].rank
                   << " nb_nodes=" << parts[ip].nb_nodes
                   << " nb_cells=" << parts[ip].nb_cells
                   << " nb_cells_effective=" << parts[ip].nb_cells_effective << endl;

      std::vector<trustIdType> node_offsets(parts.size(), 0);
      std::vector<trustIdType> elem_offsets(parts.size(), 0);
      std::vector<trustIdType> face_offsets(parts.size(), 0);

      trustIdType total_nb_nodes = 0;
      trustIdType total_nb_elem = 0;
      trustIdType total_nb_faces = 0;

      for (size_t ip = 0; ip < parts.size(); ip++)
        {
          node_offsets[ip] = total_nb_nodes;
          elem_offsets[ip] = total_nb_elem;
          face_offsets[ip] = total_nb_faces;

          total_nb_nodes += parts[ip].nb_nodes;
          total_nb_elem += parts[ip].nb_cells_effective;

          ZoneConnectivityData tmp_zc;
          const trustIdType tmp_nb_elem = read_zone_connectivity(fn, ibase, parts[ip].izone, cell_dim,
              parts[ip].zonename.c_str(), parts[ip].nb_nodes, tmp_zc);

          if (tmp_nb_elem > 0 && tmp_zc.faces.dimension(0) > 0)
            total_nb_faces += tmp_zc.faces.dimension(0);
        }

      if (total_nb_nodes == 0)
        {
          Journal() << "cgns_reader: base " << basename_c << " is empty after filtering zones ... skipping geometry" << endl;
          continue;
        }

      LataDBGeometry geom;
      geom.name_ = geom_name;
      geom.timestep_ = tstep_geom;
      lata_db.add_geometry(geom);

      BigFloatTab merged_nodes;
      merged_nodes.resize(total_nb_nodes, phys_dim);

      BigTIDTab merged_elems, merged_faces, merged_elem_faces;
      trustIdType elem_write_pos = 0, face_write_pos = 0;
      int merged_nb_comp = 0, merged_faces_nb_comp = 0, merged_elem_faces_nb_comp = 0;
      Nom merged_elem_type;
      bool found_first_part = false;

      /* 1) Pre-scan : trouver le type et les largeurs max */
      for (size_t ip = 0; ip < parts.size(); ip++)
        {
          const ZonePartInfo& p = parts[ip];
          if (!p.has_main_section || p.nb_nodes == 0)
            continue;

          ZoneConnectivityData local_zc;
          const trustIdType local_nb_elem = read_zone_connectivity(fn, ibase, p.izone, cell_dim,
                                                                   p.zonename.c_str(), p.nb_nodes, local_zc);
          if (local_nb_elem == 0)
            continue;

          if (!found_first_part)
            {
              merged_elem_type = local_zc.lata_elem_type;
              found_first_part = true;
            }
          else if (local_zc.lata_elem_type != merged_elem_type)
            {
              cerr << "cgns_reader: inconsistent element type while merging base " << basename_c << endl;
              throw LataDBError(LataDBError::READ_ERROR);
            }

          merged_nb_comp = std::max(merged_nb_comp, (int)local_zc.elements.dimension(1));
          merged_faces_nb_comp = std::max(merged_faces_nb_comp, (int)local_zc.faces.dimension(1));
          merged_elem_faces_nb_comp = std::max(merged_elem_faces_nb_comp, (int)local_zc.elem_faces.dimension(1));

         /* Journal() << "merge poly part " << p.zonename
                    << " type=" << local_zc.lata_elem_type
                    << " elem_w=" << local_zc.elements.dimension(1)
                    << " faces_w=" << local_zc.faces.dimension(1)
                    << " elem_faces_w=" << local_zc.elem_faces.dimension(1)
                    << endl; */
        }

      if (!found_first_part || merged_nb_comp <= 0)
        {
          Journal() << "cgns_reader: merged base " << basename_c
                    << " has no supported element section, skipping" << endl;
          continue;
        }


      /* 2) Allocation avec les largeurs max */
      merged_elems.resize(total_nb_elem, merged_nb_comp);
      merged_elems = -1;

      if (total_nb_faces > 0 && merged_faces_nb_comp > 0)
        {
          merged_faces.resize(total_nb_faces, merged_faces_nb_comp);
          merged_faces = -1;
        }

      if (total_nb_elem > 0 && merged_elem_faces_nb_comp > 0)
        {
          merged_elem_faces.resize(total_nb_elem, merged_elem_faces_nb_comp);
          merged_elem_faces = -1;
        }

      /* 3) Merge effectif */
      for (size_t ip = 0; ip < parts.size(); ip++)
        {
          const ZonePartInfo& p = parts[ip];

          if (!p.has_main_section || p.nb_nodes == 0)
            continue;

          BigFloatTab local_nodes;
          read_zone_coordinates(fn, ibase, p.izone, phys_dim, p.nb_nodes, local_nodes);

          for (trustIdType i = 0; i < p.nb_nodes; i++)
            for (int j = 0; j < phys_dim; j++)
              merged_nodes(node_offsets[ip] + i, j) = local_nodes(i, j);

          ZoneConnectivityData local_zc;
          const trustIdType local_nb_elem = read_zone_connectivity(fn, ibase, p.izone, cell_dim,
                                                                   p.zonename.c_str(), p.nb_nodes, local_zc);

          if (local_nb_elem == 0)
            continue;

          if (local_zc.lata_elem_type != merged_elem_type)
            {
              cerr << "cgns_reader: inconsistent element type while merging base " << basename_c << endl;
              throw LataDBError(LataDBError::READ_ERROR);
            }

          for (trustIdType i = 0; i < local_nb_elem; i++)
            for (int j = 0; j < (int)local_zc.elements.dimension(1); j++)
              {
                const trustIdType v = local_zc.elements(i, j);
                merged_elems(elem_write_pos + i, j) = (v >= 0) ? (v + node_offsets[ip]) : -1;
              }

          if (local_zc.faces.dimension(0) > 0)
            {
              const trustIdType local_nb_faces = local_zc.faces.dimension(0);

              for (trustIdType i = 0; i < local_nb_faces; i++)
                for (int j = 0; j < (int)local_zc.faces.dimension(1); j++)
                  {
                    const trustIdType v = local_zc.faces(i, j);
                    merged_faces(face_write_pos + i, j) = (v >= 0) ? (v + node_offsets[ip]) : -1;
                  }

              for (trustIdType i = 0; i < local_nb_elem; i++)
                for (int j = 0; j < (int)local_zc.elem_faces.dimension(1); j++)
                  {
                    const trustIdType f = local_zc.elem_faces(i, j);
                    merged_elem_faces(elem_write_pos + i, j) = (f >= 0) ? (f + face_offsets[ip]) : -1;
                  }

              face_write_pos += local_nb_faces;
            }

          elem_write_pos += local_nb_elem;
        }

      add_sommets_field(filename_in_master_file, geom.name_, tstep_geom, file_offset,
                        total_nb_nodes, phys_dim, merged_nodes, data_filename, lata_db);

      if (merged_nb_comp <= 0 || merged_elem_type == "" || total_nb_elem <= 0)
        {
          Journal() << "cgns_reader: merged base " << basename_c << " has no supported element section, skipping" << endl;
          continue;
        }

      add_elements_field(filename_in_master_file, geom.name_, tstep_geom, file_offset,
                         merged_elem_type, total_nb_elem, merged_nb_comp, merged_elems,
                         data_filename, lata_db);

      if (merged_faces.dimension(0) > 0 && merged_faces.dimension(1) > 0)
        add_faces_field(filename_in_master_file, geom.name_, tstep_geom, file_offset,
                        merged_faces.dimension(0), (int)merged_faces.dimension(1), merged_faces,
                        data_filename, lata_db);

      if (merged_elem_faces.dimension(0) > 0 && merged_elem_faces.dimension(1) > 0)
        add_elem_faces_field(filename_in_master_file, geom.name_, tstep_geom, file_offset,
                             merged_elem_faces.dimension(0), (int)merged_elem_faces.dimension(1), merged_elem_faces,
                             data_filename, lata_db);

      int iref = -1;
      for (size_t ip = 0; ip < parts.size(); ip++)
        if (parts[ip].has_main_section)
          {
            iref = (int)ip;
            break;
          }

      if (iref >= 0)
        {
          int nsols = 0;
          cgns_check(cg_nsols(fn, ibase, parts[iref].izone, &nsols), "cg_nsols(merged.first_zone)");

          for (int isol = 1; isol <= nsols; isol++)
            {
              const int tstep_field = isol;
              read_merged_zone_solution_fields_and_add(fn, ibase, parts,
                                                       filename_in_master_file, geom.name_,
                                                       tstep_field, file_offset,
                                                       node_offsets, elem_offsets,
                                                       total_nb_nodes, total_nb_elem,
                                                       data_filename, lata_db);
            }
        }
    }


  cgns_check(cg_close(fn), "cg_close");

#endif /* HAS_CGNS */
}

void cgns_to_lata(const char *cgns_name, const char *latafilename, bool ascii, bool fortran_blocs, bool fortran_ordering, bool fortran_indexing)
{
  std::string lata_name(latafilename);

  Motcle motcle_nom_fic(latafilename);

  if (!motcle_nom_fic.finit_par(".lata"))
    lata_name += ".lata";

  Journal() << "cgns_to_lata " << cgns_name << " -> " << lata_name << endl;
  LataDB lata_db;
  Nom dest_prefix, dest_name;
  LataOptions::extract_path_basename(lata_name.c_str(), dest_prefix, dest_name);

  // Nom du fichier .data a ecrire (sans le chemin)
  Nom datafile(dest_name);
  datafile += ".lata_single";
  lata_db.set_path_prefix(dest_prefix);

  // Nom complet du fichier cgns a lire
  LataDBDataType type;
  if (ascii)
    type.msb_ = LataDBDataType::ASCII;
  else
    type.msb_ = LataDBDataType::machine_msb_;

  type.type_ = LataDBDataType::INT64;
  type.array_index_ = fortran_indexing ? LataDBDataType::F_INDEXING : LataDBDataType::C_INDEXING;
  type.data_ordering_ = fortran_ordering ? LataDBDataType::F_ORDERING : LataDBDataType::C_ORDERING;
  type.fortran_bloc_markers_ = fortran_blocs ? LataDBDataType::BLOC_MARKERS_SINGLE_WRITE : LataDBDataType::NO_BLOC_MARKER;
  type.bloc_marker_type_ = LataDBDataType::INT64;
  type.file_offset_ = 0;
  lata_db.default_type_int_ = type;
  lata_db.default_float_type_ = LataDBDataType::REAL32;

  cgns_reader(cgns_name, datafile, lata_db);
  Journal() << "cgns_to_lata writing single_lata master file" << endl;
  lata_db.write_master_file(lata_name.c_str());
}

#endif /* WITH_CGNSLOADER */
