
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

#pragma GCC diagnostic push

#if __GNUC__ < 9
#pragma GCC diagnostic ignored "-Wsuggest-override"
#endif

#include <gtest/gtest.h>

#pragma GCC diagnostic pop

#include <verifie_pere.h>

#include <Objet_U_With_Params.h>
#include <Param.h>
#include <EChaineJDD.h>
#include <Interprete_bloc.h>



template<typename T>
bool contains(const std::map<std::string, T> map, const std::string& key) {
    return (map.find(key) != map.end());
}

class Named_Object : public Objet_U_With_Params
{
  Declare_instanciable_with_param(Named_Object);
public:



  // Implementing naming logic from Objet_U
  void nommer(const Nom& n) override { name= n; }
  const Nom& le_nom() const override { return name;}
  Nom name;


    std::string value = "";


};

Implemente_instanciable(Named_Object, "Named_Object", Objet_U_With_Params);


Sortie& Named_Object::printOn(Sortie& os) const
{
  return os;
}
void Named_Object::set_param(Param& param) const
{
  param.ajouter("value", &value, Param::OPTIONAL);
}

class Named_Object_Specialized : public Named_Object
{
  Declare_instanciable_with_param(Named_Object_Specialized);
public:
    std::string other_value = "";
};

Implemente_instanciable(Named_Object_Specialized, "Named_Object_Specialized", Named_Object);


Sortie& Named_Object_Specialized::printOn(Sortie& os) const
{
  return os;
}
void Named_Object_Specialized::set_param(Param& param) const
{
    Named_Object::set_param(param);
    param.ajouter("other_value", &other_value, Param::OPTIONAL);
}





class Test_Param : public Objet_U_With_Params
{
  Declare_instanciable_with_param(Test_Param);
public:
  void validate_params() const override
  {

  }


  std::vector<int> vec_ints;
  std::vector<double> vec_dbls;
  std::vector<std::string> vec_strs;

  std::vector<OWN_PTR(Nom)> vec_objs;

  std::map<std::string, int> map_ints;
  std::map<std::string, double> map_dbls;
  std::map<std::string, std::string> map_strs;

  std::map<std::string, OWN_PTR(Nom)> map_objs;

  std::vector<int> vec_fixed_size;

  // for testing naming of objects
  std::map<std::string, OWN_PTR(Named_Object)> map_named;


};

Implemente_instanciable(Test_Param, "Test_Param", Objet_U_With_Params);


Sortie& Test_Param::printOn(Sortie& os) const
{
  return os;
}
void Test_Param::set_param(Param& param) const
{
  param.ajouter("vec_fixed_size", &vec_fixed_size, Param::OPTIONAL,  4); // must read 4 components

  param.ajouter("vec_ints", &vec_ints, Param::OPTIONAL);
  param.ajouter("vec_dbls", &vec_dbls, Param::OPTIONAL);
  param.ajouter("vec_strs", &vec_strs, Param::OPTIONAL);
  param.ajouter("vec_objs", &vec_objs, Param::OPTIONAL);

  param.ajouter("map_ints", &map_ints, Param::OPTIONAL);
  param.ajouter("map_dbls", &map_dbls, Param::OPTIONAL);
  param.ajouter("map_strs", &map_strs, Param::OPTIONAL);
  param.ajouter("map_objs", &map_objs, Param::OPTIONAL);

  // for testing naming of objects
  param.ajouter("map_named", &map_named, Param::OPTIONAL);
}


static EChaineJDD lire_inline_jdd(const  std::string & jdd){

    return EChaineJDD(jdd.c_str());
}
static Interprete_bloc interprete_inline_jdd(const  std::string & jdd){

    EChaineJDD e = lire_inline_jdd(jdd);

    Interprete_bloc ib;

    ib.interpreter_bloc(e, Interprete_bloc::FIN , 0);

return ib;
}


TEST(Param, Vector_JSON_Syntax) {

    // the comma at the end of vec_strs must be kept for testing
    std::string jdd =  R"(

Test_Param Vector_JSON_Syntax
Lire Vector_JSON_Syntax
{
    vec_ints [ 1 , 2 ]
    vec_dbls [ 1.2 , -2.333 , 1e6 ]
    vec_strs [ azdafa , 2adfzaa , azefaz , ]
    vec_objs [ Nom val1 ,  Nom val2 ,  Motcle subtypes_allowed ]
}
End

)";

     Interprete_bloc ib = interprete_inline_jdd( jdd);

    ASSERT_TRUE(ib.objet_local_existant("Vector_JSON_Syntax"));

    auto& obj = ib.objet_local("Vector_JSON_Syntax");

    EXPECT_EQ(obj.le_type(), Nom("Test_Param"));

    Test_Param& test_obj = ref_cast(Test_Param, obj);


    EXPECT_TRUE(test_obj.map_ints.empty());
    EXPECT_TRUE(test_obj.map_dbls.empty());
    EXPECT_TRUE(test_obj.map_strs.empty());
    EXPECT_TRUE(test_obj.map_objs.empty());


    ASSERT_EQ(test_obj.vec_ints.size(), 2ul);
    ASSERT_EQ(test_obj.vec_dbls.size(), 3ul);
    ASSERT_EQ(test_obj.vec_strs.size(), 3ul);
    ASSERT_EQ(test_obj.vec_objs.size(), 3ul);

    EXPECT_EQ(test_obj.vec_ints[0], 1);
    EXPECT_EQ(test_obj.vec_ints[1], 2);

    EXPECT_EQ(test_obj.vec_dbls[0], 1.2);
    EXPECT_EQ(test_obj.vec_dbls[1], -2.333);
    EXPECT_EQ(test_obj.vec_dbls[2], 1e6);

    EXPECT_EQ(test_obj.vec_strs[0], "azdafa");
    EXPECT_EQ(test_obj.vec_strs[1], "2adfzaa");
    EXPECT_EQ(test_obj.vec_strs[2], "azefaz");

    EXPECT_EQ(std::string(test_obj.vec_objs[0]->le_type()), "Nom");
    EXPECT_EQ(std::string(test_obj.vec_objs[1]->le_type()), "Nom");
    EXPECT_EQ(std::string(test_obj.vec_objs[2]->le_type()), "Motcle");
    EXPECT_EQ(test_obj.vec_objs[0]->getString(), "val1");
    EXPECT_EQ(test_obj.vec_objs[1]->getString(), "val2");
    EXPECT_EQ(test_obj.vec_objs[2]->getString(), "SUBTYPES_ALLOWED"); // Motcle is case insensitive

}

TEST(Param, Vector_TRUST_Syntax) {
    std::string jdd =  R"(

Test_Param Vector_TRUST_Syntax
Lire Vector_TRUST_Syntax
{
    vec_ints 2 1 2
    vec_dbls 3 1.2  -2.333  1e6 
    vec_strs 3 azdafa  2adfzaa  azefaz
    vec_objs 3 Nom val1   Nom val2   Motcle subtypes_allowed

}
End

)";

     Interprete_bloc ib = interprete_inline_jdd( jdd);

    ASSERT_TRUE(ib.objet_local_existant("Vector_TRUST_Syntax"));

    auto& obj = ib.objet_local("Vector_TRUST_Syntax");

    EXPECT_EQ(obj.le_type(), Nom("Test_Param"));

    Test_Param& test_obj = ref_cast(Test_Param, obj);


    EXPECT_TRUE(test_obj.map_ints.empty());
    EXPECT_TRUE(test_obj.map_dbls.empty());
    EXPECT_TRUE(test_obj.map_strs.empty());
    EXPECT_TRUE(test_obj.map_objs.empty());


    ASSERT_EQ(test_obj.vec_ints.size(), 2ul);
    ASSERT_EQ(test_obj.vec_dbls.size(), 3ul);
    ASSERT_EQ(test_obj.vec_strs.size(), 3ul);
    ASSERT_EQ(test_obj.vec_objs.size(), 3ul);

    EXPECT_EQ(test_obj.vec_ints[0], 1);
    EXPECT_EQ(test_obj.vec_ints[1], 2);

    EXPECT_EQ(test_obj.vec_dbls[0], 1.2);
    EXPECT_EQ(test_obj.vec_dbls[1], -2.333);
    EXPECT_EQ(test_obj.vec_dbls[2], 1e6);

    EXPECT_EQ(test_obj.vec_strs[0], "azdafa");
    EXPECT_EQ(test_obj.vec_strs[1], "2adfzaa");
    EXPECT_EQ(test_obj.vec_strs[2], "azefaz");

    EXPECT_EQ(std::string(test_obj.vec_objs[0]->le_type()), "Nom");
    EXPECT_EQ(std::string(test_obj.vec_objs[1]->le_type()), "Nom");
    EXPECT_EQ(std::string(test_obj.vec_objs[2]->le_type()), "Motcle");
    EXPECT_EQ(test_obj.vec_objs[0]->getString(), "val1");
    EXPECT_EQ(test_obj.vec_objs[1]->getString(), "val2");
    EXPECT_EQ(test_obj.vec_objs[2]->getString(), "SUBTYPES_ALLOWED"); // Motcle is case insensitive

}



TEST(Param, Vector_Fixed_JSON_Syntax) {

    // the comma at the end of vec_strs must be kept for testing
    std::string jdd =  R"(

Test_Param Vector_Fixed_JSON_Syntax
Lire Vector_Fixed_JSON_Syntax
{
    vec_fixed_size [ 1 , 2 , 3 , 4 ]
}
End

)";

     Interprete_bloc ib = interprete_inline_jdd( jdd);

    ASSERT_TRUE(ib.objet_local_existant("Vector_Fixed_JSON_Syntax"));

    auto& obj = ib.objet_local("Vector_Fixed_JSON_Syntax");

    ASSERT_EQ(obj.le_type(), Nom("Test_Param"));

    Test_Param& test_obj = ref_cast(Test_Param, obj);



    ASSERT_EQ(test_obj.vec_fixed_size.size(), 4ul);
    EXPECT_EQ(test_obj.vec_fixed_size[0], 1);
    EXPECT_EQ(test_obj.vec_fixed_size[1], 2);
    EXPECT_EQ(test_obj.vec_fixed_size[2], 3);
    EXPECT_EQ(test_obj.vec_fixed_size[3], 4);

}

TEST(Param, Vector_Fixed_TRUST_Syntax) {

    // the comma at the end of vec_strs must be kept for testing
    std::string jdd =  R"(

Test_Param Vector_Fixed_TRUST_Syntax
Lire Vector_Fixed_TRUST_Syntax
{
    vec_fixed_size 4 1  2  3  4 
}
End

)";

     Interprete_bloc ib = interprete_inline_jdd( jdd);

    ASSERT_TRUE(ib.objet_local_existant("Vector_Fixed_TRUST_Syntax"));

    auto& obj = ib.objet_local("Vector_Fixed_TRUST_Syntax");

    ASSERT_EQ(obj.le_type(), Nom("Test_Param"));

    Test_Param& test_obj = ref_cast(Test_Param, obj);


    EXPECT_EQ(test_obj.vec_fixed_size.size(), 4ul);
    EXPECT_EQ(test_obj.vec_fixed_size[0], 1);
    EXPECT_EQ(test_obj.vec_fixed_size[1], 2);
    EXPECT_EQ(test_obj.vec_fixed_size[2], 3);
    EXPECT_EQ(test_obj.vec_fixed_size[3], 4);

}

TEST(Param, Map_JSON_Syntax) {

    // the comma at the end of map_strs must be kept for testing
    std::string jdd =  R"(

Test_Param Map_JSON_Syntax
Lire Map_JSON_Syntax
{
    map_ints { a : 1 , b : 2 }
    map_dbls { a : 1.2 , b : -2.333 , c : 1e6 }
    map_strs { a : azdafa , b : 2adfzaa ,
        c : azefaz , note_that : trailing_commas_are_allowed! ,  }
    map_objs { key1 : Nom val1 , key2 : Nom val2 , key3 : Motcle subtypes_allowed }
}
End

)";

     Interprete_bloc ib = interprete_inline_jdd( jdd);

    ASSERT_TRUE(ib.objet_local_existant("Map_JSON_Syntax"));

    auto& obj = ib.objet_local("Map_JSON_Syntax");

    ASSERT_EQ(obj.le_type(), Nom("Test_Param"));

    Test_Param& test_obj = ref_cast(Test_Param, obj);

    EXPECT_TRUE(test_obj.vec_ints.empty());
    EXPECT_TRUE(test_obj.vec_dbls.empty());
    EXPECT_TRUE(test_obj.vec_strs.empty());
    EXPECT_TRUE(test_obj.vec_objs.empty());

    ASSERT_TRUE(contains(test_obj.map_ints, "a"));
    ASSERT_TRUE(contains(test_obj.map_ints, "b"));
    EXPECT_EQ(test_obj.map_ints.at("a"), 1);
    EXPECT_EQ(test_obj.map_ints.at("b"), 2);

    ASSERT_TRUE(contains(test_obj.map_dbls, "a"));
    ASSERT_TRUE(contains(test_obj.map_dbls, "b"));
    ASSERT_TRUE(contains(test_obj.map_dbls, "c"));
    EXPECT_EQ(test_obj.map_dbls.at("a"), 1.2);
    EXPECT_EQ(test_obj.map_dbls.at("b"), -2.333);
    EXPECT_EQ(test_obj.map_dbls.at("c"), 1e6);

    ASSERT_TRUE(contains(test_obj.map_strs, "a"));
    ASSERT_TRUE(contains(test_obj.map_strs, "b"));
    ASSERT_TRUE(contains(test_obj.map_strs, "c"));
    EXPECT_EQ(test_obj.map_strs.at("a"), "azdafa");
    EXPECT_EQ(test_obj.map_strs.at("b"), "2adfzaa");
    EXPECT_EQ(test_obj.map_strs.at("c"), "azefaz");


    ASSERT_TRUE(contains(test_obj.map_objs, "key1"));
    ASSERT_TRUE(contains(test_obj.map_objs, "key2"));
    ASSERT_TRUE(contains(test_obj.map_objs, "key3"));
    EXPECT_EQ(std::string(test_obj.map_objs.at("key1")->le_type()), "Nom");
    EXPECT_EQ(std::string(test_obj.map_objs.at("key2")->le_type()), "Nom");
    EXPECT_EQ(std::string(test_obj.map_objs.at("key3")->le_type()), "Motcle");
    EXPECT_EQ(test_obj.map_objs.at("key1")->getString(), "val1");
    EXPECT_EQ(test_obj.map_objs.at("key2")->getString(), "val2");
    EXPECT_EQ(test_obj.map_objs.at("key3")->getString(), "SUBTYPES_ALLOWED"); // Motcle is case insensitive

}

TEST(Param, Map_TRUST_Syntax) {
     std::string  jdd =  R"(

Test_Param Map_TRUST_Syntax
Lire Map_TRUST_Syntax
{
    map_ints { a  1  b  2 }
    map_dbls { a 1.2  b  -2.333  b2  1   c  1e6 }
    map_strs { a  azdafa  b 2adfzaa
        c  azefaz  }

    map_objs { key1  Nom val1 key2  Nom val2 key3  Motcle subtypes_allowed }
}
End

)";

    Interprete_bloc ib = interprete_inline_jdd( jdd);

    ASSERT_TRUE(ib.objet_local_existant("Map_TRUST_Syntax"));

    auto& obj = ib.objet_local("Map_TRUST_Syntax");

    ASSERT_EQ(obj.le_type(), Nom("Test_Param"));

    Test_Param& test_obj = ref_cast(Test_Param, obj);

    EXPECT_TRUE(test_obj.vec_ints.empty());
    EXPECT_TRUE(test_obj.vec_dbls.empty());
    EXPECT_TRUE(test_obj.vec_strs.empty());
    EXPECT_TRUE(test_obj.vec_objs.empty());

    EXPECT_TRUE(contains(test_obj.map_ints, "a"));
    EXPECT_TRUE(contains(test_obj.map_ints, "b"));
    EXPECT_EQ(test_obj.map_ints.at("a"), 1);
    EXPECT_EQ(test_obj.map_ints.at("b"), 2);

    EXPECT_TRUE(contains(test_obj.map_dbls, "a"));
    EXPECT_TRUE(contains(test_obj.map_dbls, "b"));
    EXPECT_TRUE(contains(test_obj.map_dbls, "c"));
    EXPECT_EQ(test_obj.map_dbls.at("a"), 1.2);
    EXPECT_EQ(test_obj.map_dbls.at("b"), -2.333);
    EXPECT_EQ(test_obj.map_dbls.at("c"), 1e6);

    EXPECT_TRUE(contains(test_obj.map_strs, "a"));
    EXPECT_TRUE(contains(test_obj.map_strs, "b"));
    EXPECT_TRUE(contains(test_obj.map_strs, "c"));
    EXPECT_EQ(test_obj.map_strs.at("a"), "azdafa");
    EXPECT_EQ(test_obj.map_strs.at("b"), "2adfzaa");
    EXPECT_EQ(test_obj.map_strs.at("c"), "azefaz");


    EXPECT_TRUE(contains(test_obj.map_objs, "key1"));
    EXPECT_TRUE(contains(test_obj.map_objs, "key2"));
    EXPECT_TRUE(contains(test_obj.map_objs, "key3"));
    EXPECT_EQ(std::string(test_obj.map_objs.at("key1")->le_type()), "Nom");
    EXPECT_EQ(std::string(test_obj.map_objs.at("key2")->le_type()), "Nom");
    EXPECT_EQ(std::string(test_obj.map_objs.at("key3")->le_type()), "Motcle");
    EXPECT_EQ(test_obj.map_objs.at("key1")->getString(), "val1");
    EXPECT_EQ(test_obj.map_objs.at("key2")->getString(), "val2");
    EXPECT_EQ(test_obj.map_objs.at("key3")->getString(), "SUBTYPES_ALLOWED"); // Motcle is case insensitive

}

TEST(Param, Map_Object_Naming) {

    // the comma at the end of map_strs must be kept for testing
    std::string jdd =  R"(

Test_Param Map_Object_Naming
Lire Map_Object_Naming
{
    map_named {
        k1 : Named_Object { value foo } ,
        k2 : Named_Object { value bar } ,
        k3 : Named_Object_Specialized { value bar other_value foobar } ,
    }
}
End

)";

     Interprete_bloc ib = interprete_inline_jdd( jdd);

    ASSERT_TRUE(ib.objet_local_existant("Map_Object_Naming"));

    auto& obj = ib.objet_local("Map_Object_Naming");

    ASSERT_EQ(obj.le_type(), Nom("Test_Param"));

    Test_Param& test_obj = ref_cast(Test_Param, obj);


    ASSERT_TRUE(contains(test_obj.map_named,"k1"));
    EXPECT_EQ(std::string(test_obj.map_named.at("k1")->le_type()), "Named_Object");
    EXPECT_EQ(std::string(test_obj.map_named.at("k1")->le_nom()), "k1");
    EXPECT_EQ(test_obj.map_named.at("k1")->value, "foo");

    ASSERT_TRUE(contains(test_obj.map_named, "k2"));
    EXPECT_EQ(std::string(test_obj.map_named.at("k2")->le_type()), "Named_Object");
    EXPECT_EQ(std::string(test_obj.map_named.at("k2")->le_nom()), "k2");
    EXPECT_EQ(test_obj.map_named.at("k2")->value, "bar");

    ASSERT_TRUE(contains(test_obj.map_named, "k3"));
    ASSERT_EQ(std::string(test_obj.map_named.at("k3")->le_type()), "Named_Object_Specialized"); // assert because refcast later
    EXPECT_EQ(std::string(test_obj.map_named.at("k3")->le_nom()), "k3");
    EXPECT_EQ(test_obj.map_named.at("k3")->value, "bar");
    EXPECT_EQ(ref_cast(Named_Object_Specialized, test_obj.map_named.at("k3").valeur()).other_value, "foobar");

}
// ===== Vector edge cases =====

TEST(Param, Vector_Empty_JSON_Syntax)
{
    std::string jdd = R"(
Test_Param t
Lire t { vec_ints [ ] }
End
)";
    Interprete_bloc ib = interprete_inline_jdd(jdd);
    Test_Param& obj = ref_cast(Test_Param, ib.objet_local("t"));
    EXPECT_TRUE(obj.vec_ints.empty());
}

TEST(Param, Vector_Empty_TRUST_Syntax)
{
    std::string jdd = R"(
Test_Param t
Lire t { vec_ints 0 }
End
)";
    Interprete_bloc ib = interprete_inline_jdd(jdd);
    Test_Param& obj = ref_cast(Test_Param, ib.objet_local("t"));
    EXPECT_TRUE(obj.vec_ints.empty());
}

TEST(Param, Vector_Single_Element_JSON_Syntax)
{
    std::string jdd = R"(
Test_Param t
Lire t { vec_ints [ 42 ] }
End
)";
    Interprete_bloc ib = interprete_inline_jdd(jdd);
    Test_Param& obj = ref_cast(Test_Param, ib.objet_local("t"));
    ASSERT_EQ(obj.vec_ints.size(), 1ul);
    EXPECT_EQ(obj.vec_ints[0], 42);
}

TEST(Param, Vector_Single_Element_TRUST_Syntax)
{
    std::string jdd = R"(
Test_Param t
Lire t { vec_ints 1 42 }
End
)";
    Interprete_bloc ib = interprete_inline_jdd(jdd);
    Test_Param& obj = ref_cast(Test_Param, ib.objet_local("t"));
    ASSERT_EQ(obj.vec_ints.size(), 1ul);
    EXPECT_EQ(obj.vec_ints[0], 42);
}


TEST(Param, Vector_Int_Boundary_Values)
{
    std::string jdd = R"(
Test_Param t
Lire t { vec_ints [ 2147483647 , -2147483648 ] }
End
)";
    Interprete_bloc ib = interprete_inline_jdd(jdd);
    Test_Param& obj = ref_cast(Test_Param, ib.objet_local("t"));
    ASSERT_EQ(obj.vec_ints.size(), 2ul);
    EXPECT_EQ(obj.vec_ints[0], INT_MAX);
    EXPECT_EQ(obj.vec_ints[1], INT_MIN);
}

TEST(Param, Vector_Double_Boundary_Values)
{
    std::string jdd = R"(
Test_Param t
Lire t { vec_dbls [ 1e300 , 1e-300 , -1e300 ] }
End
)";
    Interprete_bloc ib = interprete_inline_jdd(jdd);
    Test_Param& obj = ref_cast(Test_Param, ib.objet_local("t"));
    ASSERT_EQ(obj.vec_dbls.size(), 3ul);
    EXPECT_DOUBLE_EQ(obj.vec_dbls[0], 1e300);
    EXPECT_DOUBLE_EQ(obj.vec_dbls[1], 1e-300);
    EXPECT_DOUBLE_EQ(obj.vec_dbls[2], -1e300);
}

TEST(Param, Vector_Trailing_Comma_Single_Element)
{
    // Trailing comma with a single element - exercises the [ val , ] path
    std::string jdd = R"(
Test_Param t
Lire t { vec_ints [ 1 , ] }
End
)";
    Interprete_bloc ib = interprete_inline_jdd(jdd);
    Test_Param& obj = ref_cast(Test_Param, ib.objet_local("t"));
    ASSERT_EQ(obj.vec_ints.size(), 1ul);
    EXPECT_EQ(obj.vec_ints[0], 1);
}

// ===== Map edge cases =====

TEST(Param, Map_Empty_JSON_Syntax)
{
    std::string jdd = R"(
Test_Param t
Lire t { map_ints { } }
End
)";
    Interprete_bloc ib = interprete_inline_jdd(jdd);
    Test_Param& obj = ref_cast(Test_Param, ib.objet_local("t"));
    EXPECT_TRUE(obj.map_ints.empty());
}

TEST(Param, Map_Empty_TRUST_Syntax)
{
    std::string jdd = R"(
Test_Param t
Lire t { map_ints { } }
End
)";
    Interprete_bloc ib = interprete_inline_jdd(jdd);
    Test_Param& obj = ref_cast(Test_Param, ib.objet_local("t"));
    EXPECT_TRUE(obj.map_ints.empty());
}

TEST(Param, Map_Single_Entry_JSON_Syntax)
{
    std::string jdd = R"(
Test_Param t
Lire t { map_ints { a : 42 } }
End
)";
    Interprete_bloc ib = interprete_inline_jdd(jdd);
    Test_Param& obj = ref_cast(Test_Param, ib.objet_local("t"));
    ASSERT_EQ(obj.map_ints.size(), 1ul);
    EXPECT_EQ(obj.map_ints.at("a"), 42);
}

TEST(Param, Map_Single_Entry_TRUST_Syntax)
{
    std::string jdd = R"(
Test_Param t
Lire t { map_ints { a 42 } }
End
)";
    Interprete_bloc ib = interprete_inline_jdd(jdd);
    Test_Param& obj = ref_cast(Test_Param, ib.objet_local("t"));
    ASSERT_EQ(obj.map_ints.size(), 1ul);
    EXPECT_EQ(obj.map_ints.at("a"), 42);
}

TEST(Param, Map_Trailing_Comma_Single_Entry)
{
    std::string jdd = R"(
Test_Param t
Lire t { map_ints { a : 42 , } }
End
)";
    Interprete_bloc ib = interprete_inline_jdd(jdd);
    Test_Param& obj = ref_cast(Test_Param, ib.objet_local("t"));
    ASSERT_EQ(obj.map_ints.size(), 1ul);
    EXPECT_EQ(obj.map_ints.at("a"), 42);
}

TEST(Param, Map_Double_Boundary_Values)
{
    std::string jdd = R"(
Test_Param t
Lire t { map_dbls { big : 1e300 , small : 1e-300 , neg : -1e300 } }
End
)";
    Interprete_bloc ib = interprete_inline_jdd(jdd);
    Test_Param& obj = ref_cast(Test_Param, ib.objet_local("t"));
    ASSERT_EQ(obj.map_dbls.size(), 3ul);
    EXPECT_DOUBLE_EQ(obj.map_dbls.at("big"),   1e300);
    EXPECT_DOUBLE_EQ(obj.map_dbls.at("small"), 1e-300);
    EXPECT_DOUBLE_EQ(obj.map_dbls.at("neg"),  -1e300);
}







TEST(Param_Death, Vector_Fixed_Too_Few_Elements)
{
    std::string jdd = R"(
Test_Param t
Lire t { vec_fixed_size [ 1 , 2 , 3 ] }
End
)";
    EXPECT_DEATH(interprete_inline_jdd(jdd), "");
}

TEST(Param_Death, Vector_Fixed_Too_Many_Elements)
{
    std::string jdd = R"(
Test_Param t
Lire t { vec_fixed_size [ 1 , 2 , 3 , 4 , 5 ] }
End
)";
    EXPECT_DEATH(interprete_inline_jdd(jdd), "");
}

TEST(Param_Death, Vector_Fixed_TRUST_Declared_Size_Mismatch)
{
    // declared size in stream (3) != expected_size (4)
    std::string jdd = R"(
Test_Param t
Lire t { vec_fixed_size 3 1 2 3 }
End
)";
    EXPECT_DEATH(interprete_inline_jdd(jdd), "");
}

TEST(Param_Death, Vector_Fixed_TRUST_Too_Many_Declared)
{
    // declared size in stream (5) != expected_size (4)
    std::string jdd = R"(
Test_Param t
Lire t { vec_fixed_size 5 1 2 3 4 5 }
End
)";
    EXPECT_DEATH(interprete_inline_jdd(jdd), "");
}

TEST(Param_Death, Vector_Invalid_Prefix_Token)
{
    std::string jdd = R"(
Test_Param t
Lire t { vec_ints notanumber }
End
)";
    EXPECT_DEATH(interprete_inline_jdd(jdd), "");
}

TEST(Param_Death, Vector_Invalid_Separator)
{
    std::string jdd = R"(
Test_Param t
Lire t { vec_ints [ 1 ; 2 ] }
End
)";
    EXPECT_DEATH(interprete_inline_jdd(jdd), "");
}

TEST(Param_Death, Map_Duplicate_Key_JSON_Syntax)
{
    std::string jdd = R"(
Test_Param t
Lire t { map_ints { a : 1 , a : 2 } }
End
)";
    EXPECT_DEATH(interprete_inline_jdd(jdd), "");
}

TEST(Param_Death, Map_Duplicate_Key_TRUST_Syntax)
{
    std::string jdd = R"(
Test_Param t
Lire t { map_ints { a 1 a 2 } }
End
)";
    EXPECT_DEATH(interprete_inline_jdd(jdd), "");
}

TEST(Param_Death, Map_Missing_Colon_JSON_Syntax)
{
    std::string jdd = R"(
Test_Param t
Lire t { map_ints { a 1 , b : 2 } }
End
)";
    // Once the first pair is read without colon (TRUST mode),
    // the colon-and-comma path is not set, so "," is then misread as a key.
    // Either way, this should fail.
    EXPECT_DEATH(interprete_inline_jdd(jdd), "");
}

TEST(Param_Death, Map_Invalid_Separator_JSON_Syntax)
{
    std::string jdd = R"(
Test_Param t
Lire t { map_ints { a : 1 ; b : 2 } }
End
)";
    EXPECT_DEATH(interprete_inline_jdd(jdd), "");
}
