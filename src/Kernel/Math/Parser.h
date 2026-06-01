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

#ifndef Parser_included
#define Parser_included

#include <UnaryFunction.h>
#include <TRUST_Deriv.h>
#include <TRUST_List.h>
#include <Constante.h>
#include <algorithm>
#include <Stack.h>
#include <math.h>
#include <string>
#include <Noms.h>
#include <TRUSTArray.h>

class StringTokenizer;
// FUNCTION start at 1, cause 0 is needed
enum class FUNCTION { SIN=1, ASIN, COS, ACOS, TAN, ATAN, LN, EXP, SQRT, ENT, ERF, RND, COSH,  SINH, TANH, ATANH, NOT, ABS, SGN };
/*! @brief Representation of data for the Parser class.
 *
 * @sa =
 */

class Parser
{
public :

  /**
   * Initialises the parser with the string "0": serves no purpose!!
   */
  Parser();
  Parser(const Parser&);
  virtual ~Parser();
  void set(const Parser&);

  /**
   * Constructs a Parser object with a specified string and a maximum number of variables to be indicated with the addVar method.
   */
  Parser(std::string&,int n =1);

  /**
   *
   */
  void init_parser();

  /**
   * Builds the tree corresponding to the character string. This tree must be built once and the string is evaluated by traversing this tree via the eval() method as many times as desired.
   */
  virtual void parseString();

  /**
   * Used to evaluate the mathematical expression corresponding to the character string. To do so, the tree must first be built using the parseString() method.
   */
  inline double eval() { return eval(PNodes[0]); }

  /**
   * Sets the value of the variable represented by the string sv.
   */
  inline void setVar(const char* sv, double val) { setVar(searchVar(sv),val); }

  /**
  * Sets the value of the variable represented by v.
  */
  inline void setVar(const std::string& v, double val) { setVar(searchVar(v),val); }

  /**
  * Sets the value of the variable with the specified number. This number corresponds to the order in which variables were added via the addVar() method.
  */
  inline void setVar(int i, double val)
  {
    assert(i>-1 && i<ivar);
    les_var[i] = val;
  }

  /**
   * Sets the maximum number of variables to be specified with the addVar method.
   */
  virtual void setNbVar(int nvar);

  /**
   * Allows adding a variable by specifying its representative string (e.g.: x, y1, etc.)
   */
  void addVar(const char *);

  inline int getmaxVar()     { return maxvar; }
  inline int getNbVar()      { return ivar; }
  inline std::string& getString() { return str; }
  inline void setString(const std::string& s)
  {
    str = s;
  }
  void addCst(const Constante& cst);

  /**
   * Sets the initial time and period of the impulse function.
   */
  inline void setImpulsion(double tinit, double periode)
  {
    impuls_t0 = tinit;
    impuls_T = periode;
    impuls_tn = tinit-periode;
    impuls_tempo = impuls_t0;
  }

protected:
  int test_op_binaire(int type);

  static int precedence(int);
  inline double eval(const PNodePod& node)
  {
    double x,y;
    switch(node.type)
      {
      case 1 :
        x = (node.left  != -1 ? eval(PNodes[node.left])  : 0);
        y = (node.right != -1 ? eval(PNodes[node.right]) : 0);
        return evalOp(node, x, y);  // PNode_type::OP
      case 2 :
        return node.nvalue;  // PNode_type::VALUE
      case 3 :
        x = node.value<=0 ? eval(PNodes[node.left]) : 0;
        return evalFunc(node, x);// PNode_type::FUNCTION
      case 4 :
        return les_var[node.value]; // PNode_type::VAR
      default:
        Process::exit("method eval : Unknown type for this node !!!");
        return 0;
      }
  }
  KOKKOS_INLINE_FUNCTION double evalOp(const PNodePod& node, double x, double y);
  KOKKOS_INLINE_FUNCTION double evalFunc(const PNodePod& node, double x);
  void parserState0(StringTokenizer*,PSTACK(PNode)* ,STACK(int)*);
  void parserState1(StringTokenizer*,PSTACK(PNode)* ,STACK(int)*);
  void parserState2(StringTokenizer*,PSTACK(PNode)* ,STACK(int)*);
  inline int searchVar(const std::string& s) { return searchVar(s.c_str()); }
  inline int searchVar(const char*);
  int searchCst(const std::string& v);
  int searchFunc(const std::string& v);
  int state;
  Constante c_pi ;
  double impuls_T;
  double impuls_t0;
  double impuls_tn;
  double impuls_tempo;

  PNode* root;                      // Linked list of PNode
  std::vector<PNodePod> PNodes;     // Vector of PNodePod
  std::string str;
  ArrOfDouble les_var;
  Noms les_var_names;
  LIST(Constante) les_cst;
  std::map<std::string, int> map_function_;
  int maxvar,ivar;
};

inline int Parser::searchVar(const char * sv)
{
  std::string s(sv);
  std::transform(s.begin(), s.end(), s.begin(), ::toupper);
  return les_var_names.rang(s.c_str());
}

KOKKOS_INLINE_FUNCTION double Parser::evalFunc(const PNodePod& node, double x)
{
  /* OC: New version: */
  if (node.value<=0)
    {
      int unary_function = -node.value-1;  // OC note: in node->value the stored value is the opposite of the function index in the list
      // in order to distinguish binary operators (>0) from unary functions (<0)
      // It is therefore necessary to take -node->value here to reference an element of the list
      // Furthermore, +1 is added because zero must not be used for functions
      switch (unary_function)
        {
        case static_cast<int>(FUNCTION::SIN):
          return sin(x);
        case static_cast<int>(FUNCTION::ASIN):
          return asin(x);
        case static_cast<int>(FUNCTION::COS):
          return cos(x);
        case static_cast<int>(FUNCTION::ACOS):
          return acos(x);
        case static_cast<int>(FUNCTION::TAN):
          return tan(x);
        case static_cast<int>(FUNCTION::ATAN):
          return atan(x);
        case static_cast<int>(FUNCTION::LN):
          if (x <= 0)
            Process::Kokkos_exit("Negative value x for LN(x) function used.\nCheck your data file.");
          return log(x);
        case static_cast<int>(FUNCTION::EXP):
          return exp(x);
        case static_cast<int>(FUNCTION::SQRT):
          if (x < 0)
            Process::Kokkos_exit("Negative value x for SQRT(x) function used.\nCheck your data file.");
          return sqrt(x);
        case static_cast<int>(FUNCTION::ENT):
          return (int) x;
        case static_cast<int>(FUNCTION::ERF):
#ifndef MICROSOFT
          return erf(x);
#else
          Process::exit("erf(x) fonction not implemented on Windows version.");
          return eval(x);
#endif
        case static_cast<int>(FUNCTION::RND):
#ifdef TRUST_USE_GPU
          // ToDo Kokkos: cuRAND and hipRAND
          Process::Kokkos_exit("Rnd function is not available yet for TRUST GPU version.");
          return 0;
#else
          return x*drand48();
#endif
        case static_cast<int>(FUNCTION::COSH):
          return cosh(x);
        case static_cast<int>(FUNCTION::SINH):
          return sinh(x);
        case static_cast<int>(FUNCTION::TANH):
          return tanh(x);
        case static_cast<int>(FUNCTION::ATANH):
          return atanh(x);
        case static_cast<int>(FUNCTION::NOT):
          if (x == 0) return 1;
          else return 0;
        case static_cast<int>(FUNCTION::ABS):
          return std::fabs(x);
        case static_cast<int>(FUNCTION::SGN):
          return (x > 0) - (x < 0);
        default:
          Process::Kokkos_exit("method evalFunc : Unknown function for this node !!!");
          return 0;
        }
    }
  else
    {
      Process::Kokkos_exit("method evalFunc : Unknown func !!!");
      return -1;
    }
}

// Do not inline, otherwise Parser::eval(PNode* node), which is called even more frequently,
// might not be inlined either...
KOKKOS_INLINE_FUNCTION
double Parser::evalOp(const PNodePod& node, double x, double y)
{
  switch (node.value)
    {
    case 0: // ADD
      return x + y;
    case 1: // SUBTRACT
      return  x - y;
    case 2: // MULTIPLY
      return x * y;
    case 3: // DIVIDE
      if (y==0)
        {
          Process::Kokkos_exit("Error in the Parser: x/y calculated with y equals 0. You are using a formulae with a division per 0.");
        }
      return x/y;
    case 4: // POWER
      if (y != (int)(y) && x<0)
        {
#ifdef TRUST_USE_GPU
          Process::Kokkos_exit("Error in the Parser: x^y calculated with negative value for x !");
#else
          Cerr << "Error in the Parser: x^y calculated with negative value for x (x = " << x << ") and y real y (y = " << y << " )" << finl;
          Process::exit();
#endif
        }
      return pow(x,y);
    case 5: // LT
      return (x<y)?1:0;
    case 6: // GT
      return (x>y)?1:0;
    case 7: // LE
      return (x<=y)?1:0;
    case 8: // GE
      return (x>=y)?1:0;
    case 9: // MOD
      return ((int)(x))%((int)(y));
    case 10: // MAX
      return (x>y)?x:y;
    case 11: // MIN
      return (x<y)?x:y;
    case 12: // AND
      return x && y;
    case 13: // OR
      return x || y;
    case 14: // EQ
      return (x == y);
    case 15: // NEQ
      return (x != y);
    default:
#ifdef TRUST_USE_GPU
      Process::Kokkos_exit("Method evalOp : Unknown operation during expression parsing!");
#else
      Cerr << "Method evalOp : Unknown op " << (int)node.value << "!!!" << finl;
      Process::exit();
#endif
      return 0;
    }
}
#endif
