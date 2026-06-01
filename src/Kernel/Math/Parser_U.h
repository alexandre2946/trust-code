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

#ifndef Parser_U_included
#define Parser_U_included

#include <Parser.h>



/*! @brief class Parser_U Version of the Parser class, deriving from Objet_U.
 *
 *      It allows during its use to benefit from the memory management specific to Objet_U
 *      unlike the Math/Parser class
 *
 *
 * @sa Parser
 */
class Parser_U : public Objet_U
{
  Declare_instanciable_sans_constructeur_ni_destructeur(Parser_U);

public :

  Parser_U();
  Parser_U(const Parser_U&);
  ~Parser_U() override;

  const Parser_U& operator=(const Parser_U&);



  /**
   * Constructs the tree corresponding to the character string. This tree must be constructed only once and the character string is evaluated by traversing this tree using the eval() method as many times as desired.
   */
  inline void parseString();


  /**
   * Used to evaluate the mathematical expression corresponding to the character string. For this, you must first construct the tree using the parseString() method.
   */
  inline double eval();

  /**
   * Sets the value of the variable represented by a string sv.
   */
  inline void setVar(const char* sv, double val);

  /**
   * Returns the maximum number of fixed variables.
   */
  inline int getmaxVar();

  /**
   * Returns the number of registered variables.
   */
  inline int getNbVar();

  /**
   * Sets the value of the variable represented by the string v.
   */
  inline void setVar(const std::string& v, double val);

  /**
   * Sets the value of the variable with the specified number. This number corresponds to the order in which variables were added via the addVar() method.
   */
  inline void setVar(int i, double val);
#ifdef INT_is_64_
  inline void setVar(trustIdType i, double val)
  {
    parser->setVar((int)i, val);
  };
#endif
  /**
   * Sets the maximum number of variables to be specified with the addVar method.
   */
  inline void setNbVar(int nvar);

  /**
   * Allows adding a variable by specifying its representative string (e.g.: x, y1, etc.)
   */
  inline void addVar(const char *v);

  inline std::string& getString() ;
  inline void setString(const std::string& s) ;
  inline void setString(const Nom& nom) ;
  inline void addCst(const Constante& cst);
  inline void setImpulsion(double tinit, double periode);
  Parser& getParser() { return *parser; }



private :

  Parser *parser;

};



/**
 * Builds the tree corresponding to the character string. This tree must be built once and the string is evaluated by traversing this tree via the eval() method as many times as desired.
 */
inline void Parser_U::parseString()
{
  parser->parseString();
}


/**
 * Used to evaluate the mathematical expression corresponding to the character string. For this, the tree must first be built using the parseString() method.
 */
inline double Parser_U::eval()
{
  return parser->eval();
}

/**
 * Returns the maximum number of fixed variables.
 */
inline int Parser_U::getmaxVar()
{
  return parser->getmaxVar();
}

/**
 * Returns the number of registered variables.
 */
inline int Parser_U::getNbVar()
{
  return parser->getNbVar();
}

/**
 * Sets the value of the variable represented by the string sv.
 */
inline void Parser_U::setVar(const char* sv, double val)
{
  parser->setVar(sv, val);
}

/**
 * Sets the value of the variable represented by the string v.
 */
inline void Parser_U::setVar(const std::string& v, double val)
{
  parser->setVar(v, val);
}

/**
 * Sets the value of the variable with the specified number. This number corresponds to the order in which variables were added via the addVar() method.
 */
inline void Parser_U::setVar(int i, double val)
{
  parser->setVar(i, val);
}


/**
 * Sets the maximum number of variables to be specified with the addVar method.
 */
inline void Parser_U::setNbVar(int nvar)
{
  parser->setNbVar(nvar);
}


/**
 * Allows adding a variable by specifying its representative string (e.g.: x, y1, etc.)
 */
inline void Parser_U::addVar(const char *v)
{
  parser->addVar(v);
}

inline std::string& Parser_U::getString()
{
  return parser->getString();
}


inline void Parser_U::setString(const std::string& s)
{
  parser->setString(s);
}

inline void Parser_U::setString(const Nom& nom)
{
  const char *s =  nom.getChar();
  std::string ss(s);
  std::transform(ss.begin(), ss.end(), ss.begin(), ::toupper);
  setString(ss);
}

inline void Parser_U::addCst(const Constante& cst)
{
  parser->addCst(cst);
}

inline void Parser_U::setImpulsion(double tinit, double periode)
{
  parser->setImpulsion(tinit, periode);
}





#endif
