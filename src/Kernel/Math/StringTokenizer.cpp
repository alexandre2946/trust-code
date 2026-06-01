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

#include <StringTokenizer.h>
#include <Char_ptr.h>
#include <Objet_U.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <algorithm>

//using namespace std;
using std::stringstream;

const int StringTokenizer::NUMBER=-2;
const int StringTokenizer::STRING=-3;
const int StringTokenizer::EOS=-1;


// Operator identifiers are defined only once, shared between the Parser and StringTokenizer classes.
// They are placed here and removed from the Parser class; they could be moved to a dedicated class if clearer.
/* Note: the following constants must be greater than 0 to avoid conflicts with unary function identifiers.
 */
const int StringTokenizer::ADD = 0;
const int StringTokenizer::SUBTRACT = 1;
const int StringTokenizer::MULTIPLY = 2;
const int StringTokenizer::DIVIDE = 3;
const int StringTokenizer::POWER = 4;
const int StringTokenizer::LT = 5;
const int StringTokenizer::GT = 6;
const int StringTokenizer::LE = 7;
const int StringTokenizer::GE = 8;
const int StringTokenizer::MOD = 9;
const int StringTokenizer::MAX = 10;
const int StringTokenizer::MIN = 11;
const int StringTokenizer::AND = 12;
const int StringTokenizer::OR = 13;
const int StringTokenizer::EQ = 14;
const int StringTokenizer::NEQ = 15;
// Parentheses:
const int StringTokenizer::GRP = 1000;
const int StringTokenizer::ENDGRP = 1001;


// Number of operators and keywords (parentheses "(" and ")" are not included).
// The order of operators below must match the values of the static constants ADD, SUB, etc. defined above.
//
const int StringTokenizer::nb_op=16;
const int StringTokenizer::nb_op_bis=11;
const char StringTokenizer::keyword_op[][10] = { "ADD", "SUB", "MUL", "DIV" , "POW", "LT", "GT", "LE", "GE", "MOD", "MAX", "MIN", "AND", "OR", "EQ", "NEQ" };
const char StringTokenizer::keyword_op_bis[][10] = { "+", "-", "*", "/", "^", "<", ">", "[", "]", "%", "$" };


StringTokenizer::StringTokenizer()
{
  str = std::string("0");
  init_keyword_op();
  reste = &str[0];
}

StringTokenizer::StringTokenizer(std::string& s)
{
  str = s;
  init_keyword_op();
  reste = &str[0];
}

StringTokenizer::StringTokenizer(std::string s, std::string sep)
{
  str = s;
  init_keyword_op();
  reste = &str[0];
}

StringTokenizer::~StringTokenizer()
{

  for (int i=0; i<StringTokenizer::nb_op; i++)
    {
      delete[] op_sep[i];
    }
  delete[] op_sep;

}



int StringTokenizer::check_GRP()
{
  const char* ch = str.c_str();
  int nb_o=0;
  int nb_f=0;
  int sz = (int)strlen(ch);
  for (int i=0; i<sz; i++)
    {
      if (*(ch+i)=='(')
        nb_o++;
      else if (*(ch+i)==')')
        nb_f++;
    }
  return (nb_o==nb_f);
}


int StringTokenizer::nextToken()
{
  char *tmp;
  int type_sep, length;
  tmp=find_sep(reste, type_sep, length);
  if (tmp == nullptr)
    {
      //Cout << "Fin ? " << reste << finl;
      if (reste[0] == '\0')
        {
          type = EOS;
        }
      else
        {
          if ( ((reste[0] >= '0') && ( reste[0] <= '9')) || (reste[0] == '.')  )
            {
              // Added by OC to raise an error when input is like "2x"
              int ind=(int)strlen(reste)-1;
              if ( ((reste[ind] < '0') || (reste[ind] > '9')) && (reste[ind] != '.') )
                {
                  Cerr << "The syntax " << reste << " is not allowed." << finl;
                  Process::exit();
                }
              type = NUMBER;
              stringstream stream;
              stream << reste ;
              stream >> nval;
              if (!stream.eof() )
                {
                  Cerr<<"Error conversion "<<reste<<" in number"<<finl;
                  Process::exit();
                }
            }
          else
            {
              type=STRING;
              sval = reste;
              std::transform(sval.begin(), sval.end(), sval.begin(), ::toupper);
            }
          while ((*reste++) != '\0') ;
          reste--;
        }
    }
  else if (tmp==reste)
    {
      type = type_sep;
      reste=reste+length;
    }
  else
    {
      Char_ptr token;
      token.allocate((int)(tmp-reste));

      char* tok = token.getChar();;
      int j=0;
      for (int i=0; i<tmp-reste; i++)
        {
          if ( reste[i] != ' ' )        tok[j++] = reste[i];
        }
      tok[j] = '\0';
      if (((tok[0] >= '0') && ( tok[0] <= '9')) || (tok[0] == '.')  )
        {
          int ind;
          double nval_tmp;

          type = NUMBER;
          stringstream stream;
          ind = (int)(tmp-reste-1l);
          reste = tmp;
          if (tok[ind] == 'e' || tok[ind] == 'E' )
            {
              tok[ind]=' ';
              stream << tok ;
              stream >> nval_tmp;
              char c;
              stream >> c;

              if (!stream.eof() )
                {
                  Cerr<<"Error conversion "<<tok<<" in number"<<finl;
                  Process::exit();
                }
              nextToken();
              if (type == SUBTRACT)
                {
                  nextToken();
                  if (type != NUMBER)
                    {
                      Cerr << "Error while interpreting the string " << str << finl;
                      Process::exit();
                    }
                  nval_tmp*=pow(10,-nval);
                }
              else if (type == ADD)
                {
                  nextToken();
                  if (type != NUMBER)
                    {
                      Cerr << "Error while interpreting the string " << str << finl;
                      Process::exit();
                    }
                  nval_tmp*=pow(10.,nval);
                }
              else if (type == NUMBER)
                {
                  // GF: depending on the version, the exponent is either parsed as below
                  nval_tmp*=pow(10.,nval);
                  // or as:
                  //nval_tmp = nval;
                  // Disallowed for now:
                  //Cerr << "Possible error while interpreting the string " << str << finl;
                  //Process::exit();
                }
              else
                {
                  Cerr << "Error while interpreting the string " << str << finl;
                  Process::exit();
                }
            }
          else if ( ((tok[ind] < '0') || (tok[ind] > '9')) && (tok[ind] != '.') )
            {
              Cerr << "  The syntax " << tok << " is not allowed." << finl;
              Process::exit();
              throw;
            }
          else
            {
              stream << tok ;
              stream >> nval_tmp;
              if (!stream.eof() )
                {

                  Cerr<<"Error conversion "<<tok<<" in number"<<finl;
                  Process::exit();
                }
            }
          nval = nval_tmp;
        }
      else
        {
          type=STRING;
          sval = tok;
          std::transform(sval.begin(), sval.end(), sval.begin(), ::toupper);
          reste=tmp;
        }
      // delete[] tok;
    }
  return type;
}



// Private methods:

void StringTokenizer::init_keyword_op()
{
  op_sep = new char*[nb_op];
  // Each operator keyword is surrounded by "_":
  // LT =>  _LT_  etc...
  for (int i=0; i<nb_op; i++)
    {
      const char* blanc="_";
      op_sep[i] = new char[strlen(keyword_op[i])+3];
      strcpy(op_sep[i],blanc);
      strcat(op_sep[i], keyword_op[i]);
      strcat(op_sep[i], blanc);
    }
}


/**
 * Searches the string "ch" for the next occurrence of a separator.
 * Returns the separator type ((, ), +, -, etc.) in the "type_sep" parameter.
 * Also returns the length of the string corresponding to the separator found.
 */
char* StringTokenizer::find_sep(char* ch, int& type_sep, int& length)
{
  char * trouve=nullptr;
  char * trouve_tmp;
  type_sep=-1;
  trouve_tmp = strstr(ch, "(");
  int pos=100000;
  if (trouve_tmp != nullptr)
    {
      pos = (int)(trouve_tmp-ch);
      type_sep=GRP;
      length=1;
      trouve=trouve_tmp;
    }
  trouve_tmp = strstr(ch, ")");
  if ((trouve_tmp != nullptr) && (trouve_tmp-ch<pos))
    {
      pos = (int)(trouve_tmp-ch);
      type_sep=ENDGRP;
      length=1;
      trouve=trouve_tmp;
    }
  for (int i=0; i<StringTokenizer::nb_op; i++)
    {
      trouve_tmp = strstr(ch, op_sep[i]);
      if ((trouve_tmp != nullptr) && (trouve_tmp-ch<pos))
        {
          pos = (int)(trouve_tmp-ch);
          type_sep=i;
          length=(int)strlen(op_sep[i]);
          trouve=trouve_tmp;
        }
    }
  for (int i=0; i<StringTokenizer::nb_op_bis; i++)
    {
      trouve_tmp = strstr(ch, keyword_op_bis[i]);
      if ((trouve_tmp != nullptr) && ((int)(trouve_tmp-ch)<pos))
        {
          pos = (int)(trouve_tmp-ch);
          type_sep=i;
          length=(int)strlen(keyword_op_bis[i]);
          trouve=trouve_tmp;
        }
    }

  return trouve;
}



/*
  int main()
  {
  //String s("23+34*12-13+COS ( 12 )* 2^3");
  std::string s("2+3");
  StringTokenizer tk(s);
  Cout << StringTokenizer::NUMBER << finl;
  Cout << StringTokenizer::EOS << finl;
  Cout << StringTokenizer::STRING << finl;

  while (tk.nextToken()!=StringTokenizer::EOS)
  {
  if (tk.type == StringTokenizer::STRING)
  {
  Cout << "String = " << tk.getSValue() << finl;
  }
  else if (tk.type == StringTokenizer::NUMBER)
  {
  Cout << "Value = " << tk.getNValue() << finl;
  }
  else
  {
  Cout << "Operator = " << (char) tk.type << finl;
  }
  }
  }
*/
