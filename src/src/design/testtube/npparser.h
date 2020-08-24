/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_YY_NPPARSER_H_INCLUDED
# define YY_YY_NPPARSER_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif
/* "%code requires" blocks.  */
#line 14 "tubeparser.y" /* yacc.c:1915  */

#define YYSTYPE value_struc_t*
#ifndef YY_TYPEDEF_YY_SCANNER_T
#define YY_TYPEDEF_YY_SCANNER_T
typedef void* yyscan_t;
#endif

#line 52 "npparser.h" /* yacc.c:1915  */

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    NEWLINE = 258,
    TOK_PERCENT = 259,
    EQUALS = 260,
    PERIOD = 261,
    LBRACKET = 262,
    RBRACKET = 263,
    PLUS = 264,
    LPAREN = 265,
    RPAREN = 266,
    DUPLEX = 267,
    UNPAIRED = 268,
    TOK_NAME = 269,
    TOK_RESERVED = 270,
    TOK_STRUCTURE = 271,
    TOK_DOMAIN = 272,
    TOK_STRAND = 273,
    TOK_TUBE = 274,
    TOK_STAR = 275,
    TOK_INTEGER = 276,
    TOK_FLOAT = 277,
    TOK_NUC = 278,
    TOK_SEED = 279,
    TOK_MAXSIZE = 280,
    TOK_MATERIAL = 281,
    TOK_SODIUM = 282,
    TOK_MAGNESIUM = 283,
    TOK_TEMPERATURE = 284,
    TOK_MUNFAVORABLE = 285,
    TOK_MLEAFOPT = 286,
    TOK_MRESEED = 287,
    TOK_FSPLIT = 288,
    TOK_NSPLIT = 289,
    TOK_HSPLIT = 290,
    TOK_DANGLES = 291,
    TOK_SEQ = 292,
    TOK_CONCDEF = 293,
    TOK_STOPDEF = 294,
    TOK_TRUE = 295,
    TOK_FALSE = 296,
    TOK_PRINTLEAVES = 297,
    TOK_PRINTSTEPS = 298,
    TOK_INCLUDE_ALL = 299,
    TOK_FREDECOMP = 300,
    TOK_FREFOCUS = 301,
    TOK_FSTRINGENT = 302,
    TOK_FPASSIVE = 303,
    TOK_GC_INIT = 304,
    TOK_ALLOWWOBBLE = 305,
    TOK_ALLOWMISMATCH = 306,
    TOK_DISABLEMUTWEIGHTS = 307,
    TOK_MINPAIR = 308,
    TOK_TRIALS = 309,
    TOK_OPTTIME = 310,
    TOK_SINGLE_DECOMP = 311,
    TOK_DGCLAMP = 312
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif

/* Location type.  */
#if ! defined YYLTYPE && ! defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE YYLTYPE;
struct YYLTYPE
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
};
# define YYLTYPE_IS_DECLARED 1
# define YYLTYPE_IS_TRIVIAL 1
#endif



int yyparse (value_struc_t ** expression, yyscan_t scanner);

#endif /* !YY_YY_NPPARSER_H_INCLUDED  */
