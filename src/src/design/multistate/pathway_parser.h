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

#ifndef YY_YY_PATHWAY_PARSER_H_INCLUDED
# define YY_YY_PATHWAY_PARSER_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif
/* "%code requires" blocks.  */
#line 10 "pathway_parser.y" /* yacc.c:1915  */

#ifndef YY_TYPEDEF_YY_SCANNER_T
#define YY_TYPEDEF_YY_SCANNER_T
typedef void* yyscan_t;
#endif

#line 51 "pathway_parser.h" /* yacc.c:1915  */

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    NEWLINE = 258,
    TOK_PERCENT = 259,
    EQUALS = 260,
    COMMA = 261,
    PERIOD = 262,
    LBRACKET = 263,
    RBRACKET = 264,
    LBRACE = 265,
    RBRACE = 266,
    PLUS = 267,
    MINUS = 268,
    LPAREN = 269,
    RPAREN = 270,
    DUPLEX = 271,
    UNPAIRED = 272,
    TOK_NAME = 273,
    TOK_RESERVED = 274,
    TOK_STRUCTURE = 275,
    TOK_DOMAIN = 276,
    TOK_STRAND = 277,
    TOK_TUBE = 278,
    TOK_STAR = 279,
    TOK_INTEGER = 280,
    TOK_FLOAT = 281,
    TOK_NUC = 282,
    TOK_SEED = 283,
    TOK_MAXSIZE = 284,
    TOK_MATERIAL = 285,
    TOK_SODIUM = 286,
    TOK_MAGNESIUM = 287,
    TOK_TEMPERATURE = 288,
    TOK_MBAD = 289,
    TOK_MREOPT = 290,
    TOK_MRESEED = 291,
    TOK_FSPLIT = 292,
    TOK_NSPLIT = 293,
    TOK_HSPLIT = 294,
    TOK_DANGLES = 295,
    TOK_DGCLAMP = 296,
    TOK_CONCDEF = 297,
    TOK_STOPDEF = 298,
    TOK_GLOBAL_STOPDEF = 299,
    TOK_TRUE = 300,
    TOK_FALSE = 301,
    TOK_PRINTLEAVES = 302,
    TOK_PRINTSTEPS = 303,
    TOK_REDECOMPOSE = 304,
    TOK_FSTRINGENT = 305,
    TOK_FREDECOMP = 306,
    TOK_FREFOCUS = 307,
    TOK_FPASSIVE = 308,
    TOK_GC_INIT = 309,
    TOK_ALLOWWOBBLE = 310,
    TOK_ALLOWMISMATCH = 311,
    TOK_DISABLEMUTWEIGHTS = 312,
    TOK_MINPAIR = 313,
    TOK_TRIALS = 314,
    TOK_OPTTIME = 315,
    TOK_PREVENT = 316,
    TOK_LIBRARY = 317,
    TOK_LIBSEQ = 318,
    TOK_POPULATION = 319,
    TOK_SYMMETRY_MIN = 320,
    TOK_WORDSIZE = 321,
    TOK_SOURCE = 322,
    TOK_WINDOW = 323,
    TOK_EXCLUDE = 324,
    TOK_ACCESSIBLE = 325,
    TOK_SIMILARITY = 326,
    TOK_COMPLEMENTARY = 327,
    TOK_IDENTICAL = 328,
    TOK_MATCH = 329,
    TOK_MATCHRANGE = 330,
    TOK_WEIGHT = 331,
    TOK_COMPLEX = 332,
    TOK_OFFTARGETS = 333
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef value_struc_t * YYSTYPE;
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

#endif /* !YY_YY_PATHWAY_PARSER_H_INCLUDED  */
