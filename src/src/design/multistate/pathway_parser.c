/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 2

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 1 "pathway_parser.y" /* yacc.c:339  */

#include "parsestruc.h"
#include "pathway_parser.h" 
#include "pathway_lexer.h"

int LINE_NUMBER = 1;
int yyerror(YYLTYPE * yyllocp, value_struc_t ** val, yyscan_t scanner, const char * msg);

#line 75 "pathway_parser.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 1
#endif

/* In a future release of Bison, this section will be replaced
   by #include "pathway_parser.h".  */
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
#line 10 "pathway_parser.y" /* yacc.c:355  */

#ifndef YY_TYPEDEF_YY_SCANNER_T
#define YY_TYPEDEF_YY_SCANNER_T
typedef void* yyscan_t;
#endif

#line 112 "pathway_parser.c" /* yacc.c:355  */

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

/* Copy the second part of user declarations.  */

#line 227 "pathway_parser.c" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if 1

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
# define YYCOPY_NEEDED 1
#endif


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL \
             && defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
  YYLTYPE yyls_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE) + sizeof (YYLTYPE)) \
      + 2 * YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  84
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   324

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  79
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  62
/* YYNRULES -- Number of rules.  */
#define YYNRULES  156
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  240

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   333

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   114,   114,   122,   126,   130,   134,   141,   142,   143,
     144,   145,   146,   147,   148,   149,   150,   151,   155,   159,
     166,   167,   171,   177,   181,   185,   189,   193,   200,   207,
     211,   212,   216,   222,   227,   234,   244,   250,   253,   254,
     258,   259,   260,   262,   267,   278,   294,   305,   313,   317,
     326,   327,   341,   342,   349,   350,   357,   365,   366,   373,
     374,   375,   379,   388,   389,   393,   394,   402,   407,   415,
     417,   419,   423,   429,   430,   434,   442,   449,   454,   459,
     467,   469,   471,   473,   479,   480,   487,   488,   495,   496,
     503,   508,   510,   512,   514,   521,   522,   523,   524,   528,
     536,   538,   543,   545,   550,   555,   565,   567,   569,   571,
     573,   575,   577,   579,   581,   583,   585,   587,   589,   591,
     593,   595,   597,   599,   601,   603,   605,   607,   612,   613,
     620,   626,   628,   630,   632,   634,   639,   641,   646,   648,
     650,   656,   661,   669,   671,   676,   682,   687,   694,   703,
     708,   709,   711,   716,   721,   725,   734
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NEWLINE", "TOK_PERCENT", "EQUALS",
  "COMMA", "PERIOD", "LBRACKET", "RBRACKET", "LBRACE", "RBRACE", "PLUS",
  "MINUS", "LPAREN", "RPAREN", "DUPLEX", "UNPAIRED", "TOK_NAME",
  "TOK_RESERVED", "TOK_STRUCTURE", "TOK_DOMAIN", "TOK_STRAND", "TOK_TUBE",
  "TOK_STAR", "TOK_INTEGER", "TOK_FLOAT", "TOK_NUC", "TOK_SEED",
  "TOK_MAXSIZE", "TOK_MATERIAL", "TOK_SODIUM", "TOK_MAGNESIUM",
  "TOK_TEMPERATURE", "TOK_MBAD", "TOK_MREOPT", "TOK_MRESEED",
  "TOK_FSPLIT", "TOK_NSPLIT", "TOK_HSPLIT", "TOK_DANGLES", "TOK_DGCLAMP",
  "TOK_CONCDEF", "TOK_STOPDEF", "TOK_GLOBAL_STOPDEF", "TOK_TRUE",
  "TOK_FALSE", "TOK_PRINTLEAVES", "TOK_PRINTSTEPS", "TOK_REDECOMPOSE",
  "TOK_FSTRINGENT", "TOK_FREDECOMP", "TOK_FREFOCUS", "TOK_FPASSIVE",
  "TOK_GC_INIT", "TOK_ALLOWWOBBLE", "TOK_ALLOWMISMATCH",
  "TOK_DISABLEMUTWEIGHTS", "TOK_MINPAIR", "TOK_TRIALS", "TOK_OPTTIME",
  "TOK_PREVENT", "TOK_LIBRARY", "TOK_LIBSEQ", "TOK_POPULATION",
  "TOK_SYMMETRY_MIN", "TOK_WORDSIZE", "TOK_SOURCE", "TOK_WINDOW",
  "TOK_EXCLUDE", "TOK_ACCESSIBLE", "TOK_SIMILARITY", "TOK_COMPLEMENTARY",
  "TOK_IDENTICAL", "TOK_MATCH", "TOK_MATCHRANGE", "TOK_WEIGHT",
  "TOK_COMPLEX", "TOK_OFFTARGETS", "$accept", "input", "definitions",
  "definition", "structure_def", "structure", "dpp_struc_list",
  "dpp_struc_el", "dup_struc_list", "empt_dup_struc_list", "dup_struc",
  "dup_single_el", "dup_inside_duplex", "unpaired_el", "duplex_el",
  "library_def", "domain_def", "domain", "domain_el", "domain_list",
  "domain_name", "tube_def", "name_list", "dot_prop_def",
  "rangelist_property", "name_or_blank", "rangelist", "range",
  "rangelist_property_name", "namelist_property", "spacelist_property",
  "commalist_property", "spacelist_property_name",
  "commalist_property_name", "numerical_property",
  "numerical_property_name", "dot_namelist", "namelist", "commalist",
  "domainlist_define_type", "domainlist_type", "global_prop_def",
  "namelist_global", "namelist_global_name", "seqlist_global",
  "seqlist_global_name", "numerical_global", "numerical_global_name",
  "seqlist", "boolean_global", "boolean_global_name", "boolean", "units",
  "source_def", "namelist_relation_def", "namelist_relation_type",
  "match_def", "offtarget_def", "superlist", "set", "offtargets",
  "maxsize", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333
};
# endif

#define YYPACT_NINF -188

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-188)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     187,  -188,    -6,     2,     2,  -188,     7,  -188,  -188,  -188,
    -188,  -188,  -188,  -188,  -188,  -188,  -188,  -188,  -188,  -188,
    -188,  -188,  -188,  -188,  -188,  -188,  -188,  -188,  -188,  -188,
    -188,  -188,  -188,  -188,  -188,  -188,    19,  -188,  -188,    22,
    -188,    30,  -188,  -188,  -188,    56,   247,    70,  -188,  -188,
    -188,    68,  -188,  -188,  -188,  -188,  -188,  -188,  -188,    20,
    -188,     2,  -188,  -188,    75,  -188,    10,  -188,    19,  -188,
      80,  -188,  -188,    19,  -188,  -188,  -188,    83,    90,    93,
      12,     2,    99,   104,  -188,  -188,    79,  -188,   -10,  -188,
    -188,  -188,  -188,  -188,  -188,  -188,  -188,  -188,  -188,   109,
      19,    19,   112,    19,   113,     2,    92,  -188,    28,   115,
      26,     2,    85,    92,     2,   114,   116,   117,    92,    92,
    -188,   122,   124,  -188,     2,   126,     2,   127,     2,  -188,
       2,    96,    92,  -188,   128,    92,     2,   107,  -188,  -188,
    -188,    14,  -188,  -188,  -188,  -188,   110,   118,  -188,    43,
    -188,   125,    77,  -188,  -188,  -188,    92,  -188,     2,  -188,
    -188,    92,    92,    92,    85,   129,  -188,   133,     2,   134,
     119,     2,  -188,  -188,  -188,  -188,    92,   128,  -188,  -188,
       2,    91,  -188,  -188,    77,  -188,  -188,   128,  -188,     8,
     131,  -188,   100,     5,     2,  -188,  -188,    92,     2,  -188,
      35,  -188,  -188,  -188,    77,   137,     2,    54,   136,   136,
     136,   121,  -188,     5,  -188,    77,   135,    62,   123,     2,
    -188,     2,  -188,   138,  -188,   142,  -188,    77,  -188,    77,
     141,     2,   136,   130,    94,  -188,  -188,   144,  -188,  -188
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     3,    54,     0,     0,    92,     0,   110,   100,   108,
     106,   109,   111,   112,   113,   114,   115,   116,   101,   117,
     107,   131,   127,   132,   118,   119,   120,   121,   122,   133,
     134,   135,   123,   124,   125,   104,   138,   126,    91,     0,
      93,     0,   143,   144,    94,     0,     2,     0,     7,    12,
       8,     0,     9,    11,    61,    60,    73,    74,    59,     0,
      13,     0,    10,    97,     0,    98,     0,    95,   138,    96,
       0,    14,    15,   138,    16,    17,    55,     0,     0,     0,
       0,     0,     0,     0,     1,     6,     0,     4,    84,    81,
      80,    82,    77,    72,    78,    69,    70,    71,    83,     0,
     138,   138,     0,   138,     0,     0,     0,    88,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       5,     0,     0,    85,    64,     0,     0,     0,     0,    86,
      99,    50,   128,    48,   102,     0,     0,     0,   136,   137,
     130,     0,    27,    24,    25,    26,     0,     0,    18,    20,
      23,    21,    29,    37,    39,    38,    47,    57,    56,   140,
     139,     0,   141,   145,     0,     0,    63,     0,     0,    76,
       0,    90,    52,    87,    51,    49,     0,   103,    89,   105,
       0,     0,    44,    22,     0,    36,    58,    46,    19,     0,
     154,   146,   150,     0,    75,    79,    53,   129,   142,    42,
       0,    45,    41,    40,    28,     0,   147,     0,     0,     0,
       0,     0,    68,    62,    65,    31,     0,    30,     0,     0,
     149,     0,   155,   151,   152,     0,    66,    32,    43,    33,
       0,   148,     0,     0,    34,   156,   153,     0,    35,    67
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -188,  -188,  -188,   108,  -188,    -9,  -188,    11,  -188,  -188,
    -176,  -150,  -188,   -24,   -22,  -188,  -188,  -107,  -127,  -188,
       0,  -188,  -188,  -188,  -188,  -188,  -188,   -52,  -188,  -188,
    -188,  -188,  -188,  -188,  -188,  -188,  -188,  -104,    36,  -188,
    -188,  -188,  -188,  -188,  -188,  -188,  -188,  -188,  -118,  -188,
    -188,  -188,   -59,  -188,  -188,  -188,  -188,  -188,  -188,  -187,
    -188,  -188
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    45,    46,    47,    48,   148,   149,   150,   151,   216,
     152,   153,   201,   154,   155,    49,    50,   132,   133,   171,
     129,    52,   158,    53,    54,   167,   213,   214,   100,    55,
      56,    57,   101,   102,    58,   103,    59,   206,   108,    60,
      61,    62,    63,    64,    65,    66,    67,    68,   134,    69,
      70,   140,    81,    71,    72,    73,    74,    75,   207,   190,
     191,   192
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
      51,   130,   185,    77,    78,   175,   156,   141,   204,   109,
     121,   162,   163,   211,   111,   106,   115,   177,    76,   180,
       2,   222,   223,   224,   217,    79,     2,    80,     2,   175,
     116,   212,     2,   135,   136,   175,   175,   205,     2,   227,
      82,   124,   125,   187,   127,   236,    51,   215,    83,    89,
     142,   146,   147,   234,   185,   143,    84,   144,   145,    99,
     219,   104,    90,    91,   194,   220,   107,   185,   122,   197,
     175,   138,   139,    87,   229,    88,   198,   185,   146,   147,
     105,   117,   120,    92,   185,   110,    93,    94,   112,    95,
      96,    97,   142,   146,   147,   113,    98,   143,   114,   144,
     145,   146,   147,   199,   118,   200,   238,   146,   147,   119,
     146,   147,   209,   210,   157,   231,   123,   126,   128,   131,
     137,   174,   161,   159,   166,   160,   107,   164,   172,   165,
     173,   168,   170,   179,   176,   181,   178,   184,   193,   189,
     136,   173,   218,   182,   208,   195,   221,   225,   233,   230,
     228,   232,   235,   239,    86,   188,   237,   202,   186,   203,
     183,   226,   169,     0,     0,     0,     0,     0,     0,     0,
       0,   196,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       1,     0,     0,     0,   173,     0,     0,     0,   173,     0,
       0,     0,     0,     0,     0,     2,   173,     3,     4,     5,
       6,     0,     0,     0,     0,     7,     0,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,     0,
      20,   173,     0,     0,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      85,    37,    38,     0,    39,    40,     0,     0,    41,    42,
       0,    43,     0,     0,    44,     2,     0,     3,     4,     5,
       6,     0,     0,     0,     0,     7,     0,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,     0,
      20,     0,     0,     0,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
       0,    37,    38,     0,    39,    40,     0,     0,    41,    42,
       0,    43,     0,     0,    44
};

static const yytype_int16 yycheck[] =
{
       0,   105,   152,     3,     4,   132,   113,   111,   184,    68,
      20,   118,   119,     8,    73,     5,     4,   135,    24,     5,
      18,   208,   209,   210,   200,    18,    18,     8,    18,   156,
      18,    26,    18,     5,     6,   162,   163,    29,    18,   215,
      18,   100,   101,   161,   103,   232,    46,    12,    18,    29,
       7,    16,    17,   229,   204,    12,     0,    14,    15,    59,
       6,    61,    42,    43,   168,    11,    66,   217,    78,   176,
     197,    45,    46,     3,    12,     7,   180,   227,    16,    17,
       5,    81,     3,    63,   234,     5,    66,    67,     5,    69,
      70,    71,     7,    16,    17,     5,    76,    12,     5,    14,
      15,    16,    17,    12,     5,    14,    12,    16,    17,     5,
      16,    17,    12,    13,   114,   219,     7,     5,     5,    27,
       5,    25,     5,     9,   124,     9,   126,     5,   128,     5,
     130,     5,     5,    26,     6,    25,   136,    12,     5,    10,
       6,   141,     5,    25,    13,    26,    10,    26,     6,    26,
      15,    13,    11,     9,    46,   164,    26,   181,   158,   181,
     149,   213,   126,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   171,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
       3,    -1,    -1,    -1,   194,    -1,    -1,    -1,   198,    -1,
      -1,    -1,    -1,    -1,    -1,    18,   206,    20,    21,    22,
      23,    -1,    -1,    -1,    -1,    28,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    -1,
      43,   231,    -1,    -1,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    56,    57,    58,    59,    60,    61,    62,
       3,    64,    65,    -1,    67,    68,    -1,    -1,    71,    72,
      -1,    74,    -1,    -1,    77,    18,    -1,    20,    21,    22,
      23,    -1,    -1,    -1,    -1,    28,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    -1,
      43,    -1,    -1,    -1,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    56,    57,    58,    59,    60,    61,    62,
      -1,    64,    65,    -1,    67,    68,    -1,    -1,    71,    72,
      -1,    74,    -1,    -1,    77
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,    18,    20,    21,    22,    23,    28,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      43,    47,    48,    49,    50,    51,    52,    53,    54,    55,
      56,    57,    58,    59,    60,    61,    62,    64,    65,    67,
      68,    71,    72,    74,    77,    80,    81,    82,    83,    94,
      95,    99,   100,   102,   103,   108,   109,   110,   113,   115,
     118,   119,   120,   121,   122,   123,   124,   125,   126,   128,
     129,   132,   133,   134,   135,   136,    24,    99,    99,    18,
       8,   131,    18,    18,     0,     3,    82,     3,     7,    29,
      42,    43,    63,    66,    67,    69,    70,    71,    76,    99,
     107,   111,   112,   114,    99,     5,     5,    99,   117,   131,
       5,   131,     5,     5,     5,     4,    18,    99,     5,     5,
       3,    20,    78,     7,   131,   131,     5,   131,     5,    99,
     116,    27,    96,    97,   127,     5,     6,     5,    45,    46,
     130,   116,     7,    12,    14,    15,    16,    17,    84,    85,
      86,    87,    89,    90,    92,    93,    96,    99,   101,     9,
       9,     5,    96,    96,     5,     5,    99,   104,     5,   117,
       5,    98,    99,    99,    25,    97,     6,   127,    99,    26,
       5,    25,    25,    86,    12,    90,    99,   127,    84,    10,
     138,   139,   140,     5,   116,    26,    99,    96,   116,    12,
      14,    91,    92,    93,    89,    29,   116,   137,    13,    12,
      13,     8,    26,   105,   106,    12,    88,    89,     5,     6,
      11,    10,   138,   138,   138,    26,   106,    89,    15,    12,
      26,   116,    13,     6,    89,    11,   138,    26,    12,     9
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    79,    80,    81,    81,    81,    81,    82,    82,    82,
      82,    82,    82,    82,    82,    82,    82,    82,    83,    83,
      84,    84,    85,    85,    86,    86,    86,    86,    87,    87,
      88,    88,    88,    88,    88,    88,    89,    89,    90,    90,
      91,    91,    91,    91,    92,    93,    94,    95,    96,    96,
      97,    97,    98,    98,    99,    99,   100,   101,   101,   102,
     102,   102,   103,   104,   104,   105,   105,   106,   106,   107,
     107,   107,   107,   108,   108,   109,   110,   111,   112,   113,
     114,   114,   114,   114,   115,   115,   116,   116,   117,   117,
     118,   119,   119,   119,   119,   120,   120,   120,   120,   121,
     122,   122,   123,   123,   124,   125,   126,   126,   126,   126,
     126,   126,   126,   126,   126,   126,   126,   126,   126,   126,
     126,   126,   126,   126,   126,   126,   126,   126,   127,   127,
     128,   129,   129,   129,   129,   129,   130,   130,   131,   131,
     131,   132,   133,   134,   134,   135,   136,   137,   137,   138,
     139,   139,   139,   139,   139,   139,   140
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     2,     3,     2,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     4,     5,
       1,     1,     2,     1,     1,     1,     1,     1,     3,     1,
       1,     1,     2,     2,     3,     4,     2,     1,     1,     1,
       1,     1,     1,     3,     2,     3,     5,     4,     1,     2,
       1,     2,     1,     2,     1,     2,     4,     1,     2,     1,
       1,     1,     6,     1,     0,     1,     2,     5,     1,     1,
       1,     1,     1,     1,     1,     5,     4,     1,     1,     5,
       1,     1,     1,     1,     2,     3,     1,     2,     1,     3,
       4,     1,     1,     1,     1,     1,     1,     1,     1,     3,
       1,     1,     3,     4,     1,     4,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     3,
       3,     1,     1,     1,     1,     1,     1,     1,     0,     3,
       3,     4,     5,     1,     1,     4,     5,     1,     3,     3,
       1,     3,     3,     5,     1,     3,     5
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      YY_LAC_DISCARD ("YYBACKUP");                              \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (&yylloc, expression, scanner, YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)                                \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;        \
          (Current).first_column = YYRHSLOC (Rhs, 1).first_column;      \
          (Current).last_line    = YYRHSLOC (Rhs, N).last_line;         \
          (Current).last_column  = YYRHSLOC (Rhs, N).last_column;       \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).first_line   = (Current).last_line   =              \
            YYRHSLOC (Rhs, 0).last_line;                                \
          (Current).first_column = (Current).last_column =              \
            YYRHSLOC (Rhs, 0).last_column;                              \
        }                                                               \
    while (0)
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K])


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL

/* Print *YYLOCP on YYO.  Private, do not rely on its existence. */

YY_ATTRIBUTE_UNUSED
static unsigned
yy_location_print_ (FILE *yyo, YYLTYPE const * const yylocp)
{
  unsigned res = 0;
  int end_col = 0 != yylocp->last_column ? yylocp->last_column - 1 : 0;
  if (0 <= yylocp->first_line)
    {
      res += YYFPRINTF (yyo, "%d", yylocp->first_line);
      if (0 <= yylocp->first_column)
        res += YYFPRINTF (yyo, ".%d", yylocp->first_column);
    }
  if (0 <= yylocp->last_line)
    {
      if (yylocp->first_line < yylocp->last_line)
        {
          res += YYFPRINTF (yyo, "-%d", yylocp->last_line);
          if (0 <= end_col)
            res += YYFPRINTF (yyo, ".%d", end_col);
        }
      else if (0 <= end_col && yylocp->first_column < end_col)
        res += YYFPRINTF (yyo, "-%d", end_col);
    }
  return res;
 }

#  define YY_LOCATION_PRINT(File, Loc)          \
  yy_location_print_ (File, &(Loc))

# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value, Location, expression, scanner); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, value_struc_t ** expression, yyscan_t scanner)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (yylocationp);
  YYUSE (expression);
  YYUSE (scanner);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, value_struc_t ** expression, yyscan_t scanner)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  YY_LOCATION_PRINT (yyoutput, *yylocationp);
  YYFPRINTF (yyoutput, ": ");
  yy_symbol_value_print (yyoutput, yytype, yyvaluep, yylocationp, expression, scanner);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, YYLTYPE *yylsp, int yyrule, value_struc_t ** expression, yyscan_t scanner)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                       , &(yylsp[(yyi + 1) - (yynrhs)])                       , expression, scanner);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, yylsp, Rule, expression, scanner); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif

/* Given a state stack such that *YYBOTTOM is its bottom, such that
   *YYTOP is either its top or is YYTOP_EMPTY to indicate an empty
   stack, and such that *YYCAPACITY is the maximum number of elements it
   can hold without a reallocation, make sure there is enough room to
   store YYADD more elements.  If not, allocate a new stack using
   YYSTACK_ALLOC, copy the existing elements, and adjust *YYBOTTOM,
   *YYTOP, and *YYCAPACITY to reflect the new capacity and memory
   location.  If *YYBOTTOM != YYBOTTOM_NO_FREE, then free the old stack
   using YYSTACK_FREE.  Return 0 if successful or if no reallocation is
   required.  Return 1 if memory is exhausted.  */
static int
yy_lac_stack_realloc (YYSIZE_T *yycapacity, YYSIZE_T yyadd,
#if YYDEBUG
                      char const *yydebug_prefix,
                      char const *yydebug_suffix,
#endif
                      yytype_int16 **yybottom,
                      yytype_int16 *yybottom_no_free,
                      yytype_int16 **yytop, yytype_int16 *yytop_empty)
{
  YYSIZE_T yysize_old =
    *yytop == yytop_empty ? 0 : *yytop - *yybottom + 1;
  YYSIZE_T yysize_new = yysize_old + yyadd;
  if (*yycapacity < yysize_new)
    {
      YYSIZE_T yyalloc = 2 * yysize_new;
      yytype_int16 *yybottom_new;
      /* Use YYMAXDEPTH for maximum stack size given that the stack
         should never need to grow larger than the main state stack
         needs to grow without LAC.  */
      if (YYMAXDEPTH < yysize_new)
        {
          YYDPRINTF ((stderr, "%smax size exceeded%s", yydebug_prefix,
                      yydebug_suffix));
          return 1;
        }
      if (YYMAXDEPTH < yyalloc)
        yyalloc = YYMAXDEPTH;
      yybottom_new =
        (yytype_int16*) YYSTACK_ALLOC (yyalloc * sizeof *yybottom_new);
      if (!yybottom_new)
        {
          YYDPRINTF ((stderr, "%srealloc failed%s", yydebug_prefix,
                      yydebug_suffix));
          return 1;
        }
      if (*yytop != yytop_empty)
        {
          YYCOPY (yybottom_new, *yybottom, yysize_old);
          *yytop = yybottom_new + (yysize_old - 1);
        }
      if (*yybottom != yybottom_no_free)
        YYSTACK_FREE (*yybottom);
      *yybottom = yybottom_new;
      *yycapacity = yyalloc;
    }
  return 0;
}

/* Establish the initial context for the current lookahead if no initial
   context is currently established.

   We define a context as a snapshot of the parser stacks.  We define
   the initial context for a lookahead as the context in which the
   parser initially examines that lookahead in order to select a
   syntactic action.  Thus, if the lookahead eventually proves
   syntactically unacceptable (possibly in a later context reached via a
   series of reductions), the initial context can be used to determine
   the exact set of tokens that would be syntactically acceptable in the
   lookahead's place.  Moreover, it is the context after which any
   further semantic actions would be erroneous because they would be
   determined by a syntactically unacceptable token.

   YY_LAC_ESTABLISH should be invoked when a reduction is about to be
   performed in an inconsistent state (which, for the purposes of LAC,
   includes consistent states that don't know they're consistent because
   their default reductions have been disabled).  Iff there is a
   lookahead token, it should also be invoked before reporting a syntax
   error.  This latter case is for the sake of the debugging output.

   For parse.lac=full, the implementation of YY_LAC_ESTABLISH is as
   follows.  If no initial context is currently established for the
   current lookahead, then check if that lookahead can eventually be
   shifted if syntactic actions continue from the current context.
   Report a syntax error if it cannot.  */
#define YY_LAC_ESTABLISH                                         \
do {                                                             \
  if (!yy_lac_established)                                       \
    {                                                            \
      YYDPRINTF ((stderr,                                        \
                  "LAC: initial context established for %s\n",   \
                  yytname[yytoken]));                            \
      yy_lac_established = 1;                                    \
      {                                                          \
        int yy_lac_status =                                      \
          yy_lac (yyesa, &yyes, &yyes_capacity, yyssp, yytoken); \
        if (yy_lac_status == 2)                                  \
          goto yyexhaustedlab;                                   \
        if (yy_lac_status == 1)                                  \
          goto yyerrlab;                                         \
      }                                                          \
    }                                                            \
} while (0)

/* Discard any previous initial lookahead context because of Event,
   which may be a lookahead change or an invalidation of the currently
   established initial context for the current lookahead.

   The most common example of a lookahead change is a shift.  An example
   of both cases is syntax error recovery.  That is, a syntax error
   occurs when the lookahead is syntactically erroneous for the
   currently established initial context, so error recovery manipulates
   the parser stacks to try to find a new initial context in which the
   current lookahead is syntactically acceptable.  If it fails to find
   such a context, it discards the lookahead.  */
#if YYDEBUG
# define YY_LAC_DISCARD(Event)                                           \
do {                                                                     \
  if (yy_lac_established)                                                \
    {                                                                    \
      if (yydebug)                                                       \
        YYFPRINTF (stderr, "LAC: initial context discarded due to "      \
                   Event "\n");                                          \
      yy_lac_established = 0;                                            \
    }                                                                    \
} while (0)
#else
# define YY_LAC_DISCARD(Event) yy_lac_established = 0
#endif

/* Given the stack whose top is *YYSSP, return 0 iff YYTOKEN can
   eventually (after perhaps some reductions) be shifted, return 1 if
   not, or return 2 if memory is exhausted.  As preconditions and
   postconditions: *YYES_CAPACITY is the allocated size of the array to
   which *YYES points, and either *YYES = YYESA or *YYES points to an
   array allocated with YYSTACK_ALLOC.  yy_lac may overwrite the
   contents of either array, alter *YYES and *YYES_CAPACITY, and free
   any old *YYES other than YYESA.  */
static int
yy_lac (yytype_int16 *yyesa, yytype_int16 **yyes,
        YYSIZE_T *yyes_capacity, yytype_int16 *yyssp, int yytoken)
{
  yytype_int16 *yyes_prev = yyssp;
  yytype_int16 *yyesp = yyes_prev;
  YYDPRINTF ((stderr, "LAC: checking lookahead %s:", yytname[yytoken]));
  if (yytoken == YYUNDEFTOK)
    {
      YYDPRINTF ((stderr, " Always Err\n"));
      return 1;
    }
  while (1)
    {
      int yyrule = yypact[*yyesp];
      if (yypact_value_is_default (yyrule)
          || (yyrule += yytoken) < 0 || YYLAST < yyrule
          || yycheck[yyrule] != yytoken)
        {
          yyrule = yydefact[*yyesp];
          if (yyrule == 0)
            {
              YYDPRINTF ((stderr, " Err\n"));
              return 1;
            }
        }
      else
        {
          yyrule = yytable[yyrule];
          if (yytable_value_is_error (yyrule))
            {
              YYDPRINTF ((stderr, " Err\n"));
              return 1;
            }
          if (0 < yyrule)
            {
              YYDPRINTF ((stderr, " S%d\n", yyrule));
              return 0;
            }
          yyrule = -yyrule;
        }
      {
        YYSIZE_T yylen = yyr2[yyrule];
        YYDPRINTF ((stderr, " R%d", yyrule - 1));
        if (yyesp != yyes_prev)
          {
            YYSIZE_T yysize = yyesp - *yyes + 1;
            if (yylen < yysize)
              {
                yyesp -= yylen;
                yylen = 0;
              }
            else
              {
                yylen -= yysize;
                yyesp = yyes_prev;
              }
          }
        if (yylen)
          yyesp = yyes_prev -= yylen;
      }
      {
        int yystate;
        {
          int yylhs = yyr1[yyrule] - YYNTOKENS;
          yystate = yypgoto[yylhs] + *yyesp;
          if (yystate < 0 || YYLAST < yystate
              || yycheck[yystate] != *yyesp)
            yystate = yydefgoto[yylhs];
          else
            yystate = yytable[yystate];
        }
        if (yyesp == yyes_prev)
          {
            yyesp = *yyes;
            *yyesp = yystate;
          }
        else
          {
            if (yy_lac_stack_realloc (yyes_capacity, 1,
#if YYDEBUG
                                      " (", ")",
#endif
                                      yyes, yyesa, &yyesp, yyes_prev))
              {
                YYDPRINTF ((stderr, "\n"));
                return 2;
              }
            *++yyesp = yystate;
          }
        YYDPRINTF ((stderr, " G%d", yystate));
      }
    }
}


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.  In order to see if a particular token T is a
   valid looakhead, invoke yy_lac (YYESA, YYES, YYES_CAPACITY, YYSSP, T).

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store or if
   yy_lac returned 2.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyesa, yytype_int16 **yyes,
                YYSIZE_T *yyes_capacity, yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
       In the first two cases, it might appear that the current syntax
       error should have been detected in the previous state when yy_lac
       was invoked.  However, at that time, there might have been a
       different syntax error that discarded a different initial context
       during error recovery, leaving behind the current lookahead.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      YYDPRINTF ((stderr, "Constructing syntax error message\n"));
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          int yyx;

          for (yyx = 0; yyx < YYNTOKENS; ++yyx)
            if (yyx != YYTERROR && yyx != YYUNDEFTOK)
              {
                {
                  int yy_lac_status = yy_lac (yyesa, yyes, yyes_capacity,
                                              yyssp, yyx);
                  if (yy_lac_status == 2)
                    return 2;
                  if (yy_lac_status == 1)
                    continue;
                }
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
# if YYDEBUG
      else if (yydebug)
        YYFPRINTF (stderr, "No expected tokens.\n");
# endif
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, YYLTYPE *yylocationp, value_struc_t ** expression, yyscan_t scanner)
{
  YYUSE (yyvaluep);
  YYUSE (yylocationp);
  YYUSE (expression);
  YYUSE (scanner);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/*----------.
| yyparse.  |
`----------*/

int
yyparse (value_struc_t ** expression, yyscan_t scanner)
{
/* The lookahead symbol.  */
int yychar;


/* The semantic value of the lookahead symbol.  */
/* Default value used for initialization, for pacifying older GCCs
   or non-GCC compilers.  */
YY_INITIAL_VALUE (static YYSTYPE yyval_default;)
YYSTYPE yylval YY_INITIAL_VALUE (= yyval_default);

/* Location data for the lookahead symbol.  */
static YYLTYPE yyloc_default
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
  = { 1, 1, 1, 1 }
# endif
;
YYLTYPE yylloc = yyloc_default;

    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.
       'yyls': related to locations.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    /* The location stack.  */
    YYLTYPE yylsa[YYINITDEPTH];
    YYLTYPE *yyls;
    YYLTYPE *yylsp;

    /* The locations where the error started and ended.  */
    YYLTYPE yyerror_range[3];

    YYSIZE_T yystacksize;

    yytype_int16 yyesa[20];
    yytype_int16 *yyes;
    YYSIZE_T yyes_capacity;

  int yy_lac_established = 0;
  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
  YYLTYPE yyloc;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N), yylsp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yylsp = yyls = yylsa;
  yystacksize = YYINITDEPTH;

  yyes = yyesa;
  yyes_capacity = sizeof yyesa / sizeof *yyes;
  if (YYMAXDEPTH < yyes_capacity)
    yyes_capacity = YYMAXDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  yylsp[0] = yylloc;
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;
        YYLTYPE *yyls1 = yyls;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yyls1, yysize * sizeof (*yylsp),
                    &yystacksize);

        yyls = yyls1;
        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
        YYSTACK_RELOCATE (yyls_alloc, yyls);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
      yylsp = yyls + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex (&yylval, &yylloc, scanner);
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    {
      YY_LAC_ESTABLISH;
      goto yydefault;
    }
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      YY_LAC_ESTABLISH;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
  YY_LAC_DISCARD ("shift");

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END
  *++yylsp = yylloc;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];

  /* Default location.  */
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
  YY_REDUCE_PRINT (yyn);
  {
    int yychar_backup = yychar;
    switch (yyn)
      {
          case 2:
#line 115 "pathway_parser.y" /* yacc.c:1661  */
    {
        *expression = (yyvsp[0]);
        (yyval) = (yyvsp[0]);
    }
#line 1871 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 3:
#line 122 "pathway_parser.y" /* yacc.c:1661  */
    {
        LINE_NUMBER += 1;
        (yyval) = np_tt_make_string("");
    }
#line 1880 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 4:
#line 126 "pathway_parser.y" /* yacc.c:1661  */
    {
        LINE_NUMBER += 1; 
        (yyval) = (yyvsp[-1]);
    }
#line 1889 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 5:
#line 130 "pathway_parser.y" /* yacc.c:1661  */
    {
        LINE_NUMBER += 1;
        (yyval) = np_tt_append_value((yyvsp[-2]), (yyvsp[-1]));
    }
#line 1898 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 6:
#line 134 "pathway_parser.y" /* yacc.c:1661  */
    {
        LINE_NUMBER += 1;
        (yyval) = np_tt_append_value((yyvsp[-1]), np_tt_make_string(""));
    }
#line 1907 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 18:
#line 156 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_definition(np_tt_make_int(TOK_STRUCTURE), (yyvsp[-2]), (yyvsp[0]));
    }
#line 1915 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 19:
#line 160 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_definition(np_tt_make_int(TOK_STRUCTURE), (yyvsp[-4]), (yyvsp[0]));
    }
#line 1923 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 22:
#line 172 "pathway_parser.y" /* yacc.c:1661  */
    {
        np_tt_append_string((yyvsp[-1]), (yyvsp[0])->strval);
        np_tt_destroy_value_struc((yyvsp[0]));
        (yyval) = (yyvsp[-1]);
    }
#line 1933 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 24:
#line 182 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_string("+");
    }
#line 1941 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 25:
#line 186 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_string("(");
    }
#line 1949 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 26:
#line 190 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_string(")");
    }
#line 1957 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 27:
#line 194 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_string(".");
    }
#line 1965 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 28:
#line 201 "pathway_parser.y" /* yacc.c:1661  */
    {
        np_tt_append_string((yyvsp[-2]), "+");
        np_tt_append_string((yyvsp[-2]), (yyvsp[0])->strval);
        np_tt_destroy_value_struc((yyvsp[0]));
        (yyval) = (yyvsp[-2]);
    }
#line 1976 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 31:
#line 213 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_string("+");
    }
#line 1984 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 32:
#line 217 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_string("+");
        np_tt_append_string((yyval), (yyvsp[0])->strval);
        np_tt_destroy_value_struc((yyvsp[0]));
    }
#line 1994 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 33:
#line 223 "pathway_parser.y" /* yacc.c:1661  */
    {
        np_tt_append_string((yyvsp[-1]), "+");
        (yyval) = (yyvsp[-1]);
    }
#line 2003 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 34:
#line 228 "pathway_parser.y" /* yacc.c:1661  */
    {
        np_tt_append_string((yyvsp[-2]), "+");
        np_tt_append_string((yyvsp[-2]), (yyvsp[0])->strval);
        np_tt_destroy_value_struc((yyvsp[0]));
        (yyval) = (yyvsp[-2]);
    }
#line 2014 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 35:
#line 235 "pathway_parser.y" /* yacc.c:1661  */
    {
        yyerror(&yylloc, NULL, scanner, 
            "Disconnected structure");
        YYERROR;
    }
#line 2024 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 36:
#line 245 "pathway_parser.y" /* yacc.c:1661  */
    {
        np_tt_append_string((yyvsp[-1]), (yyvsp[0])->strval);
        np_tt_destroy_value_struc((yyvsp[0]));
        (yyval) = (yyvsp[-1]);
    }
#line 2034 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 42:
#line 261 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_string("+"); }
#line 2040 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 43:
#line 263 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = (yyvsp[-1]); }
#line 2046 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 44:
#line 268 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_string("");
        for (int el_count = 0; el_count < (yyvsp[0])->intval; el_count++) {
            np_tt_append_string((yyval), ".");
        }
        np_tt_destroy_value_struc((yyvsp[0]));
    }
#line 2058 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 45:
#line 279 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_string("");
        for (int el_count = 0; el_count < (yyvsp[-1])->intval; el_count++) {
            np_tt_append_string((yyval), "(");
        }
        np_tt_append_string((yyval), (yyvsp[0])->strval);
        for (int el_count = 0; el_count < (yyvsp[-1])->intval; el_count++) {
            np_tt_append_string((yyval), ")");
        }
        
        np_tt_destroy_value_struc((yyvsp[-1]));
        np_tt_destroy_value_struc((yyvsp[0]));
    }
#line 2076 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 46:
#line 295 "pathway_parser.y" /* yacc.c:1661  */
    {

        value_struc_t * lib_int = np_tt_make_int(TOK_LIBRARY);
        np_tt_append_value(lib_int, (yyvsp[-3]));

        (yyval) = np_tt_make_definition(lib_int, (yyvsp[-2]), (yyvsp[0]));
    }
#line 2088 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 47:
#line 306 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_definition(np_tt_make_int(TOK_DOMAIN),
                (yyvsp[-2]), (yyvsp[0]));
    }
#line 2097 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 48:
#line 314 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = (yyvsp[0]);
    }
#line 2105 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 49:
#line 318 "pathway_parser.y" /* yacc.c:1661  */
    {
        np_tt_append_string((yyvsp[-1]), (yyvsp[0])->strval);
        np_tt_destroy_value_struc((yyvsp[0]));
        (yyval) = (yyvsp[-1]);
    }
#line 2115 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 51:
#line 328 "pathway_parser.y" /* yacc.c:1661  */
    {
        int repeat = (yyvsp[0])->intval;
        np_tt_destroy_value_struc((yyvsp[0]));
        char * curstr = np_tt_strdup((yyvsp[-1])->strval);
        for (int i = 0; i < repeat - 1; i++) {
            np_tt_append_string((yyvsp[-1]), curstr);
        }
        free(curstr);
        (yyval) = (yyvsp[-1]);
    }
#line 2130 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 53:
#line 343 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_append_value((yyvsp[-1]), (yyvsp[0]));
    }
#line 2138 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 55:
#line 351 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_append_string((yyvsp[-1]), "*");
    }
#line 2146 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 56:
#line 358 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_definition(np_tt_make_int(TOK_TUBE),
                (yyvsp[-2]), (yyvsp[0]));
    }
#line 2155 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 58:
#line 367 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_append_value((yyvsp[-1]), (yyvsp[0]));
    }
#line 2163 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 62:
#line 380 "pathway_parser.y" /* yacc.c:1661  */
    { 
        np_tt_append_value((yyvsp[-4]), (yyvsp[-3]));
        np_tt_append_value((yyvsp[-4]), (yyvsp[-2]));
        (yyval) = np_tt_make_definition((yyvsp[-4]), (yyvsp[-5]), (yyvsp[0]));
    }
#line 2173 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 64:
#line 390 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = NULL; }
#line 2179 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 66:
#line 395 "pathway_parser.y" /* yacc.c:1661  */
    { 
        np_tt_append_value((yyvsp[-1]), (yyvsp[0])); 
        (yyval) = (yyvsp[-1]);
    }
#line 2188 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 67:
#line 403 "pathway_parser.y" /* yacc.c:1661  */
    {
        np_tt_append_value((yyvsp[-3]), (yyvsp[-1]));
        (yyval) = (yyvsp[-3]);
    }
#line 2197 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 68:
#line 408 "pathway_parser.y" /* yacc.c:1661  */
    {
        value_struc_t * tmp_value = np_tt_make_double(((yyvsp[0]))->doubleval);
        np_tt_append_value((yyvsp[0]), tmp_value);
        (yyval) = (yyvsp[0]);
    }
#line 2207 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 69:
#line 416 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_EXCLUDE); }
#line 2213 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 70:
#line 418 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_ACCESSIBLE); }
#line 2219 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 71:
#line 420 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_MATCHRANGE); }
#line 2225 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 72:
#line 424 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_WORDSIZE); }
#line 2231 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 75:
#line 435 "pathway_parser.y" /* yacc.c:1661  */
    {
        np_tt_append_value((yyvsp[-3]), (yyvsp[-2]));
        (yyval) = np_tt_make_definition((yyvsp[-3]), (yyvsp[-4]), (yyvsp[0]));
    }
#line 2240 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 76:
#line 443 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_definition((yyvsp[-2]), (yyvsp[-3]), (yyvsp[0]));
    }
#line 2248 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 77:
#line 450 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_LIBSEQ); }
#line 2254 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 78:
#line 455 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_SOURCE); }
#line 2260 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 79:
#line 460 "pathway_parser.y" /* yacc.c:1661  */
    {
        np_tt_append_value((yyvsp[-3]), (yyvsp[-2]));
        (yyval) = np_tt_make_definition((yyvsp[-3]), (yyvsp[-4]), (yyvsp[0]));
    }
#line 2269 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 80:
#line 468 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_CONCDEF); }
#line 2275 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 81:
#line 470 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_MAXSIZE); }
#line 2281 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 82:
#line 472 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_STOPDEF); }
#line 2287 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 83:
#line 474 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_WEIGHT); }
#line 2293 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 85:
#line 481 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_append_value((yyvsp[-2]), (yyvsp[0]));
    }
#line 2301 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 87:
#line 489 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_append_value((yyvsp[-1]), (yyvsp[0]));
    }
#line 2309 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 89:
#line 497 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_append_value((yyvsp[-2]), (yyvsp[0]));
    }
#line 2317 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 90:
#line 504 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_definition((yyvsp[-3]), (yyvsp[-2]), (yyvsp[0])); }
#line 2323 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 91:
#line 509 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_SYMMETRY_MIN); }
#line 2329 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 92:
#line 511 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_STRAND); }
#line 2335 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 93:
#line 513 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_WINDOW); }
#line 2341 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 94:
#line 515 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_COMPLEX); }
#line 2347 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 99:
#line 529 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_definition((yyvsp[-2]),
                NULL, (yyvsp[0]));
    }
#line 2356 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 100:
#line 537 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_MATERIAL); }
#line 2362 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 101:
#line 539 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_DANGLES); }
#line 2368 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 102:
#line 544 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_definition((yyvsp[-2]), NULL, (yyvsp[0])); }
#line 2374 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 103:
#line 546 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_definition((yyvsp[-3]), (yyvsp[-2]), (yyvsp[0])); }
#line 2380 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 104:
#line 551 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_PREVENT); }
#line 2386 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 105:
#line 556 "pathway_parser.y" /* yacc.c:1661  */
    {
        
        np_tt_append_value((yyvsp[-3]), (yyvsp[-2]));
        (yyval) = np_tt_make_definition((yyvsp[-3]),
                NULL, (yyvsp[0]));
    }
#line 2397 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 106:
#line 566 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_MAGNESIUM); }
#line 2403 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 107:
#line 568 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_GLOBAL_STOPDEF); }
#line 2409 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 108:
#line 570 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_SODIUM); }
#line 2415 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 109:
#line 572 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_TEMPERATURE); }
#line 2421 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 110:
#line 574 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_SEED); }
#line 2427 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 111:
#line 576 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_MBAD); }
#line 2433 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 112:
#line 578 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_MREOPT); }
#line 2439 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 113:
#line 580 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_MRESEED); }
#line 2445 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 114:
#line 582 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_FSPLIT); }
#line 2451 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 115:
#line 584 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_NSPLIT); }
#line 2457 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 116:
#line 586 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_HSPLIT); }
#line 2463 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 117:
#line 588 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_DGCLAMP); }
#line 2469 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 118:
#line 590 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_FSTRINGENT); }
#line 2475 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 119:
#line 592 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_FREDECOMP); }
#line 2481 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 120:
#line 594 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_FREFOCUS); }
#line 2487 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 121:
#line 596 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_FPASSIVE); }
#line 2493 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 122:
#line 598 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_GC_INIT); }
#line 2499 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 123:
#line 600 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_MINPAIR); }
#line 2505 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 124:
#line 602 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_TRIALS); }
#line 2511 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 125:
#line 604 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_OPTTIME); }
#line 2517 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 126:
#line 606 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_POPULATION); }
#line 2523 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 127:
#line 608 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_PRINTSTEPS); }
#line 2529 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 129:
#line 614 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_append_value((yyvsp[-2]), (yyvsp[0]));
    }
#line 2537 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 130:
#line 621 "pathway_parser.y" /* yacc.c:1661  */
    {
        (yyval) = np_tt_make_definition((yyvsp[-2]), NULL, (yyvsp[0]));
    }
#line 2545 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 131:
#line 627 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_PRINTLEAVES); }
#line 2551 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 132:
#line 629 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_REDECOMPOSE); }
#line 2557 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 133:
#line 631 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_ALLOWWOBBLE); }
#line 2563 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 134:
#line 633 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_ALLOWMISMATCH); }
#line 2569 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 135:
#line 635 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_DISABLEMUTWEIGHTS); }
#line 2575 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 136:
#line 640 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(1); }
#line 2581 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 137:
#line 642 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(0); }
#line 2587 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 138:
#line 647 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = NULL; }
#line 2593 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 139:
#line 649 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = (yyvsp[-1]) ; }
#line 2599 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 140:
#line 651 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_string("%"); }
#line 2605 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 141:
#line 657 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_definition(np_tt_make_int(TOK_SOURCE), (yyvsp[-2]) , (yyvsp[0])); }
#line 2611 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 142:
#line 662 "pathway_parser.y" /* yacc.c:1661  */
    { 
        np_tt_append_value((yyvsp[-4]), (yyvsp[-3]));
        (yyval) = np_tt_make_definition((yyvsp[-4]), (yyvsp[-2]), (yyvsp[0])); 
    }
#line 2620 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 143:
#line 670 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_COMPLEMENTARY); }
#line 2626 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 144:
#line 672 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_int(TOK_IDENTICAL); }
#line 2632 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 145:
#line 677 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_definition(np_tt_make_int(TOK_MATCH), (yyvsp[-2]), (yyvsp[0])); }
#line 2638 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 146:
#line 683 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_make_definition(np_tt_make_int(TOK_OFFTARGETS), (yyvsp[-4]), (yyvsp[0])); }
#line 2644 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 147:
#line 688 "pathway_parser.y" /* yacc.c:1661  */
    { 
      value_struc_t * inner = np_tt_make_list();
      value_struc_t * outer = np_tt_make_list();
      np_tt_list_append(inner, (yyvsp[0]));
      (yyval) = np_tt_list_append(outer, inner);
    }
#line 2655 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 148:
#line 695 "pathway_parser.y" /* yacc.c:1661  */
    { 
      value_struc_t * inner = np_tt_make_list();
      np_tt_list_append(inner, (yyvsp[0]));
      (yyval) = np_tt_list_append((yyvsp[-2]), inner);
    }
#line 2665 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 149:
#line 704 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = (yyvsp[-1]); }
#line 2671 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 151:
#line 710 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = np_tt_append_value((yyvsp[-2]), (yyvsp[0])); }
#line 2677 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 152:
#line 712 "pathway_parser.y" /* yacc.c:1661  */
    { 
      np_tt_append_value((yyvsp[-2]), np_tt_make_list());
      (yyval) = np_tt_append_value((yyvsp[-2]), (yyvsp[0]));
    }
#line 2686 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 153:
#line 717 "pathway_parser.y" /* yacc.c:1661  */
    { 
      np_tt_append_value((yyvsp[-4]), (yyvsp[-2]));
      (yyval) = np_tt_append_value((yyvsp[-4]), (yyvsp[0]));
    }
#line 2695 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 154:
#line 722 "pathway_parser.y" /* yacc.c:1661  */
    {
      (yyval) = np_tt_append_value(np_tt_make_int(0), (yyvsp[0]));
    }
#line 2703 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 155:
#line 726 "pathway_parser.y" /* yacc.c:1661  */
    {
      value_struc_t * tmp = np_tt_make_int(0);
      np_tt_append_value(tmp, (yyvsp[-2]));
      (yyval) = np_tt_append_value(tmp, (yyvsp[0]));
    }
#line 2713 "pathway_parser.c" /* yacc.c:1661  */
    break;

  case 156:
#line 735 "pathway_parser.y" /* yacc.c:1661  */
    { (yyval) = (yyvsp[-1]); }
#line 2719 "pathway_parser.c" /* yacc.c:1661  */
    break;


#line 2723 "pathway_parser.c" /* yacc.c:1661  */
        default: break;
      }
    if (yychar_backup != yychar)
      YY_LAC_DISCARD ("yychar change");
  }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;
  *++yylsp = yyloc;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (&yylloc, expression, scanner, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyesa, &yyes, &yyes_capacity, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        if (yychar != YYEMPTY)
          YY_LAC_ESTABLISH;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (&yylloc, expression, scanner, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }

  yyerror_range[1] = yylloc;

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, &yylloc, expression, scanner);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  yyerror_range[1] = yylsp[1-yylen];
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;

      yyerror_range[1] = *yylsp;
      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp, yylsp, expression, scanner);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  /* If the stack popping above didn't lose the initial context for the
     current lookahead token, the shift below will for sure.  */
  YY_LAC_DISCARD ("error recovery");

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  yyerror_range[2] = yylloc;
  /* Using YYLLOC is tempting, but would change the location of
     the lookahead.  YYLOC is available though.  */
  YYLLOC_DEFAULT (yyloc, yyerror_range, 2);
  *++yylsp = yyloc;

  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if 1
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (&yylloc, expression, scanner, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, &yylloc, expression, scanner);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, yylsp, expression, scanner);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  if (yyes != yyesa)
    YYSTACK_FREE (yyes);
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 738 "pathway_parser.y" /* yacc.c:1906  */


int yyerror(YYLTYPE * yyllocp, value_struc_t ** val, yyscan_t scanner, const char * msg) {
  if (yyllocp->first_line) {
    fprintf(stderr, "Line %d\n", LINE_NUMBER);
    fprintf(stderr, "%d.%d-%d.%d: error: ", 
      yyllocp->first_line, yyllocp->first_column,
      yyllocp->last_line, yyllocp->last_column);
  } else {
    fprintf(stderr, "Error: %s", msg);
  }
  fprintf(stderr, "%s", msg);
  fprintf(stderr, "\n");
  return 0;
} 

static value_struc_t * make_val_struc(void) {
  value_struc_t * val = malloc(sizeof(value_struc_t));
  val->type = NP_TT_INTEGER;
  val->intval = 0;
  val->doubleval = 0.0;
  val->list = NULL;
  val->strval = NULL;
  val->strlen = 0;
  val->strcap = 0;
  val->stype = NULL;
  val->id = NULL;
  val->def = NULL;
  val->next = NULL;
  return val;
}

char* np_tt_strdup(const char* str)
{
  size_t len;
  char* copy;

  len = strlen(str) + 1;
  copy = (char*) malloc(len);
  if (!copy)  return 0;
  memcpy(copy, str, len);
  return copy;
}

value_struc_t * np_tt_make_int(int input) {
  value_struc_t * val = make_val_struc();
  val->strval = NULL;
  val->type = NP_TT_INTEGER;
  val->intval = input;
  val->doubleval = (double) input;
  return val;
}

value_struc_t * np_tt_make_double(double input) {
  value_struc_t * val = make_val_struc();
  val->strval = NULL;
  val->type = NP_TT_FLOAT;
  val->doubleval = input;
  val->intval = (int) input;
  return val;
}

value_struc_t * np_tt_make_string(char * input) {
  value_struc_t * val = make_val_struc();
  int l = strlen(input);
  val->strval = (char*) malloc((l + 1) * sizeof(char));
  strncpy(val->strval, input, l + 1);
  val->type = NP_TT_STRING;
  val->strlen = l + 1;
  val->strcap = l + 1;
  return val;
}

value_struc_t * np_tt_append_string(value_struc_t * cur, 
    char * next) {
  int l = strlen(next);
  char * tempchar = cur->strval;
  char * nextdup = np_tt_strdup(next);
  if (l + cur->strlen + 1 > cur->strcap) {
    tempchar = (char*) realloc(cur->strval, 
          2 * (l + cur->strlen + 1));
    check_mem(tempchar);
    cur->strval = tempchar;
    cur->strcap = 2 * (l + cur->strlen + 1);
  }
  strncat(tempchar, nextdup, l+1);
  cur->strlen = l + cur->strlen + 1;
  free(nextdup);

  return cur;
error:
  cur->strval = NULL;
  return NULL;
}

value_struc_t * np_tt_make_list() {
  value_struc_t * val = make_val_struc();
  val->type = NP_TT_LIST;  
  return val;
}

value_struc_t * np_tt_make_definition(value_struc_t * stype, value_struc_t * id,
    value_struc_t * def) {
  value_struc_t * val = make_val_struc();
  val->type = NP_TT_DEFINITION;
  val->stype = stype;
  val->id = id;
  val->def = def;
  return val;
}

value_struc_t * np_tt_append_value(value_struc_t * cur, value_struc_t * app) {
  value_struc_t * cv = cur;
  while (cv->next) {cv = cv->next; }
  cv->next = app;
  return cur;
}

value_struc_t * np_tt_list_append(value_struc_t * cur, value_struc_t * app) {
  value_struc_t * tmp = cur->list;
  if (!tmp) {
    cur->list = app;
  } else {
    while (tmp->next) tmp = tmp->next;
    tmp->next = app;
  }
  return cur;  
}


void np_tt_destroy_value_struc(value_struc_t * val) {
  value_struc_t * cur = val;
  value_struc_t * prev = cur;

  while (cur) {
    np_tt_destroy_value_struc(cur->list);
    np_tt_destroy_value_struc(cur->stype);
    np_tt_destroy_value_struc(cur->id);
    np_tt_destroy_value_struc(cur->def);
    if (cur->strval && cur->strcap > 0) {
      free(cur->strval);
    }
    prev = cur;
    cur = cur->next;
    free(prev);
  }
}

void np_tt_print_ast(value_struc_t * root, int indent) {
  value_struc_t * cur = root;
  
  char * indent_string = (char *) malloc((indent + 1) * sizeof(char));
  for (int i = 0; i < indent; i++) {
    indent_string[i] = ' ';
  }
  indent_string[indent] = '\0';
  
  while (cur) {
    if (cur->type == NP_TT_DEFINITION) {
      printf("%sDEFINITION:\n", indent_string);
      printf("%s  stype:\n", indent_string);
      np_tt_print_ast(cur->stype, indent + 4);
      printf("%s  id:\n", indent_string);
      np_tt_print_ast(cur->id, indent + 4);
      printf("%s  def:\n", indent_string);
      np_tt_print_ast(cur->def, indent + 4);
    } else if (cur->type == NP_TT_STRING) {
      printf("%s%s\n", indent_string, cur->strval);
    } else if (cur->type == NP_TT_FLOAT) {
      printf("%s%lf\n", indent_string, cur->doubleval);
    } else if (cur->type == NP_TT_INTEGER) {
      printf("%s%i\n", indent_string, cur->intval);
    }
    cur = cur->next;
  }
  free(indent_string);
}


value_struc_t * get_ast(const char * parse_string) {
  value_struc_t * vals;
  yyscan_t scanner;
  YY_BUFFER_STATE state;

  if (yylex_init(&scanner)) {
    return NULL;
  }

  state = yy_scan_string(parse_string, scanner);

  vals = NULL;
  if (yyparse(&vals, scanner)) {
    return NULL;
  }

  yy_delete_buffer(state, scanner);

  yylex_destroy(scanner);

  return vals;
}
