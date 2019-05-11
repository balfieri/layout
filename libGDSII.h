/* Copyright (C) 2005-2017 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

/*
 * libGDSII is a C++ library for working with GDSII data files.
 * Homer Reid   11/2017
 */
#ifndef LIBGDSII_H
#define LIBGDSII_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include <string>
#include <vector>
#include <set>
#include <sstream>

using namespace std;

/***************************************************************/
/* convenient shorthand typedefs *******************************/
/***************************************************************/
#ifndef iVec
  typedef vector<int>    iVec;
#endif
#ifndef dVec
  typedef vector<double> dVec;
#endif
#ifndef bVec
  typedef vector<bool> bVec;
#endif
#ifndef sVec
  typedef vector<char *> sVec;
#endif
#ifndef strVec
  typedef vector<std::string> strVec;
#endif

/****************************************************************************************/
/* A PolygonList is a collection of polygons living in the XY plane.                    */
/* PolygonList.size() is the number of polygons in the list.                            */
/* PolygonList[np].size()/2 is the number of vertices in polygon #np.                   */
/* PolygonList[np][2*nv+0, 2*nv+1] are the x,y coordinates of vertex #nv in polygon #np.*/
/****************************************************************************************/
typedef vector<dVec> PolygonList;

typedef struct { char *Text; dVec XY; int Layer; } TextString;
typedef vector<TextString> TextStringList;

/***************************************************************/
/* Data structures used to process GDSII files.                */
/*  (a) GDSIIElement and GDSIIStruct are used to store info    */
/*      on the geometry as described in the GDSII file, with   */
/*      nesting hierarchy intact.                              */
/*  (b) Flattening the hierarchy (i.e. eliminating all SREFS   */
/*      and AREFS to instantiate all objects and text directly)*/
/*      yields a table of Entity structures, organized by the  */
/*      layer on which they appear. An Entity is simply just   */ 
/*      a polygon (collection of vertices, with an optional    */ 
/*      label) or a text string (with a single vertex as       */ 
/*      reference point/location). An EntityList is a          */ 
/*      collection of Entities, in no particular order, all on */ 
/*      the same layer. An EntityTable is a collection of      */ 
/*      (LayerIndex, EntityList) pairs.                        */
/*                                                             */
/*      Note that, whereas GDSIIElements and GDSIIStructs      */
/*      represent vertices as pairs of integers (multiples of  */
/*      the GDSII database unit), Entities represent vertices  */
/*      by pairs of doubles (real-valued, continuous physical  */
/*      coordinates). The default length unit for the          */
/*      coordinates of Entity vertices is 1 micron, but this   */
/*      may be changed by specifying a nonzero value for the   */
/*      CoordinateLengthUnit argument to Flatten(), or by      */
/*      setting the environment variable LIBGDSII_LENGTH_UNIT, */
/*      to the desired length unit in meters (default=1e-6).   */
/*      Thus, to output vertex coordinates in units of         */
/*      millimeters, set LengthUnit or LIBGDSII_LENGTH_UNIT to */
/*      1.0e-3.                                                */
/***************************************************************/
enum ElementType { BOUNDARY, PATH, SREF, AREF, TEXT, NODE, BOX };

typedef struct GDSIIElement
 { 
   ElementType Type;
   int Layer, DataType, TextType, PathType;
   iVec XY;
   std::string *SName;
   int Width, Columns, Rows;
   int nsRef;
   std::string *Text;
   bool Refl, AbsMag, AbsAngle;
   double Mag, Angle;
   iVec PropAttrs;
   strVec PropValues;
 } GDSIIElement;

typedef struct GDSIIStruct
 { 
   vector<GDSIIElement *> Elements;
   bool IsPCell;
   bool IsReferenced;
   std::string *Name;

 } GDSIIStruct;

typedef struct Entity
 { char *Text;   // if NULL, the entity is a polygon; otherwise it is a text string
   dVec XY;      // vertex coordinates: 2 for a text string, 2N for an N-gon
   bool Closed;  // true if there exists an edge connecting the last to the first vertex
   char *Label;  // optional descriptive text, may be present or absent for polygons and texts
 } Entity;

typedef vector<Entity>     EntityList;
typedef vector<EntityList> EntityTable;

/**********************************************************************/
/* GDSIIData is the main class that reads and stores a GDSII geometry.*/
/**********************************************************************/
namespace libGDSII
{
  /***************************************************************/
  /* GDSIIData describes the content of a single GDSII file. *****/
  /***************************************************************/
  class GDSIIData
   {
     /*--------------------------------------------------------*/
     /*- API methods                                           */
     /*--------------------------------------------------------*/
     public:
      
       // construct from a binary GDSII file 
       GDSIIData(const std::string FileName);
       ~GDSIIData();

       void WriteDescription(const char *FileName=0);

       // list of layer indices
       iVec GetLayers();

       // get all polygons on layer Layer that contain the reference point of
       // a GDSII text element matching Text (which must also lie on layer Layer).
       // If Layer==-1, search all layers.
       // If Text==NULL, return a list of all polygons on the given layer.
       PolygonList GetPolygons(const char *Text, int Layer=-1);
       PolygonList GetPolygons(int Layer=-1);
       TextStringList GetTextStrings(int Layer=-1);

     /*--------------------------------------------------------*/
     /* API data fields                                        */
     /*--------------------------------------------------------*/
      std::string *ErrMsg; // non-null upon failure of constructor or other API routine

     /*--------------------------------------------------------*/
     /* methods intended for internal use                      */
     /*--------------------------------------------------------*/
// private:
    // constructor helper methods
      void ReadGDSIIFile(const std::string FileName, double CoordinateLengthUnit=0.0);
      int GetStructByName(std::string Name);
      void Flatten(double CoordinateLengthUnit=0.0);

     /*--------------------------------------------------------*/
     /* variables intended for internal use                    */
     /*--------------------------------------------------------*/

     // general info on the GDSII file
     std::string *LibName;
     std::string *GDSIIFileName;
     double FileUnits[2], UnitInMeters;
     set<int> LayerSet; 
     iVec Layers;

     // list of structures (hierarchical, i.e. pre-flattening)
     vector<GDSIIStruct *> Structs;

     // table of entities (flattened)
     EntityTable ETable; // ETable[nl][ne] = #neth entity on layer Layers[nl]

     /*--------------------------------------------------------*/
     /*- utility routines -------------------------------------*/
     /*--------------------------------------------------------*/
     static bool Verbose;
     static char *LogFileName;
     static void Log(const char *format, ...);
     static void ErrExit(const char *format, ...);
     static void Warn(const char *format, ...);
     static char *vstrappend(char *s, const char *format, ...);
     static char *vstrdup(const char *format, ...);
   };

/***************************************************************/
/* non-class-method geometric primitives ***********************/
/***************************************************************/
bool PointInPolygon(dVec Vertices, double X, double Y);

/***********************************************************************/
/* the next few routines implement a caching mechanism by which an API */
/* code can make multiple calls to e.g. GetPolygons() for a given GDSII*/
/* file without requiring the API code to keep track of a GDSIIData    */
/* instance, but also without re-reading the file each time.           */
/* After the final such call the API code may call ClearGDSIICache()   */
/* to free memory allocated for the cache.                             */
/***********************************************************************/
iVec GetLayers(const char *GDSIIFile);
PolygonList GetPolygons(const char *GDSIIFile, const char *Text, int Layer=-1);
PolygonList GetPolygons(const char *GDSIIFile, int Layer=-1);
TextStringList GetTextStrings(const char *GDSIIFile, int Layer=-1);
void ClearGDSIICache();

/***************************************************************/
/* non-class method utility routines                           */
/***************************************************************/
bool DumpGDSIIFile(const char *FileName);
void WriteGMSHEntity(Entity E, int Layer, const char *geoFileName, FILE **pgeoFile,
                     const char *ppFileName=0, FILE **pppFile=0);
void WriteGMSHFile(EntityTable ETable, iVec Layers, char *FileBase, bool SeparateLayers=false);

} /* namespace libGDSII */

/***************************************************************/
/* crutch function to play nice with autotools *****************/
/***************************************************************/
extern "C" {
  void libGDSIIExists();
}

#endif /* LIBGDSII_H*/


/* Copyright (C) 2005-2017 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

/*
 * ReadGDSIIFile.cc -- GDSII file reader for libGDSII
 * Homer Reid   11/2017
 */

namespace libGDSII {

typedef unsigned char  BYTE;
typedef unsigned short WORD;
typedef unsigned long  DWORD;

/***************************************************************/
/* some data structures used in this file only *****************/
/***************************************************************/

/*--------------------------------------------------------------*/
/* storage for a single data record in the GDSII file           */
/*--------------------------------------------------------------*/
typedef struct GDSIIRecord
 {
   BYTE RType; // record type

   // could use a union for the following, but I don't bother
   bool Bits[16];
   iVec iVal;
   dVec dVal;
   string *sVal;
   size_t NumVals;

 } GDSIIRecord;

/*--------------------------------------------------------------*/
/*- 'ParseState' data structure maintained while reading .GDSII */
/*- file, updated after each record is read                     */
/*--------------------------------------------------------------*/
class GDSIIData; // forward reference 
typedef struct ParseState 
 { 
   GDSIIData *Data;
   int NumRecords;
   enum { INITIAL,
          INHEADER,  INLIB,  INSTRUCT, INELEMENT,
          DONE
        } Status;
   GDSIIStruct *CurrentStruct;
   GDSIIElement *CurrentElement;

 } ParseState;

typedef string *(*RecordHandler)(GDSIIRecord Record, ParseState *PState);

const char *ElTypeNames[]=
 {"BOUNDARY", "PATH", "SREF", "AREF", "TEXT", "NODE", "BOX"};

/***************************************************************/
/* Handlers for specific types of data records in GDSII files. */
/***************************************************************/
string *handleHEADER(GDSIIRecord Record, ParseState *PState)
{
  (void) Record;
  if (PState->Status!=ParseState::INITIAL)
   return new string("unexpected record before HEADER");
  PState->Status=ParseState::INHEADER;
  return 0;
}

string *handleBGNLIB(GDSIIRecord Record, ParseState *PState)
{
  (void) Record;
  if (PState->Status!=ParseState::INHEADER)
   return new string("unexpected record BGNLIB");
  PState->Status=ParseState::INLIB;
  return 0;
}

string *handleLIBNAME(GDSIIRecord Record, ParseState *PState)
{ 
  if (PState->Status!=ParseState::INLIB)
   return new string("unexpected record LIBNAME");
  PState->Data->LibName = new string( *(Record.sVal) );
  return 0;
}

string *handleUNITS(GDSIIRecord Record, ParseState *PState)
{ 
  PState->Data->FileUnits[0] = Record.dVal[0];
  PState->Data->FileUnits[1] = Record.dVal[1];
  PState->Data->UnitInMeters =
   PState->Data->FileUnits[1] / PState->Data->FileUnits[0];
  return 0;
}

string *handleENDLIB(GDSIIRecord Record, ParseState *PState)
{ 
  (void) Record;
  if (PState->Status!=ParseState::INLIB)
   return new string("unexpected record ENDLIB");
  PState->Status=ParseState::DONE;
  return 0;
}

string *handleBGNSTR(GDSIIRecord Record, ParseState *PState)
{
  (void) Record;
  if (PState->Status!=ParseState::INLIB)
   return new string("unexpected record BGNSTR");

  // add a new structure
  GDSIIStruct *s  = new GDSIIStruct;
  s->IsReferenced = false;
  s->IsPCell      = false;
  PState->CurrentStruct = s;
  PState->Data->Structs.push_back(s);

  PState->Status=ParseState::INSTRUCT;

  return 0;
}

string *handleSTRNAME(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INSTRUCT)
   return new string("unexpected record STRNAME");
  PState->CurrentStruct->Name = new string( *(Record.sVal) );
  if( strcasestr( Record.sVal->c_str(), "CONTEXT_INFO") )
   PState->CurrentStruct->IsPCell=true;
  return 0;
}

string *handleENDSTR(GDSIIRecord Record, ParseState *PState)
{
  (void) Record;
  if (PState->Status!=ParseState::INSTRUCT)
   return new string("unexpected record ENDSTR");
  PState->Status=ParseState::INLIB;
  return 0;
}

string *handleElement(GDSIIRecord Record, ParseState *PState, ElementType ElType)
{
  (void) Record;
  if (PState->Status!=ParseState::INSTRUCT)
   return new string(   string("unexpected record") + ElTypeNames[ElType] );
  
  // add a new element
  GDSIIElement *e = new GDSIIElement;
  e->Type     = ElType;
  e->Layer    = 0;
  e->DataType = 0;
  e->TextType = 0;
  e->PathType = 0;
  e->SName    = 0;
  e->Width    = 0;
  e->Columns  = 0;
  e->Rows     = 0;
  e->Text     = 0;
  e->Refl     = false;
  e->AbsMag   = false;
  e->AbsAngle = false;
  e->Mag      = 1.0;
  e->Angle    = 0.0;
  e->nsRef    = -1;
  PState->CurrentElement = e;
  PState->CurrentStruct->Elements.push_back(e);

  PState->Status=ParseState::INELEMENT;
  return 0;
}

string *handleBOUNDARY(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, BOUNDARY); }

string *handlePATH(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, PATH); }

string *handleSREF(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, SREF); }

string *handleAREF(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, AREF); }

string *handleTEXT(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, TEXT); }

string *handleNODE(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, NODE); }

string *handleBOX(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, BOX); }

string *handleLAYER(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new string("unexpected record LAYER");
  PState->CurrentElement->Layer = Record.iVal[0];
  PState->Data->LayerSet.insert(Record.iVal[0]);
  
  return 0;
}

string *handleDATATYPE(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new string("unexpected record DATATYPE");
  PState->CurrentElement->DataType = Record.iVal[0];
  return 0;
}

string *handleTEXTTYPE(GDSIIRecord Record, ParseState *PState)
{ 
  if (    PState->Status!=ParseState::INELEMENT
       || PState->CurrentElement->Type!=TEXT
     )
   return new string("unexpected record TEXTTYPE");
  PState->CurrentElement->TextType = Record.iVal[0];
  return 0;
}

string *handlePATHTYPE(GDSIIRecord Record, ParseState *PState)
{ 
  if (    PState->Status!=ParseState::INELEMENT )
   return new string("unexpected record PATHTYPE");
  PState->CurrentElement->PathType = Record.iVal[0];
  return 0;
}

string *handleSTRANS(GDSIIRecord Record, ParseState *PState)
{ 
  if ( PState->Status!=ParseState::INELEMENT )
   return new string("unexpected record STRANS");
  PState->CurrentElement->Refl     = Record.Bits[0];
  PState->CurrentElement->AbsMag   = Record.Bits[13];
  PState->CurrentElement->AbsAngle = Record.Bits[14];
  return 0;
}

string *handleMAG(GDSIIRecord Record, ParseState *PState)
{ 
  if ( PState->Status!=ParseState::INELEMENT )
   return new string("unexpected record MAG");
  PState->CurrentElement->Mag = Record.dVal[0];
  return 0;
}

string *handleANGLE(GDSIIRecord Record, ParseState *PState)
{ 
  if ( PState->Status!=ParseState::INELEMENT )
   return new string("unexpected record ANGLE");
  PState->CurrentElement->Angle = Record.dVal[0];
  return 0;
}

string *handlePROPATTR(GDSIIRecord Record, ParseState *PState)
{ 
  if ( PState->Status!=ParseState::INELEMENT )
   return new string("unexpected record PROPATTR");
  GDSIIElement *e=PState->CurrentElement;
  e->PropAttrs.push_back(Record.iVal[0]);
  e->PropValues.push_back("");
  return 0;
}

string *handlePROPVALUE(GDSIIRecord Record, ParseState *PState)
{
  if ( PState->Status!=ParseState::INELEMENT )
   return new string("unexpected record PROPVALUE");
  GDSIIElement *e=PState->CurrentElement;
  int n=e->PropAttrs.size();
  if (n==0)
   return new string("PROPVALUE without PROPATTR");
  e->PropValues[n-1]=string( *(Record.sVal) );

  if( strcasestr( Record.sVal->c_str(), "CONTEXT_INFO") )
   PState->CurrentStruct->IsPCell=true;

  return 0;
}

string *handleXY(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new string("unexpected record XY");
  PState->CurrentElement->XY.reserve(Record.NumVals);
  for(size_t n=0; n<Record.NumVals; n++)
   PState->CurrentElement->XY.push_back(Record.iVal[n]);
  return 0;
}

string *handleSNAME(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new string("unexpected record SNAME");
  PState->CurrentElement->SName = new string( *(Record.sVal) );
  return 0;
}

string *handleSTRING(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new string("unexpected record STRING");
  PState->CurrentElement->Text = new string( *(Record.sVal) );
  return 0;
}

string *handleCOLROW(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new string("unexpected record COLROW");
  PState->CurrentElement->Columns = Record.iVal[0];
  PState->CurrentElement->Rows    = Record.iVal[1];
  return 0;
}

string *handleWIDTH(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new string("unexpected record Width");
  PState->CurrentElement->Width   = Record.iVal[0];
  return 0;
}

string *handleENDEL(GDSIIRecord Record, ParseState *PState)
{
  (void) Record;
  if (PState->Status!=ParseState::INELEMENT)
   return new string("unexpected record ENDEL");
  PState->Status = ParseState::INSTRUCT;
  return 0;
}

/***************************************************************/
/* table of GDSII data types ***********************************/
/***************************************************************/
enum DataType{ NO_DATA,    // 0x00
               BITARRAY,   // 0x01
               INTEGER_2,  // 0x02
               INTEGER_4,  // 0x03
               REAL_4,     // 0x04
               REAL_8,     // 0x05
               STRING      // 0x06
              };

/***************************************************************/
/* table of GDS record types, gleeped directly from the text of*/
/* the buchanan email                                          */
/***************************************************************/
typedef struct RecordType
 { const char    *Name;
   DataType       DType;
   RecordHandler  Handler;
 } RecordType;

const static RecordType RecordTypes[]={
 /*0x00*/  {"HEADER",       INTEGER_2,   handleHEADER},
 /*0x01*/  {"BGNLIB",       INTEGER_2,   handleBGNLIB},
 /*0x02*/  {"LIBNAME",      STRING,      handleLIBNAME},
 /*0x03*/  {"UNITS",        REAL_8,      handleUNITS},
 /*0x04*/  {"ENDLIB",       NO_DATA,     handleENDLIB},
 /*0x05*/  {"BGNSTR",       INTEGER_2,   handleBGNSTR},
 /*0x06*/  {"STRNAME",      STRING,      handleSTRNAME},
 /*0x07*/  {"ENDSTR",       NO_DATA,     handleENDSTR},
 /*0x08*/  {"BOUNDARY",     NO_DATA,     handleBOUNDARY},
 /*0x09*/  {"PATH",         NO_DATA,     handlePATH},
 /*0x0a*/  {"SREF",         NO_DATA,     handleSREF},
 /*0x0b*/  {"AREF",         NO_DATA,     handleAREF},
 /*0x0c*/  {"TEXT",         NO_DATA,     handleTEXT},
 /*0x0d*/  {"LAYER",        INTEGER_2,   handleLAYER},
 /*0x0e*/  {"DATATYPE",     INTEGER_2,   handleDATATYPE},
 /*0x0f*/  {"WIDTH",        INTEGER_4,   handleWIDTH},
 /*0x10*/  {"XY",           INTEGER_4,   handleXY},
 /*0x11*/  {"ENDEL",        NO_DATA,     handleENDEL},
 /*0x12*/  {"SNAME",        STRING,      handleSNAME},
 /*0x13*/  {"COLROW",       INTEGER_2,   handleCOLROW},
 /*0x14*/  {"TEXTNODE",     NO_DATA,     0},
 /*0x15*/  {"NODE",         NO_DATA,     0},
 /*0x16*/  {"TEXTTYPE",     INTEGER_2,   handleTEXTTYPE},
 /*0x17*/  {"PRESENTATION", BITARRAY,    0},
 /*0x18*/  {"UNUSED",       NO_DATA,     0},
 /*0x19*/  {"STRING",       STRING,      handleSTRING},
 /*0x1a*/  {"STRANS",       BITARRAY,    handleSTRANS},
 /*0x1b*/  {"MAG",          REAL_8,      handleMAG},
 /*0x1c*/  {"ANGLE",        REAL_8,      handleANGLE},
 /*0x1d*/  {"UNUSED",       NO_DATA,     0},
 /*0x1e*/  {"UNUSED",       NO_DATA,     0},
 /*0x1f*/  {"REFLIBS",      STRING,      0},
 /*0x20*/  {"FONTS",        STRING,      0},
 /*0x21*/  {"PATHTYPE",     INTEGER_2,   handlePATHTYPE},
 /*0x22*/  {"GENERATIONS",  INTEGER_2,   0},
 /*0x23*/  {"ATTRTABLE",    STRING,      0},
 /*0x24*/  {"STYPTABLE",    STRING,      0},
 /*0x25*/  {"STRTYPE",      INTEGER_2,   0},
 /*0x26*/  {"ELFLAGS",      BITARRAY,    0},
 /*0x27*/  {"ELKEY",        INTEGER_4,   0},
 /*0x1d*/  {"LINKTYPE",     NO_DATA,     0},
 /*0x1e*/  {"LINKKEYS",     NO_DATA,     0},
 /*0x2a*/  {"NODETYPE",     INTEGER_2,   0},
 /*0x2b*/  {"PROPATTR",     INTEGER_2,   handlePROPATTR},
 /*0x2c*/  {"PROPVALUE",    STRING,      handlePROPVALUE},
 /*0x2d*/  {"BOX",          NO_DATA,     0},
 /*0x2e*/  {"BOXTYPE",      INTEGER_2,   0},
 /*0x2f*/  {"PLEX",         INTEGER_4,   0},
 /*0x30*/  {"BGNEXTN",      INTEGER_4,   0},
 /*0x31*/  {"ENDTEXTN",     INTEGER_4,   0},
 /*0x32*/  {"TAPENUM",      INTEGER_2,   0},
 /*0x33*/  {"TAPECODE",     INTEGER_2,   0},
 /*0x34*/  {"STRCLASS",     BITARRAY,    0},
 /*0x35*/  {"RESERVED",     INTEGER_4,   0},
 /*0x36*/  {"FORMAT",       INTEGER_2,   0},
 /*0x37*/  {"MASK",         STRING,      0},
 /*0x38*/  {"ENDMASKS",     NO_DATA,     0},
 /*0x39*/  {"LIBDIRSIZE",   INTEGER_2,   0},
 /*0x3a*/  {"SRFNAME",      STRING,      0},
 /*0x3b*/  {"LIBSECUR",     INTEGER_2,   0}
};

#define RTYPE_HEADER		0x00
#define RTYPE_BGNLIB		0x01
#define RTYPE_LIBNAME		0x02
#define RTYPE_UNITS		0x03
#define RTYPE_ENDLIB		0x04
#define RTYPE_BGNSTR		0x05
#define RTYPE_STRNAME		0x06
#define RTYPE_ENDSTR		0x07
#define RTYPE_BOUNDARY		0x08
#define RTYPE_PATH		0x09
#define RTYPE_SREF		0x0a
#define RTYPE_AREF		0x0b
#define RTYPE_TEXT		0x0c
#define RTYPE_LAYER		0x0d
#define RTYPE_DATATYPE		0x0e
#define RTYPE_WIDTH		0x0f
#define RTYPE_XY		0x10
#define RTYPE_ENDEL		0x11
#define RTYPE_SNAME		0x12
#define RTYPE_COLROW		0x13
#define RTYPE_TEXTNODE		0x14
#define RTYPE_NODE		0x15
#define RTYPE_TEXTTYPE		0x16
#define RTYPE_PRESENTATION	0x17
#define RTYPE_UNUSED		0x18
#define RTYPE_STRING		0x19
#define RTYPE_STRANS		0x1a
#define RTYPE_MAG		0x1b
#define RTYPE_ANGLE		0x1c
#define RTYPE_UNUSED2		0x1d
#define RTYPE_UNUSED3		0x1e
#define RTYPE_REFLIBS		0x1f
#define RTYPE_FONTS		0x20
#define RTYPE_PATHTYPE		0x21
#define RTYPE_GENERATIONS	0x22
#define RTYPE_ATTRTABLE		0x23
#define RTYPE_STYPTABLE		0x24
#define RTYPE_STRTYPE		0x25
#define RTYPE_ELFLAGS		0x26
#define RTYPE_ELKEY		0x27
#define RTYPE_LINKTYPE		0x1d
#define RTYPE_LINKKEYS		0x1e
#define RTYPE_NODETYPE		0x2a
#define RTYPE_PROPATTR		0x2b
#define RTYPE_PROPVALUE		0x2c
#define RTYPE_BOX		0x2d
#define RTYPE_BOXTYPE		0x2e
#define RTYPE_PLEX		0x2f
#define RTYPE_BGNEXTN		0x30
#define RTYPE_ENDTEXTN		0x31
#define RTYPE_TAPENUM		0x32
#define RTYPE_TAPECODE		0x33
#define RTYPE_STRCLASS		0x34
#define RTYPE_RESERVED		0x35
#define RTYPE_FORMAT		0x36
#define RTYPE_MASK		0x37
#define RTYPE_ENDMASKS		0x38
#define RTYPE_LIBDIRSIZE	0x39
#define RTYPE_SRFNAME		0x3a
#define RTYPE_LIBSECUR		0x3b
#define MAX_RTYPE     		0x3b

/***************************************************************/
/***************************************************************/
/***************************************************************/
int ConvertInt(BYTE *Bytes, DataType DType)
{ 
  unsigned long long i = Bytes[0]*256 + Bytes[1];
  if (DType==INTEGER_4)
   i = i*256*256 + Bytes[2]*256 + Bytes[3];
  if (Bytes[0] & 0x80) // sign bit
   return -1*( (DType==INTEGER_2 ? 0x010000 : 0x100000000) - i );
  return i;
}

double ConvertReal(BYTE *Bytes, DataType DType)
{ 
  double Sign  = (Bytes[0] & 0x80) ? -1.0 : +1.0;
  int Exponent = (Bytes[0] & 0x7F) - 64;
  int NumMantissaBytes = (DType==REAL_4 ? 3 : 7);
  int NumMantissaBits  = 8*NumMantissaBytes;
  double Mantissa=0.0;
  for(int n=0; n<NumMantissaBytes; n++)
   Mantissa = Mantissa*256 + ((double)(Bytes[1+n]));
  return Sign * Mantissa * pow(2.0, 4*Exponent - NumMantissaBits);
}


// The allowed characters are all ASCII-printable characters, including space, except comma (,) and double quote (").
// Non-allowed characters at the end of the string are removed.
// Non-allowed characters not at the end of the string are converted to underscores.
bool IsAllowedChar(char c)
{ return isprint(c) && c!='"' && c!=','; }

string *MakeGDSIIString(char *Original, int Size)
{ 
  if (Size==0) return new string("");

  if (Size>32) Size=32;
  char RawString[33];
  strncpy(RawString, Original, Size);
  RawString[Size]=0;
  int L = strlen(RawString);
  while ( L>0 && !IsAllowedChar(RawString[L-1]) )
   RawString[--L] = 0;
  for(int n=0; n<L; n++) 
   if (!IsAllowedChar(RawString[n])) RawString[n]='_';
  return new string(RawString);
}

/***************************************************************/
/* read a single GDSII data record from the current file position */
/***************************************************************/
GDSIIRecord ReadGDSIIRecord(FILE *f, string **ErrMsg)
{
  /*--------------------------------------------------------------*/
  /* read the 4-byte file header and check that the data type     */
  /* agrees with what it should be based on the record type       */
  /*--------------------------------------------------------------*/
  BYTE Header[4];
  if ( 4 != fread(Header, 1, 4, f) )
   { *ErrMsg = new string("unexpected end of file");
     return GDSIIRecord(); // end of file
   }

  size_t RecordSize = Header[0]*256 + Header[1];
  BYTE RType        = Header[2];
  BYTE DType        = Header[3];
  
  if (RType > MAX_RTYPE)
   { *ErrMsg = new string("unknown record type");
     return GDSIIRecord();
   }
    
  if ( DType != RecordTypes[RType].DType )
   { ostringstream ss;
     ss << RecordTypes[RType].Name
        << ": data type disagrees with record type ("
        << DType
        << " != "
        << RecordTypes[RType].DType
        << ")";
     *ErrMsg = new string(ss.str());
     return GDSIIRecord();
   }

  /*--------------------------------------------------------------*/
  /*- attempt to read payload ------------------------------------*/
  /*--------------------------------------------------------------*/
  size_t PayloadSize = RecordSize - 4;
  BYTE *Payload=0;
  if (PayloadSize>0)
   { Payload = new BYTE[PayloadSize];
     if (Payload==0)
      { *ErrMsg = new string("out of memory");
        return GDSIIRecord();
      }
     if ( PayloadSize != fread((void *)Payload, 1, PayloadSize, f) )
      { delete[] Payload;
        *ErrMsg = new string("unexpected end of file");
        return GDSIIRecord();
      }
   }
 
  /*--------------------------------------------------------------*/
  /* allocate space for the record and process payload data       */
  /*--------------------------------------------------------------*/
  GDSIIRecord Record;
  Record.RType   = RType;
  Record.NumVals = 0;
  Record.sVal    = 0;

  switch(DType)
   { case NO_DATA:
       break;

     case BITARRAY:
      { Record.NumVals=1;
        WORD W = *(WORD *)Payload;
        for(unsigned nf=0, Flag=1; nf<16; nf++, Flag*=2)
         Record.Bits[nf] = (W & Flag);
      };
     break;

     case STRING:
      Record.NumVals=1;
      Record.sVal = MakeGDSIIString( (char *)Payload, PayloadSize );
      break;

     case INTEGER_2:
     case INTEGER_4:
      { size_t DataSize = (DType==INTEGER_2) ? 2 : 4;
        Record.NumVals  = PayloadSize / DataSize;
        BYTE *B=(BYTE *)Payload; 
        for(size_t nv=0; nv<Record.NumVals; nv++, B+=DataSize)
         Record.iVal.push_back( ConvertInt(B, RecordTypes[RType].DType) );
      };
     break;

     case REAL_4:
     case REAL_8:
      { size_t DataSize  = (DType==REAL_4) ? 4 : 8;
        Record.NumVals   = PayloadSize / DataSize;
        BYTE *B=(BYTE *)Payload; 
        for(size_t nv=0; nv<Record.NumVals; nv++, B+=DataSize)
         Record.dVal.push_back(ConvertReal(B, RecordTypes[RType].DType));
      };
     break;

     default:
       *ErrMsg = new string("unknown data type " + std::to_string(DType));
       return GDSIIRecord();
   };

  // success 
  *ErrMsg=0;
  delete[] Payload;
  return Record;

}

/***************************************************************/
/* get string description of GDSII record   ********************/
/***************************************************************/
string *GetRecordDescription(GDSIIRecord Record, bool Verbose=true)
{
  char Name[15];
  sprintf(Name,"%12s",RecordTypes[Record.RType].Name);
  ostringstream ss;
  ss << Name;

  if (Record.NumVals>0)
   ss << " ( " << Record.NumVals << ") ";
  
  if (!Verbose)
   return new string(ss.str());

  ss << " = ";
  switch(RecordTypes[Record.RType].DType)
   { 
     case INTEGER_2:    
     case INTEGER_4:    
      for(size_t nv=0; nv<Record.NumVals; nv++)
       ss << Record.iVal[nv] << " ";
      break;

     case REAL_4:
     case REAL_8:
      for(size_t nv=0; nv<Record.NumVals; nv++)
       ss << Record.dVal[nv] << " ";
      break;

     case BITARRAY:
      for(size_t n=0; n<16; n++)
       ss << (Record.Bits[n] ? '1' : '0');
      break;

     case STRING:
      if (Record.sVal)
       ss << *(Record.sVal);
      else
       ss << "(null)";
      break; 

     case NO_DATA:
     default:
      break; 
   };

  return new string(ss.str());

}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void InitializeParseState(ParseState *PState, GDSIIData *Data)
{
  PState->Data           = Data;
  PState->NumRecords     = 0;
  PState->CurrentStruct  = 0;
  PState->CurrentElement = 0;
  PState->Status         = ParseState::INITIAL;
}

/*--------------------------------------------------------------*/
/*- If CoordinateLengthUnit is nonzero, it sets the desired     */
/*- output unit (in meters) for vertex coordinates.             */
/*--------------------------------------------------------------*/
void GDSIIData::ReadGDSIIFile(const string FileName, double CoordinateLengthUnit)
 {
   ErrMsg=0;

   /*--------------------------------------------------------------*/
   /*- try to open the file ---------------------------------------*/
   /*--------------------------------------------------------------*/
   FILE *f=fopen(FileName.c_str(),"r");
   if (!f)
    { ErrMsg = new string("could not open " + FileName);
      return;
    }

   /*--------------------------------------------------------------*/
   /*- read records one at a time until we hit ENDLIB              */
   /*--------------------------------------------------------------*/
   ParseState PState;
   InitializeParseState(&PState, this);
   while( PState.Status != ParseState::DONE && !ErrMsg )
    { 
      // try to read the record
      GDSIIRecord Record=ReadGDSIIRecord(f, &ErrMsg);
      if (ErrMsg)
       return;

      // try to process the record if a handler is present
      PState.NumRecords++;
      RecordHandler Handler = RecordTypes[Record.RType].Handler;
      if ( Handler )
       ErrMsg = Handler(Record, &PState);
      else 
       Warn("ignoring unsupported record %s",RecordTypes[Record.RType].Name);
    }
   fclose(f);
   if (ErrMsg) return;
 
   // convert layer set to vector
   for(set<int>::iterator it=LayerSet.begin(); it!=LayerSet.end(); it++)
    Layers.push_back(*it);

   /*--------------------------------------------------------------*/
   /*- Go back through the hierarchy to note which structures are  */
   /*- referenced by others.                                       */
   /*--------------------------------------------------------------*/
   for(size_t ns=0; ns<Structs.size(); ns++)
    for(size_t ne=0; ne<Structs[ns]->Elements.size(); ne++)
     { GDSIIElement *e=Structs[ns]->Elements[ne];
       if(e->Type==SREF || e->Type==AREF)
        { e->nsRef = GetStructByName( *(e->SName) );
          if (e->nsRef==-1)
           Warn("reference to unknown struct %s ",e->SName->c_str());
          else 
           Structs[e->nsRef]->IsReferenced=true;
        }
     }

   /*--------------------------------------------------------------*/
   /*- Flatten hierarchy to obtain simple unstructured lists       */
   /*- of polygons and text labels on each layer.                  */
   /*--------------------------------------------------------------*/
   Flatten(CoordinateLengthUnit);
}

/***************************************************************/
/* Write text description of GDSII file to FileName.           */
/***************************************************************/
void GDSIIData::WriteDescription(const char *FileName)
{
  FILE *f = (FileName == 0 ? stdout : fopen(FileName,"w") );

  fprintf(f,"*\n");
  fprintf(f,"* File %s: \n",GDSIIFileName->c_str());
  if (LibName)
   fprintf(f,"* Library %s: \n",LibName->c_str());
  fprintf(f,"* Unit=%e meters (file units = {%e,%e})\n",UnitInMeters,FileUnits[0],FileUnits[1]);
  fprintf(f,"*\n");

  fprintf(f,"**************************************************\n");

  fprintf(f,"** Library %s:\n",LibName->c_str());
  fprintf(f,"**************************************************\n");
  for(size_t ns=0; ns<Structs.size(); ns++)
   { 
     GDSIIStruct *s=Structs[ns];
     fprintf(f,"--------------------------------------------------\n");
     fprintf(f,"** Struct %i: %s\n",(int )ns,s->Name->c_str());
     fprintf(f,"--------------------------------------------------\n");

    for(size_t ne=0; ne<s->Elements.size(); ne++)
     { GDSIIElement *e=s->Elements[ne];
       fprintf(f,"  Element %i: %s (layer %i, datatype %i)\n",
                    (int )ne, ElTypeNames[e->Type], e->Layer, e->DataType);
       if (e->Type==PATH || e->Type==TEXT)
        fprintf(f,"    (width %i, pathtype %i)\n",e->Width, e->PathType);
       if (e->Text)
        fprintf(f,"    (text %s)\n",e->Text->c_str());
       if (e->SName)
        fprintf(f,"    (structure %s)\n",e->SName->c_str());
       if (e->Mag!=1.0 || e->Angle!=0.0)
        fprintf(f,"    (mag %g, angle %g)\n",e->Mag,e->Angle);
       if (e->Columns!=0 || e->Rows!=0)          
        fprintf(f,"    (%i x %i array)\n",e->Columns,e->Rows);
       for(size_t n=0; n<e->PropAttrs.size(); n++)
        fprintf(f,"    (attribute %i: %s)\n",e->PropAttrs[n],e->PropValues[n].c_str());
       fprintf(f,"     XY: ");

       for(size_t nxy=0; nxy<e->XY.size(); nxy++)
        fprintf(f,"%i ",e->XY[nxy]);
       fprintf(f,"\n\n");
     }
   }
  if (FileName)
   fclose(f);
}

/***************************************************************/
/* non-class utility method to print a raw dump of all data    */
/* records in a GDSII file                                     */
/***************************************************************/
bool DumpGDSIIFile(const char *GDSIIFileName)
{
  FILE *f=fopen(GDSIIFileName,"r");
  if (!f)
   { fprintf(stderr,"error: could not open %s (aborting)\n",GDSIIFileName);
     return false;
   };

  /*--------------------------------------------------------------*/
  /*- read records one at a time ---------------------------------*/
  /*--------------------------------------------------------------*/
  int NumRecords=0;
  bool Done=false;
  while(!Done)
   { 
     string *ErrMsg=0;
     GDSIIRecord Record=ReadGDSIIRecord(f, &ErrMsg);
     if (ErrMsg)
      { fprintf(stderr,"error: %s (aborting)\n",ErrMsg->c_str());
        return false;
      }

     string *RStr = GetRecordDescription(Record);
     printf("Record %i: %s\n",NumRecords++,RStr->c_str());
     delete RStr;

     if (Record.RType==RTYPE_ENDLIB)
      Done=true;
   }
  fclose(f);

  printf("Read %i data records from file %s.\n",NumRecords,GDSIIFileName);
  return true;
}

} // namespace libGSDII


/* Copyright (C) 2005-2017 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

/*
 * Flatten.cc  -- flatten hierarchical GDSII data set
 * Homer Reid   12/2017-4/2018
 */

namespace libGDSII{

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct GTransform
 { double X0, Y0;
   double CosTheta, SinTheta;
   double Mag;
   bool Refl;
 } GTransform;

typedef vector<GTransform> GTVec;

static void ApplyGTransform(GTransform GT, double X, double Y, double *XP, double *YP)
{
   X *= GT.Mag;
   Y *= (GT.Refl ? -1.0 : 1.0)*GT.Mag;
  double NewX = GT.X0 + GT.CosTheta*X - GT.SinTheta*Y;
  double NewY = GT.Y0 + GT.SinTheta*X + GT.CosTheta*Y;
  *XP = NewX;
  *YP = NewY;
}

/***************************************************************/
/* data structure that keeps track of the status of the GMSH   */
/* file output process                                         */
/***************************************************************/
typedef struct StatusData
{ int CurrentLayer;
  double IJ2XY; // scale factor converting GDSII integer-value vertex indices to real-valued coordinates in the chosen length units
  EntityList EntitiesThisLayer;
  GTVec GTStack;
  int RefDepth;
} StatusData;

static void InitStatusData(StatusData *SD, double CoordinateLengthUnit, double PixelLengthUnit)
{ SD->CurrentLayer=-1;
  SD->IJ2XY = PixelLengthUnit / CoordinateLengthUnit;
  SD->RefDepth=0;
}

//FIXME 
static void GetPhysicalXY(StatusData *SD, double X, double Y, double *pXP, double *pYP)
{
  double XP=X, YP=Y;
  for(unsigned n=SD->GTStack.size(); n>0; n--)
   { ApplyGTransform(SD->GTStack[n-1],X,Y,&XP,&YP);
     X=XP;
     Y=YP;
   }
  *pXP = SD->IJ2XY * XP;
  *pYP = SD->IJ2XY * YP;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddBoundary(StatusData *SD, GDSIIData *Data, int ns, int ne)
{
  GDSIIStruct *s  = Data->Structs[ns];
  GDSIIElement *e = s->Elements[ne];
  if (SD->CurrentLayer!=e->Layer) return;

  vector<int> IXY = e->XY;
  int NXY         = IXY.size() / 2;

  char Label[1000];
  snprintf(Label,1000,"Struct %s element #%i (boundary)",s->Name->c_str(),ne);

  Entity E;
  E.XY.resize(IXY.size() - 2);
  E.Text   = 0;
  E.Label  = strdup(Label);
  E.Closed = true;
  for(int n=0; n<NXY-1; n++)
   GetPhysicalXY(SD, IXY[2*n+0], IXY[2*n+1], &(E.XY[2*n]), &(E.XY[2*n+1]));

  SD->EntitiesThisLayer.push_back(E);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddPath(StatusData *SD, GDSIIData *Data, int ns, int ne)
{
  GDSIIStruct *s  = Data->Structs[ns];
  GDSIIElement *e = s->Elements[ne];

  if (SD->CurrentLayer!=e->Layer) return;
  char Label[1000];
  snprintf(Label,1000,"Struct %s element #%i (path)",s->Name->c_str(),ne);

  vector<int> IXY = e->XY;
  int NXY         = IXY.size() / 2;

  double IJ2XY    = SD->IJ2XY;
  double W        = e->Width*IJ2XY;

  Entity E;
  E.Text   = 0;
  E.Label  = strdup(Label);
  E.Closed = (W!=0.0);
  int NumNodes = (W==0.0 ? NXY : 2*NXY);
  E.XY.resize(2*NumNodes); 

  if (W==0.0)
   for(int n=0; n<NXY; n++)
    GetPhysicalXY(SD, IXY[2*n+0], IXY[2*n+1], &(E.XY[2*n+0]),&(E.XY[2*n+1]));
  else
   for(int n=0; n<NXY-1; n++)
    { 
      double X1, Y1, X2, Y2;
      GetPhysicalXY(SD, IXY[2*n+0],      IXY[2*n+1],     &X1, &Y1);
      GetPhysicalXY(SD, IXY[2*(n+1)+0],  IXY[2*(n+1)+1], &X2, &Y2);

      // unit vector in width direction
      double DX = X2-X1, DY=Y2-Y1, DNorm = sqrt(DX*DX + DY*DY);
      if (DNorm==0.0) DNorm=1.0;
      double XHat = +1.0*DY / DNorm;
      double YHat = -1.0*DX / DNorm;

      E.XY[2*n+0]  = X1-0.5*W*XHat;  E.XY[2*n+1]  = Y1-0.5*W*YHat;
      int nn = 2*NXY-1-n;
      E.XY[2*nn+0] = X1+0.5*W*XHat;  E.XY[2*nn+1] = Y1+0.5*W*YHat;

      if (n==NXY-2)
       { nn=NXY-1;
         E.XY[2*nn+0] = X2-0.5*W*XHat;  E.XY[2*nn+1] = Y2-0.5*W*YHat;
         E.XY[2*nn+2] = X2+0.5*W*XHat;  E.XY[2*nn+3] = Y2+0.5*W*YHat;
       }
    }
  SD->EntitiesThisLayer.push_back(E);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddText(StatusData *SD, GDSIIData *Data, int ns, int ne)
{  
  GDSIIStruct *s  = Data->Structs[ns];
  GDSIIElement *e = s->Elements[ne];
  if (SD->CurrentLayer!=e->Layer) return;

  char Label[1000];
  snprintf(Label,1000,"Struct %s element #%i (texttype %i)",s->Name->c_str(),ne,e->TextType);

  vector<int> IXY  = e->XY;
    
  double X, Y;
  GetPhysicalXY(SD, IXY[0], IXY[1], &X, &Y);

  Entity E;
  E.XY.push_back(X);
  E.XY.push_back(Y);
  E.Text   = strdup(e->Text->c_str());
  E.Label  = strdup(Label);
  E.Closed = false;
  SD->EntitiesThisLayer.push_back(E);
}

void AddStruct(StatusData *SD, GDSIIData *Data, int ns, bool ASRef=false);

void AddASRef(StatusData *SD, GDSIIData *Data, int ns, int ne)
{
  SD->RefDepth++;

  GDSIIStruct *s   = Data->Structs[ns];
  GDSIIElement *e  = s->Elements[ne];
  vector<int> IXY  = e->XY;

  int nsRef = e->nsRef;
  if ( nsRef==-1 || nsRef>=((int)(Data->Structs.size())) )
   GDSIIData::ErrExit("structure %i (%s), element %i: REF to unknown structure %s",ns,s->Name,ne,e->SName);
    
  double Mag   = (e->Type==SREF) ? e->Mag   : 1.0;
  double Angle = (e->Type==SREF) ? e->Angle : 0.0;
  bool   Refl  = (e->Type==SREF) ? e->Refl  : false;

  GTransform GT;
  GT.CosTheta=cos(Angle*M_PI/180.0);
  GT.SinTheta=sin(Angle*M_PI/180.0);
  GT.Mag=Mag;
  GT.Refl=Refl;
  SD->GTStack.push_back(GT);
  int CurrentGT = SD->GTStack.size()-1;

  int NC=1, NR=1;
  double XYCenter[2], DeltaXYC[2]={0,0}, DeltaXYR[2]={0,0};
  XYCenter[0] = (double)IXY[0];
  XYCenter[1] = (double)IXY[1];
  if (e->Type == AREF)
   { 
     NC = e->Columns;
     NR = e->Rows;
     DeltaXYC[0] = ((double)IXY[2] - XYCenter[0]) / NC;
     DeltaXYC[1] = ((double)IXY[3] - XYCenter[1]) / NC;
     DeltaXYR[0] = ((double)IXY[4] - XYCenter[0]) / NR;
     DeltaXYR[1] = ((double)IXY[5] - XYCenter[1]) / NR;
   }

  for(int nc=0; nc<NC; nc++)
   for(int nr=0; nr<NR; nr++)
    { 
      SD->GTStack[CurrentGT].X0 = XYCenter[0] + nc*DeltaXYC[0] + nr*DeltaXYR[0];
      SD->GTStack[CurrentGT].Y0 = XYCenter[1] + nc*DeltaXYC[1] + nr*DeltaXYR[1];
      AddStruct(SD, Data, nsRef, true);
    }

  SD->GTStack.pop_back();
  SD->RefDepth--;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddElement(StatusData *SD, GDSIIData *Data, int ns, int ne)
{
  GDSIIStruct *s  = Data->Structs[ns];
  GDSIIElement *e = s->Elements[ne];

  switch(e->Type)
   { 
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case BOUNDARY:
      AddBoundary(SD, Data, ns, ne);
      break;

     case PATH:
      AddPath(SD, Data, ns, ne);
      break;

     case SREF:
     case AREF:
      AddASRef(SD, Data, ns, ne);
      break;

     case TEXT:
      AddText(SD, Data, ns, ne);
      break;

     default:
      ; // all other element types ignored for now
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddStruct(StatusData *SD, GDSIIData *Data, int ns, bool ASRef)
{
  GDSIIStruct *s=Data->Structs[ns];

  if (s->IsPCell) return;
  if (ASRef==false && s->IsReferenced) return;
  
  for(size_t ne=0; ne<s->Elements.size(); ne++)
   AddElement(SD, Data, ns, ne);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GDSIIData::Flatten(double CoordinateLengthUnit)
{
  if (CoordinateLengthUnit==0.0)
   { CoordinateLengthUnit = 1.0e-6;
     char *s=getenv("LIBGDSII_LENGTH_UNIT");
     if (s && 1==sscanf(s,"%le",&CoordinateLengthUnit))
      Log("Setting libGDSII length unit to %g meters.\n",CoordinateLengthUnit);
   }

  StatusData SD;
  InitStatusData(&SD, CoordinateLengthUnit, FileUnits[1]);
  
  for(size_t nl=0; nl<Layers.size(); nl++)
   { SD.CurrentLayer = Layers[nl];
     SD.EntitiesThisLayer.clear();
     for(size_t ns=0; ns<Structs.size(); ns++)
      AddStruct(&SD,this,ns,false);
     ETable.push_back(SD.EntitiesThisLayer);
   }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteGMSHEntity(Entity E, int Layer,
                     const char *geoFileName, FILE **pgeoFile,
                     const char *ppFileName, FILE **pppFile)
{   if ( (E.Text && !ppFileName) || (!E.Text && !geoFileName) ) return;
  
  if (E.Text)
   { FILE *ppFile = *pppFile;
     if (!ppFile)
      ppFile = *pppFile = fopen(ppFileName,"w");
     fprintf(ppFile,"View \"Layer %i %s\" {\n",Layer,E.Label);
     fprintf(ppFile,"T3 (%e,%e,%e,0) {\"%s\"};\n",E.XY[0],E.XY[1],0.0,E.Text);
     fprintf(ppFile,"};\n");
   }
  else
   { FILE *geoFile = *pgeoFile;
     if (!geoFile)
      geoFile = *pgeoFile = fopen(geoFileName,"w");
     fprintf(geoFile,"// Layer %i %s \n",Layer,E.Label);
     if (!geoFile) { fprintf(stderr,"could not open file %s (aborting)\n",geoFileName); exit(1); }

     static int NumLines=0, NumSurfaces=0, NumNodes=0;

     int Node0 = NumNodes, Line0=NumLines, NXY = E.XY.size() / 2;
    
     for(int n=0; n<NXY; n++)
      fprintf(geoFile,"Point(%i)={%e,%e,%e};\n",NumNodes++,E.XY[2*n+0],E.XY[2*n+1],0.0);
     for(int n=0; n<NXY-1; n++)
      fprintf(geoFile,"Line(%i)={%i,%i};\n",NumLines++,Node0+n,Node0+((n+1)%NXY));

     if (E.Closed)
      { fprintf(geoFile,"Line(%i)={%i,%i};\n",NumLines++,Node0+NXY-1,Node0);
        fprintf(geoFile,"Line Loop(%i)={",NumSurfaces++);
        for(int n=0; n<NXY; n++)
         fprintf(geoFile,"%i%s",Line0+n,(n==NXY-1) ? "};\n" : ",");
        fprintf(geoFile,"Plane Surface(%i)={%i};\n",NumSurfaces-1,NumSurfaces-1);
      }
     fprintf(geoFile,"\n");
   }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteGMSHFile(EntityTable ETable, iVec Layers, const char *FileBase, bool SeparateLayers)
{
  char ppFileName[100];
  snprintf(ppFileName,100,"%s.pp",FileBase);
  FILE  *ppFile = 0;

  char geoFileName[100];
  FILE  *geoFile = 0;
  if (!SeparateLayers)
   snprintf(geoFileName,100,"%s.geo",FileBase);

  for(size_t nl=0; nl<Layers.size(); nl++)
   { 
     int Layer = Layers[nl];

     if (SeparateLayers)
      snprintf(geoFileName,100,"%s.Layer%i.geo",FileBase,Layer);

     for(size_t ne=0; ne<ETable[nl].size(); ne++)
      WriteGMSHEntity(ETable[nl][ne], Layer, geoFileName, &geoFile, ppFileName, &ppFile);

     if (SeparateLayers && geoFile)
      { fclose(geoFile);
        geoFile=0;
        printf("Wrote GMSH geometry file for layer %i to %s.\n",Layer,geoFileName);
      }
   }
 
  if (geoFile)
   { fclose(geoFile);
     printf("Wrote GMSH geometry file to %s.\n",geoFileName);
   }
  if (ppFile)
   { fclose(ppFile);
     printf("Wrote GMSH post-processing file to %s.\n",ppFileName);
   }
  printf("Thank you for your support.\n");
  
}

} // namespace libGDSII


/* Copyright (C) 2005-2017 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

/*
 * libGDSII.cc -- main source file for libGDSII
 * Homer Reid   11/2017
 */
namespace libGDSII {

/***************************************************************/
/* GDSIIData constructor: create a new GDSIIData instance from */
/* a binary GDSII file.                                        */
/***************************************************************/
GDSIIData::GDSIIData(const string FileName)
{ 
  // initialize class data
  LibName       = 0;
  FileUnits[0]  = 1.0e-3; // these seem to be the default for GDSII files
  FileUnits[1]  = 1.0e-9;
  UnitInMeters  = 1.0e-6;
  GDSIIFileName = new string(FileName);
  ReadGDSIIFile(FileName);

  // at this point ErrMsg is non-null if an error occurred
  if (ErrMsg) return;
}

GDSIIData::~GDSIIData()
{
  if (GDSIIFileName) delete GDSIIFileName;
  if (ErrMsg) delete ErrMsg;
  for(size_t ns=0; ns<Structs.size(); ns++)
   { for(size_t ne=0; ne<Structs[ns]->Elements.size(); ne++)
      { if (Structs[ns]->Elements[ne]->SName) delete Structs[ns]->Elements[ne]->SName;
        if (Structs[ns]->Elements[ne]->Text)  delete Structs[ns]->Elements[ne]->Text;
        delete Structs[ns]->Elements[ne];
      }
     if (Structs[ns]->Name) delete Structs[ns]->Name;
     delete Structs[ns];
   }

  for(size_t nl=0; nl<ETable.size(); nl++)
   for(size_t ne=0; ne<ETable[nl].size(); ne++)
    { if (ETable[nl][ne].Text) free(ETable[nl][ne].Text);
      if (ETable[nl][ne].Label) free(ETable[nl][ne].Label);
    }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int GDSIIData::GetStructByName(string Name)
{ for(size_t ns=0; ns<Structs.size(); ns++)
   if ( Name == *(Structs[ns]->Name) )
    return ns;
  return -1;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
iVec GDSIIData::GetLayers()
{ return Layers; }

PolygonList GDSIIData::GetPolygons(const char *Text, int Layer)
{
  PolygonList Polygons;
  
  // first pass to find text strings matching Text, if it is non-NULL
  int TextLayer=-1;
  double TextXY[2]={HUGE_VAL, HUGE_VAL};
  if (Text)
   { for(size_t nl=0; nl<Layers.size() && TextLayer==-1; nl++)
      { if (Layer!=-1 && Layers[nl]!=Layer) continue;
        for(size_t ne=0; ne<ETable[nl].size() && TextLayer==-1; ne++)
         if ( ETable[nl][ne].Text && !strcmp(ETable[nl][ne].Text,Text) )
          { TextLayer  = Layers[nl];
            TextXY[0]  = ETable[nl][ne].XY[0];
            TextXY[1]  = ETable[nl][ne].XY[1];
          }
      }
     if (TextLayer==-1) return Polygons; // text string not found, return empty list
   }

  if (TextLayer!=-1) Layer=TextLayer;

  // second pass to find matching polygons
  for(size_t nl=0; nl<Layers.size(); nl++)
   { if (Layer!=-1 && Layers[nl]!=Layer) continue;
     for(size_t ne=0; ne<ETable[nl].size(); ne++)
      { if (ETable[nl][ne].Text!=0) continue; // we want only polygons here
        if (TextLayer==-1 || PointInPolygon(ETable[nl][ne].XY, TextXY[0], TextXY[1]))
         Polygons.push_back(ETable[nl][ne].XY);
      }
   }
  return Polygons;
}

PolygonList GDSIIData::GetPolygons(int Layer) 
 { return GetPolygons(0,Layer); }

TextString NewTextString(Entity E, int Layer)
{ TextString TS;
  TS.Text  = E.Text;
  TS.XY    = dVec(E.XY);
  TS.Layer = Layer;
  return TS;
}

TextStringList GDSIIData::GetTextStrings(int Layer)
{ 
  TextStringList TextStrings;
  for(size_t nl=0; nl<Layers.size(); nl++)
   { if (Layer!=-1 && Layers[nl]!=Layer) continue;
     for(size_t ne=0; ne<ETable[nl].size(); ne++)
      if ( ETable[nl][ne].Text )
       TextStrings.push_back( NewTextString( ETable[nl][ne], Layers[nl] ) );
   }
  return TextStrings;
}

/***************************************************************/
/* the next few routines implement a mechanism by which an API */
/* code can make multiple calls to GetPolygons() for a given   */
/* GDSII file without requiring the API code to keep track of  */
/* an instance of GDSIIData, but also without re-reading the   */
/* file each time;                                             */
/***************************************************************/
static GDSIIData *CachedGDSIIData=0;

void ClearGDSIICache()
{ if (CachedGDSIIData) delete CachedGDSIIData;
  CachedGDSIIData=0;
}

void OpenGDSIIFile(const char *GDSIIFileName)
{ 
  if (CachedGDSIIData && !strcmp(CachedGDSIIData->GDSIIFileName->c_str(),GDSIIFileName) )
   return;
  else if (CachedGDSIIData) 
   ClearGDSIICache();
  CachedGDSIIData = new GDSIIData(GDSIIFileName);
  if (CachedGDSIIData->ErrMsg)
   GDSIIData::ErrExit(CachedGDSIIData->ErrMsg->c_str());
}

iVec GetLayers(const char *GDSIIFile)
{ OpenGDSIIFile(GDSIIFile);
  return CachedGDSIIData->GetLayers();
}
  
PolygonList GetPolygons(const char *GDSIIFile, const char *Label, int Layer)
{ OpenGDSIIFile(GDSIIFile);
  return CachedGDSIIData->GetPolygons(Label,Layer);
}

PolygonList GetPolygons(const char *GDSIIFile, int Layer)
 { return GetPolygons(GDSIIFile, 0, Layer); }

TextStringList GetTextStrings(const char *GDSIIFile, int Layer)
{ OpenGDSIIFile(GDSIIFile);
  return CachedGDSIIData->GetTextStrings(Layer);
}

/***************************************************************/
/* find the value of s at which the line p+s*d intersects the  */
/* line segment connecting v1 to v2 (in 2 dimensions)          */
/* algorithm: solve the 2x2 linear system p+s*d = a+t*b        */
/* where s,t are scalars and p,d,a,b are 2-vectors with        */
/* a=v1, b=v2-v1                                               */
/***************************************************************/
bool intersect_line_with_segment(double px, double py, double dx, double dy,
                                 double *v1, double *v2, double *s)
{
  double ax = v1[0],        ay  = v1[1];
  double bx = v2[0]-v1[0],  by  = v2[1]-v1[1];
  double M00  = dx,         M10 = dy;
  double M01  = -1.0*bx,    M11 = -1.0*by;
  double RHSx = ax - px,    RHSy = ay - py;
  double DetM = M00*M11 - M01*M10;
  double L2 = bx*bx + by*by; // squared length of edge
  if ( fabs(DetM) < 1.0e-10*L2 ) // d zero or nearly parallel to edge-->no intersection
   return false;

  double t = (M00*RHSy-M10*RHSx)/DetM;
  if (t<0.0 || t>1.0) // intersection of lines does not lie between vertices
   return false;

  if (s) *s = (M11*RHSx-M01*RHSy)/DetM;
  return true;
}

// like the previous routine, but only count intersections if s>=0
bool intersect_ray_with_segment(double px, double py, double dx, double dy,
                                double *v1, double *v2, double *s)
{ double ss=0.0; if (s==0) s=&ss;
  return (intersect_line_with_segment(px,py,dx,dy,v1,v2,s) && *s>0.0);
}

/***************************************************************/
/* 2D point-in-polygon test: return 1 if the point lies within */
/* the polygon with the given vertices, 0 otherwise.           */
// method: cast a plumb line in the negative y direction from  */
/* p to infinity and count the number of edges intersected;    */
/* point lies in polygon iff this is number is odd.            */
/***************************************************************/
bool PointInPolygon(dVec Vertices, double X, double Y)
{
  size_t NV = Vertices.size() / 2;
  if (NV<3) return false;
  int num_side_intersections=0;
  for(size_t nv=0; nv<NV; nv++)
   { int nvp1 = (nv+1)%NV;
     double v1[2], v2[2];
     v1[0] = Vertices[2*nv+0  ]; v1[1]=Vertices[2*nv+  1];
     v2[0] = Vertices[2*nvp1+0]; v2[1]=Vertices[2*nvp1+1];
     if (intersect_ray_with_segment(X, Y, 0, -1.0, v1, v2, 0))
      num_side_intersections++;
   }
  return (num_side_intersections%2)==1;
}

/***************************************************************/
/* utility routines from libhrutil, duplicated here to avoid   */
/* that dependency                                             */
/***************************************************************/
bool GDSIIData::Verbose=false;
char *GDSIIData::LogFileName=0;
#define MAXSTR 1000
void GDSIIData::Log(const char *format, ...)
{
  va_list ap;
  va_start(ap,format);
  char buffer[MAXSTR];
  vsnprintf(buffer,MAXSTR,format,ap);
  va_end(ap);

  FILE *f=0;
  if (LogFileName && !strcmp(LogFileName,"stderr"))
   f=stderr;
  else if (LogFileName && !strcmp(LogFileName,"stdout"))
   f=stdout;
  else if (LogFileName)
   f=fopen(LogFileName,"a");
  if (!f) return;

  time_t MyTime;
  struct tm *MyTm;
  MyTime=time(0);
  MyTm=localtime(&MyTime);
  char TimeString[30];
  strftime(TimeString,30,"%D::%T",MyTm);
  fprintf(f,"%s: %s\n",TimeString,buffer);

  if (f!=stderr && f!=stdout) fclose(f);
}

void GDSIIData::ErrExit(const char *format, ...)
{
  va_list ap; 
  char buffer[MAXSTR];

  va_start(ap,format);
  vsnprintf(buffer,MAXSTR,format,ap);
  va_end(ap);

  fprintf(stderr,"error: %s (aborting)\n",buffer);
  Log("error: %s (aborting)",buffer);

  exit(1);
}

void GDSIIData::Warn(const char *format, ...)
{
  va_list ap;
  char buffer[MAXSTR];

  va_start(ap,format);
  vsnprintf(buffer,MAXSTR,format,ap);
  va_end(ap);

  if (Verbose) 
   fprintf(stderr,"**warning: %s \n",buffer);
  Log("warning: %s \n",buffer);

}

char *GDSIIData::vstrappend(char *s, const char *format, ...)
{
  va_list ap;
  char buffer[MAXSTR];

  va_start(ap,format);
  vsnprintf(buffer,MAXSTR,format,ap);
  va_end(ap);

  if (s==0)
   s=strdup(buffer);
  else
   { int NS=strlen(s), NB=strlen(buffer);
     s = (char *)realloc(s, NS+NB+1);
     strcpy(s + NS, buffer);
   }
  return s;
}

char *GDSIIData::vstrdup(const char *format, ...)
{
  va_list ap;
  char buffer[MAXSTR];

  va_start(ap,format);
  vsnprintf(buffer,MAXSTR,format,ap);
  va_end(ap);
  return strdup(buffer);
}

} // namespace libGSDII

/***************************************************************/
/* crutch function to play nice with autotools *****************/
/***************************************************************/
extern "C" {
void libGDSIIExists(){}
}
