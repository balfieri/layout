// Copyright (c) 2017-2019 Robert A. Alfieri
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
// 
// Layout.h - 3D Layout reader/writer in one .h file
//
// How to use it:
//
//     0) Clone this repository
//
//     1) #include "Layout.h"
//
//        Layout * layout = new Layout( "my_chip.aedt" );
//        if ( !layout->is_good ) {
//            std::cout << "Layout load failed with error: " << layout->error_msg << "\n";
//            exit( 1 );
//        }
//        layout->write( "my_chip.layout" );    // will write out the self-contained binary layout layout
//
//     2) After that, you can quickly read in the single binary layout file using:
//
//        Layout * layout = new Layout( "my_chip.layout" );  
//        if ( !layout->is_good ) {
//            std::cout << "Layout load failed with error: " << layout->error_msg << "\n";
//            exit( 1 );
//        }
//
//     3) You can also write out (export) other types of files:
//      
//        layout->write( "new_chip.aedt" );     // writes out an .aedt files
//        layout->write( "new_chip.gds" );      // writes out a .gds II file
//        layout->write( "new_chip.lst" );      // writes out a .lst file for FastCap2
//        layout->write( "new_chip.henry" );    // writes out a FastHenry2 files
//
// How it works:
//
//     1) Allocate large virtual memory 1D arrays for records and other structures.
//        These are allocated on a page boundary to make uncompressed writes faster.
//        These arrays are dynamically resized.
//     2) Read entire .aedt or other file into memory (the o/s should effectively make this work like an mmap).
//     3) Parse .aedt file using custom parser that goes character-by-character and does its own number conversions. 
//     4) Add elements to 1D arrays.  
//     5) Write to  uncompressed file is fast because all structures are aligned on a page boundary in mem and in file.
//     6) Read from uncompressed file is fast because all structures are aligned on a page boundary in mem and in file.
//        The O/S will do the equivalent of an mmap() for each array.
//
// Note: this code was taken from Bob's 3D graphics model reader/writer https://github.com/balfieri/gfx3d.
//
#ifndef _Layout_h
#define _Layout_h

#include <cstdint>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

#include <errno.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <zlib.h>

class Layout
{
public:
    typedef uint32_t uint;                  // by default, we use 32-bit integers
    typedef int32_t  _int;                  // by default, we use 32-bit integers
    typedef double   real;

    // static functions for manipulating file paths
    static void        dissect_path( std::string path, std::string& dir_name, std::string& base_name, std::string& ext_name ); // utility
    static std::string path_without_ext( std::string path, std::string * file_ext=nullptr );  // returns path without file extension part

    // constructor reads in file
    Layout( std::string top_file );
    ~Layout(); 

    bool write( std::string file_path );

    class Header                            // header (of future binary file)
    {
    public:
        uint        version;                // version
        uint        byte_cnt;               // total in-memory bytes including this header

        uint        char_cnt;               // in strings array
        uint        material_cnt;           // in materials array
        uint        layer_cnt;              // in layers array
        uint        node_cnt;               // in nodes array  
        uint        structure_cnt;          // in structures array
        uint        instance_cnt;           // in instances array
        uint        root_i;                 // index of root node in nodes array
    };

    // returns index into strings[] for s
    // creates new index if s is not already in strings[]
    uint str_get( std::string s );

    class Material
    {
    public:
        uint        name_i;                 // index of material name in strings[]
        real        relative_permittivity;           
        real        permeability;
        real        conductivity;
        real        thermal_conductivity;
        real        mass_density;
        real        specific_heat;
        real        youngs_modulus;
        real        poissons_ratio;
        real        thermal_expansion_coefficient;
    };

    // these return index in materials[] array or uint(-1) when failure or not found
    // set() will override material properties if name already exists
    uint        material_set( std::string name, const Material& material );
    uint        material_get( std::string name ) const;         

    class Layer                             // layer mapping
    {
    public:
        uint        name_i;                 // index of layer name in strings[]
        uint        gdsii_num;              // layer number of main material in GDSII file
        uint        dielectric_gdsii_num;   // layer number of dielectric material in GDSII file
        real        thickness;              // thickness in um
        uint        material_i;             // index of main material in materials[]
        uint        dielectric_material_i;  // index of dielectric material in materials[]
    };

    // these return index in layers[] array or uint(-1) when failure or not found
    // set() will override layer properties if name already exists
    uint        layer_set( uint layer_i, const Layer& layer );
    uint        layer_get( std::string name ) const;         

    enum class NODE_KIND
    {
        STR,                                // scalars
        BOOL,
        INT,
        UINT,
        REAL,
        ID,

        ASSIGN,                             // child 0 is lhs (usually an id), child 1 is rhs
        HIER,                               // child 0 is id, other children are normal children
        CALL,                               // child 0 is id, other children are args
        SLICE,                              // child 0 is id, child 1 is index before the ':', other children are other args

        GDSII_HEADER = 0x100,               // this is a GDSII HEADER node; other GDSII_KINDs follow sequentially
    };
            
    static std::string  str( NODE_KIND kind );
    
    class Node
    {
    public:
        NODE_KIND   kind;                   // kind of node
        union {
            uint        s_i;                // STR or ID - index into std::string array of std::string value
            bool        b;                  // BOOL
            _int        i;                  // INT
            uint        u;                  // UINT
            real        r;                  // REAL
            uint        child_first_i;      // non-scalars - index into nodes[] array of first child on list
        } u;
        uint        sibling_i;              // index in nodes array of sibling on list, else uint(-1)
    };

    // attempt to find name for a node
    std::string name( const Node& node ) const;

    enum class GDSII_KIND                   // these are in the order defined by the GDSII spec
    {
        HEADER,
        BGNLIB,
        LIBNAME,
        UNITS,
        ENDLIB,
        BGNSTR,
        STRNAME,
        ENDSTR,
        BOUNDARY,
        PATH,
        SREF,
        AREF,
        TEXT,
        LAYER,
        DATATYPE,
        WIDTH,
        XY,
        ENDEL,
        SNAME,
        COLROW,
        TEXTNODE,
        NODE,
        TEXTTYPE,
        PRESENTATION,
        UNUSED,
        STRING,
        STRANS,
        MAG,
        ANGLE,
        UNUSED2,
        UNUSED3,
        REFLIBS,
        FONTS,
        PATHTYPE,
        GENERATIONS,
        ATTRTABLE,
        STYPTABLE,
        STRTYPE,
        ELFLAGS,
        ELKEY,
        LINKTYPE,
        LINKKEYS,
        NODETYPE,
        PROPATTR,
        PROPVALUE,
        BOX,
        BOXTYPE,
        PLEX,
        BGNEXTN,
        ENDEXTN,
        TAPENUM,
        TAPECODE,
        STRCLASS,
        RESERVED,
        FORMAT,
        MASK,
        ENDMASKS,
        LIBDIRSIZE,
        SRFNAME,
        LIBSECUR,
    };

    static constexpr uint32_t GDSII_KIND_CNT = 0x3c;
    static std::string str( GDSII_KIND kind );

    enum class GDSII_DATATYPE
    {
        NO_DATA,
        BITARRAY,
        INTEGER_2,
        INTEGER_4,
        REAL_4,
        REAL_8,
        STRING,
    };

    static GDSII_DATATYPE kind_to_datatype( GDSII_KIND kind );
    static std::string    str( GDSII_DATATYPE datatype );

    // global scalars
    static const uint   VERSION = 0xB0BA1f01; // current version 
    bool                is_good;            // set to true if constructor succeeds
    std::string         error_msg;          // if !is_good

    // structs
    std::string         file_path;          // pathname of file passed to Layout()
    uint8_t *           mapped_region;      // != nullptr means the whole file was sucked in by read_uncompressed()
    Header *            hdr;
    Header *            max;                // holds max lengths of currently allocated arrays 

    // arrays
    char *              strings;
    Material *          materials;
    Layer *             layers;
    Node *              nodes;
    uint *              structures;         // node indexes of structures
    uint *              instances;          // node indexes of instances

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//
// IMPLEMENTATION
//
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
private:
    // maps of strings to array indexes
    std::map<std::string, uint>         str_to_str_i;           // maps std::string to unique location in strings[]

    // state used during reading of files
    std::string ext_name;
    uint8_t * nnn_start;
    uint8_t * nnn_end;
    uint8_t * nnn;
    uint      line_num;

    // state used during reading and and writing of GDSII files
    uint       gdsii_rec_cnt;
    GDSII_KIND gdsii_last_kind;
    int        gdsii_fd;
    uint8_t *  gdsii_buff;
    uint       gdsii_buff_byte_cnt;

    // state used during reading and and writing of AEDT files
    uint aedt_begin_str_i;              // these are to make it easier to compare
    uint aedt_end_str_i;
    uint true_str_i;
    uint false_str_i;

    bool layout_read( std::string file_path );          // .layout
    bool layout_write( std::string file_path );         

    bool gdsii_is_hier( GDSII_KIND kind ) const;
    GDSII_KIND gdsii_hier_end_kind( GDSII_KIND kind ) const;
    bool gdsii_is_name( GDSII_KIND kind ) const;
    bool gdsii_read( std::string file_path );           // .gds
    bool gdsii_read_record( uint& node_i );
    bool gdsii_write( std::string file );
    void gdsii_write_record( uint node_i );
    void gdsii_write_number( uint8_t * bytes, uint& byte_cnt, uint ni, GDSII_DATATYPE datatype );
    void gdsii_write_bytes( const uint8_t * bytes, uint byte_cnt );
    void gdsii_flush( void );
    
    bool aedt_read( std::string file );                 // .aedt
    bool aedt_read_expr( uint& node_i );
    bool aedt_write( std::string file );
    void aedt_write_expr( std::ofstream& out, uint node_i, std::string indent_str );

    bool cmd( std::string c );
    bool open_and_read( std::string file_name, uint8_t *& start, uint8_t *& end );

    void skip_whitespace_to_eol( uint8_t *& xxx, uint8_t *& xxx_end );  // on this line only
    void skip_whitespace( uint8_t *& xxx, uint8_t *& xxx_end );
    void skip_to_eol( uint8_t *& xxx, uint8_t *& xxx_end );
    bool eol( uint8_t *& xxx, uint8_t *& xxx_end );
    bool expect_char( char ch, uint8_t *& xxx, uint8_t * xxx_end, bool skip_whitespace_first=false );
    bool parse_number( uint node_i, uint8_t *& xxx, uint8_t *& xxx_end );
    bool parse_string( std::string& s, uint8_t *& xxx, uint8_t *& xxx_end );
    bool parse_string_i( uint& s, uint8_t *& xxx, uint8_t *& xxx_end );
    bool parse_id( uint& id_i, uint8_t *& xxx, uint8_t *& xxx_end );
    bool peek_id( uint& id_i, uint8_t *& xxx, uint8_t *& xxx_end );
    bool parse_real( real& r, uint8_t *& xxx, uint8_t *& xxx_end );
    bool parse_int( _int& i, uint8_t *& xxx, uint8_t *& xxx_end );
    bool parse_uint( uint& u, uint8_t *& xxx, uint8_t *& xxx_end );
    bool parse_bool( bool& b, uint8_t *& xxx, uint8_t *& xxx_end );
    std::string surrounding_lines( uint8_t *& xxx, uint8_t *& xxx_end );

    // allocates an array of T on a page boundary
    template<typename T>
    T * aligned_alloc( size_t cnt );

    // reallocate array if we are about to exceed its current size
    template<typename T>
    inline void perhaps_realloc( T *& array, const uint  & hdr_cnt, uint  & max_cnt, uint   add_cnt );
};

#ifdef LAYOUT_DEBUG
#define ldout if ( true )  std::cout 
#else
#define ldout if ( false ) std::cout
#endif

// these are done as macros to avoid evaluating msg (it makes a big difference)
#include <assert.h>
#define rtn_assert(  bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << msg << "\n"; assert( false ); return false; }
#define rtnn_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << msg << "\n"; assert( false );               }
#define node_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << msg << "\n"; assert( false ); goto error;   }

std::string Layout::str( Layout::GDSII_KIND kind )
{
    #define gcase( kind ) case Layout::GDSII_KIND::kind: return #kind;
    switch( kind )
    {
        gcase( HEADER )
        gcase( BGNLIB )
        gcase( LIBNAME )
        gcase( UNITS )
        gcase( ENDLIB )
        gcase( BGNSTR )
        gcase( STRNAME )
        gcase( ENDSTR )
        gcase( BOUNDARY )
        gcase( PATH )
        gcase( SREF )
        gcase( AREF )
        gcase( TEXT )
        gcase( LAYER )
        gcase( DATATYPE )
        gcase( WIDTH )
        gcase( XY )
        gcase( ENDEL )
        gcase( SNAME )
        gcase( COLROW )
        gcase( TEXTNODE )
        gcase( NODE )
        gcase( TEXTTYPE )
        gcase( PRESENTATION )
        gcase( UNUSED )
        gcase( STRING )
        gcase( STRANS )
        gcase( MAG )
        gcase( ANGLE )
        gcase( UNUSED2 )
        gcase( UNUSED3 )
        gcase( REFLIBS )
        gcase( FONTS )
        gcase( PATHTYPE )
        gcase( GENERATIONS )
        gcase( ATTRTABLE )
        gcase( STYPTABLE )
        gcase( STRTYPE )
        gcase( ELFLAGS )
        gcase( ELKEY )
        gcase( LINKTYPE )
        gcase( LINKKEYS )
        gcase( NODETYPE )
        gcase( PROPATTR )
        gcase( PROPVALUE )
        gcase( BOX )
        gcase( BOXTYPE )
        gcase( PLEX )
        gcase( BGNEXTN )
        gcase( ENDEXTN )
        gcase( TAPENUM )
        gcase( TAPECODE )
        gcase( STRCLASS )
        gcase( RESERVED )
        gcase( FORMAT )
        gcase( MASK )
        gcase( ENDMASKS )
        gcase( LIBDIRSIZE )
        gcase( SRFNAME )
        gcase( LIBSECUR )
        default: return "<unknown>"; 
    }
}

std::string Layout::str( Layout::GDSII_DATATYPE datatype )
{
    #define dcase( type ) case GDSII_DATATYPE::type: return #type;
    switch( datatype ) 
    {
        dcase( NO_DATA )
        dcase( BITARRAY )
        dcase( INTEGER_2 )
        dcase( INTEGER_4 )
        dcase( REAL_4 )
        dcase( REAL_8 )
        dcase( STRING )
        default: return "<unknown>";
    }
}

std::string Layout::str( Layout::NODE_KIND kind )
{
    #define ncase( kind ) case Layout::NODE_KIND::kind: return #kind; 
    switch( kind )
    {
        ncase( STR )
        ncase( BOOL )
        ncase( INT )
        ncase( UINT )
        ncase( REAL )
        ncase( ID )
        ncase( HIER )
        ncase( ASSIGN )
        ncase( CALL )
        ncase( SLICE )
        default: 
        {
            if ( kind >= Layout::NODE_KIND::GDSII_HEADER ) {
                Layout::GDSII_KIND gkind = Layout::GDSII_KIND(int(kind) - int(Layout::NODE_KIND::GDSII_HEADER));
                return "GDSII_" + Layout::str(gkind);
            } else {
                return "<unknown>"; 
            }
        }
    }
}

inline std::ostream& operator << ( std::ostream& os, const Layout::NODE_KIND& kind ) 
{
    os << Layout::str( kind );
    return os;
}

Layout::GDSII_DATATYPE Layout::kind_to_datatype( Layout::GDSII_KIND kind )
{
    #define kdcase( kind, type ) case GDSII_KIND::kind: return GDSII_DATATYPE::type;
    switch( kind )
    {
        kdcase( HEADER,       INTEGER_2 )
        kdcase( BGNLIB,       INTEGER_2 )
        kdcase( LIBNAME,      STRING )
        kdcase( UNITS,        REAL_8 )
        kdcase( ENDLIB,       NO_DATA )
        kdcase( BGNSTR,       INTEGER_2 )
        kdcase( STRNAME,      STRING )
        kdcase( ENDSTR,       NO_DATA )
        kdcase( BOUNDARY,     NO_DATA )
        kdcase( PATH,         NO_DATA )
        kdcase( SREF,         NO_DATA )
        kdcase( AREF,         NO_DATA )
        kdcase( TEXT,         NO_DATA )
        kdcase( LAYER,        INTEGER_2 )
        kdcase( DATATYPE,     INTEGER_2 )
        kdcase( WIDTH,        INTEGER_4 )
        kdcase( XY,           INTEGER_4 )
        kdcase( ENDEL,        NO_DATA )
        kdcase( SNAME,        STRING )
        kdcase( COLROW,       INTEGER_2 )
        kdcase( TEXTNODE,     NO_DATA )
        kdcase( NODE,         NO_DATA )
        kdcase( TEXTTYPE,     INTEGER_2 )
        kdcase( PRESENTATION, BITARRAY )
        kdcase( UNUSED,       NO_DATA )
        kdcase( STRING,       STRING )
        kdcase( STRANS,       BITARRAY )
        kdcase( MAG,          REAL_8 )
        kdcase( ANGLE,        REAL_8 )
        kdcase( UNUSED2,      NO_DATA )
        kdcase( UNUSED3,      NO_DATA )
        kdcase( REFLIBS,      STRING )
        kdcase( FONTS,        STRING )
        kdcase( PATHTYPE,     INTEGER_2 )
        kdcase( GENERATIONS,  INTEGER_2 )
        kdcase( ATTRTABLE,    STRING )
        kdcase( STYPTABLE,    STRING )
        kdcase( STRTYPE,      INTEGER_2 )
        kdcase( ELFLAGS,      BITARRAY )
        kdcase( ELKEY,        INTEGER_4 )
        kdcase( LINKTYPE,     NO_DATA )
        kdcase( LINKKEYS,     NO_DATA )
        kdcase( NODETYPE,     INTEGER_2 )
        kdcase( PROPATTR,     INTEGER_2 )
        kdcase( PROPVALUE,    STRING )
        kdcase( BOX,          NO_DATA )
        kdcase( BOXTYPE,      INTEGER_2 )
        kdcase( PLEX,         INTEGER_4 )
        kdcase( BGNEXTN,      INTEGER_4 )
        kdcase( ENDEXTN,      INTEGER_4 )
        kdcase( TAPENUM,      INTEGER_2 )
        kdcase( TAPECODE,     INTEGER_2 )
        kdcase( STRCLASS,     BITARRAY )
        kdcase( RESERVED,     INTEGER_4 )
        kdcase( FORMAT,       INTEGER_2 )
        kdcase( MASK,         STRING )
        kdcase( ENDMASKS,     NO_DATA )
        kdcase( LIBDIRSIZE,   INTEGER_2 )
        kdcase( SRFNAME,      STRING )
        kdcase( LIBSECUR,     INTEGER_2 )
        default: return GDSII_DATATYPE::NO_DATA;
    }
}

Layout::Layout( std::string top_file )
{
    is_good = false;
    error_msg = "<unknown error>";
    mapped_region = nullptr;
    nnn = nullptr;
    line_num = 1;
    file_path = top_file;

    hdr = aligned_alloc<Header>( 1 );
    memset( hdr, 0, sizeof( Header ) );
    hdr->version = VERSION;

    //------------------------------------------------------------
    // Initial lengths of arrays are large in virtual memory
    //------------------------------------------------------------
    max = aligned_alloc<Header>( 1 );
    max->node_cnt =  1024;
    max->char_cnt = max->node_cnt * 128;
    max->structure_cnt = 64;
    max->instance_cnt = 64;

    //------------------------------------------------------------
    // Read depends on file ext_name
    //------------------------------------------------------------
    std::string dir_name;
    std::string base_name;
    dissect_path( top_file, dir_name, base_name, ext_name );
    if ( ext_name == std::string( ".layout" ) ) {
        //------------------------------------------------------------
        // Read uncompressed .layout
        //------------------------------------------------------------
        if ( !layout_read( top_file ) ) return;
    } else {
        //------------------------------------------------------------
        // Allocate initial arrays
        //------------------------------------------------------------
        strings    = aligned_alloc<char>( max->char_cnt );
        nodes      = aligned_alloc<Node>( max->node_cnt );
        structures = aligned_alloc<uint>( max->structure_cnt );
        instances  = aligned_alloc<uint>( max->instance_cnt );

        if ( ext_name == std::string( ".gds" ) ) {
            if ( !gdsii_read( top_file ) ) return;
        } else if ( ext_name == std::string( ".aedt" ) ) {
            if ( !aedt_read( top_file ) ) return;
        } else {
            error_msg = "unknown top file ext_name: " + ext_name;
            return;
        }

        //------------------------------------------------------------
        // Add up byte count.
        //------------------------------------------------------------
        hdr->byte_cnt = uint  ( 1                 ) * sizeof( hdr ) +
                        uint  ( hdr->node_cnt     ) * sizeof( nodes[0] ) +
                        uint  ( hdr->char_cnt     ) * sizeof( strings[0] );

        is_good = true;
    }
}

Layout::~Layout()
{
    if ( mapped_region != nullptr ) {
        delete mapped_region;
        mapped_region = nullptr;
    } else {
        delete nodes;
        delete strings;
    }
}

bool Layout::write( std::string top_file )
{
    //------------------------------------------------------------
    // Write depends on file ext_name
    //------------------------------------------------------------
    std::string dir_name;
    std::string base_name;
    dissect_path( top_file, dir_name, base_name, ext_name );
    if ( ext_name == std::string( ".layout" ) ) {
        //------------------------------------------------------------
        // Write uncompressed .layout
        //------------------------------------------------------------
        return layout_write( top_file );
    } else {
        if ( ext_name == std::string( ".aedt" ) ) {
            return aedt_write( top_file );
        } else if ( ext_name == std::string( ".gds" ) ) {
            return gdsii_write( top_file );
        } else {
            error_msg = "unknown file ext_name: " + ext_name;
            return false;
        }
    }
}

std::string Layout::name( const Layout::Node& node ) const
{
    switch( node.kind ) 
    {
        case NODE_KIND::HIER:
        {
            break;
        }

        default:
        {
            if ( int(node.kind) >= int(NODE_KIND::GDSII_HEADER) ) {
                GDSII_KIND gkind = GDSII_KIND( int(node.kind) - int(NODE_KIND::GDSII_HEADER) );
                if ( gdsii_is_name( gkind ) ) return std::string( &strings[node.u.s_i] );
                if ( gdsii_is_hier( gkind ) ) {
                    for( uint child_i = node.u.child_first_i; child_i != uint(-1); child_i = nodes[child_i].sibling_i ) 
                    {
                        std::string n = name( nodes[child_i] );
                        if ( n != "" ) return n;
                    }
                }
            } 
            break;
        }
    }

    return "";
}

// returns array of T on a page boundary
template<typename T>
inline T * Layout::aligned_alloc( size_t cnt )
{
    void * mem = nullptr;
    posix_memalign( &mem, getpagesize(), cnt*sizeof(T) );
    return reinterpret_cast<T *>( mem );
}

// reallocate array if we are about to exceed its current size
template<typename T>
inline void Layout::perhaps_realloc( T *& array, const Layout::uint& hdr_cnt, Layout::uint& max_cnt, Layout::uint add_cnt )
{
    while( (hdr_cnt + add_cnt) > max_cnt ) {
        void * mem = nullptr;
        uint   old_max_cnt = max_cnt;
        max_cnt *= 2;
        if ( max_cnt < old_max_cnt ) {
            assert( old_max_cnt != uint(-1) );
            max_cnt = uint(-1);
        }
        T * new_array = aligned_alloc<T>( max_cnt );
        memcpy( new_array, array, hdr_cnt*sizeof(T) );
        delete array;
        array = new_array;
    }
}

bool Layout::layout_read( std::string layout_path )
{
    uint8_t * start;
    uint8_t * end;
    if ( !open_and_read( layout_path, start, end ) ) return false;
    mapped_region = start;

    //------------------------------------------------------------
    // Read header then individual arrays.
    // Each is padded out to a page boundary in the file.
    //------------------------------------------------------------
    uint8_t * _addr = start;
    size_t page_size = getpagesize();

    #define _uread( array, type, cnt ) \
        if ( (cnt) == 0 ) { \
            array = nullptr; \
        } else { \
            array = reinterpret_cast<type *>( _addr ); \
            size_t _byte_cnt = (cnt)*sizeof(type); \
            _byte_cnt += _byte_cnt % page_size; \
            _addr += _byte_cnt; \
        } \

    _uread( hdr,         Header,   1 );
    if ( hdr->version != VERSION ) {
        rtn_assert( 0, "hdr->version does not match VERSION=" + std::to_string(VERSION) + ", got " + std::to_string(hdr->version) );
    }
    max = aligned_alloc<Header>( 1 );
    memcpy( max, hdr, sizeof( Header ) );
    _uread( strings,     char,        hdr->char_cnt );
    _uread( nodes,       Node,        hdr->node_cnt );
    _uread( structures,  uint,        hdr->structure_cnt );
    _uread( instances,   uint,        hdr->instance_cnt );

    is_good = true;

    return true;
}

bool Layout::layout_write( std::string layout_path ) 
{
    cmd( "rm -f " + layout_path );
    int fd = open( layout_path.c_str(), O_CREAT|O_WRONLY|O_TRUNC|O_SYNC|S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP );
    if ( fd < 0 ) std::cout << "open() for write error: " << strerror( errno ) << "\n";
    rtn_assert( fd >= 0, "could not open() file " + layout_path + " for writing - open() error: " + strerror( errno ) );

    //------------------------------------------------------------
    // Write out header then individual arrays.
    // Each is padded out to a page boundary in the file.
    //------------------------------------------------------------
    size_t page_size = getpagesize();

    #define _uwrite( addr, byte_cnt ) \
    { \
        size_t _byte_cnt = byte_cnt; \
        _byte_cnt += _byte_cnt % page_size; \
        uint8_t * _addr = reinterpret_cast<uint8_t *>( addr ); \
        for( ; _byte_cnt != 0;  ) \
        { \
            uint _this_byte_cnt = 1024*1024*1024; \
            if ( _byte_cnt < _this_byte_cnt ) _this_byte_cnt = _byte_cnt; \
            if ( ::write( fd, _addr, _this_byte_cnt ) <= 0 ) { \
                close( fd ); \
                rtn_assert( 0, "could not write() file " + layout_path + " - write() error: " + strerror( errno ) ); \
            } \
            _byte_cnt -= _this_byte_cnt; \
            _addr     += _this_byte_cnt; \
        } \
    } \

    _uwrite( hdr,         1                  * sizeof(hdr[0]) );
    _uwrite( strings,     hdr->char_cnt      * sizeof(strings[0]) );
    _uwrite( nodes,       hdr->node_cnt      * sizeof(nodes[0]) );
    _uwrite( structures,  hdr->structure_cnt * sizeof(structures[0]) );
    _uwrite( instances,   hdr->instance_cnt  * sizeof(instances[0]) );

    fsync( fd ); // flush
    close( fd );
    cmd( "chmod +rw " + layout_path );

    return true;
}

bool Layout::gdsii_is_hier( GDSII_KIND kind ) const
{
    switch( kind )
    {
        case GDSII_KIND::BGNLIB:
        case GDSII_KIND::BGNSTR:
        case GDSII_KIND::BOUNDARY:
        case GDSII_KIND::PATH:
        case GDSII_KIND::SREF:
        case GDSII_KIND::AREF:
        case GDSII_KIND::TEXT:
        case GDSII_KIND::NODE:
            return true;

        default:
            return false;
    }
}

Layout::GDSII_KIND Layout::gdsii_hier_end_kind( Layout::GDSII_KIND kind ) const
{
    switch( kind )
    {
        case GDSII_KIND::BGNLIB:                
            return GDSII_KIND::ENDLIB;

        case GDSII_KIND::BGNSTR:
            return GDSII_KIND::ENDSTR;

        case GDSII_KIND::BOUNDARY:
        case GDSII_KIND::PATH:
        case GDSII_KIND::SREF:
        case GDSII_KIND::AREF:
        case GDSII_KIND::TEXT:
        case GDSII_KIND::NODE:
            return GDSII_KIND::ENDEL;

        default:
            assert( false );
            return GDSII_KIND::ENDEL;  // for compiler
    }
}

bool Layout::gdsii_is_name( Layout::GDSII_KIND kind ) const
{
    switch( kind )
    {
        case GDSII_KIND::LIBNAME:
        case GDSII_KIND::STRNAME:
        case GDSII_KIND::SNAME:
        case GDSII_KIND::SRFNAME:
            return true;

        default:
            return false;
    }
}

bool Layout::gdsii_read( std::string file )
{
    //------------------------------------------------------------
    // Map in file
    //------------------------------------------------------------
    line_num = 1;
    if ( !open_and_read( file, nnn_start, nnn_end ) ) return false;
    nnn = nnn_start;

    //------------------------------------------------------------
    // Create a file-level HIER node and read in all of them.
    //------------------------------------------------------------
    perhaps_realloc( nodes, hdr->node_cnt, max->node_cnt, 1 );
    uint ni = hdr->node_cnt++;
    nodes[ni].kind = NODE_KIND::HIER;
    nodes[ni].u.child_first_i = uint(-1);
    nodes[ni].sibling_i = uint(-1);
    hdr->root_i = ni;
    uint prev_i = uint(-1);
    gdsii_rec_cnt = 0;
    while( nnn < nnn_end )
    {
        uint child_i;
        if ( !gdsii_read_record( child_i ) ) return false;

        if ( prev_i == uint(-1) ) {
            nodes[ni].u.child_first_i = child_i;
        } else {
            nodes[prev_i].sibling_i = child_i;
        }

        if ( gdsii_last_kind == GDSII_KIND::ENDLIB ) break;

        prev_i = child_i;
    }
    return true;
}

static inline bool is_gdsii_allowed_char( char c ) { return isprint( c ) && c != '"' && c != ','; }

bool Layout::gdsii_read_record( uint& ni )
{
    //------------------------------------------------------------
    // Parse record header.
    //------------------------------------------------------------
    rtn_assert( (nnn + 4) <= nnn_end, "unexpected end of gdsii file rec_cnt=" + std::to_string(gdsii_rec_cnt) );
    uint32_t       byte_cnt = ( nnn[0] << 8 ) | nnn[1];
    GDSII_KIND     kind     = GDSII_KIND( nnn[2] );
    GDSII_DATATYPE datatype = GDSII_DATATYPE( nnn[3] );
    rtn_assert( byte_cnt >= 4, std::to_string(gdsii_rec_cnt) + ": gdsii record byte_cnt must be at least 4, byte_cnt=" + std::to_string(byte_cnt) + " kind=" + str(kind) );
    ldout << str(kind) << " " << str(datatype) << " byte_cnt=" << std::to_string(byte_cnt) << "\n";
    byte_cnt -= 4;
    nnn += 4;
    rtn_assert( uint32_t(kind) < GDSII_KIND_CNT, std::to_string(gdsii_rec_cnt) + ": bad gdsii record kind " + std::to_string(uint32_t(kind)) );
    rtn_assert( kind_to_datatype( kind ) == datatype, 
                std::to_string(gdsii_rec_cnt) + ": datatype=" + str(datatype) + " does not match expected datatype=" + 
                str( kind_to_datatype( kind ) ) + " for record kind " + str(kind) );
    rtn_assert( (nnn + byte_cnt) <= nnn_end, "unexpected end of gdsii file" );
  
    //------------------------------------------------------------
    // Create a GDSII node.
    //------------------------------------------------------------
    perhaps_realloc( nodes, hdr->node_cnt, max->node_cnt, 1 );
    ni = hdr->node_cnt++;
    nodes[ni].kind = NODE_KIND( int(NODE_KIND::GDSII_HEADER) + int(kind) );
    nodes[ni].sibling_i = uint(-1);

    //------------------------------------------------------------
    // Parse payload.
    //------------------------------------------------------------
    bool is_hier = gdsii_is_hier( kind );
    uint prev_i = uint(-1);
    switch( datatype )
    {
        case GDSII_DATATYPE::NO_DATA:
        {
            rtn_assert( byte_cnt == 0, "NO_DATA gdsii datatype should have no payload" );
            break;
        }

        case GDSII_DATATYPE::BITARRAY:
        {
            assert( !is_hier );
            rtn_assert( byte_cnt == 2, "BITARRAY gdsii datatype should have 2-byte payload" );
            nodes[ni].u.u = (nnn[1] << 8) | nnn[0];
            break;
        }

        case GDSII_DATATYPE::STRING:
        {
            assert( !is_hier );
            char c[1024];
            rtn_assert( byte_cnt <= (sizeof(c)-1), "STRING too big byte_cnt=" + std::to_string(byte_cnt) );
            if ( byte_cnt > 0 ) memcpy( c, nnn, byte_cnt );
            c[byte_cnt] = '\0';
            for( int i = byte_cnt-1; i >= 0; i-- ) 
            {
                if ( is_gdsii_allowed_char( c[i] ) ) break;
                c[i] = '\0';
            }
            std::string s = std::string( c );
            ldout << "    " << s << "\n";
            nodes[ni].u.s_i = str_get( s );
            break;
        }

        case GDSII_DATATYPE::INTEGER_2:
        case GDSII_DATATYPE::INTEGER_4:
        case GDSII_DATATYPE::REAL_4:
        case GDSII_DATATYPE::REAL_8:
        {
            uint datum_byte_cnt = (datatype == GDSII_DATATYPE::INTEGER_2) ? 2 :
                                  (datatype == GDSII_DATATYPE::REAL_8)    ? 8 : 4;
            uint cnt = byte_cnt / datum_byte_cnt;
            rtn_assert( (cnt*datum_byte_cnt) == byte_cnt, "datum_byte_cnt does not divide evenly" );
            uint8_t * uuu = nnn;
            for( uint i = 0; i < cnt; i++, uuu += datum_byte_cnt )
            {
                perhaps_realloc( nodes, hdr->node_cnt, max->node_cnt, 1 );
                uint child_i = hdr->node_cnt++;
                nodes[child_i].sibling_i = uint(-1);
                if ( i == 0 ) {
                    nodes[ni].u.child_first_i = child_i;
                } else {
                    nodes[prev_i].sibling_i = child_i;
                }
                if ( datatype == GDSII_DATATYPE::INTEGER_2 || datatype == GDSII_DATATYPE::INTEGER_4 ) {
                    int64_t vi = (uuu[0] << 8) | uuu[1];
                    if ( datatype == GDSII_DATATYPE::INTEGER_4 ) {
                        vi = (vi << 16) | (uuu[2] << 8) | uuu[3];
                    }
                    if ( (vi & 0x80000000) != 0 ) {
                        vi = (datatype == GDSII_DATATYPE::INTEGER_2) ? -(0x10000 - vi) : -(0x100000000LL - vi);
                    }
                    nodes[child_i].kind = NODE_KIND::INT;
                    nodes[child_i].u.i = vi;
                    ldout << "    " << vi << "\n";
                } else {
                    real sign = (uuu[0] & 0x80) ? -1.0 : 1.0;
                    real exp  = (uuu[0] & 0x7f);
                    int64_t ifrac = 0.0;
                    for( uint j = 0; j < datum_byte_cnt; j++ )
                    {
                        ldout << "bytes[" << j << "]=" << int(uuu[j]) << "\n";
                        if ( j != 0 ) ifrac = (ifrac << 8) | uuu[j];
                    }
                    nodes[child_i].kind = NODE_KIND::REAL;
                    real rexp = 4.0*(exp-64) - 8*(datum_byte_cnt-1);
                    nodes[child_i].u.r = sign * double(ifrac) * std::pow( 2.0, rexp );
                    ldout << "sign=" << ((sign < 0.0) ? "1" : "0") << " exp=" << exp << " rexp=" << rexp << " ifrac=" << ifrac << " r=" << nodes[child_i].u.r << "\n";
                }
                prev_i = child_i;
            }
            break;
        }

        default:
        {
            rtn_assert( false, "something is wrong" );
            break;
        }
    }

    gdsii_rec_cnt++;
    gdsii_last_kind = kind;
    nnn += byte_cnt;

    if ( gdsii_is_hier( kind ) ) {
        if ( kind == GDSII_KIND::BGNSTR ) {
            // record in structures[] array
            perhaps_realloc( structures, hdr->structure_cnt, max->structure_cnt, 1 );
            structures[hdr->structure_cnt++] = ni;
        } else if ( kind == GDSII_KIND::SREF || kind == GDSII_KIND::AREF ) {
            // record in instances[] array
            perhaps_realloc( instances, hdr->instance_cnt, max->instance_cnt, 1 );
            instances[hdr->instance_cnt++] = ni;
        }

        // recurse for other children
        for( ;; ) 
        {
            uint child_i;
            if ( !gdsii_read_record( child_i ) ) return false;
            GDSII_KIND gkind = GDSII_KIND( int(nodes[child_i].kind) - int(NODE_KIND::GDSII_HEADER) );
            if ( gkind == GDSII_KIND::ENDEL || gkind == GDSII_KIND::ENDSTR || gkind == GDSII_KIND::ENDLIB ) break;
            if ( prev_i == uint(-1) ) {
                nodes[ni].u.child_first_i = child_i;
            } else {
                nodes[prev_i].sibling_i = child_i;
            }
            prev_i = child_i;
        }
    }

    return true;
}

constexpr uint gdsii_buff_alloc_byte_cnt = 128*1024*1024;      // hide overhead of write() calls

bool Layout::gdsii_write( std::string gdsii_path )
{
    gdsii_buff = aligned_alloc<uint8_t>( gdsii_buff_alloc_byte_cnt );
    gdsii_buff_byte_cnt = 0;

    cmd( "rm -f " + gdsii_path );
    gdsii_fd = open( gdsii_path.c_str(), O_CREAT|O_WRONLY|O_TRUNC|O_SYNC|S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP );
    if ( gdsii_fd < 0 ) ldout << "open() for write error: " << strerror( errno ) << "\n";

    gdsii_write_record( hdr->root_i );

    gdsii_flush();
    fsync( gdsii_fd ); // flush
    close( gdsii_fd );
    cmd( "chmod +rw " + gdsii_path );

    delete gdsii_buff;
    gdsii_buff = nullptr;

    return true;
}

void Layout::gdsii_write_record( uint ni )
{
    const Node& node = nodes[ni];
    if ( int(node.kind) >= int(NODE_KIND::GDSII_HEADER) ) {
        uint8_t bytes[64*1024];
        uint    byte_cnt = 2;   // fill in byte_cnt later

        GDSII_KIND gkind = GDSII_KIND( int(node.kind) - int(NODE_KIND::GDSII_HEADER) );
        GDSII_DATATYPE datatype = kind_to_datatype( gkind );
        ldout << str(gkind) << " " << str(datatype) << "\n";
        bytes[byte_cnt++] = int(gkind);
        bytes[byte_cnt++] = int(datatype);

        uint child_i = uint(-1);

        switch( datatype )
        {
            case GDSII_DATATYPE::NO_DATA:
            {
                break;
            }

            case GDSII_DATATYPE::BITARRAY:
            {
                ldout << "    " << node.u.u << "\n";
                bytes[byte_cnt++] = node.u.u & 0xff;
                bytes[byte_cnt++] = (node.u.u >> 8) & 0xff;
                break;
            }

            case GDSII_DATATYPE::STRING:
            {
                ldout << "    " << std::string(&strings[node.u.s_i]) << "\n";
                uint len = strlen( &strings[node.u.s_i] );
                memcpy( &bytes[byte_cnt], &strings[node.u.s_i], len );
                byte_cnt += len;
                if ( byte_cnt & 1 ) bytes[byte_cnt++] = '\0';  // must have even number of bytes
                break;
            }

            case GDSII_DATATYPE::INTEGER_2:
            case GDSII_DATATYPE::INTEGER_4:
            case GDSII_DATATYPE::REAL_4:
            case GDSII_DATATYPE::REAL_8:
            {
                for( child_i = node.u.child_first_i; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
                {
                    if ( nodes[child_i].kind != NODE_KIND::INT && nodes[child_i].kind != NODE_KIND::REAL ) break; // skip GDSII children
                    gdsii_write_number( bytes, byte_cnt, child_i, datatype );
                }
                break;
            }

            default:
            {
                assert( false );
                break;
            }
        }

        // record byte count
        bytes[0] = (byte_cnt >> 8) & 0xff;
        bytes[1] = byte_cnt & 0xff;

        // transfer bytes to buffer
        gdsii_write_bytes( bytes, byte_cnt );

        if ( gdsii_is_hier( gkind ) ) {
            // recurse for rest of children
            if ( child_i == uint(-1) ) child_i = node.u.child_first_i;
            for( ; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
            {
                gdsii_write_record( child_i );
            }

            bytes[0] = 0;
            bytes[1] = 4;
            bytes[2] = uint8_t(gdsii_hier_end_kind(gkind));
            bytes[3] = uint8_t(GDSII_DATATYPE::NO_DATA);
            gdsii_write_bytes( bytes, 4 );
        }

    } else if ( node.kind == NODE_KIND::HIER ) {
        // assume file wrapper, just loop through children
        for( uint child_i = node.u.child_first_i; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
        {
            gdsii_write_record( child_i );
        }

    } else {
        rtnn_assert( false, "ignoring node kind " + str(node.kind) );
    }
}

void Layout::gdsii_write_number( uint8_t * bytes, uint& byte_cnt, uint ni, GDSII_DATATYPE datatype )
{
    switch( datatype )
    {
        case GDSII_DATATYPE::INTEGER_2:
        case GDSII_DATATYPE::INTEGER_4:
        {
            int32_t i = nodes[ni].u.i;
            ldout << "    " << i << "\n";
            uint32_t vu = (i >= 0) ? i : ((datatype == GDSII_DATATYPE::INTEGER_2) ? (0x10000 + i) : (0x100000000LL + i));
            if ( datatype == GDSII_DATATYPE::INTEGER_4 ) {
                bytes[byte_cnt++] = (vu >> 24) & 0xff;
                bytes[byte_cnt++] = (vu >> 16) & 0xff;
            }
            bytes[byte_cnt++] = (vu >> 8) & 0xff;
            bytes[byte_cnt++] = (vu >> 0) & 0xff;
            break;
        }

        case GDSII_DATATYPE::REAL_4:
        case GDSII_DATATYPE::REAL_8:
        {
            ldout << "    " << nodes[ni].u.r << "\n";
            const uint64_t du = *reinterpret_cast<uint64_t *>(&nodes[ni].u.r);
            uint64_t sign  = (du >> 63) & 1LL;
            int64_t  rexp  = (du >> 52) & 0x7ffLL;
            uint64_t ifrac = du & 0xfffffffffffffLL;
            if ( rexp != 0 ) ifrac |= 1L << 52;     // implied 1. for normalized numbers
            rexp -= (1 << 10)-1;                    // subtract double exp bias
            int datum_byte_cnt = (datatype == GDSII_DATATYPE::REAL_4) ? 4 : 8;
            if ( datum_byte_cnt == 4 ) {
                // fit into 24 bits
                ifrac >>= (53-24);
                rexp  += 53-24;
            } else {
                // fit into 56 bits
                ifrac <<= (56-53);
                rexp  -= 56-53;
            }
            ldout << "mid: rexp " << rexp << " ifrac=" << ifrac << "\n";
            while( (rexp & 0x3LL) != 0 )
            {
                rexp++;
                ifrac >>= 1;
            }
            uint8_t exp = (rexp >> 2) + 65;
            bytes[byte_cnt++] = (sign << 7) | exp;
            for( int j = datum_byte_cnt-2; j >= 0; j-- )
            {
                bytes[byte_cnt++] = (ifrac >> (8*j)) & 0xff;
            }
            for( int j = 0; j < datum_byte_cnt; j++ )
            {
                ldout << "bytes[" << j << "]=" << int(bytes[byte_cnt-datum_byte_cnt+j]) << "\n";
            }
            ldout << "sign=" << sign << " exp=" << int(exp) << " rexp=" << rexp << " ifrac=" << ifrac << " r=" << nodes[ni].u.r << "\n";
            break;
        }

        default:
        {
            assert( false );
        }
    }
}

void Layout::gdsii_write_bytes( const uint8_t * bytes, uint byte_cnt )
{
    while( byte_cnt != 0 )
    {
        uint32_t this_byte_cnt = byte_cnt;
        if ( (gdsii_buff_byte_cnt+this_byte_cnt) > gdsii_buff_alloc_byte_cnt ) {
            this_byte_cnt -= (gdsii_buff_byte_cnt+this_byte_cnt) - gdsii_buff_alloc_byte_cnt;
        }
        memcpy( &gdsii_buff[gdsii_buff_byte_cnt], bytes, this_byte_cnt );
        bytes               += this_byte_cnt;
        byte_cnt            -= this_byte_cnt;
        gdsii_buff_byte_cnt += this_byte_cnt;
        if ( gdsii_buff_byte_cnt == gdsii_buff_alloc_byte_cnt ) {
            gdsii_flush();
        }
    }
}

void Layout::gdsii_flush( void )
{
    if ( gdsii_buff_byte_cnt == 0 ) return;

    if ( ::write( gdsii_fd, gdsii_buff, gdsii_buff_byte_cnt ) <= 0 ) { 
        close( gdsii_fd ); 
        rtnn_assert( false, std::string("could not write() gdsii file - write() error: ") + strerror( errno ) );
    }
    gdsii_buff_byte_cnt = 0;
}

bool Layout::aedt_read( std::string file )
{
    //------------------------------------------------------------
    // Map in file
    //------------------------------------------------------------
    line_num = 1;
    if ( !open_and_read( file, nnn_start, nnn_end ) ) return false;
    nnn = nnn_start;

    //------------------------------------------------------------
    // Parse .aedt file contents
    // Parse the first $begin .. $end.
    //------------------------------------------------------------
    aedt_begin_str_i = str_get( "$begin" );
    aedt_end_str_i   = str_get( "$end" );
    true_str_i       = str_get( "true" );
    false_str_i      = str_get( "false" );

    if ( !aedt_read_expr( hdr->root_i ) ) return false;
    assert( nodes[hdr->root_i].kind == NODE_KIND::HIER );
    return true;
}

bool Layout::aedt_read_expr( uint& ni )
{
    perhaps_realloc( nodes, hdr->node_cnt, max->node_cnt, 1 );
    ni = hdr->node_cnt++;
    nodes[ni].sibling_i = uint(-1);

    skip_whitespace( nnn, nnn_end );
    char ch = *nnn;
    if ( ch == '\'' ) {
        ldout << "STR\n";
        nodes[ni].kind = NODE_KIND::STR;
        if ( !parse_string_i( nodes[ni].u.s_i, nnn, nnn_end ) ) return false;
    } else if ( ch == '-' || (ch >= '0' && ch <= '9') ) {
        ldout << "NUMBER\n";
        return parse_number( ni, nnn, nnn_end );
    } else {
        uint id_i;
        uint8_t * nnn_save = nnn;
        if ( !parse_id( id_i, nnn, nnn_end ) ) {
            rtn_assert( 0, "unable to parse an expression: std::string, number, or id " + surrounding_lines( nnn_save, nnn_end ) );
        }
        ldout << "ID START " << std::string(&strings[id_i]) << "\n";
        if ( id_i == aedt_begin_str_i ) {
            uint name_i;
            if ( !aedt_read_expr( name_i ) ) return false;             // STR node
            rtn_assert( nodes[name_i].kind == NODE_KIND::STR, "$begin not followed by std::string" );
            ldout << "BEGIN " << std::string(&strings[nodes[name_i].u.s_i]) << "\n";

            nodes[ni].kind = NODE_KIND::HIER;
            nodes[ni].u.child_first_i = name_i;
            uint prev_i = name_i;
            for( ;; )
            {
                skip_whitespace( nnn, nnn_end );
                uint id_i;
                if ( peek_id( id_i, nnn, nnn_end ) ) {
                    if ( id_i == aedt_end_str_i ) {
                        parse_id( id_i, nnn, nnn_end );
                        uint end_str_i;
                        if ( !parse_string_i( end_str_i, nnn, nnn_end ) ) return false;
                        ldout << "END " << std::string(&strings[end_str_i]) << "\n";
                        rtn_assert( end_str_i == nodes[name_i].u.s_i, "$end id does not match $begin id " + surrounding_lines( nnn, nnn_end ) );
                        break;
                    }
                }

                uint child_i;
                if ( !aedt_read_expr( child_i ) ) return false;
                nodes[prev_i].sibling_i = child_i;
                prev_i = child_i;
            }
        } else if ( id_i == true_str_i || id_i == false_str_i ) {
            ldout << "BOOL\n";
            nodes[ni].kind = NODE_KIND::BOOL;
            nodes[ni].u.b = id_i == true_str_i;
        } else {
            ldout << "USER ID\n";
            nodes[ni].kind = NODE_KIND::ID;
            nodes[ni].u.s_i = id_i;
        }
    }

    skip_whitespace( nnn, nnn_end );
    ch = *nnn;
    if ( ch == '=' ) {
        ldout << "ASSIGN\n";
        expect_char( ch, nnn, nnn_end );

        uint rhs_i;
        if ( !aedt_read_expr( rhs_i ) ) return false;

        perhaps_realloc( nodes, hdr->node_cnt, max->node_cnt, 1 );
        uint ai = hdr->node_cnt++;
        nodes[ai].kind = NODE_KIND::ASSIGN;
        nodes[ai].u.child_first_i = ni;
        nodes[ni].sibling_i = rhs_i;
        nodes[ai].sibling_i = uint(-1);
        ni = ai;

    } else if ( ch == '(' || ch == '[' ) {
        ldout << std::string( (ch == '(') ? "CALL\n" : "SLICE\n" );
        uint id_i = ni;
        perhaps_realloc( nodes, hdr->node_cnt, max->node_cnt, 1 );
        ni = hdr->node_cnt++;
        nodes[ni].kind = (ch == '(') ? NODE_KIND::CALL : NODE_KIND::SLICE;
        nodes[ni].u.child_first_i = id_i;
        nodes[ni].sibling_i = uint(-1);
        uint prev_i = id_i;
        expect_char( ch, nnn, nnn_end );
        for( bool have_one=false; ; have_one=true )
        {
            skip_whitespace( nnn, nnn_end );
            ch = *nnn;
            if ( (ch == ')' && nodes[ni].kind == NODE_KIND::CALL) || (ch == ']' && nodes[ni].kind == NODE_KIND::SLICE) ) {
                expect_char( ch, nnn, nnn_end );
                break;
            }
            if ( have_one ) {
                skip_whitespace( nnn, nnn_end );
                if ( ch == ':' && nodes[ni].kind == NODE_KIND::SLICE ) {
                    // skip this, it's implied
                    //
                    expect_char( ch, nnn, nnn_end );
                    continue;
                } else if ( ch == ',' ) {
                    expect_char( ',', nnn, nnn_end );
                }
            }

            uint arg_i;
            if ( !aedt_read_expr( arg_i ) ) return false;
            nodes[prev_i].sibling_i = arg_i;
            prev_i = arg_i;
        }
    }

    return true;
}

bool Layout::aedt_write( std::string file )
{
    std::ofstream out( file, std::ofstream::out );
    aedt_write_expr( out, hdr->root_i, "\n" );
    out.close();
    return false;
}

void Layout::aedt_write_expr( std::ofstream& out, uint ni, std::string indent_str )
{
    assert( ni != uint(-1) );
    const Node& node = nodes[ni];
    out << indent_str;
    switch( node.kind ) 
    {
        case NODE_KIND::STR:
            out << "'" << std::string(&strings[node.u.s_i]) << "'";
            break;

        case NODE_KIND::BOOL:
            out << (node.u.b ? "true" : "false");
            break;

        case NODE_KIND::INT:
            out << node.u.i;
            break;
            
        case NODE_KIND::UINT:
            out << node.u.u;
            break;
            
        case NODE_KIND::REAL:
            out << node.u.r;
            break;

        case NODE_KIND::ID:
            out << std::string(&strings[node.u.s_i]);
            break;

        case NODE_KIND::ASSIGN:
        {
            uint child_i = node.u.child_first_i;
            aedt_write_expr( out, child_i, "" );
            out << "=";
            child_i = nodes[child_i].sibling_i;
            aedt_write_expr( out, child_i, "" );
            break;
        }

        case NODE_KIND::CALL:
        {
            uint child_i = node.u.child_first_i;
            aedt_write_expr( out, child_i, "" );
            out << "(";
            bool have_one = false;
            for( child_i = nodes[child_i].sibling_i; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
            {
                if ( have_one ) out << ", ";
                aedt_write_expr( out, child_i, "" );
                have_one = true;
            }
            out << ")";
            break;
        }

        case NODE_KIND::SLICE:
        {
            uint child_i = node.u.child_first_i;
            aedt_write_expr( out, child_i, "" );
            out << "[";
            child_i = nodes[child_i].sibling_i;
            if ( child_i != uint(-1) ) {
                aedt_write_expr( out, child_i, "" );
                out << ":";
                bool have_one = false;
                for( child_i = nodes[child_i].sibling_i; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
                {
                    out << (have_one ? ", " : " ");
                    aedt_write_expr( out, child_i, "" );
                    have_one = true;
                }
            }
            out << "]";
            break;
        }

        case NODE_KIND::HIER:
        {
            uint id_i = node.u.child_first_i;
            uint child_i;
            if ( id_i != uint(-1) && nodes[id_i].kind == NODE_KIND::STR ) {
                out << "$begin '" << std::string(&strings[nodes[id_i].u.s_i]) << "'";
                child_i = nodes[id_i].sibling_i;
            } else {
                out << "$begin 'FILE'";
                child_i = id_i;
                id_i = uint(-1);
            }
            for( ; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
            {
                aedt_write_expr( out, child_i, indent_str + "\t" );
            }
            out << indent_str;
            if ( id_i != uint(-1) ) {
                out << "$end '" << std::string(&strings[nodes[id_i].u.s_i]) << "'";
            } else {
                out << "$end 'FILE'\n";
            }
            break;
        }

        default:
        {
            if ( int(node.kind) >= int(NODE_KIND::GDSII_HEADER) ) {
                GDSII_KIND gkind = GDSII_KIND( int(node.kind) - int(NODE_KIND::GDSII_HEADER) );
                GDSII_DATATYPE datatype = kind_to_datatype( gkind );
                if ( gdsii_is_hier( gkind ) ) {
                    out << "$begin '" << node.kind << "'";
                }
                uint child_i = uint(-1);
                switch( datatype )
                {
                    case GDSII_DATATYPE::NO_DATA:
                    {
                        break;
                    }

                    case GDSII_DATATYPE::BITARRAY:
                    {
                        out << node.kind << "(" << node.u.u << ")";
                        break;
                    }

                    case GDSII_DATATYPE::STRING: {
                        out << node.kind << "('" << std::string(&strings[node.u.s_i]) << "')";
                        break;
                    }

                    case GDSII_DATATYPE::INTEGER_2:
                    case GDSII_DATATYPE::INTEGER_4:
                    case GDSII_DATATYPE::REAL_4:
                    case GDSII_DATATYPE::REAL_8:
                    {
                        std::string vals = "";
                        for( child_i = node.u.child_first_i; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
                        {
                            if ( nodes[child_i].kind != NODE_KIND::INT && nodes[child_i].kind != NODE_KIND::REAL ) break;
                            if ( vals != "" ) vals += ", ";
                            vals += (nodes[child_i].kind == NODE_KIND::INT) ? std::to_string(nodes[child_i].u.i) : std::to_string(nodes[child_i].u.r);
                        }
                        if ( gdsii_is_hier( gkind ) ) {
                            out << indent_str << "\tARGS(" << vals << ")";
                        } else {
                            out << node.kind << "(" << vals << ")";
                        }
                        break;
                    }

                    default:
                    {
                        ldout << "ERROR: unknown GDSII_DATATYPE " << int(datatype) << "\n";
                        exit( 1 );
                    }
                }

                if ( gdsii_is_hier( gkind ) ) {
                    if ( child_i == uint(-1) ) child_i = node.u.child_first_i;
                    for( ; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
                    {
                        aedt_write_expr( out, child_i, indent_str + "\t" );
                    }
                    out << indent_str;
                    out << "$end '" << node.kind << "'";
                    break;
                }
            } else {
                ldout << "ERROR: unknown NODE_KIND " << int(node.kind) << "\n";
                exit( 1 );
            }
        }
        break;
    }
}

void Layout::dissect_path( std::string path, std::string& dir_name, std::string& base_name, std::string& ext_name ) 
{
    // in case they are not found:
    dir_name = "";
    base_name = "";
    ext_name = "";

    // ext_name 
    const int len = path.length();
    int pos = len - 1;
    for( ; pos >= 0; pos-- )
    {
        if ( path[pos] == '.' ) {
            ext_name = path.substr( pos );
            pos--;
            break;
        }
    } 
    if ( pos < 0 ) pos = len - 1;  // no ext_name, so reset for base_name
    
    // base_name
    int base_len = 0;
    for( ; pos >= 0; pos--, base_len++ )
    {
        if ( path[pos] == '/' ) {
            if ( base_len != 0 ) base_name = path.substr( pos+1, base_len );
            pos--;
            break;
        }
    }
    if ( pos < 0 ) pos = len - 1;  // no base_name, so reset for dir_name

    // dir_name is whatever's left
    if ( pos >= 0 ) dir_name = path.substr( 0, pos+1 );
}

std::string Layout::path_without_ext( std::string path, std::string * file_ext )
{
    char * path_c = strdup( path.c_str() );
    for( uint i = path.length() - 1; i != 0; i-- )
    {
        if ( path_c[i] == '.' ) {
            if ( file_ext != nullptr ) *file_ext = std::string( &path_c[i] );
            path_c[i] = '\0';
            break;
        }
    }
    return std::string( path_c );
}

bool Layout::cmd( std::string s )
{
    return system( s.c_str() ) == 0;
}

bool Layout::open_and_read( std::string file_path, uint8_t *& start, uint8_t *& end )
{
    const char * fname = file_path.c_str();
    int fd = open( fname, O_RDONLY );
    if ( fd < 0 ) ldout << "open_and_read() error reading " << file_path << ": " << strerror( errno ) << "\n";
    rtn_assert( fd >= 0, "could not open file " + file_path + " - open() error: " + strerror( errno ) );

    struct stat file_stat;
    int status = fstat( fd, &file_stat );
    if ( status < 0 ) {
        close( fd );
        rtn_assert( 0, "could not stat file " + std::string(fname) + " - stat() error: " + strerror( errno ) );
    }
    size_t size = file_stat.st_size;

    // this large read should behave like an mmap() inside the o/s kernel and be as fast
    start = aligned_alloc<uint8_t>( size );
    if ( start == nullptr ) {
        close( fd );
        rtn_assert( 0, "could not read file " + std::string(fname) + " - malloc() error: " + strerror( errno ) );
    }
    end = start + size;

    uint8_t * addr = start;
    while( size != 0 ) 
    {
        size_t _this_size = 1024*1024*1024;
        if ( size < _this_size ) _this_size = size;
        if ( ::read( fd, addr, _this_size ) <= 0 ) {
            close( fd );
            rtn_assert( 0, "could not read() file " + std::string(fname) + " - read error: " + std::string( strerror( errno ) ) );
        }
        size -= _this_size;
        addr += _this_size;
    }
    close( fd );
    return true;
}

inline void Layout::skip_whitespace( uint8_t *& xxx, uint8_t *& xxx_end )
{
    bool in_comment = false;
    for( ;; )
    {
        if ( xxx == xxx_end ) break;

        char ch = *xxx;
        if ( ch == '#' ) in_comment = true;
        if ( !in_comment && ch != ' ' && ch != '\n' && ch != '\r' && ch != '\t' ) break;

        if ( ch == '\n' || ch == '\r' ) {
            if ( ch == '\n' && xxx == nnn ) line_num++;
            in_comment = false;
        }
        xxx++;
    }
}

inline void Layout::skip_whitespace_to_eol( uint8_t *& xxx, uint8_t *& xxx_end )
{
    bool in_comment = false;
    for( ;; )
    {
        if ( xxx == xxx_end ) break;

        char ch = *xxx;
        if ( ch == '#' ) in_comment = true;
        if ( !in_comment && ch != ' ' && ch != '\n' && ch != '\r' && ch != '\t' ) break;

        if ( ch == '\n' || ch == '\r' ) {
            if ( ch == '\n' && xxx == nnn ) line_num++;
            break;
        }
        xxx++;
    }
}

inline void Layout::skip_to_eol( uint8_t *& xxx, uint8_t *& xxx_end )
{
    if ( !eol( xxx, xxx_end ) ) {
        while( xxx != xxx_end )
        {
            char ch = *xxx;
            if ( ch == '\n' || ch == '\r' ) break;
            xxx++;
        }
    }
}

inline bool Layout::eol( uint8_t *& xxx, uint8_t *& xxx_end )
{
    skip_whitespace_to_eol( xxx, xxx_end );

    if ( xxx == xxx_end || *xxx == '\n' || *xxx == '\r' ) {
        if ( xxx != xxx_end ) {
            if ( *xxx == '\n' && xxx == nnn ) line_num++;
            xxx++;
        }
        ldout << "at eol\n";
        return true;
    } else {
        ldout << "not at eol, char='" + std::string( 1, *xxx ) + "'\n";
        return false;
    }
}

inline bool Layout::expect_char( char ch, uint8_t *& xxx, uint8_t * xxx_end, bool skip_whitespace_first )
{
    if ( skip_whitespace_first ) skip_whitespace( xxx, xxx_end );
    rtn_assert( xxx != xxx_end, "premature end of file" );
    rtn_assert( *xxx == ch, "expected character '" + std::string(1, ch) + "' got '" + std::string( 1, *xxx ) + "' " + surrounding_lines( xxx, xxx_end ) );
    xxx++;
    return true;
}

inline uint Layout::str_get( std::string s )
{
    auto it = str_to_str_i.find( s );
    if ( it != str_to_str_i.end() ) return it->second;
        
    uint s_len = s.length();
    perhaps_realloc( strings, hdr->char_cnt, max->char_cnt, s_len+1 );
    uint s_i = hdr->char_cnt;
    str_to_str_i[s] = s_i;
    char * to_s = &strings[s_i];
    hdr->char_cnt += s_len + 1;
    memcpy( to_s, s.c_str(), s_len+1 );
    ldout << "str_i[" + s + "]=" + std::to_string(s_i) + " strings[]=" + std::string(&strings[s_i]) << "\n";
    return s_i;
}

inline bool Layout::parse_number( uint node_i, uint8_t *& xxx, uint8_t *& xxx_end )
{
    uint8_t * xxx_orig = xxx;

    _int i;
    if ( !parse_int( i, xxx, xxx_end ) ) return false;
    char ch = *xxx;
    if ( ch != '.' && ch != 'e' && ch != 'E' ) {
        if ( i < 0 ) {
            nodes[node_i].kind = NODE_KIND::INT;
            nodes[node_i].u.i  = i;
        } else {
            nodes[node_i].kind = NODE_KIND::UINT;
            nodes[node_i].u.u  = i;
        }
        return true;
    } else {
        xxx = xxx_orig;
        nodes[node_i].kind = NODE_KIND::REAL;
        return parse_real( nodes[node_i].u.r, xxx, xxx_end );
    }
}

inline bool Layout::parse_string( std::string& s, uint8_t *& xxx, uint8_t *& xxx_end )
{
    if ( !expect_char( '\'', xxx, xxx_end, true ) ) return false;
    s = "";
    for( ;; ) 
    {
        rtn_assert( xxx != xxx_end, "no terminating \" for std::string" );
        if ( *xxx == '\'' ) {
            xxx++;
            return true;
        }
        s += *xxx;
        xxx++;
    }
}

inline bool Layout::parse_string_i( uint& s_i, uint8_t *& xxx, uint8_t *& xxx_end )
{
    std::string s;
    if ( !parse_string( s, xxx, xxx_end ) ) return false;
    s_i = str_get( s );
    return true;
}

inline bool Layout::parse_id( uint& id_i, uint8_t *& xxx, uint8_t *& xxx_end )
{
    skip_whitespace( xxx, xxx_end );

    std::string id = "";
    while( xxx != xxx_end )
    {
        char ch = *xxx;
        if ( !( ch == '$' || ch == '_' || (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') || (id != "" && ch >= '0' && ch <= '9')) ) break;

        id += std::string( 1, ch );
        xxx++;
    }
    if ( id == "" ) return false;
    id_i = str_get( id );
    return true;
}

inline bool Layout::peek_id( uint& id_i, uint8_t *& xxx_orig, uint8_t *& xxx_end )
{
    uint8_t * xxx = xxx_orig;
    return parse_id( id_i, xxx, xxx_end );
}

inline bool Layout::parse_real( Layout::real& r, uint8_t *& xxx, uint8_t *& xxx_end )
{
    skip_whitespace( xxx, xxx_end );   
    std::string s = "";
    bool in_frac = false;
    bool has_exp = false;
    while( xxx != xxx_end )
    {
        char ch = *xxx;

        if ( ch == 'n' || ch == 'N' ) {
            // better be a NaN
            xxx++;
            if ( xxx == xxx_end || (*xxx != 'a' && *xxx != 'A') ) return false;
            xxx++;
            if ( xxx == xxx_end || (*xxx != 'n' && *xxx != 'N') ) return false;
            xxx++;
            //r = std::nan( "1" );
            r = 0.0;                    // make them zeros
            return true;
        }

        if ( ch == '-' && !in_frac ) {
            s += "-";
            xxx++;
            continue;
        }

        if ( ch == '.' && !in_frac ) {
            s += ".";
            in_frac = true;
            xxx++;
            continue;
        }

        if ( ch == 'e' || ch == 'E' ) {
            rtn_assert( !has_exp, "real has more than one 'e' exponent" );
            has_exp = true;
            s += std::string( 1, ch );
            xxx++;
            _int e10;
            if ( !parse_int( e10, xxx, xxx_end ) ) return false;
            s += std::to_string( e10 );
            continue;
        }

        if ( ch < '0' || ch > '9' ) break;

        s += std::string( 1, ch );
        xxx++;
    }

    rtn_assert( s.length() != 0, "unable to parse real in " + ext_name + " file " + surrounding_lines( xxx, xxx_end ) );

    r = std::atof( s.c_str() );
    ldout << "real=" + std::to_string( r ) << "\n";
    return true;
}

inline bool Layout::parse_int( _int& i, uint8_t *& xxx, uint8_t *& xxx_end )
{
    bool vld = false;
    i = 0;
    while( xxx != xxx_end && (*xxx == ' ' || *xxx == '\t') ) xxx++;  // skip leading spaces

    bool is_neg = false;
    while( xxx != xxx_end )
    {
        char ch = *xxx;
        if ( ch == '+' ) {
            rtn_assert( !is_neg, "-+ not allowed for an int" );
            xxx++;
            continue;
        }    
        if ( ch == '-' ) {
            rtn_assert( !is_neg, "too many minus signs" );
            is_neg = true;
            xxx++;
            continue;
        }

        if ( ch < '0' || ch > '9' ) break;
        xxx++;

        i = i*10 + (ch - '0');
        vld = true;
    }

    if ( is_neg ) i = -i;
    rtn_assert( vld, "unable to parse int in " + ext_name + " file " + surrounding_lines( xxx, xxx_end ) );
    return true;
}

inline bool Layout::parse_uint( uint& u, uint8_t *& xxx, uint8_t *& xxx_end )
{
    _int i;
    if ( !parse_int( i, xxx, xxx_end ) ) return false;
    rtn_assert( i >= 0, "parse_uint encountered negative integer" );
    u = i;
    return true;
}

inline bool Layout::parse_bool( bool& b, uint8_t *& xxx, uint8_t *& xxx_end )
{
    uint id_i;
    if ( !parse_id( id_i, xxx, xxx_end ) ) return false;
    b = id_i == true_str_i;
    return b || id_i == false_str_i;
}

std::string Layout::surrounding_lines( uint8_t *& xxx, uint8_t *& xxx_end )
{
    uint eol_cnt = 0;
    std::string s = "line " + std::to_string(line_num) + " ";
    while( eol_cnt != 10 && xxx != xxx_end )
    {
        s += std::string( 1, *xxx ) ;
        if ( *xxx == '\n' ) eol_cnt++;
        xxx++;
    }
    return s;
}

#endif
