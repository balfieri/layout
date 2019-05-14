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

    Layout( std::string top_file );
    ~Layout(); 

    bool write( std::string file_path );

    static void        dissect_path( std::string path, std::string& dir_name, std::string& base_name, std::string& ext_name ); // utility
    static std::string path_without_ext( std::string path, std::string * file_ext=nullptr );  // returns path without file extension part

    static const uint VERSION = 0xB0BA1f01; // current version 

    bool                is_good;            // set to true if constructor succeeds
    std::string         error_msg;          // if !is_good

    class real3
    {
    public:
        real c[3];
        
        real3( void ) {}
        real3( real c0, real c1, real c2 ) { c[0] = c0; c[1] = c1; c[2] = c2; }

        real   dot( const real3 &v2 ) const;
        real3  cross( const real3 &v2 ) const;
        real   length( void ) const;
        real   length_sqr( void ) const ;
        real3& normalize( void );
        real3  normalized( void ) const;
        real3  operator + ( const real3& v ) const;
        real3  operator - ( const real3& v ) const;
        real3  operator * ( const real3& v ) const;
        real3  operator * ( real s ) const;
        real3  operator / ( const real3& v ) const;
        real3  operator / ( real s ) const;
        real3& operator += ( const real3 &v2 );
        real3& operator -= ( const real3 &v2 );
        real3& operator *= ( const real3 &v2 );
        real3& operator *= ( const real s );
        real3& operator /= ( const real3 &v2 );
        real3& operator /= ( const real s );
    };

    class AABB                              // axis aligned bounding box
    {
    public:
        real3           min;                // bounding box min
        real3           max;                // bounding box max

        AABB( void ) {}
        AABB( const real3& p );             // init with one point
        AABB( const real3& p0, const real3& p1, const real3& p2 );

        void pad( real p );
        void expand( const AABB& other );
        void expand( const real3& p );
        bool encloses( const AABB& other ) const;
        bool hit( const real3& origin, const real3& direction, const real3& direction_inv, real tmin, real tmax ) const; 
    };

    class Header                            // header (of future binary file)
    {
    public:
        uint        version;                // version
        uint        byte_cnt;               // total in-memory bytes including this header
        uint        char_cnt;               // in std::string array

        uint        node_cnt;               // in nodes array  
        uint        root_i;                 // index of root node in nodes array
    };

    enum class NODE_KIND
    {
        STR,                                // scalars
        BOOL,
        INT,
        UINT,
        REAL,
        ID,

        CALL,                               // child 0 is id, other children are args
        SLICE,                              // child 0 is id, child 1 is index before the ':', other children are other args
        ASSIGN,                             // child 0 is lhs, child 1 is rhs
        HIER,                               // child 0 is id, other children are normal children
    };
            
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

    enum class INSTANCE_KIND
    {
        LAYOUT,                             // instance is another Layout (by name)
        LAYOUT_PTR,                         // instance is another Layout (read in and with a resolved pointer)
        NODE,                               // instance is a node within this Layout
    };

    // structs
    char *              mapped_region;      // != nullptr means the whole file was sucked in by read_uncompressed()
    Header *            hdr;
    Header *            max;                // holds max lengths of currently allocated arrays 

    // arrays
    char *              strings;
    Node *              nodes;

    // maps of names to array indexes
    std::map<std::string, uint>         str_to_str_i;          // maps std::string to unique location in strings[]

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
    std::string ext_name;
    char * nnn_start;
    char * nnn_end;
    char * nnn;
    uint   line_num;

    uint aedt_begin_str_i;              // these are to make it easier to compare
    uint aedt_end_str_i;
    uint true_str_i;
    uint false_str_i;

    bool read_layout( std::string file_path );          // .layout
    bool write_layout( std::string file_path );         

    bool read_gdsii( std::string file_path );           // .gds
    bool parse_gdsii_record( uint& node_i );
    bool write_gdsii( std::string file );
    void write_gdsii_record( std::ofstream& out, uint node_i );
    
    bool read_aedt( std::string file );                 // .aedt
    bool parse_aedt_expr( uint& node_i );
    bool write_aedt( std::string file );
    void write_aedt_expr( std::ofstream& out, uint node_i, std::string indent_str );

    bool open_and_read( std::string file_name, char *& start, char *& end );

    void skip_whitespace_to_eol( char *& xxx, char *& xxx_end );  // on this line only
    void skip_whitespace( char *& xxx, char *& xxx_end );
    void skip_to_eol( char *& xxx, char *& xxx_end );
    bool eol( char *& xxx, char *& xxx_end );
    bool expect_char( char ch, char *& xxx, char* xxx_end, bool skip_whitespace_first=false );
    uint get_str_i( std::string s );
    bool parse_number( uint node_i, char *& xxx, char *& xxx_end );
    bool parse_string( std::string& s, char *& xxx, char *& xxx_end );
    bool parse_string_i( uint& s, char *& xxx, char *& xxx_end );
    bool parse_name( char *& name, char *& xxx, char *& xxx_end );
    bool parse_id( uint& id_i, char *& xxx, char *& xxx_end );
    bool peek_id( uint& id_i, char *& xxx, char *& xxx_end );
    bool parse_real3( real3& r3, char *& xxx, char *& xxx_end, bool has_brackets=false );
    bool parse_real( real& r, char *& xxx, char *& xxx_end );
    bool parse_int( _int& i, char *& xxx, char *& xxx_end );
    bool parse_uint( uint& u, char *& xxx, char *& xxx_end );
    bool parse_bool( bool& b, char *& xxx, char *& xxx_end );
    std::string surrounding_lines( char *& xxx, char *& xxx_end );

    // allocates an array of T on a page boundary
    template<typename T>
    T * aligned_alloc( uint   cnt );

    // reallocate array if we are about to exceed its current size
    template<typename T>
    inline void perhaps_realloc( T *& array, const uint  & hdr_cnt, uint  & max_cnt, uint   add_cnt );

    // GDSII
    enum class GDSII_KIND
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
        ENDTEXTN,
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
    static std::string    str( GDSII_KIND kind );
    static std::string    str( GDSII_DATATYPE datatype );
};

#ifdef LAYOUT_DEBUG
#define dprint( msg ) std::cout << (msg) << "\n"
#else
#define dprint( msg )
#endif

// these are done as macros to avoid evaluating msg (it makes a big difference)
#include <assert.h>
#define rtn_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << msg << "\n"; assert( false ); return false; }
#define node_assert( bool, msg ) if ( !(bool) ) { error_msg = std::string(msg); std::cout << msg << "\n"; assert( false ); goto error;   }

inline std::istream& operator >> ( std::istream& is, Layout::real3& v ) 
{
    is >> v.c[0] >> v.c[1] >> v.c[2];
    return is;
}

inline std::ostream& operator << ( std::ostream& os, const Layout::real3& v ) 
{
    os << "[" << v.c[0] << "," << v.c[1] << "," << v.c[2] << "]";
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Layout::AABB& box ) 
{
    os << box.min << ".." << box.max;
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Layout::NODE_KIND& kind ) 
{
    #define ncase( kind ) case Layout::NODE_KIND::kind: os << #kind; break;
    switch( kind )
    {
        ncase( STR )
        ncase( BOOL )
        ncase( INT )
        ncase( UINT )
        ncase( REAL )
        ncase( ID )
        ncase( CALL )
        ncase( SLICE )
        ncase( ASSIGN )
        ncase( HIER )
        default: os << "<unknown>"; break;
    }
    return os;
}

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
        gcase( ENDTEXTN )
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
        kdcase( ENDTEXTN,     INTEGER_4 )
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

    hdr = aligned_alloc<Header>( 1 );
    memset( hdr, 0, sizeof( Header ) );
    hdr->version = VERSION;

    //------------------------------------------------------------
    // Initial lengths of arrays are large in virtual memory
    //------------------------------------------------------------
    max = aligned_alloc<Header>( 1 );
    max->node_cnt =  1024;
    max->char_cnt = max->node_cnt * 128;

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
        if ( !read_layout( top_file ) ) return;
    } else {
        //------------------------------------------------------------
        // Allocate initial arrays
        //------------------------------------------------------------
        strings = aligned_alloc<char>( max->char_cnt );
        nodes   = aligned_alloc<Node>( max->node_cnt );

        if ( ext_name == std::string( ".gds" ) ) {
            if ( !read_gdsii( top_file ) ) return;
        } else if ( ext_name == std::string( ".aedt" ) ) {
            if ( !read_aedt( top_file ) ) return;
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
        return write_layout( top_file );
    } else {
        if ( ext_name == std::string( ".aedt" ) ) {
            return write_aedt( top_file );
        } else {
            error_msg = "unknown file ext_name: " + ext_name;
            return false;
        }
    }
}

// returns array of T on a page boundary
template<typename T>
T * Layout::aligned_alloc( Layout::uint   cnt )
{
    void * mem = nullptr;
    posix_memalign( &mem, getpagesize(), cnt*sizeof(T) );
    return reinterpret_cast<T *>( mem );
}

// reallocate array if we are about to exceed its current size
template<typename T>
inline void Layout::perhaps_realloc( T *& array, const Layout::uint  & hdr_cnt, Layout::uint  & max_cnt, Layout::uint   add_cnt )
{
    while( (hdr_cnt + add_cnt) > max_cnt ) {
        void * mem = nullptr;
        uint   old_max_cnt = max_cnt;
        max_cnt *= 2;
        if ( max_cnt < old_max_cnt ) {
            assert( old_max_cnt != uint(-1) );
            max_cnt = uint(-1);
        }
        posix_memalign( &mem, getpagesize(), max_cnt*sizeof(T) );
        memcpy( mem, array, hdr_cnt*sizeof(T) );
        delete array;
        array = reinterpret_cast<T *>( mem );
    }
}

bool Layout::read_layout( std::string layout_path )
{
    char * start;
    char * end;
    if ( !open_and_read( layout_path, start, end ) ) return false;
    mapped_region = start;

    //------------------------------------------------------------
    // Write out header than individual arrays.
    // Each is padded out to a page boundary in the file.
    //------------------------------------------------------------
    char * _addr = start;
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

    is_good = true;

    return true;
}

bool Layout::write_layout( std::string layout_path ) 
{
    int fd = open( layout_path.c_str(), O_CREAT|O_WRONLY|O_TRUNC|O_SYNC|S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP );
    if ( fd < 0 ) std::cout << "open() for write error: " << strerror( errno ) << "\n";
    rtn_assert( fd >= 0, "could not open() file " + layout_path + " for writing - open() error: " + strerror( errno ) );

    //------------------------------------------------------------
    // Write out header than individual arrays.
    // Each is padded out to a page boundary in the file.
    //------------------------------------------------------------
    size_t page_size = getpagesize();

    #define _uwrite( addr, byte_cnt ) \
    { \
        size_t _byte_cnt = byte_cnt; \
        _byte_cnt += _byte_cnt % page_size; \
        char * _addr = reinterpret_cast<char *>( addr ); \
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

    fsync( fd ); // flush
    close( fd );

    return true;
}

bool Layout::read_aedt( std::string file )
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
    aedt_begin_str_i = get_str_i( "$begin" );
    aedt_end_str_i   = get_str_i( "$end" );
    true_str_i       = get_str_i( "true" );
    false_str_i      = get_str_i( "false" );

    if ( !parse_aedt_expr( hdr->root_i ) ) return false;
    assert( nodes[hdr->root_i].kind == NODE_KIND::HIER );
    return true;
}

bool Layout::parse_aedt_expr( uint& ni )
{
    perhaps_realloc( nodes, hdr->node_cnt, max->node_cnt, 1 );
    ni = hdr->node_cnt++;
    nodes[ni].sibling_i = uint(-1);

    skip_whitespace( nnn, nnn_end );
    char ch = *nnn;
    if ( ch == '\'' ) {
        dprint( "STR" );
        nodes[ni].kind = NODE_KIND::STR;
        if ( !parse_string_i( nodes[ni].u.s_i, nnn, nnn_end ) ) return false;
    } else if ( ch == '-' || (ch >= '0' && ch <= '9') ) {
        dprint( "NUMBER" );
        return parse_number( ni, nnn, nnn_end );
    } else {
        uint id_i;
        char * nnn_save = nnn;
        if ( !parse_id( id_i, nnn, nnn_end ) ) {
            rtn_assert( 0, "unable to parse an expression: std::string, number, or id " + surrounding_lines( nnn_save, nnn_end ) );
        }
        dprint( "ID START " + std::string(&strings[id_i]) );
        if ( id_i == aedt_begin_str_i ) {
            uint name_i;
            if ( !parse_aedt_expr( name_i ) ) return false;             // STR node
            rtn_assert( nodes[name_i].kind == NODE_KIND::STR, "$begin not followed by std::string" );
            dprint( "BEGIN " + std::string(&strings[nodes[name_i].u.s_i]) );

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
                        dprint( "END " + std::string(&strings[end_str_i]) );
                        rtn_assert( end_str_i == nodes[name_i].u.s_i, "$end id does not match $begin id " + surrounding_lines( nnn, nnn_end ) );
                        break;
                    }
                }

                uint child_i;
                if ( !parse_aedt_expr( child_i ) ) return false;
                nodes[prev_i].sibling_i = child_i;
                prev_i = child_i;
            }
        } else if ( id_i == true_str_i || id_i == false_str_i ) {
            dprint( "BOOL" );
            nodes[ni].kind = NODE_KIND::BOOL;
            nodes[ni].u.b = id_i == true_str_i;
        } else {
            dprint( "USER ID" );
            nodes[ni].kind = NODE_KIND::ID;
            nodes[ni].u.s_i = id_i;
        }
    }

    skip_whitespace( nnn, nnn_end );
    ch = *nnn;
    if ( ch == '=' ) {
        dprint( "ASSIGN" );
        expect_char( ch, nnn, nnn_end );

        uint rhs_i;
        if ( !parse_aedt_expr( rhs_i ) ) return false;

        perhaps_realloc( nodes, hdr->node_cnt, max->node_cnt, 1 );
        uint ai = hdr->node_cnt++;
        nodes[ai].kind = NODE_KIND::ASSIGN;
        nodes[ai].u.child_first_i = ni;
        nodes[ni].sibling_i = rhs_i;
        nodes[ai].sibling_i = uint(-1);
        ni = ai;

    } else if ( ch == '(' || ch == '[' ) {
        dprint( (ch == '(') ? "CALL" : "SLICE" );
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
            if ( !parse_aedt_expr( arg_i ) ) return false;
            nodes[prev_i].sibling_i = arg_i;
            prev_i = arg_i;
        }
    }

    return true;
}

bool Layout::write_aedt( std::string file )
{
    std::ofstream out( file, std::ofstream::out );
    write_aedt_expr( out, hdr->root_i, "\n" );
    out.close();
    return false;
}

void Layout::write_aedt_expr( std::ofstream& out, uint ni, std::string indent_str )
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
            write_aedt_expr( out, child_i, "" );
            out << "=";
            child_i = nodes[child_i].sibling_i;
            write_aedt_expr( out, child_i, "" );
            break;
        }

        case NODE_KIND::CALL:
        {
            uint child_i = node.u.child_first_i;
            write_aedt_expr( out, child_i, "" );
            out << "(";
            bool have_one = false;
            for( child_i = nodes[child_i].sibling_i; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
            {
                if ( have_one ) out << ", ";
                write_aedt_expr( out, child_i, "" );
                have_one = true;
            }
            out << ")";
            break;
        }

        case NODE_KIND::SLICE:
        {
            uint child_i = node.u.child_first_i;
            write_aedt_expr( out, child_i, "" );
            out << "[";
            child_i = nodes[child_i].sibling_i;
            if ( child_i != uint(-1) ) {
                write_aedt_expr( out, child_i, "" );
                out << ":";
                bool have_one = false;
                for( child_i = nodes[child_i].sibling_i; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
                {
                    out << (have_one ? ", " : " ");
                    write_aedt_expr( out, child_i, "" );
                    have_one = true;
                }
            }
            out << "]";
            break;
        }

        case NODE_KIND::HIER:
        {
            uint id_i = node.u.child_first_i;
            out << "$begin '" << std::string(&strings[nodes[id_i].u.s_i]) << "'";
            for( uint child_i = nodes[id_i].sibling_i; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
            {
                write_aedt_expr( out, child_i, indent_str + "\t" );
            }
            out << indent_str;
            out << "$end '" << std::string(&strings[nodes[id_i].u.s_i]) << "'";
            break;
        }

        default:
        {
            std::cout << "ERROR: unknown NODE_KIND\n";
            exit( 1 );
            break;
        }
    }
}

bool Layout::read_gdsii( std::string file )
{
    //------------------------------------------------------------
    // Map in file
    //------------------------------------------------------------
    line_num = 1;
    if ( !open_and_read( file, nnn_start, nnn_end ) ) return false;
    nnn = nnn_start;

    if ( !parse_gdsii_record( hdr->root_i ) ) return false;
    assert( nodes[hdr->root_i].kind == NODE_KIND::HIER );
    return true;
}

bool Layout::parse_gdsii_record( uint& ni )
{
    perhaps_realloc( nodes, hdr->node_cnt, max->node_cnt, 1 );
    ni = hdr->node_cnt++;
    nodes[ni].sibling_i = uint(-1);

    skip_whitespace( nnn, nnn_end );
    char ch = *nnn;
    if ( ch == '\'' ) {
        dprint( "STR" );
        nodes[ni].kind = NODE_KIND::STR;
        if ( !parse_string_i( nodes[ni].u.s_i, nnn, nnn_end ) ) return false;
    } else if ( ch == '-' || (ch >= '0' && ch <= '9') ) {
        dprint( "NUMBER" );
        return parse_number( ni, nnn, nnn_end );
    } else {
        uint id_i;
        char * nnn_save = nnn;
        if ( !parse_id( id_i, nnn, nnn_end ) ) {
            rtn_assert( 0, "unable to parse an expression: std::string, number, or id " + surrounding_lines( nnn_save, nnn_end ) );
        }
        dprint( "ID START " + std::string(&strings[id_i]) );
        if ( id_i == 0 ) {
            uint name_i;
            if ( !parse_gdsii_record( name_i ) ) return false;             // STR node
            rtn_assert( nodes[name_i].kind == NODE_KIND::STR, "$begin not followed by std::string" );
            dprint( "BEGIN " + std::string(&strings[nodes[name_i].u.s_i]) );

            nodes[ni].kind = NODE_KIND::HIER;
            nodes[ni].u.child_first_i = name_i;
            uint prev_i = name_i;
            for( ;; )
            {
                skip_whitespace( nnn, nnn_end );
                uint id_i;
                if ( peek_id( id_i, nnn, nnn_end ) ) {
                    if ( id_i == 0 ) {
                        parse_id( id_i, nnn, nnn_end );
                        uint end_str_i;
                        if ( !parse_string_i( end_str_i, nnn, nnn_end ) ) return false;
                        dprint( "END " + std::string(&strings[end_str_i]) );
                        rtn_assert( end_str_i == nodes[name_i].u.s_i, "$end id does not match $begin id " + surrounding_lines( nnn, nnn_end ) );
                        break;
                    }
                }

                uint child_i;
                if ( !parse_gdsii_record( child_i ) ) return false;
                nodes[prev_i].sibling_i = child_i;
                prev_i = child_i;
            }
        } else if ( id_i == true_str_i || id_i == false_str_i ) {
            dprint( "BOOL" );
            nodes[ni].kind = NODE_KIND::BOOL;
            nodes[ni].u.b = id_i == true_str_i;
        } else {
            dprint( "USER ID" );
            nodes[ni].kind = NODE_KIND::ID;
            nodes[ni].u.s_i = id_i;
        }
    }

    skip_whitespace( nnn, nnn_end );
    ch = *nnn;
    if ( ch == '=' ) {
        dprint( "ASSIGN" );
        expect_char( ch, nnn, nnn_end );

        uint rhs_i;
        if ( !parse_gdsii_record( rhs_i ) ) return false;

        perhaps_realloc( nodes, hdr->node_cnt, max->node_cnt, 1 );
        uint ai = hdr->node_cnt++;
        nodes[ai].kind = NODE_KIND::ASSIGN;
        nodes[ai].u.child_first_i = ni;
        nodes[ni].sibling_i = rhs_i;
        nodes[ai].sibling_i = uint(-1);
        ni = ai;

    } else if ( ch == '(' || ch == '[' ) {
        dprint( (ch == '(') ? "CALL" : "SLICE" );
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
            if ( !parse_gdsii_record( arg_i ) ) return false;
            nodes[prev_i].sibling_i = arg_i;
            prev_i = arg_i;
        }
    }

    return true;
}

bool Layout::write_gdsii( std::string file )
{
    std::ofstream out( file, std::ofstream::out );
    write_gdsii_record( out, hdr->root_i );
    out.close();
    return false;
}

void Layout::write_gdsii_record( std::ofstream& out, uint ni )
{
    assert( ni != uint(-1) );
    const Node& node = nodes[ni];
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
            write_gdsii_record( out, child_i );
            out << "=";
            child_i = nodes[child_i].sibling_i;
            write_gdsii_record( out, child_i );
            break;
        }

        case NODE_KIND::CALL:
        {
            uint child_i = node.u.child_first_i;
            write_gdsii_record( out, child_i );
            out << "(";
            bool have_one = false;
            for( child_i = nodes[child_i].sibling_i; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
            {
                if ( have_one ) out << ", ";
                write_gdsii_record( out, child_i );
                have_one = true;
            }
            out << ")";
            break;
        }

        case NODE_KIND::SLICE:
        {
            uint child_i = node.u.child_first_i;
            write_gdsii_record( out, child_i );
            out << "[";
            child_i = nodes[child_i].sibling_i;
            if ( child_i != uint(-1) ) {
                write_gdsii_record( out, child_i );
                out << ":";
                bool have_one = false;
                for( child_i = nodes[child_i].sibling_i; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
                {
                    out << (have_one ? ", " : " ");
                    write_gdsii_record( out, child_i );
                    have_one = true;
                }
            }
            out << "]";
            break;
        }

        case NODE_KIND::HIER:
        {
            uint id_i = node.u.child_first_i;
            out << "$begin '" << std::string(&strings[nodes[id_i].u.s_i]) << "'";
            for( uint child_i = nodes[id_i].sibling_i; child_i != uint(-1); child_i = nodes[child_i].sibling_i )
            {
                write_gdsii_record( out, child_i );
            }
            out << "$end '" << std::string(&strings[nodes[id_i].u.s_i]) << "'";
            break;
        }

        default:
        {
            std::cout << "ERROR: unknown NODE_KIND\n";
            exit( 1 );
            break;
        }
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

bool Layout::open_and_read( std::string file_path, char *& start, char *& end )
{
    const char * fname = file_path.c_str();
    int fd = open( fname, O_RDONLY );
    if ( fd < 0 ) std::cout << "open_and_read() error reading " << file_path << ": " << strerror( errno ) << "\n";
    rtn_assert( fd >= 0, "could not open file " + file_path + " - open() error: " + strerror( errno ) );

    struct stat file_stat;
    int status = fstat( fd, &file_stat );
    if ( status < 0 ) {
        close( fd );
        rtn_assert( 0, "could not stat file " + std::string(fname) + " - stat() error: " + strerror( errno ) );
    }
    size_t size = file_stat.st_size;

    // this large read should behave like an mmap() inside the o/s kernel and be as fast
    start = aligned_alloc<char>( size );
    if ( start == nullptr ) {
        close( fd );
        rtn_assert( 0, "could not read file " + std::string(fname) + " - malloc() error: " + strerror( errno ) );
    }
    end = start + size;

    char * addr = start;
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

inline void Layout::skip_whitespace( char *& xxx, char *& xxx_end )
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

inline void Layout::skip_whitespace_to_eol( char *& xxx, char *& xxx_end )
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

inline void Layout::skip_to_eol( char *& xxx, char *& xxx_end )
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

inline bool Layout::eol( char *& xxx, char *& xxx_end )
{
    skip_whitespace_to_eol( xxx, xxx_end );

    if ( xxx == xxx_end || *xxx == '\n' || *xxx == '\r' ) {
        if ( xxx != xxx_end ) {
            if ( *xxx == '\n' && xxx == nnn ) line_num++;
            xxx++;
        }
        dprint( "at eol" );
        return true;
    } else {
        dprint( "not at eol, char='" + std::string( 1, *xxx ) + "'" );
        return false;
    }
}

inline bool Layout::expect_char( char ch, char *& xxx, char* xxx_end, bool skip_whitespace_first )
{
    if ( skip_whitespace_first ) skip_whitespace( xxx, xxx_end );
    rtn_assert( xxx != xxx_end, "premature end of file" );
    rtn_assert( *xxx == ch, "expected character '" + std::string(1, ch) + "' got '" + std::string( 1, *xxx ) + "' " + surrounding_lines( xxx, xxx_end ) );
    xxx++;
    return true;
}

inline uint Layout::get_str_i( std::string s )
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
    dprint( "str_i[" + s + "]=" + std::to_string(s_i) + " strings[]=" + std::string(&strings[s_i]) );
    return s_i;
}

inline bool Layout::parse_number( uint node_i, char *& xxx, char *& xxx_end )
{
    char * xxx_orig = xxx;

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

inline bool Layout::parse_string( std::string& s, char *& xxx, char *& xxx_end )
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

inline bool Layout::parse_string_i( uint& s_i, char *& xxx, char *& xxx_end )
{
    std::string s;
    if ( !parse_string( s, xxx, xxx_end ) ) return false;
    s_i = get_str_i( s );
    return true;
}

inline bool Layout::parse_name( char *& name, char *& xxx, char *& xxx_end )
{
    bool vld = false;
    perhaps_realloc( strings, hdr->char_cnt, max->char_cnt, 1024 );
    name = &strings[hdr->char_cnt];

    while( xxx != xxx_end && (*xxx == ' ' || *xxx == '\t') ) xxx++;  // skip leading spaces

    uint len = 0;
    while( xxx != xxx_end )
    {
        char ch = *xxx;
        if ( ch == '\n' || ch == '\r' ) break;

        rtn_assert( len < 1024, "string is larger than 1024 characters" );
        name[len++] = ch;
        vld = true;
        xxx++;
    }

    if ( vld ) {
        name[len] = '\0';
        hdr->char_cnt += len+1;
        char * ptr;
        for( ptr = &name[len-1]; ptr != name; ptr-- )
        {
            // skip trailing spaces
            if ( *ptr != ' ' && *ptr != '\t' ) break;
            *ptr = '\0';
        }

        return *ptr != '\0';
    }

    rtn_assert( 0, "could not parse name: " + surrounding_lines( xxx, xxx_end ) );
}

inline bool Layout::parse_id( uint& id_i, char *& xxx, char *& xxx_end )
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
    id_i = get_str_i( id );
    return true;
}

inline bool Layout::peek_id( uint& id_i, char *& xxx_orig, char *& xxx_end )
{
    char * xxx = xxx_orig;
    return parse_id( id_i, xxx, xxx_end );
}

inline bool Layout::parse_real3( Layout::real3& r3, char *& xxx, char *& xxx_end, bool has_brackets )
{
    return (!has_brackets || expect_char( '[', xxx, xxx_end, true )) &&
           parse_real( r3.c[0], xxx, xxx_end ) && 
           (!has_brackets || expect_char( ',', xxx, xxx_end, true )) &&
           parse_real( r3.c[1], xxx, xxx_end ) && 
           (!has_brackets || expect_char( ',', xxx, xxx_end, true )) &&
           parse_real( r3.c[2], xxx, xxx_end ) &&
           (!has_brackets || expect_char( ']', xxx, xxx_end, true ));
}

inline bool Layout::parse_real( Layout::real& r, char *& xxx, char *& xxx_end )
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
    dprint( "real=" + std::to_string( r ) );
    return true;
}

inline bool Layout::parse_int( _int& i, char *& xxx, char *& xxx_end )
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

inline bool Layout::parse_uint( uint& u, char *& xxx, char *& xxx_end )
{
    _int i;
    if ( !parse_int( i, xxx, xxx_end ) ) return false;
    rtn_assert( i >= 0, "parse_uint encountered negative integer" );
    u = i;
    return true;
}

inline bool Layout::parse_bool( bool& b, char *& xxx, char *& xxx_end )
{
    uint id_i;
    if ( !parse_id( id_i, xxx, xxx_end ) ) return false;
    b = id_i == true_str_i;
    return b || id_i == false_str_i;
}

std::string Layout::surrounding_lines( char *& xxx, char *& xxx_end )
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

inline Layout::AABB::AABB( const Layout::real3& p )
{
    min = p;
    max = p;
}  

inline Layout::AABB::AABB( const Layout::real3& p0, const Layout::real3& p1, const Layout::real3& p2 ) 
{
    min = p0;
    max = p0;
    expand( p1 );
    expand( p2 );
}  

inline Layout::real Layout::real3::dot( const Layout::real3 &v2 ) const
{
    return c[0] * v2.c[0] + c[1] * v2.c[1] + c[2] * v2.c[2];
}

inline Layout::real3 Layout::real3::cross( const Layout::real3 &v2 ) const
{
    return real3( (c[1]*v2.c[2]   - c[2]*v2.c[1]),
                  (-(c[0]*v2.c[2] - c[2]*v2.c[0])),
                  (c[0]*v2.c[1]   - c[1]*v2.c[0]) );
}

inline Layout::real Layout::real3::length( void ) const
{ 
    return std::sqrt( c[0]*c[0] + c[1]*c[1] + c[2]*c[2] ); 
}

inline Layout::real Layout::real3::length_sqr( void ) const 
{ 
    return c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
}

inline Layout::real3& Layout::real3::normalize( void )
{
    *this /= length();
    return *this;
}

inline Layout::real3 Layout::real3::normalized( void ) const
{
    return *this / length();
}

inline Layout::real3 Layout::real3::operator + ( const Layout::real3& v2 ) const
{
    real3 r;
    r.c[0] = c[0] + v2.c[0];
    r.c[1] = c[1] + v2.c[1];
    r.c[2] = c[2] + v2.c[2];
    return r;
}

inline Layout::real3 Layout::real3::operator - ( const Layout::real3& v2 ) const
{
    real3 r;
    r.c[0] = c[0] - v2.c[0];
    r.c[1] = c[1] - v2.c[1];
    r.c[2] = c[2] - v2.c[2];
    return r;
}

inline Layout::real3 Layout::real3::operator * ( const Layout::real3& v2 ) const
{
    real3 r;
    r.c[0] = c[0] * v2.c[0];
    r.c[1] = c[1] * v2.c[1];
    r.c[2] = c[2] * v2.c[2];
    return r;
}

inline Layout::real3 operator * ( Layout::real s, const Layout::real3& v ) 
{
    return Layout::real3( s*v.c[0], s*v.c[1], s*v.c[2] );
}

inline Layout::real3 Layout::real3::operator * ( Layout::real s ) const
{
    real3 r;
    r.c[0] = c[0] * s;
    r.c[1] = c[1] * s;
    r.c[2] = c[2] * s;
    return r;
}

inline Layout::real3 Layout::real3::operator / ( const Layout::real3& v2 ) const
{
    real3 r;
    r.c[0] = c[0] / v2.c[0];
    r.c[1] = c[1] / v2.c[1];
    r.c[2] = c[2] / v2.c[2];
    return r;
}

inline Layout::real3 Layout::real3::operator / ( Layout::real s ) const
{
    real3 r;
    r.c[0] = c[0] / s;
    r.c[1] = c[1] / s;
    r.c[2] = c[2] / s;
    return r;
}

inline Layout::real3& Layout::real3::operator += ( const Layout::real3 &v2 )
{
    c[0] += v2.c[0];
    c[1] += v2.c[1];
    c[2] += v2.c[2];
    return *this;
}

inline Layout::real3& Layout::real3::operator -= ( const Layout::real3 &v2 )
{
    c[0] -= v2.c[0];
    c[1] -= v2.c[1];
    c[2] -= v2.c[2];
    return *this;
}

inline Layout::real3& Layout::real3::operator *= ( const Layout::real3 &v2 )
{
    c[0] *= v2.c[0];
    c[1] *= v2.c[1];
    c[2] *= v2.c[2];
    return *this;
}

inline Layout::real3& Layout::real3::operator *= ( const Layout::real s )
{
    c[0] *= s;
    c[1] *= s;
    c[2] *= s;
    return *this;
}

inline Layout::real3& Layout::real3::operator /= ( const Layout::real3 &v2 )
{
    c[0] /= v2.c[0];
    c[1] /= v2.c[1];
    c[2] /= v2.c[2];
    return *this;
}

inline Layout::real3& Layout::real3::operator /= ( const Layout::real s )
{
    c[0] /= s;
    c[1] /= s;
    c[2] /= s;
    return *this;
}

inline void Layout::AABB::pad( Layout::real p ) 
{
    min -= real3( p, p, p );
    max += real3( p, p, p );
}

inline void Layout::AABB::expand( const Layout::AABB& other )
{
    for( uint i = 0; i < 3; i++ )
    {
        if ( other.min.c[i] < min.c[i] ) min.c[i] = other.min.c[i];
        if ( other.max.c[i] > max.c[i] ) max.c[i] = other.max.c[i];
    }
}

inline void Layout::AABB::expand( const Layout::real3& p ) 
{
    if ( p.c[0] < min.c[0] ) min.c[0] = p.c[0];
    if ( p.c[1] < min.c[1] ) min.c[1] = p.c[1];
    if ( p.c[2] < min.c[2] ) min.c[2] = p.c[2];
    if ( p.c[0] > max.c[0] ) max.c[0] = p.c[0];
    if ( p.c[1] > max.c[1] ) max.c[1] = p.c[1];
    if ( p.c[2] > max.c[2] ) max.c[2] = p.c[2];
}

inline bool Layout::AABB::encloses( const AABB& other ) const
{
    return min.c[0] <= other.min.c[0] &&
           min.c[1] <= other.min.c[1] &&
           min.c[2] <= other.min.c[2] &&
           max.c[0] >= other.max.c[0] &&
           max.c[1] >= other.max.c[1] &&
           max.c[2] >= other.max.c[2];
}

inline bool Layout::AABB::hit( const Layout::real3& origin, const Layout::real3& direction, const Layout::real3& direction_inv, 
                              Layout::real tmin, Layout::real tmax ) const 
{
    (void)direction;
    for( uint a = 0; a < 3; a++ ) 
    {
        real dir_inv = direction_inv.c[a];
        real v0 = (min.c[a] - origin.c[a]) * dir_inv;
        real v1 = (max.c[a] - origin.c[a]) * dir_inv;
        tmin = std::fmax( tmin, std::fmin( v0, v1 ) );
        tmax = std::fmin( tmax, std::fmax( v0, v1 ) );
    }
    bool r = tmax >= std::fmax( tmin, real(0.0) );
    return r;
}

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
/***************************************************************/
/* convenient shorthand typedefs *******************************/
/***************************************************************/
#ifndef iVec
  typedef std::vector<int>    iVec;
#endif
#ifndef dVec
  typedef std::vector<double> dVec;
#endif
#ifndef bVec
  typedef std::vector<bool> bVec;
#endif
#ifndef sVec
  typedef std::vector<char *> sVec;
#endif
#ifndef strVec
  typedef std::vector<std::string> strVec;
#endif

/****************************************************************************************/
/* A PolygonList is a collection of polygons living in the XY plane.                    */
/* PolygonList.size() is the number of polygons in the list.                            */
/* PolygonList[np].size()/2 is the number of vertices in polygon #np.                   */
/* PolygonList[np][2*nv+0, 2*nv+1] are the x,y coordinates of vertex #nv in polygon #np.*/
/****************************************************************************************/
typedef std::vector<dVec> PolygonList;

typedef struct { char *Text; dVec XY; int Layer; } TextString;
typedef std::vector<TextString> TextStringList;

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
/*      label) or a text std::string (with a single vertex as       */ 
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
   std::vector<GDSIIElement *> Elements;
   bool IsPCell;
   bool IsReferenced;
   std::string *Name;

 } GDSIIStruct;

typedef struct Entity
 { char *Text;   // if NULL, the entity is a polygon; otherwise it is a text std::string
   dVec XY;      // vertex coordinates: 2 for a text std::string, 2N for an N-gon
   bool Closed;  // true if there exists an edge connecting the last to the first vertex
   char *Label;  // optional descriptive text, may be present or absent for polygons and texts
 } Entity;

typedef std::vector<Entity>     EntityList;
typedef std::vector<EntityList> EntityTable;

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
    void ReadGDSIIFile(const std::string FileName, double CoordinateLengthUnit=0.0);
    int GetStructByName(std::string Name);
    void Flatten(double CoordinateLengthUnit=0.0);

    // general info on the GDSII file
    std::string *LibName;
    std::string *GDSIIFileName;
    double FileUnits[2], UnitInMeters;
    std::set<int> LayerSet; 
    iVec Layers;

    // list of structures (hierarchical, i.e. pre-flattening)
    std::vector<GDSIIStruct *> Structs;

    // table of entities (flattened)
    EntityTable ETable; // ETable[nl][ne] = #neth entity on layer Layers[nl]

    static bool Verbose;
    static char *LogFileName;
    static void Log(const char *format, ...);
    static void ErrExit(const char *format, ...);
    static void Warn(const char *format, ...);
    static char *vstrappend(char *s, const char *format, ...);
    static char *vstrdup(const char *format, ...);
};

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
   std::string *sVal;
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

typedef std::string *(*RecordHandler)(GDSIIRecord Record, ParseState *PState);

const char *ElTypeNames[]=
 {"BOUNDARY", "PATH", "SREF", "AREF", "TEXT", "NODE", "BOX"};

/***************************************************************/
/* Handlers for specific types of data records in GDSII files. */
/***************************************************************/
std::string *handleHEADER(GDSIIRecord Record, ParseState *PState)
{
  (void) Record;
  if (PState->Status!=ParseState::INITIAL)
   return new std::string("unexpected record before HEADER");
  PState->Status=ParseState::INHEADER;
  return 0;
}

std::string *handleBGNLIB(GDSIIRecord Record, ParseState *PState)
{
  (void) Record;
  if (PState->Status!=ParseState::INHEADER)
   return new std::string("unexpected record BGNLIB");
  PState->Status=ParseState::INLIB;
  return 0;
}

std::string *handleLIBNAME(GDSIIRecord Record, ParseState *PState)
{ 
  if (PState->Status!=ParseState::INLIB)
   return new std::string("unexpected record LIBNAME");
  PState->Data->LibName = new std::string( *(Record.sVal) );
  return 0;
}

std::string *handleUNITS(GDSIIRecord Record, ParseState *PState)
{ 
  PState->Data->FileUnits[0] = Record.dVal[0];
  PState->Data->FileUnits[1] = Record.dVal[1];
  PState->Data->UnitInMeters =
   PState->Data->FileUnits[1] / PState->Data->FileUnits[0];
  return 0;
}

std::string *handleENDLIB(GDSIIRecord Record, ParseState *PState)
{ 
  (void) Record;
  if (PState->Status!=ParseState::INLIB)
   return new std::string("unexpected record ENDLIB");
  PState->Status=ParseState::DONE;
  return 0;
}

std::string *handleBGNSTR(GDSIIRecord Record, ParseState *PState)
{
  (void) Record;
  if (PState->Status!=ParseState::INLIB)
   return new std::string("unexpected record BGNSTR");

  // add a new structure
  GDSIIStruct *s  = new GDSIIStruct;
  s->IsReferenced = false;
  s->IsPCell      = false;
  PState->CurrentStruct = s;
  PState->Data->Structs.push_back(s);

  PState->Status=ParseState::INSTRUCT;

  return 0;
}

std::string *handleSTRNAME(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INSTRUCT)
   return new std::string("unexpected record STRNAME");
  PState->CurrentStruct->Name = new std::string( *(Record.sVal) );
  if( strcasestr( Record.sVal->c_str(), "CONTEXT_INFO") )
   PState->CurrentStruct->IsPCell=true;
  return 0;
}

std::string *handleENDSTR(GDSIIRecord Record, ParseState *PState)
{
  (void) Record;
  if (PState->Status!=ParseState::INSTRUCT)
   return new std::string("unexpected record ENDSTR");
  PState->Status=ParseState::INLIB;
  return 0;
}

std::string *handleElement(GDSIIRecord Record, ParseState *PState, ElementType ElType)
{
  (void) Record;
  if (PState->Status!=ParseState::INSTRUCT)
   return new std::string(   std::string("unexpected record") + ElTypeNames[ElType] );
  
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

std::string *handleBOUNDARY(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, BOUNDARY); }

std::string *handlePATH(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, PATH); }

std::string *handleSREF(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, SREF); }

std::string *handleAREF(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, AREF); }

std::string *handleTEXT(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, TEXT); }

std::string *handleNODE(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, NODE); }

std::string *handleBOX(GDSIIRecord Record, ParseState *PState)
{ return handleElement(Record, PState, BOX); }

std::string *handleLAYER(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new std::string("unexpected record LAYER");
  PState->CurrentElement->Layer = Record.iVal[0];
  PState->Data->LayerSet.insert(Record.iVal[0]);
  
  return 0;
}

std::string *handleDATATYPE(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new std::string("unexpected record DATATYPE");
  PState->CurrentElement->DataType = Record.iVal[0];
  return 0;
}

std::string *handleTEXTTYPE(GDSIIRecord Record, ParseState *PState)
{ 
  if (    PState->Status!=ParseState::INELEMENT
       || PState->CurrentElement->Type!=TEXT
     )
   return new std::string("unexpected record TEXTTYPE");
  PState->CurrentElement->TextType = Record.iVal[0];
  return 0;
}

std::string *handlePATHTYPE(GDSIIRecord Record, ParseState *PState)
{ 
  if (    PState->Status!=ParseState::INELEMENT )
   return new std::string("unexpected record PATHTYPE");
  PState->CurrentElement->PathType = Record.iVal[0];
  return 0;
}

std::string *handleSTRANS(GDSIIRecord Record, ParseState *PState)
{ 
  if ( PState->Status!=ParseState::INELEMENT )
   return new std::string("unexpected record STRANS");
  PState->CurrentElement->Refl     = Record.Bits[0];
  PState->CurrentElement->AbsMag   = Record.Bits[13];
  PState->CurrentElement->AbsAngle = Record.Bits[14];
  return 0;
}

std::string *handleMAG(GDSIIRecord Record, ParseState *PState)
{ 
  if ( PState->Status!=ParseState::INELEMENT )
   return new std::string("unexpected record MAG");
  PState->CurrentElement->Mag = Record.dVal[0];
  return 0;
}

std::string *handleANGLE(GDSIIRecord Record, ParseState *PState)
{ 
  if ( PState->Status!=ParseState::INELEMENT )
   return new std::string("unexpected record ANGLE");
  PState->CurrentElement->Angle = Record.dVal[0];
  return 0;
}

std::string *handlePROPATTR(GDSIIRecord Record, ParseState *PState)
{ 
  if ( PState->Status!=ParseState::INELEMENT )
   return new std::string("unexpected record PROPATTR");
  GDSIIElement *e=PState->CurrentElement;
  e->PropAttrs.push_back(Record.iVal[0]);
  e->PropValues.push_back("");
  return 0;
}

std::string *handlePROPVALUE(GDSIIRecord Record, ParseState *PState)
{
  if ( PState->Status!=ParseState::INELEMENT )
   return new std::string("unexpected record PROPVALUE");
  GDSIIElement *e=PState->CurrentElement;
  int n=e->PropAttrs.size();
  if (n==0)
   return new std::string("PROPVALUE without PROPATTR");
  e->PropValues[n-1]=std::string( *(Record.sVal) );

  if( strcasestr( Record.sVal->c_str(), "CONTEXT_INFO") )
   PState->CurrentStruct->IsPCell=true;

  return 0;
}

std::string *handleXY(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new std::string("unexpected record XY");
  PState->CurrentElement->XY.reserve(Record.NumVals);
  for(size_t n=0; n<Record.NumVals; n++)
   PState->CurrentElement->XY.push_back(Record.iVal[n]);
  return 0;
}

std::string *handleSNAME(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new std::string("unexpected record SNAME");
  PState->CurrentElement->SName = new std::string( *(Record.sVal) );
  return 0;
}

std::string *handleSTRING(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new std::string("unexpected record STRING");
  PState->CurrentElement->Text = new std::string( *(Record.sVal) );
  return 0;
}

std::string *handleCOLROW(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new std::string("unexpected record COLROW");
  PState->CurrentElement->Columns = Record.iVal[0];
  PState->CurrentElement->Rows    = Record.iVal[1];
  return 0;
}

std::string *handleWIDTH(GDSIIRecord Record, ParseState *PState)
{
  if (PState->Status!=ParseState::INELEMENT)
   return new std::string("unexpected record Width");
  PState->CurrentElement->Width   = Record.iVal[0];
  return 0;
}

std::string *handleENDEL(GDSIIRecord Record, ParseState *PState)
{
  (void) Record;
  if (PState->Status!=ParseState::INELEMENT)
   return new std::string("unexpected record ENDEL");
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
// Non-allowed characters at the end of the std::string are removed.
// Non-allowed characters not at the end of the std::string are converted to underscores.
bool IsAllowedChar(char c)
{ return isprint(c) && c!='"' && c!=','; }

std::string *MakeGDSIIString(char *Original, int Size)
{ 
  if (Size==0) return new std::string("");

  if (Size>32) Size=32;
  char RawString[33];
  strncpy(RawString, Original, Size);
  RawString[Size]=0;
  int L = strlen(RawString);
  while ( L>0 && !IsAllowedChar(RawString[L-1]) )
   RawString[--L] = 0;
  for(int n=0; n<L; n++) 
   if (!IsAllowedChar(RawString[n])) RawString[n]='_';
  return new std::string(RawString);
}

/***************************************************************/
/* read a single GDSII data record from the current file position */
/***************************************************************/
GDSIIRecord ReadGDSIIRecord(FILE *f, std::string **ErrMsg)
{
  /*--------------------------------------------------------------*/
  /* read the 4-byte file header and check that the data type     */
  /* agrees with what it should be based on the record type       */
  /*--------------------------------------------------------------*/
  BYTE Header[4];
  if ( 4 != fread(Header, 1, 4, f) )
   { *ErrMsg = new std::string("unexpected end of file");
     return GDSIIRecord(); // end of file
   }

  size_t RecordSize = Header[0]*256 + Header[1];
  BYTE RType        = Header[2];
  BYTE DType        = Header[3];
  
  if (RType > MAX_RTYPE)
   { *ErrMsg = new std::string("unknown record type");
     return GDSIIRecord();
   }
    
  if ( DType != RecordTypes[RType].DType )
   { std::ostringstream ss;
     ss << RecordTypes[RType].Name
        << ": data type disagrees with record type ("
        << DType
        << " != "
        << RecordTypes[RType].DType
        << ")";
     *ErrMsg = new std::string(ss.str());
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
      { *ErrMsg = new std::string("out of memory");
        return GDSIIRecord();
      }
     if ( PayloadSize != fread((void *)Payload, 1, PayloadSize, f) )
      { delete[] Payload;
        *ErrMsg = new std::string("unexpected end of file");
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
       *ErrMsg = new std::string("unknown data type " + std::to_string(DType));
       return GDSIIRecord();
   };

  // success 
  *ErrMsg=0;
  delete[] Payload;
  return Record;

}

/***************************************************************/
/* GDSIIData constructor: create a new GDSIIData instance from */
/* a binary GDSII file.                                        */
/***************************************************************/
GDSIIData::GDSIIData(const std::string FileName)
{ 
  // initialize class data
  LibName       = 0;
  FileUnits[0]  = 1.0e-3; // these seem to be the default for GDSII files
  FileUnits[1]  = 1.0e-9;
  UnitInMeters  = 1.0e-6;
  GDSIIFileName = new std::string(FileName);
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
int GDSIIData::GetStructByName(std::string Name)
{ for(size_t ns=0; ns<Structs.size(); ns++)
   if ( Name == *(Structs[ns]->Name) )
    return ns;
  return -1;
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
void GDSIIData::ReadGDSIIFile(const std::string FileName, double CoordinateLengthUnit)
 {
   ErrMsg=0;

   /*--------------------------------------------------------------*/
   /*- try to open the file ---------------------------------------*/
   /*--------------------------------------------------------------*/
   FILE *f=fopen(FileName.c_str(),"r");
   if (!f)
    { ErrMsg = new std::string("could not open " + FileName);
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
   for(std::set<int>::iterator it=LayerSet.begin(); it!=LayerSet.end(); it++)
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
   //Flatten(CoordinateLengthUnit);
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

#endif
