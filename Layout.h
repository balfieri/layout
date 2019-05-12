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
#include <map>
#include <iostream>
#include <fstream>

#include <errno.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <zlib.h>

#include "libGDSII.h"

class Layout
{
public:
    typedef uint32_t uint;                  // by default, we use 32-bit integers
    typedef int32_t  _int;                  // by default, we use 32-bit integers
    typedef double   real;

    Layout( std::string top_file );
    ~Layout(); 

    bool write( std::string file_path );

    static void dissect_path( std::string path, std::string& dir_name, std::string& base_name, std::string& ext_name ); // utility

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
        uint        char_cnt;               // in strings array

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
            uint        s_i;                // STR or ID - index into strings array of string value
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
    std::map<std::string, uint>         str_to_str_i;          // maps string to unique location in strings[]

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

    bool read_layout( std::string file_path );
    bool write_layout( std::string file_path );         // .layout

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
    bool parse_real( real& r, char *& xxx, char *& xxx_end, bool skip_whitespace_first=true );
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
    switch( kind )
    {
        case Layout::NODE_KIND::STR:            os << "STR";       break;
        case Layout::NODE_KIND::BOOL:           os << "BOOL";      break;
        case Layout::NODE_KIND::INT:            os << "INT";       break;
        case Layout::NODE_KIND::UINT:           os << "UINT";      break;
        case Layout::NODE_KIND::REAL:           os << "REAL";      break;
        case Layout::NODE_KIND::ID:             os << "ID";        break;
        case Layout::NODE_KIND::CALL:           os << "CALL";      break;
        case Layout::NODE_KIND::SLICE:          os << "SLICE";     break;
        case Layout::NODE_KIND::ASSIGN:         os << "ASSIGN";    break;
        case Layout::NODE_KIND::HIER:           os << "HIER";      break;
        default:                                os << "<unknown>"; break;
    }
    return os;
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

        if ( ext_name == std::string( ".aedt" ) ) {
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
            rtn_assert( 0, "unable to parse an expression: string, number, or id " + surrounding_lines( nnn_save, nnn_end ) );
        }
        dprint( "ID START " + std::string(&strings[id_i]) );
        if ( id_i == aedt_begin_str_i ) {
            uint name_i;
            if ( !parse_aedt_expr( name_i ) ) return false;             // STR node
            rtn_assert( nodes[name_i].kind == NODE_KIND::STR, "$begin not followed by string" );
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
    if ( *xxx != '.' ) {
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
        rtn_assert( xxx != xxx_end, "no terminating \" for string" );
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
           parse_real( r3.c[0], xxx, xxx_end, has_brackets ) && 
           (!has_brackets || expect_char( ',', xxx, xxx_end, true )) &&
           parse_real( r3.c[1], xxx, xxx_end, has_brackets ) && 
           (!has_brackets || expect_char( ',', xxx, xxx_end, true )) &&
           parse_real( r3.c[2], xxx, xxx_end, has_brackets ) &&
           (!has_brackets || expect_char( ']', xxx, xxx_end, true ));
}

inline bool Layout::parse_real( Layout::real& r, char *& xxx, char *& xxx_end, bool skip_whitespace_first )
{
    if ( skip_whitespace_first ) skip_whitespace( xxx, xxx_end );   // can span lines unlike below
    bool vld = false;
    bool is_neg = false;
    bool in_frac = false;
    bool has_exp = false;
    uint u = 0;     // integer part
    uint f = 0;     // frac part before divide
    _int e10 = 0;
    double f_factor = 1.0;
    dprint( "parse_real *xxx=" + std::string( 1, *xxx ) );
    while( xxx != xxx_end && (*xxx == ' ' || *xxx == '\t') ) xxx++;  // skip leading spaces

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
            is_neg = true;
            xxx++;
            continue;
        }

        if ( ch == '.' && !in_frac ) {
            in_frac = true;
            xxx++;
            continue;
        }

        if ( ch == 'e' || ch == 'E' ) {
            rtn_assert( !has_exp, "real has more than one 'e' exponent" );
            has_exp = true;
            xxx++;
            if ( !parse_int( e10, xxx, xxx_end ) ) return false;
            continue;
        }

        if ( ch < '0' || ch > '9' ) break;

        uint digit = ch - '0';
        if ( in_frac ) {
            f = 10*f + digit;
            f_factor *= 10.0;
        } else {
            u = 10*u + digit;
        }

        vld = true;
        xxx++;
    }

    r = double(u) + f/f_factor;
    if ( is_neg ) r = -r;
    if ( e10 != 0 ) r *= pow( 10.0, e10 );
    dprint( "real=" + std::to_string( r ) );
    rtn_assert( vld, "unable to parse real in " + ext_name + " file " + surrounding_lines( xxx, xxx_end ) );
    return vld;
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

#endif
