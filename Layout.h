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
//                                              // default is uncompressed, which is larger but faster; 
//                                              // pass true as second argument to get compressed output
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
    typedef uint64_t uint64;                // by default, we use 64-bit indexes for texels
    typedef double   real;

    Layout( std::string top_file );
    Layout( std::string layout_file, bool is_compressed );
    ~Layout(); 

    bool write( std::string file_path, bool is_compressed=false ); 

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
        uint64      byte_cnt;               // total in-memory bytes including this header
        uint64      char_cnt;               // in strings array

        uint64      node_cnt;               // in nodes array  
    };

    enum class NODE_KIND
    {
        STR,
        BOOL,
        INT,
        UINT,
        REAL,
        HIER,                          
        CALL,
    };
            
    class Node
    {
    public:
        NODE_KIND   kind;                   // kind of node
        uint        name_i;                 // index of node name in strings array of node name, if any, else uint(-1)
        union {
            uint        s_i;                // KIND_STR - index into strings array of string value
            bool        b;                  // KIND_BOOL
            _int        i;                  // KIND_INT
            uint        u;                  // KIND_UINT
            real        r;                  // KIND_REAL
            uint        child_first_i;      // KIND_HIER - index into nodes[] array of first child
            uint        arg_first_i;        // KIND_CALL - index into nodes[] array of first arg to call
        } u;
        uint        sibling_i;              // index in nodes array of sibling on list, else uint(-1)
    };

    enum class INSTANCE_KIND
    {
        LAYOUT,                             // instance is another Layout (by name)
        LAYOUT_PTR,                         // instance is another Layout (read in and with a resolved pointer)
        NODE,                               // instance is a node within this Layout
    };

    class Instance      
    {
    public:
        INSTANCE_KIND   kind;               // see above
        uint            name_i;             // index in strings array of instance name
        uint            layout_name_i;      // index in strings array of name of instanced Layout 
        uint            layout_file_name_i; // index in strings array of file name of instanced Layout 
        real3           translation;        // 3D translation of instance
        union {
            Layout *    layout_ptr;         // resolved Layout pointer for LAYOUT_PTR
            uint        node_i;             // index of Node in nodes[] array
        } u;

        bool bounding_box( const Layout * layout, AABB& box, real padding=0 ) const;
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
    char * node_start;
    char * node_end;
    char * node_c;
    uint   line_num;

    uint aedt_begin_str_i;              // these are to make it easier to compare
    uint aedt_end_str_i;
    uint true_str_i;
    uint false_str_i;

    bool load_aedt( std::string node_file, std::string dir_name );        // parse .aedt
    bool parse_aedt_node( uint& node_i, uint id_i );

    bool open_and_read( std::string file_name, char *& start, char *& end );

    bool skip_whitespace_to_eol( char *& xxx, char *& xxx_end );  // on this line only
    bool skip_whitespace( char *& xxx, char *& xxx_end );
    bool skip_to_eol( char *& xxx, char *& xxx_end );
    bool eol( char *& xxx, char *& xxx_end );
    bool expect_char( char ch, char *& xxx, char* xxx_end, bool skip_whitespace_first=false );
    uint get_str_i( std::string s );
    bool parse_expr( uint node_i, char *& xxx, char *& xxx_end );
    bool parse_number( uint node_i, char *& xxx, char *& xxx_end );
    bool parse_string( std::string& s, char *& xxx, char *& xxx_end );
    bool parse_string_i( uint& s, char *& xxx, char *& xxx_end );
    bool parse_name( char *& name, char *& xxx, char *& xxx_end );
    bool parse_id( uint& id_i, char *& xxx, char *& xxx_end );
    bool parse_real3( real3& r3, char *& xxx, char *& xxx_end, bool has_brackets=false );
    bool parse_real( real& r, char *& xxx, char *& xxx_end, bool skip_whitespace_first=true );
    bool parse_int( _int& i, char *& xxx, char *& xxx_end );
    bool parse_uint( uint& u, char *& xxx, char *& xxx_end );
    bool parse_bool( bool& b, char *& xxx, char *& xxx_end );
    std::string surrounding_lines( char *& xxx, char *& xxx_end );

    // allocates an array of T on a page boundary
    template<typename T>
    T * aligned_alloc( uint64 cnt );

    // reallocate array if we are about to exceed its current size
    template<typename T>
    inline void perhaps_realloc( T *& array, const uint64& hdr_cnt, uint64& max_cnt, uint64 add_cnt );

    bool write_uncompressed( std::string file_path );
    bool read_uncompressed( std::string file_path );
};

#define dprint( msg )
//#define dprint( msg ) std::cout << (msg) << "\n"

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

Layout::Layout( std::string top_file )
{
    is_good = false;
    error_msg = "<unknown error>";
    mapped_region = nullptr;
    node_c = nullptr;
    line_num = 1;

    hdr = aligned_alloc<Header>( 1 );
    memset( hdr, 0, sizeof( Header ) );

    //------------------------------------------------------------
    // Initial lengths of arrays are large in virtual memory
    //------------------------------------------------------------
    max = aligned_alloc<Header>( 1 );
    max->node_cnt     =  1024;
    max->char_cnt    = max->node_cnt * 128;

    //------------------------------------------------------------
    // Allocate arrays
    //------------------------------------------------------------
    strings         = aligned_alloc<char>(     max->char_cnt );
    nodes         = aligned_alloc<Node>(   max->node_cnt );

    //------------------------------------------------------------
    // Load node depending on file ext_name
    //------------------------------------------------------------
    std::string dir_name;
    std::string base_name;
    dissect_path( top_file, dir_name, base_name, ext_name );
    if ( ext_name == std::string( ".aedt" ) ) {
        if ( !load_aedt( top_file, dir_name ) ) return;
    } else {
        error_msg = "unknown top file ext_name: " + ext_name;
        return;
    }

    //------------------------------------------------------------
    // Add up byte count.
    //------------------------------------------------------------
    hdr->byte_cnt = uint64( 1                 ) * sizeof( hdr ) +
                    uint64( hdr->node_cnt      ) * sizeof( nodes[0] ) +
                    uint64( hdr->char_cnt     ) * sizeof( strings[0] );

    is_good = true;
}

Layout::Layout( std::string layout_path, bool is_compressed )
{
    is_good = false;
    mapped_region = nullptr;
    if ( !is_compressed ) {
        read_uncompressed( layout_path );
        return;
    }

    gzFile fd = gzopen( layout_path.c_str(), "r" );
    if ( fd == Z_NULL ) {
        "Could not gzopen() file " + layout_path + " for reading - gzopen() error: " + strerror( errno );
        return;
    }

    //------------------------------------------------------------
    // Reader in header then individual arrays.
    //------------------------------------------------------------
    #define _read( array, type, cnt ) \
        if ( cnt == 0 ) { \
            array = nullptr; \
        } else { \
            array = aligned_alloc<type>( cnt ); \
            if ( array == nullptr ) { \
                gzclose( fd ); \
                error_msg = "could not allocate " #array " array"; \
                assert( 0 ); \
                return; \
            } \
            char * _addr = reinterpret_cast<char *>( array ); \
            for( uint64 _byte_cnt = (cnt)*sizeof(type); _byte_cnt != 0;  ) \
            { \
                uint64 _this_byte_cnt = 1024*1024*1024; \
                if ( _byte_cnt < _this_byte_cnt ) _this_byte_cnt = _byte_cnt; \
                if ( gzread( fd, _addr, _this_byte_cnt ) <= 0 ) { \
                    gzclose( fd ); \
                    error_msg = "could not gzread() file " + layout_path + " - gzread() error: " + strerror( errno ); \
                    assert( 0 ); \
                    return; \
                } \
                _byte_cnt -= _this_byte_cnt; \
                _addr     += _this_byte_cnt; \
            } \
        } \

    _read( hdr,         Header,   1 );
    if ( hdr->version != VERSION ) {
        gzclose( fd );
        error_msg = "hdr->version does not match VERSION=" + std::to_string(VERSION) + ", got " + std::to_string(hdr->version);
        assert( 0 ); \
        return;
    }
    max = aligned_alloc<Header>( 1 );
    memcpy( max, hdr, sizeof( Header ) );
    _read( strings,     char,          hdr->char_cnt );
    _read( nodes,     Node,        hdr->node_cnt );

    gzclose( fd );

    is_good = true;
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

bool Layout::write( std::string layout_path, bool is_compressed ) 
{
    if ( !is_compressed ) return write_uncompressed( layout_path );

    gzFile fd = gzopen( layout_path.c_str(), "w" );
    rtn_assert( fd != Z_NULL, "could not gzopen() file " + layout_path + " for writing - gzopen() error: " + strerror( errno ) );

    //------------------------------------------------------------
    // Write out header than individual arrays.
    //------------------------------------------------------------
    #define _write( addr, byte_cnt ) \
    { \
        char * _addr = reinterpret_cast<char *>( addr ); \
        for( uint64 _byte_cnt = byte_cnt; _byte_cnt != 0;  ) \
        { \
            uint64 _this_byte_cnt = 1024*1024*1024; \
            if ( _byte_cnt < _this_byte_cnt ) _this_byte_cnt = _byte_cnt; \
            if ( gzwrite( fd, _addr, _this_byte_cnt ) <= 0 ) { \
                gzclose( fd ); \
                rtn_assert( 0, "could not gzwrite() file " + layout_path + " - gzwrite() error: " + strerror( errno ) ); \
            } \
            _byte_cnt -= _this_byte_cnt; \
            _addr     += _this_byte_cnt; \
        } \
    } \

    _write( hdr,         1                  * sizeof(hdr[0]) );
    _write( strings,     hdr->char_cnt      * sizeof(strings[0]) );
    _write( nodes,     hdr->node_cnt       * sizeof(nodes[0]) );

    gzclose( fd );
    return true;
}

// returns array of T on a page boundary
template<typename T>
T * Layout::aligned_alloc( Layout::uint64 cnt )
{
    void * mem = nullptr;
    posix_memalign( &mem, getpagesize(), cnt*sizeof(T) );
    return reinterpret_cast<T *>( mem );
}

// reallocate array if we are about to exceed its current size
template<typename T>
inline void Layout::perhaps_realloc( T *& array, const Layout::uint64& hdr_cnt, Layout::uint64& max_cnt, Layout::uint64 add_cnt )
{
    while( (hdr_cnt + add_cnt) > max_cnt ) {
        void * mem = nullptr;
        uint64 old_max_cnt = max_cnt;
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

bool Layout::write_uncompressed( std::string layout_path ) 
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
    _uwrite( nodes,     hdr->node_cnt       * sizeof(nodes[0]) );

    fsync( fd ); // flush
    close( fd );

    return true;
}

bool Layout::read_uncompressed( std::string layout_path )
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
    _uread( strings,     char,          hdr->char_cnt );
    _uread( nodes,     Node,        hdr->node_cnt );

    is_good = true;

    return true;
}

bool Layout::load_aedt( std::string node_file, std::string dir_name )
{
    (void)dir_name;

    //------------------------------------------------------------
    // Map in file
    //------------------------------------------------------------
    line_num = 1;
    if ( !open_and_read( node_file, node_start, node_end ) ) return false;
    node_c = node_start;

    //------------------------------------------------------------
    // Parse .aedt file contents
    // Parse the first $begin .. $end.
    //------------------------------------------------------------
    aedt_begin_str_i = get_str_i( "$begin" );
    aedt_end_str_i   = get_str_i( "$end" );
    true_str_i       = get_str_i( "true" );
    false_str_i      = get_str_i( "false" );
    uint id_i;
    uint root_i;
    if ( !parse_id( id_i, node_c, node_end ) ) return false;
    if ( !parse_aedt_node( root_i, id_i ) ) return false;
    assert( root_i == 0 );
    assert( nodes[root_i].kind == NODE_KIND::HIER );
    return true;
}

bool Layout::parse_aedt_node( uint& node_i, uint id_i )
{
    printf( "node id_i=%d id=%s\n", id_i, &strings[id_i] );
    perhaps_realloc( nodes, hdr->node_cnt, max->node_cnt, 1 );
    uint ni = hdr->node_cnt++;
    if ( id_i == aedt_begin_str_i ) {
        // HIER => parse child nodes
        uint str_i;
        if ( !parse_string_i( str_i, node_c, node_end ) ) return false;
        std::cout << "$begin " << std::string(&strings[str_i]) << "\n";
        nodes[ni].kind = NODE_KIND::HIER;
        nodes[ni].name_i = str_i;
        nodes[ni].u.child_first_i = uint(-1);
        nodes[ni].sibling_i = uint(-1);
        uint prev_i = uint(-1);
        for( ;; )
        {
            uint child_id_i;
            if ( !parse_id( child_id_i, node_c, node_end ) ) return false;
            if ( child_id_i == aedt_end_str_i ) {
                uint end_str_i;
                if ( !parse_string_i( end_str_i, node_c, node_end ) ) return false;
                rtn_assert( end_str_i == str_i, "$end id does not match $begin id " + surrounding_lines( node_c, node_end ) );
                break;
            }

            uint child_i;
            if ( !parse_aedt_node( child_i, child_id_i ) ) return false;
            if ( nodes[ni].u.child_first_i == uint(-1) ) {
                nodes[ni].u.child_first_i = child_i;
            } else {
                nodes[prev_i].sibling_i = child_i;
            }
            prev_i = child_i;
        }
    } else {
        nodes[ni].name_i = id_i;
        char ch = *node_c;
        if ( ch == '=' ) {
            // assignment
            if ( !expect_char( ch, node_c, node_end ) ) return false;
            if ( !parse_expr( ni, node_c, node_end ) ) return false;
        } else if ( ch == '(' ) {
            // CALL => parse arg list
            nodes[ni].kind = NODE_KIND::CALL;
            nodes[ni].u.child_first_i = uint(-1);
            nodes[ni].sibling_i = uint(-1);
            uint prev_i = uint(-1);
            if ( !expect_char( ch, node_c, node_end ) ) return false;
            for( bool have_one=false; ; have_one=true )
            {
                if ( !skip_whitespace( node_c, node_end ) ) return false;
                ch = *node_c;
                if ( ch == ')' ) {
                    if ( !expect_char( ch, node_c, node_end ) ) return false;
                    break;
                }

                if ( have_one ) {
                    if ( !expect_char( ',', node_c, node_end ) ) return false;
                    if ( !skip_whitespace( node_c, node_end ) ) return false;
                }
                perhaps_realloc( nodes, hdr->node_cnt, max->node_cnt, 1 );
                uint arg_i = hdr->node_cnt++;
                nodes[arg_i].name_i = uint(-1);
                nodes[arg_i].u.child_first_i = uint(-1);
                nodes[arg_i].sibling_i = uint(-1);
                if ( nodes[ni].u.child_first_i == uint(-1) ) {
                    nodes[ni].u.child_first_i = arg_i;
                } else {
                    nodes[prev_i].sibling_i = arg_i;
                }
                prev_i = arg_i;
                if ( !parse_expr( arg_i, node_c, node_end ) ) return false;
            }
            std::cout << "END CALL\n";
        } else if ( ch == '\'' ) {
            // STR
            uint str_i;
            if ( !parse_string_i( str_i, node_c, node_end ) ) return false;
        } else {
            rtn_assert( false, "unknown .aedt node id=" + std::string(&strings[id_i]) + surrounding_lines( node_c, node_end ) );
        }
    }
    return true;
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

inline bool Layout::skip_whitespace( char *& xxx, char *& xxx_end )
{
    bool in_comment = false;
    for( ;; )
    {
        if ( xxx == xxx_end ) break;

        char ch = *xxx;
        if ( ch == '#' ) in_comment = true;
        if ( !in_comment && ch != ' ' && ch != '\n' && ch != '\r' && ch != '\t' ) break;

        if ( ch == '\n' || ch == '\r' ) {
            if ( ch == '\n' && (xxx == node_c) ) line_num++;
            in_comment = false;
        }
        xxx++;
    }
    return true;
}

inline bool Layout::skip_whitespace_to_eol( char *& xxx, char *& xxx_end )
{
    bool in_comment = false;
    for( ;; )
    {
        if ( xxx == xxx_end ) break;

        char ch = *xxx;
        if ( ch == '#' ) in_comment = true;
        if ( !in_comment && ch != ' ' && ch != '\n' && ch != '\r' && ch != '\t' ) break;

        if ( ch == '\n' || ch == '\r' ) {
            if ( ch == '\n' && (xxx == node_c) ) line_num++;
            break;
        }
        xxx++;
    }
    return true;
}

inline bool Layout::skip_to_eol( char *& xxx, char *& xxx_end )
{
    if ( !eol( xxx, xxx_end ) ) {
        while( xxx != xxx_end )
        {
            char ch = *xxx;
            if ( ch == '\n' || ch == '\r' ) break;
            xxx++;
        }
    }
    return true;
}

inline bool Layout::eol( char *& xxx, char *& xxx_end )
{
    skip_whitespace_to_eol( xxx, xxx_end );

    if ( xxx == xxx_end || *xxx == '\n' || *xxx == '\r' ) {
        if ( xxx != xxx_end ) {
            if ( *xxx == '\n' && (xxx == node_c) ) line_num++;
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
    std::cout << "str_i[" << s << "]=" << s_i << " strings[]=" << std::string(&strings[s_i]) << "\n";
    return s_i;
}

inline bool Layout::parse_expr( uint node_i, char *& xxx, char *& xxx_end )
{
    char ch = *xxx;
    if ( ch == '\'' ) {
        nodes[node_i].kind = NODE_KIND::STR;
        return parse_string_i( nodes[node_i].u.s_i, xxx, xxx_end );
    } else if ( ch == '-' || (ch >= '0' && ch <= '9') ) {
        return parse_number( node_i, xxx, xxx_end );
    } else {
        nodes[node_i].kind = NODE_KIND::BOOL;
        return parse_bool( nodes[node_i].u.b, xxx, xxx_end );
    }
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
    rtn_assert( id != "", "no id found at " + surrounding_lines( xxx, xxx_end ) );
    id_i = get_str_i( id );
    return true;
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
    std::string s = "";
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
