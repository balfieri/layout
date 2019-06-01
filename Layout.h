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

    const uint NULL_I = -1;                 // null index into an array

    // FILE PATHS
    // 
    static void        dissect_path( std::string path, std::string& dir_name, std::string& base_name, std::string& ext_name ); // utility
    static std::string path_without_ext( std::string path, std::string * file_ext=nullptr );  // returns path without file extension part

    // LAYOUTS
    //
    Layout( void );
    Layout( std::string top_file, bool count_only=false );
    ~Layout(); 

    // output file format depends on file extension:
    bool write( std::string file_path );

    // write out layer info used by viewer programs (file extension determines format)
    bool write_layer_info( std::string file_path );

    // current layout dimensions
    real width( void ) const;
    real length( void ) const;
    real height( void ) const;

    // start of a new library; Layout must be empty
    uint start_library( std::string libname, real units_user=0.001, real units_meters=1e-9 );

    // instancing of other layouts
    uint inst_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, real x, real y, std::string name );
    uint inst_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, real x, real y, 
                      uint dst_layer_first, uint dst_layer_last, std::string name );
    void finalize_top_struct( uint parent_i, uint last_i, std::string top_name );              // use to create top-level struct of all insts
    uint flatten_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, std::string dst_struct_name="" );  // straight copy with flattening
    
    // fill of dielectrics or arbitrary material
    void fill_dielectrics( void );
    void fill_material( uint material_i, real x, real y, real z, real w, real l, real h );

    // PUBLIC DATA STRUCTURES
    //
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // NOTE: THERE ARE NO POINTERS IN THESE DATA STRUCTURES.
    //       WE STORE INDEXES INTO ARRAYS INSTEAD.
    //       ARRAYS ARE ALL PAGE-ALIGNED SO THAT FILE WRITE/READ
    //       TURNS INTO MMAP() EFFECTIVELY.  
    //       INDEXES ALLOW ARRAYS TO BE EASILY RELOCATED DURING
    //       RESIZING OR READING IN FROM A FILE.
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    class Header                            // header (of future binary file)
    {
    public:
        uint        version;                // version

        uint        char_cnt;               // in strings array
        uint        material_cnt;           // in materials array
        uint        layer_cnt;              // in layers array
        uint        node_cnt;               // in nodes array  
        uint        top_inst_cnt;           // in top_insts array
        uint        root_i;                 // index of root node in nodes array
    };

    // returns index into strings[] for s;
    // creates new index if s is not already in strings[]
    //
    uint str_get( std::string s );
    uint str_find( std::string s ) const;

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

    // this will clear materials and set up default set of common materials
    //
    void        materials_init( void );

    // these return index in materials[] array or NULL_I when failure or not found
    // set() will override material properties if name already exists
    //
    uint        material_set( std::string name, const Material& material );
    uint        material_get( std::string name );

    class Layer                             // layer stackup (one layer)
    {
    public:
        uint        name_i;                 	// index of layer name in strings[]
        uint        gdsii_num;              	// layer number of main material in GDSII file
        uint        dielectric_gdsii_num;   	// layer number of dielectric material in GDSII file
        uint        gdsii_datatype;         	// which datatype to use - NULL_I means all
        uint        dielectric_gdsii_datatype; 	// which datatype to use - NULL_I means all
        bool        same_zoffset_as_prev;   	// starts at same zoffset as previous layer in stackkup?
        real        thickness;              	// thickness in um
        uint        material_i;             	// index of main material in materials[]
        uint        dielectric_material_i;  	// index of dielectric material in materials[]
        uint        material_rgba;              // material color in RGBA8 format
    };

    static uint color( real r, real g, real b, real a=1.0 );
    static uint color( std::string name, real a=1.0 );

    // these return index in layers[] array or NULL_I when failure or not found
    // set() will override layer properties if name already exists
    //
    void        layer_set( uint layer_i, const Layer& layer );
    uint        layer_get( std::string name );

    enum class NODE_KIND
    {
        // these match the exact order given in the GDSII spec, so that they have int values that match
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
        // END OF GDSII KINDS

        // scalars
        STR,                                
        BOOL,
        INT,
        UINT,
        REAL,
        ID,

        // non-scalars
        ASSIGN,                                 // child 0 is lhs (usually an id), child 1 is rhs
        HIER,                                   // child 0 is id, other children are normal children
        CALL,                                   // child 0 is id, other children are args
        SLICE,                                  // child 0 is id, child 1 is index before the ':', other children are other args
    };
            
    static constexpr uint32_t GDSII_KIND_CNT = uint(NODE_KIND::LIBSECUR) + 1;

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
        uint        sibling_i;              // index in nodes array of sibling on list, else NULL_I

    };

    bool        node_is_header_footer( const Node& node ) const;// return true if node is a HEADER, BGNLIB, LIBNAME, UNITS, or ENDLIB
    bool        node_is_gdsii( const Node& node ) const;        // return true if node is a GDSII node
    bool        node_is_scalar( const Node& node ) const;       // return true if node is a scalar 
    bool        node_is_string( const Node& node ) const;       // return true if node is a (scalar) string
    bool        node_is_parent( const Node& node ) const;       // return true if node is not a scalar (i.e., could have children)
    bool        node_is_element( const Node& node ) const;      // return true if node is a GDSII element
    bool        node_is_name( const Node& node ) const;         // return true if node is a name node
    bool        node_is_hier( const Node& node ) const;         // return true if node is a hierarchy
    bool        node_is_ref( const Node& node ) const;          // return true if node is an AREF or SREF
    uint        node_last_scalar_i( const Node& node ) const;   // find node index of last scalar child, else NULL_I if none
    uint        node_name_i( const Node& node ) const;          // find name for node but return strings[] index
    std::string node_name( const Node& node ) const;            // find name for node
    uint        node_layer( const Node& node ) const;           // find LAYER value for node (an element)
    uint        node_xy_i( const Node& node ) const;            // find index of XY node within node

    static std::string  str( NODE_KIND kind );
    NODE_KIND           hier_end_kind( NODE_KIND kind ) const;      // returns corresponding end kind for hier kind
    
    enum class COPY_KIND
    {
        ONE,                                // copy only the one source node, no children
        SCALAR_CHILDREN,                    // copy scalar children only (INT, REAL, etc.)
        DEEP,                               // copy all children and descendents
        FLATTEN,                            // copy all children but flatten all REFs
    };
    uint        node_alloc( NODE_KIND kind );                   // allocate a node of the given kind
    uint        node_copy( uint parent_i, uint last_i, const Layout * src_layout, uint src_i, COPY_KIND kind, real x=0, real y=0, real angle=0, bool reflect=false, bool in_flatten=false );

    static GDSII_DATATYPE kind_to_datatype( NODE_KIND kind );
    static std::string    str( GDSII_DATATYPE datatype );

    // global scalars
    static const uint   VERSION = 0xB0BA1f01; // current version 

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




//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//
//  ___                 _                           _        _   _
// |_ _|_ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
//  | || '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
//  | || | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
// |___|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
//               |_|
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
    uint        gdsii_rec_cnt;
    NODE_KIND   gdsii_last_kind;
    int         gdsii_fd;
    uint8_t *   gdsii_buff;
    uint        gdsii_buff_byte_cnt;
    real        gdsii_units_user;
    real        gdsii_units_meters;

    // state used during reading and and writing of AEDT files
    uint aedt_begin_str_i;              // these are to make it easier to compare
    uint aedt_end_str_i;
    uint true_str_i;
    uint false_str_i;

    // struct info
    std::map< uint, uint >                      name_i_to_struct_i;
    std::map< uint, std::map<uint, bool> * > *  struct_i_to_has_layer;

    void init( bool alloc_arrays );

    bool layout_read( std::string file_path );          // .layout
    bool layout_write( std::string file_path );         

    struct TopInstInfo
    {
        uint            struct_i;
        real            x;
        real            y;
    };
    TopInstInfo *       top_insts;

    using has_layer_cache_t = std::map< uint, std::map<uint, bool>* >;
    bool node_has_layer( uint ni, uint layer_num, has_layer_cache_t * cache, std::string indent_str ) const;
    uint inst_layout_node( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, uint src_i, uint src_layer_num, 
                           uint dst_layer_num, has_layer_cache_t * cache, std::string name, std::string indent_str="" );

    void node_timestamp( Node& node );                  // adds timestamp fields

    bool gdsii_read( std::string file_path, bool count_only ); // .gds
    bool gdsii_read_record( uint& node_i, uint curr_struct_i, bool count_only );
    bool gdsii_write( std::string file );
    void gdsii_write_record( uint node_i, std::string indent_str="" );
    void gdsii_write_number( uint8_t * bytes, uint& byte_cnt, uint ni, GDSII_DATATYPE datatype, std::string indent_str );
    void gdsii_write_bytes( const uint8_t * bytes, uint byte_cnt );
    void gdsii_flush( bool for_end_of_file=false );
    
    bool aedt_read( std::string file );                 // .aedt
    bool aedt_read_expr( uint& node_i );
    bool aedt_write( std::string file );
    void aedt_write_expr( std::ofstream& out, uint node_i, std::string indent_str );

    bool gds3d_write_layer_info( std::string file );
    bool vbs_write_layer_info( std::string file );

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
#define lassert( bool, msg ) if ( !(bool) ) { std::cout << "ERROR: " << std::string(msg) << "\n"; exit( 1 ); }

std::string Layout::str( Layout::NODE_KIND kind )
{
    #define ncase( kind ) case Layout::NODE_KIND::kind: return #kind;
    switch( kind )
    {
        ncase( HEADER )
        ncase( BGNLIB )
        ncase( LIBNAME )
        ncase( UNITS )
        ncase( ENDLIB )
        ncase( BGNSTR )
        ncase( STRNAME )
        ncase( ENDSTR )
        ncase( BOUNDARY )
        ncase( PATH )
        ncase( SREF )
        ncase( AREF )
        ncase( TEXT )
        ncase( LAYER )
        ncase( DATATYPE )
        ncase( WIDTH )
        ncase( XY )
        ncase( ENDEL )
        ncase( SNAME )
        ncase( COLROW )
        ncase( TEXTNODE )
        ncase( NODE )
        ncase( TEXTTYPE )
        ncase( PRESENTATION )
        ncase( UNUSED )
        ncase( STRING )
        ncase( STRANS )
        ncase( MAG )
        ncase( ANGLE )
        ncase( UNUSED2 )
        ncase( UNUSED3 )
        ncase( REFLIBS )
        ncase( FONTS )
        ncase( PATHTYPE )
        ncase( GENERATIONS )
        ncase( ATTRTABLE )
        ncase( STYPTABLE )
        ncase( STRTYPE )
        ncase( ELFLAGS )
        ncase( ELKEY )
        ncase( LINKTYPE )
        ncase( LINKKEYS )
        ncase( NODETYPE )
        ncase( PROPATTR )
        ncase( PROPVALUE )
        ncase( BOX )
        ncase( BOXTYPE )
        ncase( PLEX )
        ncase( BGNEXTN )
        ncase( ENDEXTN )
        ncase( TAPENUM )
        ncase( TAPECODE )
        ncase( STRCLASS )
        ncase( RESERVED )
        ncase( FORMAT )
        ncase( MASK )
        ncase( ENDMASKS )
        ncase( LIBDIRSIZE )
        ncase( SRFNAME )
        ncase( LIBSECUR )

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

inline std::ostream& operator << ( std::ostream& os, const Layout::NODE_KIND& kind ) 
{
    os << Layout::str( kind );
    return os;
}

Layout::GDSII_DATATYPE Layout::kind_to_datatype( Layout::NODE_KIND kind )
{
    #define kdcase( kind, type ) case NODE_KIND::kind: return GDSII_DATATYPE::type;
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

void Layout::init( bool alloc_arrays )
{
    mapped_region = nullptr;
    nnn = nullptr;
    line_num = 1;
    file_path = "";

    if ( alloc_arrays ) {
        //------------------------------------------------------------
        // Initial lengths of arrays are large in virtual memory
        //------------------------------------------------------------
        hdr = aligned_alloc<Header>( 1 );
        memset( hdr, 0, sizeof( Header ) );
        hdr->version = VERSION;
        hdr->root_i = NULL_I;

        max = aligned_alloc<Header>( 1 );
        max->node_cnt =  1024;
        max->char_cnt = max->node_cnt * 128;
        max->material_cnt = 1024;
        max->layer_cnt = 1024;
        max->top_inst_cnt = 1024;

        //------------------------------------------------------------
        // Allocate initial arrays
        //------------------------------------------------------------
        strings    = aligned_alloc<char>( max->char_cnt );
        materials  = aligned_alloc<Material>( max->material_cnt );
        layers     = aligned_alloc<Layer>( max->layer_cnt );
        nodes      = aligned_alloc<Node>( max->node_cnt );
        top_insts  = aligned_alloc<TopInstInfo>( max->top_inst_cnt );

        materials_init();
    }
}

Layout::Layout( void )
{
    init( true );
}

Layout::Layout( std::string top_file, bool count_only ) 
{
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
        init( false );
        if ( !layout_read( top_file ) ) return;
    } else {
        //------------------------------------------------------------
        // Parse .gds or .aedt
        //------------------------------------------------------------
        init( true );
        if ( ext_name == std::string( ".gds" ) ) {
            if ( !gdsii_read( top_file, count_only ) ) return;
        } else if ( ext_name == std::string( ".aedt" ) ) {
            if ( !aedt_read( top_file ) ) return;
        } else {
            lassert( false, "unknown top file ext_name: " + ext_name );
        }

    }
}

Layout::~Layout()
{
    if ( mapped_region != nullptr ) {
        delete mapped_region;
        mapped_region = nullptr;
    } else {
        delete strings;
        delete materials;
        delete layers;
        delete nodes;
        delete top_insts;
        strings = nullptr;
        materials = nullptr;
        layers = nullptr;
        nodes = nullptr;
        top_insts = nullptr;
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
            lassert( false, "unknown file ext_name: " + ext_name );
            return false;
        }
    }
}

bool Layout::write_layer_info( std::string file_path )
{
    //------------------------------------------------------------
    // Write depends on file ext_name
    //------------------------------------------------------------
    std::string dir_name;
    std::string base_name;
    dissect_path( file_path, dir_name, base_name, ext_name );
    if ( ext_name == std::string( ".gds3d" ) ) {
        return gds3d_write_layer_info( file_path );
    } else if ( ext_name == std::string( ".vbs" ) ) {
        return vbs_write_layer_info( file_path );
    } else {
        lassert( false, "write_layer_info: unknown file ext_name: " + ext_name );
        return false;
    }
}

uint Layout::color( real r, real g, real b, real a )
{
    uint ru = uint( r * 255.0 + 0.5 );
    uint gu = uint( g * 255.0 + 0.5 );
    uint bu = uint( b * 255.0 + 0.5 );
    uint au = uint( a * 255.0 + 0.5 );

    if ( ru > 255 ) ru = 255;
    if ( gu > 255 ) gu = 255;
    if ( bu > 255 ) bu = 255;
    if ( au > 255 ) au = 255;

    return (ru << 24) | (gu << 16) | (bu << 8) | (au << 0);
}

class ColorInfo
{
public:
    const char * name;
    uint         rgb;
};

static const ColorInfo color_info[] = 
{
    { "maroon", 		0x800000 },
    { "dark red", 		0x8B0000 },
    { "brown", 			0xA52A2A },
    { "firebrick", 		0xB22222 },
    { "crimson", 		0xDC143C },
    { "red", 			0xFF0000 },
    { "r", 			0xFF0000 },
    { "tomato", 		0xFF6347 },
    { "coral", 			0xFF7F50 },
    { "indian red", 		0xCD5C5C },
    { "light coral", 		0xF08080 },
    { "dark salmon", 		0xE9967A },
    { "salmon", 		0xFA8072 },
    { "light salmon", 		0xFFA07A },
    { "orange red", 		0xFF4500 },
    { "dark orange", 		0xFF8C00 },
    { "orange", 		0xFFA500 },
    { "o", 		        0xFFA500 },
    { "gold", 			0xFFD700 },
    { "dark golden rod", 	0xB8860B },
    { "golden rod", 		0xDAA520 },
    { "pale golden rod", 	0xEEE8AA },
    { "dark khaki", 		0xBDB76B },
    { "khaki", 			0xF0E68C },
    { "olive", 			0x808000 },
    { "yellow", 		0xFFFF00 },
    { "y", 		        0xFFFF00 },
    { "yellow green", 		0x9ACD32 },
    { "dark olive green",	0x556B2F },
    { "olive drab", 		0x6B8E23 },
    { "lawn green", 		0x7CFC00 },
    { "chart reuse", 		0x7FFF00 },
    { "green yellow", 		0xADFF2F },
    { "dark green", 		0x006400 },
    { "green", 			0x00FF00 },
    { "g", 			0x00FF00 },
    { "forest green", 		0x228B22 },
    { "lime", 			0x00FF00 },
    { "lime green", 		0x32CD32 },
    { "light green", 		0x90EE90 },
    { "pale green", 		0x98FB98 },
    { "dark sea green", 	0x8FBC8F },
    { "medium spring green", 	0x00FA9A },
    { "spring green", 		0x00FF7F },
    { "sea green", 		0x2E8B57 },
    { "medium aqua marine", 	0x66CDAA },
    { "medium sea green", 	0x3CB371 },
    { "light sea green", 	0x20B2AA },
    { "dark slate gray", 	0x2F4F4F },
    { "teal", 			0x008080 },
    { "dark cyan", 		0x008B8B },
    { "aqua", 			0x00FFFF },
    { "cyan", 			0x00FFFF },
    { "c", 			0x00FFFF },
    { "light cyan", 		0xE0FFFF },
    { "dark turquoise", 	0x00CED1 },
    { "turquoise", 		0x40E0D0 },
    { "medium turquoise", 	0x48D1CC },
    { "pale turquoise", 	0xAFEEEE },
    { "aqua marine", 		0x7FFFD4 },
    { "powder blue", 		0xB0E0E6 },
    { "cadet blue", 		0x5F9EA0 },
    { "steel blue", 		0x4682B4 },
    { "corn flower blue", 	0x6495ED },
    { "deep sky blue", 		0x00BFFF },
    { "dodger blue", 		0x1E90FF },
    { "light blue", 		0xADD8E6 },
    { "sky blue", 		0x87CEEB },
    { "light sky blue", 	0x87CEFA },
    { "midnight blue", 		0x191970 },
    { "navy", 			0x000080 },
    { "dark blue", 		0x00008B },
    { "medium blue", 		0x0000CD },
    { "blue", 			0x0000FF },
    { "b", 			0x0000FF },
    { "royal blue", 		0x4169E1 },
    { "blue violet", 		0x8A2BE2 },
    { "indigo", 		0x4B0082 },
    { "dark slate blue", 	0x483D8B },
    { "slate blue", 		0x6A5ACD },
    { "medium slate blue", 	0x7B68EE },
    { "medium purple", 		0x9370DB },
    { "dark magenta", 		0x8B008B },
    { "dark violet", 		0x9400D3 },
    { "dark orchid", 		0x9932CC },
    { "medium orchid", 		0xBA55D3 },
    { "purple", 		0x800080 },
    { "p", 		        0x800080 },
    { "thistle", 		0xD8BFD8 },
    { "plum", 			0xDDA0DD },
    { "violet", 		0xEE82EE },
    { "magenta",                0xFF00FF },
    { "m",                      0xFF00FF },
    { "fuchsia", 	        0xFF00FF },
    { "orchid", 		0xDA70D6 },
    { "medium violet red", 	0xC71585 },
    { "pale violet red", 	0xDB7093 },
    { "deep pink", 		0xFF1493 },
    { "hot pink", 		0xFF69B4 },
    { "light pink", 		0xFFB6C1 },
    { "pink", 			0xFFC0CB },
    { "antique white", 		0xFAEBD7 },
    { "beige", 			0xF5F5DC },
    { "bisque", 		0xFFE4C4 },
    { "blanched almond", 	0xFFEBCD },
    { "wheat", 			0xF5DEB3 },
    { "corn silk", 		0xFFF8DC },
    { "lemon chiffon", 		0xFFFACD },
    { "light golden rod yellow",0xFAFAD2 },
    { "light yellow", 		0xFFFFE0 },
    { "saddle brown", 		0x8B4513 },
    { "sienna", 		0xA0522D },
    { "chocolate", 		0xD2691E },
    { "peru", 			0xCD853F },
    { "sandy brown", 		0xF4A460 },
    { "burly wood", 		0xDEB887 },
    { "tan", 			0xD2B48C },
    { "rosy brown", 		0xBC8F8F },
    { "moccasin", 		0xFFE4B5 },
    { "navajo white", 		0xFFDEAD },
    { "peach puff", 		0xFFDAB9 },
    { "misty rose", 		0xFFE4E1 },
    { "lavender blush",		0xFFF0F5 },
    { "linen", 			0xFAF0E6 },
    { "old lace", 		0xFDF5E6 },
    { "papaya whip", 		0xFFEFD5 },
    { "sea shell", 		0xFFF5EE },
    { "mint cream", 		0xF5FFFA },
    { "slate gray", 		0x708090 },
    { "light slate gray", 	0x778899 },
    { "light steel blue", 	0xB0C4DE },
    { "lavender", 		0xE6E6FA },
    { "floral white", 		0xFFFAF0 },
    { "alice blue", 		0xF0F8FF },
    { "ghost white", 		0xF8F8FF },
    { "honeydew", 		0xF0FFF0 },
    { "ivory", 			0xFFFFF0 },
    { "azure", 			0xF0FFFF },
    { "snow", 			0xFFFAFA },
    { "black", 			0x000000 },
    { "dim gray",               0x696969 },
    { "gray", 		        0x808080 },
    { "dark gray", 	        0xA9A9A9 },
    { "silver", 		0xC0C0C0 },
    { "light gray", 	        0xD3D3D3 },
    { "gainsboro", 		0xDCDCDC },
    { "white smoke", 		0xF5F5F5 },
    { "white", 			0xFFFFFF },
    { "copper",                 0xB87333 },
    { "nvidia green",           (118 << 16) | (185 << 8) | (0 << 0) },
};

static constexpr size_t color_info_cnt = sizeof( color_info ) / sizeof( color_info[0] );

uint Layout::color( std::string name, real a )
{
    uint au = uint( a * 255.0 + 0.5 );
    if ( au > 255 ) au = 255;
    const char * name_c = name.c_str();
    for( size_t i = 0; i < color_info_cnt; i++ )
    {
        if ( strcmp( color_info[i].name, name_c ) == 0 ) return (color_info[i].rgb << 8) | (au << 0);
    }
    lassert( false, "unknown color: " + name );
    return 0;
}

void Layout::materials_init( void )
{
    // common materials
    //
    hdr->material_cnt = 0;
    uint mi = hdr->material_cnt++;
    materials[mi] = Material{ str_get( "vacuum" ),
                              1.0,              // relative_permittivity
                              1.0,              // permeability
                              0.0,              // conductivity
                              0.0,              // thermal_conductivity
                              0.0,              // mass_density
                              0.0,              // specific_heat
                              0.0,              // youngs_modulus
                              0.0,              // poissons_ratio
                              0.0 };            // thermal_expansion_coefficient
    mi = hdr->material_cnt++;
    materials[mi] = Material{ str_get( "air" ),
                              1.0006,           // relative_permittivity
                              1.0000004,        // permeability
                              0.0,              // conductivity
                              0.026,            // thermal_conductivity
                              1.1614,           // mass_density
                              1007,             // specific_heat
                              0.0,              // youngs_modulus
                              0.0,              // poissons_ratio
                              0.0 };            // thermal_expansion_coefficient
    mi = hdr->material_cnt++;
    materials[mi] = Material{ str_get( "SiO2" ),
                              4.0,              // relative_permittivity
                              1.0000010,        // permeability
                              0.0,              // conductivity
                              1.5,              // thermal_conductivity
                              2220,             // mass_density
                              745,              // specific_heat
                              0.0,              // youngs_modulus
                              0.0,              // poissons_ratio
                              4.5e-06 };        // thermal_expansion_coefficient
    mi = hdr->material_cnt++;
    materials[mi] = Material{ str_get( "Si" ),
                              11.9,             // relative_permittivity
                              1.0000010,        // permeability
                              0.0,              // conductivity
                              148,              // thermal_conductivity
                              2330,             // mass_density
                              712,              // specific_heat
                              135000000000,     // youngs_modulus
                              0.25,             // poissons_ratio
                              2.54e-06 };       // thermal_expansion_coefficient
    mi = hdr->material_cnt++;
    materials[mi] = Material{ str_get( "Si20ohm" ),
                              11.9,             // relative_permittivity
                              1.0000020,        // permeability
                              5.0,              // conductivity
                              148,              // thermal_conductivity
                              2330,             // mass_density
                              712,              // specific_heat
                              135000000000,     // youngs_modulus
                              0.25,             // poissons_ratio
                              2.54e-06 };       // thermal_expansion_coefficient
    mi = hdr->material_cnt++;
    materials[mi] = Material{ str_get( "Al" ),
                              10000.0,          // relative_permittivity
                              1.000021,         // permeability
                              38000000,         // conductivity
                              237.5,            // thermal_conductivity
                              2689,             // mass_density
                              951,              // specific_heat
                              69000000000,      // youngs_modulus
                              0.31,             // poissons_ratio
                              2.33e-05 };       // thermal_expansion_coefficient
    mi = hdr->material_cnt++;
    materials[mi] = Material{ str_get( "Cu" ),
                              10000.0,          // relative_permittivity
                              0.999991,         // permeability
                              58000000,         // conductivity
                              400,              // thermal_conductivity
                              8933,             // mass_density
                              385,              // specific_heat
                              120000000000,     // youngs_modulus
                              0.38,             // poissons_ratio
                              1.77e-05 };       // thermal_expansion_coefficient
}

uint Layout::material_set( std::string name, const Material& material )
{
    uint mi = material_get( name );
    if ( mi == NULL_I ) {
        perhaps_realloc( materials, hdr->material_cnt, max->material_cnt, 1 );
        mi = hdr->material_cnt++;
    }
    materials[mi] = material;
    return mi;
}

uint Layout::material_get( std::string name ) 
{
    uint name_i = str_get( name );

    for( uint mi = 0; mi < hdr->material_cnt; mi++ )
    {
        if ( materials[mi].name_i == name_i ) return mi;
    }

    return NULL_I;
}

void Layout::layer_set( uint layer_i, const Layer& layer )
{
    lassert( layer_i != NULL_I, "layer_i is NULL_I" );
    lassert( layer_i <= hdr->layer_cnt, "layer_i is out of range" );
    if( layer_i == hdr->layer_cnt ) {
        perhaps_realloc( layers, hdr->material_cnt, max->material_cnt, 1 );
        hdr->layer_cnt++;
    }
    layers[layer_i] = layer;
    ldout << "Setting layer_i=" << layer_i << " gdsii_layer=" << layer.gdsii_num << "\n";
}

uint Layout::layer_get( std::string name )
{
    uint name_i = str_get( name );

    for( uint li = 0; li < hdr->layer_cnt; li++ )
    {
        if ( layers[li].name_i == name_i ) return li;
    }

    return NULL_I;
}

inline bool Layout::node_is_header_footer( const Node& node ) const
{
    switch( node.kind ) 
    {
        case NODE_KIND::HEADER:
        case NODE_KIND::BGNLIB: 
        case NODE_KIND::LIBNAME: 
        case NODE_KIND::UNITS:
        case NODE_KIND::ENDLIB:
            return true;

        default:
            return false;
    }
}

inline bool Layout::node_is_gdsii( const Node& node ) const
{
    return int(node.kind) >= 0 && int(node.kind) < GDSII_KIND_CNT;
}

inline bool Layout::node_is_scalar( const Node& node ) const
{
    switch( node.kind ) 
    {
        case NODE_KIND::STR:
        case NODE_KIND::BOOL:
        case NODE_KIND::INT:
        case NODE_KIND::UINT:
        case NODE_KIND::REAL:
        case NODE_KIND::ID:
            return true;

        default:
            return false;
    }
}

inline bool Layout::node_is_string( const Node& node ) const
{
    switch( node.kind ) 
    {
        case NODE_KIND::STR:
        case NODE_KIND::ID:
            return true;

        default:
            return kind_to_datatype( node.kind ) == GDSII_DATATYPE::STRING;
    }
}

inline bool Layout::node_is_parent( const Node& node ) const
{
    if ( node_is_scalar( node ) ) return false;

    GDSII_DATATYPE datatype = kind_to_datatype( node.kind );
    return datatype != GDSII_DATATYPE::STRING && datatype != GDSII_DATATYPE::BITARRAY;
}

inline bool Layout::node_is_hier( const Node& node ) const
{
    switch( node.kind )
    {
        case NODE_KIND::BGNLIB:
        case NODE_KIND::BGNSTR:
        case NODE_KIND::BOUNDARY:
        case NODE_KIND::PATH:
        case NODE_KIND::SREF:
        case NODE_KIND::AREF:
        case NODE_KIND::TEXT:
        case NODE_KIND::NODE:
        case NODE_KIND::HIER:
            return true;

        default:
            return false;
    }
}

inline bool Layout::node_is_element( const Node& node ) const
{
    switch( node.kind )
    {
        case NODE_KIND::BOUNDARY:
        case NODE_KIND::PATH:
        case NODE_KIND::SREF:
        case NODE_KIND::AREF:
        case NODE_KIND::TEXT:
        case NODE_KIND::NODE:
            return true;

        default:
            return false;
    }
}

Layout::NODE_KIND Layout::hier_end_kind( Layout::NODE_KIND kind ) const
{
    switch( kind )
    {
        case NODE_KIND::BGNLIB:                
            return NODE_KIND::ENDLIB;

        case NODE_KIND::BGNSTR:
            return NODE_KIND::ENDSTR;

        case NODE_KIND::BOUNDARY:
        case NODE_KIND::PATH:
        case NODE_KIND::SREF:
        case NODE_KIND::AREF:
        case NODE_KIND::TEXT:
        case NODE_KIND::NODE:
            return NODE_KIND::ENDEL;

        default:
            lassert( false, "bad kind to hier_end_kind()" );
            return NODE_KIND::ENDEL;  // for compiler
    }
}

bool Layout::node_is_name( const Node& node ) const
{
    switch( node.kind )
    {
        case NODE_KIND::LIBNAME:
        case NODE_KIND::STRNAME:
        case NODE_KIND::SNAME:
        case NODE_KIND::SRFNAME:
            return true;

        default:
            return false;
    }
}

inline bool Layout::node_is_ref( const Node& node ) const
{
    switch( node.kind )
    {
        case NODE_KIND::AREF:
        case NODE_KIND::SREF:
            return true;

        default:
            return false;
    }
}

inline uint Layout::node_last_scalar_i( const Node& node ) const
{
    lassert( node_is_parent( node ), "node_last_scalar_i: node is not a parent" );    
    uint last_i = NULL_I;
    for( uint child_i = node.u.child_first_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
    {
        if ( !node_is_scalar( nodes[child_i] ) ) break;

        last_i = child_i;
    }
    return last_i;
}

inline uint Layout::node_name_i( const Node& node ) const
{
    if ( node_is_name( node ) ) {
        return node.u.s_i;
    } else if ( node_is_hier( node ) ) {
        for( uint child_i = node.u.child_first_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i ) 
        {
            uint name_i = node_name_i( nodes[child_i] );
            if ( name_i != NULL_I ) return name_i;
        }
    }
    return NULL_I;
}

inline std::string Layout::node_name( const Node& node ) const
{
    uint name_i = node_name_i( node );
    if ( name_i != NULL_I ) return &strings[name_i];
    return "";
}

inline uint Layout::node_layer( const Node& node ) const
{
    if ( node.kind == NODE_KIND::LAYER ) {
        uint int_node_i = node.u.child_first_i;
        lassert( int_node_i != NULL_I && nodes[int_node_i].kind == NODE_KIND::INT && nodes[int_node_i].sibling_i == NULL_I, "bad LAYER node" );
        return nodes[int_node_i].u.i;
    } else {
        lassert( node_is_element( node ), "node_layer: node is not a LAYER or an ELEMENT" );
        for( uint child_i = node.u.child_first_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i ) 
        {
            if ( nodes[child_i].kind == NODE_KIND::LAYER ) {
                uint int_node_i = nodes[child_i].u.child_first_i;
                lassert( int_node_i != NULL_I && nodes[int_node_i].kind == NODE_KIND::INT && nodes[int_node_i].sibling_i == NULL_I, "bad LAYER node" );
                return nodes[int_node_i].u.i;
            }
        }
        lassert( false, "could not find layer for node" );
        return NULL_I;
    }
}

bool Layout::node_has_layer( uint ni, uint layer_num, has_layer_cache_t * cache, std::string indent_str ) const
{
    lassert( ni < hdr->node_cnt, "node index ni is out of range in node_has_layer" );
    const Node& node = nodes[ni];
    if ( node.kind == NODE_KIND::LAYER ) {
        //------------------------------------------------------------
        // Finally hit a LAYER.
        //------------------------------------------------------------
        uint li = node_layer( node );
        //ldout << indent_str << "node_has_layer: LAYER=" << li << " and want layer_num=" << layer_num << "\n";
        return li == layer_num;

    } else if ( node.kind == NODE_KIND::SREF || node.kind == NODE_KIND::AREF ) {
        //------------------------------------------------------------
        // Get the SNAME and see if that struct has the layer.
        // We may have already computed this.
        //------------------------------------------------------------
        uint name_i = node_name_i( node );
        if ( name_i == NULL_I ) return false;
        auto sit = name_i_to_struct_i.find( name_i );
        lassert( sit != name_i_to_struct_i.end(), "could not find structure with name " + std::string(&strings[name_i]) );
        uint struct_i = sit->second;
        return node_has_layer( struct_i, layer_num, cache, indent_str );

    } else if ( node_is_hier( node ) ) {
        //------------------------------------------------------------
        // If this is a struct, then we may have already computed the answer.
        //------------------------------------------------------------
        if ( node.kind == NODE_KIND::BGNSTR ) {
            auto it = cache->find( ni );
            if ( it == cache->end() ) {
                (*cache)[ni] = new std::map<uint, bool>;
                it = cache->find( ni );
            }
            std::map<uint, bool>& layer_exists = *it->second;
            auto eit = layer_exists.find( layer_num );
            if ( eit != layer_exists.end() ) return eit->second;
        }

        //------------------------------------------------------------
        // Check children.
        //------------------------------------------------------------
        bool has_layer = false;
        for( uint child_i = node.u.child_first_i; !has_layer && child_i != NULL_I; child_i = nodes[child_i].sibling_i ) 
        {
            has_layer |= node_has_layer( child_i, layer_num, cache, indent_str );
        }

        if ( node.kind == NODE_KIND::BGNSTR ) {
            //------------------------------------------------------------
            // Save the answer for this struct.
            //------------------------------------------------------------
            std::map<uint, bool>& layer_exists = *(*cache)[ni];
            layer_exists[layer_num] = has_layer;
            //ldout << indent_str << "node_has_layer: struct " << node_name( nodes[ni] ) << " layer_exists[" << layer_num << "]=" << has_layer << "\n";
        }

        return has_layer;
    }
    return false;
}

inline uint Layout::node_xy_i( const Node& node ) const
{
    for( uint child_i = node.u.child_first_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i ) 
    {
        if ( nodes[child_i].kind == NODE_KIND::XY ) return child_i;
    }
    return NULL_I;
}

inline uint Layout::node_alloc( NODE_KIND kind )
{
    perhaps_realloc( nodes, hdr->node_cnt, max->node_cnt, 1 );
    uint ni = hdr->node_cnt++;
    nodes[ni].kind = kind;
    nodes[ni].sibling_i = NULL_I;
    nodes[ni].u.child_first_i = NULL_I;
    return ni;
}

inline uint Layout::node_copy( uint parent_i, uint last_i, const Layout * src_layout, uint src_i, COPY_KIND kind, real x, real y, real angle, bool reflect, bool in_flatten )
{
    const Node& src_node = src_layout->nodes[src_i];
    if ( kind == COPY_KIND::FLATTEN ) {
        //-----------------------------------------------------
        // Flattening handles certain nodes differently than hier copies.
        //-----------------------------------------------------
        switch( src_node.kind )
        {
            case NODE_KIND::STRNAME:
            {
                //-----------------------------------------------------
                // Blow off this node if we're in a struct that's
                // being flattened.
                //-----------------------------------------------------
                if ( in_flatten ) return NULL_I;
                break;
            }

            case NODE_KIND::SREF:
            case NODE_KIND::AREF:
            {
                //-----------------------------------------------------
                // First extract the instancing params.
                //-----------------------------------------------------
                const Node * src_nodes = src_layout->nodes;
                uint sname_i = NULL_I;
                uint struct_i = NULL_I;
                uint strans = 0;
                real sangle = 0.0;
                uint col_cnt = 1;
                uint row_cnt = 1;
                real xy[3][2] = { {0, 0}, {0, 0}, {0, 0} };
                for( uint src_i = src_node.u.child_first_i; src_i != NULL_I; src_i = src_nodes[src_i].sibling_i )
                {
                    const Node& child = src_nodes[src_i];
                    switch( child.kind )
                    {
                        case NODE_KIND::SNAME:
                            lassert( sname_i == NULL_I, "SREF/AREF has duplicate SNAME child" );
                            sname_i = child.u.s_i;
                            auto it = src_layout->name_i_to_struct_i.find( sname_i );
                            lassert( it != src_layout->name_i_to_struct_i.end(), "SREF/AREF SNAME " + std::string(&src_layout->strings[sname_i]) + " does not denote a known struct" );
                            struct_i = it->second;
                            break;

                        case NODE_KIND::STRANS:
                            strans = child.u.u;
                            break;

                        case NODE_KIND::ANGLE:
                            sangle = child.u.r;
                            break;
                            
                        case NODE_KIND::COLROW:
                        {
                            lassert( src_node.kind == NODE_KIND::AREF, "COLROW not allowed for an SREF" );
                            uint gchild = child.u.child_first_i;
                            lassert( gchild != NULL_I, "COLROW has no COL value" );
                            lassert( src_nodes[gchild].kind == NODE_KIND::INT, "COLROW COL is not an INT" );
                            col_cnt = src_nodes[gchild].u.i;
                            lassert( col_cnt > 0, "COLROW COL must be non-zero" );

                            gchild = child.sibling_i;
                            lassert( gchild != NULL_I, "COLROW has no ROW value" );
                            lassert( src_nodes[gchild].kind == NODE_KIND::INT, "COLROW ROW is not an INT" );
                            row_cnt = src_nodes[gchild].u.i;
                            lassert( row_cnt > 0, "COLROW ROW must be non-zero" );
                            break;
                        }
                             
                        case NODE_KIND::XY:
                        {
                            uint i = 0;
                            for( uint gchild_i = child.u.child_first_i; gchild_i != NULL_I; gchild_i = src_nodes[gchild_i].sibling_i )
                            {
                                lassert( i == 0 || src_node.kind == NODE_KIND::AREF, "SREF may not have more than 2 XY coords" );
                                lassert( i < 6, "AREF may not have more than 6 XY coords" );
                                xy[i>>1][i&1] = real(src_nodes[gchild_i].u.i) / src_layout->gdsii_units_user + 0.5; 
                                i++;
                            }
                            lassert( i == 2 || i == 6, "wrong number of XY coords for " + str(src_node.kind) );
                            break;
                        }

                        default:
                        {
                            lassert( false, "unexpected SREF/AREF child node: " + str(child.kind) );
                            break;
                        }
                    }
                }

                //-----------------------------------------------------
                // Now the fun part.
                //-----------------------------------------------------
                in_flatten = true;
                lassert( struct_i != NULL_I, "SREF/AREF has no SNAME" );
                for( uint r = 0; r < row_cnt; r++ )
                {
                    for( uint c = 0; c < col_cnt; c++ )
                    {
                        //-----------------------------------------------------
                        // Copy the structure's children with new transformation parameters.
                        //-----------------------------------------------------
                        real inst_x = x + xy[0][0];  
                        real inst_y = y + xy[1][0];
                        real inst_angle = angle + sangle;
                        real inst_reflect = reflect;
                        for( uint child_i = src_nodes[struct_i].u.child_first_i; child_i != NULL_I; child_i = src_nodes[child_i].sibling_i )
                        {
                        }
                    }
                }
                return NULL_I;
            }

            default:
            {
                //-----------------------------------------------------
                // We can copy this node.
                //-----------------------------------------------------
                break;
            }
        }
    }

    //-----------------------------------------------------
    // Copy the main node.
    //-----------------------------------------------------
    uint dst_i = node_alloc( src_node.kind );
    bool src_is_parent = src_layout->node_is_parent( src_node );
    if ( src_is_parent ) {
        nodes[dst_i].u.child_first_i = NULL_I;
        if ( kind != COPY_KIND::ONE ) {
            //-----------------------------------------------------
            // Copy children.
            //-----------------------------------------------------
            uint dst_prev_i = NULL_I;
            for( src_i = src_node.u.child_first_i; src_i != NULL_I; src_i = src_layout->nodes[src_i].sibling_i )
            {
                if ( kind != COPY_KIND::DEEP && kind != COPY_KIND::FLATTEN && !node_is_scalar( src_layout->nodes[src_i] ) ) break;

                dst_prev_i = node_copy( dst_i, dst_prev_i, src_layout, src_i, kind, x, y, angle, reflect );
            }
        }
    } else {
        if ( node_is_string( src_node ) ) {
            nodes[dst_i].u.s_i = str_get( &src_layout->strings[src_node.u.s_i] );
        } else {
            nodes[dst_i].u = src_node.u;
        }
    }
    
    if ( last_i != NULL_I ) {
        lassert( nodes[last_i].sibling_i == NULL_I, "node_copy: nodes[last_i] sibling_i is already set" ); 
        nodes[last_i].sibling_i = dst_i;
    } else if ( parent_i != NULL_I ) {
        lassert( nodes[parent_i].u.child_first_i == NULL_I, "node_copy: nodes[parent_i] child_first_i is already set but last_i is NULL_I" ); 
        nodes[parent_i].u.child_first_i = dst_i;
    }
    return dst_i;
}

void Layout::node_timestamp( Node& node )
{
    lassert( node.kind == NODE_KIND::BGNLIB || node.kind == NODE_KIND::BGNSTR, "node_timestamp: node must be BGNLIB or BGNSTR" );
    lassert( node.u.child_first_i == NULL_I, "node_timestamp found node with children already" );

    time_t t;
    time( &t );
    struct tm * tm = gmtime( &t );

    uint prev_i = NULL_I;
    for( uint i = 0; i < 2; i++ ) 
    {
        uint ni = node_alloc( NODE_KIND::INT );
        nodes[ni].u.i = tm->tm_year;
        if ( i == 0 ) {
            node.u.child_first_i = ni;
        } else {
            nodes[prev_i].sibling_i = ni;
        }
        prev_i = ni;

        ni = node_alloc( NODE_KIND::INT );
        nodes[ni].u.i = tm->tm_mon;
        nodes[prev_i].sibling_i = ni;
        prev_i = ni;

        ni = node_alloc( NODE_KIND::INT );
        nodes[ni].u.i = tm->tm_mday;
        nodes[prev_i].sibling_i = ni;
        prev_i = ni;

        ni = node_alloc( NODE_KIND::INT );
        nodes[ni].u.i = tm->tm_hour;
        nodes[prev_i].sibling_i = ni;
        prev_i = ni;

        ni = node_alloc( NODE_KIND::INT );
        nodes[ni].u.i = tm->tm_min;
        nodes[prev_i].sibling_i = ni;
        prev_i = ni;

        ni = node_alloc( NODE_KIND::INT );
        nodes[ni].u.i = tm->tm_sec;
        nodes[prev_i].sibling_i = ni;
        prev_i = ni;
    }
}

inline Layout::real Layout::width( void ) const
{
    return 0;
}

inline Layout::real Layout::length( void ) const
{
    return 0;
}

inline Layout::real Layout::height( void ) const
{
    return 0;
}

uint Layout::start_library( std::string libname, real units_user, real units_meters )
{
    lassert( hdr->root_i == NULL_I, "starting a library when layout is not empty" );
    
    uint ni = node_alloc( Layout::NODE_KIND::HIER );   // for file level (one node)
    hdr->root_i = ni;
    uint prev_i = ni;

    ni = node_alloc( Layout::NODE_KIND::HEADER );
    nodes[prev_i].u.child_first_i = ni;
    uint ni2 = node_alloc( Layout::NODE_KIND::INT );
    nodes[ni].u.child_first_i = ni2;
    nodes[ni2].u.i = 5;
    prev_i = ni;

    ni = node_alloc( Layout::NODE_KIND::BGNLIB );
    node_timestamp( nodes[ni] );
    nodes[prev_i].sibling_i = ni;
    prev_i = node_last_scalar_i( nodes[ni] );
    
    ni = node_alloc( Layout::NODE_KIND::LIBNAME );
    nodes[prev_i].sibling_i = ni;
    nodes[ni].u.s_i = str_get( libname );
    prev_i = ni;

    ni = node_alloc( Layout::NODE_KIND::UNITS );
    gdsii_units_user = units_user;
    gdsii_units_meters = units_meters;
    nodes[prev_i].sibling_i = ni;
    ni2 = node_alloc( NODE_KIND::REAL );
    nodes[ni].u.child_first_i = ni2;
    nodes[ni2].u.r = units_user;
    prev_i = ni2;
    ni2 = node_alloc( NODE_KIND::REAL );
    nodes[prev_i].sibling_i = ni2;
    nodes[ni2].u.r = units_meters;
    prev_i = ni;

    return prev_i;
}

uint Layout::inst_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, real x, real y, std::string name )
{
    return inst_layout( parent_i, last_i, src_layout, src_struct_name, x, y, 0, hdr->layer_cnt-1, name );
}

uint Layout::inst_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, real x, real y, uint dst_layer_first, uint dst_layer_last, std::string name )
{
    //-----------------------------------------------------
    // Initialize our cache.
    //-----------------------------------------------------
    auto cache = new has_layer_cache_t;

    //-----------------------------------------------------
    // Do the instancing separately for each dst_layer.
    // Use a unique name for each.
    //-----------------------------------------------------
    for( uint i = dst_layer_first; i <= dst_layer_last; i++ )
    {
        //-----------------------------------------------------
        // Recursively find and copy all STRuctures, but 
        // only ones that have the desired src_layer_num.
        //-----------------------------------------------------
        std::string inst_name = std::to_string( i ) + "_" + name;
        uint src_layer_num = layers[i].gdsii_num;
        ldout << "inst_layout: dst_layer=" << i << " src_layer=" << layers[i].gdsii_num << " x_off=" << x << " y_off=" << y << " inst_name=" << inst_name << "\n";
        uint inst_last_i = inst_layout_node( parent_i, last_i, src_layout, src_struct_name, src_layout->hdr->root_i, src_layer_num, i, cache, inst_name );
        if ( inst_last_i != NULL_I ) {
            last_i = inst_last_i;
        }

        //-----------------------------------------------------
        // Find the top struct in the new layout.
        // It must exist.
        // Add top struct to the list of top-level structures.
        // Translate all XY nodes only in this structure.
        //-----------------------------------------------------
        std::string dst_top_struct_name = inst_name + "_" + src_struct_name;
        uint dst_top_struct_name_i = str_get( dst_top_struct_name );
        lassert( name_i_to_struct_i.find( dst_top_struct_name_i ) != name_i_to_struct_i.end(), "could not find top struct with name " + dst_top_struct_name );
        uint dst_top_struct_i = name_i_to_struct_i[dst_top_struct_name_i];
        Node& dst_top_node = nodes[dst_top_struct_i];

        perhaps_realloc( top_insts, hdr->top_inst_cnt, max->top_inst_cnt, 1 );
        uint ii = hdr->top_inst_cnt++;
        top_insts[ii] = TopInstInfo{ dst_top_struct_i, x, y };
    }
    return last_i;
}

uint Layout::flatten_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, std::string dst_struct_name )
{
    //------------------------------------------------------------
    // Locate the the src struct.
    //------------------------------------------------------------
    uint src_struct_name_i = src_layout->str_find( src_struct_name );
    auto it = src_layout->name_i_to_struct_i.find( src_struct_name_i );
    lassert( it != src_layout->name_i_to_struct_i.end(), "no src struct with the name " + src_struct_name );
    uint src_struct_i = it->second;

    //------------------------------------------------------------
    // Recursively copy the src_struct.
    //------------------------------------------------------------
    return node_copy( parent_i, last_i, src_layout, src_struct_i, COPY_KIND::FLATTEN );
}

uint Layout::inst_layout_node( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, uint src_i, uint src_layer_num, uint dst_layer_num, 
                               has_layer_cache_t * cache, std::string name, std::string indent_str )
{
    //-----------------------------------------------------
    // See if node should be copied.
    // We do this proactively.
    //-----------------------------------------------------
    bool do_copy = true;
    const Node& src_node = src_layout->nodes[src_i];
    NODE_KIND src_kind = src_node.kind;
    if ( (src_kind == NODE_KIND::BGNSTR || src_kind == NODE_KIND::SREF || src_kind == NODE_KIND::AREF) ) {
        //-----------------------------------------------------
        // If struct or ref does not have desired layer, then bail.
        //-----------------------------------------------------
        if ( !src_layout->node_has_layer( src_i, src_layer_num, cache, indent_str ) ) {
            //ldout << indent_str << str(src_node.kind) << " does NOT use src_layer=" << std::to_string(src_layer_num) << "\n";
            return NULL_I;
        }
        //ldout << indent_str << str(src_node.kind) << " uses src_layer=" << std::to_string(src_layer_num) << "\n";
    }

    if ( !src_layout->node_is_ref( src_node ) && 
          src_layout->node_is_element( src_node ) && 
          src_layout->node_layer( src_node ) != src_layer_num ) {
        do_copy = false;
        ldout << indent_str << "    src_layer=" << src_layout->node_layer( src_node ) << " desired_src_layer=" << src_layer_num << "\n";

    } else if ( src_layout->node_is_header_footer( src_node ) && src_kind != NODE_KIND::BGNLIB ) {
        //-----------------------------------------------------
        // Blow this off.
        //-----------------------------------------------------
        do_copy = false;
    }

    ldout << indent_str << str(src_kind) << ": do_copy=" << do_copy << "\n";
    if ( !do_copy ) return NULL_I;

    if ( src_kind != NODE_KIND::BGNLIB && src_kind != NODE_KIND::HIER ) {
        //-----------------------------------------------------
        // Copy this node and its scalar children.
        //-----------------------------------------------------
        uint dst_i = node_copy( parent_i, last_i, src_layout, src_i, COPY_KIND::SCALAR_CHILDREN );    // don't copy non-children yet
        lassert( dst_i != NULL_I, "node_copy should have returned a non-null node index" );
        last_i = dst_i;
        if ( src_kind == NODE_KIND::LAYER ) {
            //-----------------------------------------------------
            // Change LAYER from src_layer_num to dst_layer_num.
            //-----------------------------------------------------
            uint dst_int_node_i = nodes[dst_i].u.child_first_i;
            lassert( dst_int_node_i != NULL_I, "LAYER node has no layer number" ); 
            Node& dst_int_node = nodes[dst_int_node_i];
            lassert( dst_int_node.kind == NODE_KIND::INT && dst_int_node.u.i == src_layer_num && dst_int_node.sibling_i == NULL_I, "unexpected layer_num in LAYER node" );
            dst_int_node.u.i = dst_layer_num;
            ldout << indent_str << "    layer change: " << src_layer_num << " => " << dst_layer_num << "\n";

        } else if ( src_kind == NODE_KIND::STRNAME || src_kind == NODE_KIND::SNAME ) {
            //-----------------------------------------------------
            // Prepend "name" to struct name to uniquify.
            //-----------------------------------------------------
            std::string str_name = name + std::string("_") + &strings[nodes[dst_i].u.s_i];
            nodes[dst_i].u.s_i = str_get( str_name );
        }

        if ( src_layout->node_is_parent( src_node ) ) {
            //-----------------------------------------------------
            // Recursively copy children.
            //-----------------------------------------------------
            uint src_prev_i = src_layout->node_last_scalar_i( src_node );
            uint dst_prev_i = node_last_scalar_i( nodes[dst_i] );
            uint src_child_i = (src_prev_i == NULL_I) ? src_node.u.child_first_i : src_layout->nodes[src_prev_i].sibling_i; 
            for( ; src_child_i != NULL_I; src_child_i = src_layout->nodes[src_child_i].sibling_i )
            {
                uint dst_child_i = inst_layout_node( src_i, dst_prev_i, src_layout, src_struct_name, src_child_i, src_layer_num, dst_layer_num, cache, name, indent_str + "    " );
                if ( dst_child_i != NULL_I ) {
                    if ( dst_prev_i == NULL_I ) nodes[dst_i].u.child_first_i = dst_child_i;
                    dst_prev_i = dst_child_i;
                }
            }
        }

        if ( src_kind == NODE_KIND::BGNSTR ) {
            //-----------------------------------------------------
            // Record struct by name.
            //-----------------------------------------------------
            uint struct_name_i = node_name_i( nodes[dst_i] );
            lassert( struct_name_i != NULL_I, "could not find struct name"  );
            name_i_to_struct_i[struct_name_i] = dst_i;
        }
    } else {
        //-----------------------------------------------------
        // Skip BGNLIB/HIER and process non-scalar children.
        //-----------------------------------------------------
        for( uint src_child_i = src_node.u.child_first_i; src_child_i != NULL_I; src_child_i = src_layout->nodes[src_child_i].sibling_i )
        {
            if ( !node_is_scalar( src_layout->nodes[src_child_i] ) ) {
                uint dst_child_i = inst_layout_node( src_i, last_i, src_layout, src_struct_name, src_child_i, src_layer_num, dst_layer_num, cache, name, indent_str + "    " );
                if ( dst_child_i != NULL_I ) last_i = dst_child_i;
            } else {
                ldout << indent_str << "    " << "skipping " << str(src_layout->nodes[src_child_i].kind) << "\n";
            }
        }
    }

    return last_i;
}

void Layout::fill_dielectrics( void )
{
    //-----------------------------------------------------
    // For each layer:
    //-----------------------------------------------------
    for( uint32_t i = 0; i < hdr->layer_cnt; i++ )
    {
        //-----------------------------------------------------
        // Build a 2D kd-tree, which is a space-divided binary tree.
        // It will make it very easy to figure out where there
        // are empty spaces.
        //-----------------------------------------------------

        //-----------------------------------------------------
        // Walk the kd-tree down to the leaf level.
        // For leaves that have a triangle or rectangle, it's
        // trivial to determine the remaining empty space.
        // For empty leaves, it's obviously the entire rectangle.
        //-----------------------------------------------------
    }
}

void Layout::fill_material( uint material_i, real x, real y, real z, real w, real l, real h )
{
    //-----------------------------------------------------
    // For now, we require that the space consumed by
    // the fill material is outside the space consumed by
    // the layout.  Thus we need only add a box of the 
    // material at origin [x,y,z] and dimensions [w,l,h].
    // Simple.
    //-----------------------------------------------------
}

void Layout::finalize_top_struct( uint parent_i, uint last_i, std::string top_name )
{
    //-----------------------------------------------------
    // Create one final struct that instantiates all top_insts.
    //-----------------------------------------------------
    lassert( parent_i == NULL_I, "finalize parent_i must be NULL_I for now" );
    uint bgnstr_i = node_alloc( NODE_KIND::BGNSTR );
    nodes[last_i].sibling_i = bgnstr_i;
    node_timestamp( nodes[bgnstr_i] );
    uint prev_i = node_last_scalar_i( nodes[bgnstr_i] );

    uint strname_i = node_alloc( NODE_KIND::STRNAME );
    nodes[prev_i].sibling_i = strname_i;
    nodes[strname_i].u.s_i = str_get( top_name );
    prev_i = strname_i;

    for( size_t i = 0; i < hdr->top_inst_cnt; i++ )
    {
        uint sref_i = node_alloc( NODE_KIND::SREF );
        nodes[prev_i].sibling_i = sref_i;
        prev_i = sref_i;

        const TopInstInfo& info = top_insts[i];

        uint sname_i = node_alloc( NODE_KIND::SNAME );
        nodes[sname_i].u.s_i = node_name_i( nodes[info.struct_i] );
        nodes[sref_i].u.child_first_i = sname_i;

        uint xy_i = node_alloc( NODE_KIND::XY );
        nodes[sname_i].sibling_i = xy_i;

        uint x_i = node_alloc( NODE_KIND::INT );
        nodes[xy_i].u.child_first_i = x_i;
        nodes[x_i].u.i = int( info.x / gdsii_units_user );
        uint y_i = node_alloc( NODE_KIND::INT );
        nodes[x_i].sibling_i = y_i;
        nodes[y_i].u.i = int( info.y / gdsii_units_user );
    }
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
            lassert( old_max_cnt != NULL_I, "old_max_cnt should have been reasonable" );
            max_cnt = NULL_I;
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
        lassert( 0, "hdr->version does not match VERSION=" + std::to_string(VERSION) + ", got " + std::to_string(hdr->version) );
        return false;
    }
    max = aligned_alloc<Header>( 1 );
    memcpy( max, hdr, sizeof( Header ) );
    _uread( strings,     char,        hdr->char_cnt );
    _uread( materials,   Material,    hdr->material_cnt );
    _uread( layers,      Layer,       hdr->layer_cnt );
    _uread( nodes,       Node,        hdr->node_cnt );
    _uread( top_insts,   TopInstInfo, hdr->top_inst_cnt );

    return true;
}

bool Layout::layout_write( std::string layout_path ) 
{
    cmd( "rm -f " + layout_path );
    int fd = open( layout_path.c_str(), O_CREAT|O_WRONLY|O_TRUNC|O_SYNC|S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP );
    if ( fd < 0 ) std::cout << "open() for write error: " << strerror( errno ) << "\n";
    lassert( fd >= 0, "could not open() file " + layout_path + " for writing - open() error: " + strerror( errno ) );

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
                lassert( 0, "could not write() file " + layout_path + " - write() error: " + strerror( errno ) ); \
            } \
            _byte_cnt -= _this_byte_cnt; \
            _addr     += _this_byte_cnt; \
        } \
    } \

    _uwrite( hdr,         1                  * sizeof(hdr[0]) );
    _uwrite( strings,     hdr->char_cnt      * sizeof(strings[0]) );
    _uwrite( materials,   hdr->material_cnt  * sizeof(materials[0]) );
    _uwrite( layers,      hdr->layer_cnt     * sizeof(layers[0]) );
    _uwrite( nodes,       hdr->node_cnt      * sizeof(nodes[0]) );
    _uwrite( top_insts,   hdr->top_inst_cnt  * sizeof(top_insts[0] ));

    fsync( fd ); // flush
    close( fd );
    cmd( "chmod +rw " + layout_path );

    return true;
}

bool Layout::gdsii_read( std::string file, bool count_only )
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
    uint ni = NULL_I;
    if ( !count_only ) {
        ni = node_alloc( NODE_KIND::HIER );
        hdr->root_i = ni;
    }
    uint prev_i = NULL_I;
    gdsii_rec_cnt = 0;
    uint struct_i = NULL_I;
    while( nnn < nnn_end )
    {
        uint child_i;
        if ( !gdsii_read_record( child_i, struct_i, count_only ) ) return false;

        if ( !count_only ) {
            if ( prev_i == NULL_I ) {
                nodes[ni].u.child_first_i = child_i;
            } else {
                nodes[prev_i].sibling_i = child_i;
            }
            prev_i = child_i;
        }

        if ( gdsii_last_kind == NODE_KIND::ENDLIB ) break;
    }

    if ( count_only ) {
        std::cout << "gdsii_rec_cnt=" << gdsii_rec_cnt << "\n";
    }
    return true;
}

static inline bool is_gdsii_allowed_char( char c ) { return isprint( c ) && c != '"' && c != ','; }

bool Layout::gdsii_read_record( uint& ni, uint struct_i, bool count_only )
{
    //------------------------------------------------------------
    // Parse record header.
    //------------------------------------------------------------
    lassert( (nnn + 4) <= nnn_end, "unexpected end of gdsii file rec_cnt=" + std::to_string(gdsii_rec_cnt) );
    uint32_t       byte_cnt = ( nnn[0] << 8 ) | nnn[1];
    NODE_KIND      kind     = NODE_KIND( nnn[2] );
    if ( false && (gdsii_rec_cnt % 1000000) == 0 ) std::cout << gdsii_rec_cnt << ": " << kind << "\n";
    GDSII_DATATYPE datatype = GDSII_DATATYPE( nnn[3] );
    lassert( byte_cnt >= 4, std::to_string(gdsii_rec_cnt) + ": gdsii record byte_cnt must be at least 4, byte_cnt=" + std::to_string(byte_cnt) + " kind=" + str(kind) );
    ldout << hdr->node_cnt << ": " << str(kind) << " " << str(datatype) << " byte_cnt=" << std::to_string(byte_cnt) << "\n";
    byte_cnt -= 4;
    nnn += 4;
    lassert( uint32_t(kind) < GDSII_KIND_CNT, std::to_string(gdsii_rec_cnt) + ": bad gdsii record kind " + std::to_string(uint32_t(kind)) );
    lassert( kind_to_datatype( kind ) == datatype, 
             std::to_string(gdsii_rec_cnt) + ": datatype=" + str(datatype) + " does not match expected datatype=" + 
             str( kind_to_datatype( kind ) ) + " for record kind " + str(kind) );
    lassert( (nnn + byte_cnt) <= nnn_end, "unexpected end of gdsii file" );
  
    //------------------------------------------------------------
    // Create a GDSII node.
    //------------------------------------------------------------
    if ( !count_only ) {
        ni = node_alloc( kind );
    } else {
        ni = 0;
    }
    nodes[ni].kind = kind;

    //------------------------------------------------------------
    // Parse payload.
    //------------------------------------------------------------
    bool is_hier = node_is_hier( nodes[ni] );
    uint prev_i = NULL_I;
    switch( datatype )
    {
        case GDSII_DATATYPE::NO_DATA:
        {
            lassert( byte_cnt == 0, "NO_DATA gdsii datatype should have no payload" );
            break;
        }

        case GDSII_DATATYPE::BITARRAY:
        {
            lassert( !is_hier, "BITARRAY not allowed for hier nodes" );
            lassert( byte_cnt == 2, "BITARRAY gdsii datatype should have 2-byte payload" );
            if ( !count_only ) nodes[ni].u.u = (nnn[1] << 8) | nnn[0];
            break;
        }

        case GDSII_DATATYPE::STRING:
        {
            lassert( !is_hier, "STRING not allowed for hier nodes" );
            char c[1024];
            lassert( byte_cnt <= (sizeof(c)-1), "STRING too big byte_cnt=" + std::to_string(byte_cnt) );
            if ( byte_cnt > 0 ) memcpy( c, nnn, byte_cnt );
            c[byte_cnt] = '\0';
            for( int i = byte_cnt-1; i >= 0; i-- ) 
            {
                if ( is_gdsii_allowed_char( c[i] ) ) break;
                c[i] = '\0';
            }
            std::string s = std::string( c );
            ldout << "    " << s << "\n";
            if ( !count_only ) nodes[ni].u.s_i = str_get( s );
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
            lassert( (cnt*datum_byte_cnt) == byte_cnt, "datum_byte_cnt does not divide evenly" );
            uint8_t * uuu = nnn;
            for( uint i = 0; i < cnt; i++, uuu += datum_byte_cnt )
            {
                uint child_i = NULL_I;
                if ( !count_only ) {
                    child_i = node_alloc( NODE_KIND::INT );
                    if ( i == 0 ) {
                        nodes[ni].u.child_first_i = child_i;
                    } else {
                        nodes[prev_i].sibling_i = child_i;
                    }
                }
                if ( datatype == GDSII_DATATYPE::INTEGER_2 || datatype == GDSII_DATATYPE::INTEGER_4 ) {
                    int64_t vi = (uuu[0] << 8) | uuu[1];
                    if ( datatype == GDSII_DATATYPE::INTEGER_4 ) {
                        vi = (vi << 16) | (uuu[2] << 8) | uuu[3];
                    }
                    if ( (vi & 0x80000000) != 0 ) {
                        vi = (datatype == GDSII_DATATYPE::INTEGER_2) ? -(0x10000 - vi) : -(0x100000000LL - vi);
                    }
                    if ( !count_only ) {
                        nodes[child_i].kind = NODE_KIND::INT;
                        nodes[child_i].u.i = vi;
                    }
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
                    if ( !count_only ) {
                        nodes[child_i].kind = NODE_KIND::REAL;
                        real rexp = 4.0*(exp-64) - 8*(datum_byte_cnt-1);
                        nodes[child_i].u.r = sign * double(ifrac) * std::pow( 2.0, rexp );
                        ldout << "sign=" << ((sign < 0.0) ? "1" : "0") << " exp=" << exp << " rexp=" << rexp << " ifrac=" << ifrac << " r=" << nodes[child_i].u.r << "\n";
                        if ( kind == NODE_KIND::UNITS ) {
                            if ( i == 0 ) gdsii_units_user   = nodes[child_i].u.r;
                            if ( i == 1 ) gdsii_units_meters = nodes[child_i].u.r;
                        }
                    }
                }
                prev_i = child_i;
            }
            break;
        }

        default:
        {
            lassert( false, "something is wrong" );
            break;
        }
    }

    gdsii_rec_cnt++;
    gdsii_last_kind = kind;
    nnn += byte_cnt;

    if ( !count_only && kind == NODE_KIND::STRNAME ) {
        // record name -> struct mapping
        uint name_i = node_name_i( nodes[ni] );
        lassert( name_i != NULL_I, "STRNAME should have had a string" );
        lassert( struct_i != NULL_I, "STRNAME found outside a structure" );
        name_i_to_struct_i[name_i] = struct_i;
    }

    if ( is_hier ) {
        if ( !count_only && kind == NODE_KIND::BGNSTR ) {
            // now inside a structure
            lassert( struct_i == NULL_I, "nested BGNSTR is not allowed" );
            struct_i = ni;
        }

        // recurse for other children
        for( ;; ) 
        {
            uint child_i;
            if ( !gdsii_read_record( child_i, struct_i, count_only ) ) return false;
            NODE_KIND kind = nodes[child_i].kind;
            if ( kind == NODE_KIND::ENDEL || kind == NODE_KIND::ENDSTR || kind == NODE_KIND::ENDLIB ) break;
            if ( !count_only ) {
                if ( prev_i == NULL_I ) {
                    nodes[ni].u.child_first_i = child_i;
                } else {
                    nodes[prev_i].sibling_i = child_i;
                }
                prev_i = child_i;
            }
        }
    }

    if ( count_only ) nodes[ni].kind = kind;
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

    ldout << "Writing " << gdsii_path << " root_i=" << hdr->root_i << "\n";
    if ( hdr->root_i != NULL_I ) gdsii_write_record( hdr->root_i );

    gdsii_flush( true );
    fsync( gdsii_fd ); // flush
    close( gdsii_fd );
    cmd( "chmod +rw " + gdsii_path );

    delete gdsii_buff;
    gdsii_buff = nullptr;

    return true;
}

void Layout::gdsii_write_record( uint ni, std::string indent_str )
{
    const Node& node = nodes[ni];
    ldout << indent_str << str(node.kind) << " ni=" << ni << "\n";
    if ( node_is_gdsii( node ) ) {
        uint8_t bytes[64*1024];
        uint    byte_cnt = 2;   // fill in byte_cnt later

        GDSII_DATATYPE datatype = kind_to_datatype( node.kind );
        bytes[byte_cnt++] = int(node.kind);
        bytes[byte_cnt++] = int(datatype);

        uint child_i = NULL_I;

        switch( datatype )
        {
            case GDSII_DATATYPE::NO_DATA:
            {
                break;
            }

            case GDSII_DATATYPE::BITARRAY:
            {
                ldout << indent_str << "    " << node.u.u << "\n";
                bytes[byte_cnt++] = node.u.u & 0xff;
                bytes[byte_cnt++] = (node.u.u >> 8) & 0xff;
                break;
            }

            case GDSII_DATATYPE::STRING:
            {
                ldout << indent_str << "    " << std::string(&strings[node.u.s_i]) << "\n";
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
                for( child_i = node.u.child_first_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
                {
                    if ( nodes[child_i].kind != NODE_KIND::INT && nodes[child_i].kind != NODE_KIND::REAL ) break; // skip GDSII children
                    gdsii_write_number( bytes, byte_cnt, child_i, datatype, indent_str + "    " );
                }
                break;
            }

            default:
            {
                lassert( false, "bad datatype in gdsii_write_record" );
                break;
            }
        }

        // record byte count
        bytes[0] = (byte_cnt >> 8) & 0xff;
        bytes[1] = byte_cnt & 0xff;

        // transfer bytes to buffer
        gdsii_write_bytes( bytes, byte_cnt );

        if ( node_is_hier( node ) ) {
            // recurse for rest of children
            if ( child_i == NULL_I ) child_i = node.u.child_first_i;
            for( ; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
            {
                gdsii_write_record( child_i, indent_str + "    " );
            }

            NODE_KIND end_kind = hier_end_kind( node.kind );
            ldout << indent_str << str(end_kind) << "\n";
            bytes[0] = 0;
            bytes[1] = 4;
            bytes[2] = uint8_t(end_kind);
            bytes[3] = uint8_t(GDSII_DATATYPE::NO_DATA);
            gdsii_write_bytes( bytes, 4 );
        }

    } else if ( node.kind == NODE_KIND::HIER ) {
        // assume file wrapper, just loop through children
        for( uint child_i = node.u.child_first_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
        {
            gdsii_write_record( child_i, indent_str + "    " );
        }

    } else {
        lassert( false, "ignoring node kind " + str(node.kind) );
    }
}

void Layout::gdsii_write_number( uint8_t * bytes, uint& byte_cnt, uint ni, GDSII_DATATYPE datatype, std::string indent_str )
{
    switch( datatype )
    {
        case GDSII_DATATYPE::INTEGER_2:
        case GDSII_DATATYPE::INTEGER_4:
        {
            int32_t i = nodes[ni].u.i;
            ldout << indent_str << i << "\n";
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
            ldout << indent_str << nodes[ni].u.r << "\n";
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
            //ldout << indent_str << "mid: rexp " << rexp << " ifrac=" << ifrac << "\n";
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
                //ldout << indent_str << "bytes[" << j << "]=" << int(bytes[byte_cnt-datum_byte_cnt+j]) << "\n";
            }
            ldout << indent_str << "    sign=" << sign << " exp=" << int(exp) << " rexp=" << rexp << " ifrac=" << ifrac << " r=" << nodes[ni].u.r << "\n";
            break;
        }

        default:
        {
            lassert( false, "bad datatype in gdsii_write_number()" );
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

void Layout::gdsii_flush( bool for_end_of_file )
{
    if ( gdsii_buff_byte_cnt == 0 ) return;

    if ( for_end_of_file ) {
        // pad to multiple of 2048 bytes to keep some software happy
        //
        uint page_bytes_left = 2048 - (gdsii_buff_byte_cnt % 2048);
        for( uint i = 0; i < page_bytes_left; i++ )
        {
            gdsii_buff[gdsii_buff_byte_cnt++] = 0;
        }
    }

    if ( ::write( gdsii_fd, gdsii_buff, gdsii_buff_byte_cnt ) <= 0 ) { 
        close( gdsii_fd ); 
        lassert( false, std::string("could not write() gdsii file - write() error: ") + strerror( errno ) );
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
    lassert( nodes[hdr->root_i].kind == NODE_KIND::HIER, "aedt_read() root node should have been a HIER node, got " + str(nodes[hdr->root_i].kind) );
    return true;
}

bool Layout::aedt_read_expr( uint& ni )
{
    ni = node_alloc( NODE_KIND::HIER ); // will change

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
            lassert( 0, "unable to parse an expression: std::string, number, or id " + surrounding_lines( nnn_save, nnn_end ) );
        }
        ldout << "ID START " << std::string(&strings[id_i]) << "\n";
        if ( id_i == aedt_begin_str_i ) {
            uint name_i;
            if ( !aedt_read_expr( name_i ) ) return false;             // STR node
            lassert( nodes[name_i].kind == NODE_KIND::STR, "$begin not followed by std::string" );
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
                        lassert( end_str_i == nodes[name_i].u.s_i, "$end id does not match $begin id " + surrounding_lines( nnn, nnn_end ) );
                        break;
                    }
                }

                uint child_i;
                if ( !aedt_read_expr( child_i ) ) return false;
                lassert( child_i != 0, "aedt_read_expr() could not read expression within BEGIN" );
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

        uint ai = node_alloc( NODE_KIND::ASSIGN );
        nodes[ai].u.child_first_i = ni;
        nodes[ni].sibling_i = rhs_i;
        ni = ai;

    } else if ( ch == '(' || ch == '[' ) {
        ldout << std::string( (ch == '(') ? "CALL\n" : "SLICE\n" );
        uint id_i = ni;
        ni = node_alloc( (ch == '(') ? NODE_KIND::CALL : NODE_KIND::SLICE );
        nodes[ni].u.child_first_i = id_i;
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
            lassert( arg_i != 0, "could not read expression in CALL or SLICE" );
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
    return true;
}

void Layout::aedt_write_expr( std::ofstream& out, uint ni, std::string indent_str )
{
    lassert( ni != NULL_I, "aedt_write_expr: bad node index" );
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
            for( child_i = nodes[child_i].sibling_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
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
            if ( child_i != NULL_I ) {
                aedt_write_expr( out, child_i, "" );
                out << ":";
                bool have_one = false;
                for( child_i = nodes[child_i].sibling_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
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
            if ( id_i != NULL_I && nodes[id_i].kind == NODE_KIND::STR ) {
                out << "$begin '" << std::string(&strings[nodes[id_i].u.s_i]) << "'";
                child_i = nodes[id_i].sibling_i;
            } else {
                out << "$begin 'FILE'";
                child_i = id_i;
                id_i = NULL_I;
            }
            for( ; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
            {
                aedt_write_expr( out, child_i, indent_str + "\t" );
            }
            out << indent_str;
            if ( id_i != NULL_I ) {
                out << "$end '" << std::string(&strings[nodes[id_i].u.s_i]) << "'";
            } else {
                out << "$end 'FILE'\n";
            }
            break;
        }

        default:
        {
            if ( node_is_gdsii( node ) ) {
                GDSII_DATATYPE datatype = kind_to_datatype( node.kind );
                if ( node_is_hier( node ) ) {
                    out << "$begin '" << node.kind << "'";
                }
                uint child_i = NULL_I;
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
                        for( child_i = node.u.child_first_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
                        {
                            if ( nodes[child_i].kind != NODE_KIND::INT && nodes[child_i].kind != NODE_KIND::REAL ) break;
                            if ( vals != "" ) vals += ", ";
                            vals += (nodes[child_i].kind == NODE_KIND::INT) ? std::to_string(nodes[child_i].u.i) : std::to_string(nodes[child_i].u.r);
                        }
                        if ( node_is_hier( node ) ) {
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

                if ( node_is_hier( node ) ) {
                    if ( child_i == NULL_I ) child_i = node.u.child_first_i;
                    for( ; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
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

bool Layout::gds3d_write_layer_info( std::string file )
{
    std::ofstream out( file, std::ofstream::out );

    real height = 0.0;
    for( uint i = 0; i < hdr->layer_cnt; i++ )
    {
        const Layer& layer = layers[i];

        out << "LayerStart: " << &strings[layer.name_i] << "\n";
        out << "Layer:      " << i << "\n";
        if ( layer.gdsii_datatype != NULL_I ) out << "Datatype:   " << layer.gdsii_datatype << "\n";
        out << "Height:     " << height << "\n";
        real thickness = layer.thickness / gdsii_units_user;
        out << "Thickness:  " << thickness << "\n";
        out << "Red:        " << (real((layer.material_rgba >> 24) & 0xff) / 255.0) << "\n";
        out << "Green:      " << (real((layer.material_rgba >> 16) & 0xff) / 255.0) << "\n";
        out << "Blue:       " << (real((layer.material_rgba >>  8) & 0xff) / 255.0) << "\n";
        if ( false ) out << "Filter:     " << (1.0 - (real((layer.material_rgba >>  0) & 0xff) / 255.0)) << "\n";
        if ( i <= 9 ) out << "Shortkey:   " << i << "\n";
        out << "Show:       " << 1 << "\n";
        out << "LayerEnd\n\n";

        if ( i != (hdr->layer_cnt-1) && !layers[i+1].same_zoffset_as_prev ) height += thickness;
    }

    out.close();
    return true;
}

bool Layout::vbs_write_layer_info( std::string file )
{
    return false;
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
    lassert( fd >= 0, "could not open file " + file_path + " - open() error: " + strerror( errno ) );

    struct stat file_stat;
    int status = fstat( fd, &file_stat );
    if ( status < 0 ) {
        close( fd );
        lassert( 0, "could not stat file " + std::string(fname) + " - stat() error: " + strerror( errno ) );
    }
    size_t size = file_stat.st_size;

    // this large read should behave like an mmap() inside the o/s kernel and be as fast
    start = aligned_alloc<uint8_t>( size );
    if ( start == nullptr ) {
        close( fd );
        lassert( 0, "could not read file " + std::string(fname) + " - malloc() error: " + strerror( errno ) );
    }
    end = start + size;

    uint8_t * addr = start;
    while( size != 0 ) 
    {
        size_t _this_size = 1024*1024*1024;
        if ( size < _this_size ) _this_size = size;
        if ( ::read( fd, addr, _this_size ) <= 0 ) {
            close( fd );
            lassert( 0, "could not read() file " + std::string(fname) + " - read error: " + std::string( strerror( errno ) ) );
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
    lassert( xxx != xxx_end, "premature end of file" );
    lassert( *xxx == ch, "expected character '" + std::string(1, ch) + "' got '" + std::string( 1, *xxx ) + "' " + surrounding_lines( xxx, xxx_end ) );
    xxx++;
    return true;
}

uint Layout::str_get( std::string s )
{
    auto it = str_to_str_i.find( s );
    if ( it != str_to_str_i.end() ) return it->second;
        
    uint s_len = s.length();
    perhaps_realloc( strings, hdr->char_cnt, max->char_cnt, s_len+1 );
    uint s_i = hdr->char_cnt;
    char * to_s = &strings[s_i];
    hdr->char_cnt += s_len + 1;
    memcpy( to_s, s.c_str(), s_len+1 );
    str_to_str_i[s] = s_i;
    ldout << "str_i[" + s + "]=" + std::to_string(s_i) << "\n";
    return s_i;
}

uint Layout::str_find( std::string s ) const
{
    auto it = str_to_str_i.find( s );
    if ( it != str_to_str_i.end() ) return it->second;
    return NULL_I;
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
        lassert( xxx != xxx_end, "no terminating \" for std::string" );
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
            lassert( !has_exp, "real has more than one 'e' exponent" );
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

    lassert( s.length() != 0, "unable to parse real in " + ext_name + " file " + surrounding_lines( xxx, xxx_end ) );

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
            lassert( !is_neg, "-+ not allowed for an int" );
            xxx++;
            continue;
        }    
        if ( ch == '-' ) {
            lassert( !is_neg, "too many minus signs" );
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
    lassert( vld, "unable to parse int in " + ext_name + " file " + surrounding_lines( xxx, xxx_end ) );
    return true;
}

inline bool Layout::parse_uint( uint& u, uint8_t *& xxx, uint8_t *& xxx_end )
{
    _int i;
    if ( !parse_int( i, xxx, xxx_end ) ) return false;
    lassert( i >= 0, "parse_uint encountered negative integer" );
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
