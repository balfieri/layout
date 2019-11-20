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
#include <cstdio>

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
    static const uint NULL_I = -1;                 // null index into an array

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

    // write out material info used by viewer programs (file extension determines format)
    bool write_material_info( std::string file_path );

    // start of a new library; Layout must be empty
    uint start_library( std::string libname, real units_user=0.001, real units_meters=1e-9 );

    class real4
    {
    public:
        real c[4];
        
        real4( void )                               { c[0] = 0;  c[1] = 0;  c[2] = 0;  c[3] = 0;  }
        real4( real c0, real c1, real c2, real c3 ) { c[0] = c0; c[1] = c1; c[2] = c2; c[3] = c3; }

        real   dot( const real4 &v2 ) const;
        real   length( void ) const;
        real   length_sqr( void ) const ;
        real4& normalize( void );
        real4  normalized( void ) const;
        real4  operator + ( const real4& v ) const;
        real4  operator - ( const real4& v ) const;
        real4  operator * ( const real4& v ) const;
        real4  operator * ( real s ) const;
        real4  operator / ( const real4& v ) const;
        real4  operator / ( real s ) const;
        real4& operator += ( const real4 &v2 );
        real4& operator -= ( const real4 &v2 );
        real4& operator *= ( const real4 &v2 );
        real4& operator *= ( const real s );
        real4& operator /= ( const real4 &v2 );
        real4& operator /= ( const real s );
        std::string str( void ) const;
    };

    class real3
    {
    public:
        real c[3];
        
        real3( void )                      { c[0] = 0;  c[1] = 0;  c[2] = 0;  }
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
        std::string str( void ) const;
    };

    class real2
    {
    public:
        real c[2];

        real2( void )             { c[0] = 0;  c[1] = 0;  }
        real2( real c0, real c1 ) { c[0] = c0; c[1] = c1; }

        real   dot( const real2 &v2 ) const;
        real   length( void ) const;
        real   length_sqr( void ) const;
        real2& normalize( void );
        real2  normalized( void ) const;
        bool   colinear_is_same_dir( const real2& v2 ) const;
        bool   is_on_segment( const real2& p2, const real2& p3, bool include_endpoints, real epsilon=0.000001 ) const;  // returns true if this point is on line segment (p2, p3)
        bool   is_left_of_segment( const real2& p2, const real2& p3, real epsilon=0.000001 ) const;  // returns true if this point is to the left  of segment (p2, p3)
        bool   is_right_of_segment( const real2& p2, const real2& p3, real epsilon=0.000001 ) const; // returns true if this point is to the right of segment (p2, p3)
        int    orientation( const real2& p2, const real2& p3 ) const;         // returns: 0=colinear, 1=clockwise, 2=counterclockwise
        bool   lines_intersection( const real2& p2, const real2& p3, const real2& p4, real2& ip, real epsilon=0.000001 ) const;
        bool   segments_intersect( const real2& p2, const real2& p3, const real2& p4, bool include_endpoints ) const;
        bool   segments_intersection( const real2& p2, const real2& p3, const real2& p4, real2& ip, 
                                      bool include_p1_p2, bool include_p3_p4, bool include_colinear, bool& colinear, bool& reverse ) const;
        void   pad_segment( const real2& p2, real pad, real2& p1_new, real2& p2_new ) const;       // (p1_new, p2_new) are new endpoints
        void   perpendicular_segment( const real2& p2, real length, real2& p3, real2& p4 ) const;  // (p2, p3) will pass through p1
        void   parallel_segment( const real2& p2, real dist, real2& p1_new, real2& p2_new ) const; // (p1_new, p2_new) are new endpoints
        real2  operator + ( const real2& v ) const;
        real2  operator - ( const real2& v ) const;
        real2  operator * ( const real2& v ) const;
        real2  operator * ( real s ) const;
        real2  operator / ( const real2& v ) const;
        real2  operator / ( real s ) const;
        real2& operator += ( const real2 &v2 );
        real2& operator -= ( const real2 &v2 );
        real2& operator *= ( const real2 &v2 );
        real2& operator *= ( const real s );
        real2& operator /= ( const real2 &v2 );
        real2& operator /= ( const real s );
        bool   operator == ( const real2 &v2 ) const; 
        bool   operator != ( const real2 &v2 ) const; 
        bool   nearly_equal( const real2 &v2, real epsilon=0.000001 ) const;
        std::string str( void ) const;
    };

    // Axis-Aligned Bounding Rectangle (2D)
    //
    struct AABR
    {
        real2   min;
        real2   max;

        AABR( void ) {}
        AABR( const real2& p );                                 // init with one point

        void pad( real p );
        void expand( const AABR& other );
        void expand( const real2& p );
        void intersect( const AABR& other );
        bool encloses( const AABR& other ) const;
        bool encloses( const real2& p ) const;
        bool intersects( const AABR& other ) const;
        AABR quadrant( uint i, uint j ) const;
        AABR enclosing_square( real scale_factor=1.0 ) const;
    };

    class Matrix                                        // used for instancing to transform
    {
    public:
        real            m[4][4];

        Matrix( void )          { identity(); }

        void   identity(  void );                       // make this the identity matrix
        void   translate( const real3& translation );   // translate this matrix by a real3
        void   scale(     const real3& scaling );       // scale this matrix by a real3
        void   rotate_xz( double radians );             // rotate by radians in xz plane (yaw)
        void   rotate_yz( double radians );             // rotate by radians in yz plane (pitch)
        void   rotate_xy( double radians );             // rotate by radians in xy plane (roll)

        Matrix operator + ( const Matrix& m ) const;    // add two matrices
        Matrix operator - ( const Matrix& m ) const;    // subtract two matrices
        bool   operator == ( const Matrix& m ) const;   // return true if matrices are equal
        bool   operator != ( const Matrix& m ) const;   // return true if matrices are unequal

        void   multiply(    double s );                 // multiply this matrix by scalar

        real4  row( uint r ) const;                     // returns row r as a vector
        real4  column( uint c ) const;                  // returns column c as a vector
        void   transform( const real4& v, real4& r ) const; // r = *this * v
        void   transform( const real3& v, real3& r, bool div_by_w=false ) const; // r = *this * v (and optional divide by w)
        void   transform( const Matrix& M2, Matrix& M3 ) const; // M3 = *this * M2
        void   transpose( Matrix& mt ) const;           // return the transpose of this matrix 
    };

    enum class CONFLICT_POLICY
    {
        MERGE_NONE_ALLOW_NONE,                          // abort if any conflict is encountered
        MERGE_EXACT_ONLY,                               // merge exact overlaps only; abort partial overlaps
        MERGE_ALL,                                      // merge any conflict, including partial overlaps
    };

    // instancing of other layouts
    uint inst_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, real x, real y, std::string name );
    uint inst_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, real x, real y, 
                      uint dst_layer_first, uint dst_layer_last, std::string name );
    uint inst_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, const Matrix& M,
                      uint dst_layer_first, uint dst_layer_last, std::string name );
    void finalize_top_struct( uint parent_i, uint last_i, std::string top_name );              // use to create top-level struct of all insts
    uint flatten_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string struct_name, 
                         CONFLICT_POLICY conflict_policy, const Matrix& M = Matrix() );  
    AABR bounding_rect( uint layer_first=0, uint layer_last=0xffffffff ) const;
     
    // fill of dielectrics or arbitrary material
    uint fill_dielectrics( uint parent_i, const AABR& brect, uint layer_first=0, uint layer_last=0xffffffff );

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
        uint        leaf_node_cnt;          // in leaf_nodes array
        uint        bah_node_cnt;           // in bah_nodes array
        uint        root_i;                 // index of root node in nodes array
        uint *      layer_bah_root_i;       // index of per-layer root node in bah_nodes array
        AABR *      layer_bah_brect;        // overall bounding rectangle of per-layer BAH
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
        uint        gdsii_num;              	// layer number of main material in (input) GDSII file
        uint        dielectric_gdsii_num;   	// layer number of dielectric material in (output) GDSII file
        uint        gdsii_datatype;         	// which datatype to use - NULL_I means all
        uint        dielectric_gdsii_datatype; 	// which datatype to use - NULL_I means all
        bool        same_zoffset_as_prev;   	// starts at same zoffset as previous layer in stackkup?
        real        thickness;              	// thickness in um
        uint        material_i;             	// index of main material in materials[]
        uint        dielectric_material_i;  	// index of dielectric material in materials[]
        uint        material_rgba;              // material color in RGBA8 format
        uint        bah_root_i;                 // index in bah_nodes[] of this layer's BAH root node (if !same_zoffset_as_prev)
        AABR        bah_brect;                  // bounding rectangle of this layer's BAH             (if !same_zoffset_as_prev)
    };

    static uint color( real r, real g, real b, real a=1.0 );
    static uint color( std::string name, real a=1.0 );
    static std::string color_name( uint color );                // returns closest name for color

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
            
    static std::string  str( NODE_KIND kind );
    NODE_KIND           hier_end_kind( NODE_KIND kind ) const;      // returns corresponding end kind for hier kind
    
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

    // NODE CREATION
    uint        node_alloc( NODE_KIND kind );                   // allocate a node of the given kind
    uint        node_alloc_int( _int i );                       // allocate INT node
    uint        node_alloc_real( real r );                      // allocate REAL node
    uint        node_alloc_str( std::string s );                // allocate STR node
    uint        node_alloc_boundary( uint parent_i, uint last_i, uint gdsii_num, uint datatype, const real2 * vtx, uint vtx_cnt );
    void        node_alloc_timestamp( Node& node );             // adds timestamp fields to node

    // NODE QUERIES
    bool        node_is_header_footer( const Node& node ) const;// return true if node is a HEADER, BGNLIB, LIBNAME, UNITS, or ENDLIB
    bool        node_is_gdsii( const Node& node ) const;        // return true if node is a GDSII node
    bool        node_is_scalar( const Node& node ) const;       // return true if node is a scalar 
    bool        node_is_string( const Node& node ) const;       // return true if node is a (scalar) string
    bool        node_is_parent( const Node& node ) const;       // return true if node is not a scalar (i.e., could have children)
    bool        node_is_element( const Node& node ) const;      // return true if node is a GDSII element
    bool        node_is_name( const Node& node ) const;         // return true if node is a name node
    bool        node_is_hier( const Node& node ) const;         // return true if node is a hierarchy
    bool        node_is_ref( const Node& node ) const;          // return true if node is an AREF or SREF
    uint        node_last_i( const Node& node ) const;          // find node index of last child, else NULL_I if none
    uint        node_last_scalar_i( const Node& node ) const;   // find node index of last scalar child, else NULL_I if none
    uint        node_name_i( const Node& node ) const;          // find name for node but return strings[] index
    std::string node_name( const Node& node ) const;            // find name for node
    uint        node_layer( const Node& node ) const;           // find LAYER value for node (an element)
    uint        node_width( const Node& node ) const;           // find WIDTH value for node 
    uint        node_bah_layer( const Node& node ) const;       // find LAYER value for flattened node and get layer to use for BAH
    uint        node_xy_i( const Node& node ) const;            // find index of XY node within node
    void        node_xy_replace_polygon( Node& node, const real2 * vtx, uint vtx_cnt ); // replace X,Y children from new set of vertices
    uint        node_pathtype( const Node& node ) const;        // find PATHTYPE for node and return number (default: 0)
    uint        node_datatype( const Node& node ) const;        // find DATATYPE for node and return number (default: 0)

    // NODE COPIES
    enum class COPY_KIND
    {
        ONE,                                // copy only the one source node, no children
        SCALAR_CHILDREN,                    // copy scalar children only (INT, REAL, etc.)
        DEEP,                               // copy all children and descendents
        FLATTEN,                            // copy all children but flatten all REFs
    };

    uint        node_copy( uint parent_i, uint last_i, const Layout * src_layout, uint src_i, COPY_KIND copy_kind, 
                           CONFLICT_POLICY conflict_policy=CONFLICT_POLICY::MERGE_NONE_ALLOW_NONE, const Matrix& M = Matrix(), bool in_flatten=false );

    // RAW POLYGONS
    real2 *     polygon_alloc( uint vtx_cnt ) const;
    real2 *     polygon_alloc( const Node& xy_node, uint& vtx_cnt ) const;
    real2 *     polygon_alloc( const AABR& brect, uint& vtx_cnt ) const;
    real2 *     polygon_copy( const real2 * other, uint other_cnt ) const;
    real2 *     polygon_reversed( const real2 * other, uint other_cnt ) const;
    void        polygon_dealloc( real2 * vtx_array ) const;
    real2 *     polygon_merge_or_intersect( bool do_merge, const real2 * vtx1, uint vtx1_cnt, const real2 * vtx2, uint vtx2_cnt, uint& vtx_cnt ) const;
    AABR        polygon_brect( const real2 * vtx, uint vtx_cnt ) const;
    std::string polygon_str( const real2 * vtx, uint vtx_cnt, std::string color, real xy_scale=4.0, real x_off=0.0, real y_off=0.0 ) const;
    bool        polygon_eq( const real2 * vtx1, uint vtx1_cnt, const real2 * vtx2, uint vtx2_cnt ) const;
    bool        polygon_encloses( const real2 * vtx, uint vtx_cnt, const real2& p ) const;
    bool        polygon_includes( const real2 * vtx, uint vtx_cnt, const real2& v ) const;  // v is already in the vtx list?
    bool        polygon_is_ccw( const real2 * vtx, uint vtx_cnt ) const;        // are vertices in counterclockwise winding order?
    real2 *     polygon_ccw( const real2 * vtx, uint vtx_cnt ) const;           // reverses vertices if they are in clockwise winding order, else just a copy

    struct TopInstInfo
    {
        uint            struct_i;
        Matrix          M;
    };

    // Bounding Area Hierarchy
    //
    struct Leaf_Node
    {
        uint            node_i;             // in nodes[] array
        AABR            brect;              // bounding rectangle
    };

    struct BAH_Node                         // quadtree node, bounding rect is implicit based on position in tree
    {
        uint            child_i[2][2];      // bah_nodes[] or leaf_nodes[] index for quad space
        bool            child_is_leaf[2][2];// true=Leaf_Node, false=is BAH_Node
    };

    static GDSII_DATATYPE kind_to_datatype( NODE_KIND kind );
    static std::string    str( GDSII_DATATYPE datatype );

    // global scalars
    static const uint   VERSION = 0xB0BA1f01; // current version 

    // structs
    std::string         file_path;          // pathname of file passed to Layout()
    uint8_t *           mapped_region;      // != nullptr means the whole file was sucked in by read_uncompressed()
    Header *            hdr;
    Header *            max;                // holds max lengths of currently allocated arrays 

    // arrays                               // these are relocatable and will all get written out to .layout file
    char *              strings;
    Material *          materials;
    Layer *             layers;
    Node *              nodes;
    TopInstInfo *       top_insts;
    Leaf_Node *         leaf_nodes;
    BAH_Node *          bah_nodes;


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
    void init( bool alloc_arrays );

    // STRUCTS
    std::map<std::string, uint>                 str_to_str_i;           // maps std::string to unique location in strings[]
    std::map< uint, uint >                      name_i_to_struct_i;
    std::string all_struct_names( std::string delim = "\n    " ) const;

    // INSTANCING AND NODE COPYING
    using has_layer_cache_t = std::map< uint, std::map<uint, bool>* >;
    bool node_has_layer( uint ni, uint layer_num, has_layer_cache_t * cache, std::string indent_str ) const;
    uint inst_layout_node( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, uint src_i, uint src_layer_num, 
                           uint dst_layer_num, has_layer_cache_t * cache, std::string name, std::string indent_str="" );
    uint node_flatten_ref( uint parent_i, uint last_i, const Layout * src_layout, uint src_i, CONFLICT_POLICY conflict_policy, const Matrix& M );
    uint node_convert_path_to_boundary( uint parent_i, uint last_i, const Layout * src_layout, uint src_i,
                           CONFLICT_POLICY conflict_policy, const Matrix& M );
    void node_transform_xy( uint parent_i, uint xy_i, COPY_KIND copy_kind, CONFLICT_POLICY conflict_policy, const Matrix& M );

    // BAH 
    uint bah_node_alloc( void );
    void bah_add( uint leaf_i, uint layer, CONFLICT_POLICY conflict_policy );
    void bah_insert( uint bah_i, const AABR& brect, uint leaf_i, const AABR& leaf_brect, CONFLICT_POLICY conflict_policy, std::string indent_str="" );
    bool bah_leaf_nodes_intersect( const AABR& quadrant_brect, uint li1, uint li2, bool& is_exact ); 

    // FILL
    uint fill_dielectric_bah( uint parent_i, uint last_i, uint layer_i, uint bah_i, const AABR& rect );
    uint fill_dielectric_leaf( uint parent_i, uint last_i, uint layer_i, uint bah_i, uint x, uint y, const AABR& rect );
    uint fill_dielectric_rect( uint parent_i, uint last_i, uint layer_i, const AABR& rect );
    uint fill_dielectric_polygon( uint parent_i, uint last_i, uint layer_i, uint datatype, const real2 * vtx, uint vtx_cnt );

    // LAYOUT I/O
    bool layout_read( std::string file_path );          // .layout
    bool layout_write( std::string file_path );         

    // GDSII I/O
    uint        gdsii_rec_cnt;
    NODE_KIND   gdsii_last_kind;
    int         gdsii_fd;
    uint8_t *   gdsii_buff;
    uint        gdsii_buff_byte_cnt;
    real        gdsii_units_user;
    real        gdsii_units_meters;

    bool gdsii_read( std::string file_path, bool count_only ); // .gds
    bool gdsii_read_record( uint& node_i, uint curr_struct_i, bool count_only );
    bool gdsii_write( std::string file );
    void gdsii_write_record( uint node_i, std::string indent_str="" );
    void gdsii_write_number( uint8_t * bytes, uint& byte_cnt, uint ni, GDSII_DATATYPE datatype, std::string indent_str );
    void gdsii_write_bytes( const uint8_t * bytes, uint byte_cnt );
    void gdsii_flush( bool for_end_of_file=false );
    
    // AEDT I/O
    uint aedt_begin_str_i;              // these are to make it easier to compare
    uint aedt_end_str_i;
    uint true_str_i;
    uint false_str_i;

    bool aedt_read( std::string file );                 // .aedt
    bool aedt_read_expr( uint& node_i );
    bool aedt_write( std::string file, bool for_raw );
    void aedt_write_expr( std::ofstream& out, uint node_i, std::string indent_str, bool for_raw );

    bool gds3d_write_layer_info( std::string file );
    bool hfss_write_layer_info( std::string file );

    bool hfss_write_material_info( std::string file );

    // PARSING
    std::string ext_name;
    uint8_t * nnn_start;
    uint8_t * nnn_end;
    uint8_t * nnn;
    uint      line_num;

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

    // SYSTEM 
    bool cmd( std::string c );

    // allocates an array of T on a page boundary
    template<typename T>
    T * aligned_alloc( size_t cnt );

    // reallocate array if we are about to exceed its current size
    template<typename T>
    void perhaps_realloc( T *& array, const uint  & hdr_cnt, uint  & max_cnt, uint   add_cnt );

    // printf-style formatting for ostream output
    static std::string putf( const char * fmt, ... );
};

#ifdef LAYOUT_DEBUG
#define ldout if ( true )  std::cout 
#else
#define ldout if ( false ) std::cout
#endif

// these are done as macros to avoid evaluating msg (it makes a big difference)
#define lassert( bool, msg ) if ( !(bool) ) { std::cout << "ERROR: " << std::string(msg) << "\n"; exit( 1 ); }

static inline std::string str( Layout::real r, int n = 30 )
{
    std::ostringstream out;
    out.precision( n );
    out << std::fixed << r;
    return out.str();
}

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
        max->leaf_node_cnt = 128;
        max->bah_node_cnt = 128;

        //------------------------------------------------------------
        // Allocate initial arrays
        //------------------------------------------------------------
        strings    = aligned_alloc<char>( max->char_cnt );
        materials  = aligned_alloc<Material>( max->material_cnt );
        layers     = aligned_alloc<Layer>( max->layer_cnt );
        nodes      = aligned_alloc<Node>( max->node_cnt );
        top_insts  = aligned_alloc<TopInstInfo>( max->top_inst_cnt );
        leaf_nodes = aligned_alloc<Leaf_Node>( max->leaf_node_cnt );
        bah_nodes  = aligned_alloc<BAH_Node>( max->bah_node_cnt );

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
        delete leaf_nodes;
        delete bah_nodes;
        strings = nullptr;
        materials = nullptr;
        layers = nullptr;
        nodes = nullptr;
        top_insts = nullptr;
        leaf_nodes = nullptr;
        bah_nodes = nullptr;
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
            return aedt_write( top_file, false );
        } else if ( ext_name == std::string( ".raw" ) ) {
            return aedt_write( top_file, true );
        } else if ( ext_name == std::string( ".gds" ) ) {
            return gdsii_write( top_file );
        } else {
            lassert( false, "unknown file ext_name: " + ext_name );
            return false;
        }
    }
}

bool Layout::write_layer_info( std::string file_name )
{
    //------------------------------------------------------------
    // Write depends on file ext_name
    //------------------------------------------------------------
    std::string dir_name;
    std::string base_name;
    dissect_path( file_name, dir_name, base_name, ext_name );
    if ( ext_name == std::string( ".gds3d" ) ) {
        return gds3d_write_layer_info( file_name );
    } else if ( ext_name == std::string( ".tech" ) || ext_name == std::string( ".hfss" ) ) {
        return hfss_write_layer_info( file_name );
    } else {
        lassert( false, "write_layer_info: unknown file ext_name: " + ext_name );
        return false;
    }
}

bool Layout::write_material_info( std::string file_name )
{
    //------------------------------------------------------------
    // Write depends on file ext_name
    //------------------------------------------------------------
    std::string dir_name;
    std::string base_name;
    dissect_path( file_name, dir_name, base_name, ext_name );
    if ( ext_name == std::string( ".amat" ) || ext_name == std::string( ".hfss" ) ) {
        return hfss_write_material_info( file_name );
    } else {
        lassert( false, "write_material_info: unknown file ext_name: " + ext_name );
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
    { "dark_red", 		0x8B0000 },
    { "brown", 			0xA52A2A },
    { "firebrick", 		0xB22222 },
    { "crimson", 		0xDC143C },
    { "red", 			0xFF0000 },
    { "r", 			0xFF0000 },
    { "tomato", 		0xFF6347 },
    { "coral", 			0xFF7F50 },
    { "indian_red", 		0xCD5C5C },
    { "light_coral", 		0xF08080 },
    { "dark_salmon", 		0xE9967A },
    { "salmon", 		0xFA8072 },
    { "light_salmon", 		0xFFA07A },
    { "orange_red", 		0xFF4500 },
    { "dark_orange", 		0xFF8C00 },
    { "orange", 		0xFFA500 },
    { "o", 		        0xFFA500 },
    { "gold", 			0xFFD700 },
    { "dark_golden_rod", 	0xB8860B },
    { "golden_rod", 		0xDAA520 },
    { "pale_golden_rod", 	0xEEE8AA },
    { "dark_khaki", 		0xBDB76B },
    { "khaki", 			0xF0E68C },
    { "olive", 			0x808000 },
    { "yellow", 		0xFFFF00 },
    { "y", 		        0xFFFF00 },
    { "yellow_green", 		0x9ACD32 },
    { "dark_olive_green",	0x556B2F },
    { "olive_drab", 		0x6B8E23 },
    { "lawn_green", 		0x7CFC00 },
    { "chartreuse", 		0x7FFF00 },
    { "green_yellow", 		0xADFF2F },
    { "dark_green", 		0x006400 },
    { "green", 			0x00FF00 },
    { "g", 			0x00FF00 },
    { "forest_green", 		0x228B22 },
    { "lime", 			0x00FF00 },
    { "lime_green", 		0x32CD32 },
    { "light_green", 		0x90EE90 },
    { "pale_green", 		0x98FB98 },
    { "dark_sea_green", 	0x8FBC8F },
    { "medium_spring_green", 	0x00FA9A },
    { "spring_green", 		0x00FF7F },
    { "sea_green", 		0x2E8B57 },
    { "medium_aqua_marine", 	0x66CDAA },
    { "medium_sea_green", 	0x3CB371 },
    { "light_sea_green", 	0x20B2AA },
    { "dark_slate_gray", 	0x2F4F4F },
    { "teal", 			0x008080 },
    { "dark_cyan", 		0x008B8B },
    { "cyan", 			0x00FFFF },
    { "c", 			0x00FFFF },
    { "aqua", 			0x00FFFF },
    { "light_cyan", 		0xE0FFFF },
    { "dark_turquoise", 	0x00CED1 },
    { "turquoise", 		0x40E0D0 },
    { "medium_turquoise", 	0x48D1CC },
    { "pale_turquoise", 	0xAFEEEE },
    { "aqua_marine", 		0x7FFFD4 },
    { "powder_blue", 		0xB0E0E6 },
    { "cadet_blue", 		0x5F9EA0 },
    { "steel_blue", 		0x4682B4 },
    { "corn_flowerNblue", 	0x6495ED },
    { "deep_sky_blue", 		0x00BFFF },
    { "dodger_blue", 		0x1E90FF },
    { "light_blue", 		0xADD8E6 },
    { "sky_blue", 		0x87CEEB },
    { "light_sky_blue", 	0x87CEFA },
    { "midnight_blue", 		0x191970 },
    { "navy", 			0x000080 },
    { "dark_blue", 		0x00008B },
    { "medium_blue", 		0x0000CD },
    { "blue", 			0x0000FF },
    { "b", 			0x0000FF },
    { "royal_blue", 		0x4169E1 },
    { "blue_violet", 		0x8A2BE2 },
    { "indigo", 		0x4B0082 },
    { "dark_slate_blue", 	0x483D8B },
    { "slate_blue", 		0x6A5ACD },
    { "medium_slate_blue", 	0x7B68EE },
    { "medium_purple", 		0x9370DB },
    { "dark_magenta", 		0x8B008B },
    { "dark_violet", 		0x9400D3 },
    { "dark_orchid", 		0x9932CC },
    { "medium_orchid", 		0xBA55D3 },
    { "purple", 		0x800080 },
    { "p", 		        0x800080 },
    { "thistle", 		0xD8BFD8 },
    { "plum", 			0xDDA0DD },
    { "violet", 		0xEE82EE },
    { "magenta",                0xFF00FF },
    { "m",                      0xFF00FF },
    { "fuchsia", 	        0xFF00FF },
    { "orchid", 		0xDA70D6 },
    { "medium_violet_red", 	0xC71585 },
    { "pale_violetNred", 	0xDB7093 },
    { "deep_pink", 		0xFF1493 },
    { "hot_pink", 		0xFF69B4 },
    { "light_pink", 		0xFFB6C1 },
    { "pink", 			0xFFC0CB },
    { "antique_white", 		0xFAEBD7 },
    { "beige", 			0xF5F5DC },
    { "bisque", 		0xFFE4C4 },
    { "blanched_almond", 	0xFFEBCD },
    { "wheat", 			0xF5DEB3 },
    { "corn_silk", 		0xFFF8DC },
    { "lemon_chiffon", 		0xFFFACD },
    { "light_golden_rod_yellow",0xFAFAD2 },
    { "light_yellow", 		0xFFFFE0 },
    { "saddle_brown", 		0x8B4513 },
    { "sienna", 		0xA0522D },
    { "chocolate", 		0xD2691E },
    { "peru", 			0xCD853F },
    { "sandy_brown", 		0xF4A460 },
    { "burly_wood", 		0xDEB887 },
    { "tan", 			0xD2B48C },
    { "rosy_brown", 		0xBC8F8F },
    { "moccasin", 		0xFFE4B5 },
    { "navajo_white", 		0xFFDEAD },
    { "peach_puff", 		0xFFDAB9 },
    { "misty_rose", 		0xFFE4E1 },
    { "lavender_blush",		0xFFF0F5 },
    { "linen", 			0xFAF0E6 },
    { "old_lace", 		0xFDF5E6 },
    { "papaya_whip", 		0xFFEFD5 },
    { "seashell", 		0xFFF5EE },
    { "mint_cream", 		0xF5FFFA },
    { "slate_gray", 		0x708090 },
    { "light_slate_gray", 	0x778899 },
    { "light_steel_blue", 	0xB0C4DE },
    { "lavender", 		0xE6E6FA },
    { "floral_white", 		0xFFFAF0 },
    { "alice_blue", 		0xF0F8FF },
    { "ghost_white", 		0xF8F8FF },
    { "honeydew", 		0xF0FFF0 },
    { "ivory", 			0xFFFFF0 },
    { "azure", 			0xF0FFFF },
    { "snow", 			0xFFFAFA },
    { "black", 			0x000000 },
    { "dim_gray",               0x696969 },
    { "gray", 		        0x808080 },
    { "dark_gray", 	        0xA9A9A9 },
    { "silver", 		0xC0C0C0 },
    { "light_gray", 	        0xD3D3D3 },
    { "gainsboro", 		0xDCDCDC },
    { "white_smoke", 		0xF5F5F5 },
    { "w", 			0xFFFFFF },
    { "white", 			0xFFFFFF },
    { "copper",                 0xB87333 },
    { "nvidia_green",           (118 << 16) | (185 << 8) | (0 << 0) },
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

std::string Layout::color_name( uint color )
{
    // return closest color's name
    _int rdiff_best = 1000000;
    _int gdiff_best = 1000000;
    _int bdiff_best = 1000000;

    _int r = (color >> 24) & 0xff;
    _int g = (color >> 16) & 0xff;
    _int b = (color >>  8) & 0xff;

    size_t best_i = 0;
    for( size_t i = 0; i < color_info_cnt; i++ )
    {
        uint rgb = color_info[i].rgb;
        _int rr = (rgb >> 16) & 0xff;
        _int gg = (rgb >>  8) & 0xff;
        _int bb = (rgb >>  0) & 0xff;

        _int rdiff = std::abs(r - rr);
        _int gdiff = std::abs(g - gg);
        _int bdiff = std::abs(b - bb);

        if ( (rdiff+gdiff+bdiff) < (rdiff_best+gdiff_best+bdiff_best) ) {
            best_i = i;        
            if ( (rdiff+gdiff+bdiff) == 0 ) break;
        }
    }

    return color_info[best_i].name;
}

void Layout::materials_init( void )
{
    // common materials
    //
    // permittivity   - amount of charge needed to generate one unit of electric flux in a given medium
    // relativity_permittivity (aka dialectric_constant) - permittivity/permittivity_of_vacuum
    // permeability   - measure of the ability of a material to support the formation of a magnetic field 
    //                  within itself, otherwise known as distributed inductance 
    // conductivity   - 1/resistivity 
    // thermal_conductivity - measure of ability to conduct heat
    // mass_density   - mass/volume
    // specific_heat  - Joule/Kelvin
    // youngs_modulus - uniaxial_stress/strain
    // poissons_ratio - transverse_strain/axial_strain
    // thermal_expansion_coefficient - 1/V * dV/dT (pressure held constant)
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
    mi = hdr->material_cnt++;
    materials[mi] = Material{ str_get( "port" ),// fake material representing a source/sink point in a simulation
                              1.0,              // relative_permittivity
                              1.0,              // permeability
                              0.0,              // conductivity
                              0.0,              // thermal_conductivity
                              0.0,              // mass_density
                              0.0,              // specific_heat
                              0.0,              // youngs_modulus
                              0.0,              // poissons_ratio
                              0.0 };            // thermal_expansion_coefficient
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

inline std::string Layout::putf( const char * fmt, ... ) 
{
    char buff[1024];
    va_list args;
    va_start( args, fmt );
    vsnprintf( buff, sizeof(buff), fmt , args );
    va_end( args );
    return buff;
}

inline std::ostream& operator << ( std::ostream& os, const Layout::real4& v ) 
{
    os << "[" << v.c[0] << "," << v.c[1] << "," << v.c[2] << "," << v.c[3] << "]";
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Layout::real3& v ) 
{
    os << "[" << v.c[0] << "," << v.c[1] << "," << v.c[2] << "]";
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Layout::real2& v ) 
{
    os << "[" << v.c[0] << "," << v.c[1] << "]";
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Layout::Matrix& m ) 
{
    os << "[ [" << m.m[0][0] << "," << m.m[0][1] << "," << m.m[0][2] << "," << m.m[0][3] << "],\n" << 
          "  [" << m.m[1][0] << "," << m.m[1][1] << "," << m.m[1][2] << "," << m.m[1][3] << "],\n" << 
          "  [" << m.m[2][0] << "," << m.m[2][1] << "," << m.m[2][2] << "," << m.m[2][3] << "],\n" << 
          "  [" << m.m[3][0] << "," << m.m[3][1] << "," << m.m[3][2] << "," << m.m[3][3] << "] ]\n";
    return os;
}

inline std::ostream& operator << ( std::ostream& os, const Layout::AABR r )
{
    os << r.min << " .. " << r.max;
    return os;
}

inline Layout::real Layout::real4::dot( const Layout::real4 &v2 ) const
{
    return c[0] * v2.c[0] + c[1] * v2.c[1] + c[2] * v2.c[2] + c[3] * v2.c[3];
}

inline Layout::real Layout::real4::length( void ) const
{ 
    return std::sqrt( c[0]*c[0] + c[1]*c[1] + c[2]*c[2] + c[3]*c[3] ); 
}

inline Layout::real Layout::real4::length_sqr( void ) const 
{ 
    return c[0]*c[0] + c[1]*c[1] + c[2]*c[2] + c[3]*c[3];
}

inline Layout::real4& Layout::real4::normalize( void )
{
    *this /= length();
    return *this;
}

inline Layout::real4 Layout::real4::normalized( void ) const
{
    return *this / length();
}

inline Layout::real4 Layout::real4::operator + ( const Layout::real4& v2 ) const
{
    real4 r;
    r.c[0] = c[0] + v2.c[0];
    r.c[1] = c[1] + v2.c[1];
    r.c[2] = c[2] + v2.c[2];
    r.c[3] = c[3] + v2.c[3];
    return r;
}

inline Layout::real4 Layout::real4::operator - ( const Layout::real4& v2 ) const
{
    real4 r;
    r.c[0] = c[0] - v2.c[0];
    r.c[1] = c[1] - v2.c[1];
    r.c[2] = c[2] - v2.c[2];
    r.c[3] = c[3] - v2.c[3];
    return r;
}

inline Layout::real4 Layout::real4::operator * ( const Layout::real4& v2 ) const
{
    real4 r;
    r.c[0] = c[0] * v2.c[0];
    r.c[1] = c[1] * v2.c[1];
    r.c[2] = c[2] * v2.c[2];
    r.c[3] = c[3] * v2.c[3];
    return r;
}

inline Layout::real4 operator * ( Layout::real s, const Layout::real4& v ) 
{
    return Layout::real4( s*v.c[0], s*v.c[1], s*v.c[2], s*v.c[3] );
}

inline Layout::real4 Layout::real4::operator * ( Layout::real s ) const
{
    real4 r;
    r.c[0] = c[0] * s;
    r.c[1] = c[1] * s;
    r.c[2] = c[2] * s;
    r.c[3] = c[3] * s;
    return r;
}

inline Layout::real4 Layout::real4::operator / ( const Layout::real4& v2 ) const
{
    real4 r;
    r.c[0] = c[0] / v2.c[0];
    r.c[1] = c[1] / v2.c[1];
    r.c[2] = c[2] / v2.c[2];
    r.c[3] = c[3] / v2.c[3];
    return r;
}

inline Layout::real4 Layout::real4::operator / ( Layout::real s ) const
{
    real4 r;
    r.c[0] = c[0] / s;
    r.c[1] = c[1] / s;
    r.c[2] = c[2] / s;
    r.c[3] = c[3] / s;
    return r;
}

inline Layout::real4& Layout::real4::operator += ( const Layout::real4 &v2 )
{
    c[0] += v2.c[0];
    c[1] += v2.c[1];
    c[2] += v2.c[2];
    c[3] += v2.c[3];
    return *this;
}

inline Layout::real4& Layout::real4::operator -= ( const Layout::real4 &v2 )
{
    c[0] -= v2.c[0];
    c[1] -= v2.c[1];
    c[2] -= v2.c[2];
    c[3] -= v2.c[3];
    return *this;
}

inline Layout::real4& Layout::real4::operator *= ( const Layout::real4 &v2 )
{
    c[0] *= v2.c[0];
    c[1] *= v2.c[1];
    c[2] *= v2.c[2];
    c[3] *= v2.c[3];
    return *this;
}

inline Layout::real4& Layout::real4::operator *= ( const Layout::real s )
{
    c[0] *= s;
    c[1] *= s;
    c[2] *= s;
    c[3] *= s;
    return *this;
}

inline Layout::real4& Layout::real4::operator /= ( const Layout::real4 &v2 )
{
    c[0] /= v2.c[0];
    c[1] /= v2.c[1];
    c[2] /= v2.c[2];
    c[3] /= v2.c[3];
    return *this;
}

inline Layout::real4& Layout::real4::operator /= ( const Layout::real s )
{
    c[0] /= s;
    c[1] /= s;
    c[2] /= s;
    c[3] /= s;
    return *this;
}

inline std::string Layout::real4::str( void ) const
{
    return "[" + ::str(c[0]) + "," + ::str(c[1]) + "," + ::str(c[2]) + "," + ::str(c[3]) + "]";
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

inline std::string Layout::real3::str( void ) const
{
    return "[" + ::str(c[0]) + "," + ::str(c[1]) + "," + ::str(c[2]) + "]";
}

inline Layout::real Layout::real2::dot( const Layout::real2 &v2 ) const
{
    return c[0] * v2.c[0] + c[1] * v2.c[1];
}

inline Layout::real Layout::real2::length( void ) const
{ 
    return std::sqrt( c[0]*c[0] + c[1]*c[1] );
}

inline Layout::real Layout::real2::length_sqr( void ) const 
{ 
    return c[0]*c[0] + c[1]*c[1];
}

inline Layout::real2& Layout::real2::normalize( void )
{
    *this /= length();
    return *this;
}

inline Layout::real2 Layout::real2::normalized( void ) const
{
    return *this / length();
}

inline bool Layout::real2::colinear_is_same_dir( const Layout::real2& v2 ) const
{
    return std::signbit( c[0] ) == std::signbit( v2.c[0] ) &&
           std::signbit( c[1] ) == std::signbit( v2.c[1] );
}

inline bool Layout::real2::is_on_segment( const real2& p2, const real2& p3, bool include_endpoints, real epsilon ) const
{
    const real2& p1 = *this;

    real2 _min( std::min( p2.c[0], p3.c[0] ), std::min( p2.c[1], p3.c[1] ) );
    real2 _max( std::max( p2.c[0], p3.c[0] ), std::max( p2.c[1], p3.c[1] ) );
    bool on_segment = (p1.c[0]+epsilon) >= _min.c[0] && (p1.c[0]-epsilon) <= _max.c[0] &&
                      (p1.c[1]+epsilon) >= _min.c[1] && (p1.c[1]-epsilon) <= _max.c[1];
    ldout << ", _min=" << _min << " _max=" << _max << " on_segment=" << on_segment;
    return on_segment && (include_endpoints || (!p1.nearly_equal( p2 ) && !p1.nearly_equal( p3 )));
}

inline bool Layout::real2::is_left_of_segment( const real2& p2, const real2& p3, real epsilon ) const
{
    return (p3.c[0] - p2.c[0])*(c[1] - p2.c[1]) > (p3.c[1] - p2.c[1])*(c[0] - p2.c[0]) &&
           !is_on_segment( p2, p3, false, epsilon );
}

inline bool Layout::real2::is_right_of_segment( const real2& p2, const real2& p3, real epsilon ) const
{
    return (p3.c[0] - p2.c[0])*(c[1] - p2.c[1]) < (p3.c[1] - p2.c[1])*(c[0] - p2.c[0]) &&
           !is_on_segment( p2, p3, false, epsilon );
}

inline int Layout::real2::orientation( const real2& p2, const real2& p3 ) const
{
    const real2& p1 = *this;

    // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
    // for details of below formula.
    //
    real val = (p2.c[1] - p1.c[1]) * (p3.c[0] - p2.c[0]) -
               (p2.c[0] - p1.c[0]) * (p3.c[1] - p2.c[1]);

    return (val == 0) ? 0 :     // colinear
           (val >  0) ? 1 : 2;  // clockwise or counterclockwise
}

inline bool Layout::real2::segments_intersect( const real2& p2, const real2& p3, const real2& p4, bool include_endpoints ) const
{
    const real2& p1 = *this;

    // Find the four orientations needed for general and special cases
    int o1 = p1.orientation( p2, p3 );
    int o2 = p1.orientation( p2, p4 );
    int o3 = p3.orientation( p4, p1 );
    int o4 = p3.orientation( p4, p2 );

    // General case
    if ( o1 != o2 && o3 != o4 ) return true;

    // Special Cases
    // p1, p2 and p3 are colinear and p3 lies on segment p1p2
    if ( o1 == 0 && p3.is_on_segment( p1, p2, include_endpoints ) ) return true;

    // p1, p2 and p4 are colinear and p4 lies on segment p1p2
    if ( o2 == 0 && p4.is_on_segment( p1, p2, include_endpoints ) ) return true;

    // p3, p4 and p1 are colinear and p1 lies on segment p3p4
    if ( o3 == 0 && p1.is_on_segment( p3, p4, include_endpoints ) ) return true; 

    // p3, p4 and p2 are colinear and p2 lies on segment p3p4
    if ( o4 == 0 && p2.is_on_segment( p3, p4, include_endpoints ) ) return true;

    return false; 
}

inline bool Layout::real2::lines_intersection( const real2& p2, const real2& p3, const real2& p4, real2& ip, real epsilon ) const
{
    const real2& p1 = *this;

    // Line (p1, p2) is represented as a1x + b1y = c1 
    real  a1 = p2.c[1] - p1.c[1];
    real  b1 = p1.c[0] - p2.c[0];
    real  c1 = a1*p1.c[0] + b1*p1.c[1];

    // Line (p3, p4) is represented as a2x + b2y = c2 
    real  a2 = p4.c[1] - p3.c[1];
    real  b2 = p3.c[0] - p4.c[0];
    real  c2 = a2*p3.c[0] + b2*p3.c[1];
  
    real determinant = a1*b2 - a2*b1; 
  
    if ( determinant >= -epsilon && determinant <= epsilon ) {
        ldout << " lines_intersection=NONE";
        return false;   // parallel
    } else {
        ip.c[0] = (b2*c1 - b1*c2) / determinant; 
        ip.c[1] = (a1*c2 - a2*c1) / determinant; 
        ldout << " lines_intersection ip=" << ip.str() << " determinant=" << ::str(determinant);
        return true;
    } 
}

inline bool Layout::real2::segments_intersection( const real2& p2, const real2& p3, const real2& p4, real2& ip, 
                                                  bool include_p1_p2, bool include_p3_p4, bool include_colinear, bool& colinear, bool& reverse ) const
{
    const real2& p1 = *this;

    ldout << " checking intersection of segment p1=[" << p1 << ", p2=" << p2 << "] against p3=[" << p3 << ", p4=" << p4 << 
                        "] include_p1_p2_endpoints=" << include_p1_p2 << " include_p3_p4_endpoints=" << include_p3_p4 << 
                        " include_colinear=" << include_colinear << ": ";

    colinear = false;
    reverse = false;
    if ( p1.lines_intersection( p2, p3, p4, ip ) ) {
        if ( ip.is_on_segment( p1, p2, include_p1_p2 ) ) {
            ldout << ", on_p1_p2=true";
            if ( ip.is_on_segment( p3, p4, include_p3_p4 ) ) {
                ldout << ", on_p3_p4=true => INTERSECTION\n";
                return true;
            } else {
                ldout << ", on_p3_p4=false";
            }
        } else {
            ldout << ", on_p1_p2=false";
        }
    } else if ( include_colinear ) {
        bool p3_on_p1_p2_proper = p3.is_on_segment( p1, p2, false );
        bool p4_on_p1_p2_proper = p4.is_on_segment( p1, p2, false );
        bool p3_on_p1_p2 = p3_on_p1_p2_proper || (include_p3_p4 && p3.nearly_equal( p2 ));
        bool p4_on_p1_p2 = p4_on_p1_p2_proper || (include_p3_p4 && p4.nearly_equal( p2 ));
        ldout << ", p3_on_p1_p2=" << p3_on_p1_p2_proper << "," << p3_on_p1_p2 << 
                  " p4_on_p1_p2=" << p4_on_p1_p2_proper << "," << p4_on_p1_p2;
        if ( p3_on_p1_p2 != p4_on_p1_p2 ) {
            real2 other;
            if ( p3_on_p1_p2 ) {
                ip = p3;
                other = p4;
            } else if ( p4_on_p1_p2 ) {
                ip = p4;
                other = p3;
                reverse = true;
            }
            real2 p1_p2_dir    = p2 - p1;
            real2 ip_other_dir = other - ip;
            if ( p1_p2_dir.colinear_is_same_dir( ip_other_dir ) ) {
                ldout << ", reverse=" << reverse << " => INTERSECTION (colinear same dir)\n";
                colinear = true;
                return true;
            } else {
                ldout << ", but wrong dir";
            }
        } else {
            ldout << ", p3_on_p1_p2=" << p3_on_p1_p2 << " p4_on_p1_p2=" << p4_on_p1_p2;
        }
    }
    ldout << "\n";
    return false;
}

void Layout::real2::pad_segment( const real2& p2, real pad, real2& p1_new, real2& p2_new ) const
{
    const real2& p1 = *this;

    real2 dxy = p1 - p2;
    dxy.normalize();
    p1_new = p1 + real2(  pad*dxy.c[0],  pad*dxy.c[1] );
    p2_new = p2 + real2( -pad*dxy.c[0], -pad*dxy.c[1] );
}

void Layout::real2::perpendicular_segment( const real2& p2, real length, real2& p3, real2& p4 ) const
{
    const real2& p1 = *this;

    real2 dxy = p1 - p2;
    dxy.normalize();
    real l2 = length / 2.0;
    p3 = p1 + real2(  l2*dxy.c[1], -l2*dxy.c[0] );
    p4 = p1 + real2( -l2*dxy.c[1],  l2*dxy.c[0] );
}

void Layout::real2::parallel_segment( const real2& p2, real dist, real2& p1_new, real2& p2_new ) const
{
    const real2& p1 = *this;
    const real2 v12 = p2 - p1;

    real2 p3;
    real2 p4;
    perpendicular_segment( p2, std::abs(dist) * 2.0, p3, p4 );
    p1_new = (dist >= 0.0) ? p3 : p4;
    p2_new = p1_new + v12;
}

inline Layout::real2 Layout::real2::operator + ( const Layout::real2& v2 ) const
{
    real2 r;
    r.c[0] = c[0] + v2.c[0];
    r.c[1] = c[1] + v2.c[1];
    return r;
}

inline Layout::real2 Layout::real2::operator - ( const Layout::real2& v2 ) const
{
    real2 r;
    r.c[0] = c[0] - v2.c[0];
    r.c[1] = c[1] - v2.c[1];
    return r;
}

inline Layout::real2 Layout::real2::operator * ( const Layout::real2& v2 ) const
{
    real2 r;
    r.c[0] = c[0] * v2.c[0];
    r.c[1] = c[1] * v2.c[1];
    return r;
}

inline Layout::real2 Layout::real2::operator * ( Layout::real s ) const
{
    real2 r;
    r.c[0] = c[0] * s;
    r.c[1] = c[1] * s;
    return r;
}

inline Layout::real2 Layout::real2::operator / ( const Layout::real2& v2 ) const
{
    real2 r;
    r.c[0] = c[0] / v2.c[0];
    r.c[1] = c[1] / v2.c[1];
    return r;
}

inline Layout::real2 Layout::real2::operator / ( Layout::real s ) const
{
    real2 r;
    r.c[0] = c[0] / s;
    r.c[1] = c[1] / s;
    return r;
}

inline Layout::real2& Layout::real2::operator += ( const Layout::real2 &v2 )
{
    c[0] += v2.c[0];
    c[1] += v2.c[1];
    return *this;
}

inline Layout::real2& Layout::real2::operator -= ( const Layout::real2 &v2 )
{
    c[0] -= v2.c[0];
    c[1] -= v2.c[1];
    return *this;
}

inline Layout::real2& Layout::real2::operator *= ( const Layout::real2 &v2 )
{
    c[0] *= v2.c[0];
    c[1] *= v2.c[1];
    return *this;
}

inline Layout::real2& Layout::real2::operator *= ( const Layout::real s )
{
    c[0] *= s;
    c[1] *= s;
    return *this;
}

inline Layout::real2& Layout::real2::operator /= ( const Layout::real2 &v2 )
{
    c[0] /= v2.c[0];
    c[1] /= v2.c[1];
    return *this;
}

inline Layout::real2& Layout::real2::operator /= ( const Layout::real s )
{
    c[0] /= s;
    c[1] /= s;
    return *this;
}

inline bool Layout::real2::operator == ( const Layout::real2 &v2 ) const
{
    return c[0] == v2.c[0] && c[1] == v2.c[1];  
}

inline bool Layout::real2::operator != ( const Layout::real2 &v2 ) const
{
    return c[0] != v2.c[0] || c[1] != v2.c[1];  
}

bool Layout::real2::nearly_equal( const real2 &v2, real epsilon ) const
{
    real2 diff = *this - v2;
    return diff.c[0] >= -epsilon && diff.c[0] <= epsilon && diff.c[1] >= -epsilon && diff.c[1] <= epsilon;
}

inline std::string Layout::real2::str( void ) const
{
    return "[" + ::str(c[0]) + "," + ::str(c[1]) + "]";
}

inline void Layout::Matrix::identity( void )
{
    for( uint i = 0; i < 4; i++ )
    {
        for( uint j = 0; j < 4; j++ )
        {
            m[i][j] = (i == j) ? 1 : 0;
        }
    }
}

inline void Layout::Matrix::translate( const real3& translation )
{
    m[0][3] += translation.c[0];
    m[1][3] += translation.c[1];
    m[2][3] += translation.c[2];
}

inline void Layout::Matrix::scale( const real3& scaling )
{
    m[0][0] *= scaling.c[0];
    m[1][1] *= scaling.c[1];
    m[2][2] *= scaling.c[2];
}

void Layout::Matrix::rotate_xy( double radians )
{
    if ( radians == 0.0 ) return;
    double c = cos( radians );
    double s = sin( radians );
    Matrix M2;
    M2.m[0][0] = c;
    M2.m[0][1] = s;
    M2.m[1][0] = -s;
    M2.m[1][1] = c;

    // order: *this = *this * M2  
    Matrix M1 = *this;
    M1.transform( M2, *this );
}

void Layout::Matrix::rotate_xz( double radians )
{
    if ( radians == 0.0 ) return;
    double c = cos( radians );
    double s = sin( radians );
    Matrix M2;
    M2.m[0][0] = c;
    M2.m[0][2] = s;
    M2.m[2][0] = -s;
    M2.m[2][2] = c;

    // order: *this = *this * M2  
    Matrix M1 = *this;
    M1.transform( M2, *this );
}

void Layout::Matrix::rotate_yz( double radians )
{
    if ( radians == 0.0 ) return;
    double c = cos( radians );
    double s = sin( radians );
    Matrix M2;
    M2.m[1][1] = c;
    M2.m[1][2] = -s;
    M2.m[2][1] = s;
    M2.m[2][2] = c;

    // order: *this = *this * M2  
    Matrix M1 = *this;
    M1.transform( M2, *this );
}

inline Layout::Matrix Layout::Matrix::operator + ( const Matrix& M ) const
{
    Matrix r;
    for( uint i = 0; i < 4; i++ )
    {
        for( uint j = 0; j < 4; j++ )
        {
            r.m[i][j] = this->m[i][j] + M.m[i][j];
        }
    }
    return r;
}

inline Layout::Matrix Layout::Matrix::operator - ( const Matrix& M ) const
{
    Matrix r;
    for( uint i = 0; i < 4; i++ )
    {
        for( uint j = 0; j < 4; j++ )
        {
            r.m[i][j] = this->m[i][j] - M.m[i][j];
        }
    }
    return r;
}

inline bool Layout::Matrix::operator == ( const Matrix& M ) const
{
    for( uint i = 0; i < 4; i++ )
    {
        for( uint j = 0; j < 4; j++ )
        {
            if ( this->m[i][j] != M.m[i][j] ) return false;
        }
    }
    return true;
}

inline bool Layout::Matrix::operator != ( const Matrix& M ) const
{
    for( uint i = 0; i < 4; i++ )
    {
        for( uint j = 0; j < 4; j++ )
        {
            if ( this->m[i][j] != M.m[i][j] ) return true;
        }
    }
    return false;
}

inline void Layout::Matrix::multiply( double s ) 
{
    for( uint i = 0; i < 4; i++ )
    {
        for( uint j = 0; j < 4; j++ )
        {
            m[i][j] *= s;
        }
    }
}

Layout::real4 Layout::Matrix::row( uint r ) const
{
    real4 v;
    for( uint32_t c = 0; c < 4; c++ ) 
    {
        v.c[c] = m[r][c];
    }
    return v;
}

Layout::real4 Layout::Matrix::column( uint c ) const
{
    real4 v;
    for( uint32_t r = 0; r < 4; r++ ) 
    {
        v.c[r] = m[r][c];
    }
    return v;
}

void Layout::Matrix::transform( const real4& v, real4& r ) const
{
    // order: r = *this * v
    for( uint i = 0; i < 4; i++ )
    {
        double sum = 0.0;               // use higher-precision here
        for( uint j = 0; j < 4; j++ )
        {
            double partial = m[i][j];
            partial *= v.c[j];
            sum += partial;
        }
        r.c[i] = sum;
    }
}

void Layout::Matrix::transform( const real3& v, real3& r, bool div_by_w ) const
{
    // order: r = *this * v
    if ( div_by_w ) {
        real4 v4 = real4( v.c[0], v.c[1], v.c[2], 1.0 );
        real4 r4;
        transform( v4, r4 );
        r.c[0] = r4.c[0];
        r.c[1] = r4.c[1];
        r.c[2] = r4.c[2];
        r /= r4.c[3];                   // w
    } else {
        for( uint i = 0; i < 3; i++ )
        {
            double sum = 0.0;               // use higher-precision here
            for( uint j = 0; j < 3; j++ )
            {
                double partial = m[i][j];
                partial *= v.c[j];
                sum += partial;
            }
            r.c[i] = sum;
        }
    }
}

void Layout::Matrix::transform( const Matrix& M2, Matrix& M3 ) const
{
    // order: M3 = *this * M2
    for( uint r = 0; r < 4; r++ )
    {
        for( uint c = 0; c < 4; c++ )
        {
            double sum = 0.0;
            for( int k = 0; k < 4; k++ )
            {
                double partial = m[r][k];
                partial *= M2.m[k][c];
                sum += partial;
            }
            M3.m[r][c] = sum;
        }
    }
}

void Layout::Matrix::transpose( Layout::Matrix& mt ) const
{
    for( int i = 0; i < 4; i++ )
    {
        for( int j = 0; j < 4; j++ )
        {
            mt.m[j][i] = m[i][j];
        }
    }
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
    return int(node.kind) >= 0 && uint32_t(node.kind) < GDSII_KIND_CNT;
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

inline uint Layout::node_last_i( const Node& node ) const
{
    lassert( node_is_parent( node ), "node_last_i: node is not a parent" );    
    uint last_i = NULL_I;
    for( uint child_i = node.u.child_first_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
    {
        if ( nodes[child_i].sibling_i == NULL_I ) break;

        last_i = child_i;
    }
    return last_i;
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

inline uint Layout::node_width( const Node& node ) const
{
    if ( node.kind == NODE_KIND::WIDTH ) {
        uint int_node_i = node.u.child_first_i;
        lassert( int_node_i != NULL_I && nodes[int_node_i].kind == NODE_KIND::INT && nodes[int_node_i].sibling_i == NULL_I, "bad WIDTH node" );
        return nodes[int_node_i].u.i;
    } else {
        lassert( node_is_element( node ), "node_width: node is not a WIDTH or an ELEMENT" );
        for( uint child_i = node.u.child_first_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i ) 
        {
            if ( nodes[child_i].kind == NODE_KIND::WIDTH ) {
                uint int_node_i = nodes[child_i].u.child_first_i;
                lassert( int_node_i != NULL_I && nodes[int_node_i].kind == NODE_KIND::INT && nodes[int_node_i].sibling_i == NULL_I, "bad WIDTH node" );
                return nodes[int_node_i].u.i;
            }
        }
        lassert( false, "could not find width for node" );
        return NULL_I;
    }
}

inline uint Layout::node_bah_layer( const Node& node ) const
{
    uint layer_i = node_layer( node );
    lassert( layer_i < hdr->layer_cnt, "node_layer() is out of range of current layer_cnt" );
    while( layers[layer_i].same_zoffset_as_prev ) 
    {
        //------------------------------------------------------------
        // Move to previous.
        //------------------------------------------------------------
        lassert( layer_i != 0, "layers[0].same_zoffset_as_prev should be false" );
        layer_i--;
    }
    return layer_i;
}

std::string Layout::all_struct_names( std::string delim ) const
{
    std::string names = "";
    for( auto it = name_i_to_struct_i.begin(); it != name_i_to_struct_i.end(); it++ )
    {
        uint name_i = it->first;    
        std::string name = std::string( &strings[name_i] );
        if ( names != "" ) names += delim;
        names += name;
    }
    return names;
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
        lassert( sit != name_i_to_struct_i.end(), "could not find structure with name " + std::string(&strings[name_i]) + 
                                                  ", available structures:\n    " + all_struct_names() );
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

void Layout::node_xy_replace_polygon( Node& node, const real2 * vtx, uint vtx_cnt ) // replace X,Y children from new set of vertices
{
    lassert( node.kind == NODE_KIND::XY, "node_xy_replace_polygon() not called with XY node" );
    uint dst_xy_prev_i = NULL_I;
    for( uint v = 0; v < vtx_cnt; v++ )
    {
        uint dst_x_i = node_alloc_int( vtx[v].c[0] / gdsii_units_user );        
        if ( dst_xy_prev_i == NULL_I ) {
            node.u.child_first_i = dst_x_i;
        } else {
            nodes[dst_xy_prev_i].sibling_i = dst_x_i;
        } 
        
        uint dst_y_i = node_alloc_int( vtx[v].c[1] / gdsii_units_user );
        nodes[dst_x_i].sibling_i = dst_y_i;
        dst_xy_prev_i = dst_y_i;
    }
}

inline uint Layout::node_pathtype( const Node& node ) const
{
    if ( node.kind == NODE_KIND::PATHTYPE ) {
        uint child_i = node.u.child_first_i;
        lassert( child_i != NULL_I, "PATHTYPE node has no child INT" );
        lassert( nodes[child_i].kind == NODE_KIND::INT, "PATHTYPE child node is not an INT" );
        return nodes[child_i].u.i;
    } else {    
        lassert( node.kind == NODE_KIND::PATH, "pathtype() should have been called on a PATH node" );
        for( uint child_i = node.u.child_first_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i ) 
        {
            if ( nodes[child_i].kind == NODE_KIND::PATHTYPE ) return nodes[child_i].u.i;
        }
    }
    return 0;   // default
}

inline uint Layout::node_datatype( const Node& node ) const
{
    if ( node.kind == NODE_KIND::DATATYPE ) {
        uint child_i = node.u.child_first_i;
        lassert( child_i != NULL_I, "DATATYPE node has no child INT" );
        lassert( nodes[child_i].kind == NODE_KIND::INT, "DATATYPE child node is not an INT" );
        return nodes[child_i].u.i;
    } else {    
        for( uint child_i = node.u.child_first_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i ) 
        {
            if ( nodes[child_i].kind == NODE_KIND::DATATYPE ) return nodes[child_i].u.i;
        }
    }
    return 0;   // default
}

inline Layout::AABR::AABR( const Layout::real2& p )
{
    min = p;
    max = p;
}  

inline void Layout::AABR::pad( Layout::real p ) 
{
    min -= real2( p, p );
    max += real2( p, p );
}

inline void Layout::AABR::expand( const Layout::AABR& other )
{
    for( uint i = 0; i < 2; i++ )
    {
        if ( other.min.c[i] < min.c[i] ) min.c[i] = other.min.c[i];
        if ( other.max.c[i] > max.c[i] ) max.c[i] = other.max.c[i];
    }
}

inline void Layout::AABR::expand( const Layout::real2& p ) 
{
    if ( p.c[0] < min.c[0] ) min.c[0] = p.c[0];
    if ( p.c[1] < min.c[1] ) min.c[1] = p.c[1];
    if ( p.c[0] > max.c[0] ) max.c[0] = p.c[0];
    if ( p.c[1] > max.c[1] ) max.c[1] = p.c[1];
}

inline bool Layout::AABR::encloses( const AABR& other ) const
{
    return min.c[0] <= other.min.c[0] &&
           min.c[1] <= other.min.c[1] &&
           max.c[0] >= other.max.c[0] &&
           max.c[1] >= other.max.c[1];
}

inline bool Layout::AABR::encloses( const real2& p ) const
{
    return min.c[0] <= p.c[0] &&
           min.c[1] <= p.c[1] &&
           max.c[0] >= p.c[0] &&
           max.c[1] >= p.c[1];
}

inline void Layout::AABR::intersect( const Layout::AABR& other )
{
    for( uint i = 0; i < 2; i++ )
    {
        if ( other.min.c[i] > min.c[i] ) min.c[i] = other.min.c[i];
        if ( other.max.c[i] < max.c[i] ) max.c[i] = other.max.c[i];
    }
}

inline bool Layout::AABR::intersects( const AABR& other ) const
{
    return !( min.c[0] >= other.max.c[0] ||
              min.c[1] >= other.max.c[1] ||
              max.c[0] <= other.min.c[0] ||
              max.c[1] <= other.min.c[1] );
}

inline Layout::AABR Layout::AABR::quadrant( uint i, uint j ) const
{
    AABR quad = *this;
    real2 half;
    half.c[0] = (quad.max.c[0] - quad.min.c[0]) / 2.0;
    half.c[1] = (quad.max.c[1] - quad.min.c[1]) / 2.0;
    if ( j == 0 ) {
        quad.max.c[0] -= half.c[0];
    } else {
        quad.min.c[0] += half.c[0];
    }
    if ( i == 0 ) {
        quad.max.c[1] -= half.c[1];
    } else {
        quad.min.c[1] += half.c[1];
    }
    return quad;
}

inline Layout::AABR Layout::AABR::enclosing_square( real scale_factor ) const
{
    AABR square = *this;
    real w  = max.c[0] - min.c[0];
    real h  = max.c[1] - min.c[1];
    real wh = std::max( w, h ) * scale_factor;  // side of square, scaled
    square.max.c[0] = min.c[0] + wh;
    square.max.c[1] = min.c[1] + wh;
    return square;
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
    node_alloc_timestamp( nodes[ni] );
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
    ni2 = node_alloc_real( units_user );
    nodes[ni].u.child_first_i = ni2;
    prev_i = ni2;
    ni2 = node_alloc_real( units_meters );
    nodes[prev_i].sibling_i = ni2;
    prev_i = ni;

    return prev_i;
}

uint Layout::inst_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, real x, real y, std::string name )
{
    return inst_layout( parent_i, last_i, src_layout, src_struct_name, x, y, 0, hdr->layer_cnt-1, name );
}

uint Layout::inst_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, 
                          real x, real y, uint dst_layer_first, uint dst_layer_last, std::string name )
{
    Matrix M;
    M.translate( real3( x, y, 0 ) );
    return inst_layout( parent_i, last_i, src_layout, src_struct_name, M, dst_layer_first, dst_layer_last, name );
}

uint Layout::inst_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string src_struct_name, 
                          const Matrix& M, uint dst_layer_first, uint dst_layer_last, std::string name )
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
        ldout << "inst_layout: dst_layer=" << i << " src_layer=" << layers[i].gdsii_num << " inst_name=" << inst_name << "\n";
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
        lassert( name_i_to_struct_i.find( dst_top_struct_name_i ) != name_i_to_struct_i.end(), "could not find top struct with name " + dst_top_struct_name + 
                                                                                               ", available structures:\n    " + all_struct_names() );

        uint dst_top_struct_i = name_i_to_struct_i[dst_top_struct_name_i];
        Node& dst_top_node = nodes[dst_top_struct_i];

        perhaps_realloc( top_insts, hdr->top_inst_cnt, max->top_inst_cnt, 1 );
        uint ii = hdr->top_inst_cnt++;
        top_insts[ii] = TopInstInfo{ dst_top_struct_i, M };
    }
    return last_i;
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
            lassert( dst_int_node.kind == NODE_KIND::INT && dst_int_node.u.i == int(src_layer_num) && dst_int_node.sibling_i == NULL_I, "unexpected layer_num in LAYER node" );
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
                uint dst_child_i = inst_layout_node( dst_i, dst_prev_i, src_layout, src_struct_name, src_child_i, src_layer_num, dst_layer_num, cache, name, indent_str + "    " );
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
                uint dst_child_i = inst_layout_node( parent_i, last_i, src_layout, src_struct_name, src_child_i, src_layer_num, dst_layer_num, cache, name, indent_str + "    " );
                if ( dst_child_i != NULL_I ) last_i = dst_child_i;
            } else {
                ldout << indent_str << "    " << "skipping " << str(src_layout->nodes[src_child_i].kind) << "\n";
            }
        }
    }

    return last_i;
}

void Layout::finalize_top_struct( uint parent_i, uint last_i, std::string top_name )
{
    //-----------------------------------------------------
    // Create one final struct that instantiates all top_insts.
    //-----------------------------------------------------
    lassert( parent_i == NULL_I, "finalize parent_i must be NULL_I for now" );
    uint bgnstr_i = node_alloc( NODE_KIND::BGNSTR );
    nodes[last_i].sibling_i = bgnstr_i;
    node_alloc_timestamp( nodes[bgnstr_i] );
    uint prev_i = node_last_scalar_i( nodes[bgnstr_i] );

    uint strname_i = node_alloc( NODE_KIND::STRNAME );
    nodes[prev_i].sibling_i = strname_i;
    nodes[strname_i].u.s_i = str_get( top_name );
    prev_i = strname_i;

    name_i_to_struct_i[nodes[strname_i].u.s_i] = bgnstr_i;

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

        real x = info.M.m[0][3];
        real y = info.M.m[1][3];
        Matrix Mt;
        Mt.translate( real3( x, y, 0 ) );
        lassert( info.M == Mt, "top-level instances must have only an XY translation for now; no other transformations at the top" );
        uint x_i = node_alloc( NODE_KIND::INT );
        nodes[xy_i].u.child_first_i = x_i;
        nodes[x_i].u.i = int( x / gdsii_units_user );      
        uint y_i = node_alloc( NODE_KIND::INT );
        nodes[x_i].sibling_i = y_i;
        nodes[y_i].u.i = int( y / gdsii_units_user );
    }
}

uint Layout::flatten_layout( uint parent_i, uint last_i, const Layout * src_layout, std::string struct_name, 
                             CONFLICT_POLICY conflict_policy, const Matrix& M )
{
    //------------------------------------------------------------
    // Locate the the src struct.
    //------------------------------------------------------------
    uint src_struct_name_i = src_layout->str_find( struct_name );
    auto it = src_layout->name_i_to_struct_i.find( src_struct_name_i );
    lassert( it != src_layout->name_i_to_struct_i.end(), "no src struct with the name " + struct_name + ", available structures:\n    " + src_layout->all_struct_names() );
    uint src_struct_i = it->second;

    //------------------------------------------------------------
    // Recursively copy the src_struct.
    //------------------------------------------------------------
    uint dst_struct_i = node_copy( parent_i, last_i, src_layout, src_struct_i, COPY_KIND::FLATTEN, conflict_policy, M );
    lassert( nodes[dst_struct_i].kind == NODE_KIND::BGNSTR, "should have gotten back a BGNSTR node after flatten-copy" );

    uint dst_name_i = str_get( struct_name );
    lassert( name_i_to_struct_i.find( dst_name_i ) == name_i_to_struct_i.end(), "struct name " + struct_name + " should not exist yet in flattened layout" );
    name_i_to_struct_i[dst_name_i] = dst_struct_i;

    return dst_struct_i;
}

Layout::AABR Layout::bounding_rect( uint layer_first, uint layer_last ) const
{
    //-----------------------------------------------------
    // Take union of all layers.
    //-----------------------------------------------------
    layer_last = std::min( layer_last, hdr->layer_cnt-1 );
    AABR brect;
    bool have_one = false;
    for( uint i = layer_first; i <= layer_last; i++ )
    {
        if ( layers[i].bah_root_i == NULL_I ) continue;

        if ( !have_one ) {
            brect = layers[i].bah_brect;
            have_one = true;
        } else {
            brect.expand( layers[i].bah_brect );
        }
    }
    lassert( have_one, "found no bounding area hierarchy (BAH) so could not determine bounding rectangle" );
    return brect;
}

uint Layout::fill_dielectrics( uint parent_i, const AABR& brect, uint layer_first, uint layer_last )
{
    //-----------------------------------------------------
    // For each layer that is not part of previous layer:
    //-----------------------------------------------------
    lassert( nodes[parent_i].kind == NODE_KIND::BGNSTR, "fill_dielectrics first argument should be index of BGNSTR to fill" );
    uint last_i = node_last_i( nodes[parent_i] );
    layer_last = std::min( layer_last, hdr->layer_cnt-1 );
    for( uint32_t i = layer_first; i <= layer_last; i++ )
    {
        if ( layers[i].same_zoffset_as_prev ) continue;
        lassert( layers[i].bah_root_i != NULL_I, "layer " + std::to_string(i) + " has no bounding area hierarchy (BAH); is it flattened?" );

        //-----------------------------------------------------
        // First, fill any space that's inside the overall brect but
        // outside the layer's bah_brect.
        //-----------------------------------------------------
        if ( brect.min.c[0] < layers[i].bah_brect.min.c[0] ) {
            // left
            AABR rect = brect;
            rect.max.c[0] = layers[i].bah_brect.min.c[0];
            last_i = fill_dielectric_rect( parent_i, last_i, i, rect );
        }
        if ( brect.max.c[0] > layers[i].bah_brect.max.c[0] ) {
            // right
            AABR rect = brect;
            rect.min.c[0] = layers[i].bah_brect.max.c[0];
            last_i = fill_dielectric_rect( parent_i, last_i, i, rect );
        }
        if ( brect.min.c[1] < layers[i].bah_brect.min.c[1] ) {
            // top
            AABR rect = layers[i].bah_brect;
            rect.min.c[1] = brect.min.c[1];
            rect.max.c[1] = layers[i].bah_brect.min.c[1];
            last_i = fill_dielectric_rect( parent_i, last_i, i, rect );
        }
        if ( brect.max.c[1] > layers[i].bah_brect.max.c[1] ) {
            // bottom
            AABR rect = layers[i].bah_brect;
            rect.min.c[1] = layers[i].bah_brect.max.c[1];
            rect.max.c[1] = brect.max.c[1];
            last_i = fill_dielectric_rect( parent_i, last_i, i, rect );
        }

        //-----------------------------------------------------
        // Next, walk the entire BAH for this layer.
        // We'll fill empty space at each leaf.
        //-----------------------------------------------------
        last_i = fill_dielectric_bah( parent_i, last_i, i, layers[i].bah_root_i, layers[i].bah_brect );
    }
    return last_i;
}

uint Layout::fill_dielectric_bah( uint parent_i, uint last_i, uint layer_i, uint bah_i, const AABR& rect )
{
    //-----------------------------------------------------
    // Fill each quadrant.
    //-----------------------------------------------------
    lassert( bah_i != NULL_I, "fill_dielectric_bah: bah_i should not be NULL_I" );
    const BAH_Node& node = bah_nodes[bah_i];
    AABR crect;
    real2 mid( (rect.min.c[0] + rect.max.c[0]) / 2.0, (rect.min.c[1] + rect.max.c[1]) / 2.0 );
    for( uint x = 0; x < 2; x++ ) 
    {
        crect.min.c[0] = (x == 0) ? rect.min.c[0] : mid.c[0];
        crect.max.c[0] = (x == 0) ? mid.c[0]      : rect.max.c[0];
        for( uint y = 0; y < 2; y++ ) 
        {
            crect.min.c[1] = (y == 0) ? rect.min.c[1] : mid.c[1];
            crect.max.c[1] = (y == 0) ? mid.c[1]      : rect.max.c[1];
            if ( node.child_is_leaf[x][y] ) {
                // leaf could be empty or nonempty
                last_i = fill_dielectric_leaf( parent_i, last_i, layer_i, bah_i, x, y, crect );
            } else {
                // recurse to child bah node
                last_i = fill_dielectric_bah( parent_i, last_i, layer_i, node.child_i[x][y], crect );
            }
        }
    }
    return last_i;
}

uint Layout::fill_dielectric_leaf( uint parent_i, uint last_i, uint layer_i, uint bah_i, uint x, uint y, const AABR& rect )
{
    uint li = bah_nodes[bah_i].child_i[x][y];
    if ( li == NULL_I ) {
        //-----------------------------------------------------
        // Add new rectangular element.
        //-----------------------------------------------------
        last_i = fill_dielectric_rect( parent_i, last_i, layer_i, rect );
    } else {
        //-----------------------------------------------------
        // Look up XY node in leaf's node.
        // Convert it to simple polygon format.
        //-----------------------------------------------------
        Leaf_Node& leaf = leaf_nodes[li];
        const Node& node = nodes[leaf.node_i];
        lassert( node_is_element( node ), "leaf is not an element" );    
        uint xy_i = node_xy_i( node );
        lassert( xy_i != NULL_I, "leaf has no XY node" );
        Node& xy = nodes[xy_i];
        uint vtx_cnt;
        real2 * vtx = polygon_alloc( xy, vtx_cnt );

        //-----------------------------------------------------
        // Find the vertex in the polygon that is closest to 
        // our rect.min corner.  
        //-----------------------------------------------------
        uint best_v = 0;
        real best_dist_sqr = 1e100;
        for( uint v = 0; v < (vtx_cnt-1); v++ ) 
        {
            real this_dist_sqr = (rect.min - vtx[v]).length_sqr();
            if ( this_dist_sqr < best_dist_sqr ) {
                best_v = v;
                best_dist_sqr = this_dist_sqr;
            }
        }

        //-----------------------------------------------------
        // Start with a line segment from the min corner to 
        // the closest vertex.  Then go all the way around the polygon
        // back to the same vertex, then back to the min corner.
        //-----------------------------------------------------
        uint vtx2_cnt = vtx_cnt + 5;  
        real2 * vtx2 = polygon_alloc( vtx2_cnt );
        uint v2 = 0;
        vtx2[v2++] = rect.min;
        for( uint i = 0; i < (vtx_cnt-1); i++ )
        {
            uint v3 = (best_v + i) % (vtx_cnt-1);
            vtx2[v2++] = vtx[v3];
        }
        vtx2[v2++] = rect.min;

        delete[] vtx;

        //-----------------------------------------------------
        // Now trace the outer rectangle starting with rect.min.
        //-----------------------------------------------------
        vtx2[v2++] = rect.min;
        vtx2[v2++] = real2( rect.min.c[0], rect.max.c[1] );
        vtx2[v2++] = rect.max;
        vtx2[v2++] = real2( rect.max.c[0], rect.min.c[1] );
        vtx2[v2++] = rect.min;
        lassert( v2 == vtx2_cnt, "something is wrong" );

        //-----------------------------------------------------
        // Add new polygonal element.
        //-----------------------------------------------------
        uint datatype = node_datatype( node );
        last_i = fill_dielectric_polygon( parent_i, last_i, layer_i, datatype, vtx2, vtx2_cnt );

        delete[] vtx2;
    }

    return last_i;
}

uint Layout::fill_dielectric_rect( uint parent_i, uint last_i, uint layer_i, const AABR& rect )
{
    //-----------------------------------------------------
    // Add a dielectric rectangle on the given layer.
    // The caller must ensure that it doesn't overlap anything.
    //-----------------------------------------------------
    const real2 vtx[] = { rect.min, 
                          real2( rect.min.c[0], rect.max.c[1] ),
                          rect.max,
                          real2( rect.max.c[0], rect.min.c[1] ),
                          rect.min };
    last_i = fill_dielectric_polygon( parent_i, last_i, layer_i, NULL_I, vtx, 5 ); 
    return last_i;
}

uint Layout::fill_dielectric_polygon( uint parent_i, uint last_i, uint layer_i, uint datatype, const real2 * vtx, uint vtx_cnt )
{
    //-----------------------------------------------------
    // Create new BOUNDARY node.
    // LAYER is dielectric layer.
    // DATATYPE is passed in.
    // XY contains polygon vertices.
    //
    // There's already a utility routine to do this.
    //-----------------------------------------------------
    last_i = node_alloc_boundary( parent_i, last_i, layers[layer_i].dielectric_gdsii_num, datatype, vtx, vtx_cnt );
    return last_i;
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

inline uint Layout::node_alloc_int( _int i )
{
    uint ni = node_alloc( NODE_KIND::INT );
    nodes[ni].u.i = i;
    return ni;
}

inline uint Layout::node_alloc_real( real r )
{
    uint ni = node_alloc( NODE_KIND::REAL );
    nodes[ni].u.r = r;
    return ni;
}

inline uint Layout::node_alloc_str( std::string s )
{
    uint ni = node_alloc( NODE_KIND::STR );
    nodes[ni].u.s_i = str_get( s );
    return ni;
}

uint Layout::node_alloc_boundary( uint parent_i, uint last_i, uint gdsii_num, uint datatype, const real2 * vtx, uint vtx_cnt )
{
    // BOUNDARY
    uint boundary_i = node_alloc( NODE_KIND::BOUNDARY );
    if ( last_i != NULL_I ) {
        lassert( nodes[last_i].sibling_i == NULL_I, "node_alloc_boundary last_i sibling_i should not be set" );
        nodes[last_i].sibling_i = boundary_i;
    }

    // LAYER
    lassert( gdsii_num != NULL_I, "BOUNDARY must have a non-null LAYER gdsii_num" );
    uint layer_i = node_alloc( NODE_KIND::LAYER );
    nodes[boundary_i].u.child_first_i = layer_i;
    nodes[layer_i].u.child_first_i = node_alloc_int( gdsii_num );
    last_i = layer_i;

    if ( datatype != NULL_I ) {
        // DATATYPE
        uint datatype_i = node_alloc( NODE_KIND::DATATYPE );
        nodes[last_i].sibling_i = datatype_i;
        nodes[datatype_i].u.child_first_i = node_alloc_int( datatype );    
        last_i = datatype_i;
    }

    // XY
    uint xy_i = node_alloc( NODE_KIND::XY );
    nodes[last_i].sibling_i = xy_i;

    uint xy_prev_i = NULL_I;
    for( uint i = 0; i < vtx_cnt; i++ )
    {
        uint x_i = node_alloc_int( vtx[i].c[0] / gdsii_units_user );        
        if ( xy_prev_i == NULL_I ) {
            nodes[xy_i].u.child_first_i = x_i;
        } else {
            nodes[xy_prev_i].sibling_i = x_i;
        } 
        
        uint y_i = node_alloc_int( vtx[i].c[1] / gdsii_units_user );
        nodes[x_i].sibling_i = y_i;
        xy_prev_i = y_i;
    }

    return boundary_i;
}

void Layout::node_alloc_timestamp( Node& node )
{
    lassert( node.kind == NODE_KIND::BGNLIB || node.kind == NODE_KIND::BGNSTR, "node_alloc_timestamp: node must be BGNLIB or BGNSTR" );
    lassert( node.u.child_first_i == NULL_I, "node_alloc_timestamp found node with children already" );

    time_t t;
    time( &t );
    struct tm * tm = gmtime( &t );

    uint prev_i = NULL_I;
    for( uint i = 0; i < 2; i++ ) 
    {
        uint ni = node_alloc_int( tm->tm_year );
        if ( i == 0 ) {
            node.u.child_first_i = ni;
        } else {
            nodes[prev_i].sibling_i = ni;
        }
        prev_i = ni;

        ni = node_alloc_int( tm->tm_mon );
        nodes[prev_i].sibling_i = ni;
        prev_i = ni;

        ni = node_alloc_int( tm->tm_mday );
        nodes[prev_i].sibling_i = ni;
        prev_i = ni;

        ni = node_alloc_int( tm->tm_hour );
        nodes[prev_i].sibling_i = ni;
        prev_i = ni;

        ni = node_alloc_int( tm->tm_min );
        nodes[prev_i].sibling_i = ni;
        prev_i = ni;

        ni = node_alloc_int( tm->tm_sec );
        nodes[prev_i].sibling_i = ni;
        prev_i = ni;
    }
}

uint Layout::node_copy( uint parent_i, uint last_i, const Layout * src_layout, uint src_i, COPY_KIND copy_kind, 
                        CONFLICT_POLICY conflict_policy, const Matrix& M, bool in_flatten )
{
    const Node& src_node = src_layout->nodes[src_i];
    if ( copy_kind == COPY_KIND::FLATTEN ) {
        //-----------------------------------------------------
        // Flattening handles certain nodes differently than hier copies.
        //-----------------------------------------------------
        switch( src_node.kind )
        {
            case NODE_KIND::SREF:
            case NODE_KIND::AREF:
                return node_flatten_ref( parent_i, last_i, src_layout, src_i, conflict_policy, M );

            case NODE_KIND::STRNAME:
                if ( in_flatten ) return NULL_I;  // blow this off
                break;

            default:
                break;
        }
    }

    uint dst_first_i;
    if ( in_flatten && src_node.kind == NODE_KIND::PATH ) {
        //-----------------------------------------------------
        // Convert PATH to BOUNDARY     
        //-----------------------------------------------------
        dst_first_i = node_convert_path_to_boundary( parent_i, last_i, src_layout, src_i, conflict_policy, M );
    } else {
        //-----------------------------------------------------
        // Normal Copy
        //-----------------------------------------------------
        dst_first_i = node_alloc( src_node.kind );
        if ( src_layout->node_is_parent( src_node ) ) { 
            if ( copy_kind != COPY_KIND::ONE ) {
                //-----------------------------------------------------
                // Copy children.
                //-----------------------------------------------------
                uint dst_prev_i = NULL_I;
                for( src_i = src_node.u.child_first_i; src_i != NULL_I; src_i = src_layout->nodes[src_i].sibling_i )
                {
                    const Node& src = src_layout->nodes[src_i];

                    if ( copy_kind == COPY_KIND::DEEP || copy_kind == COPY_KIND::FLATTEN || node_is_scalar( src ) ) {
                        dst_prev_i = node_copy( dst_first_i, dst_prev_i, src_layout, src_i, copy_kind, conflict_policy, M, in_flatten );
                    } else {
                        break;
                    }
                }
            }
        } else if ( node_is_string( src_node ) ) {
            // STR
            nodes[dst_first_i].u.s_i = str_get( &src_layout->strings[src_node.u.s_i] );

        } else {
            // INT, REAL
            nodes[dst_first_i].u = src_node.u;
        }
    }

    //-----------------------------------------------------
    // For each dst_i (could be more than one for PATH to BOUNDARY conversion):
    //-----------------------------------------------------
    for( uint dst_i = dst_first_i; dst_i != NULL_I; dst_i = nodes[dst_i].sibling_i )
    {
        //-----------------------------------------------------
        // Connect the new dst_i node to parent_i or last_i.
        //-----------------------------------------------------
        if ( last_i != NULL_I ) {
            lassert( nodes[last_i].sibling_i == NULL_I || nodes[last_i].sibling_i == dst_i, "node_copy: nodes[last_i] sibling_i is already set" ); 
            lassert( parent_i == NULL_I || nodes[parent_i].u.child_first_i != NULL_I, "node_copy: nodes[parent_i] is not set but last_i is set" );
            nodes[last_i].sibling_i = dst_i;
        } else if ( parent_i != NULL_I ) {
            lassert( nodes[parent_i].u.child_first_i == NULL_I, "node_copy: nodes[parent_i] child_first_i is already set but last_i is NULL_I" ); 
            nodes[parent_i].u.child_first_i = dst_i;
        }
        last_i = dst_i;

        if ( nodes[dst_i].kind == NODE_KIND::XY ) {
            //-----------------------------------------------------
            // Transform each X,Y pair by M.
            // Also keep track of the bounding rectangle if we're flattening.
            //-----------------------------------------------------
            node_transform_xy( parent_i, dst_i, copy_kind, conflict_policy, M );
        }
    }

    return last_i;
}

uint Layout::node_flatten_ref( uint parent_i, uint last_i, const Layout * src_layout, uint src_i, CONFLICT_POLICY conflict_policy, const Matrix& M )
{
    //-----------------------------------------------------
    // First extract the instancing params.
    //-----------------------------------------------------
    const Node * src_nodes = src_layout->nodes;
    const Node&  src_node  = src_nodes[src_i];
    uint sname_i = NULL_I;
    uint struct_i = NULL_I;
    uint strans = 0;
    real smag = 1.0;
    real sangle = 0.0;
    bool abs_smag = false;
    bool abs_sangle = false;
    bool sreflection = false;
    uint col_cnt = 1;
    uint row_cnt = 1;
    real xy[3][2] = { {0, 0}, {0, 0}, {0, 0} };
    for( uint child_i = src_node.u.child_first_i; child_i != NULL_I; child_i = src_nodes[child_i].sibling_i )
    {
        const Node& child = src_nodes[child_i];
        switch( child.kind )
        {
            case NODE_KIND::SNAME:
            {
                lassert( sname_i == NULL_I, "SREF/AREF has duplicate SNAME child" );
                sname_i = child.u.s_i;
                auto it = src_layout->name_i_to_struct_i.find( sname_i );
                lassert( it != src_layout->name_i_to_struct_i.end(), "SREF/AREF SNAME " + std::string(&src_layout->strings[sname_i]) + 
                                                                     " does not denote a known struct; available structs are:\n" + all_struct_names() );
                struct_i = it->second;
                break;
            }

            case NODE_KIND::STRANS:
            {
                //-----------------------------------------------------
                // Bit 0 (the leftmost bit) specifies reflection. If it is set, then reflection 
                // about the X-axis is applied before angular rotation. For AREFs, the entire array lattice is reflected, 
                // with the individual array elements rigidly attached. 
                // 
                // Bit 13 flags absolute magnification. 
                // Bit 14 flags absolute angle. 
                //
                // Bit 15 (the rightmost bit) and all remaining bits are reserved for future use and must be cleared. 
                //-----------------------------------------------------
                strans = child.u.u;
                sreflection  = (strans >> 0)  & 1;
                sreflection |= (strans >> 7)  & 1;  // hack
                abs_smag     = (strans >> 13) & 1;
                abs_sangle   = (strans >> 14) & 1;
                ldout << "strans=" << strans << " sreflection=" << sreflection << 
                         " abs_smag=" << abs_smag << " abs_sangle=" << abs_sangle << "\n";
                break;
            }

            case NODE_KIND::MAG:
            {
                uint gchild_i = child.u.child_first_i;
                lassert( src_nodes[gchild_i].kind == NODE_KIND::REAL, "MAG node does not have a child that is a REAL" );
                smag = src_nodes[gchild_i].u.r;
                break;
            }

            case NODE_KIND::ANGLE:
            {
                //-----------------------------------------------------
                // For an AREF, the ANGLE rotates the entire array lattice (with the individual 
                // array elements rigidly attached) about the array reference point. 
                // If this record is omitted, an angle of zero degrees is as- sumed.
                //
                // Convert from degrees to radians.
                //-----------------------------------------------------
                uint gchild_i = child.u.child_first_i;
                lassert( src_nodes[gchild_i].kind == NODE_KIND::REAL, "ANGLE node does not have a child that is a REAL" );
                ldout << "raw angle=" << src_nodes[gchild_i].u.r << "\n";
                sangle = src_nodes[gchild_i].u.r * real(M_PI) / 180.0;
                break;
            }

            case NODE_KIND::COLROW:
            {
                //-----------------------------------------------------
                // Array dimensions
                //-----------------------------------------------------
                lassert( src_node.kind == NODE_KIND::AREF, "COLROW not allowed for an SREF" );
                uint gchild_i = child.u.child_first_i;
                lassert( gchild_i != NULL_I, "COLROW has no COL value" );
                lassert( src_nodes[gchild_i].kind == NODE_KIND::INT, "COLROW COL is not an INT" );
                col_cnt = src_nodes[gchild_i].u.i;
                lassert( col_cnt > 0, "COLROW COL must be non-zero" );

                gchild_i = src_nodes[gchild_i].sibling_i;
                lassert( gchild_i != NULL_I, "COLROW has no ROW value" );
                lassert( src_nodes[gchild_i].kind == NODE_KIND::INT, "COLROW ROW is not an INT" );
                row_cnt = src_nodes[gchild_i].u.i;
                lassert( row_cnt > 0, "COLROW ROW must be non-zero" );
                break;
            }
                 
            case NODE_KIND::XY:
            {
                //-----------------------------------------------------
                // A text or SREF element must have only one pair of coordinates.
                //
                // An AREF has exactly three pairs of coordinates, which specify the orthogonal array lattice. 
                // In an AREF the first point is the array reference point. 
                // The second point locates a position which is dis- placed from the reference point by 
                // the inter-column spacing times the number of columns. 
                // The third point locates a position which is displaced from the reference 
                // point by the inter-row spacing times the number of rows.
                //-----------------------------------------------------
                uint i = 0;
                for( uint gchild_i = child.u.child_first_i; gchild_i != NULL_I; gchild_i = src_nodes[gchild_i].sibling_i )
                {
                    lassert( i < 2 || src_node.kind == NODE_KIND::AREF, "SREF may not have more than 2 XY coords" );
                    lassert( i < 6, "AREF may not have more than 6 XY coords" );
                    xy[i>>1][i&1] = real(src_nodes[gchild_i].u.i) * src_layout->gdsii_units_user;
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
    if ( src_node.kind == NODE_KIND::AREF ) {
        smag = 1.0;
        sangle = 0.0;
        sreflection = false;
    }
    lassert( struct_i != NULL_I, "SREF/AREF has no SNAME" );

    for( uint r = 0; r < row_cnt; r++ )
    {
        for( uint c = 0; c < col_cnt; c++ )
        {
            //-----------------------------------------------------
            // Apply transformations to previous ones.
            //-----------------------------------------------------
            Matrix inst_M = M;

            if ( smag != 1.0 || sreflection ) {
                real3 sv{ smag, sreflection ? -smag : smag, 1 };
                inst_M.scale( sv );
                ldout << "scale mag=" << smag << " sreflection=" << sreflection << " vec=" << sv << " new M=\n" << inst_M << "\n";
            }

            if ( sangle != 0.0 ) {
                inst_M.rotate_xy( sangle );
                ldout << "rotate angle=" << sangle << " new M=\n" << inst_M << "\n";
            }

            if ( src_node.kind == NODE_KIND::SREF ) {
                real3 tv{ xy[0][0], xy[0][1], 0 };
                inst_M.translate( tv );
                ldout << "sref translate vec=" << tv << " M=\n" << inst_M << "\n";
            } else {
                real dxy[2][2] = { { (xy[1][0] - xy[0][0]) / real(col_cnt), (xy[1][1] - xy[0][1]) / real(col_cnt) },
                                   { (xy[2][0] - xy[0][0]) / real(row_cnt), (xy[2][1] - xy[0][1]) / real(row_cnt) } };
                real rr = r;
                real cc = c;
                real x_translate = xy[0][0] + cc*dxy[0][0] + rr*dxy[1][0];
                real y_translate = xy[0][1] + rr*dxy[1][1] + cc*dxy[0][1];
                real3 tv{ x_translate, y_translate, 0 };
                inst_M.translate( real3(x_translate, y_translate, 0) );
                ldout << "aref translate vec=" << tv << " M=" << inst_M << "\n";
            }

            //-----------------------------------------------------
            // Copy the structure's children with new transformations.
            // But first skip the timestamp which we do not need.
            //-----------------------------------------------------
            uint src_child_i = src_nodes[struct_i].u.child_first_i;
            for( uint j = 0; j < 12; j++ )
            {
                lassert( src_child_i != NULL_I && src_nodes[src_child_i].kind == NODE_KIND::INT, 
                         "structure " + std::string(&src_layout->strings[sname_i]) + " timestamp ended prematurely" );
                src_child_i = src_nodes[src_child_i].sibling_i;
            }
            for( ; src_child_i != NULL_I; src_child_i = src_nodes[src_child_i].sibling_i )
            {
                uint dst_child_i = node_copy( parent_i, last_i, src_layout, src_child_i, COPY_KIND::FLATTEN, conflict_policy, inst_M, true );
                if ( dst_child_i != NULL_I ) last_i = dst_child_i;
            }
        }
    }
    return last_i;
}

uint Layout::node_convert_path_to_boundary( uint parent_i, uint last_i, const Layout * src_layout, uint src_i,
                                            CONFLICT_POLICY conflict_policy, const Matrix& M )
{
    //-----------------------------------------------------
    // PATH -> BOUNDARY conversion
    // 
    // First record the width, pathtype, etc.
    //-----------------------------------------------------
    lassert( src_layout->nodes[src_i].kind == NODE_KIND::PATH, "node_convert_path_to_boundary not called on PATH node" );
    real width = 0;
    uint pathtype = 0;
    uint datatype = NULL_I;
    uint layer = NULL_I;
    uint dst_first_i = NULL_I;
    uint dst_prev_i = NULL_I;
    for( src_i = src_layout->nodes[src_i].u.child_first_i; src_i != NULL_I; src_i = src_layout->nodes[src_i].sibling_i )
    {
        const Node& src = src_layout->nodes[src_i];

        switch( src.kind )
        {
            case NODE_KIND::WIDTH:              width    = src_layout->node_width( src ) * gdsii_units_user; break;
            case NODE_KIND::PATHTYPE:           pathtype = src_layout->node_pathtype( src );    break;
            case NODE_KIND::DATATYPE:           datatype = src_layout->node_datatype( src );    break;
            case NODE_KIND::LAYER:              layer    = src_layout->node_layer( src );       break;
            case NODE_KIND::XY:
            {
                //-----------------------------------------------------
                // Add each line segment separately.
                //
                // For pathtype 0:
                //     Get parallel segments width/2 from the src segment on each side.
                //     Add a rectangle through these points.
                //
                // For pathtype 1:
                //     Pad the src segment by width/2 on each end.
                //     Goto pathtype 1.
                //
                // For pathtype 2:
                //     Get parallel segments as for pathtype 0.
                //     Then extend the src segment to get the two arc endpoints.
                //     Add 4 arcs from arc endpoints to parallel segment endpoints.
                //     The arcs are tesselated.
                //
                // Remember to add a final segment from the last point back to the first of the new BOUNDARY.
                //-----------------------------------------------------
                real2 p0;
                real2 p1;
                bool have_x = false;
                bool have_p0 = false;
                for( uint src_xy_i = src.u.child_first_i; src_xy_i != NULL_I; src_xy_i = src_layout->nodes[src_xy_i].sibling_i )
                {
                    const Node& src_xy = src_layout->nodes[src_xy_i];
                    lassert( src_xy.kind == NODE_KIND::INT, "XY coordinate should be an INT" );
                    real xy_r = real(src_xy.u.i) * src_layout->gdsii_units_user;
                    if ( !have_x ) {
                        p1.c[0] = xy_r;
                        have_x = true;
                    } else {    
                        p1.c[1] = xy_r;
                        if ( have_p0 ) {
                            // see comments above 
                            lassert( pathtype >= 0 && pathtype <= 2, "pathtype is out of range 0 .. 2" );      
                            lassert( width >= 0.0, "PATH WIDTH should be >= 0.0" );
                            real2 s0 = p0;
                            real2 s1 = p1;
                            if ( pathtype == 1 ) p0.pad_segment( p1, width/2.0, s0, s1 );        
                            real2 s2, s3, s4, s5;
                            s0.parallel_segment( s1,  width/2.0, s2, s3 );
                            s0.parallel_segment( s1, -width/2.0, s4, s5 );
                            if ( pathtype == 2 ) p0.pad_segment( p1, width/2.0, s0, s1 );        

                            real2 vertex[14];
                            uint  c = 0;
                            vertex[c++] = s2;
                            vertex[c++] = s3;
                            if ( pathtype == 2 ) vertex[c++] = s1;
                            vertex[c++] = s5;
                            vertex[c++] = s4;
                            if ( pathtype == 2 ) vertex[c++] = s0;
                            vertex[c++] = s2;  // close the loop
                            lassert( c <= 14, "oops" );

                            //-----------------------------------------------------
                            // We have what we need to make the per-segment BOUNDARY,
                            // so do it.
                            //-----------------------------------------------------
                            uint dst_i = node_alloc_boundary( parent_i, dst_prev_i, layer, datatype, vertex, c );
                            if ( dst_first_i == NULL_I ) dst_first_i = dst_i;
                            dst_prev_i = dst_i;
                            uint dst_xy_i = node_xy_i( nodes[dst_i] );

                            //-----------------------------------------------------
                            // Transform XY pairs using M.
                            // Also add to BAH.
                            //-----------------------------------------------------
                            node_transform_xy( dst_i, dst_xy_i, COPY_KIND::FLATTEN, conflict_policy, M );
                        }

                        p0 = p1;
                        have_x = false;
                        have_p0 = true; 
                    } 
                }
                break;
            }

            default:
            {
                lassert( false, "unexpected PATH child kind " + str(src.kind) );
                break;
            }
        }
    }

    return dst_first_i;
}

void Layout::node_transform_xy( uint parent_i, uint xy_i, COPY_KIND copy_kind, CONFLICT_POLICY conflict_policy, const Matrix& M )
{
    lassert( nodes[xy_i].kind == NODE_KIND::XY, "node_transform_xy is not called on an XY node" );      
    bool is_y = false;
    real3 v;
    v.c[2] = 1.0;
    uint i = 0;
    uint prev_i = NULL_I;
    AABR brect;
    for( uint child_i = nodes[xy_i].u.child_first_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
    {
        //---------------------------------------------------------
        // Transform next X,Y pair
        //---------------------------------------------------------
        v.c[i&1] = real(nodes[child_i].u.i) * gdsii_units_user;
        if ( i&1 ) {
            real3 r;
            M.transform( v, r, true ); // divide by w
            if ( copy_kind == COPY_KIND::FLATTEN ) {
                ldout << "M=\n" << M << "\n";
                ldout << "v=" << v << "\n";
                ldout << "r=" << r << "\n";
            }

            real2 p( r.c[0] / gdsii_units_user,    // new X
                     r.c[1] / gdsii_units_user );  // new Y
            nodes[prev_i].u.i  = p.c[0];
            nodes[child_i].u.i = p.c[1];

            if ( copy_kind == COPY_KIND::FLATTEN ) {
                real2 r2( r.c[0], r.c[1] );
                if ( i == 1 ) {
                    brect = AABR( r2 );
                } else {
                    brect.expand( r2 );    
                }
                ldout << "brect=" << brect << "\n";     
            }
        }
        i++;
        prev_i = child_i;
    }

    if ( copy_kind == COPY_KIND::FLATTEN ) {
        //-----------------------------------------------------
        // The parent should be an expanded element.
        // Add it to leaf_nodes[] which has the computed brect.
        // 
        // Add the element (parent) to the BAH.
        //-----------------------------------------------------
        lassert( parent_i != uint(-1) && node_is_element( nodes[parent_i] ), "XY parent should be an element" );
        uint bah_layer_i = node_bah_layer( nodes[parent_i] );

        perhaps_realloc( leaf_nodes, hdr->leaf_node_cnt, max->leaf_node_cnt, 1 );
        uint leaf_i = hdr->leaf_node_cnt++;
        leaf_nodes[leaf_i].node_i = parent_i;
        leaf_nodes[leaf_i].brect  = brect;

        bah_add( bah_layer_i, leaf_i, conflict_policy );
    }
}

uint Layout::bah_node_alloc( void )
{
    perhaps_realloc( bah_nodes, hdr->bah_node_cnt, max->bah_node_cnt, 1 );
    uint bi = hdr->bah_node_cnt++;
    for( uint i = 0; i < 2; i++ )
    {
        for( uint j = 0; j < 2; j++ )
        {
            bah_nodes[bi].child_i[i][j]       = NULL_I;
            bah_nodes[bi].child_is_leaf[i][j] = true;
        }
    }
    return bi;
}

void Layout::bah_add( uint bah_layer_i, uint leaf_i, CONFLICT_POLICY conflict_policy )
{
    if ( layers[bah_layer_i].bah_root_i == NULL_I ) {
        //------------------------------------------------------------
        // This is the first leaf node.
        // Create the first BAH node with the leaf in quadrant 0,0.
        // But make sure the bounding rectangle is actually square in shape and
        // scaled up by 2x to make sense with the child in quadrant 0,0.
        //------------------------------------------------------------
        uint bi = bah_node_alloc();
        bah_nodes[bi].child_i[0][0] = leaf_i;

        layers[bah_layer_i].bah_root_i = bi;
        layers[bah_layer_i].bah_brect  = leaf_nodes[leaf_i].brect.enclosing_square( 2.0 );

    } else {
        const AABR& leaf_brect = leaf_nodes[leaf_i].brect;
        while( !layers[bah_layer_i].bah_brect.encloses( leaf_brect ) )
        {
            //------------------------------------------------------------
            // The entire BAH needs to be expanded because it won't contain
            // the new leaf node.  Create a new BAH_Node above the current root
            // and put the current root in the quadrant that allows for expansion 
            // in the right direction.
            //------------------------------------------------------------
            AABR bah_brect = layers[bah_layer_i].bah_brect;
            real w = bah_brect.max.c[0] - bah_brect.min.c[0];
            real h = bah_brect.max.c[1] - bah_brect.min.c[1];
            uint i, j;
            if ( bah_brect.min.c[0] < leaf_brect.min.c[0] ) {
                j = 0;
                bah_brect.max.c[0] += w;
            } else {
                j = 1;
                bah_brect.min.c[0] -= w;
            }
            if ( bah_brect.min.c[1] < leaf_brect.min.c[1] ) {
                i = 0;
                bah_brect.max.c[1] += h;
            } else {
                i = 1;
                bah_brect.min.c[1] -= h;
            }
            
            uint bi = bah_node_alloc();
            bah_nodes[bi].child_i[i][j]       = layers[bah_layer_i].bah_root_i;
            bah_nodes[bi].child_is_leaf[i][j] = false;

            layers[bah_layer_i].bah_root_i = bi;
            layers[bah_layer_i].bah_brect  = bah_brect;
        }

        //------------------------------------------------------------
        // Insert leaf recursively starting at the root.
        //------------------------------------------------------------
        bah_insert( layers[bah_layer_i].bah_root_i, layers[bah_layer_i].bah_brect, leaf_i, leaf_nodes[leaf_i].brect, conflict_policy );
    }
}

void Layout::bah_insert( uint bi, const AABR& bah_brect, uint li, const AABR& leaf_brect, CONFLICT_POLICY conflict_policy, std::string indent_str )
{
    //------------------------------------------------------------
    // Insert into each quadrant that intersects with the new leaf.
    //------------------------------------------------------------
    std::cout << indent_str << "bah_insert bi=" << bi << " bah_brect=" << bah_brect << " li=" << li << " leaf_brect=" << leaf_brect << "\n";
    indent_str += "  ";
    std::string indent_str2 = indent_str + "  ";
    for( uint i = 0; i < 2; i++ )
    {
        for( uint j = 0; j < 2; j++ )
        {
            AABR quadrant_brect = bah_brect.quadrant( i, j );
            if ( quadrant_brect.intersects( leaf_brect ) ) {
                std::cout << indent_str << "quadrant[" << i << "," << j << "]=" << quadrant_brect << " intersects with new leaf\n";
                if ( bah_nodes[bi].child_i[i][j] == NULL_I ) {
                    //------------------------------------------------------------
                    // No child.  Make this leaf the child.
                    //------------------------------------------------------------
                    std::cout << indent_str2 << "empty child, make this leaf the child\n";
                    bah_nodes[bi].child_i[i][j] = li;
                    bah_nodes[bi].child_is_leaf[i][j] = true;

                } else if ( !bah_nodes[bi].child_is_leaf[i][j] ) {
                    //------------------------------------------------------------
                    // Call this recursively on the child BAH node.
                    // But first calculate the child's bounding rectangle based on i,j.
                    //------------------------------------------------------------
                    std::cout << indent_str2 << "child is not a leaf, call recursively...\n";
                    AABR quadrant_leaf_brect = leaf_brect.quadrant( i, j );
                    bah_insert( bah_nodes[bi].child_i[i][j], quadrant_brect, li, quadrant_leaf_brect, conflict_policy, indent_str2 );

                } else {
                    //------------------------------------------------------------
                    // Child is a leaf.
                    // First, check for overlap with the leaf that is already there.
                    //------------------------------------------------------------
                    uint li2 = bah_nodes[bi].child_i[i][j];
                    bool is_exact;
                    if ( bah_leaf_nodes_intersect( quadrant_brect, li, li2, is_exact ) ) {
                        std::cout << indent_str2 << "child is a leaf, and intersects with new leaf\n";
                        if ( is_exact ) {
                            //------------------------------------------------------------
                            // Exact duplicates are not added to the BAH.
                            // If they are allowed, discard this new node.
                            //------------------------------------------------------------
                            lassert( conflict_policy != CONFLICT_POLICY::MERGE_NONE_ALLOW_NONE, 
                                     "overlapping elements are not allowed by confict_policy MERGE_NONE_ALLOW_NONE" );
                        } else {
                            //------------------------------------------------------------
                            // Partial intersections are allowed only with MERGE_ALL
                            // and the merge was already done by bah_leaf_nodes_intersect().
                            //------------------------------------------------------------
                            lassert( conflict_policy == CONFLICT_POLICY::MERGE_ALL,
                                     "merging partially overlapping elements is allowed only by conflict_policy MERGE_ALL" );
                        }
                    } else {
                        //------------------------------------------------------------
                        // We have to add a BAH_Node as the new child, then
                        // recursively add both leaves to it.  At some point, they won't conflict.
                        //------------------------------------------------------------
                        std::cout << indent_str2 << "child is a non-intersecting leaf, replace with new BAH node and recurse for both leaves...\n";
                        uint cbi = bah_node_alloc();
                        bah_nodes[bi].child_i[i][j] = cbi;
                        bah_nodes[bi].child_is_leaf[i][j] = false;
                        for( uint ii = 0; ii < 2; ii++ )
                        {
                            for( uint jj = 0; jj < 2; jj++ )
                            {
                                bah_nodes[cbi].child_i[ii][jj]     = NULL_I;
                                bah_nodes[cbi].child_is_leaf[i][j] = true;
                            }
                        }
                        AABR quadrant_leaf_brect  = leaf_brect.quadrant( i, j );
                        AABR quadrant_leaf2_brect = leaf_nodes[li2].brect.quadrant( i, j );
                        bah_insert( cbi, quadrant_brect, li,  quadrant_leaf_brect,  conflict_policy, indent_str2 );
                        bah_insert( cbi, quadrant_brect, li2, quadrant_leaf2_brect, conflict_policy, indent_str2 );
                    }
                }
            }
        }
    }
}

bool Layout::bah_leaf_nodes_intersect( const AABR& quadrant_brect, uint li1, uint li2, bool& is_exact ) 
{
    //------------------------------------------------------------
    // First, clamp the leaf brects by the current quadrant_brect.      
    // Return false if resulting leaf brects do not intersect.
    //------------------------------------------------------------
    Leaf_Node& leaf1 = leaf_nodes[li1];
    Leaf_Node& leaf2 = leaf_nodes[li2];
    AABR leaf1_brect = leaf1.brect;
    AABR leaf2_brect = leaf2.brect;
    leaf1_brect.intersect( quadrant_brect );     
    leaf2_brect.intersect( quadrant_brect );     
    if ( !leaf1_brect.intersects( leaf2_brect ) ) return false;

    //------------------------------------------------------------
    // See if there is an exact intersection by comparing XY nodes.
    //------------------------------------------------------------
    const Node& node1 = nodes[leaf1.node_i];
    const Node& node2 = nodes[leaf2.node_i];
    lassert( node_is_element( node1 ), "leaf1 is not an element" );    
    lassert( node_is_element( node2 ), "leaf2 is not an element" );    
    lassert( node1.kind == node2.kind, "ERROR: leaf1 and leaf2 must have same kind if bounding rects overlap, got " +
                                       str(node1.kind) + " and " + str(node2.kind) )  ;
    uint xy1_i = node_xy_i( node1 );
    uint xy2_i = node_xy_i( node2 );
    lassert( xy1_i != NULL_I, "leaf1 has no XY node" );
    lassert( xy2_i != NULL_I, "leaf2 has no XY node" );
    const Node& xy1 = nodes[xy1_i];
    Node& xy2 = nodes[xy2_i];

    is_exact = true;  // for now
    uint c1_i, c2_i;
    for( c1_i = xy1.u.child_first_i, c2_i = xy2.u.child_first_i;
         c1_i != NULL_I || c2_i != NULL_I;
         c1_i = nodes[c1_i].sibling_i, c2_i = nodes[c2_i].sibling_i )
    {
        lassert( c1_i == NULL_I || nodes[c1_i].kind == NODE_KIND::INT, "leaf1 XY coord is not an INT" );
        lassert( c2_i == NULL_I || nodes[c2_i].kind == NODE_KIND::INT, "leaf2 XY coord is not an INT" );
        if ( c1_i == NULL_I || c2_i == NULL_I || nodes[c1_i].u.i != nodes[c2_i].u.i ) {
            is_exact = false;
            break;
        }
    }
    if ( is_exact ) return true;

    //------------------------------------------------------------
    // Complicated Case
    //
    // We probably have a partial intersection, but there's
    // a chance there is no intersection at all or one polygon completely covers the other.   
    // Attempt to intersect and merge them.
    //
    // First, make this easier by copying the vertexes from
    // both polygons into arrays that we can process more easily.
    // Also create a square from the bah_brect.
    //------------------------------------------------------------
    uint vtx1_cnt;
    uint vtx2_cnt;
    uint vtxr_cnt;
    real2 * vtx1 = polygon_alloc( xy1,            vtx1_cnt );
    real2 * vtx2 = polygon_alloc( xy2,            vtx2_cnt );
    real2 * vtxr = polygon_alloc( quadrant_brect, vtxr_cnt );
    ldout << "xy1, xy2, qbrect:\n";
    ldout << polygon_str( vtx1, vtx1_cnt, "red" ) << "\n";
    ldout << polygon_str( vtx2, vtx2_cnt, "green" ) << "\n";
    ldout << polygon_str( vtxr, vtxr_cnt, "yellow" ) << "\n";

    //------------------------------------------------------------
    // Intersect both leaf polygons with the quadrant brect polygon.
    //------------------------------------------------------------
    uint vtx1r_cnt;
    uint vtx2r_cnt;
    real2 * vtx1r = polygon_merge_or_intersect( false, vtx1, vtx1_cnt, vtxr, vtxr_cnt, vtx1r_cnt );
    real2 * vtx2r = polygon_merge_or_intersect( false, vtx2, vtx2_cnt, vtxr, vtxr_cnt, vtx2r_cnt );
    ldout << "xy1r, xy2r:\n";
    ldout << polygon_str( vtx1r, vtx1r_cnt, "red" ) << "\n";
    ldout << polygon_str( vtx2r, vtx2r_cnt, "green" ) << "\n";

    //------------------------------------------------------------
    // Merge the resultant polygons to get the final merged polygon.  
    // If nullptr is returned, that means that the two didn't really intersect
    // at all in this quadrant and nothing should be changed.
    //------------------------------------------------------------
    uint vtx_cnt;
    real2 * vtx = polygon_merge_or_intersect( true, vtx1r, vtx1r_cnt, vtx2r, vtx2r_cnt, vtx_cnt );
    if ( vtx != nullptr ) {
        //------------------------------------------------------------
        // Replace the 2nd leaf node's polygon with the merged polygon.
        //------------------------------------------------------------
        node_xy_replace_polygon( xy2, vtx, vtx_cnt );
    }

    //------------------------------------------------------------
    // Deallocate vtx arrays.
    //------------------------------------------------------------
    if ( vtx != nullptr ) polygon_dealloc( vtx );
    polygon_dealloc( vtx1r );
    polygon_dealloc( vtx2r );
    polygon_dealloc( vtx1 );
    polygon_dealloc( vtx2 );
    polygon_dealloc( vtxr );
    return vtx != nullptr; 
}

Layout::real2 * Layout::polygon_alloc( uint vtx_cnt ) const
{
    //------------------------------------------------------------
    // Later, we'll cache these.
    //------------------------------------------------------------
    return new real2[vtx_cnt];
}

Layout::real2 * Layout::polygon_alloc( const Node& xy_node, uint& vtx_cnt ) const
{
    //------------------------------------------------------------
    // First count number of vertices.
    //------------------------------------------------------------
    lassert( xy_node.kind == NODE_KIND::XY, "expected XY, got " + str(xy_node.kind) );
    vtx_cnt = 0;
    bool have_x = false;                
    for( uint ci = xy_node.u.child_first_i; ci != NULL_I; ci = nodes[ci].sibling_i )
    {
        if ( have_x ) vtx_cnt++;
        have_x = !have_x;
    }
    lassert( !have_x, "no Y coord" );

    //------------------------------------------------------------
    // Now allocate that vtx_array and copy them out.
    //------------------------------------------------------------
    real2 * vtx_array = polygon_alloc( vtx_cnt );
    uint i = 0;
    for( uint ci = xy_node.u.child_first_i; ci != NULL_I; ci = nodes[ci].sibling_i )
    {
        if ( !have_x ) {
            vtx_array[i].c[0] = real(nodes[ci].u.i) * gdsii_units_user;
        } else {
            vtx_array[i++].c[1] = real(nodes[ci].u.i) * gdsii_units_user;
        }
        have_x = !have_x;
    }
    return vtx_array;
}

Layout::real2 * Layout::polygon_alloc( const AABR& brect, uint& vtx_cnt ) const
{
    vtx_cnt = 5;        
    real2 * vtx_array = new real2[vtx_cnt];
    vtx_array[0] = brect.min;
    vtx_array[1].c[0] = brect.max.c[0]; 
    vtx_array[1].c[1] = brect.min.c[1]; 
    vtx_array[2] = brect.max;
    vtx_array[3].c[0] = brect.min.c[0]; 
    vtx_array[3].c[1] = brect.max.c[1]; 
    vtx_array[4] = vtx_array[0];
    return vtx_array;
}

Layout::real2 * Layout::polygon_copy( const real2 * other, uint other_cnt ) const
{
    real2 * vtx = polygon_alloc( other_cnt );
    memcpy( vtx, other, other_cnt * sizeof( real2 ) );
    return vtx;
}

Layout::real2 * Layout::polygon_reversed( const real2 * other, uint other_cnt ) const
{
    real2 * vtx = polygon_alloc( other_cnt );
    for( uint i = 0; i < other_cnt; i++ )
    {
        vtx[i] = other[other_cnt-1-i];
    }
    return vtx;
}

void Layout::polygon_dealloc( real2 * vtx_array ) const
{
    delete[] vtx_array;
}

Layout::real2 * Layout::polygon_merge_or_intersect( bool do_merge, const real2 * _vtx1, uint vtx1_cnt, const real2 * _vtx2, uint vtx2_cnt, 
                                                    uint& vtx_cnt ) const
{
    //------------------------------------------------------------
    // To simplify things, make sure poly1 and poly2 are in counterclockwise winding order.
    //------------------------------------------------------------
    real2 * vtx1 = polygon_ccw( _vtx1, vtx1_cnt );
    real2 * vtx2 = polygon_ccw( _vtx2, vtx2_cnt );
    ldout << "\npolygon_merge_or_intersect: do_merge=" << do_merge << " poly1_was_ccw=" << polygon_is_ccw( _vtx1, vtx1_cnt ) << " ccw_poly1=" << polygon_str( vtx1, vtx1_cnt, "red" ) <<
                                                                      " poly2_was_ccw=" << polygon_is_ccw( _vtx2, vtx2_cnt ) << " ccw_poly2=" << polygon_str( vtx2, vtx2_cnt, "green" ) << "\n";

    //------------------------------------------------------------
    // First find any intersection point between the two polygons, 
    // excluding segment endpoints.
    // TODO: not correct, can still have polygons that overlap and meet only at common endpoints
    //------------------------------------------------------------
    uint i, i2, j, j2;
    real2 ip;
    j = vtx2_cnt;
    for( i = 0; i < vtx1_cnt; i++ )
    {
        i2 = (i + 1) % vtx1_cnt;
        for( j = 0; j < vtx2_cnt; j++ )
        {
            j2 = (j + 1) % vtx2_cnt;
            bool colinear_unused;
            bool reverse_unused;
            if ( vtx1[i].segments_intersection( vtx1[i2], vtx2[j], vtx2[j2], ip, false, false, false, colinear_unused, reverse_unused ) ) break;
        }
        if ( j != vtx2_cnt ) break;
    }   

    real2 * vtx;
    if ( i == vtx1_cnt ) {
        ldout << "no proper intersection, checking for complete containment...\n";

        //------------------------------------------------------------
        // NO INTERSECTION
        //
        // Now, we need to check to see if one polygon
        // completely contains the other.  We can use the bounding boxes 
        // to determine this.  One will contain the other.
        //------------------------------------------------------------
        lassert( j == vtx2_cnt, "j should be vtx2_cnt at this point" );
        AABR brect1 = polygon_brect( vtx1, vtx1_cnt ); 
        AABR brect2 = polygon_brect( vtx2, vtx2_cnt ); 
        bool brect1_encloses = brect1.encloses( brect2 );
        bool brect2_encloses = !brect1_encloses && brect2.encloses( brect1 );
        if ( brect1_encloses || brect2_encloses ) {
            if ( brect1_encloses ) {
                if ( do_merge ) {
                    vtx_cnt = vtx1_cnt;
                    vtx = polygon_copy( vtx1, vtx1_cnt );
                    ldout << "using polygon1 for merge\n";
                } else {
                    vtx_cnt = vtx2_cnt;
                    vtx = polygon_copy( vtx2, vtx2_cnt );
                    ldout << "using polygon2 for intersection\n";
                }
            } else {
                if ( do_merge ) {
                    vtx_cnt = vtx2_cnt;
                    vtx = polygon_copy( vtx2, vtx2_cnt );
                    ldout << "using polygon2 for merge\n";
                } else {
                    vtx_cnt = vtx1_cnt;
                    vtx = polygon_copy( vtx1, vtx1_cnt );
                    ldout << "using polygon1 for intersection\n";
                }
            }
        } else {
            // they don't overlap at all 
            ldout << "neither contains the other\n";
            vtx = nullptr;  
            vtx_cnt = 0;
        }

        ldout << "polygon_merge_or_intersect: result=" << polygon_str( vtx, vtx_cnt, "blue" ) << "\n";
        polygon_dealloc( vtx1 );
        polygon_dealloc( vtx2 );
        return vtx;
    }

    //------------------------------------------------------------
    // INTERSECTION
    //
    // Allocate a vtx array for the resultant polygon.
    // Start with the intersection point from above.
    //------------------------------------------------------------
    uint vtx_cnt_max = 2*(vtx1_cnt + vtx2_cnt);
    vtx = new real2[vtx_cnt_max];
    vtx_cnt = 0;                

    uint          curr           = 0;                       // current polygon index
    uint          other          = 1;                       // other   polygon index
    const real2 * vtxn[2]        = { vtx1, vtx2 };
    const uint    vtxn_cnt[2]    = { vtx1_cnt, vtx2_cnt };  
    bool          have_ip        = true;
    bool          have_colinear_ip = false;
    uint          curr_s0_i      = i;
    uint          curr_s1_i      = i2;
    uint          other_s0_i     = j;
    uint          other_s1_i     = j2;
    bool          is_ccw     = true;
    for( ;; ) 
    {
        lassert( vtx_cnt != vtx_cnt_max, "vtx array grew bigger than expected" );
        if ( have_colinear_ip ) {
            //------------------------------------------------------------
            // Record the new intersection point, ip.
            //------------------------------------------------------------
            ldout << "\nsave colinear intersection point as next vertex in result: " << ip.str() << " is_ccw=" << is_ccw << "\n\n";
            vtx[vtx_cnt++] = ip;
            ldout << "partial result=" << polygon_str( vtx, vtx_cnt, "blue" );

            //------------------------------------------------------------
            // If we're back at the first intersection point, we are done.
            // However, make sure the last vertex is an exact copy of the first.
            //------------------------------------------------------------
            if ( vtx_cnt != 1 && vtx[vtx_cnt-1].nearly_equal( vtx[0] ) ) { 
                vtx[vtx_cnt-1] = vtx[0];  // exact copy to ensure polygon is water-tight
                ldout << "back at first intersection point, so done, vtx[0]=" << vtx[0] << " vtx[last]=" << vtx[vtx_cnt-1] << "\n";
                break;
            }

            curr_s0_i  = other_s0_i;
            curr_s1_i  = other_s1_i;
            curr       = other;
            other      = 1 - curr;

        } else if ( have_ip ) {
            const real2& curr_s0  = vtxn[curr][curr_s0_i];
            const real2& curr_s1  = vtxn[curr][curr_s1_i];
            const real2& other_s0 = vtxn[other][other_s0_i];
            const real2& other_s1 = vtxn[other][other_s1_i];
            ldout << "\n";
            ldout << "curr_seg: [ " << curr_s0.str() << ", " << ip.str() << ", " << curr_s1.str() << " ] curr=" << curr << " curr_s0_i=" << curr_s0_i << " curr_s1_i=" << curr_s1_i << " is_ccw=" << is_ccw << "\n";
            ldout << "other_seg:[ " << other_s0.str() << ", " << other_s1.str() << " ] other=" << other << " other_s0_i=" << other_s0_i << " other_s1_i=" << other_s1_i << "\n";

            //------------------------------------------------------------
            // Record the new intersection point, ip.
            //------------------------------------------------------------
            ldout << "save intersection point as next vertex in result: " << ip.str() << " is_ccw=" << is_ccw << "\n\n";
            vtx[vtx_cnt++] = ip;
            ldout << "partial result=" << polygon_str( vtx, vtx_cnt, "blue" );

            //------------------------------------------------------------
            // If we're back at the first intersection point, we are done.
            // However, make sure the last vertex is an exact copy of the first.
            //------------------------------------------------------------
            if ( vtx_cnt != 1 && vtx[vtx_cnt-1].nearly_equal( vtx[0] ) ) { 
                vtx[vtx_cnt-1] = vtx[0];  // exact copy to ensure polygon is water-tight
                ldout << "back at first intersection point, so done, vtx[0]=" << vtx[0] << " vtx[last]=" << vtx[vtx_cnt-1] << "\n";
                break;
            }

            lassert( !polygon_includes( vtx, vtx_cnt-1, ip ), "intersection point should not already be on the vtx[] list" );

            //------------------------------------------------------------
            // For merge,        head outside the current polygon.
            // For intersection, head inside  the current polygon.
            //------------------------------------------------------------
            bool other_s0_is_left  = other_s0.is_left_of_segment(  curr_s0, curr_s1 );
            bool other_s0_is_right = other_s0.is_right_of_segment( curr_s0, curr_s1 );
            bool other_s1_is_left  = other_s1.is_left_of_segment(  curr_s0, curr_s1 );
            bool other_s1_is_right = other_s1.is_right_of_segment( curr_s0, curr_s1 );
            bool use_other_s0_for_s0;
            bool use_other_s1_for_s0;
            if ( is_ccw ) {
                use_other_s0_for_s0 = do_merge ? other_s0_is_left : other_s1_is_left;
                use_other_s1_for_s0 = do_merge ? other_s1_is_left : other_s0_is_left;
            } else {
                use_other_s0_for_s0 = do_merge ? other_s0_is_right : other_s1_is_right;
                use_other_s1_for_s0 = do_merge ? other_s1_is_right : other_s0_is_right;
            }
            ldout << " other_s0_is_left=" << other_s0_is_left << " other_s0_is_right=" << other_s0_is_right <<
                     " other_s1_is_left=" << other_s1_is_left << " other_s1_is_right=" << other_s1_is_right <<
                     " is_ccw=" << is_ccw << " do_merge=" << do_merge << " use_other_s0_for_s0=" << use_other_s0_for_s0 << " use_other_s1_for_s0=" << use_other_s1_for_s0 << "\n"; 
            lassert( use_other_s0_for_s0 != use_other_s1_for_s0, "use_other_s0_for_s0 and use_other_s1_for_s0 can't both be set" );

            other_s0_i  = use_other_s0_for_s0 ? other_s0_i : other_s1_i;
            other_s1_i  = use_other_s1_for_s0 ? other_s0_i : other_s1_i;

            //------------------------------------------------------------
            // There's an ugly case we need to check when ip == curr_s1.
            // If segment [ip, other_s1] is already part of of an existing
            // segment in the result polygon, then we need to 
            // pretend that we don't have an intersection point.
            //------------------------------------------------------------
            if ( ip.nearly_equal( curr_s1 ) ) {
                ldout << " checking to see if ip=curr_s1=" << ip << ", other_s1=" << other_s1 << " is already part of an existing result segment...\n";
                for( uint k = 0; k < (vtxn_cnt[other]-1); k++ )
                {
                    bool ip_on_seg          = ip.is_on_segment( vtxn[other][k], vtxn[other][k+1], true );
                    bool other_s1_on_seg    = other_s1.is_on_segment( vtxn[other][k], vtxn[other][k+1], true );
                    bool seg_s0_on_ip_other = vtxn[other][k].is_on_segment( ip, other_s1, true );
                    bool seg_s1_on_ip_other = vtxn[other][k+1].is_on_segment( ip, other_s1, true );
                    ldout << "\n  k=" << k << " seg_s0=" << vtxn[other][k] << " seg_s1=" << vtxn[other][k+1] << 
                             " ip_on_seg=" << ip_on_seg << " other_s1_on_seg=" << other_s1_on_seg << 
                             " seg_s0_on_ip_other=" << seg_s0_on_ip_other << " seg_s1_on_ip_other=" << seg_s1_on_ip_other << "\n";
                    if ( ip_on_seg && other_s1_on_seg ) {
                        have_ip = false;
                        break;
                    }
                }
            }

            if ( have_ip ) {
                //------------------------------------------------------------
                // Set curr_s0_i = chosen other s0
                // Set curr_s1_i = chosen other s1
                // Switch to the other polygon.
                // Follow (ip, curr_s1).
                //------------------------------------------------------------
                curr_s0_i  = other_s0_i;
                curr_s1_i  = other_s1_i;
                is_ccw     = use_other_s0_for_s0;
                curr       = other;
                other      = 1 - curr;
            } else {
                //------------------------------------------------------------
                // Backed out.
                //------------------------------------------------------------
                ip = curr_s1;
                vtx[vtx_cnt-1] = ip;
                ldout << "new curr_seg was already visited, switching back to using previous endpoint: " << ip.str() << " is_ccw=" << is_ccw << "\n\n"; 
                lassert( !polygon_includes( vtx, vtx_cnt, ip ), "intersection point should not already be on the vtx[] list" );
                curr_s0_i = curr_s1_i;
                curr_s1_i = (is_ccw ? (curr_s0_i+1) : (curr_s0_i+vtxn_cnt[curr]-2)) % (vtxn_cnt[curr]-1);
            }

        } else {
            //------------------------------------------------------------
            // Record the endpoint (curr_s1_i) of the current segment.
            // Then make that the curr_s0_i and calculate the next curr_s1_i 
            // based on our current direction.
            //------------------------------------------------------------
            ip = vtxn[curr][curr_s1_i];
            ldout << "save endpoint as next vertex in result: " << ip.str() << " is_ccw=" << is_ccw << "\n\n";
            lassert( !polygon_includes( vtx, vtx_cnt, ip ), "intersection point should not already be on the vtx[] list" );
            vtx[vtx_cnt++] = ip;

            curr_s0_i = curr_s1_i;
            curr_s1_i = (is_ccw ? (curr_s0_i+1) : (curr_s0_i+vtxn_cnt[curr]-2)) % (vtxn_cnt[curr]-1);
        }
        ldout << "new curr_seg: [ " << vtxn[curr][curr_s0_i] << ", " << ip << ", " << vtxn[curr][curr_s1_i] << " ] curr=" << curr << 
                 " curr_s0_i=" << curr_s0_i << " curr_s1_i=" << curr_s1_i << " is_ccw=" << is_ccw << "\n";
        lassert( !ip.nearly_equal( vtxn[curr][curr_s1_i] ), "curr_seg is a point" );

        //------------------------------------------------------------
        // See if there's an intersection point along (ip, curr_s1) with
        // the (new) other segment.  Choose the one that is closest to ip.
        //------------------------------------------------------------
        ldout << "checking for intersection between new curr_seg and some other_seg...\n";
        uint  best_k = NULL_I;
        uint  best_k2 = NULL_I;
        bool  best_colinear = false;
        bool  best_reverse = false;
        real  best_dist = 1e100;
        real2 best_ip;
        for( uint k = 0; k < (vtxn_cnt[other]-1); k++ )
        {
            uint k2 = k + 1;
            real2 this_ip;
            const real2 curr_s1 = vtxn[curr][curr_s1_i];
            bool colinear = false;
            bool reverse = false;
            if ( ip.segments_intersection( curr_s1, vtxn[other][k], vtxn[other][k2], this_ip, true, false, true, colinear, reverse ) ) {
                if ( !this_ip.nearly_equal( ip ) ) {
                    real this_ip_dist = (this_ip - ip).length();     
                    if ( best_k == NULL_I || this_ip_dist < best_dist ) {
                        best_k       = k;
                        best_k2      = k2;
                        best_colinear = colinear;
                        best_reverse = reverse;
                        best_dist    = this_ip_dist;
                        best_ip      = this_ip;
                        ldout << " ip=" << ip << " this_ip=" << this_ip << " best_ip=" << best_ip << " best_dist=" << best_dist << 
                                 " best_colinear=" << best_colinear << " best_reverse=" << best_reverse << "\n";
                    }
                } else {
                        ldout << " REJECTED\n";
                }
            }
        }

        have_ip = best_k != NULL_I;
        if ( have_ip ) {
            //------------------------------------------------------------
            // There was an intersection.
            // Set the new ip to best_ip and go back to the top of the loop.
            //------------------------------------------------------------
            ip = best_ip;
            ldout << " closest intersection point along new curr_seg: " << ip << " colinear=" << best_colinear << " reverse=" << best_reverse << "\n";
            have_colinear_ip = best_colinear;
            other_s0_i = best_reverse ? best_k2 : best_k;
            other_s1_i = best_reverse ? best_k  : best_k2;
        } else {
            have_colinear_ip = false;
            ldout << " no intersection point before end of new curr_seg\n";
        }
    }    

    ldout << "polygon_merge_or_intersect: result=" << polygon_str( vtx, vtx_cnt, "blue" ) << "\n";
    polygon_dealloc( vtx1 );
    polygon_dealloc( vtx2 );
    lassert( vtx_cnt == 0 || vtx_cnt >= 4, "polygon_merge_or_intersect: not a polygon vtx_cnt=" + std::to_string(vtx_cnt) );
    return vtx;
}

Layout::AABR Layout::polygon_brect( const real2 * vtx, uint vtx_cnt ) const
{
    AABR brect( vtx[0] );       
    for( uint i = 1; i < vtx_cnt; i++ )
    {
        brect.expand( vtx[i] );
    }
    return brect;      
}

std::string Layout::polygon_str( const real2 * vtx, uint vtx_cnt, std::string color, real xy_scale, real x_off, real y_off ) const
{
    std::string s = "<polygon points='";
    for( uint i = 0; i < vtx_cnt; i++ )
    {
        if ( i != 0 ) s += " ";
        s += ::str( vtx[i].c[0] );
        s += ",";
        s += ::str( vtx[i].c[1] );
    }
    s += "' transform='translate(" + ::str(x_off) + " " + ::str(2000-y_off) + 
          ") scale(" + ::str(xy_scale) + " " + ::str(-xy_scale) + 
          ")' style=\"fill:" + color + ";fill-opacity:0.5;stroke:black;stroke-width:0.1\"/>";
    return s;
}

bool Layout::polygon_eq( const real2 * vtx1, uint vtx1_cnt, const real2 * vtx2, uint vtx2_cnt ) const
{
    if ( vtx1_cnt != vtx2_cnt ) return false;
    for( uint i = 0; i < vtx1_cnt; i++ )
    {
        if ( vtx1[i] != vtx2[i] ) return false;
    }
    return true;
}

bool Layout::polygon_encloses( const real2 * vtx, uint vtx_cnt, const real2& p ) const
{
    AABR brect = polygon_brect( vtx, vtx_cnt );
    return brect.encloses( p );
}

bool Layout::polygon_includes( const real2 * vtx, uint vtx_cnt, const real2& v ) const
{
    for( uint i = 0; i < vtx_cnt; i++ )
    {
        if ( v.nearly_equal( vtx[i] ) ) return true;
    }
    return false;
}

bool Layout::polygon_is_ccw( const real2 * vtx, uint vtx_cnt ) const
{
    real sum = 0.0;
    for( uint i = 0; i < (vtx_cnt-1); i++ ) 
    {
        const real2& v1 = vtx[i];
        const real2& v2 = vtx[(i+1) % vtx_cnt];
        sum += (v2.c[0] - v1.c[0]) * (v2.c[1] + v1.c[1]);
    }
    return sum < 0.0;
}

Layout::real2 * Layout::polygon_ccw( const real2 * vtx, uint vtx_cnt ) const
{
    return polygon_is_ccw( vtx, vtx_cnt ) ? polygon_copy( vtx, vtx_cnt ) : polygon_reversed( vtx, vtx_cnt );
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
    _uread( leaf_nodes,  Leaf_Node,   hdr->leaf_node_cnt );
    _uread( bah_nodes,   BAH_Node,    hdr->bah_node_cnt );

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
    _uwrite( leaf_nodes,  hdr->leaf_node_cnt * sizeof(leaf_nodes[0] ));
    _uwrite( bah_nodes,   hdr->bah_node_cnt  * sizeof(bah_nodes[0] ));

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

            if ( !count_only ) nodes[ni].u.u = *reinterpret_cast<uint16_t *>( nnn );
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
                    real exp  = (uuu[0] & 0x7f) - 64;
                    real mantissa = 0;
                    for( uint j = 0; j < datum_byte_cnt; j++ )
                    {
                        ldout << "bytes[" << j << "]=" << int(uuu[j]) << "\n";
                        if ( j != 0 ) mantissa = mantissa*256.0 + real(uuu[j]);
                    }

                    if ( !count_only ) {
                        nodes[child_i].kind = NODE_KIND::REAL;
                        real rexp = 4.0*exp - 8*(datum_byte_cnt-1);
                        nodes[child_i].u.r = sign * mantissa * std::pow( 2.0, rexp );
                        ldout << "sign=" << ((sign < 0.0) ? "1" : "0") << " exp=" << exp << " rexp=" << rexp << " mantissa=" << mantissa << " r=" << nodes[child_i].u.r << "\n";
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
            lassert( false, "unexpected datatype" );
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
            NODE_KIND ckind = nodes[child_i].kind;
            if ( ckind == NODE_KIND::ENDEL || ckind == NODE_KIND::ENDSTR || ckind == NODE_KIND::ENDLIB ) break;
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
        lassert( false, "gdsii_write_record: unexpected node kind " + str(node.kind) );
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
                uint id2_i;
                if ( peek_id( id2_i, nnn, nnn_end ) ) {
                    if ( id2_i == aedt_end_str_i ) {
                        parse_id( id2_i, nnn, nnn_end );
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

bool Layout::aedt_write( std::string file, bool for_raw )
{
    std::ofstream out( file, std::ofstream::out );
    aedt_write_expr( out, hdr->root_i, "\n", for_raw );
    out.close();
    return true;
}

void Layout::aedt_write_expr( std::ofstream& out, uint ni, std::string indent_str, bool for_raw )
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
            aedt_write_expr( out, child_i, "", for_raw );
            out << "=";
            child_i = nodes[child_i].sibling_i;
            aedt_write_expr( out, child_i, "", for_raw );
            break;
        }

        case NODE_KIND::CALL:
        {
            uint child_i = node.u.child_first_i;
            aedt_write_expr( out, child_i, "", for_raw );
            out << "(";
            bool have_one = false;
            for( child_i = nodes[child_i].sibling_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
            {
                if ( !have_one ) out << " ";
                if ( have_one ) out << ", ";
                aedt_write_expr( out, child_i, "", for_raw );
                have_one = true;
            }
            if ( have_one ) out << " ";
            out << ")";
            break;
        }

        case NODE_KIND::SLICE:
        {
            uint child_i = node.u.child_first_i;
            aedt_write_expr( out, child_i, "", for_raw );
            out << "[";
            child_i = nodes[child_i].sibling_i;
            if ( child_i != NULL_I ) {
                aedt_write_expr( out, child_i, "", for_raw );
                out << ":";
                bool have_one = false;
                for( child_i = nodes[child_i].sibling_i; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
                {
                    out << (have_one ? ", " : " ");
                    aedt_write_expr( out, child_i, "", for_raw );
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
                if ( for_raw ) {
                    out << "BGNSTR" << std::string(&strings[nodes[id_i].u.s_i]);
                } else {
                    out << "$begin '" << std::string(&strings[nodes[id_i].u.s_i]) << "'";
                }
                child_i = nodes[id_i].sibling_i;
            } else {
                if ( for_raw ) {
                    out << "BGNFILE";
                } else {
                    out << "$begin 'FILE'";
                }
                child_i = id_i;
                id_i = NULL_I;
            }
            for( ; child_i != NULL_I; child_i = nodes[child_i].sibling_i )
            {
                aedt_write_expr( out, child_i, indent_str + (for_raw ? "    " : "\t"), for_raw );
            }
            out << indent_str;
            if ( id_i != NULL_I ) {
                if ( for_raw ) {
                    out << "ENDHIER " << std::string(&strings[nodes[id_i].u.s_i]); 
                } else {
                    out << "$end '" << std::string(&strings[nodes[id_i].u.s_i]) << "'";
                }
            } else {
                if ( for_raw ) {
                    out << "ENDFILE\n";
                } else {
                    out << "$end 'FILE'\n";
                }
            }
            break;
        }

        default:
        {
            if ( node_is_gdsii( node ) ) {
                GDSII_DATATYPE datatype = kind_to_datatype( node.kind );
                if ( node_is_hier( node ) ) {
                    if ( for_raw ) {
                        out << str(node.kind);
                    } else {
                        out << "$begin '" << node.kind << "'";
                    }
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
                        out << node.kind << "( " << node.u.u << " )";
                        break;
                    }

                    case GDSII_DATATYPE::STRING: {
                        out << node.kind << "( '" << std::string(&strings[node.u.s_i]) << "' )";
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
                            std::string this_indent_str = indent_str;
                            if ( for_raw ) {
                                if ( node.kind == NODE_KIND::BGNLIB || node.kind == NODE_KIND::BGNSTR ) {
                                    this_indent_str = "";
                                } else {
                                    this_indent_str = "    ";
                                }
                            }
                            out << this_indent_str << (for_raw ? "" : "ARGS") << ((vals != "") ? "( ": "(") << vals << ((vals != "") ? " )" : ")");
                        } else {
                            if ( vals != "" ) {
                                out << node.kind << "( " << vals << " )";
                            } else {
                                out << node.kind << "()";
                            }
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
                        aedt_write_expr( out, child_i, indent_str + (for_raw ? "    " : "\t"), for_raw );
                    }
                    out << indent_str;
                    if ( for_raw ) {
                        if ( node.kind != NODE_KIND::HIER ) {
                            out << hier_end_kind(node.kind);
                        } else {
                            out << "END" << node.kind;
                        }
                    } else {
                        out << "$end '" << node.kind << "'";
                    }
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
        if ( thickness == 0 ) thickness = gdsii_units_user;
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

bool Layout::hfss_write_layer_info( std::string file )
{
    std::ofstream out( file, std::ofstream::out );

    out << "UNITS um\n";
    out << "# L#    Name                    Color             Elevation    Thickness\n";
    out << "#-----------------------------------------------------------------------\n";
    real height = 0.0;
    for( uint i = 0; i < hdr->layer_cnt; i++ )
    {
        const Layer& layer = layers[i];
        real thickness = layer.thickness;
        std::string color = color_name( layer.material_rgba );
        
        out << putf( "%-4d    %-20s    %-10s   %10.3f   %10.3f\n", i, &strings[layer.name_i], color.c_str(), height, thickness );

        if ( i != (hdr->layer_cnt-1) && !layers[i+1].same_zoffset_as_prev ) height += thickness;
    }

    out.close();
    return true;
}

bool Layout::hfss_write_material_info( std::string file )
{
    std::ofstream out( file, std::ofstream::out );

    for( uint i = 0; i < hdr->material_cnt; i++ )
    {
        const Material& material = materials[i];
        const char * name = &strings[material.name_i];
        out << "$begin '" << name << "'\n";
        out << "\t$begin 'MaterialDef'\n";
        out << "\t\t$begin '" << name << "'\n";
        out << "\t\t\tCoordinateSystemType='Cartesian'\n";
	out << "\t\t\t$begin 'AttachedData'\n";
	out << "\t\t\t$end 'AttachedData'\n";
	out << "\t\t\t$begin 'ModifierData'\n";
	out << "\t\t\t$end 'ModifierData'\n";
        out << "\t\t\tpermittivity='" << material.relative_permittivity << "'\n";
        out << "\t\t\tpermeability='" << material.permeability << "'\n";
        out << "\t\t\tthermal_conductivity='" << material.thermal_conductivity << "'\n";
        out << "\t\t\tmass_density='" << material.mass_density << "'\n";
        out << "\t\t\tspecific_heat='" << material.specific_heat << "'\n";
        out << "\t\t$end '" << name << "'\n";
        out << "\t$end 'MaterialDef'\n";
        out << "$end '" << name << "'\n\n";
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

bool Layout::open_and_read( std::string file_name, uint8_t *& start, uint8_t *& end )
{
    const char * fname = file_name.c_str();
    int fd = open( fname, O_RDONLY );
    if ( fd < 0 ) ldout << "open_and_read() error reading " << file_name << ": " << strerror( errno ) << "\n";
    lassert( fd >= 0, "could not open file " + file_name + " - open() error: " + strerror( errno ) );

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

#endif
