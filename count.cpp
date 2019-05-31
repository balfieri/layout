// count.cpp - count number of input records only
//
#include "Layout.h"

void die( std::string msg )
{
    std::cout << "ERROR: " << msg << "\n";
    exit( 1 );
}

int main( int argc, const char * argv[] )
{
    if ( argc != 2 ) die( "usage: count <file>" );
    std::string file = argv[1];

    std::string file_ext;
    std::string file_base = Layout::path_without_ext( file, &file_ext );

    std::cout << "Reading " + file + "...\n";
    Layout * layout = new Layout( file, true );

    std::cout << std::flush;
    return 0; 
}
