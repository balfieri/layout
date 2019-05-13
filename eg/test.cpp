// test.cpp - test Layout class 
//
#include "Layout.h"

void die( std::string msg )
{
    std::cout << "ERROR: " << msg << "\n";
    exit( 1 );
}

int main( int argc, const char * argv[] )
{
    if ( argc != 2 ) die( "usage: test <file>" );
    std::string file = argv[1];

    std::string file_ext;
    std::string file_base = Layout::path_without_ext( file, &file_ext );

    if ( file_ext == ".gds" ) {
        std::cout << "Reading " + file + "...\n";
        GDSIIData * gds = new GDSIIData( file );
        if ( gds->ErrMsg ) die( "could not read " + file + ": " + *gds->ErrMsg );

    } else if ( file_ext == ".aedt" ) {
        std::cout << "Reading " + file + "...\n";
        Layout * layout = new Layout( file );
        if ( !layout->is_good ) die( "could not read " + file + ": " + layout->error_msg );

        std::string layout_file = file_base + ".layout";
        std::cout << "Writing " + layout_file + "...\n";
        layout->write( layout_file );
        delete layout;

        std::cout << "Reading " + layout_file + "...\n";
        layout = new Layout( layout_file );
        if ( !layout->is_good ) die( "could not read " + layout_file + ": " + layout->error_msg );

        file = file_base + ".exported.aedt";
        std::cout << "Writing " + file + "...\n";
        layout->write( file );
        delete layout;

    } else {
        die( "unknown file extension: " + file_ext );
    }

    std::cout << std::flush;
    return 0; 
}
