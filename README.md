Single Layout.h file provides:

<ul>
<li>reading and writing .gds files (GDSII industry standard layout format)</li>
<li>reading and writing .aedt files for Ansys HFSS 3D field solver</li> 
<li>reading and writing .layout files (our Layout.h binary representation which will soon be more compressed than .gds)</li>
<li>writing fastcap2/fastercap/fftcap input files for capacitance 3D field solver</li>
<li>writing fasthenry2 input files for resistance+inductance 3D field solver</li>
<li>instancing structures or instances from one or more source files</li>
<li>instancing all instances from one or more source files</li>
<li>filling empty 3D space (within an arbitrary box) with a dielectric material</li>
<li>flattening instances into uniquified structures</li>
<li>writing out a single file with resultant edits</li>
</ul>

<pre>
How to use it:

    0) Clone this repository

    1) #include "Layout.h"

       Layout * layout = new Layout( "my_chip.aedt" );
       layout->write( "my_chip.layout" );    // will write out the self-contained binary layout layout

    2) After that, you can quickly read in the single binary layout file using:

       Layout * layout = new Layout( "my_chip.layout" );  

    3) You can also write out (export) other types of files:
     
       layout->write( "new_chip.gds" );      // writes out a .gds II file
       layout->write( "new_chip.aedt" );     // writes out an .aedt files
       layout->write( "new_chip.lst" );      // writes out a .lst file for FastCap2
       layout->write( "new_chip.henry" );    // writes out a FastHenry2 files

       layout->write_layer_info( "new_chip.gds3d" ); // writes out layer info used by GDS3D viewer app
</pre>

<p>
There are a couple example programs in this directory:
<ul>
<li>count.cpp - parses a .gds file but doesn't store; prints out number of records found; runs fast on large .gds files; doit.count builds and runs it</li>
<li>test.cpp - can read in a .gds, .aedt, or .layout file, then writes out .layout, .gds, and .aedt files; doit.test builds and runs it</li>
</ul>
