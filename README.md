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

<h2>How to Use It</h2>

0. Clone this repository.  This code should build on Linux, MacOS, and Cygwin with a C++14 compiler.</li>

```
git clone https://github.com/balfieri/layout
```

1. Include Layout.h in your C++ program.

```
#include "Layout.h"
```

2. Import some .gds file and optionally write out our binary Layout binary which requires no translation.  
By the way, you may read in multiple .gds files at the same time.

```
Layout * layout = new Layout( "my_chip.gds" );   // creates a new layout with .gds imported into Layout's internal format
layout->write( "my_chip.layout" );               // optional: writes out Layout-format file (useful for large layouts)
```

3. Later, you can quickly read in the Layout-format file.  It essentially mmap's the file into memory with no translation required.

```
Layout * layout = new Layout( "my_chip.layout" ); // optional: useful only for large layouts
```

4. You can also export other types of files:
```
layout->write( "new_chip.gds" );                 // writes out a .gds II file
layout->write( "new_chip.aedt" );                // writes out an .aedt file for HFSS (includes layer info)
layout->write( "new_chip.fastcap" );             // writes out files for FastCap2, FFTCap 
layout->write( "new_chip.fasthenry" );           // writes out FastHenry2 files
layout->write( "new_chip.obj" );                 // writes out Alias/Wavefront 3D .obj/.mtl files 
(those last 4 are not fully implemented yet)

layout->write_layer_info( "new_chip.gds3d" );    // writes out layer mapping info used by GDS3D viewer app
layout->write_layer_info( "new_chip.tech" );     // writes out layer mapping info used by HFSS

layout->write_material_info( "new_chip.amat" );  // writes out material information used by HFSS
```

<h2>Example Programs</h2>

<p>
There are some example programs in this directory:

- ```count.cpp``` - parses a .gds file but doesn't build a tree; instead, prints out the number of records found; runs fast on large .gds files; ```doit.count``` builds and runs it

- ```test.cpp``` - can read in a .gds, .aedt, or .layout file, then writes out .layout, .gds, and .aedt files; ```doit.test``` builds and runs it

I need to add a more interesting test case that is not proprietary.


Bob Alfieri<br>
Chapel Hill, NC
