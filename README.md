# Table of Contents

- [Overview](#overview)
- [Basics](#basics)
- [Instantiating Layouts](#instantiating-layouts)
- [Advanced Layout Editing](#advanced-layout-editing)
- [Example Programs](#example-programs)

# Overview

This repo contains a single Layout.h file that provides the following VLSI layout operations:

<ul>
<li>reading and writing .gds files (GDSII industry standard layout format)</li>
<li>reading and writing .aedt files for Ansys HFSS 3D field solver</li> 
<li>reading and writing .layout files (our Layout.h binary representation which will soon be more compressed than .gds)</li>
<li>writing .fastcap input files for capacitance 3D field solver</li>
<li>writing .fasthenry input files for resistance+inductance 3D field solver</li>
<li>instancing structures or instances from one or more source files</li>
<li>instancing all instances from one or more source files</li>
<li>filling empty 3D space (within an arbitrary box) with a dielectric material</li>
<li>flattening instances into uniquified structures</li>
<li>writing out a single file with resultant edits</li>
</ul>

# Basics

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
layout->write( "new_chip.gds" );                 // writes out a .gds file
layout->write( "new_chip.raw" );                 // writes out a .raw file with internal node structure in ASCII
layout->write( "new_chip.aedt" );                // writes out an .aedt file for HFSS (includes layer info)
layout->write( "new_chip.fastcap" );             // writes out files for FastCap2, FFTCap 
layout->write( "new_chip.fasthenry" );           // writes out FastHenry2 files
layout->write( "new_chip.obj" );                 // writes out Alias/Wavefront 3D .obj/.mtl files 
(those last 4 are not fully implemented yet)

layout->write_layer_info( "new_chip.gds3d" );    // writes out layer mapping info used by GDS3D viewer app
layout->write_layer_info( "new_chip.tech" );     // writes out layer mapping info used by HFSS

layout->write_material_info( "new_chip.amat" );  // writes out material information used by HFSS
```

# Instantiating Layouts

TBD 

# Advanced Layout Editing

TBD

# Example Programs

<p>
There are some example programs in this directory:

- ```count.cpp``` - parses a .gds file but doesn't build a tree; instead, prints out the number of records found; runs fast on large .gds files; ```doit.count``` builds and runs it:
```
doit.count example1.gds
```

- ```test.cpp``` - can read in a .gds, .aedt, or .layout file, then writes out .layout, .gds, and .aedt files:
```
doit.test example1.gds
doit.test example1.aedt
doit.test example1.layout
```

- ```flatten.cpp``` - reads in a .gds/.aedt/.layout file, flattens it, then writes out a .flat.<ext> output:
```
doit.flatten example1.gds                       # flattened layout written to example1.flat.gds
doit.flatten example1.gds flattened.gds         # flattened layout written instead to flattened.gds
doit.flatten example1.gds flattened.layout      # flattened layout written instead to flattened.layout
```

I need to add an editing test case that is not proprietary.


Bob Alfieri<br>
Chapel Hill, NC
