Single Layout.h file provides:

<ul>
<li>reading and writing .gds files (GDSII industry standard layout format)</li>
<li>reading and writing .aedt files for Ansys HFSS 3D field solver</li> 
<li>reading and writing .layout files (our Layout.h binary representation)</li>
<li>writing fastcap2/fastercap/fftcap input files for capacitance 3D field solver</li>
<li>writing fasthenry2 input files for resistance+inductance 3D field solver</li>
<li>instancing structures or instances from one or more source files</li>
<li>instancing all instances from one or more source files</li>
<li>filling empty 3D space (within an arbitrary box) with a dielectric material</li>
<li>flattening instances into uniquified structures</li>
<li>writing out a single file with resultant edits</li>
</ul>

<p>
See instructions at the top of Layout.h for how to include the Layout class in your program.</p>

<p>
There are a couple example programs:
<ul>
<li>count.cpp - parses a .gds file but doesn't store; prints out number of records found; runs fast on large .gds files; doit.count builds and runs it</li>
<li>test.cpp - can read in a .gds, .aedt, or .layout file, then writes out .layout, .gds, and .aedt files; doit.test builds and runs it</li>
</ul>
