Single Layout.h file for chip layout provides:

<ul>
<li>reading/writing .gds files (GDSII industry standard layout format)</li>
<li>reading/writing .aedt files for Ansys HFSS 3D field solver</li> 
<li>reading/writing .layout files (Layout.h binary representation)</li>
<li>writing fastcap2 input files for capacitance 3D field solver</li>
<li>writing fasthenry2 input files for resistance+inductance 3D field solver</li>
<li>instancing elements (structures, paths, etc.) from one or more source files</li>
<li>instancing entire source files</li>
<li>filling empty 3D space with an arbitrary material (typically a dielectric)</li>
<li>flattening instances into unique elements</li>
<li>writing out a single file with results<li>
</ul>

<p>
See instructions at the top of Layout.h for how to include the Layout class in your program.</p>

<p>
See eg/* for example programs that use Layout.h</p>
