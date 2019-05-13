Single Layout.h file for chip layout.  Layout.h provides:

<ul>
<li>efficient binary format (.layout) for both traversal and I/O</li>
<li>reading/writing .gds files (GDSII industry standard layout format)</li>
<li>reading/writing .aedt files for Ansys HFSS 3D field solver</li> 
<li>writing fastcap2 input files for capacitance 3D field solver</li>
<li>writing fasthenry2 input files for resistance+inductance 3D field solver</li>
</ul>

<p>
See instructions at the top of Layout.h for how to include the Layout class in your program.</p>

<p>
See eg/* for example programs that use Layout.h</p>
