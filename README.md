Bob's tools for chip layout using single .h files:

<ul>
<li>Layout.h - Bob's layout format, plus reader/writer for .gdsii, .aedt (for HFSS), fastcap2, and fasthenry2 3D models
<ul>
    <li>libGDSII.h - GDSII reader/writer (from https://github.com/HomerRead/libGDSII, but in one .h file) used by Layout.h
</ul>
</ul>

<p>
See instructions at the top of Layout.h for how to include the Layout class in your program.</p>
<p>
See eg/* for example programs that use Layout.h</p>
