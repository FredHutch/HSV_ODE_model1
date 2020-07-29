# HSV_ODE_model1
Initial C++ Stochastic ODE model for HSV infection (with OpenGL/GTK+ 2.0 GUI)

This repository contains a stochastic HSV model written in C++ for building on a Linux UBUNTU platform.
It requires GTKGLEXT, GTK 2.0 and OpenGL for its graphical user interface (GUI).  
This models was developed a number of years ago, so the required packages may not be available.
The only files treated specially are those from gtkglext.  The path to the installed libraries
and headers for that package must be set as the value of the GTKGLEXT environment variable prior to
running "make" in the Linux directory.

This model was used in support of these papers

    Rapid localized spread and immunologic containment define herpes simplex virus-2 reactivation in the human genital tract,Elife,Schiffer,et al,2013
    Rapid Viral Expansion and Short Drug Half-Life Explain the Incomplete Effectiveness of Current Herpes Simplex Virus 2-Directed Antiviral Agents,
    Antimicrobial Agents and Chemotherapy, Schiffer,et al,2013
    Herpes simplex virus-2 transmission probability estimates based on quantity of viral shedding,
    Journal of The Royal Society Interface, Schiffer,et al,2014
    
This model requires several libraries.  You will also need to load the Mesa module to compile and 
run it.  Depending on your UBUNTU version other modules may also be required including:

    GSL,GLib,GTK+,libGLU and Pango (version with pangox)

To help in case required versions of certain libraries are unavailable, some required libraries and header files are contained 
in the extra_libs directory of the GIT repository FredHutch/HSV_models.  The path to the extra_libs directory must be set as 
the value of the GTKGLEXT environment variable before running make in the Linux directory.

