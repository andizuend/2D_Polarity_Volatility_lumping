# 2D_Polarity_Volatility_lumping
This repository provides code-level access to a 2D framework for lumping of gas and/or aerosol components into an adjustable set of surrogate components to represent a complex system of hundreds to tens of thousands of components. The 2D space and related lumping approaches are based on the predicted polarity and volatility of the system's molecules.
This 2D framework is described in detail in a related scientific modeling article by Amaladhasan et al. (*in prep.*).
The framework includes the Aerosol Inorganic–Organic Mixtures Functional groups Activity Coefficients ([AIOMFAC](https://aiomfac.lab.mcgill.ca "AIOMFAC")) model (core code) to enable the computation of activity coefficient ratios as one option for expressing the polarity of organic molecules. The surrogate selection methods include grid-based sampling of the 2D polarity–volatility space or the use of the *k*-means clustering method. The latter is typically the recommended method; see Amaladhasan et al. (*in prep.*).  

## Dependencies
- `SMILES_to_AIOMFAC_inp` As part of the outputs from the framework, input files for the AIOMFAC(-web) model are generated for each surrogate system. This is done, in part, by calling the [SMILES_to_AIOMFAC (S2AS) tool](https://github.com/andizuend/S2AS__SMILES_to_AIOMFAC). For this, the 2D-lumping framework program expects to be located and accessible under the same parent directory as the 'SMILES_to_AIOMFAC_inp' tool folder.
- The `CustomizedPlots_Dislin` folder contains an additional Fortran program for the generation of plots using 2D framework output located in the folder `Output_lumping`. Such plots may require specific settings near the top of the `CustomizedPlotting` program .f90 file; e.g., for `lumpResChar`, `maxrows`, and the file range to be considered. Output from `CustomizedPlots_Dislin` will be located in its own `Output_Plots` subfolder.
- **Dislin** plotting library (optional). The plotting code also expects that Dislin is installed at its default location of `c:\dislin` (on a Windows machine). [Dislin download and installation information](https://www.dislin.de/index.html "Dislin") is provided on its website. For examples see also the [Dislin_x_y_plots repository](https://github.com/andizuend/Dislin_x_y_plot). If plotting with Dislin not is desired, the whole `CustomizedPlots_Dislin` folder can be removed from your local copy.

# Quick guide to running the 2D lumping framework
The 2D lumping framework can be run simply via an executable, but more conveniently also via an IDE like MS Visual Studio Community with Intel's oneAPI ifx Fortran compiler.
#### Necessary input files:
1. tbd (see details in the article by Amaladhasan et al.).

## Running a customized case
Once the input file and settings file has been prepared, you are ready to run the 2D lumping framework.
- Option 1: (MS Visual Studio on Windows)
	- TBA 
- Option 2: (via executable)
	- Open ...
