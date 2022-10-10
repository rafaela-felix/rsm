# Reciprocal space maps generator

![downloads](https://img.shields.io/github/downloads/rafaela-felix/rsm/total)

Python and MATLAB versions of a simple module to generate three-dimensional reciprocal space maps from X-ray diffraction data.

## Download

Use 

```
git clone https://github.com/rafaela-felix/rsm.git
```

for obtaining the source codes.

### Requirements

The following requirements are needed for using the Python code:

* Python (version >= 3.6)
* Python packages (pip installable):
  * numpy (>= 1.22.0)
  * matplotlib (>= 3.5.0)
  * imageio (>= 2.9.0)
  * pyvista (>= 0.34.1) or any other library for mesh visualization (Paraview, vtkplotter, ParaView, Open3D, Vedo, etc)
  
## Usage

Inside the `rsm` folder or any directory containing the codes,  

* Call the function `ReciprocalSpaceMaps` in the MATLAB Command Window.
```
ReciprocalSpaceMaps(wl, twoth_d, phi_d, D, p, Ny, Nx, n0, m0, nb, mb, folder, fnamepng, sclfactor, Az, El, th0, dth)
``` 
* Import the module `rms`, then call it on a Python terminal.

```
from rms import rms

rsm(detector_info, scan_info, path, sclfactor, fnamepng, cmap, op, md, nd, el, az)
```

## Basic input

The following input is required for both codes:

* `path`: location of image files
* `fnamepng`: output file name

* **Detector parameters**:
  * `D`: sample-detector distance (mm)
  * `p`: pixel size (mm)
  * `Nx`: number of pixels along `xd` (horizontal direction)
  * `Ny`: number of pixels along `yd` (vertical direction)
  * `m0, n0`: central pixel - hit by the direct beam (horizontal and vertical coordinates)
  
* **Scan parameters**:
  * `wl`: wavelength (angstrom)
  * `twoth_d`: elevation of detector arm (deg)
  * `phi_d`: azimuth of detector arm (deg)
  * `th0`: initial `th` value (deg)
  * `dth`: increment in the rocking curve angle th (deg)
  
* **Plot parameters**:
  * `Az, El`: Azimuthal and elevation angles of 3D plot
  * `sclfactor`: isointensity curve levels (% of max value) 
  * `md, nd`: dead pixels ((horizontal and vertical coordinates, if there is)
  
The Python code has two additional (optional) parameters: `op` and `cmap` remaining for opacity and color of each isosurface. Also, the detector and scan parameters are expected as dictionaries (`detector_info` and `scan_info`).

# Citation

If you use the codes in your research, please consider citing the following work:

```
@article{3d_reciprocal_space_maps,
  title={A simple recipe to create three-dimensional reciprocal space maps},
  author={F. S. Penacchio, Rafaela and B. Estradiote, Mauricio and Morelhao, Sergio L. and Fornari, Celso I. and Kagerer, Philipp and Dittmar, Marco and Muller, Simon and Reinert, Friedrich},
  journal={},
  pages={},
  year={2022},
  doi={}
  publisher={}
}
```
