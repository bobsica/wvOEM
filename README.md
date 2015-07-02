#wvOEM

####Table of Contents

1. [Description](#description)
2. [Setup - The basics of getting started with wvOEM](#setup)
    * [Setup requirements](#setup-requirements)
    * [Beginning with wvOEM](#beginning-with-wvOEM)
3. [Usage - Configuration options and additional functionality](#usage)
4. [Reference - An under-the-hood peek at what the module is doing and how](#reference)
    * [How to Run a Retrieval](#how-to-run-a-retrieval)
    * [Retrieved Parameters](#retrieved-parameters)
5. [Limitations - OS compatibility, etc.](#limitations)
6. [Development - Guide for contributing to the module](#development)
7. [Contributors](#contributors)
8. [Notes](#notes)

##Description

This matlab m-file and its subroutines retrieves water vapour mixing ratio, optical depth, lidar constants, Ångstrom exponent, dead time and background from RALMO lidar measurements. RALMO is described in detail in [this *Atmospheric Measurements Technique* paper.](http://www.atmos-meas-tech.net/6/1329/2013/amt-6-1329-2013.html)

##Setup
###Setup Requirements

1. **oem solver:** *atmlab-2.2.0* or higher is required for the OEM solver (*qpack*, Eriksson, P. et al. (2005), *J. of Quant Spec and Rad Trans,* 91(1), 47–64, doi:10.1016/j.jqsrt.2004.05.050). This free matlab package is available at [http://www.sat.ltu.se/arts/tools/](http://www.sat.ltu.se/arts/tools/) .
2. **additional matlab toolboxes required:** wvOEM using the function *smooth* in the curvefit toolbox. If *smooth* is not available you will have to write your own *N* point smoothing function. (/Applications/MATLAB_R2015a.app/toolbox/curvefit/curvefit/smooth.m).
3. **the following modules in the repository are from the [Mathworks File Exchange:](http://www.mathworks.com/matlabcentral/fileexchange/) curvefit.m, derivative.m, fitlinenp.m (fitline.m without plotting) and fwhmquiet.m (fwhm for return on error with a "-999" error flag).

###Beginning with wvOEM	

This code is not been cleaned for user friendliness. Most items which need to be configured are self explanatory and reside in the beginning of the script. However, there are some hard-coded quantities which would primarily be in makeQ2d0.m. For instance, nominal digital channel dead time is set at 4 ns in makeQ2d0.m. These assumptions will be updated in time.

##Usage

There are some self-explanatory logicals to set in the beginning of the file. Basically once these options are set you put data in and get retrieved parameters out.

##Reference

###How to Run a Retrieval

####Data Files

The following data files are required to live in /Users/BobSica/Dropbox/matlab/matlabWork/fromMCH/ralmodata/:

1. alphaBetayyyymmdd.mat
2. ralmoyyyymmddhhmmcombined.mat
3. sndyyyymmddhhmm.mat
4. S0S3yyyymmddhhmmchan2noE.mat

The first file is made by find_RayleighExtinct.m, a simple script to call Maxime's program that fetches the NCEP analysis for the day for a given location, then calculates the molecular extinction and backscatter coefficients. These are required to determine extinction from the RALMO aerosol scattering ratio.

The second and third files come from the Meteoswiss DataWarehouse, the former being the standard RALMO analysis with glued analog and digital channels, the latter the radiosonde measurements, typically at 0 and 12 UT.

The final file gets the raw measurements from the data files using a driver script I wrote (loopCalcSdata2chan.m) which calls Gianni's calcSdata_ral.m script, using the time interval specified in adt.conf.

###Retrieved Parameters

Information about OEM lidar retrieval for temperature can be found in Sica and Haefele (*Retrieval of temperature from a multiple-channel Rayleigh-scatter lidar using an optimal estimation method,* Appl. Opt., 54(8), 1872, doi:10.1364/AO.54.001872. This script applies the same principles to retrieval of the following quantities in the retrieval vector, *x*.

    % this version uses 4 data channels
    % retrieving 2*m + 10 parameters, where m is the length of the retrieval grid
    % x(1:m) is water vapour mixing ratio, q 
    % x(m+1:2*m) is optical depth 
    % x(end-9) is analog channel WV (water vapour) lidar constant
    % x(end-8) is analog channel N2 (nitrogen) lidar constant
    % x(end-7) is digital channel N2 lidar constant
    % x(end-6) is the Angstrom exponent 
    % x(end-5) is dead time for SH (water vapour counts)
    % x(end-4) is dead time for SN (nitrogen counts)
    % x(end-3) is H analog background 
    % x(end-2) is N analog background 
    % x(end-1) is H digital background 
    % x(end) is N digital background

##Limitations

This is not an operational code and has not be fully tested. It is completely specific to the Meteoswiss RALMO lidar. It should run on any OS that runs Matlab 2013 or higher, though it has only been tested with Matlab 2015a on Mac OSX (Yosemite).

##Development

This script is intended to be used with RALMO measurements only. To obtain RALMO measurements please contact [the Meteoswiss office in Payerne, Switzerland.](http://www.meteoswiss.admin.ch/home/about-us/contact.html)

##Contributors

R. Sica is to blame for the coding and holds copyright on all modules in this repository unless otherwise stated in the header. This project is a joint collaboration between [Alexander Haefele, Meteoswiss](http://4dweb.proclim.ch/4dcgi/proclim/en/Detail_Person?haefelea.27291) and [R. Sica, The University of Western Ontario](http://www.physics.uwo.ca/people/faculty_web_pages/sica.html).

##Notes

1 July 2015: this script (at version 2d0) is put into GIT versioning.
