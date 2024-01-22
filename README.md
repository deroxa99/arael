<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/deroxa99/arael">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a>

<h3 align="center">ARAEL</h3>

  <p align="center">
    Adaptable tRAjectory Evaluation tooL
    <br />
    <a href="https://github.com/deroxa99/arael"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/deroxa99/arael">View Demo</a>
    ·
    <a href="https://github.com/deroxa99/arael/issues">Report Bug</a>
    ·
    <a href="https://github.com/deroxa99/arael/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

[Product Name Screen Shot][product-screenshot]

'ARAEL' is a MATLAB toolbox for planetary and interplanetary trajectories simulation.

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- GETTING STARTED -->
## Getting Started

Arael is an orbit propagater that allows the user to choose between different propatgation methods, making use of various levels of approximation of the dynamics in order to allow for a trade off between high-fidelity or high-performance simulations.                                                                        
The selectable methods are:
- 'hifi'
- 'approx'
- 'averaged'
- 'full'

## Prerequisites

- **Spice**: This function makes use of the Spice toolkit developed by NASA, install spice from the [NASA](: https://naif.jpl.nasa.gov/naif/toolkit.html) website.
- **Toolboxes**: 
 - Symbolic Manipulation Toolbox

## Installation

1. Download the main folder 'arael-main.zip' and extract it into your MATLAB working directory.
2. Add the 'arael-main' folder and subfolder to the path.
3. Into your script, call the function 'arael.m' (se the 'arael.m' description for detailed input).
  > [!IMPORTANT]
  > The input is caps sensitive.
  > [!NOTE]
  > The input order is not 
4. The simulation will start and the output will be saved in the Workspace (se the 'arael.m' description for detailed output).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Modalities

### _HIFI_
- GRAVITY: implement a complete harmonic expansion for the gravity field in order to simulate the dynamics in the most advanced way possible. The integration is carried out using Gauss Equinoctial elements and the supported central bodies are:
 - 'EARTH': expansion up to order 360 (EGM96)
  - Source: https://cddis.nasa.gov/926/egm96/
 - 'MOON': expansion up to order 1200 (GRGM1200A)
  - Source: https://pgda.gsfc.nasa.gov/products/50
 - 'MARS': expansion up to order 80 (GGM1025A)
  - Source: https://pds-ppi.igpp.ucla.edu/search/view/?f=yes&id=pds://PPI/MGS-M-RSS-5-SDP-V1.0/DATA/RS_SHA/GGM1025A&o=1
 - 'VENUS': expansion up to order 180 (MGN180u)
  - Source: https://pds-geosciences.wustl.edu/mgn/mgn-v-rss-5-gravity-l2-v1/mg_5201/gravity/

- THIRD BODY: Perturbation due to third bodies can be implemented for the following bodies:
 - 'SUN','MERCURY','VENUS','EARTH','MOON','MARS BARYCENTER','JUPITER BARYCENTER','SATURN BARYCENTER', 'URANUS BARYCENTER','NEPTUNE BARYCENTER','PLUTO' BARYCENTER.

- AIR DRAG: not yet implemented

- SOLAR RADIATION PRESSURE: Perturbation due to solar radiation hitting the surface of the satellite. Umbra and penumbra regions due to the shadow of the primary attractor are also taken into account.

### _APPROX_
- GRAVITY: implement an approximated model for the gravity field, to speed up the simulation. The integration is carried out using Gauss Equinoctial elements and up to now it is possible to simulate only zonal harmonics effect for the following central bodies:
 - 'EARTH': zonal harmonics up to J6;
 - 'MOON': zonal harmonics up to J4;
 - 'MARS': zonal harmonics up to J4;
 - 'VENUS': zonal harmonics up to J4;
 
- THIRD BODY: Perturbation due to third bodies can be implemented for the following bodies:
  - 'SUN','MERCURY','VENUS','EARTH','MOON','MARS BARYCENTER','JUPITER BARYCENTER','SATURN BARYCENTER', 'URANUS BARYCENTER','NEPTUNE BARYCENTER','PLUTO' BARYCENTER.
 
- AIR DRAG: not yet implemented
 
- SOLAR RADIATION PRESSURE: Perturbation due to solar radiation hitting the surface of the satellite. Umbra and penumbra regions due to the shadow of the primary attractor are also taken into account.

### _AVERAGED_
- GRAVITY: not yet implemented

- THIRD BODY: not yet implemented

- AIR DRAG: not yet implemented

- SOLAR RADIATION PRESSURE: not yet implemented

## _FULL_
- GRAVITY: implement an N-body problem where all the attractor are treated as point masses. The integration is carried out using the cartesian state. In this case perturb.n is not used. Allowed observer bodies are:
 - 'SUN'

- THIRD BODY: In this mode, the list of third bodies given in input is used to compute the N-body dynamics. Allowed bodies are:
 - 'SUN','MERCURY','VENUS','EARTH','MOON','MARS BARYCENTER','JUPITER BARYCENTER','SATURN BARYCENTER', 'URANUS BARYCENTER','NEPTUNE BARYCENTER','PLUTO' BARYCENTER.

- AIR DRAG: not yet implemented

- SOLAR RADIATION PRESSURE: Perturbation due to solar radiation hitting the surface of the satellite. Shadow is not considered since the S/C is supposed to orbit the Sun.

<!-- USAGE EXAMPLES -->
## Usage

In the following section an example for the _'approx'_ mode is shown.

> [!WARNING]
> Make sure to place the main script into the same wotking folder as the _'arael-main'_ folder and add to path.

[Example screenshot][example-screenshot]

'''
%%% ARAEL EXAMPLE
clearvars
close all
clc

%% Data

% initial conditions
init_cond.x0 = [3.5786e04; 0; 0; 3.3374; 0; 0]; % state of the spacecraft [km,km/s]
init_cond.et = 0; % initial time in seconds after J2000 (2000/01/01 - 12:00:00.000)
init_cond.tSpan = 0:60:10*24*3600; % integration time-span

% gravity
perturb.n = 0;

% third body
perturb.TB = {'SUN','MOON'};

% solar radiation pressure
perturb.srp = 'on';

% reference system
ref_sys.inertial = 'J2000'; % inertial axes
ref_sys.obs = 'EARTH'; % center of the ref. system (observer)

% settings
settings.mode = 'approx'; % integration mode
settings.rel_tol = 1e-09;
settings.abs_tol = 1e-10;

% satellite
spacecraft.m = 850; % mass [kg]
spaceraft.A = 2; % surface [m^2]
spacecraft.cR = 1.8; % reflectivity coefficient [-]

%% integrate
[t,y] = arael(init_cond,ref_sys,perturb,settings);

%% post processing

.
.
.
'''

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Alessio Derobertis - [linkedin-shield][linkedin-url] - alessioderobertis@outlook.com

Project Link: [https://github.com/deroxa99/arael](https://github.com/deroxa99/arael)

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/deroxa99/arael.svg?style=for-the-badge
[contributors-url]: https://github.com/deroxa99/arael/graphs/contributors
[issues-shield]: https://img.shields.io/github/issues/deroxa99/arael.svg?style=for-the-badge
[issues-url]: https://github.com/deroxa99/arael/issues
[license-shield]: https://img.shields.io/github/license/deroxa99/arael.svg?style=for-the-badge
[license-url]: https://github.com/deroxa99/arael/blob/main/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/in/alessio-derobertis-4b4a831b7/
[product-screenshot]: images/arael.png
[example-screenshhot]: images/example.png
