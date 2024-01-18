%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ARAEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Adaptable tRAjectory Evaluation tooL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
% AUTHOR: Alessio Derobertis
%-------------------------------------------------------------------------%
% GITHUB: https://github.com/deroxa99/arael
%-------------------------------------------------------------------------%
% INSTALLATION: download the 'arael' folder and place it in the MATLAB 
% working directory. Call the function 'arael.m' for propagating the orbit.
% More info is contained in the description of the 'arael.m' function. Some 
% examples are also available in the 'examples' folder.
%-------------------------------------------------------------------------%
% THIS FUNCTION MAKE USE OF THE SPICE TOOLKIT DEVELOPED BY NASA, LEARN
% MORE ABOUT SPICE AT: https://naif.jpl.nasa.gov/naif/toolkit.html
%-------------------------------------------------------------------------%
% Required toolbox: - symbolic manipulation toolbox
%-------------------------------------------------------------------------%
% RELEASE NOTES: Arael is an orbit propagater that allows the user to     
% choose between different propatgation methods, making use of various 
% levels of approximation of the dynamics in order to allow for a trade off
% between high-fidelity or high-performance simulations.                          
%                                                                         
% The selectable methods are:
%  - 'hifi'
%  - 'approx'
%  - 'averaged'
%  - 'full'
%
%------------------------------- HIFI ------------------------------------%
%  - GRAVITY: implement a complete harmonic expansion for the gravity 
%    field in order to simulate the dynamics in the most advanced way 
%    possible. The integration is carried out using Gauss Equinoctial 
%    elements and the supported central bodies are:
%    - 'EARTH': expansion up to order 360 (EGM96)
%               - Source: https://cddis.nasa.gov/926/egm96/
%    - 'MOON': expansion up to order 1200 (GRGM1200A)
%               - Source: https://pgda.gsfc.nasa.gov/products/50
%    - 'MARS': expansion up to order 80 (GGM1025A)
%               - Source: https://pds-ppi.igpp.ucla.edu/search/view/
%                         ?f=yes&id=pds://PPI/MGS-M-RSS-5-SDP-V1.0/DATA/
%                         RS_SHA/GGM1025A&o=1
%    - 'VENUS': expansion up to order 180 (MGN180u)
%               - Source: https://pds-geosciences.wustl.edu/mgn/
%                         mgn-v-rss-5-gravity-l2-v1/mg_5201/gravity/
%    To speed up the computation, pre-computed polynomials are available 
%    for expansions of order n = [25,50,60,75,100,150,200,300].
%
%  - THIRD BODY: Perturbation due to third bodies can be implemented for
%    the following bodies:
%    - 'SUN','MERCURY','VENUS','EARTH','MOON','MARS BARYCENTER',
%      'JUPITER BARYCENTER','SATURN BARYCENTER', 'URANUS BARYCENTER',
%      'NEPTUNE BARYCENTER','PLUTO' BARYCENTER.
%  
%  - AIR DRAG: not yet implemented
%
%  - SOLAR RADIATION PRESSURE: not yet implemented
%
%------------------------------ APPROX -----------------------------------%
%  - GRAVITY: implement an approximated model for the gravity field, to 
%    speed up the simulation. The integration is carried out using Gauss 
%    Equinoctial elements and up to now it is possible to simulate only 
%    zonal harmonics effect for the
%    following central bodies:
%     - 'EARTH': zonal harmonics up to J6;
%     - 'MOON': zonal harmonics up to J4;
%     - 'MARS': zonal harmonics up to J4;
%     - 'VENUS': zonal harmonics up to J4;
%
%  - THIRD BODY: Perturbation due to third bodies can be implemented for
%    the following bodies:
%    - 'SUN','MERCURY','VENUS','EARTH','MOON','MARS BARYCENTER',
%      'JUPITER BARYCENTER','SATURN BARYCENTER', 'URANUS BARYCENTER',
%      'NEPTUNE BARYCENTER','PLUTO' BARYCENTER.
%
%  - AIR DRAG: not yet implemented
%
%  - SOLAR RADIATION PRESSURE: not yet implemented
%
%----------------------------- AVERAGED ----------------------------------%
%  - GRAVITY: not yet implemented
%
%  - THIRD BODY: not yet implemented
%
%  - AIR DRAG: not yet implemented
%
%  - SOLAR RADIATION PRESSURE: not yet implemented
%
%------------------------------- FULL ------------------------------------%
%  - GRAVITY: implement an N-body problem where all the attractor are
%    treated as point masses. The integration is carried out using the
%    cartesian state. In this case perturb.n is not used. Allowed observer
%    bodies are:
%    - 'SUN'
%
%  - THIRD BODY: In this mode, the list of third bodies given in input is
%    used to compute the N-body dynamics. Allowed bodies are:
%    - 'SUN','MERCURY','VENUS','EARTH','MOON','MARS BARYCENTER',
%      'JUPITER BARYCENTER','SATURN BARYCENTER', 'URANUS BARYCENTER',
%      'NEPTUNE BARYCENTER','PLUTO' BARYCENTER.
%
%  - AIR DRAG: not yet implemented
%
%  - SOLAR RADIATION PRESSURE: not yet implemented
%
%-------------------------------------------------------------------------%
% Copyright (c) 2024 Alessio Derobertis
%-------------------------------------------------------------------------%