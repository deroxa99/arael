function equi = kep2equi(kep)
    % ---------------------------------------------------------------------
    % KEP2EQUI - Convert keplerian elements into equinoctial elements
    % ---------------------------------------------------------------------
    % Convert a set of keplerian elements into one of equinoctial elements
    % ---------------------------------------------------------------------
    % INPUTS:
    % kep          : [6,1] - keplerian elements of the orbiter:
    %                        - SMA:   SEMI-MAJOR AXIS [km] 
    %                        - ECC:   ECCENTRICITY [-]
    %                        - INC:   INCLINATION [rad]
    %                        - LNODE: LONGITUDE OF THE ASCENDING NODE [rad]
    %                        - ARGP:  ARGUMENT OF PERIAPSIS [rad]
    %                        - THETA: TRUE ANOMALY [rad]
    %
    % OUTPUTS:
    % equi         : [6,1] - equinoctial elements of the orbiter. Ordered
    %                        in the following way: [p, f, g, h, k, L]
    % ---------------------------------------------------------------------
    % CONTRIBUTORS: Alessio Derobertis
    % ---------------------------------------------------------------------

    % unpack keplerian elements
    SMA = kep(1);
    ECC = kep(2);
    INC = kep(3);
    LNODE = kep(4);
    ARGP = kep(5);
    THETA = kep(6);

    % compute equinoctial elements
    omega = ARGP + LNODE;
    L = THETA + omega;
    p = SMA * (1 - ECC^2);
    f = ECC * cos(omega);
    g = ECC * sin(omega);
    tan_i2 = tan(INC/2);
    h = tan_i2 * cos(LNODE);
    k = tan_i2 * sin(LNODE);

    % assemble equinoctial elements
    equi = [p; f; g; h; k; L];
    
end
