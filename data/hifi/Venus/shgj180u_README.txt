CCSD3ZF0000100000001NJPL3KS0PDSX##mark##

/* File information */

PDS_VERSION_ID               = "PDS3"
FILE_NAME                    = "SHGJ180U.A01"

RECORD_TYPE                  = FIXED_LENGTH
RECORD_BYTES                 = 122
FILE_RECORDS                 =    16551
LABEL_RECORDS                =       79
^SHADR_HEADER_TABLE            =       80
^SHADR_COEFFICIENTS_TABLE      =       82

SPACECRAFT_NAME              = "MAGELLAN"
TARGET_NAME                  = "VENUS"
INSTRUMENT_NAME              = "RADIO SCIENCE SUBSYSTEM"
DATA_SET_ID                  = "MGN-V-RSS-5-GRAVITY-L2-V1.0"
OBSERVATION_TYPE             = "GRAVITY FIELD"
PRODUCT_ID                   = "RS-SHGJ180U.A01"
PRODUCT_RELEASE_DATE         = 1997-12-17
DESCRIPTION                  = "This file contains coefficients
and related data for the MGNP180U Spherical Harmonic model of 
the Venus gravity field.  Input data are from radio tracking of the
Magellan spacecraft.  This product is a set of two ASCII tables:
a header table and a coefficients table.  Definitions of the
tables follow.  The Magellan Venus gravity model is known as
the Spherical Harmonics Gravity ASCII Data Record (SHGADR) and
is produced by the Magellan Gravity Science Team at JPL under
the direction of W.L. Sjogren."

START_ORBIT_NUMBER           = 5758
STOP_ORBIT_NUMBER            = 15019
PRODUCT_CREATION_TIME        = 1997-12-17T00:00:00.000
PRODUCER_FULL_NAME           = "WILLIAM L. SJOGREN"
PRODUCER_INSTITUTION_NAME    = "JET PROPULSION LABORATORY"
PRODUCT_VERSION_TYPE         = "PRELIMINARY"
PRODUCER_ID                  = "MGN GRAVSCI TEAM"
SOFTWARE_NAME                = "SHAGRV;V1.0"

/* Structure Objects */

OBJECT               = SHADR_HEADER_TABLE
ROWS                       = 1
COLUMNS                    = 8
ROW_BYTES                  = 137
ROW_SUFFIX_BYTES           = 107
INTERCHANGE_FORMAT         = ASCII
DESCRIPTION                = "The SHADR header includes
descriptive information about the spherical harmonic 
coefficients which follow in SHADR_COEFFICIENTS_TABLE.  
The header consists of a single record of eight (delimited) 
data columns requiring 137 bytes, a pad of 105 unspecified 
ASCII characters, an ASCII carriage-return, and an ASCII 
line-feed."

  OBJECT                   = COLUMN
    NAME                         = "REFERENCE RADIUS"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 1
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "KILOMETER"
    DESCRIPTION                  = "The assumed reference 
    radius of the spherical planet."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "CONSTANT"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 25
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "N/A"
    DESCRIPTION                  = "For a gravity field model
    the assumed gravitational constant GM in km cubed
    per seconds squared for the planet.  For a topography
    model, set to 1."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "UNCERTAINTY IN CONSTANT"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 49
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "N/A"
    DESCRIPTION                  = "For a gravity field model
    the uncertainty in the gravitational constant GM in km
    cubed per seconds squared for the planet.  For a topography
    model, set to 0."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "DEGREE OF FIELD"
    DATA_TYPE                    = ASCII_INTEGER
    START_BYTE                   = 73
    BYTES                        = 5
    FORMAT                       = "I5"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The degree of model field."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "ORDER OF FIELD"
    DATA_TYPE                    = ASCII_INTEGER
    START_BYTE                   = 79
    BYTES                        = 5
    FORMAT                       = "I5"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The order of the model field."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "NORMALIZATION STATE"
    DATA_TYPE                    = ASCII_INTEGER
    START_BYTE                   = 85
    BYTES                        = 5
    FORMAT                       = "I5"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The normalization indicator.
    For gravity field:
        0   coefficients are unnormalized
        1   coefficients are normalized
        2   other."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "REFERENCE LONGITUDE"
    POSITIVE_LONGITUDE_DIRECTION = "EAST"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 91
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "DEGREE"
    DESCRIPTION                  = "The reference longitude for
    the spherical harmonic expansion; normally 0."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "REFERENCE LATITUDE"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 115
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "DEGREE"
    DESCRIPTION                  = "The reference latitude for
    the spherical harmonic expansion; normally 0."
  END_OBJECT               = COLUMN

END_OBJECT           = SHADR_HEADER_TABLE

OBJECT               = SHADR_COEFFICIENTS_TABLE
  ROWS                    = 16470
  COLUMNS                  = 6
  ROW_BYTES                = 107
  ROW_SUFFIX_BYTES         = 15
  INTERCHANGE_FORMAT       = ASCII
  DESCRIPTION              = "The SHADR coefficients table 
  contains the coefficients for the spherical harmonic model.
  Each row in the table contains the degree index m, the
  order index n, the coefficients Cmn and Smn, and the
  uncertainties in Cmn and Smn.  The (delimited) data 
  require 107 ASCII characters; these are followed by a pad 
  of 13 unspecified ASCII characters, an ASCII carriage-
  return, and an ASCII line-feed."

  OBJECT                   = COLUMN
    NAME                         = "COEFFICIENT DEGREE"
    DATA_TYPE                    = ASCII_INTEGER
    START_BYTE                   = 1
    BYTES                        = 5
    FORMAT                       = "I5"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The degree index m of the
    C and S coefficients in this record."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "COEFFICIENT ORDER"
    DATA_TYPE                    = ASCII_INTEGER
    START_BYTE                   = 7
    BYTES                        = 5
    FORMAT                       = "I5"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The order index n of the
    C and S coefficients in this record."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "C"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 13
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The coefficient Cmn
    for this spherical harmonic model."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "S"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 37
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The coefficient Smn
    for this spherical harmonic model."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "C UNCERTAINTY"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 61
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The uncertainty in the
    coefficient Cmn for this spherical harmonic model."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "S UNCERTAINTY"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 85
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The uncertainty in the
    coefficient Smn for this spherical harmonic model."
  END_OBJECT               = COLUMN

END_OBJECT           = SHADR_COEFFICIENTS_TABLE

END             
