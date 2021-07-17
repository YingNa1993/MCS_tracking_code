MODULE irt_parameters

INTEGER, PARAMETER    :: domainsize_x = 498!497!
INTEGER, PARAMETER    :: domainsize_y = 284!249!

LOGICAL, PARAMETER    :: llonlatgrid = .FALSE.
REAL, PARAMETER       :: unit_area = 12345 ! km^2
! only used if llonlatgrid=.TRUE., otherwise set to arbitrary value:
REAL, PARAMETER       :: lat_first = 15.11718 !20.0390625 !
REAL, PARAMETER       :: lat_inc = 0.140625
REAL, PARAMETER       :: lon_inc = 0.140625

LOGICAL, PARAMETER    :: lperiodic_x = .FALSE.
LOGICAL, PARAMETER    :: lperiodic_y = .FALSE.

INTEGER, PARAMETER    :: n_fields = 1   ! number of additional averaging fields

! bins of coarse velocity field
INTEGER, PARAMETER    :: time_steps = 720    ! total number of timesteps
INTEGER, PARAMETER    :: nt_bins = 16        ! 6 hourly
INTEGER, PARAMETER    :: nx_bins = 1
INTEGER, PARAMETER    :: ny_bins = 1

REAL, PARAMETER       :: threshold = 241.0            ! for intensity
REAL, PARAMETER       :: threshold1 = 225.0            ! for intensity
REAL, PARAMETER       :: minimum_size = 60000    !changed     ! events smaller than that will be sorted out
REAL, PARAMETER       :: threshold2 = 2.0        !changed     ! for intensity
REAL, PARAMETER       :: minimum_size2 = 3600    !changed     ! events smaller than that will be sorted out
INTEGER, PARAMETER    :: minhour = 6 

REAL, PARAMETER       :: termination_sensitivity=1.0      ! Choose value between 0.0 and 1.0

REAL, PARAMETER       :: max_velocity = 25.      !changed   25*0.14*111=388km/h=108m/s=3.5degree/h  ! adjust acordingly
                                              ! velocities>max_velocity will be ignored to remove outliers
! define a minimal number of cells required for a coarse grained coordinate to be evaluated 
! if there are less, missing value will be assigned to that coarse cell
INTEGER, PARAMETER    :: min_cells = 100          !NOT USE

INTEGER, PARAMETER    :: max_no_of_cells=5000  ! buffer size, increase if necessary
INTEGER, PARAMETER    :: max_no_of_tracks=5000    ! buffer size, increase if necessary
INTEGER, PARAMETER    :: max_length_of_track=1000  ! buffer size, increase if necessary

REAL, PARAMETER       :: miss=-9999.           ! value<miss ==> missing_value

END MODULE irt_parameters

