module 2STREAM_Mod

      USE twostream_master_m

      implicit NONE

      PUBLIC 2STREAM_Run

      !  2STREAM input structure
      TYPE 2STREAM
        integer     :: NBEAMS = 1          ! Number of solar zenith angles
        integer     :: N_USER_STREAMS = 1  ! Number of Viewing zenith angles
        integer     :: N_USER_RELAZMS = 1  ! Number of relative azimuth angles
        integer     :: N_USER_LEVELS  = 1  ! Number of user-defined vertical output levels
        integer     :: N_USER_OBSGEOMS = 1 ! Number of azimuth angles calculated by surface supplement
      END TYPE 2STREAM
