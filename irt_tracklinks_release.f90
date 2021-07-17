! needs the output of "irt objects" as input file
! writes out single cell tracks, ignoring joining and splitting
! Compile: gfortran -O0 -o irt_tracklinks_release.x irt_tracklinks_release.f90 irt_parameters.f90

PROGRAM irt_tracklinks

USE irt_parameters, ONLY: max_no_of_tracks, max_length_of_track, n_fields

IMPLICIT NONE

! track header info
INTEGER    :: track_length0(max_no_of_tracks,2)
INTEGER    :: first_timestep0(max_no_of_tracks,2)
INTEGER    :: track_length_ofmcs0(max_no_of_tracks,2)
INTEGER    :: starttimestep0(max_no_of_tracks,2)
INTEGER    :: status_beginning0(max_no_of_tracks,2)
INTEGER    :: status_end0(max_no_of_tracks,2)
INTEGER    :: track_id0(max_no_of_tracks,2)

! connecting track IDs
INTEGER    :: parent_track1(max_no_of_tracks,2)
INTEGER    :: parent_track2(max_no_of_tracks,2)
INTEGER    :: child_track1(max_no_of_tracks,2)
INTEGER    :: child_track2(max_no_of_tracks,2)

! status_beginning=0 : Track emerges by itself
! status_beginning=1 : Track is fragment of splitting event
! status_beginning=2 : Track is result of merging event
! status_beginning=-1: Track started by contact with missing values

! status_end=0 : Track dissipates
! status_end=1 : Track terminates by merging event
! status_end=2 : Track terminates by splitting event
! status_end=-1: Track terminates by contact with missing values

! track body info
INTEGER    :: track_id
INTEGER    :: cell_timestep(max_length_of_track,max_no_of_tracks,2)
INTEGER    :: cell_age1(max_length_of_track,max_no_of_tracks,2)
INTEGER    :: cell_age2(max_length_of_track,max_no_of_tracks,2)
INTEGER    :: cell_number(max_length_of_track,max_no_of_tracks,2)
INTEGER    :: cell_id(max_length_of_track,max_no_of_tracks,2)
INTEGER    :: mcsornot(max_length_of_track,max_no_of_tracks,2)
INTEGER    :: mcsornot1(max_length_of_track,max_no_of_tracks,2)
REAL       :: totarea(max_length_of_track,max_no_of_tracks,2)
REAL       :: totarea1(max_length_of_track,max_no_of_tracks,2)
REAL       :: totarea2(max_length_of_track,max_no_of_tracks,2)
REAL       :: totarea3(max_length_of_track,max_no_of_tracks,2)
REAL       :: prcp_mean1(max_length_of_track,max_no_of_tracks,2)
REAL       :: prcp_mean2(max_length_of_track,max_no_of_tracks,2)
REAL       :: prcp_mean3(max_length_of_track,max_no_of_tracks,2)
REAL       :: field_mean(max_length_of_track,n_fields+1,max_no_of_tracks,2)
REAL       :: field_min(max_length_of_track,n_fields+1,max_no_of_tracks,2)
REAL       :: field_max(max_length_of_track,n_fields+1,max_no_of_tracks,2)
INTEGER    :: xfirst(max_length_of_track,max_no_of_tracks,2)
INTEGER    :: xlast(max_length_of_track,max_no_of_tracks,2)
INTEGER    :: yfirst(max_length_of_track,max_no_of_tracks,2)
INTEGER    :: ylast(max_length_of_track,max_no_of_tracks,2)
REAL       :: center_of_mass_x(max_length_of_track,max_no_of_tracks,2)
REAL       :: center_of_mass_y(max_length_of_track,max_no_of_tracks,2)
REAL       :: xmax(max_length_of_track,max_no_of_tracks,2)
REAL       :: ymax(max_length_of_track,max_no_of_tracks,2)
REAL       :: xmin(max_length_of_track,max_no_of_tracks,2)
REAL       :: ymin(max_length_of_track,max_no_of_tracks,2)
REAL       :: velocity_x(max_length_of_track,max_no_of_tracks,2)
!REAL       :: velocity_y(max_length_of_track,max_no_of_tracks,2)
INTEGER    :: lfl(max_length_of_track,max_no_of_tracks,2)
INTEGER    :: slfl(max_length_of_track,max_no_of_tracks,2)
INTEGER    :: lbl(max_length_of_track,max_no_of_tracks,2)
INTEGER    :: slbl(max_length_of_track,max_no_of_tracks,2)

INTEGER    :: i1,i2,j
INTEGER    :: counter         ! track index of previous time step
INTEGER    :: counter_max(2)
INTEGER    :: lvl,ilvl1,ilvl2 ! time level 1 or 2 for actual time step (1-lvl for previous time step)
INTEGER    :: actual_time ! actual time step
LOGICAL    :: l_read_two_records, l_finished

CHARACTER (len=1)    :: sternchen

CHARACTER (len=15)   :: filenamestring

CHARACTER (len=150) :: input_filename 
CHARACTER (len=150) :: output_filename 

CALL getarg(1,filenamestring)

input_filename = "/disk1/nay/PE1/4_iterative_raincell_tracking/3_3_olrpermonth/gpm_1/output/irt_tracks_output_"&
                  //filenamestring//".txt"
output_filename ="/disk1/nay/PE1/4_iterative_raincell_tracking/3_3_olrpermonth/gpm_1/output/irt_tracklinks_output_"&
                  //filenamestring//".txt"  


OPEN(10,FILE=trim(input_filename),FORM='formatted', ACTION='read')
OPEN(20,FILE=trim(output_filename),FORM='formatted', ACTION='write')

counter=1 ! counter of tracks in memory
lvl=1     ! 1 or 2, level of actual record (3-lvl for precious record)
actual_time=1
l_read_two_records=.TRUE.
l_finished=.FALSE.
counter_max=0

!WRITE (*,*) "debug1", max_no_of_tracks

! beginning of main loop
DO

WRITE(*,*) "Read in ",max_no_of_tracks,"Level=",lvl

DO counter=1,max_no_of_tracks

! read in track header
READ (10,*,END=100) sternchen
READ (10,*) track_id0(counter,lvl),first_timestep0(counter,lvl),track_length0(counter,lvl), &
            starttimestep0(counter,lvl),track_length_ofmcs0(counter,lvl),&
            status_beginning0(counter,lvl),status_end0(counter,lvl)

parent_track1(counter,lvl)=0
parent_track2(counter,lvl)=0
child_track1(counter,lvl)=0
child_track2(counter,lvl)=0

! read in track record
DO j=1,track_length0(counter,lvl)
  !WRITE(*,*) j,counter,lvl
  READ(10,*) track_id,cell_timestep(j,counter,lvl),cell_id(j,counter,lvl),cell_number(j,counter,lvl), &
             cell_age1(j,counter,lvl),cell_age2(j,counter,lvl),&
             mcsornot(j,counter,lvl),mcsornot1(j,counter,lvl),&
             totarea(j,counter,lvl),totarea1(j,counter,lvl),totarea2(j,counter,lvl),totarea3(j,counter,lvl),&
             prcp_mean1(j,counter,lvl),prcp_mean2(j,counter,lvl),prcp_mean3(j,counter,lvl),&
	           field_mean(j,:,counter,lvl),field_min(j,:,counter,lvl),field_max(j,:,counter,lvl), &
             xfirst(j,counter,lvl),xlast(j,counter,lvl),yfirst(j,counter,lvl),ylast(j,counter,lvl), &
	           center_of_mass_x(j,counter,lvl),center_of_mass_y(j,counter,lvl), &
             xmax(j,counter,lvl),ymax(j,counter,lvl),&
             xmin(j,counter,lvl),ymin(j,counter,lvl),&
             velocity_x(j,counter,lvl), &
             lfl(j,counter,lvl),slfl(j,counter,lvl),lbl(j,counter,lvl),slbl(j,counter,lvl)
ENDDO

! max_no_of_tracks
ENDDO

!WRITE (*,*) "debug2"

counter_max(lvl)=max_no_of_tracks

IF (l_read_two_records) THEN
  lvl=3-lvl
  l_read_two_records = .FALSE.
  CYCLE
ENDIF

GOTO 200

100 CONTINUE
WRITE(*,*) "End of input reached, remaining no of tracks: ",counter-1
counter_max(lvl)=counter-1
l_finished=.TRUE.

200 CONTINUE

! find connection track ID numbers
!WRITE(*,*) "Establish links, lvl1:",counter_max(1),", lvl2:",counter_max(2)
!STOP

DO ilvl1=1,2
  DO i1=1,counter_max(ilvl1)
    DO ilvl2=1,2
      DO i2=1,counter_max(ilvl2)
        IF (ilvl1.EQ.ilvl2 .AND. i1.EQ.i2) CYCLE
        DO j=1,track_length0(i2,ilvl2)
          !WRITE(*,*) ilvl1,i1,ilvl1,i2,j
          IF (lbl(1,i1,ilvl1).EQ.cell_number(j,i2,ilvl2)) THEN 
            parent_track1(i1,ilvl1)=track_id0(i2,ilvl2)
          ENDIF
          IF (slbl(1,i1,ilvl1).EQ.cell_number(j,i2,ilvl2)) THEN
            parent_track2(i1,ilvl1)=track_id0(i2,ilvl2)
          ENDIF
          IF (lfl(track_length0(i1,ilvl1),i1,ilvl1).EQ.cell_number(j,i2,ilvl2)) THEN 
            child_track1(i1,ilvl1)=track_id0(i2,ilvl2)
          ENDIF
          IF (slfl(track_length0(i1,ilvl1),i1,ilvl1).EQ.cell_number(j,i2,ilvl2)) THEN 
            child_track2(i1,ilvl1)=track_id0(i2,ilvl2)
          ENDIF
        ENDDO ! j
      ENDDO ! i2
    ENDDO ! ilvl2
  ENDDO ! i1
ENDDO ! ilvl1

!WRITE(*,*) "Write out ",max_no_of_tracks,"Level=",3-lvl, counter_max(3-lvl)

DO counter=1,counter_max(3-lvl)

! write out track header
WRITE (20,*) "*"
WRITE (20,*) track_id0(counter,3-lvl),first_timestep0(counter,3-lvl),track_length0(counter,3-lvl), &
            track_length_ofmcs0(counter,3-lvl),starttimestep0(counter,3-lvl),&
            status_beginning0(counter,3-lvl),status_end0(counter,3-lvl), &
            parent_track1(counter,3-lvl),parent_track2(counter,3-lvl), &
            child_track1(counter,3-lvl),child_track2(counter,3-lvl)

! write out track record
DO j=1,track_length0(counter,3-lvl)
  WRITE(20,*) track_id0(counter,3-lvl),cell_timestep(j,counter,3-lvl), &
             cell_id(j,counter,3-lvl),cell_number(j,counter,3-lvl), &
             cell_age1(j,counter,3-lvl),cell_age2(j,counter,3-lvl),&
             mcsornot(j,counter,3-lvl),mcsornot1(j,counter,3-lvl),&
             totarea(j,counter,3-lvl),totarea1(j,counter,3-lvl),totarea2(j,counter,3-lvl),totarea3(j,counter,3-lvl), &
             prcp_mean1(j,counter,3-lvl),prcp_mean2(j,counter,3-lvl),prcp_mean3(j,counter,3-lvl),&
	           field_mean(j,:,counter,3-lvl),field_min(j,:,counter,3-lvl),field_max(j,:,counter,3-lvl), &
             xfirst(j,counter,3-lvl),xlast(j,counter,3-lvl),yfirst(j,counter,3-lvl),ylast(j,counter,3-lvl), &
	           center_of_mass_x(j,counter,3-lvl),center_of_mass_y(j,counter,3-lvl), &
             xmax(j,counter,3-lvl),ymax(j,counter,3-lvl),&
             xmin(j,counter,3-lvl),ymin(j,counter,3-lvl),&
             velocity_x(j,counter,3-lvl), &
             lfl(j,counter,3-lvl),slfl(j,counter,3-lvl),lbl(j,counter,3-lvl),slbl(j,counter,3-lvl)
ENDDO

! max_no_of_tracks
ENDDO

IF (l_finished) EXIT

! switch level
lvl=3-lvl

! end main loop
ENDDO

! write out remaining tracks

!WRITE(*,*) "Write out remaining ",counter_max,"Level=",lvl, counter_max(lvl)

DO counter=1,counter_max(lvl)

! write out track header
WRITE (20,*) "*"
WRITE (20,*) track_id0(counter,lvl),first_timestep0(counter,lvl),track_length0(counter,lvl), &
            track_length_ofmcs0(counter,lvl),starttimestep0(counter,lvl),&
            status_beginning0(counter,lvl),status_end0(counter,lvl), &
            parent_track1(counter,lvl),parent_track2(counter,lvl), &
            child_track1(counter,lvl),child_track2(counter,lvl)

! write out track record
DO j=1,track_length0(counter,lvl)
  WRITE(20,*) track_id0(counter,lvl),cell_timestep(j,counter,lvl), &
             cell_id(j,counter,lvl),cell_number(j,counter,lvl), &
             cell_age1(j,counter,lvl),cell_age2(j,counter,lvl),&
             mcsornot(j,counter,lvl),mcsornot1(j,counter,lvl),&
             totarea(j,counter,lvl),totarea1(j,counter,lvl),totarea2(j,counter,lvl),totarea3(j,counter,lvl), &
             prcp_mean1(j,counter,lvl),prcp_mean2(j,counter,lvl),prcp_mean3(j,counter,lvl),&
	           field_mean(j,:,counter,lvl),field_min(j,:,counter,lvl),field_max(j,:,counter,lvl), &
             xfirst(j,counter,lvl),xlast(j,counter,lvl),yfirst(j,counter,lvl),ylast(j,counter,lvl), &
	           center_of_mass_x(j,counter,lvl),center_of_mass_y(j,counter,lvl), &
             xmax(j,counter,lvl),ymax(j,counter,lvl),&
             xmin(j,counter,lvl),ymin(j,counter,lvl),&
             velocity_x(j,counter,lvl), &
             lfl(j,counter,lvl),slfl(j,counter,lvl),lbl(j,counter,lvl),slbl(j,counter,lvl)
ENDDO

! max_no_of_tracks
ENDDO


!200 CONTINUE

 CLOSE(20)
 CLOSE(10)

END PROGRAM irt_tracklinks


