! identifies objects and establishes links
! needs input data as SRV file
! Compile: ifort -no-wrap-margin -o irt_objects_release.x irt_objects_release.f90 irt_parameters.f90

PROGRAM irt_objects

USE irt_parameters, ONLY: domainsize_x, domainsize_y, lperiodic_x, lperiodic_y, &
    n_fields, time_steps, nt_bins, nx_bins, ny_bins, threshold, threshold1, threshold2,&
    minimum_size,minimum_size2, max_no_of_cells, miss,max_velocity

IMPLICIT NONE

INTEGER              :: domsize_x, domsize_y

INTEGER              :: ii, ij,ix,iy,idx,idy
REAL, ALLOCATABLE    :: input_field(:,:,:)
REAL, ALLOCATABLE    :: overlay_field(:,:)
INTEGER, ALLOCATABLE :: event_number(:,:,:)
LOGICAL, ALLOCATABLE :: occupied(:,:)
LOGICAL, ALLOCATABLE :: occupied_consecuprcp(:,:)
LOGICAL, ALLOCATABLE :: occupied_consecuprcp1(:,:)
INTEGER              :: i,j,a,b,c,d,e,f,t,it,fileid
INTEGER              :: xxx,yyy,xxx1,yyy1,maxl(2)
INTEGER              :: counter_actual, counter_previous, counter_mixed
INTEGER              :: counter_total_actual, counter_total_previous
LOGICAL              :: delete_cell
REAL                 :: totarea_sum
REAL                 :: xx1,yy1,xx2,yy2

! for cell statistics
REAL                 :: totarea(max_no_of_cells,2)      !CCS area
REAL                 :: totarea1(max_no_of_cells,2)     !changed area of prcp>0 in CCS
REAL                 :: totarea2(max_no_of_cells,2)     !changed area of prcp>=2mm/h in CCS
REAL                 :: totarea3(max_no_of_cells,2)     !changed area of consecutive prcp>=2mm/h from maxloc
!REAL                 :: totarea3(max_no_of_cells,2)    !changed area of consecutive prcp>=2mm/h from COM

REAL                 :: prcp_mean1(max_no_of_cells,2)   !changed prcp>0 in CCS
REAL                 :: prcp_mean2(max_no_of_cells,2)   !changed prcp>=2mm/h in CCS
REAL                 :: prcp_mean3(max_no_of_cells,2)   !changed consecutive prcp>=2mm/h
REAL                 :: pmax(max_no_of_cells,2), tbmin(max_no_of_cells,2)
REAL                 :: field_mean(max_no_of_cells,2,n_fields+1)
REAL                 :: field_min(max_no_of_cells,2,n_fields+1)
REAL                 :: field_max(max_no_of_cells,2,n_fields+1)
REAL                 :: center_of_mass_x(max_no_of_cells,2),center_of_mass_y(max_no_of_cells,2) ! of prcp>0
INTEGER              :: cell_age1(max_no_of_cells,2), cell_age2(max_no_of_cells,2)
INTEGER              :: first_point_x(max_no_of_cells,2),first_point_y(max_no_of_cells,2)
INTEGER              :: first_point_x1(max_no_of_cells,2),first_point_y1(max_no_of_cells,2)
INTEGER              :: first_point_x2(max_no_of_cells,2),first_point_y2(max_no_of_cells,2)
INTEGER              :: xfirst(max_no_of_cells,2), xlast(max_no_of_cells,2)
INTEGER              :: yfirst(max_no_of_cells,2), ylast(max_no_of_cells,2)
INTEGER              :: xmax(max_no_of_cells,2), ymax(max_no_of_cells,2)
INTEGER              :: xmin(max_no_of_cells,2), ymin(max_no_of_cells,2)

! variables for the overlay field
INTEGER              :: xfirst_mixed,xlast_mixed,yfirst_mixed,ylast_mixed

! for link statistics
INTEGER              :: largest_forward_link(max_no_of_cells)
REAL                 :: largest_forward_link_size(max_no_of_cells)
INTEGER              :: second_largest_forward_link(max_no_of_cells)
REAL                 :: second_largest_forward_link_size(max_no_of_cells)
INTEGER              :: largest_backward_link(max_no_of_cells)
REAL                 :: largest_backward_link_size(max_no_of_cells)
INTEGER              :: second_largest_backward_link(max_no_of_cells)
REAL                 :: second_largest_backward_link_size(max_no_of_cells)

! diagnosed cell velocities
REAL                 :: velocity_x(max_no_of_cells)
REAL                 :: velocity_y(max_no_of_cells)

! dummy variables for velocity diagnosis
REAL                 :: area_weight,dx,dy

INTEGER              :: srv_header_input(8)
INTEGER              :: srv_header_coarse(8)

! file names
CHARACTER (len=45)   :: filepath
CHARACTER (len=150)   :: input_filename(n_fields+1)
CHARACTER (len=150)   :: output_filename
CHARACTER (len=150)   :: mask_filename
CHARACTER (len=150)   :: coarsevel_filename
CHARACTER (len=15)    :: filenamestring !nicam_cn_200403 

CHARACTER (len=1)    :: iteration_str
INTEGER              :: headlen

! time handling
INTEGER              :: n_actual,n_previous,n_previous_it
LOGICAL              :: previous_exists
INTEGER              :: previousdate, previoustime, previoustimestep

! which iteration?
INTEGER              :: iteration = 1   ! 1: first iteration
                                        ! 2: second iteration

! for advection velocity field
REAL                 :: coarse_vel_x(nx_bins,ny_bins,nt_bins)
REAL                 :: coarse_vel_y(nx_bins,ny_bins,nt_bins)
REAL                 :: coarse_dummy(nx_bins,ny_bins)
REAL                 :: vx,vy,vv



filepath = '/disk1/nay/PE1/4_iterative_raincell_tracking/'
CALL getarg(1,filenamestring)
!filenamestring1="Tb_"//filenamestring
input_filename(1) = '/disk1/nay/GPM_data/inputsrvdata/irt_objects_input_Tb_'//filenamestring//'.srv'
input_filename(2) = '/disk1/nay/GPM_data/inputsrvdata/irt_objects_input_'//filenamestring//'.srv'

!!not use!! DO i=0,n_fields
!!not use!!  WRITE(input_filename(i+1),'(A18,I2.2,A4)") "irt_objects_input_",i,".srv"
!!not use!! ENDDO

output_filename    = filepath//'3_3_olrpermonth/gpm_1/output/irt_objects_output_'//filenamestring//'.txt'
mask_filename      = filepath//'3_3_olrpermonth/gpm_1/output/irt_objects_mask_'//filenamestring//'.srv'
coarsevel_filename = filepath//'3_3_olrpermonth/gpm_1/output/irt_advection_field_'//filenamestring//'.srv'

  iteration = 1

!CALL getarg(1,iteration_str)
!IF (iteration_str .EQ. "1") THEN
!   iteration = 1
!ELSEIF (iteration_str .EQ. "2") THEN
!   iteration = 2
!ELSE
!  WRITE (*,*) "ERROR: 1 for first iteration, 2 for subsequent iterations"
!  STOP
!ENDIF


DO fileid=1,n_fields+1
  OPEN(fileid,FILE=trim(input_filename(fileid)),FORM='unformatted', ACTION='read')
ENDDO

OPEN(20,FILE=trim(output_filename),FORM='formatted', ACTION='write')
OPEN(30,FILE=trim(mask_filename),FORM='unformatted', ACTION='write')

! if not first iteration, read in the coarse velocity field
! deleted


n_actual   = 1
n_previous = 2

previous_exists = .FALSE.  ! because it's the first time step
counter_actual = 0
counter_total_actual = 0
counter_total_previous = 0
!previoustimestep=-1

! If periodic boundary conditions are switched off, a 1-gridbix thick 
! frame of missing values will be laid around the field, and the domainsize
! has to be increased by 2 in both dimensions
domsize_x = domainsize_x+2
domsize_y = domainsize_y+2

ALLOCATE(input_field(domsize_x,domsize_y,n_fields+1))
ALLOCATE(overlay_field(domsize_x,domsize_y))
ALLOCATE(event_number(domsize_x,domsize_y,4))
ALLOCATE(occupied(domsize_x,domsize_y))
ALLOCATE(occupied_consecuprcp(domsize_x,domsize_y))
ALLOCATE(occupied_consecuprcp1(domsize_x,domsize_y))

!!!!! beginning of main loop !!!!!!!!

DO previoustimestep=0,time_steps-1

!write(*,*) previoustimestep
  previousdate=srv_header_input(3)  ! header(3) : Date in the form: YYYYMMDDhhmmss
  previoustime=srv_header_input(4)  ! header(4) : Time increment
  !previoustimestep=previoustimestep+1

  DO fileid=1,n_fields+1
    READ (fileid,END=200) srv_header_input
    input_field(:,:,fileid) = miss-1.
    READ (fileid) input_field(2:domsize_x-1,2:domsize_y-1,fileid)
  ENDDO

  occupied(:,:)=.FALSE.
  occupied_consecuprcp(:,:)=.FALSE.
  occupied_consecuprcp1(:,:)=.FALSE.

  event_number(:,:,n_actual)=0    !!! n_actual   = 1
  counter_previous = MAX(1,counter_actual)   !!! counter_actual = 0
  counter_actual=1  !nevent

  ! identification of patches
  DO iy=1, domsize_y
    DO ix=1, domsize_x
      IF (input_field(ix,iy,1).LT.threshold .AND. input_field(ix,iy,1).GT.(miss+100) .AND. .NOT. occupied(ix,iy)) THEN
        ii = ix
        ij = iy
        totarea(counter_actual,n_actual)=0
        totarea1(counter_actual,n_actual)=0
        totarea2(counter_actual,n_actual)=0
        totarea3(counter_actual,n_actual)=0
        prcp_mean1(counter_actual,n_actual)=0
        prcp_mean2(counter_actual,n_actual)=0
        prcp_mean3(counter_actual,n_actual)=0        
        field_mean(counter_actual,n_actual,:)=0
        field_min(counter_actual,n_actual,:)=1E+9
        field_max(counter_actual,n_actual,:)=0
        center_of_mass_x(counter_actual,n_actual)= 0
        center_of_mass_y(counter_actual,n_actual)= 0
        xfirst(counter_actual,n_actual)=ii
        xlast(counter_actual,n_actual)=ii
        yfirst(counter_actual,n_actual)=ij
        ylast(counter_actual,n_actual)=ij
        xmax(counter_actual,n_actual)= 0 
        ymax(counter_actual,n_actual)= 0 
        pmax(counter_actual,n_actual)=0 
        xmin(counter_actual,n_actual)= 0 
        ymin(counter_actual,n_actual)= 0 
        tbmin(counter_actual,n_actual)=225

        delete_cell = .FALSE.

        CALL area(ii, ij, domsize_x,domsize_y,input_field,occupied, counter_actual, &
          event_number(:,:,n_actual),totarea(counter_actual,n_actual),&
          totarea1(counter_actual,n_actual),totarea2(counter_actual,n_actual),&
          prcp_mean1(counter_actual,n_actual),prcp_mean2(counter_actual,n_actual),&
          field_mean(counter_actual,n_actual,:),field_min(counter_actual,n_actual,:),&
          field_max(counter_actual,n_actual,:), &
          center_of_mass_x(counter_actual,n_actual),center_of_mass_y(counter_actual,n_actual), &
          delete_cell,xfirst(counter_actual,n_actual),xlast(counter_actual,n_actual), &
          yfirst(counter_actual,n_actual),ylast(counter_actual,n_actual),&
          xmax(counter_actual,n_actual),ymax(counter_actual,n_actual),pmax(counter_actual,n_actual),&
          xmin(counter_actual,n_actual),ymin(counter_actual,n_actual),tbmin(counter_actual,n_actual))


        IF (delete_cell) THEN
        ! delete this cell by overwriting it with -1 in event_number
 	        CALL set_event_number_to_value(domsize_x,domsize_y,-1,event_number(:,:,n_actual),counter_actual, &
	             xfirst(counter_actual,n_actual),xlast(counter_actual,n_actual), &
				 yfirst(counter_actual,n_actual),ylast(counter_actual,n_actual))

        ELSEIF (counter_actual .GT. max_no_of_cells) THEN
          WRITE(*,*) "ERROR: number of cells >",max_no_of_cells
          STOP

        !!!!!!!!!!!!!!!!!!!!!!!!changed!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ELSEIF (field_min(counter_actual,n_actual,1) .GE. 225) THEN
        ! delete this cell by overwriting it with -1 in event_number
          CALL set_event_number_to_value(domsize_x,domsize_y,-1,event_number(:,:,n_actual),counter_actual, &
               xfirst(counter_actual,n_actual),xlast(counter_actual,n_actual), &
               yfirst(counter_actual,n_actual),ylast(counter_actual,n_actual))


        ELSEIF (totarea(counter_actual,n_actual) .GT. 0.1) THEN  !!!!!changed

          field_mean(counter_actual,n_actual,:)=field_mean(counter_actual,n_actual,:)/ &
                                                totarea(counter_actual,n_actual)
          IF (prcp_mean1(counter_actual,n_actual) .GT. 0.1) THEN                                      
          center_of_mass_x(counter_actual,n_actual)=center_of_mass_x(counter_actual,n_actual)/ &
                                                    prcp_mean1(counter_actual,n_actual)
          center_of_mass_y(counter_actual,n_actual)=center_of_mass_y(counter_actual,n_actual)/ &
                                                    prcp_mean1(counter_actual,n_actual)

          ENDIF                                                                                          
          IF (totarea1(counter_actual,n_actual) .GT. 0.1) THEN                                      
          prcp_mean1(counter_actual,n_actual)=prcp_mean1(counter_actual,n_actual)/ &
                                            totarea1(counter_actual,n_actual) 
          ENDIF  
          IF (totarea2(counter_actual,n_actual) .GT. 0.1) THEN                                 
          prcp_mean2(counter_actual,n_actual)=prcp_mean2(counter_actual,n_actual)/ &  
                                            totarea2(counter_actual,n_actual)  !!!!!!!!!!!!!!changed   
          ENDIF  

        	! take care for periodic boundary conditions
        	IF (center_of_mass_x(counter_actual,n_actual) .GE. domsize_x+1) THEN
        	  center_of_mass_x(counter_actual,n_actual)=center_of_mass_x(counter_actual,n_actual)-domsize_x
        	ENDIF
        	IF (center_of_mass_x(counter_actual,n_actual) .LE. 1) THEN
        	  center_of_mass_x(counter_actual,n_actual)=center_of_mass_x(counter_actual,n_actual)+domsize_x
        	ENDIF
        	IF (center_of_mass_y(counter_actual,n_actual) .GE. domsize_y+1) THEN
        	  center_of_mass_y(counter_actual,n_actual)=center_of_mass_y(counter_actual,n_actual)-domsize_y
        	ENDIF
        	IF (center_of_mass_y(counter_actual,n_actual) .LE. 1) THEN
        	  center_of_mass_y(counter_actual,n_actual)=center_of_mass_y(counter_actual,n_actual)+domsize_y
        	ENDIF

          IF (pmax(counter_actual,n_actual).GT.0.1) THEN
          xxx1 = xmax(counter_actual,n_actual)
          yyy1 = ymax(counter_actual,n_actual)
          CALL area_consecuprcp(xxx1,yyy1,domsize_x,domsize_y,input_field, &
                    occupied_consecuprcp1,totarea3(counter_actual,n_actual),prcp_mean3(counter_actual,n_actual),delete_cell)  

          IF (totarea3(counter_actual,n_actual) .GT. 0.1) THEN
            prcp_mean3(counter_actual,n_actual)=prcp_mean3(counter_actual,n_actual)/ &
                                                totarea3(counter_actual,n_actual) 
          ENDIF
      	  ENDIF

          !!!!changed
          first_point_x1(counter_actual,n_actual)=center_of_mass_x(counter_actual,n_actual)
          first_point_y1(counter_actual,n_actual)=center_of_mass_y(counter_actual,n_actual)
          first_point_x2(counter_actual,n_actual)=xmin(counter_actual,n_actual)
          first_point_y2(counter_actual,n_actual)=ymin(counter_actual,n_actual)         
          first_point_x(counter_actual,n_actual)=ix
          first_point_y(counter_actual,n_actual)=iy
          cell_age1(counter_actual,n_actual)=1   ! new born cell. Will be adjusted later
          cell_age2(counter_actual,n_actual)=1
          counter_actual=counter_actual+1


        ELSE   ! if cell is too small, delete it
	        CALL set_event_number_to_value(domsize_x,domsize_y,-1,event_number(:,:,n_actual),counter_actual, &
	             xfirst(counter_actual,n_actual),xlast(counter_actual,n_actual), &
				 yfirst(counter_actual,n_actual),ylast(counter_actual,n_actual))
        ENDIF

      ENDIF
    ENDDO
  ENDDO


  CALL write_srv(REAL(event_number(2:domsize_x-1,2:domsize_y-1,n_previous)),previousdate,previoustime,30)


  counter_total_previous = counter_total_actual
  counter_total_actual = counter_total_actual+counter_previous-1

  ! do we have to read in the next step also?
  IF (.NOT. previous_exists) THEN
    previous_exists = .TRUE.
    ! flip time indices and advance 1 time step
    n_actual   = 3-n_actual
    n_previous = 3-n_previous
    !WRITE (*,*) "JUMP"
    CYCLE
  ENDIF

! perform second iteration on "previous" field
! deleted
  n_previous_it = n_previous


! calculate events for the mixed field
overlay_field = 0
DO iy=1,domsize_y
  DO ix=1,domsize_x

    IF (input_field(ix,iy,1) .LE. miss) THEN
      overlay_field(ix,iy) = miss-1
    ELSE IF(event_number(ix,iy,n_actual) .NE. 0 .OR. event_number(ix,iy,n_previous_it) .NE. 0) THEN
      overlay_field(ix,iy) = 1.
    ENDIF
  ENDDO
ENDDO


occupied(:,:)=.FALSE.
occupied_consecuprcp(:,:)=.FALSE.
occupied_consecuprcp1(:,:)=.FALSE.

counter_mixed=1

!identification of patches for the overlaid field

DO iy=1, domsize_y
  DO ix=1, domsize_x
    IF (overlay_field(ix,iy) .GE. 0.5 .AND. .NOT. occupied(ix,iy)) THEN
      ii = ix
      ij = iy
      xfirst_mixed=ii
      xlast_mixed=ii
      yfirst_mixed=ij
      ylast_mixed=ij
      delete_cell = .FALSE.
      CALL overlay_area(miss, ii, ij, domsize_x,domsize_y, overlay_field, &
                       occupied,counter_mixed,event_number(:,:,3),  &
	                     delete_cell,xfirst_mixed,xlast_mixed,yfirst_mixed,ylast_mixed)

      IF (delete_cell) THEN
        CALL set_event_number_to_value(domsize_x,domsize_y,0,event_number(:,:,3),counter_mixed, &
	                               xfirst_mixed,xlast_mixed,yfirst_mixed,ylast_mixed)
      ELSEIF (counter_mixed .LE. max_no_of_cells) THEN
        counter_mixed=counter_mixed+1                !! not use after here
      ELSE
        WRITE(*,*) "ERROR: number of cells >",max_no_of_cells
        STOP
      ENDIF
    ENDIF
  ENDDO
ENDDO

!WRITE (*,*) 'timestep:',previoustimestep,', time:',srv_header_input(4)


! identify foreward links and velocity
largest_forward_link=0
largest_forward_link_size=0
second_largest_forward_link=0
second_largest_forward_link_size=0
DO i=1,counter_previous-1
  !original code:if first point does not belong to overlay area,it can not link objects to track.add center location
  !first_point_x1 : center of prcp mass x first_point_x:  x of min tb 
  a = event_number(first_point_x(i,n_previous),first_point_y(i,n_previous),3) 
  d = event_number(first_point_x1(i,n_previous),first_point_y1(i,n_previous),3) 
  e = event_number(first_point_x2(i,n_previous),first_point_y2(i,n_previous),3) 

  IF (a==0) THEN 
    largest_forward_link(i) = -1    ! track interrupted by missing values
    !WRITE(*,*) "value 0:",first_point_x(i,n_previous),first_point_y(i,n_previous)
    CYCLE
  ENDIF

  velocity_x(i) = 0.
  velocity_y(i) = 0.
  vv=0.
  !area_weight = 0.
  DO j=1,counter_actual-1
    b = event_number(first_point_x(j,n_actual),first_point_y(j,n_actual),3)
    c = event_number(first_point_x1(j,n_actual),first_point_y1(j,n_actual),3)
    f = event_number(first_point_x2(j,n_actual),first_point_y2(j,n_actual),3)

    IF (b==0) CYCLE
    IF (a==b.or.d==c.or.e==f.or.d==f.or.e==c) THEN  !!!changed


      ! take care for cyclic boundary conditions:
      !IF (dx .GT. domsize_x/2)  dx = dx-domsize_x
      !IF (dx .LT. -domsize_x/2) dx = dx+domsize_x
      !IF (dy .GT. domsize_y/2)  dy = dy-domsize_y
      !IF (dy .LT. -domsize_y/2) dy = dy+domsize_y
      !velocity_x(i) = velocity_x(i)+dx*totarea(j,n_actual)
      !velocity_y(i) = velocity_y(i)+dy*totarea(j,n_actual)
      !area_weight = area_weight+totarea(j,n_actual)
      if (center_of_mass_x(j,n_actual).gt.0 .and. center_of_mass_y(j,n_actual).gt.0) then 
      	xx1=center_of_mass_x(j,n_actual)
      	yy1=center_of_mass_y(j,n_actual)
      else 
      	xx1=xmin(j,n_actual)
      	yy1=ymin(j,n_actual)
      endif
      if (center_of_mass_x(i,n_previous).gt.0 .and. center_of_mass_y(i,n_previous).gt.0) then 
      	xx2=center_of_mass_x(i,n_previous)
      	yy2=center_of_mass_y(i,n_previous)
      else 
      	xx2=xmin(i,n_previous)
      	yy2=ymin(i,n_previous)
      endif

      dx = xx1 - xx2
      dy = yy1 - yy2
      vv = sqrt(dx*dx+dy*dy)

      IF (vv .LT. max_velocity) THEN

      IF (totarea(j,n_actual) .GT. largest_forward_link_size(i)) THEN
        second_largest_forward_link(i) = largest_forward_link(i)
	    second_largest_forward_link_size(i) = largest_forward_link_size(i)
	    largest_forward_link(i) = counter_total_actual+j
	    largest_forward_link_size(i) = totarea(j,n_actual)
	    velocity_x(i) = vv
      ELSEIF (totarea(j,n_actual) .GT. second_largest_forward_link_size(i)) THEN
		second_largest_forward_link(i) = counter_total_actual+j
		second_largest_forward_link_size(i) = totarea(j,n_actual)
      ENDIF
      ENDIF

    ENDIF
      
  ENDDO
  !IF (area_weight.EQ.0) WRITE (*,*) "area_weight=0",counter_actual,velocity_x(i)
  ! determine cell velocity if possible, otherwise assign missing value
  !IF (area_weight .GE. minimum_size) THEN
  !  velocity_x(i) = velocity_x(i)/area_weight
  !  velocity_y(i) = velocity_y(i)/area_weight
  !ELSE
  !  velocity_x(i) = miss
  !  velocity_y(i) = miss
  !ENDIF

ENDDO

! write sizes
IF (counter_previous==1) THEN
  !IF (cell_age1(1,n_previous).LT.cell_age2(1,n_previous)) THEN
  !  WRITE(*,*) "cell age(cp=1):",cell_age1(1,n_previous),cell_age2(1,n_previous)
  !ENDIF
  WRITE(20,*) previoustimestep,0,counter_total_previous+i, &
        cell_age1(1,n_previous),cell_age2(1,n_previous),totarea(1,n_previous), &
        totarea1(1,n_previous),totarea2(1,n_previous),totarea3(1,n_previous),&
        prcp_mean1(1,n_previous),prcp_mean2(1,n_previous),prcp_mean3(1,n_previous),&
	    field_mean(1,n_previous,:),field_min(1,n_previous,:),field_max(1,n_previous,:), &
	    xfirst(1,n_previous),xlast(1,n_previous),yfirst(1,n_previous),ylast(1,n_previous), &
	    center_of_mass_x(1,n_previous),center_of_mass_y(1,n_previous), &
	    xmax(1,n_previous),ymax(1,n_previous), &
	    xmin(1,n_previous),ymin(1,n_previous), &
	    velocity_x(1), &
	    largest_forward_link(1),largest_forward_link_size(1), &
	    second_largest_forward_link(1),second_largest_forward_link_size(1), &
	    largest_backward_link(1),largest_backward_link_size(1), &
	    second_largest_backward_link(1),second_largest_backward_link_size(1)
ENDIF

DO i=1,counter_previous-1
  !IF (cell_age1(i,n_previous).LT.cell_age2(i,n_previous)) THEN
  !  WRITE(*,*) "cell age:",cell_age1(i,n_previous),cell_age2(i,n_previous)
  !ENDIF
  WRITE(20,*) previoustimestep,i,counter_total_previous+i, &
          cell_age1(i,n_previous),cell_age2(i,n_previous),totarea(i,n_previous), &
          totarea1(i,n_previous),totarea2(i,n_previous), totarea3(i,n_previous), &
          prcp_mean1(i,n_previous),prcp_mean2(i,n_previous),prcp_mean3(i,n_previous),&
	      field_mean(i,n_previous,:),field_min(i,n_previous,:),field_max(i,n_previous,:), &
	      xfirst(i,n_previous),xlast(i,n_previous),yfirst(i,n_previous),ylast(i,n_previous), &
	      center_of_mass_x(i,n_previous),center_of_mass_y(i,n_previous), &
	      xmax(i,n_previous),ymax(i,n_previous), &
	      xmin(i,n_previous),ymin(i,n_previous), &
	      velocity_x(i), &
	      largest_forward_link(i),largest_forward_link_size(i), &
	      second_largest_forward_link(i),second_largest_forward_link_size(i), &
	      largest_backward_link(i),largest_backward_link_size(i), &
	      second_largest_backward_link(i),second_largest_backward_link_size(i)
ENDDO

! identify backward links
largest_backward_link=0
largest_backward_link_size=0
second_largest_backward_link=0
second_largest_backward_link_size=0
DO i=1,counter_actual-1
  a = event_number(first_point_x(i,n_actual),first_point_y(i,n_actual),3)
  d = event_number(first_point_x1(i,n_actual),first_point_y1(i,n_actual),3)
  e = event_number(first_point_x2(i,n_actual),first_point_y2(i,n_actual),3)
  IF (a==0) THEN
    largest_backward_link(i) = -1    ! track interrupted by missing values
    !WRITE(*,*) "value 0:",first_point_x(i,n_actual),first_point_y(i,n_actual)
    CYCLE
  ENDIF
  !IF (a==0) THEN
  !  WRITE(*,*) "value -1:",first_point_x(i,n_actual),first_point_y(i,n_actual)
  !ENDIF
  DO j=1,counter_previous-1
    b = event_number(first_point_x(j,n_previous),first_point_y(j,n_previous),3)
    c = event_number(first_point_x1(j,n_previous),first_point_y1(j,n_previous),3)
    f = event_number(first_point_x2(j,n_previous),first_point_y2(j,n_previous),3)
    VV = 0.
    IF (b==0) CYCLE
    IF (a==b.or.d==c.or.e==f.or.d==f.or.e==c) THEN 
      ! object inherits age from oldest predecessor
      IF (cell_age1(j,n_previous) .GE. cell_age1(i,n_actual)) THEN
        cell_age1(i,n_actual)=cell_age1(j,n_previous)+1 ! aging cell
      ENDIF

      if (center_of_mass_x(i,n_actual).gt.0 .and. center_of_mass_y(i,n_actual).gt.0) then 
      	xx1=center_of_mass_x(i,n_actual)
      	yy1=center_of_mass_y(i,n_actual)
      else 
      	xx1=xmin(i,n_actual)
      	yy1=ymin(i,n_actual)
      endif
      if (center_of_mass_x(j,n_previous).gt.0 .and. center_of_mass_y(j,n_previous).gt.0) then 
      	xx2=center_of_mass_x(j,n_previous)
      	yy2=center_of_mass_y(j,n_previous)
      else 
      	xx2=xmin(j,n_previous)
      	yy2=ymin(j,n_previous)
      endif

      dx = xx2 - xx1
      dy = yy2 - yy1
      vv = sqrt(dx*dx+dy*dy)

      !dx = center_of_mass_x(j,n_previous) - center_of_mass_x(i,n_actual)
      !dy = center_of_mass_y(j,n_previous) - center_of_mass_y(i,n_actual)
      !vv = sqrt(dx*dx+dy*dy)

      IF (vv .LT. max_velocity) THEN

      IF (totarea(j,n_previous) .GT. largest_backward_link_size(i)) THEN
        second_largest_backward_link(i) = largest_backward_link(i)
		second_largest_backward_link_size(i) = largest_backward_link_size(i)
		largest_backward_link(i) = counter_total_previous+j
		largest_backward_link_size(i) = totarea(j,n_previous)
        ! object inherits age from largest predecessor
        cell_age2(i,n_actual)=cell_age2(j,n_previous)+1 ! aging cell
      ELSEIF (totarea(j,n_previous) .GT. second_largest_backward_link_size(i)) THEN
		second_largest_backward_link(i) = counter_total_previous+j
		second_largest_backward_link_size(i) = totarea(j,n_previous)
      ENDIF
      ENDIF

    ENDIF
  ENDDO
  !IF (cell_age1(i,n_actual).LT.cell_age2(i,n_actual)) THEN
  !  WRITE(*,*) "cell age:",cell_age1(i,n_actual),cell_age2(i,n_actual)
  !ENDIF
ENDDO

! flip time indices
n_actual   = 3-n_actual
n_previous = 3-n_previous

! end main loop
ENDDO

200 CONTINUE

! close input/output files
CLOSE(20)

DO fileid=1,n_fields+1
  CLOSE(fileid)
ENDDO


END PROGRAM irt_objects


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                           SUBROUTINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE write_srv(field,date,time,file_id)

  USE irt_parameters, ONLY: domainsize_x, domainsize_y

  IMPLICIT NONE
  !INTEGER, INTENT(IN)   :: domainsize_x,domainsize_y
  REAL, INTENT(IN)      :: field(domainsize_x,domainsize_y)
  INTEGER, INTENT(IN)   :: date,time
  INTEGER, INTENT(IN)   :: file_id
  INTEGER               :: srv_header(8)

  srv_header(1) = 1	      ! Code
  srv_header(2) = 1	      ! Level
  srv_header(3) = date        ! Datum
  srv_header(4) = time        ! Zeitinkrement
  srv_header(5) = domainsize_x
  srv_header(6) = domainsize_y
  srv_header(7) = 0
  srv_header(8) = 0

  WRITE (file_id) srv_header
  WRITE (file_id) field
  
  RETURN

END SUBROUTINE write_srv


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE set_event_number_to_value(domsize_x,domsize_y,value,event_number,nevent, &
                                     xfirst,xlast,yfirst,ylast)

  IMPLICIT NONE
  INTEGER, INTENT(IN)	    :: domsize_x,domsize_y
  INTEGER, INTENT(IN)       :: value
  INTEGER, INTENT(INOUT)    :: event_number(domsize_x,domsize_y)
  INTEGER, INTENT(IN)       :: nevent
  INTEGER, INTENT(IN)       :: xfirst,xlast,yfirst,ylast
  INTEGER                   :: ix,iy,ix_mod,iy_mod
  
  DO iy=yfirst,ylast
    DO ix=xfirst,xlast
      ix_mod = MOD(ix-1+domsize_x,domsize_x)+1
      iy_mod = MOD(iy-1+domsize_y,domsize_y)+1
      IF (event_number(ix_mod,iy_mod)==nevent) event_number(ix_mod,iy_mod)=value
    ENDDO
  ENDDO
  
  RETURN

END SUBROUTINE set_event_number_to_value

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!! deleted SUBROUTINE decrease_resolution(domsize_x,domsize_y,miss,field)

!!!!!!!changed+++++++++++++++++++++++++++++++++++++++++++++++++++++++

RECURSIVE SUBROUTINE area_consecuprcp(ii,ij,domsize_x,domsize_y,input_field, &
                                      occupied_consecuprcp,totarea3,prcp_mean3,delete_cell)
  USE irt_parameters, ONLY: n_fields, threshold,threshold2, miss, llonlatgrid, &
                            unit_area, lat_first, lat_inc, lon_inc
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: ii, ij
  INTEGER, INTENT(IN)    :: domsize_x,domsize_y
  REAL, INTENT(IN)       :: input_field(domsize_x,domsize_y,n_fields+1)
  LOGICAL, INTENT(INOUT) :: occupied_consecuprcp(domsize_x,domsize_y)
  REAL, INTENT(INOUT)    :: totarea3
  REAL, INTENT(INOUT)    :: prcp_mean3
  INTEGER                :: i, ii_mod, ij_mod
  INTEGER                :: icell(4), jcell(4)
  LOGICAL, INTENT(INOUT) :: delete_cell

  REAL, PARAMETER        :: deg2rad = 3.141592654/180.
  REAL                   :: gridboxarea

  ! Indices for all 4 flow directions
  icell = (/ ii+1, ii  , ii-1, ii  /)
  jcell = (/ ij  , ij+1, ij  , ij-1 /)
  
  ii_mod = MOD(ii-1+domsize_x,domsize_x)+1
  ij_mod = MOD(ij-1+domsize_y,domsize_y)+1

  IF (.NOT. occupied_consecuprcp(ii_mod,ij_mod)) THEN
     occupied_consecuprcp(ii_mod,ij_mod) = .TRUE.

     IF (input_field(ii_mod,ij_mod,1) .LT. threshold .AND. input_field(ii_mod,ij_mod,1) .GT.(miss+100) &
        .AND. input_field(ii_mod,ij_mod,2) .GE. threshold2 ) THEN
        ! center of mass is now weighted by intensity!!!
        ! area in km^2
        gridboxarea = lon_inc*lat_inc*unit_area*COS((lat_first+ij_mod*lat_inc)*deg2rad)

        totarea3 = totarea3 + gridboxarea
        prcp_mean3 = prcp_mean3 + input_field(ii_mod,ij_mod,2)*gridboxarea

        DO i=1,4
           CALL area_consecuprcp(icell(i),jcell(i),domsize_x,domsize_y,input_field,&
                    occupied_consecuprcp,totarea3,prcp_mean3, delete_cell)
        ENDDO

     ELSEIF (input_field(ii_mod,ij_mod,1) .LE. miss) THEN
        ! if cell touches missing value, this cell will be deleted
        delete_cell = .TRUE.
     ENDIF
  ENDIF

  RETURN


END SUBROUTINE area_consecuprcp


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

RECURSIVE SUBROUTINE area(ii,ij,domsize_x,domsize_y,input_field, &
      occupied,nevent,event_number,totarea,totarea1,totarea2,prcp_mean1,prcp_mean2,&
      field_mean,field_min,field_max, &
		  COMx,COMy,delete_cell,xfirst,xlast,yfirst,ylast,xmax,ymax,pmax,xmin,ymin,tbmin)

  USE irt_parameters, ONLY: n_fields, threshold,threshold2, miss, llonlatgrid, &
                            unit_area, lat_first, lat_inc, lon_inc

  IMPLICIT NONE

  INTEGER, INTENT(IN)	 :: ii, ij
  INTEGER, INTENT(IN)	 :: domsize_x,domsize_y
  !REAL, INTENT(IN)	 :: miss
  !INTEGER, INTENT(IN)	 :: n_fields
  REAL, INTENT(IN)       :: input_field(domsize_x,domsize_y,n_fields+1)
  INTEGER,INTENT(INOUT)  :: event_number(domsize_x,domsize_y)
  LOGICAL, INTENT(INOUT) :: occupied(domsize_x,domsize_y)
  INTEGER, INTENT(INOUT) :: nevent
  REAL, INTENT(INOUT)    :: field_mean(n_fields+1),field_min(n_fields+1),field_max(n_fields+1)
  REAL, INTENT(INOUT)    :: prcp_mean1,prcp_mean2,pmax,tbmin
  REAL, INTENT(INOUT)    :: totarea,totarea1
  REAL, INTENT(INOUT)    :: totarea2  !!area of prcp>= 2mm/h
  REAL, INTENT(INOUT)	   :: COMx,COMy !!center of mass of precipitation
  INTEGER, INTENT(INOUT) :: xfirst,xlast,yfirst,ylast,xmax,ymax,xmin,ymin
  INTEGER		 :: i, ii_mod, ij_mod,fieldid
  INTEGER		 :: icell(4), jcell(4)
  !REAL, INTENT(IN)       :: threshold
  LOGICAL, INTENT(INOUT) :: delete_cell

  REAL, PARAMETER        :: deg2rad = 3.141592654/180.
  REAL                   :: gridboxarea
  
  ! Indices for all 4 flow directions
  icell = (/ ii+1, ii  , ii-1, ii  /)
  jcell = (/ ij  , ij+1, ij  , ij-1 /)
  
  ii_mod = MOD(ii-1+domsize_x,domsize_x)+1
  ij_mod = MOD(ij-1+domsize_y,domsize_y)+1
  

  IF (.NOT. occupied(ii_mod,ij_mod)) THEN
     ! take care for periodic boundary conditions:
     IF (ii .LT. xfirst) xfirst=ii
     IF (ii .GT. xlast) xlast=ii
     IF (ij .LT. yfirst) yfirst=ij
     IF (ij .GT. ylast) ylast=ij
     occupied(ii_mod,ij_mod) = .TRUE.

     IF (input_field(ii_mod,ij_mod,1) .LT. threshold .AND. input_field(ii_mod,ij_mod,1) .GT.(miss+100) ) THEN
	      ! center of mass is now weighted by intensity!!!
        ! area in km^2
        gridboxarea = lon_inc*lat_inc*unit_area*COS((lat_first+ij_mod*lat_inc)*deg2rad)
        totarea = totarea + gridboxarea

        IF (input_field(ii_mod,ij_mod,1) .LT. tbmin) THEN 
        	xmin = ii_mod
            ymin = ij_mod
            tbmin= input_field(ii_mod,ij_mod,1)
        ENDIF

        !!!!changed
        IF (input_field(ii_mod,ij_mod,2) .GT. 0.1 ) THEN    
          totarea1 = totarea1 + gridboxarea
          prcp_mean1 = prcp_mean1 + input_field(ii_mod,ij_mod,2)*gridboxarea

          COMx = COMx + ii*input_field(ii_mod,ij_mod,2)*gridboxarea
          COMy = COMy + ij*input_field(ii_mod,ij_mod,2)*gridboxarea  

          IF (input_field(ii_mod,ij_mod,2) .GT. pmax ) THEN 
            xmax = ii_mod
            ymax = ij_mod
            pmax = input_field(ii_mod,ij_mod,2)
          ENDIF

        IF (input_field(ii_mod,ij_mod,2) .GE. threshold2 ) THEN   !!! threshold2 = 2mm/h 
          totarea2 = totarea2 + gridboxarea
          prcp_mean2 = prcp_mean2 + input_field(ii_mod,ij_mod,2)*gridboxarea
        ENDIF

        ENDIF


        event_number(ii_mod,ij_mod) = nevent
        !IF (gridboxarea .LT. 0.5) WRITE (*,*) "gridboxarea=",gridboxarea
	      DO fieldid=1,n_fields+1
	        field_mean(fieldid) = field_mean(fieldid) + input_field(ii_mod,ij_mod,fieldid)*gridboxarea
	        field_min(fieldid)  = MIN(field_min(fieldid),input_field(ii_mod,ij_mod,fieldid))
	        field_max(fieldid)  = MAX(field_max(fieldid),input_field(ii_mod,ij_mod,fieldid))
	      ENDDO


	    
        DO i=1,4
           CALL area(icell(i),jcell(i),domsize_x,domsize_y,input_field,&
	              occupied,nevent,event_number,totarea,totarea1,totarea2,prcp_mean1,prcp_mean2,&
                field_mean,field_min,field_max, &
		            COMx,COMy,delete_cell,xfirst,xlast,yfirst,ylast,xmax,ymax,pmax,xmin,ymin,tbmin)
        ENDDO

      ELSEIF (input_field(ii_mod,ij_mod,1) .LE. miss) THEN
        ! if cell touches missing value, this cell will be deleted
        delete_cell = .TRUE.
      ENDIF
  ENDIF

  RETURN

END SUBROUTINE area

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

RECURSIVE SUBROUTINE overlay_area(miss, ii,ij,domsize_x,domsize_y,input_field,occupied,nevent, &
                                  event_number,delete_cell,xfirst,xlast,yfirst,ylast)
  IMPLICIT NONE
  INTEGER, INTENT(IN)	 :: ii, ij
  INTEGER, INTENT(IN)	 :: domsize_x,domsize_y
  REAL, INTENT(IN)	 :: miss
  REAL, INTENT(IN)       :: input_field(domsize_x,domsize_y)
  INTEGER,INTENT(INOUT)  :: event_number(domsize_x,domsize_y)
  LOGICAL, INTENT(INOUT) :: occupied(domsize_x,domsize_y)
  INTEGER, INTENT(INOUT) :: nevent
  INTEGER, INTENT(INOUT) :: xfirst,xlast,yfirst,ylast
  INTEGER		 :: i, ii_mod, ij_mod
  INTEGER		 :: icell(4), jcell(4)
  LOGICAL, INTENT(INOUT) :: delete_cell
  
  ! Indices for all 4 flow directions
  icell = (/ ii+1, ii  , ii-1, ii  /)
  jcell = (/ ij  , ij+1, ij  , ij-1 /)
  
  ii_mod = MOD(ii-1+domsize_x,domsize_x)+1
  ij_mod = MOD(ij-1+domsize_y,domsize_y)+1
  
  IF (.NOT. occupied(ii_mod,ij_mod)) THEN
     ! take care for periodic boundary conditions:
     IF (ii .LT. xfirst) xfirst=ii
     IF (ii .GT. xlast) xlast=ii
     IF (ij .LT. yfirst) yfirst=ij
     IF (ij .GT. ylast) ylast=ij
     occupied(ii_mod,ij_mod) = .TRUE.
     IF (input_field(ii_mod,ij_mod) .GE. 0.5) THEN
	! center of mass is now weighted by intensity!!!
        event_number(ii_mod,ij_mod) = nevent
		DO i=1,4
           CALL overlay_area(miss, icell(i),jcell(i),domsize_x,domsize_y,input_field,occupied,nevent, &
	             event_number,delete_cell,xfirst,xlast,yfirst,ylast)
        ENDDO
     ! Bugfix! Do NOT delete boundary cells in the overlaid field!!!
     !ELSEIF (input_field(ii_mod,ij_mod) .LE. miss) THEN
     !   ! if cell touches missing value, this cell will be deleted
     !   delete_cell = .TRUE.
     ENDIF
  ENDIF

  RETURN

!END SUBROUTINE overlay_area

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! linear time interpolation of velocity field

! deleted

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! bilinear spatial interpolation of coarse velocity field

! deleted

END SUBROUTINE
