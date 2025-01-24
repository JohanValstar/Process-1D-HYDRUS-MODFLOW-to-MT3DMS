!  postproces_hydrus_modflow2mt3dms.f90 
!
!  FUNCTIONS:
!  postproces_hydrus_modflow2mt3dms - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: postproces_hydrus_modflow2mt3dms
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program postproces_hydrus_modflow2mt3dms

    implicit none
    ! Variables
    character*500 filename, filename_zonation, filename_solute_details, filename_out,text,root_out,ftlfile_in,ftlfile_out,  filename_flux
    character*250 line,lineprev,lineprev_store
    character*50 zoneformat,format_out, timestepchar
    integer ncol,nrow,i,icol,irow,first,last
    integer, allocatable:: zonation(:,:)
    integer timestep,nprofile,profile_nr,ntimestep
    real time,timeprev,conc
    real, allocatable:: conc_wt(:),times(:)
    real flux(2)
    double precision vbot
    real time2
    character*15 char1
    character*42 char2,charvbot
    integer zone
    CHARACTER VERSION*11
    integer MTWEL,MTDRN,MTRCH,MTEVT,MTRIV,MTGHB,MTCHD,ISS,NPERFL,MTSTR,MTRES,MTFHB,MTDRT,MTETS,MTTLK,MTIBS,MTLAK,MTMNW,MTSWT,MTSFR,MTUZF
    real cell_area
    logical noquotes,lastloc
    character*3 nametype(100),before_this_type
    
    

    ! Body of postproces_hydrus_modflow2mt3dms
    
    ! get argument file name
    call GET_COMMAND(text)
    noquotes=.true.
    do i=1,500
        if ( text(i:i).eq.'"' .or. text(i:i).eq."'" )then
            if (noquotes)then
                 noquotes=.false.
            else
                noquotes=.true.             
            end if    
        end if     
        if ( text(i:i).eq.' '.and.noquotes)then
            first = i 
            exit
        end if   
    end do 
    write(filename,*)text(i:500)    
    filename = adjustl(filename)

    ! read argument file
    open(unit=1,file = filename,status='old')
    read (1,*) filename_zonation
    read(1,*) filename_solute_details
    read(1,*) filename_flux
    read(1,*) ncol,nrow
    read(1,*)root_out
    read(1,*)format_out
    read(1,*)ftlfile_in
    read(1,*)ftlfile_out
    read(1,*) cell_area
    close(1)
    
    ! allocate and read zonation array    
    allocate(zonation(ncol,nrow))
    open(unit=1,file=filename_zonation, status='old')
    read(1,*)
    read(1,*)
    read(1,'(a500)')text
    text=adjustl(text) 
    if (text(1:8).eq.'INTERNAL'.or.text(1:8).eq.'internal' ) then
        ! read format of data
        do i=9,500
            if ( text(i:i).eq.'(') first = i
            if ( text(i:i).eq.')')then
               last = i 
               exit
            end if
        end do    
        write(zoneformat,*)text(first:last)    
        do irow=1,nrow
            read(1,zoneformat) (zonation(icol,irow),icol=1,ncol)
        end do    
    else
        write(*,*)'third line not starting with "INTERNAL" or internal" in file ', filename_zonation  
        write(*,*)'stopping'
        stop
    end if   
    close(1)

    
    ! find number of timesteps
    
    open (unit=1,file=filename_solute_details,status='old')
    read(1,*)
    read(1,*)
    ntimestep=0
5  read(1,*,end=9) profile_nr,time,conc
    if (profile_nr.eq.1)then
        ntimestep=ntimestep+1
    end if
    goto 5
9   close(1)

    
    
    
    
    !allocate and initialize array
    nprofile = maxval(zonation)
    allocate(conc_wt(nprofile),times(ntimestep))
    conc_wt = 0.0  
    
    !read from concentration data 
    open (unit=1,file=filename_solute_details,status='old')
    read(1,*)
    read(1,*)
    timeprev=0
    timestep=0
10  read(1,*,end=99) profile_nr,time,conc
    if (profile_nr.eq.1)then
        timestep=timestep+1
        times(timestep)=time
         
    end if    
    if (profile_nr.eq.1.and.timeprev.gt.0)then
        ! all relevant concentrations of previous time step are known, process them
        write(timestepchar, *)timestep-1
        timestepchar = adjustl(timestepchar)
        filename_out = trim(root_out)//timestepchar
            
        open(unit=2,file=filename_out)
        do irow=1,nrow
            write(2,format_out) (conc_wt(zonation(icol,irow)),icol=1,ncol)
        end do    
        close(2)    
    end if
    
    conc_wt(profile_nr)= conc
    timeprev = time
    goto 10
           
99  write(timestepchar, *)timestep
    ! and also process for the last time step
    timestepchar = adjustl(timestepchar)
    filename_out = trim(root_out)//timestepchar
            
    open(unit=2,file=filename_out)
    do irow=1,nrow
        write(2,format_out) (conc_wt(zonation(icol,irow)),icol=1,ncol)
    end do    
    close(2)    
    close(1)     
    

    ! Adapt the FTL flux file  
    open(unit=1,file=ftlfile_in)
    ! first find location  where RCH flux files should be added 
12  read(1,'(a80)',end=200)line
    i=1
    if(line(2:2).eq."'")then
        write(nametype(i),*)  line(3:5)
       if (line(3:5).eq.'EVT'.or.line(3:5).eq.'ETS'.or.line(3:5).eq.'RIV'.or.line(3:5).eq.'GHB'.or.&
           line(3:5).eq.'STR'.or.line(3:5).eq.'RES'.or.line(3:5).eq.'FHB'.or.line(3:5).eq.'MNW'.or.&
           line(3:5).eq.'DRT')then
           ! if fluxes for one of these modules are present, recharge flux should be added before this 
           before_this_type = line(3:5)
           lastloc = .false. 
           goto 15 
       end if    
       if(line(3:5).eq.nametype(1))then
           ! if fluxes for one of these modules are not present, recharge flux should be added after all the fluxes for each time step 
           before_this_type = nametype(1)
           lastloc = .true. 
           goto 15
       end if
       i=i+1 
    end if
    goto 12  
15  close(1)
    
    
    ! and now insert the fluxes on the right locations
    open(unit=1,file=ftlfile_in)
    open(unit=2,file=ftlfile_out)
    open(unit=3,file= filename_flux)
    
    
    ! read and copy some header
    read(1,*) VERSION,MTWEL,MTDRN,MTRCH,MTEVT,MTRIV,MTGHB,MTCHD,ISS,NPERFL, &
       MTSTR,MTRES,MTFHB,MTDRT,MTETS,MTTLK,MTIBS,MTLAK,MTMNW,MTSWT,MTSFR,MTUZF
    MTRCH = ncol*nrow
    write(2,*)VERSION,MTWEL,MTDRN,MTRCH,MTEVT,MTRIV,MTGHB,MTCHD,ISS,NPERFL,&
         MTSTR,MTRES,MTFHB,MTDRT,MTETS,MTTLK,MTIBS,MTLAK,MTMNW,MTSWT,MTSFR,MTUZF
    
    timestep = 1
20  if(line(2:2).eq."'")then
        ! one array header line is read before the type of flux is known, store it to reuse it later
        lineprev_store = lineprev
    end if
    lineprev= line
    read(1,'(a80)',end=200)line
    if(line(3:5).eq.before_this_type) then
        call insert_recharge_flux(ncol,nrow,zonation,cell_area,timestep,times,ntimestep)
        write(2,'(a80)')lineprev_store 
    end if
    write(2,'(a80)')line 
    goto 20
    
200 close(1)
    if (lastloc) then
        write(2,'(a80)')lineprev_store  
        call insert_recharge_flux(ncol,nrow,zonation,cell_area,timestep,times,ntimestep)
    end if    
    close(2)
    end program postproces_hydrus_modflow2mt3dms

    
!   *********************************************************************************************************************    
    subroutine insert_recharge_flux(ncol,nrow,zonation,cell_area,timestep,times,ntimestep)
    implicit none
    character*250 line,lineprev
    integer ncol,nrow,icol,irow
    integer zonation(ncol,nrow)
    real flux(2)
    double precision vbot
    real time2,dt
    character*15 char1
    character*42 char2
    real cell_area
    integer timestep,ntimestep
    integer zone
    real times(ntimestep)
    
    zone=1
30  read(3,'(a15,E14.8, a43, D15.10)')  char1,time2,char2,vbot
    if (abs(time2-times(timestep)).le.1.0e-4)then
        flux(zone)=vbot
        zone=zone+1
        if (zone.eq.3)then
            timestep = timestep + 1
            zone = 1
            goto 80
        end if    
    end if 
    goto 30
      
80  write(2,'(a20, I6 )') " 'RCH             ' ",ncol*nrow
    do irow=1,nrow
        ! write  layer array of recharge layer  (always 1 )  
        write(2,'(28a2)')(' 1', icol=1,ncol)
    end do
    !write flux area, correct for area and minus sign 
    do irow=1,nrow  
          write(2,'(28E16.6 )')(-cell_area*flux(zonation(icol,irow)),icol=1,ncol)     
    end do
    return
    end subroutine