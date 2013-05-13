function read_mini, run, camcol, ra, dec, rad, type, rerun=rerun
    ;dir = '/scr/apache11/users/cpl/2009-11-16.v3/'
    ; miniSweeps-008100-star.fits.gz
    ;file = string(dir, 'miniSweeps-', run, '-', type, '.fits.gz', $
    ;              format='(a,a,i06,a,a0,a)')
    ;splog, "reading ", file

    objs = sweep_readobj(run, camcol, type=type, rerun=rerun)
    ;objs = mrdfits(file, 1)
    splog, "for run-camcol, ", run, camcol, " read ", n_elements(objs), " objects"

    ; mag limits -- require S/N of 10 in r&i
    w = where((objs.modelflux[2]*sqrt(objs.modelflux_ivar[2]) GE 10) AND $
              (objs.modelflux[3]*sqrt(objs.modelflux_ivar[3]) GE 10), cnt)
    if (cnt eq 0) then begin
       splog, "!!!!! no brighter objects in ", run, "-", camcol
       return, 0
    end
    objs = objs[w]
    splog, "    ", n_elements(objs), " objects with S/N >= 10 for ", run, "-", camcol


    dudmask = sdss_flagval('OBJECT1', 'BRIGHT') OR $ 
              sdss_flagval('OBJECT1', 'EDGE') OR $
              sdss_flagval('OBJECT1', 'BLENDED')
    dudmask2 = sdss_flagval('OBJECT2', 'INTERP_CENTER')
    ok = where(((objs.objc_flags AND dudmask) EQ 0) AND ((objs.objc_flags2 AND dudmask2) EQ 0), cnt)
    if (cnt eq 0) then begin
       splog, "!!!!! no passable objects in ", run, "-", camcol
       return, 0
    end
    objs = objs[ok]
    splog, "     kept ", n_elements(objs), " objects"

    spherematch, objs.ra, objs.dec, ra, dec, rad, m1, m2, d12, maxmatch=0, $
                 estnmatch=n_elements(objs)
    splog, "matched ", n_elements(m1), " objects"
    if (n_elements(m1) eq 1) && (m1 eq -1) then $
       return,0
    objs = objs[m1]

                                ; Remake the table:
                                ;  flux -> mags
                                ;  mags[5] -> u,g,r,i,z
    astrom_struct = {thing_id:objs[0].thing_id*0, $
                     ra:objs[0].ra*0.0, dec:objs[0].dec*0.0, $
                     objc_flags:objs[0].objc_flags*0, $
                     objc_flags2:objs[0].objc_flags2*0, $
                     modelmag_u:0.0, modelmag_g:0.0, modelmag_r:0.0, modelmag_i:0.0, modelmag_z:0.0, $
                     modelmag_ivar_u:0.0, modelmag_ivar_g:0.0, modelmag_ivar_r:0.0, $
                     modelmag_ivar_i:0.0, modelmag_ivar_z:0.0, $
                     objc_type:objs[0].objc_type*0, $
                     object_id:ulong64(objs[0].thing_id)*0ULL $
                    }

    astrom = replicate(astrom_struct, n_elements(m1))
    struct_assign, objs, astrom

    ; Make a scrutable synthetic object ID
    ; RRRR C FFFF OOOOO
    astrom.object_id = ((ulong64(objs.RUN)*10ULL + ulong64(objs.CAMCOL))*10000ULL + ulong64(objs.FIELD))*100000ULL + ulong64(objs.ID)

    astrom.modelmag_u = 22.5 - 2.5*alog10(objs.modelflux[0])
    astrom.modelmag_g = 22.5 - 2.5*alog10(objs.modelflux[1])
    astrom.modelmag_r = 22.5 - 2.5*alog10(objs.modelflux[2])
    astrom.modelmag_i = 22.5 - 2.5*alog10(objs.modelflux[3])
    astrom.modelmag_z = 22.5 - 2.5*alog10(objs.modelflux[4])
    
    astrom.modelmag_ivar_u = (alog(10.)/2.5)^2 * (objs.modelflux[0] > 0)^2 * objs.modelflux_ivar[0]
    astrom.modelmag_ivar_g = (alog(10.)/2.5)^2 * (objs.modelflux[1] > 0)^2 * objs.modelflux_ivar[1]
    astrom.modelmag_ivar_r = (alog(10.)/2.5)^2 * (objs.modelflux[2] > 0)^2 * objs.modelflux_ivar[2]
    astrom.modelmag_ivar_i = (alog(10.)/2.5)^2 * (objs.modelflux[3] > 0)^2 * objs.modelflux_ivar[3]
    astrom.modelmag_ivar_z = (alog(10.)/2.5)^2 * (objs.modelflux[4] > 0)^2 * objs.modelflux_ivar[4]

    return, astrom
 end

    
