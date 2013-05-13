pro mash_sweeps, run, objtype=objtype, rerun=rerun, outdir=outdir, $
                 ra=ra, dec=dec, matchrad=matchrad, $
                 recalibrate=recalibrate, $
                 astrom=astrom, logical_cols=logical_cols, colnames=colnames

                                ; optionally set up to recalibrate
    if keyword_set(recalibrate) then begin
       setenv,"PHOTO_CALIB="+recalibrate
       splog, 'recalibrating using PHOTO_CALIB=', getenv('PHOTO_CALIB')
    endif

    if not keyword_set(objtype) then begin
       stars = 0
       mash_sweeps, run, objtype='star', rerun=rerun, ra=ra, dec=dec, matchrad=matchrad, $
                    astrom=stars, logical_cols=logical_cols, $
                    recalibrate=recalibrate
       gals = 0
       mash_sweeps, run, objtype='gal', rerun=rerun, ra=ra, dec=dec, matchrad=matchrad, $
                    astrom=gals, logical_cols=logical_cols, $
                    recalibrate=recalibrate
       if keyword_set(stars) then begin
          if keyword_set(gals) then astrom = [stars, gals] else astrom = stars
       end else begin
          astrom = gals
       end

       if (not keyword_set(astrom)) then begin
          splog, "!!!!! no objects at all in run ", run
          return
       end
    end else begin
       allcols = 0
       for camcol=1,6 do begin
          if keyword_set(objtype) then begin
             objs = sweep_readobj(run, camcol, type=objtype, rerun=rerun)
          end else begin
             objs = [sweep_readobj(run, camcol, type='star', rerun=rerun), $
                     sweep_readobj(run, camcol, type='gal', rerun=rerun)] 
          end
          
          if not keyword_set(objs) then begin
             splog, "!!!!! no objects in ", run, "-", camcol
             continue
          end
          splog, "read ", n_elements(objs), " objects for ", run, "-", camcol
          
                                ; trim to survey primary
          w = where(objs.resolve_status AND sdss_flagval('RESOLVE_STATUS', 'SURVEY_PRIMARY'), cnt)
          if (cnt eq 0) then begin
             splog, "!!!!! no primary objects in ", run, "-", camcol
             continue
          end
          
          objs = objs[w]
          splog, "    ", n_elements(objs), " primary objects for ", run, "-", camcol
          
                                ; Optionally filter by position
          if keyword_set(matchrad) then begin
             spherematch, objs.ra, objs.dec, ra, dec, matchrad, m1, m2, d12, $
                          maxmatch=0, estnmatch=n_elements(objs)
             if (n_elements(m1) EQ 1) && (m1 eq -1) then begin
                splog, "!! no objects within ", matchrad, " deg of ", ra, ", ", dec
                continue
             end
             objs = objs[m1]
             splog, "     really kept ", n_elements(objs), " objects"
          end

                                ; optionally recalibrate
          if keyword_set(recalibrate) then begin
             objs0 = objs
             sdss_recalibrate, objs
          endif

                                ; mag limits -- require S/N of 10 in r&i
          w = where((objs.modelflux[2]*sqrt(objs.modelflux_ivar[2]) GE 10) AND $
                    (objs.modelflux[3]*sqrt(objs.modelflux_ivar[3]) GE 10), cnt)
          if (cnt eq 0) then begin
             splog, "!!!!! no brighter objects in ", run, "-", camcol
             continue
          end
          objs = objs[w]
          splog, "    ", n_elements(objs), " objects with S/N >= 10 for ", run, "-", camcol
          
                                ; 
          dudmask = sdss_flagval('OBJECT1', 'BRIGHT') OR $ 
                    sdss_flagval('OBJECT1', 'EDGE') OR $
                    sdss_flagval('OBJECT1', 'BLENDED')
          dudmask2 = sdss_flagval('OBJECT2', 'INTERP_CENTER')
          ok = where(((objs.objc_flags AND dudmask) EQ 0) AND ((objs.objc_flags2 AND dudmask2) EQ 0), cnt)
          if (cnt eq 0) then begin
             splog, "!!!!! no passable objects in ", run, "-", camcol
             continue
          end
          objs = objs[ok]
          splog, "     kept ", n_elements(objs), " objects"
          
          
          allcols = keyword_set(allcols) ? [allcols, objs] : objs
          splog, "     ", n_elements(allcols), " objects for the run"
       end
       
       if (not keyword_set(allcols)) then begin
          splog, "!!!!! no objects at all in run ", run
          return
       end
       
                                ; Remake the table:
                                ;  flux -> mags
                                ;  mags[5] -> u,g,r,i,z
       astrom_struct = {id:long64(objs[0].thing_id)*0LL, $
                        ra:objs[0].ra*0.0, dec:objs[0].dec*0.0, $
                        u:0.0, g:0.0, r:0.0, i:0.0, z:0.0, $
                        u_err:0.0, g_err:0.0, r_err:0.0, i_err:0.0, z_err:0.0, $
                        objc_flags:objs[0].objc_flags*0, $
                        objc_flags2:objs[0].objc_flags2*0, $
                        calib_status:objs[0].calib_status*0, $
                        thing_id:objs[0].thing_id*0, $
                        starnotgal:0B $
                       }
       colnames = [['id', 'id'], $
                   ['ra', 'ra'], ['dec', 'dec'], $
                   ['u', 'u'], ['g', 'g'], ['r', 'r'], ['i', 'i'], ['z', 'z'], $
                   ['u_err', 'u_err'], ['g_err', 'g_err'], ['r_err', 'r_err'], $
                   ['i_err', 'i_err'], ['z_err', 'z_err'], $
                   ['objc_flags', 'objc_flags'], ['objc_flags2', 'objc_flags2'], $
                   ['calib_status', 'calib_status'], $
                   ['thing_id', 'thing_id'], ['starnotgal', 'starnotgal']]

       astrom = replicate(astrom_struct, n_elements(allcols))
       struct_assign, allcols, astrom
       
                                ; Make a scrutable synthetic object ID
                                ; RRRR C FFFF OOOOO
       astrom.id = ((ulong64(allcols.RUN)*10LL + ulong64(allcols.CAMCOL))*10000LL + ulong64(allcols.FIELD))*100000LL + ulong64(allcols.ID)
       
       astrom.u = 22.5 - 2.5*alog10(allcols.modelflux[0])
       astrom.g = 22.5 - 2.5*alog10(allcols.modelflux[1])
       astrom.r = 22.5 - 2.5*alog10(allcols.modelflux[2])
       astrom.i = 22.5 - 2.5*alog10(allcols.modelflux[3])
       astrom.z = 22.5 - 2.5*alog10(allcols.modelflux[4])
       
       astrom.u_err = 1/sqrt((alog(10.)/2.5)^2 * (allcols.modelflux[0] > 0)^2 * allcols.modelflux_ivar[0])
       astrom.g_err = 1/sqrt((alog(10.)/2.5)^2 * (allcols.modelflux[1] > 0)^2 * allcols.modelflux_ivar[1])
       astrom.r_err = 1/sqrt((alog(10.)/2.5)^2 * (allcols.modelflux[2] > 0)^2 * allcols.modelflux_ivar[2])
       astrom.i_err = 1/sqrt((alog(10.)/2.5)^2 * (allcols.modelflux[3] > 0)^2 * allcols.modelflux_ivar[3])
       astrom.z_err = 1/sqrt((alog(10.)/2.5)^2 * (allcols.modelflux[4] > 0)^2 * allcols.modelflux_ivar[4])
       
       astrom.starnotgal = allcols.objc_type EQ 6
    end

    if keyword_set(outdir) then begin
       fxbhmake, hdr, 0
       sxaddpar, hdr, 'SWEEPDIR', getenv('PHOTO_SWEEP')
       sxaddpar, hdr, 'RESLVDIR', getenv('PHOTO_RESOLVE')
                                ; optionally recalibrate
       if keyword_set(recalibrate) then begin
          sxaddpar, hdr, 'CALIBDIR', getenv('PHOTO_CALIB')
       endif
       
       sxaddpar, hdr, 'RERUN', rerun
       typeStr = keyword_set(objtype) ? string("-", objtype) : ""
       mwrfits, astrom, string(outdir, '/astromSweeps-', run, typeStr, '.fits', $
                               format='(a,a,i04,a0,a)'), hdr, alias=colnames, /create
    endif
 end

    
