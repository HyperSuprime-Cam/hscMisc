
;setenv, 'PHOTO_RESOLVE=/u/dss/resolve/2009-11-16'
;setenv, 'PHOTO_SWEEP=/u/dss/sweeps/2009-11-16.v3'
pro z
    setenv, 'PHOTO_RESOLVE=/u/dss/resolve/2010-05-23'
    setenv, 'PHOTO_SWEEP=/u/dss/sweeps/dr8_final'
    
    ra = 15*[ten(10,00,29), ten(02,18,0), ten(23,29,08.280), ten(16,30,00.000), ten(14,17,54),    ten(10,57,0), ten(20,51,27), ten(0,22,4)]
    dec =   [ten(02,12,21), ten(-5,0,0),  ten(-03,01,58.80), ten(0,0,0),        ten(52,30,31.00), ten(57,40,0), ten(0,55,30), ten(0,-38,11)]
    rad =   [3.0, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.0, 1.0]
    name =  ["COSMOS", "SXDS", "CFHQS", "GalPlane", "CFHTLS-D3", "LockmanHole", "HM1", "HM2"]
    
    rerun = 301
    objtypes = ['star','gal']
    
    allobjs = 0
    for o_i=0,n_elements(objtypes)-1 do begin
       alltypeobjs = 0
       objtype = objtypes[o_i]
       
       for pointing=0,n_elements(ra)-1 do begin
          flist = flist_sub(ra[pointing], dec[pointing], rad[pointing])
          runids = flist.run
          ii = uniq(runids, sort(runids))
          print, "ra, dec: ", ra[pointing], dec[pointing]
          print, "runs: ", runids[ii]
          
          objlist = 0
          for i=0,n_elements(ii)-1 do begin
             run = runids[ii[i]]
             for camcol=1,6 do begin
                runobjs = read_mini(run, camcol, ra[pointing], dec[pointing], rad[pointing], objtype, rerun=rerun)
                if keyword_set(runobjs) then begin
                   objlist = keyword_set(objlist) ? [objlist, runobjs] : runobjs 
                   alltypeobjs = keyword_set(alltypeobjs) ? [alltypeobjs, runobjs] : runobjs 
                endif else begin
                   splog, "FAILED to match any objects in run ", run
                endelse
             endfor
             splog, "after ", i, " of ", ii, " runs: ", n_elements(objlist) 
          endfor
          
          fxbhmake, hdr, 0
          sxaddpar, hdr, 'PROJECT', "HSC DC2 pointings"
          sxaddpar, hdr, 'PTG_NAME', name[pointing]
          sxaddpar, hdr, 'PTG_RA',   ra[pointing]
          sxaddpar, hdr, 'PTG_DEC',  dec[pointing]
          sxaddpar, hdr, 'PTG_RAD',  rad[pointing]
          for r_i=0,n_elements(ii)-1 do begin
             sxaddpar, hdr, string('BSSRUN', r_i+1, format='(a,i0)'), runids[ii[r_i]], 'The BOSS runs we draw from'
          endfor
          
          mwrfits, objlist, string('miniSweeps-hsc_dc2-', pointing+1, '-', objtype, '.fits', $
                                   format='(a,i02,a,a0,a)'), hdr, /create
       endfor
       mwrfits, alltypeobjs, string('miniSweeps-hsc_dc2-', objtype, '.fits', $
                                    format='(a,a0,a)'), /create
       allobjs = keyword_set(allobjs) ? [allobjs, alltypeobjs] : alltypeobjs 
    endfor
    
    mwrfits, allobjs, 'miniSweeps-hsc_dc2-all.fits'
 end
