pro z4
    setenv, 'PHOTO_RESOLVE=/u/dss/resolve/2010-05-23'
    setenv, 'PHOTO_SWEEP=/u/dss/sweeps/dr8_final'
    rerun = 301
    outdir = '/scr/apache12/users/cpl/HSC/astrom/r8.3'

    ra = 15*[ten(10,00,29), ten(02,18,0), ten(23,29,08.280), ten(16,30,00.000), ten(14,17,54),    ten(10,57,0)]
    dec =   [ten(02,12,21), ten(-5,0,0),  ten(-03,01,58.80), ten(0,0,0),        ten(52,30,31.00), ten(57,40,0)]
    rad =   [3.0, 1.5, 1.5, 1.5, 1.5, 1.5]
    name =  ["COSMOS", "SXDS", "CFHQS", "GalPlane", "CFHTLS-D3", "LockmanHole"]

    allobjs = 0
    objtypes = ['star','gal']
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
             runobjs = 0
             mash_sweeps, run, objtype=objtype, rerun=rerun, $
                          ra=ra[pointing], dec=dec[pointing], matchrad=rad[pointing], $
                          astrom=runobjs, colnames=colnames
             if keyword_set(runobjs) then begin
                objlist = keyword_set(objlist) ? [objlist, runobjs] : runobjs 
                alltypeobjs = keyword_set(alltypeobjs) ? [alltypeobjs, runobjs] : runobjs 
             endif else begin
                splog, "FAILED to match any objects in run ", run
             endelse
             splog, "after ", i, " of ", n_elements(ii), " runs: ", n_elements(objlist) 
          endfor
          
          fxbhmake, hdr, 0
          sxaddpar, hdr, 'PROJECT', "HSC DC2 pointings"
          sxaddpar, hdr, 'RERUN', rerun
          sxaddpar, hdr, 'SWEEPDIR', getenv('PHOTO_SWEEP')
          sxaddpar, hdr, 'PTG_NAME', name[pointing]
          sxaddpar, hdr, 'PTG_RA',   ra[pointing]
          sxaddpar, hdr, 'PTG_DEC',  dec[pointing]
          sxaddpar, hdr, 'PTG_RAD',  rad[pointing]
          for r_i=0,n_elements(ii)-1 do begin
             sxaddpar, hdr, string('BSSRUN', r_i+1, format='(a,i0)'), runids[ii[r_i]], 'The BOSS runs we draw from'
          endfor
          
          mwrfits, objlist, string('miniSweeps-hsc_dc2-', pointing+1, '-', objtype, '.fits', $
                                   format='(a,i02,a,a0,a)'), hdr, alias=colnames, /create
       endfor
       fxbhmake, hdr, 0
       sxaddpar, hdr, 'PROJECT', "HSC DC2 pointings"
       sxaddpar, hdr, 'RERUN', rerun
       sxaddpar, hdr, 'SWEEPDIR', getenv('PHOTO_SWEEP')
       mwrfits, alltypeobjs, string('miniSweeps-hsc_dc2-', objtype, '.fits', $
                                    format='(a,a0,a)'), hdr, alias=colnames, /create
       allobjs = keyword_set(allobjs) ? [allobjs, alltypeobjs] : alltypeobjs 
    endfor
    
    mwrfits, allobjs, 'miniSweeps-hsc_dc2-all.fits', hdr, alias=colnames, /create
 end

