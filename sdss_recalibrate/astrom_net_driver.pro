; 
pro astrom_net_driver, runsfile=runsfile, recalibrate=recalibrate, $
                       rerun=rerun, outdir=outdir
    setenv, 'PHOTO_RESOLVE=/u/dss/resolve/2010-05-23'
    setenv, 'PHOTO_SWEEP=/u/dss/sweeps/dr9'
    
    rerun = keyword_set(rerun) ? rerun : 301
    outdir = keyword_set(outdir) ? outdir : '/scr/fahl2/users/cloomis/astrom/scratch'

    if keyword_set(runsfile) then begin 
       readcol, runsfile, runs, format='L'
    end else begin
       window_read, flist=flist
       all_runs = flist.run
       runs = all_runs[uniq(all_runs, sort(all_runs))]
    end

    ; Got this far...
    ;w = where(runs GE 3842)
    ;runs = runs[w]

    print,"processing runs: ", runs
    
    for i=0,n_elements(runs)-1 do begin
       run = runs[i]
       mash_sweeps, run, rerun=rerun, outdir=outdir, $
                    recalibrate=recalibrate
    endfor
 end
