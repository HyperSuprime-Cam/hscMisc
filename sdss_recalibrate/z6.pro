pro z6, filename
    setenv, 'PHOTO_RESOLVE=/u/dss/resolve/2010-05-23'
    setenv, 'PHOTO_SWEEP=/u/dss/sweeps/dr9'
    rerun = 301
    outdir = '/scr/fahl2/users/cloomis/astrom/dr9.full.2'

    if keyword_set(filename) then begin 
       readcol, filename, runs, format='L'
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
       mash_sweeps, run, rerun=rerun, outdir=outdir
    endfor
 end
