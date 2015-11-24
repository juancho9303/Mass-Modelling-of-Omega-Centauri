for (i=1; i <= 51; i+=1) {
  
    s1 = "w.cluster."//i

  if (access (s1 // ".fits")) 
  {
      fxcor ("w.cluster."//i//".fits",
"w.cluster.25.fits", apertures="*", cursor="", continuum="both", filter="none",
rebin="smallest", pixcorr=no, osample="*", rsample="*", apodize=0.2,
function="gaussian", width=INDEF, height=0., peak=no, minwidth=3.,
maxwidth=21., weights=1., background=0., window=INDEF, wincenter=INDEF,
output="vel", verbose="long", imupdate=no, graphics="stdgraph", interactive=no,
autowrite=yes, autodraw=yes, ccftype="image", observatory="ctio",
continpars="", filtpars="", keywpars="")

}
;
}
