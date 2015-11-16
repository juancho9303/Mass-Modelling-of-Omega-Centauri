for (i=1; i <= 63; i+=1) {
  
    s1 = "s_"//i

  if (access (s1 // ".fits")) 
  {
      fxcor ("s_"//i//".fits",
"sdr_1.fits", apertures="*", cursor="", continuum="both", filter="none",
rebin="smallest", pixcorr=no, osample="*", rsample="*", apodize=0.2,
function="gaussian", width=INDEF, height=0., peak=no, minwidth=3.,
maxwidth=21., weights=1., background=0., window=INDEF, wincenter=INDEF,
output="sdr_1", verbose="long", imupdate=no, graphics="stdgraph", interactive=no,
autowrite=yes, autodraw=yes, ccftype="image", observatory="ctio",
continpars="", filtpars="", keywpars="")

}
;
}
