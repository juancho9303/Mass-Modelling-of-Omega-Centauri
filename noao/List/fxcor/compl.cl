procedure completitud_cl()

int p,N

begin

noao
digiphot
daophot

unlearn daophot

datapars.scale=1.0
datapars.fwhmpsf=2.077
datapars.sigma=0.075
datapars.readnoise=5.25
datapars.epadu=2.0
centerpars.calgorithm="centroid"
fitskypars.salgorithm="mode"
fitskypars.annulus=4.0
fitskypars.dannulus=3.5
photpars.apertures=4.0
photpars.zmag=26.414
daopars.function=auto
daopars.psfrad=4.0
daopars.fitrad=20.0
daopars.recenter=yes
daopars.fitsky=no

!python randc.py
#!bl < mksf606w.bl

N = 3


for (p=1;p<N;p+=1) {
print("mksynth Random/rand"//p//".dat 1.0 syntf606w"//p//".fits rezx=5200.0 rezy=5200.0 random=NO rdnoise=5.25 backgr=0.075 epadu=111.5 zpoint=4.1e12 bias=0 photnoise=NO psftype=USER psffile=m4f606w.psf.fits", >>&"mksf606w.bl")

!bl < mksf606w.bl

imarith (operand1="syntf606w"//p//".fits", op="+", operand2="/home/daniel/Documentos/Estudio/M/ngc6121/New/ngc6121_F606W.fits", result="syntf606w"//p//".fits")

print "imagen "//p//""
print "_______________Daofind________________"
daofind (image="syntf606w"//p//".fits", output="default", verify-, interactive-, verbose-)

copy (input="syntf606w"//p//".fits.coo.1", output="Coord/syntf606w"//p//".coo")

!rm syntf606w*.fits.coo.1

print "_______________Phot___________________"
phot (image="syntf606w"//p//".fits", coords="Coord/syntf606w"//p//".coo", output="Mag/default", verify-, update-, interactive-, verbose-)


print "__________________Allstar________________"
allstar (image="syntf606w"//p//".fits", photfile="Mag/default", psfimage="ngc6121_F606W.fits.psf.1.fits", allstarfile="Als/default", rejfile="Als/default", subimage="Als/default", verify-, update-, verbose-)

}

#for (j=1;j<N;j+=1) {
#tmatch (input1="rand"//j//".dat", input2="syntf606w"//j//".fits.mag.1", output="match.mag"//#j//".dat", match1="c1,c2", match2="c7,c8", maxnorm=1.0, incol1="c1,c2,c3", incol2="c30,c31")
#}



#for (j=1;j<N;j+=1) {
#tmatch (input1="rand"//j//".dat", input2="syntf606w"//j//".fits.als.1", output="match.als"//#j//".dat", match1="c1,c2", match2="c2,c3", maxnorm=1.0, incol1="c1,c2,c3", incol2="c4,c5")
#}

#!rm -rf syntf606*.fits

#operand2="/home/esteban/Dropbox/UdeA/Students/Second_pop/ori_ima/HST_10775_64_ACS_WFC_F606W_drz.fits[1]"

end
