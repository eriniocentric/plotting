#!/usr/local/bin/perl

# program to match spatial resolution of two images of a galaxy in different filters 
# first image J, second image IRAC1

use PDL;
use lib "$ENV{HOME}/Software";
use PDL::NiceSlice;
use Astro::WCS::LibWCS qw( :functions);
use code::MatchImages qw( matchplatescale dilate);
use PDL::Graphics::PGPLOT;
use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');
set_pgplot_options('CharSize' => 2.3,'HardCH'=> 1.8,'HardLW'=>1.2,'AspectRatio'=>1,'Colour'=>'black');

$datadir = "$ENV{HOME}/Software/SINGS/data.spitzer.caltech.edu/popular/sings/20070410_enhanced_v1/";
#path to IRAC and ancillary Halpha data
$datadir2mass = "$ENV{HOME}/Software/SINGS/irsa.ipac.caltech.edu/data/LGA/images/";
$outputdir = "$ENV{HOME}/Software/SINGS/outputimages/";

open IN, "$ENV{HOME}/Software/SINGS/object.list" or die "Error opening SINGS object.list file";
open OUT, ">$ENV{HOME}/Software/SINGS/output.dat" or die "Error openeing output.dat file";

while ( <IN>) {
@names = split;
$name = $names[0];
$name2mass = $names[1];
next if ($name2mass =~ 'NA');

$fitsimageJ = "$datadir2mass".$name2mass."/".$name2mass."_mosaic_j.fits";
$fitsimageL = "$datadir".$name."/IRAC/".$name."_v7.phot.1.fits";
$fitsimageM = "$datadir".$name."/IRAC/".$name."_v7.phot.2.fits";
$fitsimageHa = "$datadir".$name."/Ancillary/".$name."_HA_SUB_dr4.fits";

$outfitsJ = "/Users/mentuch/Software/SINGS/outputimages/".$name."J.fits";
$outfitsL = "/Users/mentuch/Software/SINGS/outputimages/".$name."L.fits";
$outfitsM = "/Users/mentuch/Software/SINGS/outputimages/".$name."M.fits";
$outfitsHa = "/Users/mentuch/Software/SINGS/outputimages/".$name."Ha.fits";

#match plate scale for images, output is 500 x 500 pixel
#plate scale matched with larger image

next unless (-e $fitsimageHa);
#next unless ($fitsimageHa =~ '6946');

#matchplatescale($fitsimageJ,$fitsimageL,$fitsimageM,$fitsimageHa,$name);
#get metallicity, SFR and other quantities from SINGS website info

open SAMPLE, "$ENV{HOME}/Software/SINGS/sample.dat";

while (<SAMPLE>) {
$radius = 200.; 
next unless m/$name/;

@v = split;
#($name2,$type,$nuc,$a,$dummy,$b,$vr,$dist,$Mopt,$LIRtoOPT,$LRIR,$logMH,$SFR,$OH) = pdl(@v); #this notation doesn't work
$name2 = $v[0]; $type = $v[1]; $nuc = $v[2]; $a = $v[3]; $b = $v[4]; $SFR = $v[12]; $OH = $v[13];

$radius = 60.*($a/2.) / 1.22;# in pixels

$itype = 0 if ($type =~ m/S/);
$itype = 1 if ($type =~ m/I/);
$itype = 2 if ($type =~ m/E/);
}
close SAMPLE;

#use adaptsmooth to ensure map of S/N > 10 across image

$outfitsJsmooth = "/Users/mentuch/Software/SINGS/outputimages/".$name."Jsm.fits";
$outfitsLsmooth = "/Users/mentuch/Software/SINGS/outputimages/".$name."Lsm.fits";
$outfitsMsmooth = "/Users/mentuch/Software/SINGS/outputimages/".$name."Msm.fits";
$filemaskJ = "/Users/mentuch/Software/SINGS/outputimages/".$name."maskJ.fits";
$filemaskL = "/Users/mentuch/Software/SINGS/outputimages/".$name."maskL.fits";
$filemaskM = "/Users/mentuch/Software/SINGS/outputimages/".$name."maskM.fits";

$Jzp = 1594.; $Jmagzp =20.775899;

$J = rfits $outfitsJ;
$sigma = $J->hdr->{SIGMA};

$sigma *= $Jzp * 10**(-0.4 * $Jmagzp) * 1e-6 * 3282.85 / ((2.777777845e-4)** 2 ); #image in MJy/sr

if (0) {
system "rm $outfitsJsmooth";
system "rm $outfitsLsmooth";
system "rm $outfitsMsmooth";
system "rm $filemaskJ";
system "rm $filemaskL";
system "rm $filemaskM";

system "/usr/local/bin/adaptsmooth/adaptsmooth -g -r $sigma -s 30.0 -G 1 $outfitsJ $outfitsJsmooth $filemaskJ";
system "/usr/local/bin/adaptsmooth/adaptsmooth -g -r 0.07 -s 30.0 -G 1 $outfitsL $outfitsLsmooth $filemaskL";
system "/usr/local/bin/adaptsmooth/adaptsmooth -g -r 0.07 -s 30.0 -G 1 $outfitsM $outfitsMsmooth $filemaskM";
}
$J = rfits $outfitsJsmooth;
$L = rfits $outfitsLsmooth;
$Ha = rfits $outfitsHa;
$M = rfits $outfitsMsmooth;

$maskJ = rfits $filemaskJ;
$maskL = rfits $filemaskL;
$maskM = rfits $filemaskM;



#determine mask from 2MASS star catalog, only use point sources which have subtracted HA>3
# over 3 pixel radius aperture (indicative of true PSs due to mismatching PSFs in subtracting)

($xstar,$ystar) = rcols("$ENV{HOME}/Software/SINGS/starcat/".$name.".cat");

$i = which($xstar >= 0 & $xstar < 500 & $ystar >= 0 & $ystar < 500);
 
$xstar = floor($xstar($i));
$ystar = floor($ystar($i));
$n = nelem($i);

$starmask = ones(500,500);
$sum = zeroes($n);
$i = 0;
for ($i == 0; $i < $n; $i++) {
  $x = $xstar(($i));
  $y = $ystar(($i));

  $xs = max(pdl(0,$x-2));
  $xe = min(pdl(499,$x+2));

  $ys = max(pdl(0,$y-2));
  $ye = min(pdl(499,$y+2));
  $sub = $Ha($xs:$xe,$ys:$ye)->clump(-1);

  $sum($i) .= sum($Ha($xs:$xe,$ys:$ye));
  #only use point sources that have large negative pixel values as a 
  #result of mismatching PSFs in subtracting
  if (any($sub < -0.05) ) {
    $starmask($xs:$xe,$ys:$ye) .= 0;
  }
}

$maskJgood = $maskJ < 100;
$maskLgood = $maskL < 100;
$maskMgood = $maskM < 100;

$maskcircle = zeroes(500,500);
$maskcircle(250-$radius:250+$radius,250-$radius:250+$radius) .= 1;

$diff = $L-$J;
$div = $L/$J;

$diff2 = $M-$J;
$div2 = $M/$J;

$diff->sethdr( $J->hdr);
$div->sethdr( $J->hdr);
$diff2->sethdr( $J->hdr);
$div2->sethdr( $J->hdr);

wfits $diff, "$outputdir".$name."dif.fits";
wfits $div, "$outputdir".$name."div.fits";
wfits $diff2, "$outputdir".$name."dif2.fits";
wfits $div2, "$outputdir".$name."div2.fits";

$dbin = 0.03;
$start = -1.0; 
$end = 3.;
$numbins = int(($end - $start) / $dbin);
$hist = histogram($div->clump(-1), $dbin, $start, $numbins);

$plotbins = sequence($numbins)*$dbin + $start;

$med = median($diff);

#Take median of four patches on Ha_SUB image to get noise

$N1 = median($Ha(0:10,0:10));
$n2 = median($Ha(30:40,10:20));
$n3 = median($Ha(450:460,480:490));
$n4 = median($Ha(372:382,296:306));
$n5 = median($Ha(70:80,410:420));
$n6 = median($Ha(100:110,100:110));

$sigmaHA = abs(median(pdl($n1,$n2,$n3,$n4,$n5,$n6)));
$sigmaHA = 0.08 if ($sigmaHA <= 0);
$maskraw = $div > 0.9 & $div < 1e4;

#dilate mask by 3 pixels

if (1) {
$mask = dilate($maskraw,6);

wfits $maskraw, "/Users/mentuch/Software/SINGS/outputimages/".$name."rawmask.fits";
wfits $mask,"/Users/mentuch/Software/SINGS/outputimages/".$name."mask.fits";
} 
else {
$maskraw = rfits "/Users/mentuch/Software/SINGS/outputimages/".$name."rawmask.fits";
$mask = rfits "/Users/mentuch/Software/SINGS/outputimages/".$name."mask.fits";
}
#dev '/xs', 5, 1,{AspectRatio=>0.2,WindowWidth=>15};
dev "$ENV{HOME}/Software/SINGS/postagestamps/$name.ps/vcps", 5, 1,{AspectRatio=>0.2};
imag -log10(abs($J*$maskJgood)),-1,1,{Axis=>'Empty',Title=>'J','Color'=>black};
imag -log10(abs($L*$maskLgood)),-1,1,{Axis=>'Empty',Title=>'[3.6]'};
imag -log10(abs( $div*$maskJgood*$maskLgood)), -0.5,0.5,{Axis=>'Empty',Title=>'J-[3.6]'};
imag -log10(abs($Ha)),-log10($sigmaHA*100),-log10($sigmaHA*10),{Axis=>'Empty',Title=>'H\ga'};
hold;
cont $mask*$maskJgood*$maskLgood,{Color=>red,Ncontours=>2};
release;

$mipsfile = "$datadir".$name."/MIPS/".$name."_mips24_image_v5-0.fits";
if (-e $mipsfile) {
  $mipsremap = "$ENV{HOME}/Software/SINGS/outputimages/".$name."_mips24_remap.fits";
  ($ra,$dec) = split(' ', qx (nedpos -d $name));
  $pix = 500;
  $outputplatescale = 1.22328355; #arcsec/pix for IRAC Ch1
  system "remap -j $ra $dec -p $outputplatescale -y $pix $pix -o $mipsremap $mipsfile";
  $mips = rfits($mipsremap);
  imag -log10(abs($mips)),-2,0.5,{Axis=>'Empty', Title=> 'MIPS 24 micron'};
  hold;
  cont $mask*$maskJgood*$maskLgood,{Color=>red,Ncontours=>2};
  release;
}

close_window();
#system "convert $ENV{HOME}/Software/SINGS/postagestamps/$name.ps $ENV{HOME}/Software/SINGS/postagestamps/$name.jpg";
$Ha_mask = $Ha*$mask;
$Ha_mask->sethdr( $J->hdr);
wfits $Ha_mask, "$outputdir".$name."_Ha_mask.fits";

#add up flux in Ha and difference image intensity over mask values

$Ha_maskflux = sum($Ha*$maskraw*$maskJgood*$maskLgood*$starmask);
$J_maskflux = sum($J*$maskraw*$maskJgood*$maskLgood*$starmask);
$L_maskflux = sum($L*$maskraw*$maskJgood*$maskLgood*$starmask);
$mips_maskflux = sum($mips*$maskraw*$maskJgood*$maskLgood*$starmask);
$diff_maskflux2 = sum($diff*$maskraw*$maskJgood*$maskLgood*$starmask);
$nummask = sum($maskraw*$maskJgood*$maskLgood*$starmask);

$diff_maskflux = $L_maskflux - $J_maskflux;
$div_maskflux = $L_maskflux / $J_maskflux;

print "DIFF MASKS : $diff_maskflux and $diff_maskflux2 \n";
#measure total Ha and total HA in mask

$masknoisesub = $Ha > 5*$sigmaHA;

$Ha_insidemaskflux = sum($Ha*$masknoisesub*$mask*$maskJgood*$maskLgood*$starmask);
$Ha_flux = sum($Ha*$masknoisesub*$maskJgood*$maskLgood*$starmask);

#measure total mips flux and mips flux in mask
if (-e $mipsfile) {

  $masknoisesub = $mips > 0.1;
  $mips_insidemaskflux = sum($mips*$masknoisesub*$mask*$maskJgood*$maskLgood*$starmask);
  $mips_flux = sum($mips*$masknoisesub*$maskJgood*$maskLgood*$starmask);
}

#get metallicity, SFR and other quantities from SINGS website info

open SAMPLE, "$ENV{HOME}/Software/SINGS/sample.dat";

while (<SAMPLE>) {
next unless m/$name/;

@v = split;
#($name2,$type,$nuc,$a,$dummy,$b,$vr,$dist,$Mopt,$LIRtoOPT,$LRIR,$logMH,$SFR,$OH) = pdl(@v); #this notation doesn't work
$name2 = $v[0]; $type = $v[1]; $nuc = $v[2]; $SFR = $v[12]; $OH = $v[13];
print "$type $SFR \n";
$itype = 0 if ($type =~ m/S/);
$itype = 1 if ($type =~ m/I/);
$itype = 2 if ($type =~ m/E/);
}
close SAMPLE;

printf OUT "%10s  %3.0f  %8.6f  %8.6f  %6.0f  %8.4f   %10.6g  %8.4f  %8.4f  %8.4f  %8.4f  %10.6g  %10.6g\n", $name, $itype, $SFR, $OH, $nummask, $Ha_maskflux,$mips_maskflux,$diff_maskflux,$div_maskflux,$Ha_insidemaskflux,$Ha_flux,$mips_insidemaskflux,$mips_flux;
}
close IN;
close OUT;

($name, $itype, $SFR, $OH, $nummask, $Ha_maskflux,$mips_maskflux,$diff_maskflux,$div_maskflux,$Ha_insidemaskflux,$Ha_flux,$mips_insidemaskflux,$mips_flux) = rcols "Software/SINGS/output.dat";

$Ha_maskflux /= $nummask;
$mips_maskflux /= $nummask;
$diff_maskflux /= $nummask;

badmask($mips_maskflux,-99);
dev 'Software/SINGS/Ha_mips.ps/vcps', 1, 1, {AspectRatio=>1,Charsize=>1.4};
env -0.3,2.2,-2,0.7, {Axis=>'logXY',XTitle=>'masksum(mips)/npix',YTitle=>'masksum(Ha)/npix'};
$sel = which($mips_maskflux > 0 & $Ha_maskflux > 0 & $diff_maskflux > 0);
points log10($mips_maskflux($sel)),log10($Ha_maskflux($sel));
close_window;

dev 'Software/SINGS/J-L_mips.ps/vcps', 1, 1, {AspectRatio=>1,Charsize=>1.4};
env -0.3,2.2,-2.,0.7, {Axis=>'logXY',XTitle=>'masksum(mips)/npix',YTitle=>'masksum([J-L])/npix'};;
points log10($mips_maskflux($sel)),log10($diff_maskflux($sel));
close_window;

dev 'Software/SINGS/J-L_Ha.ps/vcps', 1, 1, {AspectRatio=>1,Charsize=>1.4};
env -2.4,0.7,-2,0.7, {Axis=>'logXY',XTitle=>'masksum(Ha)/npix',YTitle=>'masksum([J-L])/npix'};;
points log10($Ha_maskflux($sel)),log10($diff_maskflux($sel));
close_window;

$i = which($Ha_flux > 0 & $diff_maskflux > 0  );
dev 'Software/SINGS/Hamaskfraction_J-L.ps/vcps', 1, 1, {AspectRatio=>1,Charsize=>1.4};
points log10($diff_maskflux($i)), $Ha_insidemaskflux($i)/$Ha_flux($i),{XRange=>[-2.7,0.5],YRange=>[-0.1,1.1],YTitle=>'Ha Flux in mask/ Total Flux',XTitle=>'J-L Flux in mask/npix'};
hold;
points log10($diff_maskflux($i)), $Ha_insidemaskflux($i)/$Ha_flux($i),{Sympol=>circle};
close_window;

$ratio = $Ha_insidemaskflux($i)/$Ha_flux($i);
dev 'Software/SINGS/HAmask_hist.ps/vcps';
bin hist($ratio,0,1.1,0.17),{YRange=>[0,30],XTitle=>'Fraction of Halpha Flux in Mask',YTitle=>'n',Title=>'SINGS Sample = 54 Galaxies',Charsize=>1.4,XRange=>[0,1.05]};
close_window;

$i = which($mips_flux > 0);
dev 'Software/SINGS/Mipsmaskfraction_vs_MIPSflux.ps/vcps',1,1,{AspectRatio=>1};
points log10($mips_flux($i)),$mips_insidemaskflux($i)/$mips_flux($i),{XRange=>[3,6],YRange=>[-0.1,0.8],Axis=>'logX',XTitle=>'Total MIPS flux',YTitle=>'MIPS flux in mask/ Total MIPS flux',Charsize=>1.4};
close_window;

$ratio2 = $mips_insidemaskflux($i)/$mips_flux($i);

dev 'Software/SINGS/Mipsmask_hist.ps/vcps';
bin hist($ratio2,0,1.1,0.17),{YRange=>[0,10],XTitle=>'Fraction of MIPS Flux in Mask',YTitle=>'n',Title=>'SINGS Sample = 30 Galaxies',Charsize=>1.4,XRange=>[0,1.05]};
close_window;
