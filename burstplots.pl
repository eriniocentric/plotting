use PDL; 
use PDL::NiceSlice; 
use PDL::Graphics::PGPLOT::Window; 
use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');
set_pgplot_options('CharSize' => 0.8,'HardCH'=> 1.0,'HardLW'=>2.5);

chdir "$ENV{HOME}/Software/Masses/";

system "sed s/SA// <GDDSVIzK-masses.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.dat";
($id1,$z1,$K1,$massVIzK,$massVIzK_err,$ageVIzK,$chi2VIzK,$tVIzK,$tburstVIzK) = rcols("$ENV{HOME}/Software/Masses/temp.dat",0,1,2,3,4,18,20,22,24);

system "sed s/SA// <GDDSVIzK_nb-masses.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.dat";
($id1_nb,$z1_nb,$K1_nb,$massVIzK_nb,$massVIzK_err_nb,$ageVIzK_nb,$chi2VIzK_nb,$tVIzK_nb,$tburstVIzK_nb) = rcols("$ENV{HOME}/Software/Masses/temp.dat",0,1,2,3,4,18,20,22,24);

system "sed s/SA// <GDDSVIzK12-masses.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.dat";
($id2,$z2,$K2,$massVIzK12,$massVIzK12_err,$fburstVIzK12,$ageVIzK12,$chi2VIzK12,$tVIzK12,$tburstVIzK12) = rcols("$ENV{HOME}/Software/Masses/temp.dat",0,1,2,3,4,14,18,20,22,24);

system "sed s/SA// <GDDSVIzK12_nb-masses.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.dat";
($id2_nb,$z2_nb,$K2_nb,$massVIzK12_nb,$massVIzK12_err_nb,$ageVIzK12_nb,$chi2VIzK12_nb,$tVIzK12_nb,$tburstVIzK12_nb) = rcols("$ENV{HOME}/Software/Masses/temp.dat",0,1,2,3,4,18,20,22,24);

system "sed s/SA// <GDDSVIzK1234-masses.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.dat";
($id3,$z3,$K3,$massVIzK1234,$massVIzK1234_err,$ageVIzK1234,$chi2VIzK1234,$tVIzK1234,$tburstVIzK1234) = rcols("$ENV{HOME}/Software/Masses/temp.dat",0,1,2,3,4,18,20,22,24);

system "sed s/SA// <GDDSVIzK1234_nb-masses.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.dat";
($id3_nb,$z3_nb,$K3_nb,$massVIzK1234_nb,$massVIzK1234_err_nb,$ageVIzK1234_nb,$chi2VIzK1234_nb,$tVIzK1234_nb,$tburstVIzK1234_nb) = rcols("$ENV{HOME}/Software/Masses/temp.dat",0,1,2,3,4,18,20,22,24);

system "sed s/SA// <~/Software/PAH/pahratio.txt> tmp.txt";
system "sed s/-// <tmp.txt> temp.txt";
($idpah_,$zpah_,$confz_,$fratio3_,$fratio4_,$sfr2000_,$sfrOII_,$restUB_) = rcols("temp.txt",0,1,2,3,4,5,6,7);

$idpah = zeroes(nelem($id1));
$zpah  = zeroes(nelem($id1));
$fratio3 = zeroes(nelem($id1));
$fratio4 = zeroes(nelem($id1));
$sfr2000 = zeroes(nelem($id1));
$sfrOII  = zeroes(nelem($id1));
$restUB  = zeroes(nelem($id1));
$conf   = zeroes(nelem($id1));

for ($i=0;$i<=nelem($id1)-1;$i++) {
    if (any $idpah_ == $id1(($i))) {
	$tmp = which($idpah_ == $id1(($i)));
	$idpah($i) .= $idpah_(($tmp(0)));	
	$zpah($i)  .= $zpah_(($tmp(0)));
	$fratio3($i) .= $fratio3_(($tmp(0)));
	$fratio4($i) .= $fratio4_(($tmp(0)));
	$sfr2000($i) .= $sfr2000_(($tmp(0)));
	$sfrOII($i)  .= $sfrOII_(($tmp(0)));
	$restUB($i) .= $restUB_(($tmp(0)));
	$conf($i) .= $confz_(($tmp(0)));
    }
    else {
       	$idpah($i) .= -9999;
	$zpah($i)  .= -9999;
	$fratio3($i) .= -9999;
	$fratio4($i) .= -9999;
	$sfr2000($i) .= -9999;
	$sfrOII($i)  .= -9999;
	$restUB($i) .= -9999;
	$conf($i) .= -9999;
    }
}

$idxUB1 = which($restUB > 0.1  & $conf>1 & $conf<10);
$idxUB2 = which($restUB <= 0.1 & $restUB >= -0.1  & $conf>1 & $conf<10);
$idxUB3 = which($restUB < -0.1 & $conf>1 & $conf<10);

$opt = {Device => 'figures/VIzK1234_vs_VIzK1234_nb.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($massVIzK1234_nb($idxUB1),$massVIzK1234($idxUB1), {Colour=>red, Symbol=>circle, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($massVIzK1234_nb($idxUB2),$massVIzK1234($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($massVIzK1234_nb($idxUB3),$massVIzK1234($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(20),sequence(20),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 10.5,9.3,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK1234 Mass (Without Burst) ','VIzK1234 Mass (With Burst)');
$win->release; 
$win->close();

$opt = {Device => 'figures/VIzK12_vs_VIzK12_nb.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($massVIzK12_nb($idxUB1),$massVIzK12($idxUB1), {Colour=>red, Symbol=>circle, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($massVIzK12_nb($idxUB2),$massVIzK12($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($massVIzK12_nb($idxUB3),$massVIzK12($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(20),sequence(20),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 10.5,9.3,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK12 Mass (Without Burst) ','VIzK12 Mass (With Burst)');
$win->release; 
$win->close();

$opt = {Device => 'figures/ageVIzK1234_vs_ageVIzK1234_nb.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($ageVIzK1234_nb($idxUB1),$ageVIzK1234($idxUB1), {Colour=>red, Symbol=>circle}); #, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($ageVIzK1234_nb($idxUB2),$ageVIzK1234($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($ageVIzK1234_nb($idxUB3),$ageVIzK1234($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(5000),sequence(5000),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 1500,500,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK1234 Age (Without Burst) ','VIzK1234 Age (With Burst)');
$win->release; 
$win->close();

$opt3 = {Device => 'figures/VIzK_vs_VIzK_nb.ps/cps'};
$win3 = PDL::Graphics::PGPLOT::Window->new($opt3);
$win3->points($massVIzK_nb($idxUB1),$massVIzK($idxUB1), {Colour=>red, Symbol=>circle, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win3->hold;
$win3->points($massVIzK_nb($idxUB2),$massVIzK($idxUB2), {Colour=>orange,Symbol=>square});
$win3->points($massVIzK_nb($idxUB3),$massVIzK($idxUB3), {Colour=>blue,Symbol=>star});
$win3->line(sequence(20),sequence(20),{colour=>red});
$win3->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 10.5,9.3,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win3->label_axes('VIzK Mass (Without Burst) ','VIzK Mass (With Burst)');
$win3->release; 
$win3->close();

$opt = {Device => 'figures/ageVIzK_vs_ageVIzK_nb.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($ageVIzK_nb($idxUB1),$ageVIzK($idxUB1), {Colour=>red, Symbol=>circle}); #, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($ageVIzK_nb($idxUB2),$ageVIzK($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($ageVIzK_nb($idxUB3),$ageVIzK($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(5000),sequence(5000),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 1500,1000,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK Age (Without Burst) ','VIzK Age (With Burst)');
$win->release; 
$win->close();

$opt = {Device => 'figures/VIzK1234_vs_VIzK_nb.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($massVIzK_nb($idxUB1),$massVIzK1234($idxUB1), {Colour=>red, Symbol=>circle, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($massVIzK_nb($idxUB2),$massVIzK1234($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($massVIzK_nb($idxUB3),$massVIzK1234($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(20),sequence(20),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 10.5,9.3,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK Mass (Without Burst) ','VIzK1234 Mass (With Burst)');
$win->release; 
$win->close();

$opt = {Device => 'figures/VIzK1234_vs_VIzK.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($massVIzK($idxUB1),$massVIzK1234($idxUB1), {Colour=>red, Symbol=>circle, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($massVIzK($idxUB2),$massVIzK1234($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($massVIzK($idxUB3),$massVIzK1234($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(20),sequence(20),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 10.5,9.3,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK Mass (With Burst) ','VIzK1234 Mass (With Burst)');
$win->release; 
$win->close();

$opt = {Device => 'figures/VIzK12_vs_VIzK.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($massVIzK($idxUB1),$massVIzK12($idxUB1), {Colour=>red, Symbol=>circle, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($massVIzK($idxUB2),$massVIzK12($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($massVIzK($idxUB3),$massVIzK12($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(20),sequence(20),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 10.5,9.3,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK Mass (With Burst) ','VIzK12 Mass (With Burst)');
$win->release; 
$win->close();
 

$opt = {Device => 'figures/VIzK12_vs_VIzK1234.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($massVIzK1234($idxUB1),$massVIzK12($idxUB1), {Colour=>red, Symbol=>circle, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($massVIzK1234($idxUB2),$massVIzK12($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($massVIzK1234($idxUB3),$massVIzK12($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(20),sequence(20),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 10.5,9.3,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK1234 Mass (With Burst) ','VIzK12 Mass (With Burst)');
$win->release; 
$win->close();

$opt = {Device => 'figures/ageVIzK1234_vs_ageVIzK.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($ageVIzK($idxUB1),$ageVIzK1234($idxUB1), {Colour=>red, Symbol=>circle}); #, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($ageVIzK($idxUB2),$ageVIzK1234($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($ageVIzK($idxUB3),$ageVIzK1234($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(5000),sequence(5000),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 2000,1000,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK Age (With Burst) ','VIzK1234 Age (With Burst)');
$win->release; 
$win->close();

#compare to maraston models that have no burst models but are csp

system "sed s/SA// <~/Software/maraston/masses/masses_VIzK.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.dat";
($id1_mara,$z1_mara,$K1_mara,$massVIzK_mara,$massVIzK_err_mara,$ageVIzK_mara) = rcols("$ENV{HOME}/Software/Masses/temp.dat",0,1,2,3,4,18);

system "sed s/SA// <~/Software/maraston/masses/masses_VIzK12.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.dat";
($id12_mara,$z12_mara,$K12_mara,$massVIzK12_mara,$massVIzK12_err_mara,$ageVIzK12_mara) = rcols("$ENV{HOME}/Software/Masses/temp.dat",0,1,2,3,4,18);

system "sed s/SA// <~/Software/maraston/masses/masses_VIzK123.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.dat";
($id2_mara,$z2_mara,$K2_mara,$massVIzK123_mara,$massVIzK123_err_mara,$ageVIzK123_mara) = rcols("$ENV{HOME}/Software/Masses/temp.dat",0,1,2,3,4,18);

system "sed s/SA// <~/Software/maraston/masses/masses_VIzK1234.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.dat";
($id3_mara,$z3_mara,$K3_mara,$massVIzK1234_mara,$massVIzK1234_err_mara,$ageVIzK1234_mara) = rcols("$ENV{HOME}/Software/Masses/temp.dat",0,1,2,3,4,18);

$ix = zeroes(nelem($id1));

for ($i=0;$i<=nelem($id1)-1;$i++) {
	$tmp = which($id1_mara == $id1(($i))); 
	$ix(($i)) .= $tmp(0);	
}

$id1_mara = $id1_mara($ix); $z1_mara = $z1_mara($ix); $K1_mara = $K1_mara($ix); $massVIzK_mara = $massVIzK_mara($ix); $massVIzK_err_mara = $massVIzK_err_mara($ix); $ageVIzK_mara = $ageVIzK_mara($ix);

$id2_mara = $id2_mara($ix); $z2_mara = $z2_mara($ix); $K2_mara = $K2_mara($ix); $massVIzK123_mara = $massVIzK123_mara($ix); $massVIzK123_err_mara = $massVIzK123_err_mara($ix); $ageVIzK123_mara = $ageVIzK123_mara($ix);

$id3_mara = $id3_mara($ix); $z3_mara = $z3_mara($ix); $K3_mara = $K3_mara($ix); $massVIzK1234_mara = $massVIzK1234_mara($ix); $massVIzK1234_err_mara = $massVIzK1234_err_mara($ix); $ageVIzK1234_mara = $ageVIzK1234_mara($ix);


$opt = {Device => 'figures/VIzK1234_mara_vs_VIzK1234_nb.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($massVIzK1234_nb($idxUB1),$massVIzK1234_mara($idxUB1), {Colour=>red, Symbol=>circle, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($massVIzK1234_nb($idxUB2),$massVIzK1234_mara($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($massVIzK1234_nb($idxUB3),$massVIzK1234_mara($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(20),sequence(20),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 10.5,9.3,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK1234 Mass (PEGASE Without Burst) ','VIzK1234 Mass (Maraston)');
$win->release; 
$win->close();

$opt = {Device => 'figures/VIzK_mara_vs_VIzK_nb.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($massVIzK_nb($idxUB1),$massVIzK_mara($idxUB1), {Colour=>red, Symbol=>circle, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($massVIzK_nb($idxUB2),$massVIzK_mara($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($massVIzK_nb($idxUB3),$massVIzK_mara($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(20),sequence(20),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 10.5,9.3,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK Mass (PEGASE Without Burst) ','VIzK Mass (Maraston)');
$win->release; 
$win->close();

$opt = {Device => 'figures/VIzK1234_mara_vs_VIzK_mara.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($massVIzK_mara($idxUB1),$massVIzK1234_mara($idxUB1), {Colour=>red, Symbol=>circle, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($massVIzK_mara($idxUB2),$massVIzK1234_mara($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($massVIzK_mara($idxUB3),$massVIzK1234_mara($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(20),sequence(20),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 10.5,9.3,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK Mass (Maraston) ','VIzK1234 Mass (Maraston)');
$win->release; 
$win->close();

$opt = {Device => 'figures/ageVIzK1234_mara_vs_ageVIzK1234_nb.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($ageVIzK1234_nb($idxUB1),$ageVIzK1234_mara($idxUB1), {Colour=>red, Symbol=>circle}); #, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($ageVIzK1234_nb($idxUB2),$ageVIzK1234_mara($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($ageVIzK1234_nb($idxUB3),$ageVIzK1234_mara($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(5000),sequence(5000),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 1500,1000,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK1234 Age (PEGASE Without Burst) ','VIzK1234 Age (Maraston)');
$win->release; 
$win->close();

$opt = {Device => 'figures/VIzK12_mara_vs_VIzK12_nb.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($massVIzK12_nb($idxUB1),$massVIzK12_mara($idxUB1), {Colour=>red, Symbol=>circle, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($massVIzK12_nb($idxUB2),$massVIzK12_mara($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($massVIzK12_nb($idxUB3),$massVIzK12_mara($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(5000),sequence(5000),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 1500,1000,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK12 Mass (PEGASE Without Burst) ','VIzK12 Mass (Maraston)');
$win->release; 
$win->close();

$opt = {Device => 'figures/ageVIzK_mara_vs_ageVIzK_nb.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($ageVIzK_nb($idxUB1),$ageVIzK_mara($idxUB1), {Colour=>red, Symbol=>circle}); #, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($ageVIzK_nb($idxUB2),$ageVIzK_mara($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($ageVIzK_nb($idxUB3),$ageVIzK_mara($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(5000),sequence(5000),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 1500,1000,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK Age (PEGASE Without Burst) ','VIzK Age (Maraston)');
$win->release; 
$win->close();

$opt = {Device => 'figures/ageVIzK1234_mara_vs_ageVIzK_mara.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($ageVIzK_mara($idxUB1),$ageVIzK1234_mara($idxUB1), {Colour=>red, Symbol=>circle}); #, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win->hold;
$win->points($ageVIzK_mara($idxUB2),$ageVIzK1234_mara($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($ageVIzK_mara($idxUB3),$ageVIzK1234_mara($idxUB3), {Colour=>blue,Symbol=>star});
$win->line(sequence(5000),sequence(5000),{colour=>red});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 1500,1000,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('VIzK Age (Maraston) ','VIzK1234 Age (Maraston)');
$win->release; 
$win->close();

$massdif = $massVIzK12-$massVIzK;

$opt = {Device =>'figures/massdif_vs_pahratio.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$idx = which(abs($zpah-1.4) < 0.4 & $conf > 1 & $conf <5);
$win->points($massdif($idx),$fratio4($idx),{Col=>red,Symbol=>circle});
$win->hold;
$idx = which(abs($zpah-0.78) < 0.13 & $conf > 1 & $conf <5);
$win->points($massdif($idx),$fratio3($idx),{Col=>blue,Symbol=>star});
$win->label_axes('Mass (VIzK12) - Mass (VIzK)','Flux Ratio');
$win->release; 
$win->close();

$massdif = $massVIzK12-$massVIzK;

$opt = {Device => 'figures/massdif_vs_agem.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($ageVIzK12($idxUB1),$massdif($idxUB1), {Colour=>red, Symbol=>circle});
$win->hold;
$win->points($ageVIzK12($idxUB2),$massdif($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($ageVIzK12($idxUB3),$massdif($idxUB3), {Colour=>blue,Symbol=>star});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 500,0.08,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('Mass Weighted Age from VIzK12 SED Fitting','VIzK12 Mass - VIzK Mass');
$win->release; 
$win->close();

$opt = {Device => 'figures/massdif_vs_age.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($tVIzK12($idxUB1),$massdif($idxUB1), {Colour=>red, Symbol=>circle});
$win->hold;
$win->points($tVIzK12($idxUB2),$massdif($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($tVIzK12($idxUB3),$massdif($idxUB3), {Colour=>blue,Symbol=>star});
$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 500,0.08,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
#$win->label_axes('Age from VIzK12 SED Fitting','VIzK12 Mass - VIzK Mass');
$win->release; 
$win->close();

$timesinceburst = $tVIzK12-$tburstVIzK12;

$idxUB1 = which($restUB > 0.1  & $conf>1 & $conf<10 & $fburstVIzK12>-3);
$idxUB2 = which($restUB <= 0.1 & $restUB >= -0.1  & $conf>1 & $conf<10 & $fburstVIzK12>-3);
$idxUB3 = which($restUB < -0.1 & $conf>1 & $conf<10 & $fburstVIzK12>-3);

$opt = {Device => 'figures/massdif_vs_tburst.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($timesinceburst($idxUB1),$massdif($idxUB1), {Colour=>red, Symbol=>circle});
$win->hold;
$win->points($timesinceburst($idxUB2),$massdif($idxUB2), {Colour=>orange,Symbol=>square});
$win->points($timesinceburst($idxUB3),$massdif($idxUB3), {Colour=>blue,Symbol=>star});
#$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 500,0.08,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('Time since burst from VIzK12 SED Fitting','VIzK12 Mass - VIzK Mass');
$win->release; 
$win->close();


$opt = {Device =>'figures/fratio_vs_tburst.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$idx = which(abs($zpah-1.4) < 0.4 & $conf > 1 & $conf <5 & $fburstVIzK12>-3);
$win->points($timesinceburst($idx),$fratio4($idx),{Col=>red});
$win->hold;
$idx = which(abs($zpah-0.78) < 0.13 & $conf > 1 & $conf <5 & $burstVIzK12>-3);
$win->points($timesinceburst($idx),$fratio3($idx),{Col=>blue});
$win->label_axes('Time since Burst (Myr)','3.3 micron Rest Frame Flux Ratio');
$win->release; 
$win->close();

$opt = {Device => 'figures/tburst_vs_ageM.ps/cps'};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->points($ageVIzK12($idxUB1),$timesinceburst($idxUB1),{Colour=>red, Symbol=>circle,XRange=>[0,6000],YRange=>[0,1300]});
$win->hold;
$win->points($ageVIzK12($idxUB2),$timesinceburst($idxUB2),{Colour=>orange,Symbol=>square});
$win->points($ageVIzK12($idxUB3),$timesinceburst($idxUB3),{Colour=>blue,Symbol=>star});
#$win->legend(['UB rest > 0.1','0.1 > UB rest > -0.1','UB rest < -0.1'], 500,0.08,{Colour=>[red,orange,blue],Symbol=>[circle,square,star]});
$win->label_axes('Mass weighted age (Myr)','Time since burst (Myr)');
$win->release; 
$win->close();
